#!/usr/bin/env python
__version__ = '0.1.0'

import os
import sys
import csv
import time
import string
import subprocess as subP
from itertools import product
from copy import deepcopy
from multiprocessing import Pool, cpu_count, Value, Manager
from math import log, cosh, exp
from optparse import OptionParser, OptionGroup
# dependencies
from numpy import std, mean

_time_init = time.time()


"""
generate telomere counts from bam files
by 
A) aligning TTAGGGx6 or CCCTAAx6
B) counting TTAGGG (or CCCTAA the reverse-complement)
and correcting that for total read-depth and on-target rate
"""

comp = string.maketrans("ACGTNacgtn-", "TGCANtgcan-")
def revcomp(s):
    return s[::-1].translate(comp)


class MATCHBED:
    # given a bed file look for match 
    # with a target region
    def __init__ ( self, bedFile, bam_header) :
        self.bed = bedFile
        self.bin = 10**5
        self.db = dict()
        self.bam_header = bam_header
        self.is_target = self.bed_is_target
        try:
            # if bedfile is not a file
            self.parseBed()
            self.validate_contigs()
        except IOError:
            self.is_target = self.wg_is_target


    def validate_contigs( self ):
        header = set(self.bam_header).discard("*")
        if len(set(self.db.keys()) - set(self.bam_header)) == 0:
            reporting.reportInfo("contigs check passed")
            return True

        # if not; not good
        print "ERROR: bam file headers and bed file contigs do not match"
        print "CONTIGS in bedfile", set(self.db.keys())
        print "CONTIGS in bamfile", set(self.bam_header)
        print "unmatched in bamfile", set(self.db.keys()) - set(self.bam_header)
        print "ERROR: Unable to continue, fix the file headers."
        sys.exit()
    

    def parseBed ( self ) :
        with open(self.bed) as inf:
            for line in inf:
                line = line.strip().split("\t")

                try:
                    C, S, E = line[:3]
                    S = int(S)
                    E = int(E)

                    mS = S / self.bin
                    mE = E / self.bin

                    try:
                        self.db[C].setdefault(mS, list()).append((S, E))
                    except KeyError:
                        self.db[C] = {mS : [(S, E)]}

                    if mS != mE :
                        self.db[C].setdefault(mE, list()).append((S, E))
                except ValueError:
                    print "WARNING: improper bed line", line       


    def wg_is_target( self, C, S, E ):
        return False

        
    def bed_is_target ( self, C, S, E ):
        if C in self.db:
            mS = S / self.bin           
            for start, end in self.db[C].get(mS, list()):
                if start > E:
                    return False
                elif start <= S < end or start < E <= end or S <= start <= end <= E:
                    return True

            mE = E / self.bin
            if mS != mE :
                for start, end in self.db[C].get(mE, list()):
                    if start > E:
                        return False
                    elif start <= S < end or start < E <= end or S <= start <= end <= E:
                        return True        
        return False


class REPORTING:
    def __init__( self, silence ):
        self.silence = silence        
        self.start_main = time.time()
        
        try:
            rows, columns = os.popen('stty size', 'r').read().split()
        except ValueError:
            columns = "80"

        self.column = int(columns)


    def gethms( self, S):
        days, remain = divmod( S , (3600*24) )
        hours, remain = divmod( remain , 3600 )
        minutes, seconds = divmod( remain , 60 )
        return ':'.join(['%d%s'%(p,q) for p,q in zip((days, hours, minutes, seconds), ('d', 'h', 'm', 's')) if p > 0 ])


    def report_progress( self ):
        R = reads_done.value
        E = (time.time() - self.start_main)
        self.synopsis('Processed %s reads in %s [%.0f reads/sec]'%(self.splitthousands(str(R)), self.gethms(E), (R/E)))
        time.sleep(1)
        return 
           

    def synopsis( self, text ):
        if not self.silence:
            return
        
        sys.stderr.write('\r%s' % ''.join([text, ' '*(self.column-len(text))]))
        sys.stderr.flush()
        return


    def splitthousands( self, s, sep=','):  
        if len(s) <= 3:
            return s  
        return self.splitthousands(s[:-3], sep) + sep + s[-3:]


    def reportInfo( self, s):
        print >> sys.stdout, "INFO:", s


    def reportError(self, s):
        print >> sys.stdout, "ERROR:", s


class CALCTEL:
    def __init__(self, opts):
        self.MT = opts.mintel
        self.CNT = opts.count / 6.
        self.PHRED = opts.phred
        self.QUAL = opts.qual
        
        if self.MT.isdigit():
            self.MT = int(self.MT)
            self.mintel = self.mintel_abs
        else:
            self.MT = float(self.MT) / 6.
            self.mintel = self.mintel_frac

        if not self.QUAL:
            self.get_read_qual = bool
        else:
            self.get_read_qual = self.getReadQual

        if opts.dups:
            self.isDup = self.isDups
            self.corr_dups = self.corrDups
        else:
            self.isDup = self.noDups
            self.corr_dups = self.corrNoDups

        if opts.CC == '1':
            def correctionC(F_OFFT):
                # use an exponential decay correction
                return 10**(-1*F_OFFT**2)

        elif opts.CC == '2':
            def correctionC(F_OFFT):
                # use a exponential correction
                return (1 - F_OFFT)**-1

        elif opts.CC == '3':
            def correctionC(F_OFFT):
                # use an hyperbolic cosine correction
                return cosh(F_OFFT*5)

        elif opts.CC == '5':
            def correctionC(F_OFFT):
                # use a fifth order polynomial
                return 1-(10*F_OFFT**3-15*F_OFFT**4+6*F_OFFT**5)

        elif opts.CC == '6':
            def correctionC(F_OFFT):
                # use a third order polynomial
                return 1-(3*F_OFFT**2-2*F_OFFT**3)

        elif opts.CC == '7':
            def correctionC(F_OFFT):
                # exponential decay correction
                return (1.0-F_OFFT)**3

        elif opts.CC == '8':
            def correctionC(F_OFFT):
                # exponential decay correction
                return exp(-F_OFFT*3)

        else:
            # this is the default and CC=4
            def correctionC(F_OFFT):
                return 10**(F_OFFT**2)**4

        self.telcount = self.telcountA
        self.correctionC = correctionC


        if opts.type == 'adv':
            # change any two bases at a time in the first three positions
            self.Freps = set(["".join(p) for p in product(('T',), 'ACGT', 'ACGT',('GGG',))] + ["".join(p) for p in product('ACGT',('T',), 'ACGT', ('GGG',))] + ["".join(p) for p in product('ACGT', 'ACGT', ('TGGG',))])
            self.Rreps = [revcomp(p) for p in self.Freps]
            
        elif opts.type == 'complex':
            # change any one base at a time in the first three positions
            self.Freps = set(["".join(p) for p in product('ACGT', ('TAGGG',))] + ["".join(p) for p in product(('T',), 'ACGT', ('AGGG',))] + ["".join(p) for p in product(('TT',),'ACGT', ('GGG',))])
            self.Rreps = [revcomp(p) for p in self.Freps]
            
        elif opts.type == 'simple':
            # only second position variablity
            self.Freps = ['TCAGGG', 'TGAGGG', 'TTAGGG', 'TAAGGG']
            self.Rreps = [revcomp(p) for p in self.Freps]

        elif opts.type == 'basic':
            # only canonical variant sequence
            self.telcount = self.telcountB
            self.Freps = 'TTAGGG'
            self.Rreps = 'CCCTAA'


    def isDups( self, s ):
        return 'd' in s


    def noDups( self, s ):
        # nothing is a duplicate
        return 0


    def corrDups( self, Vt, Vl, Vm, D ):
        # correct for percent duplication
        P_DUPS = 1. - D
        return Vt*P_DUPS, Vl*P_DUPS, Vm*P_DUPS        


    def corrNoDups( self, Vt, Vl, Vm, D ):
        # correct for percent duplication
        return Vt, Vl, Vm       


    def getReadQual( self, qual ):
        # given a list of qual values
        # determine whether its pass or fail (30*6)
        return bool(sum([ord(p)-self.PHRED for p in qual])/(len(qual)*self.QUAL))


    def lookupRG( self, read, RG='null' ):
        for p in read[10:]:
            if p.startswith('RG:Z:'):
                return p.split('RG:Z:')[1]
        return RG


    def mintel_frac(self, read, repeat):
        # calculate relative fraction repeat based and Read-Lenght
        if read.count(repeat) >= len(read) * self.MT:
            return True
        return False
        

    def mintel_abs(self, read, repeat):
        # compare absolute repeat counts in read
        if read.count(repeat) >= self.MT:
            return True
        return False
        

    def count_tel(self, cnt, read):
        if cnt >= len(read) * self.CNT:
            return True
        return False


    def telcountA(self, Seq, Qual ):
        # count anything but basic
        # check if repeat 
        # minimum exact canonical-repeats
        # then, use a read-length correction
        for reps, R in zip((self.Freps, self.Rreps), ('TTAGGG', 'CCCTAA')):
            if self.mintel(Seq, R):
                telrepeats = dict((p, Seq.count(p)) for p in reps)
                cnt = sum(telrepeats.values())
                if self.count_tel(cnt, Seq):
                    return 1., cnt, [p for p,q in telrepeats.items() for i in xrange(q)]
        return 0., 0., list()


    def telcountB(self, Seq, Qual ):
        # count only basic
        # check if repeat 
        # use a read-length correction
        for R in (self.Freps, self.Rreps):
            cnt = Seq.count(R)
            if self.count_tel(cnt, Seq):
                return 1., cnt, [R] * cnt
        return 0., 0., list()


    def CORR(self, TELR, EXPF, ALLT, ONT, F_OFFT, CHRM, DUPS):
        CC = self.correctionC(F_OFFT)

        Vt = 10**6 * (TELR / (ALLT-ONT)) / CC
        # each repeat is 6-bases long, so you could multiply by 6
        Vl = 10**6 * (EXPF / (ALLT-ONT)) / CC
        try:
            Vm = 10**6 * (TELR / CHRM / ALLT)
        except ZeroDivisionError:
            Vm = 0

        return self.corr_dups(Vt, Vl, Vm, DUPS)

        
# get working
parser = OptionParser(usage="%prog -f [-b -t]", version=__version__ )
groups = OptionGroup(parser, "User input options")
groups.add_option("-f", dest="file",
                 help="path/to/infile.bam or comma-separated list of files to be merged and counted (only single-threaded supported for merged) [REQUIRED]")

groups.add_option("-b", dest="bedfile",
                 help="/path/to/targets.bed or whole genome 'wg' [%default]", default='wg')

groups.add_option('-t', dest='threads', type='int',
                 help="use multiprocessing based on read groups and each chr use 1 for non-mp, 0 for debug mp everthing else is multi-threaded [%default]",
                 default=8)

groups.add_option("-o", dest="outdir",
                 help="path/to/outdir/ [%default]",
                 default=os.getcwd())
parser.add_option_group(groups)

groups = OptionGroup(parser, "Advanced calibration options 'NOT recommended to alter unless you are recalibrating'")
groups.add_option("--count", dest="count",
                 metavar=' ', type='float',
                 help="number of repeats to count as a fraction of read-length [%default]",
                 default=0.8)

groups.add_option("--mintel", dest="mintel",
                 metavar=' ', type='str',
                 help="minimum fraction (ie:0.3) or count (ie:2) of canonical TTAGGG repeats in read [%default]",
                 default='2')

groups.add_option('--CC', dest='CC',
                 metavar=' ',
                 help="correction coefficient formula to use (1, 2, 3, 4, 5, 6, 7, 8) [%default]",
                 default='6')

groups.add_option('--qual', dest='qual',
                 metavar=' ', type='int',
                 help="the avg base quality score cutoff for good reads [%default]",
                 default=20)

groups.add_option('--type', dest='type', choices=('basic', 'simple', 'complex', 'adv'),
                 metavar=' ', type='choice',
                 help='counting approach; basic=TTAGGG, simple=TNAGGG, complex="TTNGGG,TNAGGG,NTAGGG" derivates total of 10", adv="TNNGGG,NTNGGG,NNAGGG" derivates total of 37" [%default]',
                 default='complex')

groups.add_option('--phred', dest='phred',
                 metavar=' ', type='int',
                 help="the quality score offset [%default]",
                 default=33)

groups.add_option('--dups', dest='dups',
                 metavar=' ', action='store_true',
                 help="correct for PCR/optical duplicates [%default]",
                 default=True)
parser.add_option_group(groups)

groups = OptionGroup(parser, "Runtime options")
groups.add_option('--verbose', dest='verbose',
                 metavar=' ', action='store_true',
                 help="print counts [%default]",
                 default=False)

groups.add_option('--samtools', dest='samtools',
                 metavar=' ',
                 help="/path/to/samtools [%default]",
                 default='samtools')

groups.add_option('--prefix', dest='prefix',
                 metavar=' ',
                 help="prefix for output name [%default]",
                 default='')
parser.add_option_group(groups)
(option, arg)= parser.parse_args()

reporting = REPORTING(option.verbose)

reporting.reportInfo("python PID %s"%os.getpid())

need_SP_merge = False
if ',' in option.file:
    option.threads = 1
    need_SP_merge = True
    
elif not option.file.startswith('ftp') and not os.path.isfile(option.file):
    reporting.reportError("need to provide an input.bam file")
    reporting.reportError("%s is not a valid file"%option.file)
    sys.exit(1)

if option.bedfile == None or option.bedfile.upper() == 'WG':
    option.bedfile = 'wg'
    reporting.reportInfo("running whole genome calculation")
elif not os.path.isfile(option.bedfile):
    reporting.reportError("bed file %s not located" %option.bedfile)
    sys.exit(1)

if not option.file.startswith('ftp'):
    # check .bam.bai file
    # rename .bai file
    if not os.path.isfile(option.file+'.bai'):
        if os.path.isfile(option.file.replace('.bam', '.bai')):
            os.rename(option.file.replace('.bam', '.bai'),
                      option.file+'.bai')
        else:
            subP.call([option.samtools, 'index', option.file],
                      stdout=subP.PIPE, bufsize=0) 

# collect contig names from bam file
Chr_list = sorted([p.strip().split('\t')[0] for p in subP.Popen([option.samtools, 'idxstats', option.file.split(',')[0]], stdout=subP.PIPE, bufsize=0).stdout])
reporting.reportInfo("Bamfile contig count is %s"%len(Chr_list))


# These are the types of information that we are collecting
keepCountsDB = {'allR' : 1.0, 'ontR' : 0.0, 'telR' : 0.0,
                'chrM' : 0.0, 'telU' : 0.0, 'telC': 0.0,
                'dups' : 0.0, 'freq' : dict()}

# 1 'allR' : cnt of all reads in the read group
# 2 'ontR' : cnt of on target reads in the read group
# 3 'telR' : cnt of telomeric reads in the read group
# 4 'chrM' : cnt of chrM reads in the read group
# 5 'telU' : cnt of telomeric repeats from the telR unmapped reads in the read group
# 6 'telC' : cnt of telomeric repeat from the telR reads in the read group
# 'dups' : cnt of duplicate reads in the read group
# 'freq' : cnt of individual telomeric hexamers in the read group

def counter_SP( merge ):
    global tabixtools

    if merge:
        flist = option.file.split(',')
        # read samtools view in tandem for all the files
        reporting.reportInfo("reading from %s files in tandem" % len(flist))
        stmt = list()
        for each in flist:
            stmt.extend([option.samtools, 'view', '-X', each, ';'])

        bamF = subP.Popen(" ".join(stmt).rstrip(";"), shell=True, stdout=subP.PIPE, bufsize=0)
        
    else:
        reporting.synopsis("starting %s" %sname[:40])
        bamF = subP.Popen([option.samtools, 'view', '-X', option.file],
                           stdout=subP.PIPE, bufsize=0)

    RGdata = dict()
    for cnt, read in enumerate(bamF.stdout, 1):
        read = read.strip().split("\t")
        Qual = read[10]

        if calctel.get_read_qual(Qual):
            Seq = read[9].upper()
            Chr = read[2]
            RG = calctel.lookupRG(read)

            # add to total read-group reads
            try:
                RGdata[RG]['allR'] += 1.
            except KeyError:
                RGdata[RG] = deepcopy(keepCountsDB)

            RGdata[RG]['dups'] += calctel.isDup(read[1])

            if Chr == '*':
                # unmapped
                tc, freq, telrepeats = calctel.telcount(Seq, Qual)
                RGdata[RG]['telR'] += tc
                RGdata[RG]['telU'] += tc
                RGdata[RG]['telC'] += freq
                for trep in telrepeats:
                    RGdata[RG]['freq'][trep] = RGdata[RG]['freq'].get(trep, 0) + 1

            else:
                # mapped
                Spos = int(read[3])
                Epos = Spos + len(Seq)

                if Chr in ('chrM', 'chrMT', 'MT', 'M'):
                    RGdata[RG]['chrM'] += 1.
                elif tabixtools.is_target(Chr, Spos, Epos):
                    # add to total on-target-reads
                    RGdata[RG]['ontR'] += 1.
                else:
                    tc, freq, telrepeats = calctel.telcount(Seq, Qual)
                    RGdata[RG]['telR'] += tc
                    RGdata[RG]['telC'] += freq
                    for trep in telrepeats:
                        RGdata[RG]['freq'][trep] = RGdata[RG]['freq'].get(trep, 0) + 1

    reporting.synopsis('Processed %s reads in %s'%(reporting.splitthousands(str(cnt)), reporting.gethms(time.time()-_time_init)))
    reads_done.value = cnt
    bamF.terminate()
    return RGdata


def counter_MP( Chr ):
    global tabixtools

    # collect read-count, on-target-reads, telomeric reads, chrM reads, telomeric reads in unmapped, tel_repeat_cnt
    reporting.synopsis("starting %s with %s" %(sname[:25], Chr))
    cnt = 0
    RGdata = dict()

    if Chr in ('chrM', 'chrMT', 'MT', 'M'):
        # chrM reads
        bamF = subP.Popen([option.samtools, 'view', '-X', 
                           option.file, Chr],
                           stdout=subP.PIPE, bufsize=0)

        for cnt, read in enumerate(bamF.stdout, 1):
            read = read.strip().split("\t")
            Qual = read[10]          
            if calctel.get_read_qual(Qual):
                RG = calctel.lookupRG(read)
                Spos = int(read[3])
                Epos = Spos + len(Qual)

                # add to total read-group reads
                try:
                    RGdata[RG]['allR'] += 1.
                except KeyError:
                    RGdata[RG] = deepcopy(keepCountsDB)

                RGdata[RG]['dups'] += calctel.isDup(read[1])

                if tabixtools.is_target(Chr, Spos, Epos):
                    # add to total on-target-reads
                    RGdata[RG]['ontR'] += 1.
                RGdata[RG]['chrM'] += 1.

    elif Chr == '*':
        # unmapped
        bamF = subP.Popen([option.samtools, 'view', '-X', 
                           '-f4', option.file],
                           stdout=subP.PIPE, bufsize=0)

        for cnt, read in enumerate(bamF.stdout, 1):
            read = read.strip().split("\t")
            Qual = read[10]          
            if calctel.get_read_qual(Qual):
                Seq = read[9].upper()
                RG = calctel.lookupRG(read)

                # add to total read-group reads
                try:
                    RGdata[RG]['allR'] += 1.
                except KeyError:
                    RGdata[RG] = deepcopy(keepCountsDB)

                RGdata[RG]['dups'] += calctel.isDup(read[1])

                # check if telomeric repeat
                tc, freq, telrepeats = calctel.telcount(Seq, Qual)
                RGdata[RG]['telR'] += tc
                RGdata[RG]['telU'] += tc
                RGdata[RG]['telC'] += freq
                for trep in telrepeats:
                    RGdata[RG]['freq'][trep] = RGdata[RG]['freq'].get(trep, 0) + 1

    else:
        # mapped
        bamF = subP.Popen([option.samtools, 'view', '-X', option.file, Chr],
                           stdout=subP.PIPE, bufsize=0)
        
        for cnt, read in enumerate(bamF.stdout, 1):
            read = read.strip().split("\t")
            Qual = read[10]
            if calctel.get_read_qual(Qual):
                Seq = read[9].upper()
                Chr = read[2]
                Spos = int(read[3])
                Epos = Spos + len(Seq)
                RG = calctel.lookupRG(read)

                # add to total read-group reads
                try:
                    RGdata[RG]['allR'] += 1.
                except KeyError:
                    RGdata[RG] = deepcopy(keepCountsDB)
    
                RGdata[RG]['dups'] += calctel.isDup(read[1])

                if tabixtools.is_target(Chr, Spos, Epos):
                    # add to total on-target-reads
                    RGdata[RG]['ontR'] += 1.
                else:
                    tc, freq, telrepeats = calctel.telcount(Seq, Qual)
                    RGdata[RG]['telR'] += tc
                    RGdata[RG]['telC'] += freq
                    for trep in telrepeats:
                        RGdata[RG]['freq'][trep] = RGdata[RG]['freq'].get(trep, 0) + 1

    bamF.terminate()
    # add to total read-group reads
    reads_done.value += cnt
    reporting.report_progress()

    return RGdata


def merge_RGdata(s):
    global header

    # first check for data type, should be a dict
    if type(s) is not dict:
        # collect all the outputs from the multiprocesses
        # they are tuples of dicts of dicts
        TelRepeats = dict()
        RGdata = dict()
        for each_dict in s:
            # there are multiple RG tags with dict of values
            for k,v in each_dict.iteritems():
                TR = v.pop('freq') # the telrepeat dict
                try:
                    for i in v.keys():
                        RGdata[k][i] += v[i]
                    for l,m in TR.items():
                        TelRepeats[k][l] = TelRepeats[k].get(l, 0) + m
                except KeyError:
                    RGdata[k] = v
                    TelRepeats[k] = TR
    else:
        RGdata = s
        TelRepeats = dict((p, RGdata[p].pop('freq')) for p in RGdata)
        
    # for each read-group in the bam file generate the fraction

    outlist = list()
    for_ont = list() # calc mean based on off-target tel read cnt
    for_chrm = list() # calc mean based on chrM ratio
    for_ontC = list() # calc mean based on off-target telomere repeat cnt
    TelRepeats_all = dict()
    dup_rates = list() # calc mean duplication rate
    ont_rates = list() # calc mean on-target rate

    for RG in RGdata:
        ALLT = RGdata[RG]['allR']
        ONT = RGdata[RG]['ontR']
        TELR = RGdata[RG]['telR']
        CHRM = RGdata[RG]['chrM']
        TELUN = RGdata[RG]['telU']
        EXPF = RGdata[RG]['telC']
        DUPS = RGdata[RG]['dups'] / ALLT
        
        # ASSUME: the better the capture efficiency
        # the worse the telomeric fraction
        # i.e. if your on-target is 0 you should assume that
        # the telomeric count represents WGS so take it as is.
        F_OFFT = ONT / ALLT

        Vt, Vl, Vm = calctel.CORR(TELR, EXPF, ALLT, ONT, F_OFFT, CHRM, DUPS)

        if option.verbose:
            print RG, RGdata[RG]   
            print Vt, Vl, Vm
        
        for_ont.append(Vt)
        for_ontC.append(Vl)
        for_chrm.append(Vm)
        dup_rates.append(DUPS)
        ont_rates.append(F_OFFT)

        TelRepeats_RG = dict()
        # combine reverse-complement sequences
        for T in TelRepeats[RG].keys():
            tcnt = TelRepeats[RG][T]
            # look for rev-comp sequences
            if T.startswith('CCC'):
                T = revcomp(T)

            TelRepeats_all[T] = TelRepeats_all.get(T, 0) + tcnt
            TelRepeats_RG[T] = TelRepeats_RG.get(T, 0) + tcnt
            
        TRcnt = float(sum(TelRepeats_RG.values()))
        TR = sorted(TelRepeats_RG.items(), key=lambda x: x[1], reverse=True)
        TRstr = "\t".join(["%s:%.4f" %(p, q/TRcnt) for p,q in TR])

        outlist.append('\t'.join([RG, '%i'%ALLT, '%i'%ONT, '%i'%TELR, '%i'%CHRM,  '%i'%TELUN, '%i'%EXPF, str(Vt), str(Vm), str(Vl), '%.2f'%DUPS, TRstr ]) + '\n')
   
    TRcnt = float(sum(TelRepeats_all.values()))
    TR = sorted(TelRepeats_all.items(), key=lambda x: x[1], reverse=True)
    TRstr = "\t".join(["%s:%.4f" %(p, q/TRcnt) for p,q in TR if q/TRcnt > 0.001])
        
    outlist.append('\n' + '\t'.join([sname, 'MEAN+-SD',
                                     'RTPMO-w/offt', 
                                     str(mean(for_ont)),
                                     str(std(for_ont)), 
                                     'RTPKM-w/chrM', 
                                     str(mean(for_chrm)),
                                     str(std(for_chrm)),
                                     'RTPMO-w/offtCNT',
                                     str(mean(for_ontC)),
                                     str(std(for_ontC)),
                                     'ONTARGET-RATE',
                                     str(mean(ont_rates)),
                                     str(std(ont_rates)),
                                     'DUP-RATE',
                                     str(mean(dup_rates)),
                                     str(std(dup_rates)),
                                     TRstr]) + '\n')
    return outlist


option.outdir = os.path.realpath(option.outdir)
if not os.path.isdir(option.outdir):
    os.makedirs(option.outdir)
    
# start processing
calctel = CALCTEL(option)
tabixtools = MATCHBED(option.bedfile, Chr_list)
sname = os.path.basename(option.file).split(',')[0].split('.bam')[0]

reads_done = Value('i', 0)
if option.threads == 1:
    reporting.reportInfo("Running single-threaded")
    D = counter_SP(need_SP_merge)
else:
    print "INFO: Running multi-threaded",
    # start with unmapped, since it may take the most
    # then, add the rest of the contigs
    if option.threads == 0:
        print "debug-mode"
        # debug mp mode
        D = [counter_MP(p) for p in Chr_list]       
    else:
        print "with %s" %option.threads
        jobP = Pool(option.threads)
        jobR = jobP.map_async(counter_MP, Chr_list)
        jobP.close()
        jobP.join()
        D = jobR.get()

output_data = merge_RGdata(D)
outname = "" if not len(option.prefix) else "%s_"%option.prefix
outname += '_'.join([sname, "telcount_by", str(option.count), option.type, str(option.qual), option.CC, str(option.mintel)]) 
outname += "" if not option.dups else "_dups"
outname += '.tsv'

with open(os.path.join(option.outdir, outname), "w") as outf:
    outf.write('\t'.join(['SAMPLE', sname]) + '\n')
    outf.write('\t'.join(['TARGET', option.bedfile]) + '\n')
    cmd = " ".join(["%s=%s"%p for p in sorted(option.__dict__.items())])
    outf.write('\t'.join(['CMD', sys.argv[0] + " %s"%cmd]) + '\n')
    outf.write('\t'.join(['RG', 'ALL', 'ONTARGET', 'TELOMERE', 
                          'CHRM', 'TEL-UNMAPPED', 'CNT_TEL_LEN', 
                          'RTPMO-w/offt',
                          'RTPKM-w/chrM',
                          'RTPMO-w/offtCNT',
                          'FRACTION_DUPS',
                          'TELOMERE_SEQ_FOUND']) + '\n')   
    outf.writelines(output_data)

CNT = reads_done.value
TIME = time.time()-_time_init

reporting.reportInfo("Process completed: %s reads in %s [%s r/s]" %(reporting.splitthousands(str(CNT)), reporting.gethms(TIME), reporting.splitthousands(str(int(float(CNT)/TIME)))))
