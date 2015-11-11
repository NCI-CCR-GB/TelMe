#!/usr/local/env python

"""
This is a cluster job submitter
can do local or swarm
"""

import sys
import os
import subprocess as SP
from itertools import product
import math
import random
import glob
import platform
from optparse import OptionParser

def updateOpts(option, opt_str, value, parser, *args, **kwargs ):
    """
    update the user-input path variables using callback
    """
    if option.dest == 'outpath':
        value = os.path.realpath(os.path.expanduser(value))
        if not os.path.isdir(value):
            os.mkdir(value)
        os.chdir(value)

    elif option.dest == 'infiles':
        value = os.path.realpath(os.path.expanduser(value))

    elif not option.dest.islower():
        pass
        
    setattr(parser.values, option.dest, value)


def submitcluster(opts):
    SP.call("~/bin/qsubstuff.sh -t %s python %s %s"%(options.threads, options.script, opts),
        shell=True)
    
def submitbiowulf(opts):
    return "python %s %s\n"%(options.script, opts)


parser=OptionParser()

parser.add_option('-q', dest='queue', default='all',
           help='designate a queue for biowulf or run all [%default]')

parser.add_option('-b', dest='bed',
           help='the bedfile for bait regions "PRESET"')

parser.add_option('-t', dest='threads', default=8, type='int',
           help='number of threads to use for each run [%default]')

parser.add_option('-d', dest='data', default='1000g', type='choice',
           choices=('1000g', 'nci60'),
           help='data file type to process [%default]')

parser.add_option('-s', dest='script', default='~/bin/telME/telME.py',
           help='/path/to/telomere_counter.py [%default]')

parser.add_option('-i', dest='infiles', default=os.path.dirname(os.getcwd()),
           help='/path/to/input_dir for bam files [%default]',
           action="callback", callback=updateOpts, type="str")

parser.add_option('-o', dest='outpath', default=os.getcwd(),
           help='/path/to/output_dir [%default]',
           action="callback", callback=updateOpts, type="str")

parser.add_option('-m', dest='multi', default=2, type='int',
           help='multiply the ccr queue odds [%default]')

parser.add_option('--SL', dest='sampleList', default='all',
           help='comma-sep list of sample names used for startswith [%default]')

parser.add_option('-C', dest='COUNT', default='0.8',
                  help='%default')

parser.add_option('--CC', dest='CC', default='4,5,6,8',
                  help='%default')

parser.add_option('-Q', dest='QUAL', default='20',
                  help='%default')

parser.add_option('-T', dest='TYPE', default='basic,complex',
                  help='%default')

parser.add_option('-D', dest='DUPS', default="--dups",
                  help='%default')

parser.add_option('-M', dest='MINTEL', default='2',
                  help='%default')
           
(options, args) = parser.parse_args()

if platform.node().split(".")[0] == "biowulf":
    cluster = 'BIOWULF'
else:
    cluster = 'linux2-local'

if options.data == '1000g':
    bamfiles  = ('*exome*.bam', '*low_c*.bam')
    if options.bed == None:
        bedfiles = ('~/beds/1000g_SureSelect_SeqCap_Exomes_merged_sorted_b37lo.bed', 'wg')
    else:
        bedfiles = (options.bed, 'wg')
else:
    bamfiles = ('*.bam',)
    bedfiles = ('~/beds/0024A_Exome_InSolution_hg19.bed',)

allVars = ('COUNT', 'CC', 'QUAL', 'TYPE', 'MINTEL', 'DUPS')
allVarValues = [tuple(getattr(options, opt).split(',')) for opt in allVars]
allJobs = [p for p in product(*allVarValues)]
allBams = len(glob.glob(os.path.join(options.infiles, "*.bam")))
print 'Running:',allVarValues

if allBams:
    len_jobs = len(allJobs)*len(bamfiles)*allBams
else:
    print "Unable to locate bam files in", options.infiles
    sys.exit(1)
print '# of jobs = %s' % (len_jobs)

### for 1000g files
runOpts = "--count=%s --CC=%s --qual=%s --type=%s --mintel=%s %s"
if cluster.upper().startswith('BIO'):
    nameCode = random.randint(100, 999)
    """
    bundle = max(1, int(math.ceil(len_jobs / 8000.))) # 4000 jobs but 2-4 jobs per nore
    print 'using -b', bundle
    userInp=raw_input('is this OK [Y] or #: ')
    if not userInp.upper().startswith('Y'):
        try:
            bundle = int(userInp)
        except ValueError:
            sys.exit()
            """
    if options.queue == 'ccr':
        Qprob = ['-q ccr']
    elif options.queue == 'norm':
        Qprob = ['']
    else:
        Qprob = list()
        # there are 14648 norm cores
        # and 4096 ccr cores
        # scale by current batch limits
        for line in SP.Popen('batchlim', stdout=SP.PIPE).stdout:
            line = line.strip().split()
            if len(line) and line[0] == 'norm':
                ODDS=round(14648.0*float(line[1])/10**6)
                print "norm queue odds:", ODDS
                Qprob.extend(['' for i in range(int(ODDS))])

            elif len(line) and line[0] == 'ccr':
                # ccr runs faster so double the odds
                ODDS=round(options.multi*4096.0*float(line[1])/10**6)
                print 'ccr queue odds:', ODDS
                Qprob.extend(['-q ccr' for i in range(int(ODDS))])
                break

        random.shuffle(Qprob)

    bundle = max(1, int(math.ceil(len_jobs / 8000.))) # 4000 jobs but 2-4 jobs per nore
    print 'using -b', bundle
    userInp=raw_input('is this OK [Y] or #: ')
    if not userInp.upper().startswith('Y'):
        try:
            bundle = int(userInp)
        except ValueError:
            sys.exit()

    fileOpts = "-b %s -f /scratch/%s" + ' -t %s ' %(options.threads)
    # BIOWULF
    for bedf, bamf in zip(bedfiles, bamfiles):
        for bamfile in glob.glob(os.path.join(options.infiles, bamf)):
            filename = os.path.basename(bamfile)
            stmt = fileOpts%(bedf, filename)
            fname = filename.split('.bam')[0]
            swarmName = '_swarm%s_%s' %(nameCode, fname)
            with open(swarmName, 'w') as outf:
                for opts in allJobs:
                    outf.write(submitbiowulf(stmt + runOpts%opts))
            
            # now swarm the job
            SP.call('swarm -b %s -t %s -R gpfs -f %s %s --prologue "clearscratch && cp %s /scratch/"'%(bundle, options.threads, swarmName, random.choice(Qprob), os.path.join(options.infiles,fname+'.b*')), shell=True)
        
else:
    # LOCAL
    userInp=raw_input('is this OK [Y/N]: ')
    if not userInp.upper().startswith('Y'):
        sys.exit()

    stmt = '-t %s '%(options.threads-1) + '-b %s -f %s ' 
    for bedf, bamf in zip(bedfiles, bamfiles):
        for opts in allJobs:
            for bamfile in glob.glob('../%s'%bamf):
                submitcluster(stmt%(bedf, bamfile)+runOpts%opts)


# touch the bam folder to make sure that the timestamp is updated
SP.call('touch %s'%options.infiles, shell=True)
