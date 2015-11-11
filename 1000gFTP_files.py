
import glob
import ftplib
from ftplib import FTP
import os
import sys
import random
from itertools import product
import multiprocessing as mp
import subprocess as subP

try:
    local=sys.argv[1]
except IndexError:
    local=False

def picard(args):
    for each in args:
        if '.mapped.' in each:
            MAPPED = os.path.basename(each)
        else:
            UNMAPPED = os.path.basename(each)

    MERGED = os.path.basename(MAPPED).replace('.mapped.', '.merged.', 1)

    # picard merge the files
    stmt = "java -jar -Xmx4g ~/bin/PICARD/picard-tools-1.110/MergeSamFiles.jar INPUT=%s INPUT=%s OUTPUT=%s ASSUME_SORTED=true USE_THREADING=true VERBOSITY=ERROR VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true COMPRESSION_LEVEL=9 && " %(MAPPED, UNMAPPED, MERGED)
    # remove the downloaded files
    stmt += "rm -f %s %s && " %(MAPPED, UNMAPPED)
    # rename the index file
    stmt += "mv %s %s \n" %(MERGED.replace('.bam','.bai',1), MERGED.replace('.bam','.bam.bai',1))

    return stmt


def samtoolsMerge(args):
    for each in args:
        if '.mapped.' in each:
            MAPPED = os.path.basename(each)
        else:
            UNMAPPED = os.path.basename(each)

    MERGED = MAPPED.replace('.mapped.', '.merged.', 1)

    stmt = "samtools cat %s %s | samtools sort -l 9 -@ 4 -m 2G - %s && samtools index %s" % (args[0], args[1], MERGED.split('.bam')[0], MERGED)
    

# collect all the list of files from FTP site
ftp = FTP('ftp-trace.ncbi.nih.gov')
ftp.login()
ftp.cwd('/1000genomes/ftp/data/')
snames = ftp.nlst()
random.shuffle(snames)
snames = random.sample(snames, 40)

"""
# Broad samples
snames = ['HG02307','HG02309','NA12399','HG01464','NA19921',
          'NA19346','HG00551','NA19474','NA19247','HG01950',
          'HG01492','NA19725','NA19092','NA19713','NA20893']
# BGI samples
snames = ["HG00657","HG01173","HG01848","HG01882",
          "HG02260","NA21143","HG01082","HG00525",
          "NA18593","HG00342","HG00422","HG01350"]
# WU samples
snames = ['NA18968','NA18969','NA18970','NA18972','NA18974',
          'NA18975','NA18976','NA18981','NA19098','NA19119',
          'NA19141','NA19143','NA19152','NA19153','NA19159',
          'NA19204','NA19206','NA19207','NA19209','NA19200',
          'NA19131']
# BCM samples
snames = ['NA18645', 'NA19096', 'NA19117', 'NA21105', 'NA19310',
          'NA21107', 'HG00355', 'HG01630', 'HG01624', 'HG00351',
          'NA18748']

# WU2 samples
snames = ['NA18972','NA19321','NA18924','NA20276','HG01874','NA19770',
          'NA20332','NA19777','NA19143','NA19144','HG00096','NA19746',
          'HG01497','NA19204','HG01494','NA19750','NA19160','NA18923',
          'NA20342','HG02470','NA20289','NA19783','NA19732','NA20346',
          'HG01498','HG02070','NA20127','HG01954','NA18969']

# BCM2 samples
snames = ['HG00235','NA19146','NA21108','NA20363','HG00349','HG00130',
          'HG02496','NA19740','HG01617','NA18527','NA19309','NA18645',
          'NA19310','NA19185','HG02285','HG02084','HG01271','NA19117',
          'HG02484','NA19149','HG01275','NA19214','NA19001','HG00129',
          'HG01618','NA19031','HG02002','NA21109','NA18643','NA21111',
          'HG02508','HG02138','NA18957','NA19028','NA12815','NA18642',
          'HG00355','HG02104','HG01625']
"""
ftppwd = "ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/data/"
files = {'G': list(), 'E': list()}

while len(snames):
    sname=snames.pop()
    try:
        G = ftp.nlst(sname + '/alignment/*mapped*.bam')
        E = ftp.nlst(sname + '/exome_alignment/*mapped*.bam')
        if len(G) == len(E) == 2 :
            files['G'].append([os.path.join(ftppwd, p) for p in G])
            files['E'].append([os.path.join(ftppwd, p) for p in E])           
    except ftplib.error_temp :
        pass
ftp.close()

if local:
    qsub = '~/bin/qsubstuff.sh -m 4 %s && rm %s'
else:
    qsub = 'swarm -t 4 -R gpfs -f %s --epilogue "rm %s"'

# download the files and fire up the merge process
joblist = list()
with open('download_file_list', 'w') as outf:
    for p in files:
        for q in files[p]:
            # easier on the FTP download speeds
            subP.call('wget %s %s'%tuple(q), shell=True)

            if os.path.isfile(os.path.basename(q[0])) and os.path.isfile(os.path.basename(q[1])):
                tname = os.path.basename(q[1]).replace('.bam', '.jobat', 1)
                with open(tname, 'w') as ff:
                    ff.write( picard(q) )

                # once the files are on disk, merge
                subP.call(qsub%(tname, tname), shell=True)

                joblist.append( picard(q) )
                outf.write("\n".join(q) + "\n")

"""
# download all the files at once
# run wget in the file
subP.call('wget -i download_file_list', shell=True)
           
if local:
    for each in joblist:
        subP.call('~/bin/qsubstuff.sh -m 4 %s'%each, shell=True)

else:
    with open("_swarm_1000g_files", 'w') as outf:
        outf.writelines(joblist)
        
"""
