#!/usr/bin/env python
__version__ = '0.0.1'
import os
import sys
import csv
import time
import glob
import random
#from scipy.stats import ttest_ind as ttest
from optparse import OptionParser
from Bio import Cluster
from copy import deepcopy

import numpy as np
import matplotlib
matplotlib.use("Agg")
# or use 'PDF'
import matplotlib.pyplot as plt
prop = matplotlib.font_manager.FontProperties(size=6) 

_time_init = time.time()

### GO INTO OPTIONS ###
parser= OptionParser(usage="%prog", version=__version__)

parser.add_option('-f', dest='indir',
                  help='/path/to/infile_directory [REQUIRED]')

parser.add_option('-d', dest='deco',
                  help='decorator string to locate files [%default]',
                  default='.tsv')

parser.add_option('--repeats', dest='repeats', type='int',
                  default=-8,
                  help='location of telomeric repeat sequences [%default]')

parser.add_option('-p', dest='plot', type='choice',
                  choices=('offt', 'chrm', 'offtcnt'),
                  default='offt',
                  help='plot value choices: offt,chrM,offtcnt  [%default]')

parser.add_option('-i', dest='ignore',
                  help='comma-separated list of sample names to ignore [%default]', default='')

parser.add_option('-n', dest='normals',
                  help='/path/to/summary_normarls file [%default]')

(option, arg)= parser.parse_args()

# collect and summarize
outprefix = os.path.join(option.indir, "_".join([p for p in ("summary", option.plot, option.deco.lstrip("_").replace(".tsv", "", 1).replace("*", "-")) if len(p)]))

infile = outprefix + '.txt'
ignore_snames = option.ignore.split(',')
flist = glob.glob(os.path.join(option.indir, "*%s"%option.deco))
if len(flist) == 0:
    print 'no files found, check deco %s' %option.deco
    sys.exit(1)

outdata = dict()
for fname in flist:
    with open(fname) as inf:
        line = inf.readline().strip().split("\t")
        sname = line[1]
        for line in inf:
            pass
        # write the last line
        outdata[sname] = line

prefix = os.path.commonprefix([p[::-1] for p in outdata.keys()])[::-1]
try:
    outdata = dict((p.split(prefix)[0],q) for p,q in outdata.items())
except ValueError:
    pass

scaleDB = {'offt' : 'RTPMO-w/offt',
           'chrm' : 'RTPKM-w/chrM',
           'offtcnt' : 'RTPMO-w/offtCNT',
           'dups' : 'DUP-RATE',
           'ont' :  'ONTARGET-RATE'}

main_values = scaleDB.values()
dataDB = dict()

with open(infile, 'w') as outf:
    # read the file in memory
    for sname, line in outdata.items():
        line = line.strip().split("\t")
        if sname not in ignore_snames:
            dataDB[sname] = dict()
            outf.write(outdata[sname])

            """
            # consume the repeats
            dataDB[sname]['repeats'] = dict()
            for i in xrange(1, len(line)):
                if ':' in line[-i]:
                    dataDB[sname]['repeats'].update(tuple(line[-i].split(":")))
                else:
                    break
                    """
            
            for i in xrange(2,len(main_values)*3,3):
                if line[i] in main_values:
                    dataDB[sname][line[i]] = [float(line[i+1]), float(line[i+2])]
                else:
                    break

            dataDB[sname]['repeats'] = dict(tuple(p.split(":")) for p in line[i+3:])

# if present; read the normals file as well
normalsDB = dict()
normal_repeatsDB = dict()
if option.normals != None and os.path.isfile(option.normals):
    with open(option.normals) as inf:
        csv_inf = csv.reader(inf, delimiter='\t')
        for line in csv_inf:
            sname = line[0].split(prefix)[0]

            normalsDB[sname] = dict()
            for i in xrange(2,len(main_values)*3,3):
                if line[i] in main_values:
                    normalsDB[sname][line[i]] = [float(line[i+1]), float(line[i+2])]
                else:
                    break

            normalsDB[sname]['repeats'] = dict(tuple(p.split(":")) for p in line[i:])
           
            for repeat in normalsDB[sname]['repeats']:
                normalsDB[sname]['repeats'][repeat] = float(normalsDB[sname]['repeats'][repeat])
                normal_repeatsDB.setdefault(repeat, list()).append(normalsDB[sname]['repeats'][repeat])

    # define normal mean and SD
    normal_data = [normalsDB[p][scaleDB[option.plot]][0] for p in normalsDB]
    normal_mean = np.mean(normal_data)
    normal_std = np.std(normal_data)

    normal_repeatsDB = dict((p,sum(normal_repeatsDB[p])/len(normalsDB)) for p in normal_repeatsDB)


# plot 1
# bar graph for sorted scale
sort_data = list()
repeatsDB = dict()
for sname in dataDB.keys():
    line = [sname] + dataDB[sname][scaleDB[option.plot]]
    sort_data.append(line)
    # also collect repeat info
    for repeat in dataDB[sname]['repeats']:
        dataDB[sname]['repeats'][repeat] = float(dataDB[sname]['repeats'][repeat])
        repeatsDB.setdefault(repeat, list()).append(dataDB[sname]['repeats'][repeat])
  
sort_data = sorted(sort_data, key=lambda x:x[1])
len_data = len(sort_data)
ind = np.arange(len_data)
Y_axis = list()
X_axis = list()
E_bars = list() # error bars
for p in sort_data:
    X_axis.append(p[0])
    Y_axis.append(p[1])
    E_bars.append(p[2])


PLOT = plt.bar(ind, Y_axis, yerr=E_bars, linewidth=0.1)
bw = PLOT.patches[0].get_width() / 2
plt.xlim(xmax=len_data)
ymin, ymax = plt.ylim()
plt.xticks(ind + bw, X_axis, rotation=90, fontsize=8, ha='center')

plot_title = "%s Values plot"%scaleDB[option.plot]

if len(normalsDB):
    # plot normals
    plt.axhspan(normal_mean - normal_std, normal_mean + normal_std, color='r', alpha=0.2)
    plt.axhline(normal_mean, color='r', alpha=0.2)

    """
    # run stats
    p_value = ttest([p for p in Y_axis if p<normal_mean], [p for p in Y_axis if p>normal_mean])[1]#, equal_var=False)
    plot_title += '\nLonger_vs_shorther than Normals p_value=%0.2e'%p_value
    """
    
plt.title(plot_title, fontsize=16, fontweight='bold')
plt.tight_layout()
plt.savefig(outprefix + '.pdf')
plt.close()

# plot 2
# distribution of all the hexamer counts
# collect the most prevelant repeats

bck_repeatsDB = deepcopy(repeatsDB)
with open(outprefix + '_cluster.txt', 'w') as outf:
    repeatsDB = dict((p,sum(repeatsDB[p])/len_data) for p in repeatsDB)
    sort_repeats = sorted(repeatsDB.items(), key=lambda x: x[1], reverse=True)
    limit_bottom = len_data * [0]
    if len(normalsDB):
        repeatsDB = normal_repeatsDB
        X_axis.append('NORMALS')
        limit_bottom.append(0)
        len_data += 1
        ind = np.arange(len_data)

    outf.write("name\t" + "\t".join(X_axis) + "\n")
    legendDB = dict()
    #all_colors = open(os.path.join(os.path.dirname(sys.argv[0]), 'colors_hex_values.txt')).read().split('\n')
    #random.shuffle(all_colors)
    all_colors = ("navy", "green", "maroon",
                  "blue", "lime", "red",
                  "teal", "olive", "purple",
                  "aqua",  "yellow", "fuchsia",
                  "FireBrick", "Indigo", "DodgerBlue", "DarkGreen")

    top_repeats = list()

    delta_repeats = dict()

    for repeat,color in zip(sort_repeats, all_colors[:15]):
        repeat = repeat[0]
        top_repeats.append(repeat)
        temp = [float(dataDB[p]['repeats'].get(repeat, '0.')) for p in X_axis if p in dataDB]
        delta_repeats[repeat] = temp
        if len(normalsDB):
            temp.append(np.mean(normal_repeatsDB.get(repeat, 0.)))
        outf.write('%s\t' %repeat + '\t'.join([str(p-repeatsDB[repeat]) for p in temp]) + '\n')    
        legendDB[repeat] = plt.bar(ind, temp, bottom=limit_bottom, color=color, ec=color, linewidth=0.1)
        limit_bottom = [p+q for p,q in zip(limit_bottom, temp)]
                
    # add-up to 1.0 with summing up all others       
    temp = [1.0-p for p in limit_bottom]
    legendDB['other'] = plt.bar(ind, temp, bottom=limit_bottom, color='gray', ec='gray', linewidth=0.1)
    top_repeats.append('other')
    #what are 'others'
    repeatsDB['other'] = sum([repeatsDB[p] for p in repeatsDB if p not in top_repeats])/len(dataDB)
    
    outf.write('%s\t' % 'other' + '\t'.join([str(p-repeatsDB["other"]) for p in temp]) + '\n')

try:
    bw = legendDB['TTAGGG'].patches[0].get_width() / 2
except KeyError:
    print >> sys.stderr, '%s ERROR: no TTAGGG for' %os.path.basename(sys.argv[0]), option.deco
    
# format the plot
plt.ylabel("fraction", fontsize=12, fontweight='bold')
plt.title("Distribution of Hexameric Repeats", fontsize=16, fontweight='bold')
plt.xticks(ind+bw, X_axis, rotation=90, size=8, ha='center')
plt.xlim(xmax=len_data)
ymin_scale = min([dataDB[p]['repeats'].get(top_repeats[0], 0) for p in X_axis if p in dataDB])*.9
plt.ylim(ymax=1, ymin=ymin_scale)

#top_repeats += ['other']
plt.figlegend([legendDB[p][0] for p in top_repeats], top_repeats,
              ncol=2, numpoints=1, prop=prop, loc=5)
plt.tight_layout()
plt.savefig(outprefix + '_hexrepeats.pdf', format="pdf", bbox_inches='tight')
plt.close()

"""
# plot 3
# group hexamers
# define groups

with open(outprefix + '_groups.txt', 'w') as outf:
    repeatsDB = dict((p,sum(repeatsDB[p])/len_data) for p in repeatsDB)
    sort_repeats = sorted(repeatsDB.items(), key=lambda x: x[1], reverse=True)
    limit_bottom = len_data * [0]
    if len(normalsDB):
        repeatsDB = normal_repeatsDB
        X_axis.append('NORMALS')
        limit_bottom.append(0)
        len_data += 1
        ind = np.arange(len_data)

    outf.write("name\t" + "\t".join(X_axis) + "\n")
    legendDB = dict()
    all_colors = ("b", "m", "g", "r", "y", "c", 'k')
    all_colors = ("navy", "green", "maroon",
                  "blue", "lime", "red",
                  "teal", "olive", "purple",
                  "aqua",  "yellow", "fuchsia",
                  "gray", "black")

    top_repeats = list()

    delta_repeats = dict()

    for repeat,color in zip(sort_repeats, all_colors):
        repeat = repeat[0]
        top_repeats.append(repeat)
        temp = [float(dataDB[p]['repeats'].get(repeat, '0.')) for p in X_axis if p in dataDB]
        delta_repeats[repeat] = temp
        if len(normalsDB):
            temp.append(np.mean(normal_repeatsDB.get(repeat, 0.)))
        outf.write('%s\t' %repeat + '\t'.join([str(p-repeatsDB[repeat]) for p in temp]) + '\n')    
        legendDB[repeat] = plt.bar(ind, temp, bottom=limit_bottom, color=color, ec=color, linewidth=0.1)
        limit_bottom = [p+q for p,q in zip(limit_bottom, temp)]
                
    # add-up to 1.0 with summing up all others       
    #temp = [1.0-p for p in limit_bottom]
    #legendDB['other'] = plt.bar(ind, temp, bottom=limit_bottom, color='silver', ec='0.75', linewidth=0.1)

    # what are 'others'
    #repeatsDB['other'] = sum([repeatsDB[p] for p in repeatsDB if p not in top_repeats])/len(dataDB)
    
    #outf.write('%s\t' % 'other' + '\t'.join([str(p-repeatsDB["other"]) for p in temp]) + '\n')



# plot 4
# line-graph for change in repeat-fraction
# use an running-average of three
AVG = 2
for repeat,color in zip(delta_repeats, all_colors):
    data = [np.mean(delta_repeats[repeat][i:i+AVG]) for i in xrange(0, len(delta_repeats[repeat])-AVG, AVG)]
    mdata = np.mean(data)
    data = [p-mdata for p in data]
    plt.plot(data, color=color, linestyle='-', label=repeat)
plt.tight_layout()
#plt.figlegend()
plt.savefig(outprefix + '_delta_hexrepeats.pdf', format="pdf", bbox_inches='tight')
plt.close()
"""

# plot 3 heatmap
# can you generate the clusters and save the files for viewing
cluster_method= "m"
"""
method== s : pairwise single-linkage clustering
method== m : pairwise maximum- (or complete-) linkage clustering
method== c : pairwise centroid-linkage clustering
method== a : pairwise average-linkage clustering
"""
cluster_dist= "c"
"""
dist== c : correlation;
dist== a : absolute value of the correlation;
dist== u : uncentered correlation;
dist== x : absolute uncentered correlation;
dist== s : Spearman s rank correlation;
dist== k : Kendall s;
dist== e : Euclidean distance;
dist== b : City-block distance.
"""


"""
if len(repeatsDB) > 1:
    record= Cluster.read(open(outprefix + '_cluster.txt'))
    genetree= record.treecluster(transpose=0, method=cluster_method, dist=cluster_dist)
    genetree.scale()
    exptree= record.treecluster(transpose=1, method=cluster_method, dist=cluster_dist)
    exptree.scale()
    record.save(outprefix + '_cluster', geneclusters=genetree,  expclusters=exptree)

else:
    print 'need adv counting for a heatmap'
"""
