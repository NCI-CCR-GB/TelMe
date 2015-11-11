#merge
import os
import sys
import glob
import numpy as np
import multiprocessing as multiP
import subprocess as subP
from itertools import product
from scipy.stats import linregress, pearsonr, spearmanr
from optparse import OptionParser

parser = OptionParser()

parser.add_option("-f", dest="inpath",
                 help="/path/to/indir/ [%default]",
                 default=os.getcwd())

parser.add_option("-d", dest="data",
                  type="choice",
                  choices=('nci60', '1000g', 'plot', 'merge'),
                  help="data type is it 'nci60, 1000g, merge' or just need 'plot' [%default]",
                  default='1000g')

parser.add_option("-v", dest="value",
                  type="choice",
                  choices=('offt', 'chrm', 'offtcnt'),
                  help="which value to calculate correlation 'offt', 'chrm', 'offtcnt' [%default]",
                  default='offt')

parser.add_option("-m", dest='merge',
                  help='need to merge first [%default]',
                  action='store_true', default=False)

parser.add_option("--corr", dest='corr',
                  help='calculate genome vs exome numbers ony [%default]',
                  action='store_true', default=False)

(option, arg)= parser.parse_args()

data_value = {'offt' : 3, 'chrm' : 6, 'offtcnt' : 9}


def collect_matrix(s):
    return set(os.path.basename(p).split('_telcount_by_')[1] for p in glob.glob(os.path.join(s, "*_telcount_by_*.tsv")))


def summary_telcnt(s):
    subP.call(['python',
               os.path.join(os.path.dirname(sys.argv[0]), 'plot_telCNT.py'), '-f' , option.inpath, '-d', s, '-p', option.value])
    return

matrix = collect_matrix(option.inpath)
plot3D = list()
if option.data == 'nci60':
    # calc R2 between MMQPCR
    data1 = {'RPMI8226':0.104903901,'SR':0.147910047,'K562':0.318848784,'HL60':0.347084249,'CCRFCEM':0.60698061,'MOLT4':0.826795241,'SF295':0.154506851,'SNB19':0.215696539,'SNB75':0.237510368,'SF539':0.481706139,'U251':0.71016758,'SF268':1.371578386,'BT549':0.116924174,'MDAMB231':0.125665072,'T47D':0.164070382,'MDAMB468':0.184282268,'MCF7':0.239252154,'HS578T':0.399737576,'COLO205':0.213311666,'SW620':0.216553238,'KM12':0.255927101,'HCC2998':0.393384364,'HT29':0.463047027,'HCT116':0.491488212,'HCT15':0.800910246,'NCIH522':0.155292589,'EKVX':0.246122888,'NCIH226':0.248969967,'NCIH322M':0.282402414,'HOP62':0.33963512,'HOP92':0.445477351,'NCIH460':0.629361409,'A549ATCC':0.685670604,'NCIH23':2.500496309,'UACC62':0.218036691,'SKMEL2':0.329932565,'UACC257':0.337038227,'MALME3M':0.373990615,'SKMEL5':0.481831025,'MDAMB435':0.905810283,'M14':1.958366049,'SKMEL28':2.458716443,'LOXIMVI':11.62419,'OVCAR5':0.270016965,'IGROV1':0.303356597,'OVCAR3':0.374108667,'OVCAR8':0.474659411,'NCIADRRES':0.720260714,'SKOV3':0.754147662,'OVCAR4':0.967697883,'PC3':0.116381103,'DU145':0.352909736,'TK10':0.096153677,'UO31':0.129639628,'7860':0.157232698,'A498':0.265178878,'ACHN':0.220651497,'RXF393':0.49838021,'CAKI1':1.02882806,'SN12C':2.417144842}

    if option.merge:
        jobP = multiP.Pool(8)
        jobR = jobP.map_async(summary_telcnt, matrix)
        jobP.close()
        jobP.join()
        jobR.get()

    def normname(outdata):
        prefix = os.path.commonprefix([p[::-1] for p in outdata.keys()])[::-1]
        outdata = dict((p.split(prefix)[0].translate(None," .-_,()\/").upper(),q) for p,q in outdata.items())

        return outdata
    

    with open("r2-corr-%s.tsv" %option.value, "w") as outf:
        outf.write('\t'.join(['file1', 'MMQPCR', 'pearsonC', 'spearmanC', 'linearC']) + '\n')
        for each in matrix:
            each = each.replace('.tsv', '.txt', 1)
            flist = sorted(glob.glob(os.path.join(option.inpath, "summary_%s_%s"%(option.value, each))))
            if len(flist) == 1:
                # calculate
                data2 = dict()
                with open(flist[0]) as inf:
                    for line in inf:
                        line = line.strip().split()
                        data2[line[0]] = float(line[data_value[option.value]])
                data2 = normname(data2)
                order = sorted(data1.items(), key=lambda x:x[1])
                data = [list(), list()]
                incomplete=False
                for K,V in order:
                    if K not in data2:
                        incomplete=True
                        print flist[0], "missing", K
                    else:
                        data[0].append(V)
                        data[1].append(data2.get(K,""))

                # calc linear correlation
                slope, intercept, r_value, p_value, std_err = linregress(data)
                r2 = r_value ** 2
                # calc pearson correlation
                pC, pp_value = pearsonr(*data)
                # calc spearman correlation
                sC, sp_value = spearmanr(*data)
                
                print pC, sC, r2, each,
                if incomplete:
                    print "INCOMPLETE"
                else:
                    print
                    
                outf.write("\t".join([flist[0], 'MMQPCR', str(pC), str(sC), str(r2), str(incomplete)]) + "\n")
                plot3D.append([each, pC])                
            else:
                print "missing file for", each, len(flist)
                outf.write("\t".join(["missing file", "_".join(each)]) + "\n")

elif option.data == '1000g':
    #calculate R2 for genome vs exome

    if option.corr:
        #MM = [p for p in matrix]
        import random
        Xdata = list() # Oe/Fe ontarget/offtarget
        Ydata = list() # Te/Tg exon_cnt/genome_cnt

        for V in matrix:
            flist = glob.glob(os.path.join(option.inpath,"*exome*%s"%V))
            flist2 = glob.glob(os.path.join(option.inpath,"*low_c*%s"%V))
            if len(flist) == len(flist2) and len(flist) > 10:
                #print 'located', V, len(flist), len(flist2)
                flist.sort()
                flist2.sort()
                for genome_file, exome_file in zip(flist2, flist):
                    if os.path.basename(genome_file).split(".")[0] != os.path.basename(exome_file).split(".")[0]:
                        print os.path.basename(genome_file).split(".")[0], os.path.basename(exome_file).split(".")[0]
                        break

                    telE = list() # relative tel counts
                    ontE = list() # on target reads fraction
                    with open(exome_file) as inf:
                        line = inf.readline()
                        while not line.startswith('RG'):
                            line = inf.readline()

                        headerDB=dict((p,q) for q,p in enumerate(line.strip().split()))
                        # data lines from now on
                        line = inf.readline().strip()
                        while len(line):
                            line = line.split('\t')
                            telE.append(float(line[headerDB['TELOMERE']])/(float(line[headerDB['ALL']])-float(line[headerDB['ONTARGET']])))
                            ontE.append(float(line[headerDB['ONTARGET']])/float(line[headerDB['ALL']]))
                            line = inf.readline().strip()

                    telG = list() # relative tel counts
                    with open(genome_file) as inf:
                        line = inf.readline()
                        while not line.startswith('RG'):
                            line = inf.readline()

                        # data lines from now on
                        line = inf.readline().strip()
                        while len(line):
                            line = line.split('\t')
                            telG.append(float(line[headerDB['TELOMERE']])/float(line[headerDB['ALL']]))
                            line = inf.readline().strip()

                    ontE = np.mean(ontE)
                    if ontE > 0.1:
                        Xdata.append(ontE)
                        Ydata.append( np.mean(telE) / np.mean(telG) )

        import matplotlib.pyplot as plt
        #print Xdata
        #print Ydata
        plt.plot(Xdata, Ydata, 'bo')
        plt.show()
        
        sys.exit()
    
    if option.merge:
        jobP = multiP.Pool(8)
        jobR = jobP.map_async(summary_telcnt, ["%s*%s"%p for p in product(('.exome.', '.low_coverage.'), matrix)])
        jobP.close()
        jobP.join()
        jobR.get()

    with open("r2-corr-%s.tsv" %option.value, "w") as outf:
        outf.write('\t'.join(['file1', 'file2', 'pearsonC', 'spearmanC', 'linearC']) + '\n')
        for each in matrix:
            each = each.replace('.tsv', '.txt', 1)
            flist = sorted(glob.glob(os.path.join(option.inpath, "summary_%s_*%s"%(option.value, each))))
            if len(flist) == 2:
                # calculate
                data1 = dict()
                with open(flist[0]) as inf:
                    for line in inf:
                        line = line.strip().split()
                        data1[line[0].split('.')[0]] = float(line[data_value[option.value]])
                data2 = dict()
                with open(flist[1]) as inf:
                    for line in inf:
                        line = line.strip().split()
                        data2[line[0].split('.')[0]] = float(line[data_value[option.value]])
                        
                order = sorted(data1.items(), key=lambda x:x[1])
                data = [list(), list()]
                incomplete=False
                for K,V in order:
                    if K not in data2:
                        incomplete=True
                        print >> sys.stderr, flist[0], "missing", K
                    else:
                        data[0].append(V)
                        data[1].append(data2.get(K,""))

                #print data
                # calc linear correlation
                slope, intercept, r_value, p_value, std_err = linregress(data)
                r2 = r_value ** 2
                # calc pearson correlation
                pC, pp_value = pearsonr(*data)
                # calc spearman correlation
                sC, sp_value = spearmanr(*data)
                
                print pC, sC, r2, each, 
                if incomplete:
                    print 'incomplete'
                else:
                    print
                                
                outf.write("\t".join([flist[0], flist[1], str(pC), str(sC), str(r2)]) + "\n")
                plot3D.append([each, pC])                
            else:
                print "missing file for", each, len(flist)
                outf.write("\t".join(["missing file", "_".join(each)]) + "\n")

elif option.data == 'merge':
    # just run merge:
    jobP = multiP.Pool(8)
    jobR = jobP.map_async(summary_telcnt, matrix)
    jobP.close()
    jobP.join()
    jobR.get()


                    
if option.data == 'plot':
    if os.path.isfile("r2-corr.tsv"):
        # read the file
        plot3D = list()
        with open("r2-corr.tsv") as inf:
            for line in inf:
                if not line.startswith('missing'):
                    line = line.strip().split("\t")
                    D = ".txt"
                    if '_dups.txt' in line[0]:
                        D = "_dups.txt"

                    V = line[0].split(D)[0].split("_")[-4:] + [D]
                    plot3D.append([V, float(line[-1])])

    from mpl_toolkits.mplot3d import Axes3D
    import matplotlib
    matplotlib.use('agg')
    import matplotlib.pyplot as plt

    # keep C, T, Q, R2
    data = dict()
    for V, R2 in plot3D:
        C, T, Q, CC, D = V
        K = "_".join([T, D.split('.txt')[0], C])
        data.setdefault(K, list([[], [], []]))
        [data[K][i].append(float(j)) for i,j in enumerate([Q,CC,R2])]

    data = sorted(data.items(), key=lambda x: x[1][-1], reverse=True)

    fig = plt.figure(figsize=(9,11))
    for i, each in enumerate(data, 1):
        ax = fig.add_subplot(4,3,i, projection='3d')
        KEY, V = each
        X, Y, Z = V
        ax.scatter(X, Y, Z)

        ax.set_xlabel('QUAL')
        ax.set_xlim((0,40))
        #ax.set_xticks(sorted(set(X)))
        #ax.set_xticklabels(sorted(set(X)))
        ax.set_ylabel('CC')
        ax.set_ylim((0,5))
        #ax.set_xticks(sorted(set(Y)))
        #ax.set_yticklabels(sorted(set(Y)))
        ax.set_zlabel('R2')
        ax.set_zticklabels((0, 0.2, 0.4, 0.6, 0.8, 1.))
        ax.set_zlim((0,1.))
        ax.set_title(KEY+"_%.3f"%(sum(Z)/len(Z)))
    plt.tight_layout()
    plt.savefig("r2-corr.pdf",format='pdf')
    #plt.show()
