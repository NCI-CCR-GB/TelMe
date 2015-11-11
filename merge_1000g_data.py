#merge
import os
import sys
import glob
import numpy as np
from math import cosh

def CORR(TELR, EXPF, ALLT, ONT, F_ONT, CC):
    if CC == '1':
        # use an exponential decay correction
        Vt = 10**6 * (TELR / (ALLT-ONT)) / 10**(-1*F_ONT**2)
        # each repeat is 6-bases long, so you could multiply by 6
        Vl = 10**6 * (EXPF / (ALLT-ONT)) / 10**(-1*F_ONT**2)
    elif CC == '2':
        # use a exponential correction
        Vt = 10**6 * (TELR / (ALLT-ONT)) * (1 - F_ONT)**-1
        # each repeat is 6-bases long, so you could multiply by 6
        Vl = 10**6 * (EXPF / (ALLT-ONT)) * (1 - F_ONT)**-1
    elif CC == '3':
        # use an hyperbolic cosine correction
        Vt = 10**6 * (TELR / (ALLT-ONT)) * cosh(F_ONT*5)
        # each repeat is 6-bases long, so you could multiply by 6
        Vl = 10**6 * (EXPF / (ALLT-ONT)) * cosh(F_ONT*5)
    elif CC == '4':
        # use an exponential correction
        Vt = 10**6 * (TELR / (ALLT-ONT)) * 10**(F_ONT**2)**4
        # each repeat is 6-bases long, so you could multiply by 6
        Vl = 10**6 * (EXPF / (ALLT-ONT)) * 10**(F_ONT**2)**4
    return Vt, Vl


plot3D = list()
if not len(sys.argv[1]):
    for mapped in glob.iglob("*.mapped.*"):
        fpath, mapped = os.path.split(mapped)
        print mapped
        COUNT, TYPE, QUAL, CC = mapped.split('_telcount_by_')[1].split('.tsv')[0].split('_')

        unmapped = os.path.join(fpath, mapped.replace('.mapped.', '.unmapped.', 1))
        merged = os.path.join(fpath, mapped.replace('.mapped.', '.merged.', 1))

        RGdata = dict()
        outlist = list()
        TelRepeats = dict()
        for fname in (mapped, unmapped):
            with open(fname) as inf:
                line = inf.readline()
                while not line.startswith('RG'):
                    outlist.append(line)
                    line = inf.readline()

                header = line.strip().split('\t')
                #RG	ALL	ONTARGET	TELOMERE	CHRM	TEL-UNMAPPED	CNT_TEL_LEN	RTPMO-w/offt	RTPKM-w/chrM	RTPMO-w/offtCNT	TELOMERE_SEQ_FOUND
                line = inf.readline().strip()
                while len(line):
                    #print line
                    line = line.split('\t')
                    RG = line[0]
                    RGdata.setdefault(RG, list()).append([float(p) for p in line[1:7]])
                    TelRepeats.setdefault(RG, dict())

                    telcnt = int(line[6])
                    for tel in line[10:]:
                        R, F = tel.split(":")
                        F = round(telcnt*float(F))
                        TelRepeats[RG][R] = TelRepeats[RG].get(R,0) + F

                    line = inf.readline().strip()

        outlist.append('\t'.join(['RG', 'ALL', 'ONTARGET', 'TELOMERE', 
                                  'CHRM', 'TEL-UNMAPPED', 'CNT_TEL_LEN', 
                                  'RTPMO-w/offt',
                                  'RTPKM-w/chrM',
                                  'RTPMO-w/offtCNT',
                                  'TELOMERE_SEQ_FOUND']) + '\n')

        for_ont = list() # calc mean based on off-target tel read cnt
        for_chrm = list() # calc mean based on chrM ratio
        for_ontC = list() # calc mean based on off-target telomere repear cnt
        TelRepeats_all = dict()

        for RG in RGdata:
            ALLT, ONT, TELR, CHRM, TELUN, EXPF = np.sum(RGdata[RG], axis=0)

            # ASSUME: the better the capture efficieny
            # the worse the telomeric fraction
            # i.e. if your on-target is 0 you should assume
            # the the telomeric count represents WGS
            # so take it as is.
            F_ONT = ONT / ALLT

            Vt, Vl = CORR(TELR, EXPF, ALLT, ONT, F_ONT, CC)

            try:
                Vm = 10**3 * (TELR / CHRM)
            except ZeroDivisionError:
                Vm = 0

            for_ont.append(Vt)
            for_ontC.append(Vl)
            for_chrm.append(Vm)

            TelRepeats_RG = dict()
            # combine reverse-complement sequences
            for T in TelRepeats[RG].keys():
                tcnt = TelRepeats[RG][T]
                # look for rev-comp sequences
                if T.startswith('CCC') and T != 'CCCGGG':
                    T = revcomp(T)

                TelRepeats_all[T] = TelRepeats_all.get(T, 0) + tcnt
                TelRepeats_RG[T] = TelRepeats_RG.get(T, 0) + tcnt

            TRcnt = float(sum(TelRepeats_RG.values()))
            TR = sorted(TelRepeats_RG.items(), key=lambda x: x[1], reverse=True)
            TRstr = "\t".join(["%s:%.4f" %(p, q/TRcnt) for p,q in TR])

            outlist.append('\t'.join([RG, '%i'%ALLT, '%i'%ONT, '%i'%TELR, '%i'%CHRM,  '%i'%TELUN, '%i'%EXPF,  str(Vt), str(Vm), str(Vl), TRstr ]) + '\n')

        TRcnt = float(sum(TelRepeats_all.values()))
        TR = sorted(TelRepeats_all.items(), key=lambda x: x[1], reverse=True)
        TRstr = "\t".join(["%s:%.4f" %(p, q/TRcnt) for p,q in TR if q/TRcnt > 0.001])

        sname = os.path.basename(merged)
        outlist.append('\n' + '\t'.join([sname, 'MEAN+-SD',
                                         'RTPMO-w/offt', 
                                         str(np.mean(for_ont)),
                                         str(np.std(for_ont)), 
                                         'RTPKM-w/chrM', 
                                         str(np.mean(for_chrm)),
                                         str(np.std(for_chrm)),
                                         'RTPMO-w/offtCNT',
                                         str(np.mean(for_ontC)),
                                         str(np.std(for_ontC)),
                                         TRstr])+'\n')
        with open(merged, 'w') as outf:
            outf.writelines(outlist)

elif sys.argv[1].upper() == 'NCI60':
    # calc R2 
    #calculate R2 for genome vs exome
    from scipy.stats import linregress
    from itertools import product

    data1 = {'RPMI8226':0.104903901,'SR':0.147910047,'K562':0.318848784,'HL60':0.347084249,'CCRFCEM':0.60698061,'MOLT4':0.826795241,'SF295':0.154506851,'SNB19':0.215696539,'SNB75':0.237510368,'SF539':0.481706139,'U251':0.71016758,'SF268':1.371578386,'BT549':0.116924174,'MDAMB231':0.125665072,'T47D':0.164070382,'MDAMB468':0.184282268,'MCF7':0.239252154,'HS578T':0.399737576,'COLO205':0.213311666,'SW620':0.216553238,'KM12':0.255927101,'HCC2998':0.393384364,'HT29':0.463047027,'HCT116':0.491488212,'HCT15':0.800910246,'NCIH522':0.155292589,'EKVX':0.246122888,'NCIH226':0.248969967,'NCIH322M':0.282402414,'HOP62':0.33963512,'HOP92':0.445477351,'NCIH460':0.629361409,'A549ATCC':0.685670604,'NCIH23':2.500496309,'UACC62':0.218036691,'SKMEL2':0.329932565,'UACC257':0.337038227,'MALME3M':0.373990615,'SKMEL5':0.481831025,'MDAMB435':0.905810283,'M14':1.958366049,'SKMEL28':2.458716443,'LOXIMVI':11.62419,'OVCAR5':0.270016965,'IGROV1':0.303356597,'OVCAR3':0.374108667,'OVCAR8':0.474659411,'NCIADRRES':0.720260714,'SKOV3':0.754147662,'OVCAR4':0.967697883,'PC3':0.116381103,'DU145':0.352909736,'TK10':0.096153677,'UO31':0.129639628,'7860':0.157232698,'A498':0.265178878,'ACHN':0.220651497,'RXF393':0.49838021,'CAKI1':1.02882806,'SN12C':2.417144842}

    with open("r2-corr.tsv", "w") as outf:
        for each in product(('0.6', '0.9'),
                            ('simple', 'adv'),
                            ('0', '15', '30'),
                            ('1', '2', '3', '4'),
                            ('.txt', '_dups.txt')):
            flist = sorted(glob.glob("summary*_%s_%s_%s_%s%s"%each))
            if len(flist) == 1:
                # calculate
                data2 = dict()
                with open(flist[0]) as inf:
                    for line in inf:
                        line = line.strip().split()
                        data2[line[0].split('_')[0]] = float(line[3])

                order = sorted(data1.items(), key=lambda x:x[1])
                data = [list(), list()]
                for K,V in order:
                    if K not in data2:
                        print flist[0], "missing", K
                    else:
                        data[0].append(V)
                        data[1].append(data2.get(K,""))

                try:
                    slope, intercept, r_value, p_value, std_err = linregress(data)
                    r2 = r_value ** 2
                    print r2, each
                    outf.write("\t".join([flist[0], 'MMQPCR', str(r2)]) + "\n")
                    plot3D.append([each, r2])
                except TypeError:
                    print flist[0]            
            else:
                print "missing file for", each
                outf.write("\t".join(["missing file", "_".join(each)]) + "\n")

elif sys.argv[1].upper() != 'PLOT':
    #calculate R2 for genome vs exome
    from scipy.stats import linregress
    from itertools import product
    with open("r2-corr.tsv", "w") as outf:
        for each in product(('0.6', '0.9'),
                            ('simple', 'adv'),
                            ('0', '15', '30'),
                            ('1', '2', '3', '4')
                            ('.txt', '_dups.txt')):
            flist = sorted(glob.glob("summary*_%s_%s_%s_%s%s"%each))

            if len(flist) == 2:
                # calculate
                data1 = dict()
                with open(flist[0]) as inf:
                    for line in inf:
                        line = line.strip().split()
                        data1[line[0].split('.')[0]] = float(line[3])
                data2 = dict()
                with open(flist[1]) as inf:
                    for line in inf:
                        line = line.strip().split()
                        data2[line[0].split('.')[0]] = float(line[3])
                        
                order = sorted(data1.items(), key=lambda x:x[1])
                data = [list(), list()]

                for K,V in order:
                    data[0].append(V)
                    data[1].append(data2.get(K,""))

                #print data
                slope, intercept, r_value, p_value, std_err = linregress(data)
                r2 = r_value ** 2
                print r2, each
                outf.write("\t".join([flist[0], flist[1], str(r2)]) + "\n")
                plot3D.append([each, r2])                
            else:
                print "missing file for", each
                outf.write("\t".join(["missing file", "_".join(each)]) + "\n")

   
if len(plot3D) or sys.argv[1].upper() == 'PLOT':
    if not len(plot3D) and os.path.isfile("r2-corr.tsv"):
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
    elif not len(plot3D):
        print "run merge first"
        sys.exit()

    from mpl_toolkits.mplot3d import Axes3D
    import matplotlib.pyplot as plt

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    # keep C, T, Q, R2
    data = {'adv': list([[], [], [], []]),
            'simple': list([[], [], [], []]),
            'simple_dups': list([[], [], [], []]),
            'adv_dups': list([[], [], [], []])}

    for V, R2 in plot3D:
        C, T, Q, CC, D = V
        [data[T+D.split('.txt')[0]][i].append(float(j)) for i,j in enumerate([C,Q,CC,R2])]

    legendDB = dict()
    for LEG, KEYS in zip((('r', 'o'), ('b', '^'), ('g','s'), ('k','p')),
        data):
        print LEG, KEYS
        X,Y,Z,S = data[KEYS]
        legendDB[ax.scatter(X, Y, Z, s=[p*100 for p in S], color=LEG[0], marker=LEG[1])] = KEYS

    ax.set_xlabel('COUNT')
    ax.set_xlim((0.6, 0.9))
    ax.set_ylabel('QUAL')
    ax.set_zlabel('CC')

    ax.legend(legendDB.items())

    plt.show()
