from time import sleep as sl

import pandas as pd
import numpy as np


import glob

infiles = glob.glob("/Users/michaelpayne/Documents/UNSW/Salmonella/Enterobacteriaceae_core/cgMLST_loci_to_use_*.txt")

for i in infiles:
    perc = i.split('/')[-1][-6:-4]
    print perc

    cgMLST_subset = open(i,"r").read()

    subset = cgMLST_subset.split('\n')
    print subset[:10]


    profiles = pd.read_csv("/Users/michaelpayne/Documents/UNSW/Salmonella/data_aquisition/enterobase_data/cgMLST_definitions/cgMLST-v2-profiles-8-6-17",sep="\t",header=0)

    new_profiles = profiles[["ST"]+subset]

    print profiles.shape

    print new_profiles.shape

    new_lists = new_profiles.values.tolist()

    new_lists = [list(new_profiles.columns.values)] + new_lists

    print len(new_lists)
    print len(new_lists[0])


    # print len(new_lists),len(new_lists[0])

    count = 1
    out_profiles = {}
    conversion = {}
    new_old = {}
    test = []

    print new_lists[0]

    for i in new_lists[1:]:
        st = i[0]
        unique = "-".join(map(str, i[1:]))
        if unique not in test:
            out_profiles[count] = [str(count)]+unique.split("-")
            conversion[unique] = count
            new_old[count] = [st]
            test.append(unique)
            count += 1
        else:
            c = conversion[unique]
            new_old[c].append(st)

    new_out = [new_lists[0]]+[out_profiles[x] for x in sorted(out_profiles.keys())]

    for i in new_out[:10]:
        print i[:10]
        sl(0.3)

    print len(new_out)
    print len(new_out[0])

    outfile = open("/Users/michaelpayne/Documents/UNSW/Salmonella/Enterobacteriaceae_core/entbact_"+perc+"_MLST_profiles.txt","w")

    for i in new_out:
        outfile.write("\t".join(i)+ "\n")

    outfile.close()

    outfile = open("/Users/michaelpayne/Documents/UNSW/Salmonella/Enterobacteriaceae_core/cgMLST_to_entbac_"+perc+"_MLST.txt","w")

    outfile.write("cgMLSTv2 SeqType\tentbacMLST SeqType\n")

    for i in new_old:
        for j in new_old[i]:
            outfile.write(str(j)+"\t" + str(i) + "\n")

    outfile.close()

annots = glob.glob("/Users/michaelpayne/Documents/UNSW/Salmonella/Enterobacteriaceae_core/cgMLST_to_entbac_conv/cgMLST_to_entbac_*_MLST.txt")
outannot = open("/Users/michaelpayne/Documents/UNSW/Salmonella/Enterobacteriaceae_core/cgMLST_to_entbac_conv/cgMLST_to_entbac_all_for_tree.txt","w")
outd = {}
header = "cgMLST_type"
for i in annots:
    perc = i.split('/')[-1][-11:-9]
    print i
    inf = open(i,"r").readlines()
    header += "\tentbac_"+perc
    for j in inf[1:]:
        col = j.strip('\n').split('\t')
        if col[0] not in outd:
            outd[col[0]] = [col[1]]
        else:
            outd[col[0]].append(col[1])

outannot.write(header + "\n")

for i in map(str,sorted(map(int,outd.keys()))):
    outannot.write(i + "\t" + "\t".join(outd[i])+'\n')

outannot.close()