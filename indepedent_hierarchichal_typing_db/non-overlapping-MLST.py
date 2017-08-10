from time import sleep as sl
import random

incgmlst_loci = "/Users/michaelpayne/Documents/UNSW/Salmonella/Enterobacteriaceae_core/cgMLST_accessions.txt"

from time import sleep as sl

import pandas as pd
import numpy as np


import glob
random.seed(1234)

def select_random(lst,num):
    pos = random.sample(xrange(len(lst)), num)
    newlis = []
    for i in pos:
        newlis.append(lst[i])
    for i in newlis:
        lst.remove(i)
    return newlis,lst

def select_genes(cgmlst_loci):
    genes7 = [x for x in open("/Users/michaelpayne/Documents/UNSW/Salmonella/Enterobacteriaceae_core/STMMW_7_gene_accs.txt").read().split('\r')]
    cgmlst_genes = [x for x in open(cgmlst_loci).read().split('\r')]
    genes=cgmlst_genes
    set_dict = {}
    print len(genes)
    for i in genes7:
        genes.remove(i)
    # set_dict["7"] = genes7
    print len(genes)
    genes20,genes = select_random(cgmlst_genes,20)
    set_dict["20"] = genes20
    print len(genes)
    genes50,genes = select_random(cgmlst_genes,50)
    set_dict["50"] = genes50
    print len(genes)
    genes100,genes = select_random(cgmlst_genes,100)
    set_dict["100"] = genes100
    print len(genes)
    genes250,genes = select_random(cgmlst_genes,250)
    set_dict["250"] = genes250
    print len(genes)
    genes500,genes = select_random(cgmlst_genes,500)
    set_dict["500"] = genes500
    print len(genes)
    genes750,genes = select_random(cgmlst_genes,750)
    set_dict["750"] = genes750
    print len(genes)
    genes1000,genes = select_random(cgmlst_genes,1000)
    set_dict["1000"] = genes1000
    print len(genes)
    return set_dict

sets_dict = select_genes(incgmlst_loci)

profiles = pd.read_csv("/Users/michaelpayne/Documents/UNSW/Salmonella/data_aquisition/enterobase_data/cgMLST_definitions/cgMLST-v2-profiles-8-6-17",sep="\t",header=0)

tot_dict = {}
header = "ST"
for s in map(str,sorted(map(int,sets_dict.keys()))):
    print s
    header += "\t"+str(s)+"_string\t"+str(s)+"_type"
    subset = sets_dict[s]

    print subset[:10]

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

    outfile = open("/Users/michaelpayne/Documents/UNSW/Salmonella/new_typing_scheme/MLST_schemes/"+s+"_gene_MLST_scheme.txt","w")

    for i in new_out:
        outfile.write("\t".join(i)+ "\n")

    outfile.close()

    outfile = open("/Users/michaelpayne/Documents/UNSW/Salmonella/new_typing_scheme/MLST_schemes/cgMLST_to_entbac_"+s+"_gene_MLST.txt","w")

    outfile.write("cgMLSTv2 SeqType\t"+s+" gene MLST SeqType\n")

    for i in new_old:
        for j in new_old[i]:
            outfile.write(str(j)+"\t" + str(i) + "\n")
            if str(j) not in tot_dict:
                tot_dict[str(j)] = [str(i),str(i)]
            else:
                next_level = tot_dict[str(j)][-2] + "-" + str(i)
                tot_dict[str(j)] += [next_level,str(i)]
    outfile.close()

outannot = open("/Users/michaelpayne/Documents/UNSW/Salmonella/new_typing_scheme/MLST_schemes/cgMLST_to_new_all_for_tree.txt","w")

outannot.write(header+'\tcgMLST_string\n')

for i in tot_dict:
    outannot.write(str(i)+'\t'+'\t'.join(tot_dict[i]) +'\t'+ tot_dict[i][-2]+'-'+i+'\n')
outannot.close()

# annots = glob.glob("/Users/michaelpayne/Documents/UNSW/Salmonella/new_typing_scheme/MLST_schemes/cgMLST_to_entbac_*_gene_MLST.txt")
#
# outd = {}
# header = "cgMLST_type"
# for i in annots:
#     perc = i.split('/')[-1][-11:-9]
#     print i
#     inf = open(i,"r").readlines()
#     header += "\tentbac_"+perc
#     for j in inf[1:]:
#         col = j.strip('\n').split('\t')
#         if col[0] not in outd:
#             outd[col[0]] = [col[1]]
#         else:
#             outd[col[0]].append(col[1])
#
# outannot.write(header + "\n")
#
# for i in map(str,sorted(map(int,outd.keys()))):
#     outannot.write(i + "\t" + "\t".join(outd[i])+'\n')
#
# outannot.close()