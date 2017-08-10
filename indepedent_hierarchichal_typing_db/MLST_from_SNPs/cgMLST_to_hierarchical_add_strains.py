from time import sleep as sl
import random
import sys
from time import sleep as sl
import glob
import pandas as pd


profiles = sys.argv[1]

inaccs = glob.glob(sys.argv[2]+"*")

outfolder = sys.argv[3]

strains = sys.argv[4]

sets_dict = {}

for i in inaccs:
    acc=i.split("/")[-1].split("_")[0]
    # print acc
    inf = open(i,"r").read().strip("\n").splitlines()
    # print len(inf)
    sets_dict[acc]=inf

# for i in sets_dict:
#     print i,sets_dict[i][:10]

profiles = pd.read_csv(profiles,sep="\t",header=0)

tot_dict = {}
header = "ST"

nl = list(sets_dict.keys())
nl = [x for x in nl if x != 'cgMLST']
schemels = map(str,sorted(map(int,nl))) + ["cgMLST"]

for s in schemels:
    print s
    header += "\t"+str(s)+"_string\t"+str(s)+"_type"
    subset = sets_dict[s]

    # print subset[:10]

    new_profiles = profiles[["ST"]+subset]

    # print profiles.shape

    # print new_profiles.shape

    ## check if profile at x gene level is identical to preexisting one.


    new_lists = new_profiles.values.tolist()

    new_lists = [list(new_profiles.columns.values)] + new_lists

    # print len(new_lists)
    # print len(new_lists[0])

    profilefile = open(outfolder + s + "_gene_profiles.txt", "r").read().strip('\n').splitlines()
    oldlists = [x.split('\t') for x in profilefile]
    conversion = {"-".join()}

    # print len(new_lists),len(new_lists[0])

    count = 1
    out_profiles = {}
    conversion = {}
    new_old = {}
    test = []

    # print new_lists[0]

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

    # for i in new_out[:10]:
    #     print i[:10]
    #     sl(0.3)

    # print len(new_out)
    # print len(new_out[0])

    outfile = open(outfolder+s+"_gene_profiles.txt","")

    for i in new_out:
        outfile.write("\t".join(i)+ "\n")

    outfile.close()

    outfile = open(outfolder+s+"_gene_annotation.txt","")

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

outannot = open(outfolder+"all_scheme_annotation.txt","w")

outannot.write(header+'\tcgMLST_string\n')

for i in tot_dict:
    outannot.write(str(i)+'\t'+'\t'.join(tot_dict[i]) +'\t'+ tot_dict[i][-2]+'-'+i+'\n')
    # outannot.write(str(i) + '\t' + '\t'.join(tot_dict[i]) + '\n')

outannot.close()

def annot_strains(profiles,strains):
    outf = strains.replace("assignments.txt","hierarchical_assignments.txt")
    outf = open(outf,"w")
    strains = open(strains,"r").read().strip('\n').splitlines()
    outf.write("Strain\t"+header[3:]+'\n')
    for i in strains[1:]:
        col = i.split("\t")
        outf.write(col[0]+'\t'+"\t".join(tot_dict[col[1]])+'\n')
    outf.close()




annot_strains(tot_dict,strains)

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
