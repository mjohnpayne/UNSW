from time import sleep as sl

import sys

from collections import Counter

'''
For each level of hST (starting at 1445) if a strain is identical at that ST but different at the next less specific one (560 in this case) :
 
find the most common 560 gene st for a given 1445 st and change all 560 gene STs in that group to the most common

Then do the same sequentially for less specific ST schemes i.e. for 560 check and modify 420, for 420 check and modify 260 etc etc

STEPS:

1 - get dict of each cgMLST ST with corresponding list of 9 sub STs
2 - iterate over 1445 STs collecting 560 ST types for each 1445 ST
3 - if an 1445 ST has > 1 560 ST and there is one in majority convert all strings to have the majority 560 ST
4 - use above output strings to feed into same process for next pair of schemes


'''

def check_and_convert_sts(stdict,lower_scheme_pos):
    types_dict = {}
    for i in stdict:
        low_st = stdict[i][lower_scheme_pos]
        high_st = stdict[i][lower_scheme_pos-1]
        if stdict[i][lower_scheme_pos] not in types_dict:
            types_dict[low_st] = [high_st]
        else:
            types_dict[low_st].append(high_st)
    conv = {}
    for i in types_dict:
        c = Counter(types_dict[i])
        countls = c.most_common()
        if len(set(types_dict[i])) == 1:
            conv[i] = None
        elif countls[0][1]-countls[1][1] >= 0:
            conv[i] = countls[0][0]
        else:
            conv[i] = None
    # for i in conv:
    #     if conv[i] != None:
    #         print i,conv[i]
    #         sl(0.2)
    outdict = {}
    for i in stdict:
        low_st = stdict[i][lower_scheme_pos]
        high_st = stdict[i][lower_scheme_pos - 1]
        if conv[low_st] == None:
            outdict[i] = stdict[i]
        else:
            if high_st == conv[low_st]:
                outdict[i] = stdict[i]
            else:
                lst = stdict[i][:lower_scheme_pos-1] + [conv[low_st]] + stdict[i][lower_scheme_pos:]
                outdict[i] = lst
    return outdict


def run_all(dict,num):
    print num
    if num == 0:
        return dict
    else:
        dict1 = check_and_convert_sts(dict,num)
        num = num-1
        return run_all(dict1,num)




hST = open(sys.argv[1],"r").read().splitlines()

stdict = {}




### from snp base scheme
for i in hST[1:]:
    col = i.split("\t")
    if col[0] != "":
        stdict[col[0]] = col[2::2]
# print stdict

### from enterobase scheme

# for i in hST[1:]:
#     col = i.split("\t")
#     if col[0] != "":
#         stdict[col[0]] = col[1].split("-")
# print stdict




final_dict = run_all(stdict,9)

outnm = sys.argv[1].strip(".txt")+"_hcorrected0.txt"

outhST = open(outnm,"w")




###snp_based

header = hST[0].split('\t')

header = ["corr_"+x for x in header[2::2]]

outhST.write("Strain\tcorr_hier_string_ST\tcorr_final4_hier_string\t"+"\t".join(header)+"\n")

### enterobase based

# header = hST[0].split('\t')
#
# header = ["corr_"+x for x in header[3:]] + ["cgMLST"]
#
# outhST.write("cgMLST_type\tcorr_hier_string_ST\tcorr_final4_hier_string\t"+"\t".join(header)+"\n")




for i in final_dict.keys():
    outhST.write(i + '\t' + "-".join(final_dict[i]) + "\t" +"-".join(final_dict[i][-3:]) +"\t" + "\t".join(final_dict[i]) + '\n')

outhST.close()


# check_and_convert_sts(stdict,8)