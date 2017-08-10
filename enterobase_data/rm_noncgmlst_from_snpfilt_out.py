from time import sleep as sl

import glob

infiles = glob.glob("/Users/michaelpayne/Documents/UNSW/Salmonella/Multiple_SNP_calls_testing/snp_filt_outputs/*.csv")

in_ref_pos = open("/Users/michaelpayne/Documents/UNSW/Salmonella/data_aquisition/enterobase_data/cgMLST_definitions/LT2_MLST_position_data/LT2_allele_locations.txt","r")

posls = []




for i in in_ref_pos:
    col = i.strip('\n').split('\t')
    posls.append([int(col[1]),int(col[2])])

print posls[:10]

def within_cg(poslist,position):
    check = 0
    for i in poslist:
        if i[0] <= int(position) <= i[1]:
            check = 1
    return check


for i in infiles:
    name = i.split("/")[-1][:-9]
    inf = open(i,"r").readlines()
    print name
    outfile = open("/".join(i.split('/')[:-1]) + '/' + name + "_cgMLST_only.txt","w")
    for j in inf:
        if j[0] == "#" or j[0] == "f":
            outfile.write(j)
        else:
            col = j.strip('\n').split(',')
            pos = col[2]
            if within_cg(posls,pos) == 1:
                outfile.write(j)
    outfile.close()



