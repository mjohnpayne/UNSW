from time import sleep as sl

import glob

infiles = glob.glob("/Users/michaelpayne/Documents/UNSW/Salmonella/Multiple_SNP_calls_testing/nullarbor_random_20_snp_files/*")

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
    outfile = open(i[:-4] + "cg_only.txt","w")
    outfile.write(inf[0])
    for j in inf[1:]:
        col = j.strip('\n').split('\t')
        pos = col[1]
        if within_cg(posls,pos) == 1:
            outfile.write(j)
    outfile.close()



