from time import sleep as sl

import glob

inls = glob.glob("/Users/michaelpayne/Documents/UNSW/Salmonella/new_typing_scheme/cgMLSTv2_stm/stm_scheme/dbstarter/hierMLST_stm_loci_accession_lists/*")

outfolder = "/Users/michaelpayne/Documents/UNSW/Salmonella/new_typing_scheme/cgMLSTv2_stm/stm_scheme/dbstarter/hierMLST_stm_hierarchical_ST/"

for i in inls:
    last = i.split("/")[-1]
    id = last[:last.find("_")]

    lis = open(i,"r").read().splitlines()
    profs = ["1"]*len(lis)

    outprof = open(outfolder+id+"_gene_profiles.txt","w")
    outannot = open(outfolder+id+"_gene_annotation.txt","w")
    outannot.write("cgMLST_type\t"+id+"_gene_MLST_seqtype\n1\t1\n")
    outprof.write("ST\t" + "\t".join(lis) + "\n" + "1\t" + "\t".join(profs) + "\n")
    outprof.close()
    outannot.close()