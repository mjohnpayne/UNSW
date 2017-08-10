from time import sleep as sl

outf = open("/Users/michaelpayne/Documents/UNSW/Salmonella/LT2_genome/LT2_cgMLST_gene_conversions.txt","w")

for i in open("/Users/michaelpayne/Documents/UNSW/Salmonella/LT2_genome/LT2_cgMLST_gene_hits.txt","r").read().splitlines():
    col = i.split("\t")
    stm = col[0].split(" ")
    stm_acc = ""
    for i in stm:
        if "locus_tag" in i:
            stm_acc = i.strip("[locus_tag=").strip("]")
    cgmlst = col[1][:col[1].rindex("_")]
    outf.write(stm_acc+'\t'+cgmlst+'\n')

outf.close()