from time import sleep as sl
from Bio import SeqIO

def cat_multifastas(f1,f2,out):
    in1 = SeqIO.parse(f1,"fasta")
    in1d = {}
    for i in in1:
        in1d[i.id] = str(i.seq)

    in2 = SeqIO.parse(f2, "fasta")
    in2d = {}
    for i in in2:
        in2d[i.id] = str(i.seq)

    outf = open(out,"w")
    for i in in1d:
        for j in in2d:
            if i == j:
                s = in1d[i] + in2d[i]
                outf.write(">"+i+'\n'+s+'\n')
    outf.close()



# cat_multifastas("/Users/michaelpayne/Documents/UNSW/Salmonella/new_typing_scheme/cgMLST/new_cgmlst_strains_concat_snps.fasta","/Users/michaelpayne/Documents/UNSW/Salmonella/new_typing_scheme/cgMLST_STM/stm_int_cgmlst_concat_snps.fasta","/Users/michaelpayne/Documents/UNSW/Salmonella/new_typing_scheme/cgMLST_STM/all_cgmlst_concat_snps.fasta")

cat_multifastas("/Users/michaelpayne/Documents/UNSW/Salmonella/new_typing_scheme/cgMLST/new_cgmlst_strains_concat_snps.fasta","/Users/michaelpayne/Documents/UNSW/Salmonella/new_typing_scheme/cgMLST_stm_genes/cgMLST_stm_genes_concat_snps.fasta","/Users/michaelpayne/Documents/UNSW/Salmonella/new_typing_scheme/cgMLST_stm_genes/stm_genes_cgmlst_concat_snps.fasta")