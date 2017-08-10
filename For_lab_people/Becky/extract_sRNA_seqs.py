from time import sleep as sl

from Bio import SeqIO
from Bio import SeqRecord

infile = SeqIO.parse("/Users/michaelpayne/Documents/UNSW/Students/Becky/GCF_000195715.1_ASM19571v1_genomic.fna","fasta")

genome = ""
for i in infile:
    genome = i

genseq = genome.seq
print len(genseq)


infile = open("/Users/michaelpayne/Documents/UNSW/Students/Becky/NC_002929_sRNA.out_annotated.txt","r").readlines()

outfile = open("/Users/michaelpayne/Documents/UNSW/Students/Becky/tohama_sRNAs.fasta","w")

for i in infile[6:]:
    col = i.split('\t')
    if "********" in i:
        break
    elif len(col) > 5:
        sRNA=genseq[int(col[9])-1:int(col[10])]
        if "<" in col[11]:
            sRNA=str(SeqRecord.SeqRecord(sRNA).reverse_complement().seq)
        sRNA = str(sRNA)
        outfile.write(">"+col[0]+"\n"+sRNA+"\n")

outfile.close()