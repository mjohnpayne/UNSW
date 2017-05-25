__author__ = 'mjohnpayne'
from Bio import SeqIO

infile = SeqIO.parse("/Users/mjohnpayne/Documents/PhD/wt_genome/pmfa1_annot_scaf.fasta","fasta")
sizes = []
for i in infile:
    print len(i.seq), i.id
    sizes.append(len(i.seq))

print sorted(sizes)
