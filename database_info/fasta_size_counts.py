__author__ = 'mjohnpayne'

from Bio import SeqIO
import sys

inf = SeqIO.parse("/Users/mjohnpayne/Documents/UNSW/Salmonella/data_aquisition/enterobase_data/concat_MLST_seqs.fasta",'fasta')

out = {}
for i in inf:
    if len(i.seq) not in out:
        out[len(i.seq)] = [1,[i.id]]
    else:
        out[len(i.seq)][0] +=1
        out[len(i.seq)][1].append(i.id)


for i in sorted(out.keys()):
    print i,out[i][0],out[i][1][:5]