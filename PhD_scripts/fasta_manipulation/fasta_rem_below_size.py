__author__ = 'mjohnpayne'

from Bio import SeqIO
import sys

inf = SeqIO.parse(sys.argv[1],'fasta')

outfolder = sys.argv[2]

cutoff = sys.argv[3]

outfile = outfolder + sys.argv[1].split('/')[-1].replace('.fasta','_rm_small.fasta')

out = []
for i in inf:
    if len(i.seq) > int(cutoff):
        out += [i]

SeqIO.write(out,outfile,'fasta')