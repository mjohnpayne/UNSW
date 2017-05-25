__author__ = 'mjohnpayne'

import sys
from Bio import SeqIO

inp = SeqIO.parse(sys.argv[1],'fasta')

for i in inp:
    outfile = sys.argv[2] + '/' + i.id + '.fasta'
    SeqIO.write(i,outfile,'fasta')
