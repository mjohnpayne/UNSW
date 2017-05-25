
from Bio import SeqIO
from Bio.Seq import Seq
import sys

inseq = SeqIO.parse(sys.argv[1],'fasta')

outseq = open(sys.argv[1][:-6] + "_dedup.fasta",'w')

# outseq = open("/".join(sys.argv[1].split('/')[:-1])+"/"+sys.argv[1].split('/')[-1][:-6] + "_dedup.fasta",'w')


seqs = {}
for i in inseq:
    seqs[i.id] = str(i.seq)

for i in seqs:
    outseq.write('>' + i + '\n' + seqs[i] + '\n')

outseq.close()
