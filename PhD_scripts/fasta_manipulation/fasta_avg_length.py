__author__ = 'mjohnpayne'


from Bio import SeqIO
import glob
import sys

infasta = glob.glob(sys.argv[1] + "/*.fasta")

def avg_fasta(fs):
    tmp = SeqIO.parse(fs,"fasta")
    ln = 0
    c = 0
    for j in tmp:
        ln += len(j.seq)
        c +=1
    return float(ln)/c

for i in infasta:
    name = i.split('/')[-1][:-6]
    print name,avg_fasta(i)