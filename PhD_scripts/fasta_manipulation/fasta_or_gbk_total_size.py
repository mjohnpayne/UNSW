import sys
from Bio import SeqIO
from Bio import GenBank
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_nucleotide
from Bio.Alphabet import IUPAC

inp = sys.argv[1]

inf = SeqIO.parse(inp, "fasta")
totlen = 0
num = 0
for i in inf:
    totlen += len(i.seq)
    num += 1
print"Total Size of fastas: " + str(totlen)
print "Average Size of fasta: " + str(float(totlen)/num)


