import sys
from Bio import SeqIO
from Bio import GenBank
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_nucleotide
from Bio.Alphabet import IUPAC


def fastasize(input):
    inf = SeqIO.parse(input, "fasta")
    for i in inf:
        print i.id + '\t' + str(len(i.seq))

def gbsize(input):
    inf = SeqIO.parse(input, "genbank")
    for i in inf:
        print i.id + '\t' + str(len(i.seq))


inp = sys.argv[1]

if '.fa' in inp or '.fasta' in inp:
    fastasize(inp)
elif '.gb' in inp or '.gbk' in inp:
    gbsize(inp)



