from Bio import SeqIO
import sys

fasta_file = sys.argv[1]#"/Users/lukemn/Documents/prc/pABH1_13Illumina/velvet/velk40/contigs.fa"
subset_file = sys.argv[2]#"/Users/lukemn/Documents/prc/pABH1_13Illumina/velvet/velk40/contigRetrieve.txt"
output = sys.argv[3]#"/Users/lukemn/Documents/prc/pABH1_13Illumina/velvet/velk40/Retrieved.fa"

wanted = set()
with open(subset_file) as f:
    for line in f:
        #print line
        line = line.strip()
        if line != "":
            #print line
            wanted.add(line)

fasta_sequences = SeqIO.parse(open(fasta_file),"fasta")
with open(output, "w") as f:
    for seq in fasta_sequences:
        if seq.id in wanted:
            SeqIO.write(seq, f, "fasta")

fasta_sequences.close()
