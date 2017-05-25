__author__ = 'mjohnpayne'


from Bio import SeqIO
import sys
import math


infile = sys.argv[1]

splitno = int(sys.argv[2])

infasta = SeqIO.parse(infile,"fasta")

fastals = []

for i in infasta:
    fastals.append(i)

#outfile = open(infile.strip('.fasta') + '_' + str(fileno) + '.fasta','w')

filelen = float(len(fastals))

st = 0
en = splitno
fileno = 1

while st < filelen:
    if en > filelen:
        SeqIO.write(fastals[st:],infile.strip('.fasta') + '_' + str(fileno) + '.fasta','fasta')
    else:
        SeqIO.write(fastals[st:en],infile.strip('.fasta') + '_' + str(fileno) + '.fasta','fasta')
    fileno +=1
    st += splitno
    en += splitno



