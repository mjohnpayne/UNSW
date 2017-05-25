__author__ = 'mjohnpayne'

from Bio import SeqIO
import sys

infasta = SeqIO.parse('/Volumes/MP_HD/Pm_Ts_Tf_Pf_comparison/Tf_databases/TF_vel_scaff_prot.fasta','fasta')

newfasta = []
for i in infasta:
    add = 6 - len(i.id)
    newid = 'TFLA_' + str(add*'0') + i.id[1:] + '0'
    print i.id
    i.id = newid
    i.name = ''
    i.description = ''
    print i
    newfasta.append(i)

outfasta = SeqIO.write(newfasta,'/Volumes/MP_HD/Pm_Ts_Tf_Pf_comparison/Tf_databases/TF_vel_scaff_prot_rename.fasta','fasta')