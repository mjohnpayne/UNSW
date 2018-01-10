from time import sleep as sl

from Bio import SeqIO
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.SearchIO import BlastIO
from Bio.Blast import NCBIXML
import time
import subprocess as sub
import glob

start_time = time.time()

cline = NcbiblastnCommandline(query="/Users/michaelpayne/Documents/UNSW/Salmonella/new_typing_scheme/allele_scheme_dev/blast_testing/LT2_genome.fasta", db="/Users/michaelpayne/Documents/UNSW/Salmonella/new_typing_scheme/allele_scheme_dev/blast_testing/hierMLST_stm_alleles.fasta",evalue=10,perc_identity=100, out="/Users/michaelpayne/Documents/UNSW/Salmonella/new_typing_scheme/allele_scheme_dev/blast_testing/LT2_blast_all_alleles.xml", outfmt=5,max_target_seqs=1000000,max_hsps=3,word_size=20)
#
print(("--- %s seconds ---" % (time.time() - start_time)))
#
stdout, stderr = cline()
#
# for i in stdout:
#     print(i)
#     sl(0.5)
#
print(("--- %s seconds ---" % (time.time() - start_time)))
#
# r_handle = open("/Users/michaelpayne/Documents/UNSW/Salmonella/new_typing_scheme/allele_scheme_dev/blast_testing/LT2_blast_all_alleles.xml")
#
# blast_records = list(NCBIXML.parse(r_handle))
#
# # for x in blast_records:
# #     print(x
# # c=0
# # for result in blast_records:
# #     print("\n\nresult")
# #     for alignment in result.alignments:
# #         for hsp in alignment.hsps:
# #             name = alignment.title.encode('ascii','ignore')
# #             # if hsp.identities == alignment.length:
# #             if "-1" in name and "--" not in name:
# #                 c+=1
# #                 # print('****Alignment****')
# #                 # print('sequence:', alignment.title)
# #                 # print('length:', alignment.length)
# #                 # print('e value:', hsp.expect)
# #                 # print('identities:', hsp.identities)
# #                 # print('positives:', hsp.positives)
# #                 # print('align_length:', hsp.align_length)
# #                 # print(hsp.query[0:75] + '...')
# #                 # print(hsp.match[0:75] + '...')
# #                 # print(hsp.sbjct[0:75] + '...')
# #                 # sl(0.5)
# #     print(c)
# hitls = []
#
# c=0
# for result in blast_records:
#     print("\n\nresult")
#     for alignment in result.alignments:
#         # name = alignment.title.encode('ascii','ignore')
#         # print(name
#         gename = alignment.title.split(" ")[-1]
#         hitls.append(gename)
#         c2 = 0
#         for hsp in alignment.hsps:
#             if hsp.identities == alignment.length:
#             # if "-1" in name and "--" not in name:
#                 c+=1
#                 c2+=1
#                 if c2 > 1:
#                     print('****Alignment****')
#                     print(('sequence:', alignment.title))
#                     print(('length:', alignment.length))
#                     print('e value:', hsp.expect)
#                     print(('identities:', hsp.identities))
#                     print('positives:', hsp.positives)
#                     print('align_length:', hsp.align_length)
#                     print(hsp.query)
#                     print(hsp.match)
#                     print(hsp.sbjct)
#     print(c)
#
#
# print(("--- %s seconds ---" % (time.time() - start_time)))
#
# allele_db = SeqIO.parse("/Users/michaelpayne/Documents/UNSW/Salmonella/new_typing_scheme/scheme_development/betas/cgMLST_b4/hierMLST_stm_alleles.fasta","fasta")
#
# full_idls = [x.id for x in allele_db]
#
# print((len(full_idls)))
# print((len(hitls)))
# nohits = set(full_idls).difference(set(hitls))
# print((len(nohits)))
# print(nohits)
# nohits = list(nohits)
#
# allele_db = SeqIO.parse("/Users/michaelpayne/Documents/UNSW/Salmonella/new_typing_scheme/scheme_development/betas/cgMLST_b4/hierMLST_stm_alleles.fasta","fasta")
#
# for i in allele_db:
#     if i.id in nohits:
#         print((i.id,len(i.seq)))

# start_time = time.time()
# allelelist = glob.glob("/Users/michaelpayne/Documents/UNSW/Salmonella/new_typing_scheme/scheme_development/betas/cgMLST_b4/hierMLST_stm_alleles/*.fasta")
#
# # make all db separately 45 secs, 800isolates 183 secs
# # make all db separately R2 813 isolates
#
# # together starter = 0.4s, 800isolates 4.2s
#
#
# def makedb(fasta,dbfolder):
#     out = fasta.split("/")[-1]
#     path = dbfolder+"/"+out
#     sub.Popen(['makeblastdb', '-dbtype','nucl','-in',fasta,"-out",path],stdout=sub.PIPE,stderr=sub.PIPE)
#
# dbpath = "/Users/michaelpayne/Documents/UNSW/Salmonella/new_typing_scheme/scheme_development/betas/cgMLST_b4/hierMLST_stm_alleles/blastdbs"
#
#
# for i in allelelist:
#     makedb(i, dbpath)
#     # output, errors = p.communicate()
#     # print(output.decode())
#
# print(("--- %s seconds ---" % (time.time() - start_time)))