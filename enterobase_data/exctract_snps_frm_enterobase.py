from time import sleep as sl

import numpy as np
import csv
from time import sleep as sl
from Bio import SeqIO
from Bio import Seq

## ST with allele numbers in cols
cgmlst_types = "/Users/michaelpayne/Documents/UNSW/Salmonella/data_aquisition/enterobase_data/cgMLST_definitions/cgMLST-v2-profiles-8-6-17"

## strain metadata
stm_seqtypes = "/Users/michaelpayne/Documents/UNSW/Salmonella/data_aquisition/enterobase_data/cgMLST_definitions/STM_cgmlst_22-5-17.txt"

## folder with fasta files containing allele seqs, one file per gene
allele_files = "/Users/michaelpayne/Documents/UNSW/Salmonella/data_aquisition/enterobase_data/cgMLST_definitions/cgMLST_v2_alleles-8-6-17/"

## reference used for analysis with SNPS
LT2_genome = "/Users/michaelpayne/Documents/UNSW/Salmonella/data_aquisition/enterobase_data/cgMLST_definitions/LT2_genome.fa"

LT2_ST = "5539"

## file containing SRA identifiers for 20 strains used in comparisons

strain_identifier_ls = "/Users/michaelpayne/Documents/UNSW/Salmonella/Multiple_SNP_calls_testing/strains_used_for_tests.txt"

#### convert metadata and ST schema into list of lists

# schema = [x.split('\t') for x in open(cgmlst_types,"r").readlines()]
# print len(schema)
#
# types = [x.split('\t') for x in open(stm_seqtypes,"r").readlines()]

# ## get LT2 ST sequences and generate fasta with all LT2 allele seqs
#
# loci = schema[0]
#
# ST_defs = {}
#
# for i in schema[1:]:
#     if i[0] == LT2_ST:
#         print schema[0]
#         print i
#         st_seqs = []
#         for j in range(1,len(i)-1):
#             gene = loci[j]
#             allele = i[j]
#             seqs = SeqIO.parse(allele_files+"cgMLST_v2-"+ gene,"fasta")
#             for s in seqs:
#                 if s.id == gene + "_" + allele:
#                     st_seqs.append(s)
#         SeqIO.write(st_seqs,"/Users/michaelpayne/Documents/UNSW/Salmonella/data_aquisition/enterobase_data/cgMLST_definitions/LT2_cgMLST_alleles.fa","fasta")
#

########################

# Run through LT2 alleles and find their position in the LT2 genome
#
# Lt2_alleles = SeqIO.parse(
#     "/Users/michaelpayne/Documents/UNSW/Salmonella/data_aquisition/enterobase_data/cgMLST_definitions/LT2_cgMLST_alleles.fa",
#     "fasta")
#
# LT2_positions = open(
#     "/Users/michaelpayne/Documents/UNSW/Salmonella/data_aquisition/enterobase_data/cgMLST_definitions/LT2_allele_locations.txt",
#     "w")
#
# LT2_gen = SeqIO.parse(LT2_genome, "fasta")
# LT2_gen_s = ""
# for i in LT2_gen:
#     LT2_gen_s = i
#
# pos_dict = {}
#
# for i in Lt2_alleles:
#     st = LT2_gen_s.seq.find(i.seq)
#     if st == -1:
#         irev = i.reverse_complement()
#         st = LT2_gen_s.seq.find(irev.seq)
#     en = st + len(i.seq)
#     pos_dict[st] = [en, i.id[:i.id.find("_", 7)]]
#
# for i in list(sorted(pos_dict.keys())):
#     LT2_positions.write(pos_dict[i][1] + '\t' + str(i + 1) + '\t' + str(pos_dict[i][0]) + '\n')
#
# LT2_positions.close()

################################

# extract alleles for each test strain

# align to corresponding LT2 allele

# identify missmatches and report position and nature

# adjust position to genome position using LT2 allele location data - take into account rev comp