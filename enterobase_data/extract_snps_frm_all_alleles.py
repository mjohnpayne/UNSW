from time import sleep as sl

import numpy as np
import csv
from time import sleep as sl
from Bio import SeqIO
from Bio import Seq
import glob
from Bio import pairwise2
import sys

################################

def revcomp(base):
    if base == "T" or base == "t":
        return "A"
    if base == "A" or base == "a":
        return "T"
    if base == "C" or base == "c":
        return "G"
    if base == "G" or base == "g":
        return "C"


def align_and_xtract_snps(ref,query):
    pos = {}
    align = pairwise2.align.globalms(ref, query, 2, -1, -10, -.1)
    r = align[0][0]
    q = align[0][1]
    num = 0
    indel = 0
    if len(ref) != len(query):
        indel = 1
    for i in range(len(r)):
        if r[i] == "-":
            continue
        elif q[i] == "-":
            num+=1
        else:
            if r[i] != q[i]:
                pos[num] = [r[i],q[i]]
            num+=1
    return pos,indel

def extract_SNPS_from_alleles(LT2_alleles,LT2_positions,alleles_fasta, outsnps):
    LT2_alleles1 = SeqIO.parse(LT2_alleles,"fasta")
    LT2_allele_no = {"_".join(i.id.split("_")[:-1]):i.id.split("_")[-1] for i in LT2_alleles1}
    # for i in LT2_allele_no:
    #     print i,LT2_allele_no[i]
    #     sl(0.2)
    LT2_alleles1 = SeqIO.parse(LT2_alleles, "fasta")
    LT2_alleles2 = {"_".join(i.id.split("_")[:-1]):i for i in LT2_alleles1}
    # for i in LT2_alleles2:
    #     print i
    #     sl(0.2)

    LT2_positions = open(LT2_positions,"r").readlines()

    LT2_pos = {}

    for i in LT2_positions:
        col = i.strip('\n').split('\t')
        LT2_pos[col[0]] = col[1:]

    # alleles_fastas = glob.glob(allele_folder + "*")

    indels = {}

    indout = outsnps.replace(".txt","_indels.txt")
    ind_out = open(
        indout,
        "w")
    ind_out.write("gene\tno. of alleles\tno of indel containing alleles\n")

    # alleles_of_interest = "/srv/scratch/z3521839/enterobase-API/cgMLST_v2_all_alleles.fasta"
    # alleles_of_interest = "/srv/scratch/z3521839/enterobase-API/alleles_of_interest.fa"
    alleles_of_interest = alleles_fasta
    allele_dict = {i.id:i for i in SeqIO.parse(alleles_of_interest,"fasta")}
    c = 1
    # snp_profiles = open("/srv/scratch/z3521839/enterobase-API/all_allele_snps/all_snp_profiles_test.txt","w")
    snp_profiles = open(outsnps, "w")
    for i in allele_dict:
        gene = "_".join(i.split("_")[:-1])
        # LT2_allelen = LT2_allele_no[gene]
        if gene in LT2_alleles2:
            LT2_allele = LT2_alleles2[gene]
            snps,ind = align_and_xtract_snps(LT2_allele.seq,allele_dict[i].seq)
            snpls = []

            if gene not in indels:
                indels[gene] = [ind]
            else:
                indels[gene] += [ind]
            if LT2_pos[gene][2] == "1":
                mod = int(LT2_pos[gene][0])
                newsnps = {}
                for j in snps:
                    newsnps[j+mod] = snps[j]
                snpls = [newsnps[t][0] + str(t) + newsnps[t][1] for t in newsnps]
            elif LT2_pos[gene][2] == "-1":
                mod = int(LT2_pos[gene][1])
                newsnps = {}
                for j in snps:
                    newsnps[mod - j] = snps[j]
                snpls = [revcomp(newsnps[t][0]) + str(t) + revcomp(newsnps[t][1]) for t in newsnps]
            if c % 100 == 0:
                print c
            c+=1
            snp_profiles.write(i + '\t' + ','.join(snpls) + '\n')
    snp_profiles.close()
    st_snps = {}

    for i in indels:
        ind_out.write(i + '\t' + str(len(indels[i])) + '\t' + str(indels[i].count(1)) + '\n')
    ind_out.close()




LT2allele = sys.argv[3]#"/srv/scratch/z3521839/enterobase-API/processing/LT2_MLST_position_data/LT2_cgMLST_alleles.fa"
LT2pos = sys.argv[4]#"/srv/scratch/z3521839/enterobase-API/processing/LT2_MLST_position_data/LT2_allele_locations.txt"
input_fasta = sys.argv[1]
outfile = sys.argv[2]
extract_SNPS_from_alleles(LT2allele,LT2pos,input_fasta,outfile)
