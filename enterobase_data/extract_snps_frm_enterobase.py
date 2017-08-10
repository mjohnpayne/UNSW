from time import sleep as sl

import numpy as np
import csv
from time import sleep as sl
from Bio import SeqIO
from Bio import Seq
import glob
from Bio import pairwise2


## ST with allele numbers in cols
cgmlst_types = "/srv/scratch/z3521839/enterobase-API/cgMLST_v2-profiles"

## strain metadata
stm_seqtypes = "/srv/scratch/z3521839/enterobase-API/processing/STM_cgmlst_22-5-17.txt"

## folder with fasta files containing allele seqs, one file per gene
allele_files = "/srv/scratch/z3521839/enterobase-API/cgMLST_v2_alleles-8-6-17/"

## reference used for analysis with SNPS
LT2_genome = "/srv/scratch/z3521839/enterobase-API/processing/LT2_genome.fa"

LT2_ST = "5539"

## file containing SRA identifiers for 20 strains used in comparisons

strain_identifier_ls = "/srv/scratch/z3521839/enterobase-API/processing/strains_used_for_tests.txt"

#### convert metadata and ST schema into list of lists

schema = [x.split('\t') for x in open(cgmlst_types,"r").read().split('\n')]


metadata = [x.split('\t') for x in open(stm_seqtypes,"r").read().split('\r')]

# ## get LT2 ST sequences and generate fasta with all LT2 allele seqs


def extract_allele(schem,ST):
    loci = schem[0]
    for i in schem[1:]:
        if i[0] == ST:
            print schem[0]
            print i
            st_seqs = []
            for j in range(1,len(i)-1):
                gene = loci[j]
                allele = i[j]
                seqs = SeqIO.parse(allele_files+"cgMLST_v2-"+ gene,"fasta")
                for s in seqs:
                    if s.id == gene + "_" + allele:
                        st_seqs.append(s)
            SeqIO.write(st_seqs,"/srv/scratch/z3521839/enterobase-API/processing/LT2_MLST_position_data/LT2_cgMLST_alleles.fa","fasta")
    return

# extract_allele(schema, LT2_ST)

########################

# Run through LT2 alleles and find their position in the LT2 genome
#

def LT2_cgmlst_pos_find(inf,outf):
    Lt2_alleles = SeqIO.parse(
        inf,
        "fasta")
    LT2_positions = open(
        outf,
        "w")

    LT2_gen = SeqIO.parse(LT2_genome, "fasta")
    LT2_gen_s = ""
    for i in LT2_gen:
        LT2_gen_s = i

    pos_dict = {}

    for i in Lt2_alleles:
        st = LT2_gen_s.seq.find(i.seq)
        if st == -1:
            irev = i.reverse_complement()
            st = LT2_gen_s.seq.find(irev.seq)
            en = st + len(i.seq)
            pos_dict[st] = [en, i.id[:i.id.find("_", 7)],"-1"]
        else:
            en = st + len(i.seq)
            pos_dict[st] = [en, i.id[:i.id.find("_", 7)], "1"]

    for i in list(sorted(pos_dict.keys())):
        LT2_positions.write(pos_dict[i][1] + '\t' + str(i + 1) + '\t' + str(pos_dict[i][0]) + '\t' + str(pos_dict[i][2]) +'\n')

    LT2_positions.close()
    return

# LT2_cgmlst_pos_find("/Users/michaelpayne/Documents/UNSW/Salmonella/data_aquisition/enterobase_data/cgMLST_definitions/LT2_MLST_position_data/LT2_cgMLST_alleles.fa","/Users/michaelpayne/Documents/UNSW/Salmonella/data_aquisition/enterobase_data/cgMLST_definitions/LT2_allele_locations.txt")

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
# extract_STs for strains of interest

def extract_sts(metdata,sch,ids):
    idlis = open(ids,"r").readlines()
    idlis = [i.replace('\n',"") for i in idlis]
    id_dict = {'ST':'name'}
    stdict = {}
    stdict[sch[0][0]] = sch[0][1:]
    for j in metdata:
        seq = j[2].split(";")[0]
        if seq in idlis:
            stdict[j[34]] = []
            id_dict[j[34]] = seq
    for j in sch:
        if j[0] in list(stdict.keys()):
            stdict[j[0]] = j[1:]
            print j[0]
    return stdict,id_dict


st_list,ids = extract_sts(metadata,schema,strain_identifier_ls)

##########################
# extract alleles for each test strain

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

def extract_SNPS_from_alleles(LT2_alleles,LT2_positions,allele_folder,outf,st_inf,id_inf):
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
    ind_out = open(
        "/srv/scratch/z3521839/enterobase-API/processing/indel_loci.txt",
        "w")
    ind_out.write("gene\tno. of alleles\tno of indel containing alleles\n")
    alleles_out = []
    alleles_of_interest = "/srv/scratch/z3521839/enterobase-API/processing/alleles_of_interest.fa"

    ## iterate over positions in ST dictionary list to get alleles for each of the STs then go back and collate each ST snp profile
    allele_dict = {}
    # for j in range(1,len(st_inf['ST'])):
    #     gene = st_inf['ST'][j]
    #     allele_list = []
    #     allele_seqs = {i.id:i for i in SeqIO.parse(allele_folder+ "cgMLST_v2-" + gene,"fasta")}
    #     for i in list(sorted(st_inf.keys())):
    #         if i == 'ST':
    #             continue
    #         else:
    #             allele = st_inf[i][j]
    #             if allele not in allele_list:
    #                 allele_list.append(allele)
    #                 if gene + "_" + allele not in allele_seqs:
    #                     print gene,allele
    #                 else:
    #                     alleles_out.append(allele_seqs[gene + "_" + allele])
    #                     allele_dict[gene+"-"+allele] = allele_seqs[gene + "_" + allele]
    #
    # SeqIO.write(alleles_out,alleles_of_interest,"fasta")
    allele_dict = {i.id:i for i in SeqIO.parse(alleles_of_interest,"fasta")}

    snp_dict = {}
    c = 1
    snp_profiles = open("/srv/scratch/z3521839/enterobase-API/processing/snp_profiles_of_interest.txt","w")
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
            snp_dict[i] = snpls
    for i in snp_dict:
        snp_profiles.write(i + '\t' + ','.join(snp_dict[i]) + '\n')
    snp_profiles.close()
    st_snps = {}

    for i in indels:
        ind_out.write(i + '\t' + str(len(indels[i])) + '\t' + str(indels[i].count(1)) + '\n')

    for i in st_inf:
        print i
        st_snps[id_inf[i]] = []
        if i != "ST":
            for j in range(len(st_inf[i])):
                gene = st_inf["ST"][j]
                allele_no = st_inf[i][j]
                key = gene+"_"+allele_no
                if key in snp_dict:
                    st_snps[id_inf[i]] += snp_dict[key]
                else:
                    print key
    for i in st_snps:
        out = open(outf + i + ".txt","w")
        out.write("Position\tReference(LT2)\tStrain\n")
        curlist = st_snps[i]
        ordered_dict = {i[1:-1]:[i[0],i[-1]] for i in curlist}
        for i in map(str,list(sorted(map(int,ordered_dict.keys())))):
            out.write(i + '\t' + ordered_dict[i][0] + '\t' + ordered_dict[i][1] + '\n')
        out.close()




LT2allele = "/srv/scratch/z3521839/enterobase-API/processing/LT2_MLST_position_data/LT2_cgMLST_alleles.fa"
LT2pos = "/srv/scratch/z3521839/enterobase-API/processing/LT2_MLST_position_data/LT2_allele_locations.txt"
ST_alleles = "/srv/scratch/z3521839/enterobase-API/processing/cgMLST_snps_"
extract_SNPS_from_alleles(LT2allele,LT2pos,allele_files,ST_alleles,st_list,ids)

#extract_SNPS_from_alleles(LT2allele,LT2pos,allele_files)

# align to corresponding LT2 allele

# identify missmatches and report position and nature

# adjust position to genome position using LT2 allele location data - take into account rev comp