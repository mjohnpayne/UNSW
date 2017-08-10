from time import sleep as sl

import sys
import os

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
# prefix = "/Users/michaelpayne/Documents/UNSW/Salmonella/new_typing_scheme/cgMLST_STM/db_starter/stm_int_cgmlst"#sys.argv[1]
# prefix = "/Users/michaelpayne/Documents/UNSW/Salmonella/new_typing_scheme/cgMLST_stm_genes/db_starter/cgMLST_stm_genes"
# prefix = "/Users/michaelpayne/Documents/UNSW/Salmonella/new_typing_scheme/cgMLST_stm_int/db_starter/cgMLST_stm_int"
# prefix="/Users/michaelpayne/Documents/UNSW/Salmonella/new_typing_scheme/cgMLSTv2/cgMLST_v2_frm_snps/dbstarter/hierMLST"
prefix="/Users/michaelpayne/Documents/UNSW/Salmonella/new_typing_scheme/cgMLSTv2_stm/stm_scheme_sep_subtypes/dbstarter/hierMLST_stm"

profiles = prefix + "_profiles.txt"

s_name = prefix.split("/")[-1]

SNPdb = prefix + "_snpdb.txt"

allele_location = prefix + "_allele_locations.txt"

# insnps = sys.argv[2]

# strain = "_".join(insnps.split("/")[-1].split("_")[:-1])

# gen = SeqIO.parse(sys.argv[3],"fasta")

gen = SeqIO.parse('/Users/michaelpayne/Documents/UNSW/Salmonella/new_typing_scheme/LT2_info/LT2_genome.fa',"fasta")
genome = ""
for i in gen:
    genome = str(i.seq)

alleles = prefix + "_alleles.fasta"

allele_folder = prefix + "_alleles/"

annots = prefix + "_assignments.txt"

loci_lists = prefix + "_gene_accessions.txt"

''' make 
profiles file 
gene accessions file
snpdb file
assignments file
alleles file
individual allele files in folder 

using allele location file and genome.fasta'''

def make_dbstart(infile,genome):

    prefix = infile.replace("_allele_locations.txt","")

    inf = open(infile,"r").read().splitlines()

    outprof = open(prefix+"_profiles.txt","w")

    outgene = open(prefix+"_gene_accessions.txt","w")

    outsnpdb = open(prefix+"_snpdb.txt","w")

    outasign = open(prefix+"_assignments.txt","w")

    outasign.write("Strain\t"+s_name+"_ST\nLT2\t1\n")

    os.mkdir(prefix + "_alleles")

    inf = [x.split('\t') for x in inf]
    l1 = ["ST"]
    l2 = ["1"]
    slist = []

    for i in inf:
        name = i[0]
        l1.append(name)
        l2.append("1")
        st = int(i[1])
        en = int(i[2])
        s = SeqRecord(Seq(genome[st-1:en-1]),id=name+"-1",description="")
        SeqIO.write(s,prefix+"_alleles/"+name+"-1.fasta","fasta")
        slist.append(s)

    SeqIO.write(slist,prefix+"_alleles.fasta","fasta")


    outsnpdb.write("\t1\t\n".join(l1[1:])+"\t1\t\n")

    outgene.write("\n".join(l1[1:])+"\n")

    outprof.write("\t".join(l1)+"\n"+"\t".join(l2)+'\n')

    outprof.close()
    outgene.close()
    outsnpdb.close()
    outasign.close()

make_dbstart(allele_location,genome)