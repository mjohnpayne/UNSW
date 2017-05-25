from time import sleep as sl
import random
import glob
from Bio import SeqIO

indata = open("/Users/mjohnpayne/Documents/UNSW/Salmonella/data_aquisition/enterobase_data/7_gene/subspecies-type1-7gene-data.txt","r").readlines()

infastas = glob.glob("/Users/mjohnpayne/Documents/UNSW/Salmonella/data_aquisition/enterobase_data/7_gene/7_gene_alleles/*.fasta")

strain_list = ["aberdeen","birkenhead","chester","infantis","muenchen","saintpaul","stanley","virchow","waycross","typhimurium","enteritidis","wangata","hvittingfoss","weltevreden","zanzibar"]

outannot = open("/Users/mjohnpayne/Documents/UNSW/Salmonella/data_aquisition/enterobase_data/7_gene/14_reduced_tree_annotation.txt","w")

def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False

straindata = {}

for i in indata[1:]:
    col = i.strip('\n').split("\t")
    if col[19] != "":
        straindata[col[0]] = [col[11],col[12],col[19],col[36]]+[col[39:]]

strains = {}

for i in indata[1:]:
    col = i.strip('\n').split("\t")
    sv = str.lower(col[19])
    if col[11] != "" and sv != "" and is_number(sv[0]) == False and "predicted" not in sv and col[34] != "NaN":
        if int(col[34]) > 0:
            sv = sv.split(" ")[0]
            if sv not in strains:
                strains[sv] = {}
                strains[sv][col[11]] = [col[0]]
            else:
                if col[11] not in strains[sv]:
                    strains[sv][col[11]] = [col[0]]
                else:
                    strains[sv][col[11]].append(col[0])

##dict(country(listofSRRs))

random.seed(234)


strainsold = strains
count = 0
for i in strains:
    for j in strains[i]:
        if str.lower(i) in strain_list:
            if len(strains[i][j]) > 30:
                count +=30
                numlis = random.sample(range(0, len(strains[i][j])), 30)
                ls = strains[i][j]
                strains[i][j] = []
                for n in numlis:
                    strains[i][j].append(ls[n])
            else:
                count += len(strains[i][j])
        else:
            if len(strains[i][j]) > 5:
                count +=5
                numlis = random.sample(range(0, len(strains[i][j])), 5)
                ls = strains[i][j]
                strains[i][j] = []
                for n in numlis:
                    strains[i][j].append(ls[n])
            else:
                count += len(strains[i][j])

##store serotype alleles as sequences in dict

alleles = {}

for i in infastas:
    gene = i[-10:-6]
    alleles[gene] = {}
    allele_f = SeqIO.parse(i,"fasta")
    for i in allele_f:
        alleles[gene][i.id[5:]] = i.seq

# outannot.write("uberstrain\tserovar\n")
#
# for i in sorted(strains.keys()):
#     for j in strains[i]:
#             for k in strains[i][j]:
#                 outannot.write(k + '\t' + i + '\n')
# outannot.close()

genes = sorted(alleles.keys())
print genes
outannot.write("uberstrain\tserovar-subset\tserovar\tpred-serovar\tcontinent\tcountry\n")

concats = []
for i in sorted(strains.keys()):
    for j in strains[i]:
            for k in strains[i][j]:
                try:
                    sequence = ""
                    for t in range(len(genes)):
                        sequence += alleles[genes[t]][straindata[k][4][t]]
                    if str.lower(straindata[k][2]) in strain_list:
                        outannot.write(k + '\t' + straindata[k][2]+ '\t' + straindata[k][2]+ '\t' + straindata[k][3]+ '\t' + straindata[k][0]+ '\t' + straindata[k][1] + '\n')
                    else:
                        outannot.write(k + '\t' + "" + '\t' + straindata[k][2] + '\t' + straindata[k][0] + '\t' + straindata[k][1] + '\n')
                    concats.append(SeqIO.SeqRecord(sequence,k,description=""))
                except:
                    pass

outannot.close()
newconcats = []
for i in concats:
    if len(i.seq) == 3336:
        newconcats.append(i)

SeqIO.write(newconcats,"/Users/mjohnpayne/Documents/UNSW/Salmonella/data_aquisition/enterobase_data/7_gene/14_concat_MLST_seqs.fasta","fasta")
