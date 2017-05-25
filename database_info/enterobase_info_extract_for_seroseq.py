
import sys
import random

from time import sleep as sl



infile = open(sys.argv[1],"r").readlines()

outfile = open(sys.argv[1][:-4] + "_seqsero_10only_accessions.txt","w")

outfile2 = open(sys.argv[1][:-4] + "_seqsero_10only_details.txt","w")

serovar_list = ["Aberdeen","Chester","Infantis","Muenchen","Saintpaul","Stanley","Virchow","Waycross","Typhimurium","Enteritidis"]

strains = {}
count = {}
countycount = {}
deets = {}

# lines into list
# count of number of strains in each country

random.seed(234)

for i in infile[1:]:
    col = i.strip('\n').split("\t")
    if col[19] != "0" and col[19] != "" and col[2].split(";")[0][:3] in ("SRR", "ERR", "DRR"):
        if col[19] not in strains:
            strains[col[19]] = [col[2].split(";")[0]]
            deets[col[19]] = ['\t'.join(col) + '\n']
        else:
            strains[col[19]].append(col[2].split(";")[0])
            deets[col[19]] += ['\t'.join(col) + '\n']

##dict(country(listofSRRs))

for i in strains:
    if i in serovar_list:
        count = 0
        if len(strains[i]) < 5:
            outfile.write("\n".join(strains[i]) +"\n")
            outfile2.write(deets[i])
        else:
            numlis = random.sample(range(0,len(strains[i])),5)
            for j in numlis:
                outfile.write(strains[i][j]+"\n")
                outfile2.write(deets[i][j])


outfile.close()
outfile2.close()