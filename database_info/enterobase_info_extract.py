
import sys
import random

from time import sleep as sl



infile = open(sys.argv[1],"r").readlines()

outfile = open(sys.argv[1][:-4] + "_limited_accessions.txt","w")
strains = {}
count = {}
countycount = {}

# lines into list
# count of number of strains in each country

random.seed(234)

for i in infile[1:]:
    col = i.strip('\n').split("\t")
    if col[11] != "" and col[12] != "" and col[2].split(";")[0][:3] in ("SRR", "ERR", "DRR"):
        if col[12] not in strains:
            strains[col[12]] = [col[2].split(";")[0]]
        else:
            strains[col[12]].append(col[2].split(";")[0])

##dict(country(listofSRRs))

for i in strains:
    count = 0
    if len(strains[i]) < 5:
        outfile.write("\n".join(strains[i]) +"\n")
    else:
        numlis = random.sample(range(0,len(strains[i])),5)
        for j in numlis:
            outfile.write(strains[i][j]+"\n")


outfile.close()