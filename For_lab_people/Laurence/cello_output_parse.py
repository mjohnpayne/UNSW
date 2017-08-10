from time import sleep as sl
import re
import sys

infile = open(sys.argv[1],"r")

id = ""

outpath = sys.argv[1].strip(".txt")+"_result.txt"

outfile = open(outpath,"w")

outfile.write("gene_id\tLocation\tPosition")

for i in infile:
    line = i.strip('\n')
    if "SeqID" in line:
        col = line.split("|")
        id = col[1]
        outfile.write('\n' + id)
    elif line[0] == "*":
        continue
    elif "*" in line:
        line = re.sub(r"\s+", '-', line)
        col = line.split("-")
        loc = col[1]
        sig = col[2]
        outfile.write('\t' +loc + '\t' + sig)

outfile.write("\n")

outfile.close()