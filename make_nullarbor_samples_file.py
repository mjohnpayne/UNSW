from time import sleep as sl

import glob
import sys

inf = sys.argv[1]

filelis = glob.glob(inf + '/*fastq*')
print filelis
outfile = open(sys.argv[2],"w")

for i in filelis:
    print i
    if "_1.fastq" in i:
        acc = i.split('/')[-1].strip("_1.fastq")
        outfile.write(acc +'\t' + i + '\t' + i.strip("_1.fastq")+'_2.fastq\n')
outfile.close()
