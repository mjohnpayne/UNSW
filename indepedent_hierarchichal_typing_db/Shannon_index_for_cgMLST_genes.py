from time import sleep as sl

import numpy as np
import pandas as pd
from csv import DictReader
import csv

#
# for each gene calculate number of alleles, number of occurences for each allele


# def shannon_index_gene(genelist):




def makedict(infile):
    columns = []
    with open(infile, 'rU') as f:
        reader = csv.reader(f,delimiter="\t")
        for row in reader:
            if columns:
                for i, value in enumerate(row):
                    columns[i].append(value)
            else:
                # first row
                columns = [[value] for value in row]
    # you now have a column-major 2D array of your file.
    as_dict = {c[0]: c[1:] for c in columns}
    return(as_dict)

def get_shannon_simp(col):
    s = len(col)
    dict = {}
    tot = 0
    simtot = 0
    for i in col:
        if i not in dict:
            dict[i] = 1
        else:
            dict[i] +=1
    for i in dict:
        stat = float(dict[i])/s
        shan = stat*np.log(stat)
        simtot += float(dict[i])*(float(dict[i])-1)
        tot += shan
    shannon = (-1*tot)
    simpson = 1-(simtot/(s*(s-1)))
    return shannon,simpson,s



def main(infile):
    indict = makedict(infile)
    outfile = open("/Users/michaelpayne/Documents/UNSW/Salmonella/new_typing_scheme/shannon_div_cgMLST_genes.txt","w")
    outfile.write("Gene\tShannon diversity\tShannon evenness\tSimpsons diversity index\n")
    for i in indict:
        shan,simp,s = get_shannon_simp(indict[i])
        shan_even = shan/s
        outfile.write(i+"\t"+str(round(shan,3))+"\t" + str(round(shan_even,10))+"\t" + str(round(simp,3))+'\n')
    outfile.close()





inf = "/Users/michaelpayne/Documents/UNSW/Salmonella/data_aquisition/enterobase_data/cgMLST_definitions/cgMLST-v2-profiles-8-6-17"

main(inf)