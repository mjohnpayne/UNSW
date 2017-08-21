from time import sleep as sl
import sys
'''

snpdb snp locations

length of each locus

total length of all loci

counts number of snps per locus (in snps/kb) over data

get ratio of mutations per locus and mutations in all loci average (in snps/kb)

'''

def main(prefix):
    locils = open(prefix+"_allele_locations.txt","r").read().splitlines()
    col = [x.split("\t") for x in locils]
    locis = {x[0]:int(x[2])-int(x[1]) for x in col}
    snpdb = open(prefix+"_snpdb.txt","r").read().splitlines()
    col2 = [x.split("\t") for x in snpdb]
    snpno = {}
    snpls = []
    totlen = 0
    for i in locis:
        snpno[i] = []
        totlen += locis[i]
    print "\n\nLoci done\n\n"
    c=0
    for i in col2:
        snps = i[2].split(",")
        if c%10000 == 0:
            print str((float(c)/2381301)*100) + "% done"
        c+=1
        for j in snps:
            if j not in snpno[i[0]]:
                snpno[i[0]].append(j)
                snpls.append(j)
        # snpls.append(i[2].split(","))
    nsnpno = {}
    for i in snpno:
        nsnpno[i] = list(set(snpno[i]))
    snpls=list(set(snpls))
    genome_val = (float(len(snpls))*1000)/totlen

    outfile = open(prefix + "_mutation_rates.txt","w")
    outfile.write("Loci\tlength(bp\tsnp number\tsnps/kb\n")

    outfile.write("Genome\t"+str(totlen)+"\t"+str(len(snpls))+"\t"+str(genome_val)+"\n")
    for i in locis:
        loci_val = float(len(nsnpno[i]))*1000/int(locis[i])
        outfile.write(i+"\t"+str(locis[i])+"\t"+str(len(nsnpno[i]))+"\t"+str(loci_val)+"\t"+str(loci_val-genome_val)+"\n")
    outfile.close()

pref = sys.argv[1]
main(pref)


