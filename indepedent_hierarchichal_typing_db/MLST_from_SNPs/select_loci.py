


import random

incgmlst_loci = "/Users/michaelpayne/Documents/UNSW/Salmonella/Enterobacteriaceae_core/cgMLST_accessions.txt"

random.seed(1234)

def select_random(lst,num):
    pos = random.sample(xrange(len(lst)), num)
    newlis = []
    for i in pos:
        newlis.append(lst[i])
    for i in newlis:
        lst.remove(i)
    return newlis,lst

def select_genes(cgmlst_loci):
    genes7 = [x for x in open("/Users/michaelpayne/Documents/UNSW/Salmonella/Enterobacteriaceae_core/STMMW_7_gene_accs.txt").read().split('\r')]
    cgmlst_genes = [x for x in open(cgmlst_loci).read().split('\r')]
    genes=cgmlst_genes
    set_dict = {}
    print len(genes)
    for i in genes7:
        genes.remove(i)
    # set_dict["7"] = genes7
    print len(genes)
    genes20,genes = select_random(genes,20)
    set_dict["20"] = genes20
    print len(genes)
    genes50,genes = select_random(genes,50)
    set_dict["50"] = genes50
    print len(genes)
    genes100,genes = select_random(genes,100)
    set_dict["100"] = genes100
    print len(genes)
    genes250,genes = select_random(genes,250)
    set_dict["250"] = genes250
    print len(genes)
    genes500,genes = select_random(genes,500)
    set_dict["500"] = genes500
    print len(genes)
    genes750,genes = select_random(genes,750)
    set_dict["750"] = genes750
    print len(genes)
    genes1000,genes = select_random(genes,1000)
    set_dict["1000"] = genes1000
    print len(cgmlst_genes)
    print len(genes)
    cgmlst = [x for x in open(cgmlst_loci).read().split('\r')]
    return set_dict,cgmlst

sets_dict,cgmlst = select_genes(incgmlst_loci)

## output gene set lists

outpath = "/Users/michaelpayne/Documents/UNSW/Salmonella/new_typing_scheme/cgMLST/new_cgmlst_loci_accession_lists/"

for i in sets_dict:
    outf = open(outpath+i+"_gene_accessions.txt","w")
    outf.write("\n".join(sets_dict[i])+"\n")
    outf.close()

## get scheme genetic length

sets_dict["3002"] = cgmlst

def print_scheme_length(genelists,lenfile):

    genes7 = [x for x in open("/Users/michaelpayne/Documents/UNSW/Salmonella/Enterobacteriaceae_core/STMMW_7_gene_accs.txt").read().split('\r')]
    genelists["7"] = genes7

    pos = open(lenfile,"r").read().strip('\n').split('\n')
    pos = [x.split('\t') for x in pos]
    pos = {x[0]:(int(x[2])-int(x[1])) for x in pos}
    for i in sorted(map(int,genelists.keys())):
        i = str(i)
        print "Scheme " + i
        sum = 0
        for j in genelists[i]:
            sum += pos[j]
        print "Scheme " + i,sum
lens = "/Users/michaelpayne/Documents/UNSW/Salmonella/new_typing_scheme/cgMLST/db_starter/new_cgmlst_allele_locations.txt"

print_scheme_length(sets_dict,lens)