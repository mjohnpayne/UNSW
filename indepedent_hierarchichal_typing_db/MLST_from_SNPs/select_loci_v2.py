from time import sleep as sl

##v2 includes vivians typable genes and variability sorting - only from 3002 enterobase core genes

import glob
import numpy as np
import operator

'''steps

1st get cgmlst and 84% overlap genes
2nd divide up by presence in core gene sets. i.e. 99 then 97 then 95 etcetc
3rd sort by decreasing diversity - simpsons index

now have ordered list

divide into a number of schemes

use schemes to estimate curve of gene number to ST number
(could also calculate gene number to time discrimination average)


use curve to calculate number of genes to give equal increases in ST number between schemes
(could also select gene number to give various time difference sensitivities)


'''

import random

incgmlst_loci = open("/Users/michaelpayne/Documents/UNSW/Salmonella/Enterobacteriaceae_core/cgMLST_accessions.txt","r").read().splitlines()
cgmlst = open("/Users/michaelpayne/Documents/UNSW/Salmonella/Enterobacteriaceae_core/cgMLST_accessions.txt","r").read().splitlines()
typeable_loci = open("/Users/michaelpayne/Documents/UNSW/Salmonella/new_typing_scheme/vivian_typable_genes.txt","r").read().splitlines()

stm_core_genes = open("/Users/michaelpayne/Documents/UNSW/Salmonella/new_typing_scheme/cgMLST_STM/stm_core_info/core_gene_list_song-total_lst.txt","r").read().splitlines()[1:]



# variability = open("/Users/michaelpayne/Documents/UNSW/Salmonella/new_typing_scheme/diversity estimates/shannon_div_cgMLST_genes.txt","r").read().splitlines()[1:]

gene_con = open("/Users/michaelpayne/Documents/UNSW/Salmonella/LT2_genome/LT2_cgMLST_gene_conversions.txt","r").read().splitlines()
gene_conv = {x.split('\t')[0]:x.split('\t')[1] for x in gene_con}
gene_convr = {x.split('\t')[1]:x.split('\t')[0] for x in gene_con}

##remove 7 MLST genes from list

genes7 = [x for x in open("/Users/michaelpayne/Documents/UNSW/Salmonella/Enterobacteriaceae_core/STMMW_7_gene_accs.txt").read().split('\r')]
print len(incgmlst_loci)

for i in genes7:
    incgmlst_loci.remove(i)

print len(incgmlst_loci)

ent_and_typ = []
ent_only = []
for i in incgmlst_loci:
    if gene_convr[i] in stm_core_genes:
        ent_and_typ.append(i)
    else:
        ent_only.append(i)

print len(incgmlst_loci),len(ent_and_typ),len(ent_only)



def sortbyvar(lst):
    invar = open(
        "/Users/michaelpayne/Documents/UNSW/Salmonella/new_typing_scheme/diversity estimates/shannon_div_cgMLST_genes.txt",
        "r").read().splitlines()[1:]
    inv = [(x.split('\t')[0],float(x.split('\t')[3])) for x in invar if x.split('\t')[0] != "ST"]
    inv.sort(key=operator.itemgetter(1),reverse=True)

    invd={x[0]:float(x[1]) for x in inv}

    nlst = [(x,invd[x]) for x in lst]
    nlst.sort(key=operator.itemgetter(1),reverse=True)

    outls = []

    for i in nlst:
        outls.append(i[0])

    return outls


def give_non_overlap(smlst,lrglst):
    outls = set(lrglst).difference(set(smlst))
    return list(outls)

out = sortbyvar(ent_and_typ)


def sort_by_cons(lst):
    # print len(lst)
    cons_files = glob.glob(
        "/Users/michaelpayne/Documents/UNSW/Salmonella/Enterobacteriaceae_core/entbacMLST_loci/cgMLST_loci_to_use_*")
    nlist = []
    for i in cons_files:
        num = int(i[-6:-4])
        nlist.append(num)
    nlist.sort(reverse=True)
    totlst = []
    outlst = []
    for j in map(str,nlist):
        # print j
        ifile = open("/Users/michaelpayne/Documents/UNSW/Salmonella/Enterobacteriaceae_core/entbacMLST_loci/cgMLST_loci_to_use_"+j+".txt","r").read().splitlines()
        # print len(ifile)
        sectionlist = give_non_overlap(totlst,ifile)
        sectionsubset = []
        totlst +=sectionsubset
        for i in lst:
            if i in sectionlist:
                sectionsubset.append(i)
                lst.remove(i)
        # print len(lst),len(sectionsubset)
        # print len(sectionsubset)
        sortedsection = sortbyvar(sectionsubset)
        # print len(sortedsection)
        outlst +=sortedsection
        # print len(outlst),"\n\n"

    # print len(lst)
    outlst += sortbyvar(lst)
    # print len(outlst)
    return outlst


orderlist = sort_by_cons(ent_and_typ)
print len(orderlist)

orderlist += sort_by_cons(ent_only)

div = [20,30,100,160,260,420,560,1445]

print len(orderlist)

def newlsts(sizes,lst):
    outd = {}
    for i in sizes:
        pos = i - 1
        outd[str(i)] = lst[:i]
        lst = lst[i:]
        print len(lst)
    return outd


gene_dict = newlsts(div,orderlist)

for i in map(str,sorted(map(int,gene_dict.keys()))):
    print i, gene_dict[i]
gene_dict["3002"] = cgmlst

genes7 = [x for x in open("/Users/michaelpayne/Documents/UNSW/Salmonella/Enterobacteriaceae_core/STMMW_7_gene_accs.txt").read().split('\r')]
gene_dict["7"] = genes7

## output gene set lists

outpath = "/Users/michaelpayne/Documents/UNSW/Salmonella/new_typing_scheme/cgMLSTv2/cgmlstv2_loci_accession_lists/"

for i in gene_dict:
    outf = open(outpath+i+"_gene_accessions.txt","w")
    outf.write("\n".join(gene_dict[i])+"\n")
    outf.close()

## get scheme genetic length



def print_scheme_length(genelists,lenfile):

    pos = open(lenfile,"r").read().strip('\n').split('\n')
    pos = [x.split('\t') for x in pos]
    pos = {x[0]:(int(x[2])-int(x[1])) for x in pos}
    for i in sorted(map(int,genelists.keys())):
        i = str(i)
        # print "Scheme " + i
        sum = 0
        for j in genelists[i]:
            sum += pos[j]
        print "Scheme " + i,sum
lens = "/Users/michaelpayne/Documents/UNSW/Salmonella/new_typing_scheme/cgMLST/db_starter/new_cgmlst_allele_locations.txt"

print_scheme_length(gene_dict,lens)
