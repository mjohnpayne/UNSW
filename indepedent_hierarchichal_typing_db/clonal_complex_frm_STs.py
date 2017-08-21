#! /Users/michaelpayne/anaconda/bin/python

from time import sleep as sl
import glob
'''
For each ST at highest level pairwise check allele by allele

if limit of CC is reached stop and move to next ST (should save time)

output dict of CC = [list of STs]

start off with limit of 10 to evaluate


TODO

add assignment to strains using _assignments file

set up different cutoff in script
'''

import sys

# ls = [1,4,2,3,5,2]
#
# def trial(ls,cutoff):
#     c=0
#     for i in ls:
#         if c > cutoff:
#             return "exceed"
#         else:
#             c+=1
#     return "finished"
#
#
# print trial(ls,3)


cutoffs = {'7':1,'20':1,'30':1,'100':1,'160':1,'260':2,'420':2,'560':2,'1445':3, "cgMLST":5,"stmcgMLST":8}

schemels = ['7','20','30','100','160','260','420','560','1445', "cgMLST","stmcgMLST"]

def check_within_cutoff(cutoff,s1,s2):
    c = 0
    for i in range(len(s1)):
        if c > cutoff:
            return False
        else:
            if s1[i] == s2[i]:
                continue
            else:
                c+=1
    return c



def main():
    prefix = sys.argv[1]
    profiless = glob.glob(prefix + "_hierarchical_ST/*_gene_profiles.txt")
    overall = {}
    for p in profiless:
        print p
        name = p.split("/")[-1].replace("_gene_profiles.txt","")
        profiles = open(p,"r").read().splitlines()
        cut = cutoffs[name]
        names = profiles[0].split('\t')[1:]
        print len(names)
        profs = {}
        cc = 1
        ccs = {}
        for i in profiles[1:]:
            lst = i.split("\t")
            if len(profs.keys()) > 0:
                t=0
                for i in profs:
                    if check_within_cutoff(cut,profs[i],lst[1:]) == False:
                        continue
                    else:
                        ccs[lst[0]] = ccs[i]
                        t=1
                if t==0:
                    ccs[lst[0]] = cc
                    cc +=1
                profs[lst[0]] = lst[1:]
            else:
                profs[lst[0]] = lst[1:]
                ccs[lst[0]] = cc
        print cc-1
        print len(profiles)
        hassign = open(prefix + "_hierarchical_assignments.txt","r").read().splitlines()
        hassign1 = hassign[0].split('\t')
        indx = hassign1.index(name + "_type")
        strain = {}
        for i in hassign[1:]:
            col = i.split('\t')
            strain[col[0]] = col[indx]
        outfile = open(prefix + "_hierarchical_ST/" + name + "_gene_" + str(cut) + "snp_cc.txt","w")

        outfile.write("Strain_ID\t"+name+"-type\t"+str(cut)+"snpCC\n")
        for i in strain:
            # print i,strain[i],ccs[strain[i]]
            # sl(0.5)
            if i not in overall:
                overall[i] = {name:ccs[strain[i]]}
            else:
                overall[i][name] = ccs[strain[i]]
            outfile.write(str(i) + '\t' + str(strain[i]) + "\t" + str(ccs[strain[i]]) + "\n")
        outfile.close()

    allout = open(prefix + "_cc_hierarchy.txt","w")
    allout.write("Strain\t" + "\t".join(["cc"+x+"["+str(cutoffs[x])+"]" for x in schemels])+"\tCC-hier\n")
    for i in overall:
        allout.write(i)
        lst = []
        for j in schemels:
            allout.write("\t" + str(overall[i][j]))
            lst.append(str(overall[i][j]))
        allout.write("\t"+"-".join(lst)+"\n")
    allout.close()




main()

