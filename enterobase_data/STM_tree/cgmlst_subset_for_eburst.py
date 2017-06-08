from time import sleep as sl
import random

present = []
outfile = open("/Users/michaelpayne/Documents/UNSW/Salmonella/STM_tree_frm_entbase/STM_cgmlst_for_distance_for_eburst_reduced.txt","w")

random.seed(234)

randlist = random.sample(xrange(0,2980),1000)

infile = open("/Users/michaelpayne/Documents/UNSW/Salmonella/STM_tree_frm_entbase/STM_cgmlst_for_distance_for_eburst_dedup.txt","r").readlines()

for i in randlist:
    outfile.write(infile[i]+'\n')

outfile.close()