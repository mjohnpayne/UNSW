from time import sleep as sl

import glob

import os

change_files = glob.glob("/Users/michaelpayne/Documents/UNSW/Salmonella/Sophie_STM_outbreak/nullarbor_snps/*_snps.tab")

table = open("/Users/michaelpayne/Documents/UNSW/Salmonella/Sophie_STM_outbreak/SraRunTable.txt","r").read().splitlines()

conv = {x.split('\t')[9]:x.split('\t')[4] for x in table}

for i in change_files:
    name = i.split("/")[-1].replace("_snps.tab","")
    newname = "/".join(i.split("/")[:-1]) + "/" + conv[name] + "_snps.tab"
    # os.rename(i,newname)
