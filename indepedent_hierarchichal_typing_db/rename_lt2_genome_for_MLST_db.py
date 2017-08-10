from time import sleep as sl

from Bio import SeqIO







outallele = "/Users/michaelpayne/Documents/UNSW/Salmonella/new_typing_scheme/cgMLST/alleles/"

outprofile = open("/Users/michaelpayne/Documents/UNSW/Salmonella/new_typing_scheme/cgMLST/new_cgmlst_profiles.txt","w")


outls = [["ST"],["1"]]
for i in SeqIO.parse("/Users/michaelpayne/Documents/UNSW/Salmonella/new_typing_scheme/MLST_schemes/LT2_cgMLST_alleles.fa","fasta"):
    name = "_".join(i.id.split("_")[:-1])
    i.id = name+"-1"
    i.description=""
    SeqIO.write(i,outallele+name+".tfa","fasta")
    outls[0].append(name)
    outls[1].append("1")


outprofile.write("\t".join(outls[0])+'\n'+"\t".join(outls[1])+"\n")

outprofile.close()





