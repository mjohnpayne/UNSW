
import sys

infile = open(sys.argv[1],"r").readlines()

outfile = open(sys.argv[2],"w")
for i in infile:
    if "Input files" in i:
        f = i.strip('\n').split('\t')[1]
    elif "antigenic profile:" in i:
        profile = i.strip('\n').split('\t')[1]
    elif "Predicted serotype" in i:
        sero = i.strip('\n').split('\t')[1]
        outfile.write("%s\t%s\t%s\n"%(f,profile,sero))
    else:
        continue

outfile.close()
