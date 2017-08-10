from time import sleep as sl

infile = open("/Users/michaelpayne/Documents/UNSW/Salmonella/LT2_genome/LT2_cgMLST_genome_hits.txt","r").read().strip('\n').splitlines()

outf = open("/Users/michaelpayne/Documents/UNSW/Salmonella/LT2_genome/LT2_cgMLST_genome_positions.txt","w")

for i in infile:
    col = i.split("\t")
    ids = "_".join(col[0].split("_")[:-1])
    genstart = col[8]
    genend = col[9]
    if int(col[9]) > int(col[8]):
        orient = "+"
        st = genstart
        en = genend
    else:
        orient = "-"
        st = genend
        en = genstart
    outf.write("{0}\t{1}\t{2}\t{3}\n".format(ids, st, en, orient))
outf.close()
