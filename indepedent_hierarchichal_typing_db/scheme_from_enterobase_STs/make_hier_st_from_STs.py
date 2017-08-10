from time import sleep as sl

import glob


ins = glob.glob("/Users/michaelpayne/Documents/UNSW/Salmonella/new_typing_scheme/cgMLSTv2/cgmlstv2_from_enterobase/*_MLST.txt")

outd = {}

div = ['7','20','30','100','160','260','420','560','1445']

for i in div:
    for j in ins:
        if "/"+i+"_" in j:
            inf = open(j,"r").read().splitlines()
            inf = [x.split('\t') for x in inf]
            for h in inf[1:]:
                if h[0] not in outd:
                    outd[h[0]] = [h[1]]
                else:
                    outd[h[0]].append(h[1])


outf = open("/Users/michaelpayne/Documents/UNSW/Salmonella/new_typing_scheme/cgMLSTv2/cgmlstv2_from_enterobase/hier_for_tree_annot.txt","w")

outf.write("cgMLST_type\thier_type\tlast4hiertype\t"+"\t".join(div)+"\n")

for i in sorted(outd.keys()):
    outf.write(i + '\t' + "-".join(outd[i]) + "-" + i + "\t" +"-".join(outd[i][-3:]) + "-" + i +"\t" + "\t".join(outd[i]) + '\n')
outf.close()