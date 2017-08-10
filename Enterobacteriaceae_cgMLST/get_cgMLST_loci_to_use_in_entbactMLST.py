from time import sleep as sl

from venn_func import venn3

import glob


## make dictionary converting between cgMLST loci and STM gene names and vice versa

blastresults = open("/Users/michaelpayne/Documents/UNSW/Salmonella/LT2_genome/LT2_cgMLST_gene_hits.txt","r")

# outfile = open("/Users/michaelpayne/Documents/UNSW/Salmonella/Enterobacteriaceae_core/cgMLST_to_LT2_STM.txt","w")
# outfile.write("cgMLST_tag\tLT2_gene(STM)\n")

g2cg = {}
cg2g = {}

cgMLST_genes = []

for i in blastresults:
    col = i.strip('\n').split('\t')
    #print col[0],col[1]
    genstm = ""
    for j in col[0].split(" "):
        if "locus_tag" in j:
            genstm = j[j.find("=")+1:-1]
    cgstm = col[1][:col[1].find("_",7)]
    #print genstm,cgstm
    # outfile.write(cgstm + '\t' + genstm + '\n')
    g2cg[genstm] = cgstm
    cg2g[cgstm] = genstm
    cgMLST_genes.append(genstm)
# outfile.close()


## get list of STM genes that are conserved in enterobacteraciae

ent_lst = glob.glob("/Users/michaelpayne/Documents/UNSW/Salmonella/Enterobacteriaceae_core/Roary_20_genome/r2/All_genus*/gene_presence_absence.csv")

for i in ent_lst:
    perc = i.split('/')[-2][-2:]

    ent_cons = open(i,"r")

    cons_list = []

    for i in ent_cons:
        col = i.split(',')
        LT2 = col[29][1:col[29].find(".")]
        # print col[3],col[4]
        if col[3]=='"20"' and col[4] =='"20"':
            cons_list.append(LT2)

    venn3([cgMLST_genes,cons_list],names=["cgMLST","entbact_"+perc+"_core"],show_plot=False,out=True,path="/Users/michaelpayne/Documents/UNSW/Salmonella/Enterobacteriaceae_core/cgMLST_entbac_"+perc+"_venn.pdf")

    cgset = set(cgMLST_genes)
    entbactset = set(cons_list)

    entbact_cg_overlap = cgset.intersection(entbactset)

    entbact_diff = entbactset.difference(cgset)
    cg_diff = cgset.difference(entbactset)

    # print len(cg_diff),len(entbact_cg_overlap),len(entbact_diff)

    stm_to_use = list(entbact_cg_overlap)
    cgMLST_to_use = [g2cg[x] for x in list(entbact_cg_overlap)]

    outfile = open("/Users/michaelpayne/Documents/UNSW/Salmonella/Enterobacteriaceae_core/cgMLST_loci_to_use_"+perc+".txt","w")

    outfile.write("\n".join(cgMLST_to_use))
    outfile.close()

    outfile = open("/Users/michaelpayne/Documents/UNSW/Salmonella/Enterobacteriaceae_core/STM_loci_to_use_"+perc+".txt","w")

    outfile.write("\n".join(stm_to_use))
    outfile.close()