from time import sleep as sl


present = []
outfile = open("/Users/michaelpayne/Documents/UNSW/Salmonella/STM_tree_frm_entbase/STM_cgmlst_for_distance_for_eburst_dedup.txt","w")
for i in open("/Users/michaelpayne/Documents/UNSW/Salmonella/STM_tree_frm_entbase/STM_cgmlst_for_distance_for_eburst.txt","r").readlines():
    col = i.split('\t')
    new = '\t'.join(col[1:])
    if new not in present:
        outfile.write(i+'\n')
        present += [new]
outfile.close()