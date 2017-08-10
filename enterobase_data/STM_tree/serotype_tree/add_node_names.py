from time import sleep as sl

# from ete3 import Tree
#
#
# tree = Tree("/Users/michaelpayne/Documents/UNSW/Salmonella/STM_tree_frm_entbase/R2/STM_cgmlst_seqtype_tree_fix.nwk")
#
# edge = 0
# for node in tree.traverse():
#    if not node.is_leaf():
#       node.name = "NODE_%d" %edge
#       edge += 1
#
# tree.write(outfile="/Users/michaelpayne/Documents/UNSW/Salmonella/STM_tree_frm_entbase/R2/STM_cgmlst_seqtype_nodenames.nwk",format=1)

infile = open("/Users/michaelpayne/Documents/UNSW/Salmonella/Enteritidis/mleecomp_tree/Enteriditis_phylip_tree.nwk",'r').read()


tree = infile.split(":")
new = []
count = 1
for i in tree:
    if i[-1].isdigit() == False:
        new.append(i+"Node_"+str(count))
        count+=1
    else:
        new.append(i)

outfile = open("/Users/michaelpayne/Documents/UNSW/Salmonella/Enteritidis/mleecomp_tree/Enteriditis_phylip_tree_nodes.nwk","w")

outfile.write(":".join(new))

outfile.close()