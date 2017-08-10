from time import sleep as sl
from ete3 import treeview
from ete3 import Tree, faces, AttrFace, TreeStyle, NodeStyle, PhyloTree, PieChartFace


t = PhyloTree(newick="/Users/michaelpayne/Documents/UNSW/Salmonella/new_typing_scheme/cgMLSTv2_stm/stm_scheme/Trees/stm_mlst_3outberak_tree.nwk",format=1)



ancestor=t.get_common_ancestor("anast_45","anast_42")

print ancestor.get_ascii(attributes=["name"], show_internal=False)

ancestor.write(format=1,outfile="/Users/michaelpayne/Documents/UNSW/Salmonella/new_typing_scheme/cgMLSTv2_stm/stm_scheme/Trees/anast_tree.nwk")