from time import sleep as sl
from ete3 import Tree, faces, AttrFace, TreeStyle, NodeStyle

t = Tree("/Users/michaelpayne/Documents/UNSW/Salmonella/Enteritidis/mleecomp_tree/Enteriditis_phylip_tree_nodes.nwk",format=1)

### make subtree from node

# Nd = "5995"
#
# for node in tree.get_descendants():
#     if node.name == "Node_"+Nd:
#         subtree = node
#         print subtree
#         subtree.write(format=1,outfile="/Users/michaelpayne/Documents/UNSW/Salmonella/Enteritidis/mleecomp_tree/Enteriditis_subclade_node_"+Nd+".nwk")



def layout(node):
    if node.is_leaf():
        N = AttrFace("name", fsize=5)
        faces.add_face_to_node(N, node, 0, position="aligned")

# Set dashed blue lines in all leaves
nst1 = NodeStyle()
nst1["bgcolor"] = "LightSteelBlue"
nst2 = NodeStyle()
nst2["bgcolor"] = "Moccasin"
nst3 = NodeStyle()
nst3["bgcolor"] = "DarkSeaGreen"
nst4 = NodeStyle()
nst4["bgcolor"] = "Khaki"

Nd = "5995"

for node in t.get_descendants():
    if node.name == "Node_5995":
        node.set_style(nst1)
    elif node.name == "Node_5743":
        node.set_style(nst2)
        # subtree.write(format=1,outfile="/Users/michaelpayne/Documents/UNSW/Salmonella/Enteritidis/mleecomp_tree/Enteriditis_subclade_node_"+Nd+".nwk")



ts = TreeStyle()
ts.layout_fn = layout
ts.show_leaf_name = False
t.render("/Users/michaelpayne/Documents/UNSW/Salmonella/Enteritidis/mleecomp_tree/ete3_test.png", w=1000, tree_style=ts)
