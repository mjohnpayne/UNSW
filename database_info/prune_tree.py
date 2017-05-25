from ete3 import Tree


def collapsed_leaf(node):
    if len(node2labels[node]) == 1:
        return True
    else:
        return False

tree = Tree("/Users/mjohnpayne/Documents/UNSW/Salmonella/data_aquisition/genometrakr_data/PDG000000002.910.reference_target.tree.newick",format=1)

node2labels = t.get_cached_content(store_attr="name")


for node in tree.get_leaves():