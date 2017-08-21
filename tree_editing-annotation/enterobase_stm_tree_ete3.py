from time import sleep as sl
from ete3 import treeview
from ete3 import Tree, faces, AttrFace, TreeStyle, NodeStyle, PhyloTree, PieChartFace

# t = PhyloTree(newick="/Users/michaelpayne/Documents/UNSW/Salmonella/STM_tree_frm_entbase/R2/STM_cgmlst_seqtype_tree_fix2.nwk",format=1)
# t.resolve_polytomy()
#t = PhyloTree(newick="/Users/michaelpayne/Documents/UNSW/Salmonella/new_typing_scheme/cgMLSTv2/cgMLST_v2_frm_snps/NZ_outbreak_only.nwk")



# ancestor=t.get_common_ancestor("42208","39901")
#
# print ancestor.get_ascii(attributes=["name"], show_internal=False)
#
# ancestor.write(format=1,outfile="/Users/michaelpayne/Documents/UNSW/Salmonella/STM_tree_frm_entbase/R2/STM_NSW_clade.nwk")

import colorsys

##gets N colours equally spaces as hex strings

def get_N_HexCol(N=5):

    HSV_tuples = [(x*1.0/N, 0.8, 0.8) for x in xrange(N)]
    hex_out = []
    for rgb in HSV_tuples:
        rgb = map(lambda x: int(x*255),colorsys.hsv_to_rgb(*rgb))
        hex_out.append("".join(map(lambda x: chr(x).encode('hex'),rgb)))

    HSV_tuples = [(x*1.0/N, 0.1, 0.9) for x in xrange(N)]
    hex2_out = []
    for rgb in HSV_tuples:
        rgb = map(lambda x: int(x*255),colorsys.hsv_to_rgb(*rgb))
        hex2_out.append("".join(map(lambda x: chr(x).encode('hex'),rgb)))

    return hex_out, hex2_out

## makes dict of strain and hier STs in shemels order
# in_hier = open("/Users/michaelpayne/Documents/UNSW/Salmonella/new_typing_scheme/cgMLSTv2/cgmlstv2_from_enterobase/hier_for_tree_anno_hcorrected0.txt","r").read().splitlines()
in_hier = open("/Users/michaelpayne/Documents/UNSW/Salmonella/new_typing_scheme/cgMLSTv2_stm/stm_scheme/hierMLST_stm_hierarchical_assignments_hcorrected0.txt","r").read().splitlines()

in_cc = open("/Users/michaelpayne/Documents/UNSW/Salmonella/new_typing_scheme/cgMLSTv2_stm/stm_scheme/hierMLST_stm_cc_hierarchy.txt","r").read().splitlines()

schemes = ['7','20','30','100','160','260','420','560','1445', "cgMLST","stmcgMLST"]

annots = {}

# t = Tree(newick="/Users/michaelpayne/Documents/UNSW/Salmonella/new_typing_scheme/cgMLSTv2_stm/stm_scheme/Trees/all_inc_soph.nwk")
t=Tree(newick="/Users/michaelpayne/Documents/UNSW/Salmonella/new_typing_scheme/cgMLSTv2_stm/stm_scheme/Trees/Westmead_strains/subclade_westmead.nwk")
stlist = []
for node in t.traverse():
    stlist.append(node.name)


cols = {}



# for i in in_hier[1:]:
#     i = i.split('\t')
#     sts = i[1].split("-")
#     annots[i[0]] = sts

####for corr hierarchical


# for i in in_hier[1:]:
#     i = i.split('\t')
#     if i[0] in stlist:
#         sts = i[1].split("-")
#         annots[i[0]] = sts

#### for clonal complex

for i in in_cc[1:]:
    i = i.split('\t')
    if i[0] in stlist:
        sts = i[-1].split("-")
        annots[i[0]] = sts

sttypes = {}

print "done1"
for i in range(len(schemes)):
    sttypes[schemes[i]] = []
    c=0
    for j in annots:
        if annots[j][i] not in sttypes[schemes[i]]:
            sttypes[schemes[i]].append(annots[j][i])
            c+=1
    print schemes[i],c

##
print "done2"

cold = {}
for i in range(len(schemes)):
    cold[schemes[i]] = {}
    ls = []
    for j in annots:
        ls.append(annots[j][i])
    st = list(set(ls))
    colist, colist2 = get_N_HexCol(len(st))
    for j in range(len(st)):
        cold[schemes[i]][st[j]] = ["#"+colist[j],"#"+colist2[j]]

print "done3"

# t = PhyloTree(newick="/Users/michaelpayne/Documents/UNSW/Salmonella/STM_tree_frm_entbase/R2/STM_cgmlst_seqtype_tree_fix2.nwk",format=1)
# t.resolve_polytomy()
#t = PhyloTree(newick="/Users/michaelpayne/Documents/UNSW/Salmonella/new_typing_scheme/cgMLSTv2/cgMLST_v2_frm_snps/NZ_outbreak_only.nwk")



# ancestor=t.get_common_ancestor("42208","39901")
#
# print ancestor.get_ascii(attributes=["name"], show_internal=False)
#
# ancestor.write(format=1,outfile="/Users/michaelpayne/Documents/UNSW/Salmonella/STM_tree_frm_entbase/R2/STM_NSW_clade.nwk")
# t = PhyloTree(newick="/Users/michaelpayne/Documents/UNSW/Salmonella/STM_tree_frm_entbase/R2/2_aus_clades.nwk",format=1)


##Add schemes at faces

#usedschemes = ["7","cgMLST"]
usedschemes = ['7','20','30','100','160','260','420','560','1445', "cgMLST","stmcgMLST"]
name = 1
c=0
for node in t.traverse():
    c+=1
    for i in schemes:
        node.__setattr__(i," x")
    if node.is_leaf():
        continue
    else:
        node.name="node_"+str(name)
        name+=1

print c," Nodes"

print "done4"

# for node in t.traverse():
#     print node.name,node.__getattribute__('7')
#     sl(0.3)


for leaf in t.traverse():
    for i in annots:
        if leaf.name == i:
            for j in range(len(schemes)):
                leaf.add_feature(schemes[j]," "+annots[i][j])

print "done4"

# for node in t.traverse():
#     print node.name
#     print node.__getattribute__("7")
#     sl(0.4)
# print t.get_ascii(attributes=["name", '7','20','30','100','160','260','420','560','1445', "cgMLST"], show_internal=False)

c = 0
nodetypes = {}
for i in usedschemes:
    print i,"\n\n\n"
    c = 0
    nodetypes[i] = {}
    for j in sttypes[i]:
        c+=1
        if c % 100 == 0:
            print c," types done"
        for node in t.get_monophyletic(values=[" "+j], target_attr=i):
            # print i, j
            # print node.get_ascii(attributes=["name", '7','20','30','100','160','260','420','560','1445', "cgMLST"], show_internal=True)
            # node.del_feature(i)
            node.add_feature(i," "+j)
            nodetypes[i][node.name] = j
            for k in node.children:
                nodetypes[i][k.name] = j


# print t.get_ascii(attributes=["name", '7','20','30','100','160','260','420','560','1445', "cgMLST"], show_internal=True)

bg_scheme = "7"
for s in usedschemes:
    def my_layout(node):
        nstyle = NodeStyle()
        nstyle["size"] = 1
        if node.is_leaf():
            node.set_style(nstyle)
            nameface = AttrFace("name", fsize=16,)
            faces.add_face_to_node(nameface, node, column=0, aligned=False)
            col=1
            for i in range(len(schemes)):
                if schemes[i] in usedschemes:
                    if node.name in annots:
                        st = annots[node.name][i]
                        stcol = cold[schemes[i]][st][0]
                        nface = AttrFace(schemes[i], fsize=18,fgcolor=stcol,)
                        faces.add_face_to_node(nface, node, column=col, aligned=True)
                        col+=1
        else:
            if node.name in nodetypes[s]:
                nst1 = NodeStyle()
                nst1["bgcolor"] = cold[s][nodetypes[s][node.name]][1]
                nst1["size"] = 0
                node.set_style(nst1)
            else:
                nst1 = NodeStyle()
                nst1["size"] = 0
                node.set_style(nst1)
            # col=0
            # for i in range(len(schemes)):
            #     if schemes[i] in usedschemes:
            #         # st = node.get_cached_content()
            #         # stcol = cold[schemes[i]][st]
            #         nface = AttrFace(schemes[i], fsize=12)#, fgcolor=stcol)
            #         faces.add_face_to_node(nface, node, column=col, position="branch-right")
            #         col+=1


    ts = TreeStyle()
    # Do not add leaf names automatically ts.show_leaf_name = False
    # Use my custom layout
    ts.show_leaf_name = False
    ts.layout_fn = my_layout
    ts.tree_width=1500
    # ts.mode = 'c'


    # t.show(tree_style=ts)

    t.render(file_name="/Users/michaelpayne/Documents/UNSW/Salmonella/new_typing_scheme/cgMLSTv2_stm/stm_scheme/Trees/Westmead_strains/westmead_sub_cc3_"+s+".pdf",tree_style=ts,dpi=200)


