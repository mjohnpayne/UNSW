from time import sleep as sl

infile = open('/Users/michaelpayne/Documents/UNSW/Salmonella/STM_tree_frm_entbase/STM_cgmlst_for_distance.txt','r').readlines()

ls = []

for i in infile:
    col = i.split('\t')
    ls += [int(col[0])]

print sorted(ls)[:10]