
innewick = open("/Users/mjohnpayne/Documents/UNSW/Salmonella/STM_tree_frm_entbase/R2/STM_cgmlst_seqtype_tree.nwk",'r').read()
outnewick = open("/Users/mjohnpayne/Documents/UNSW/Salmonella/STM_tree_frm_entbase/R2/STM_cgmlst_seqtype_tree_fix.nwk",'w')


divd = innewick.split(":-")
newtree = [divd[0]]

for i in divd[1:]:
    new = ["0.00000" + i[7:]]
    newtree += new

outnewick.write(":".join(newtree))

outnewick.close()