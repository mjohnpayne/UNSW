from time import sleep as sl

in_hier = "/Users/michaelpayne/Documents/UNSW/Salmonella/new_typing_scheme/cgMLST/new_cgmlst_hierarchical_assignments.txt"

#in_new = "/Users/michaelpayne/Documents/UNSW/Salmonella/new_typing_scheme/cgMLST_STM/stm_int_cgmlst_assignments.txt"
in_new="/Users/michaelpayne/Documents/UNSW/Salmonella/new_typing_scheme/cgMLST_stm_genes/cgMLST_stm_genes_assignments.txt"

#out = "/Users/michaelpayne/Documents/UNSW/Salmonella/new_typing_scheme/cgMLST_STM/stm_int_cgmlst_assignments_hier.txt"
out = "/Users/michaelpayne/Documents/UNSW/Salmonella/new_typing_scheme/cgMLST_stm_genes/stm_cgmlst_assignments_hier.txt"

def add_scheme_to_hier(hier,new,out):
    inh = open(hier,"r").read().splitlines()
    inn = open(new,"r").read().splitlines()
    outf = open(out,"w")
    outf.write(inh[0]+'\tstm_ST_string\tstm_ST\n')
    for i in inh[1:]:
        i = i.split('\t')
        for j in inn[1:]:
            j = j.split('\t')
            if i[0] == j[0]:
                newhier = i[-2] +"-"+ j[1]
                outf.write("\t".join(i) + '\t' + newhier+'\t' + j[1] + '\n')
    outf.close()



add_scheme_to_hier(in_hier,in_new,out)