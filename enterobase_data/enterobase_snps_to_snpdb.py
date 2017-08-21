from time import sleep as sl


inf = open("/Users/michaelpayne/Documents/UNSW/Salmonella/new_typing_scheme/cgMLSTv2_stm/stm_scheme/locus_relative_mutation_rates/enterobase_derived/all_enterobase_snps.txt","r").read().splitlines()
outf = open("/Users/michaelpayne/Documents/UNSW/Salmonella/new_typing_scheme/cgMLSTv2_stm/stm_scheme/locus_relative_mutation_rates/enterobase_derived/all_enterobase_snpdb.txt","w")
c=0
for i in inf:
    col = i.split('\t')
    c+=1
    splt = col[0].rfind("_")
    name = col[0][:splt]
    allele = col[0][splt+1:]
    outf.write(name + '\t' + allele + '\t' + col[1] + '\n')

outf.close()

