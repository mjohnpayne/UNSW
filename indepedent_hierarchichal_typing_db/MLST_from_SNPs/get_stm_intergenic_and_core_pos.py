from time import sleep as sl
from Bio import SeqIO



genels = open("/Users/michaelpayne/Documents/UNSW/Salmonella/new_typing_scheme/cgMLST_STM/stm_core_info/core_gene_list_song-total_lst.txt","r").read().splitlines()[1:]

stm_genepos = open("/Users/michaelpayne/Documents/UNSW/Salmonella/LT2_genome/LT2_genome.gb.gff","r").read().splitlines()[6:]

intergenic_ls = open("/Users/michaelpayne/Documents/UNSW/Salmonella/new_typing_scheme/cgMLST_STM/stm_core_info/core_gene_list_song-intergenic.txt","r").read().splitlines()[2:]
intergenic_pos = {x.split('\t')[1]:x.split('\t')[0].split("-") for x in intergenic_ls}

gene_con = open("/Users/michaelpayne/Documents/UNSW/Salmonella/LT2_genome/LT2_cgMLST_gene_conversions.txt","r").read().splitlines()
gene_conv = {x.split('\t')[0]:x.split('\t')[1] for x in gene_con}
gene_convr = {x.split('\t')[1]:x.split('\t')[0] for x in gene_con}
cgMLST_genes_positions = open("/Users/michaelpayne/Documents/UNSW/Salmonella/LT2_genome/LT2_cgMLST_genome_positions.txt","r").read().splitlines()
cgMLST_pos = {x.split('\t')[0]:x.split('\t')[1:] for x in cgMLST_genes_positions}



def merge_two_dicts(x, y):
    """Given two dicts, merge them into a new dict as a shallow copy."""
    z = x.copy()
    z.update(y)
    return z

#### get all stm gene positions

stm_genepos = [x.split('\t') for x in stm_genepos]
stmpos = {}
for i in stm_genepos:
    if len(i) > 1:
        if i[2] == "gene":
            name = i[8].split(";")[0][3:]
            stmpos[name] = [i[3],i[4],i[6]]
print len(stmpos)




#### need to get merged list of cgMLST genes as STMs as well as additional STM genes used - as well as their positions - if gene doesn't have cgmlst equivalent need to add it with position info

outlist = {}
for i in stmpos:
    if i in gene_conv.keys():
        outlist[gene_conv[i]] = cgMLST_pos[gene_conv[i]] + [i] + ["core"]
    elif i in genels:
        # if i not in stmpos:
        #     continue#outlist[i] = "Problematic"
        # else:
        outlist[i] = stmpos[i]+[i]+["stmcore"]
    else:
        outlist[i] = stmpos[i] + [i] + ["noncore"]

print len(outlist)
# used_stm = [x[3] for x in outlist.values()]
# print used_stm[:10]

# for i in stmpos:
#     if i not in used_stm:
#         outlist[i] = stmpos[i]+[i]+["noncore"]

#
# for i in outlist:
#     print i,outlist[i]
#     sl(0.2)





LT2_genome = SeqIO.parse("/Users/michaelpayne/Documents/UNSW/Salmonella/LT2_genome/LT2_genome.fasta","fasta")

genome = ""

# for i in cgMLST_pos:
#     print i,cgMLST_pos[i]
#     sl(0.2)

#
# for i in LT2_genome:
#     genome = i.seq
#
# print genome[:1000]
c = 0
c3 = 0
c5 = 0
old = {}


### cycle through intergenic regions and trim to edges of existing genes - if existing gene within intergenic region

#### first find "intergenic" regions that contain genes from core typing

intergenics = []

for i in intergenic_pos:
    test = False
    intergenic = [int(intergenic_pos[i][0]),int(intergenic_pos[i][1])]
    int_st = intergenic[0]
    int_en = intergenic[1]
    for j in outlist:
        gene_st = int(outlist[j][0])
        gene_en = int(outlist[j][1])
        if gene_st > int_st and gene_en < int_en:
            intergenics.append([int_st,gene_st-1])
            intergenics.append([gene_en+1, int_en])
            test = True
    if test == False:
        intergenics.append([int_st,int_en])

print len(intergenic_pos),len(intergenics)

#### second trim intergenic regions start and ends to match cgMLST scheme

newints = []
for i in intergenics:
    int_st = i[0]
    int_en = i[1]
    for j in outlist:
        gene_st = int(outlist[j][0])
        gene_en = int(outlist[j][1])
        if gene_en > int_st and gene_en < int_en:
            int_st = gene_en+1
        elif int_en > gene_st and gene_en > int_en:
            int_en = gene_st-1
    newints.append([int_st,int_en])

print len(newints)

# for i in range(len(intergenics)):
#     print intergenics[i][1]-intergenics[i][0],newints[i][1]-newints[i][0]
#     sl(0.3)

### get adjacent genes for each intergenic region


gends = sorted([int(outlist[x][1]) for x in outlist])
gstarts = sorted([int(outlist[x][0]) for x in outlist])

endsd = {int(outlist[x][1]):x for x in outlist}
startsd = {int(outlist[x][0]):x for x in outlist}

### smallest int_start - gene end = 5' (non-negative)
### smallest gene start - int end = 3' (non-negative)
intergenics = {}

for i in newints:
    st = i[0]
    en = i[1]
    # print st,en
    endsd = {st-int(outlist[x][1]): x for x in outlist}
    lst5 = [st-x for x in gends]
    pos_dists = [x for x in lst5 if x >= 0]
    smallest5=""
    smallest3=""

    if len(pos_dists) >0:
        smallest5 = min(pos_dists)
    # print smallest5, endsd[smallest5]

    startsd = {int(outlist[x][0])-en: x for x in outlist}
    lst3 = [x-en for x in gstarts]
    pos_lists3 = [x for x in lst3 if x >=0]
    if smallest5 != "":
        smallest3 = min(pos_lists3)
    # print smallest3,startsd[smallest3]
    # print '\n'
    if smallest5 == "":
        continue
    else:
        st = st - smallest5 + 1
        en = en + smallest3 - 1
        intergenics[endsd[smallest5] + "..." + startsd[smallest3]] = (st, en,"","","intergenic_core")

print len(intergenics)

intergenics["STMMW_45431...STM0001"] = (1,189,"","","intergenic_core")

# for i in intergenics:
#     print i,intergenics[i]
#     sl(0.3)
#

outdict = {}
for i in outlist:
    outdict[i] = (int(outlist[i][0]), int(outlist[i][1]), outlist[i][2],outlist[i][3],outlist[i][4])


corec = 0
stmc = 0
noncorec = 0

for i in outlist:
    if outdict[i][4] == "core":
        corec +=1
    elif outdict[i][4] == "stmcore":
        stmc +=1
    elif outdict[i][4] == "noncore":
        noncorec +=1

print corec
print stmc
print noncorec

#
# outdict = merge_two_dicts(intergenics,outdict)
#
# outboth = open("/Users/michaelpayne/Documents/UNSW/Salmonella/new_typing_scheme/cgMLST_STM/stm_intergenic_positions.txt","w")
# outinter = open("/Users/michaelpayne/Documents/UNSW/Salmonella/new_typing_scheme/cgMLST_STM/intergenic_positions.txt","w")
# outstm = open("/Users/michaelpayne/Documents/UNSW/Salmonella/new_typing_scheme/cgMLST_STM/stm_positions.txt","w")
#
# rev = {}
# posls = []
#
# for i in outdict:
#     rev[outdict[i][0]] = i
#     posls += [outdict[i][0]]
#
# posls = sorted(posls)
#
# for i in posls:
#     i = rev[i]
#     if outdict[i][4] == "intergenic_core":
#         outboth.write(i + '\t' + str(outdict[i][0]) + '\t' + str(outdict[i][1]) + '\t' + outdict[i][2] + "\n")
#         outinter.write(i + '\t' + str(outdict[i][0]) + '\t' + str(outdict[i][1]) + '\t' + outdict[i][2] + "\n")
#     elif outdict[i][4] == "stmcore":
#         outboth.write(i + '\t' + str(outdict[i][0]) + '\t' + str(outdict[i][1]) + '\t' + outdict[i][2] + "\n")
#         outstm.write(i + '\t' + str(outdict[i][0]) + '\t' + str(outdict[i][1]) + '\t' + outdict[i][2] + "\n")
#
# outboth.close()
# outinter.close()
# outstm.close()
#
# ### su