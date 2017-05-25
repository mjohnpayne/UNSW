from time import sleep as sl

innewick = open("/Users/mjohnpayne/Documents/UNSW/Salmonella/data_aquisition/genometrakr_data/9-5-17_gtrackr_tree.newick",'r').read()

print innewick[:500],"/n"

inannot = open("/Users/mjohnpayne/Documents/UNSW/Salmonella/data_aquisition/genometrakr_data/9-5-17_gtrackr_metadata.tsv",'r').readlines()

strain_list = ["aberdeen","birkenhead","chester","infantis","muenchen","saintpaul","stanley","virchow","waycross","typhimurium","enteritidis"]

#annot = open("/Users/mjohnpayne/Documents/UNSW/Salmonella/data_aquisition/genometrakr_data/aberdeen_annot.txt","w")

# def rename(intree,indict,count):
#     if count == len(indict):
#         return intree
#     else:
#         key = sorted(indict.keys())[count]
#         retree = intree.replace(key,indict[key])
#         count += 1
#         return rename(retree,indict,count)

rename_dict = {}

## replace taxon names with serovar - or append


outtree = open("/Users/mjohnpayne/Documents/UNSW/Salmonella/data_aquisition/genometrakr_data/9-5-17_gtrackr_serovar_tree.newick",'w')

for i in inannot[1:]:
    col =  i.split('\t')
    #print col[37],col[31]
    if col[31] == "":
        rename_dict[col[37]] = "ND"
    else:
        rename_dict[col[37]] =str.lower(col[31])

print len(rename_dict)

for f_key, f_value in rename_dict.items():
    if f_value == "":
        innewick = innewick.replace(f_key, "none")
    else:
        innewick = innewick.replace(f_key, f_value)
    # print f_key

finaltree = innewick
# finaltree = rename(innewick,rename_dict,0)

print finaltree[:500]
outtree.write(finaltree)
outtree.close()


## output annotation file for figtree


# annot.write('Isolate\tSerovar\n')
#
# for i in inannot[4:]:
#     if i[0] != ";":
#         col =  i[2:-2].split("|")
#         if col[5] == "":
#             annot.write(col[0] + '\t' + '' + '\n')
#         else:
#             if "salmonella" in col[5] or 'Salmonella' in col[5]:
#                 sv = str.lower(col[5])
#                 sv = sv.replace("salmonella ","")
#                 annot.write(col[0] + '\t' + sv + '\n')
#             else:
#                 annot.write("\'" + col[0] + '\'\t' + str.lower(col[5]) +'\n')
#     else:
#         break
#
# annot.close()

##output annot only for single name


# for s in strain_list:
#     annot = open("/Users/mjohnpayne/Documents/UNSW/Salmonella/data_aquisition/genometrakr_data/" + s + "_annot.txt","w")
#
#     annot.write('Isolate\tSerovar\n')
#
#     for i in inannot[4:]:
#         if i[0] != ";":
#             col =  i[2:-2].split("|")
#             sv = str.lower(col[5])
#             if s not in sv:
#                 annot.write(col[0] + '\t' + '' + '\n')
#             else:
#                 annot.write(col[0] + '\t' + s + '\n')
#         else:
#             break
#
#     annot.close()