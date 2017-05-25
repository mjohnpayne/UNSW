import numpy as np
import csv
from time import sleep as sl

cgmlst_types = "/Users/mjohnpayne/Documents/UNSW/Salmonella/STM_tree_frm_entbase/SALwgMLST.cgMLSTv1.profiles.list.txt"

stm_seqtypes = "/Users/mjohnpayne/Documents/UNSW/Salmonella/STM_tree_frm_entbase/STM_cgmlst_22-5-17.txt"



cgmlst_mat = csv.reader(open(cgmlst_types,"r"),delimiter='\t')

stm_types = csv.reader(open(stm_seqtypes,"r"),delimiter='\t')

types = list(stm_types)
schema = list(cgmlst_mat)

print len(schema[0])

# check what ST are in typhimurium strains

stm_st = []
all = 0
for i in open(stm_seqtypes,"r").readlines()[1:]:
    col = i.strip('\n').split('\t')
    all +=1
    if col[-1] not in stm_st:
        stm_st.append(col[-1])

print all,len(stm_st)

##make profiles file with only STM profiles

# stm_cgmlst_profiles = open("/Users/mjohnpayne/Documents/UNSW/Salmonella/STM_tree_frm_entbase/STM_cgmlst.profiles.list","w")
# stm_cgmlst_profiles = open("/Users/mjohnpayne/Documents/UNSW/Salmonella/STM_tree_frm_entbase/STM_cgmlst_for_distance.txt","w")

count = 0
col = []

for i in open(cgmlst_types,"r").readlines():
    # if count == 0:
    #     continue
        # stm_cgmlst_profiles.write(i)
        # stm_cgmlst_profiles.write("10804\n3003\n")
    # else:
    col = i.strip('\n').split('\t')
    #print col[0],len(col[1:])
    if col[0] in stm_st:
        # stm_cgmlst_profiles.write(i)
        count+=1
print count,col[0],len(col[1:])

# stm_cgmlst_profiles.close()

# remove genes that do not vary in any ST

# stm_cgmlst_profiles = csv.reader(open("/Users/mjohnpayne/Documents/UNSW/Salmonella/STM_tree_frm_entbase/STM_cgmlst.profiles.list","r"),delimiter="\t")
#
# stm_cgmlst_profiles = list(stm_cgmlst_profiles)
#
# print len(stm_cgmlst_profiles),len(stm_cgmlst_profiles[0])
#
# stm_arr = np.array(stm_cgmlst_profiles).astype("str")
#
# test_arr = stm_arr[1:,1:]
# #col1 = set(list(test_arr[:,2]))
#
# #print col1
#
# for i in stm_arr.T[1:]:
#     nums = list(set(i[1:]))
#     if len(nums) < 5:
#         print len(i[1:])
#         counts = []
#         for j in nums:
#             print j,list(i).count(j)
#         sl(0.4)
#
# tfstm_arr = np.all(test_arr == test_arr[0,:], axis = 0)
#
# print np.size(tfstm_arr),np.size(tfstm_arr)-np.count_nonzero(tfstm_arr)