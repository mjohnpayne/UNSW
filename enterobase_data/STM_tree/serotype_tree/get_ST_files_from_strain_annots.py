
# TODO fix dates

from time import sleep as sl
import datetime

def most_common(lst):
    return max(set(lst), key=lst.count)

## record number of strains, number in Australia % , majority country and %, earliest detection, latest detection, animal niche

in_annots = open("/Users/michaelpayne/Documents/UNSW/Salmonella/STM_tree_frm_entbase/STM_cgmlst_22-5-17.txt","r").read().split("\n")

in_rMLST = [x.split('\t') for x in open("/Users/michaelpayne/Documents/UNSW/Salmonella/STM_tree_frm_entbase/STM_rMLST_30-5-17.txt","r").read().split('\n')]
rMLST = {x[0]:x[34] for x in in_rMLST[1:] if x[34] != "NaN"}


in_7gene = [x.split('\t') for x in open("/Users/michaelpayne/Documents/UNSW/Salmonella/STM_tree_frm_entbase/STM_7gene_30-5-17.txt","r").read().split('\n')]
sevengene = {x[0]:x[34] for x in in_7gene[1:] if x[34] != "NaN"}

# for i in sevengene:
#     print i,sevengene[i]
#     sl(0.2)

ST_annots = open("/Users/michaelpayne/Documents/UNSW/Salmonella/STM_tree_frm_entbase/STM_cgmlst_22-5-17_ST_annots.txt","w")

ST = {}

print len(in_annots)

for i in in_annots[1:]:
    col = i.split('\t')
    seqtype = col[34]
    if seqtype not in ST and seqtype != "NaN":
        ST[seqtype] = {}
        ST[seqtype]['count'] = 1
        if col[0] in rMLST:
            ST[seqtype]['rmlst'] = rMLST[col[0]]
        if col[0] in sevengene:
            ST[seqtype]['7gene'] = sevengene[col[0]]
        if col[23] != "":
            ST[seqtype]['institution'] = [col[23]]
        else:
            ST[seqtype]['institution'] = []
        if col[11] != "":
            ST[seqtype]['continent'] = [col[11]]
        else:
            ST[seqtype]['continent'] = []
        if col[12] != "":
            ST[seqtype]['country'] = [col[12]]
        else:
            ST[seqtype]['country'] = [""]
        if col[12] == "Australia":
            ST[seqtype]['Australia'] = 1
        else:
            ST[seqtype]['Australia'] = 0
        if col[4] != "":
            ST[seqtype]['source_niche'] = [col[4]]
        else:
            ST[seqtype]['source_niche'] = [""]
        if col[7] != "" and col[8] != "":
            ST[seqtype]['dates'] = [datetime.date(int(col[7]),int(col[8]),1)]
        else:
            ST[seqtype]['dates'] = []
    elif seqtype != "NaN":
        ST[seqtype]['count'] += 1
        if col[0] in rMLST:
            ST[seqtype]['rmlst'] = rMLST[col[0]]
        if col[0] in sevengene:
            ST[seqtype]['7gene'] = sevengene[col[0]]
        if col[11] != "":
            ST[seqtype]['continent'] += [col[11]]
        ST[seqtype]['country'] += [col[12]]
        ST[seqtype]['source_niche'] += [col[4]]
        ST[seqtype]['institution'] += [col[23]]
        if col[12] == "Australia":
            ST[seqtype]['Australia'] += 1
        if col[7] != "" and col[8] != "":
            if "dates" in ST[seqtype]:
                ST[seqtype]['dates'] += [datetime.date(int(col[7]),int(col[8]),1)]
            else:
                ST[seqtype]['dates'] = [datetime.date(int(col[7]), int(col[8]), 1)]

ST_annots.write("Sequence type\tNumber of strains\tIn Australia?\tMajority Country\tMajority Country %\tEarliest detection\tLatest detection\tAnimal Niche\trMLST type\t7gene MLST type\tinstitution\n")
# ST_annots.write("Sequence type\tNumber of strains\tMajority Country\tMajority Country %\tEarliest detection\tLatest detection\tAnimal Niche\n")
# for i in list(sorted(map(int,ST.keys()))):
#     i = str(i)
#     if ST[i]['Australia'] > 0:
#         print i,ST[i]['country'],ST[i]['Australia'],most_common(ST[i]["country"])
#         sl(0.2)

for i in list(sorted(map(int,ST.keys()))):
    i = str(i)
    common = ""
    source = ""
    common_count_perc = 0.0
    sevgeneST = ""
    rST = ""
    instit = ""
    if '7gene' in ST[i]:
        sevgeneST = ST[i]['7gene']
    if 'rmlst' in ST[i]:
        rST = ST[i]['rmlst']
    if len(set(ST[i]['source_niche'])) > 1:
        source = "mixed"
    elif len(set(ST[i]['source_niche'])) == 0:
        source = ""
    else:
        source = ST[i]['source_niche'][0]
    if len(ST[i]["dates"]) > 0:
        mindate = str(min(ST[i]["dates"]))
        maxdate = str(max(ST[i]["dates"]))
    else:
        mindate = ""
        maxdate = ""
    if len(set(ST[i]['country'])) > 1:
        common = most_common(ST[i]["country"])
        common_count_perc = (float(ST[i]["country"].count(common))/float(ST[i]["count"]))*100
    elif len(set(ST[i]['country'])) == 1:
        common = ST[i]["country"][0]
        common_count_perc = 100
    else:
        common = ""
        common_count_perc = ""
    if len(set(ST[i]['institution'])) > 1:
        instit = "mixed"
    elif len(set(ST[i]['institution'])) == 0:
        instit = ""
    else:
        instit = ST[i]['institution'][0]
    # if ST[i]['Australia'] > 0:
    #     print i,ST[i]['country'],ST[i]['Australia'],common
    #     sl(0.2)
    ST_annots.write(i+'\t'+str(ST[i]["count"])+'\t'+str(ST[i]["Australia"])+'\t'+common+'\t'+str(common_count_perc)+'\t'+mindate+'\t'+maxdate+'\t'+source+'\t'+ rST + '\t' + sevgeneST +'\t'+ instit +'\n')
    # ST_annots.write(i + '\t' + str(ST[i]["count"]) + '\t' + common + '\t' + str(common_count_perc) + '\t' + mindate + '\t' + maxdate + '\t' + source + '\n')

ST_annots.close()