
# TODO fix dates

from time import sleep as sl
import datetime

def most_common(lst):
    return max(set(lst), key=lst.count)

## record number of strains, number in Australia % , majority country and %, earliest detection, latest detection, animal niche

in_annots = open("/Users/michaelpayne/Documents/UNSW/Salmonella/Enteritidis/MLST_definitions-annotations/Enteritidis_cgMLST.txt","r").read().split("\n")

in_rMLST = [x.split('\t') for x in open("/Users/michaelpayne/Documents/UNSW/Salmonella/Enteritidis/MLST_definitions-annotations/Enteritidis_rMLST.txt","r").read().split('\n')]
rMLST = {x[0]:x[34] for x in in_rMLST[1:] if x[34] != "NaN"}


in_7gene = [x.split('\t') for x in open("/Users/michaelpayne/Documents/UNSW/Salmonella/Enteritidis/MLST_definitions-annotations/Enteritidis_7gene.txt","r").read().split('\n')]
sevengene = {x[0]:x[34] for x in in_7gene[1:] if x[34] != "NaN"}

in70entbactMLST = [x.split("\t") for x in open("/Users/michaelpayne/Documents/UNSW/Salmonella/Enterobacteriaceae_core/70-perc-blast/cgMLST_to_entbacMLST.txt","r").read().split('\n')]

in70entbactMLST = {x[0]:x[1] for x in in70entbactMLST}

# for i in sevengene:
#     print i,sevengene[i]
#     sl(0.2)

ST_annots = open("/Users/michaelpayne/Documents/UNSW/Salmonella/Enteritidis/mleecomp_tree/Enteriditis_tree_annots_26-6-17.txt","w")

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
        if col[0] in in70entbactMLST:
            ST[seqtype]['entbac70MLST'] = in70entbactMLST[col[0]]
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
        elif col[7] != "" and col[8] == "":
            ST[seqtype]['dates'] = [datetime.date(int(col[7]), 1, 1)]
        else:
            ST[seqtype]['dates'] = []
        if col[24] != "":
            ST[seqtype]['Phage-T'] = [col[24]]
    elif seqtype != "NaN":
        ST[seqtype]['count'] += 1
        if col[0] in rMLST:
            ST[seqtype]['rmlst'] = rMLST[col[0]]
        if col[0] in sevengene:
            ST[seqtype]['7gene'] = sevengene[col[0]]
        if col[0] in in70entbactMLST:
            ST[seqtype]['entbac70MLST'] = in70entbactMLST[col[0]]
        if col[11] != "":
            ST[seqtype]['continent'] += [col[11]]
        ST[seqtype]['country'] += [col[12]]
        ST[seqtype]['source_niche'] += [col[4]]
        ST[seqtype]['institution'] += [col[23]]
        if 'Phage-T' in ST[seqtype]:
            ST[seqtype]['Phage-T'] += [col[24]]
        if col[12] == "Australia":
            ST[seqtype]['Australia'] += 1
        if "dates" in ST[seqtype]:
            if col[7] != "" and col[8] != "":
                ST[seqtype]['dates'] += [datetime.date(int(col[7]),int(col[8]),1)]
            elif col[7] != "" and col[8] == "":
                ST[seqtype]['dates'] += [datetime.date(int(col[7]), 1, 1)]
        else:
            if col[7] != "" and col[8] != "":
                ST[seqtype]['dates'] == [datetime.date(int(col[7]),int(col[8]),1)]
            elif col[7] != "" and col[8] == "":
                ST[seqtype]['dates'] == [datetime.date(int(col[7]), 1, 1)]

print len(ST)

ST_annots.write("Sequence type\tNumber of strains\tIn Australia?\tMajority Country\tMajority Country %\tEarliest detection\tLatest detection\tAnimal Niche\trMLST type\t7gene MLST type\tinstitution\tPhage Type\tAus only\tminyear\tmaxyear\tenbact70\n")

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
    PT = ""
    Aus_only = ""
    eb70MLST = ""
    if '7gene' in ST[i]:
        sevgeneST = ST[i]['7gene']
    if 'Phage-T' in ST[i]:
        PT = most_common(ST[i]['Phage-T'])
    if 'rmlst' in ST[i]:
        rST = ST[i]['rmlst']
    if 'entbac70MLST' in ST[i]:
        eb70MLST = ST[i]["entbac70MLST"]
    if len(set(ST[i]['source_niche'])) > 1:
        source = "mixed"
    elif len(set(ST[i]['source_niche'])) == 0:
        source = ""
    else:
        source = ST[i]['source_niche'][0]
    if len(ST[i]["dates"]) > 0:
        mindate = min(ST[i]["dates"])
        maxdate = max(ST[i]["dates"])
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
    if ST[i]["Australia"] > 0:
        Aus_only = "Australia"
    # if ST[i]['Australia'] > 0:
    #     print i,ST[i]['country'],ST[i]['Australia'],common
    #     sl(0.2)
    minyear = ""
    maxyear = ""
    if mindate != "":
        minyear = mindate.year
    if maxdate != "":
        maxyear = maxdate.year


    ST_annots.write(i+'\t'+str(ST[i]["count"])+'\t'+str(ST[i]["Australia"])+'\t'+common+'\t'+str(common_count_perc)+'\t'+str(mindate)+'\t'+str(maxdate)+'\t'+source+'\t'+ rST + '\t' + sevgeneST +'\t'+ instit +'\t'+ PT +'\t'+ Aus_only + '\t' +str(minyear) + '\t' + str(maxyear) + '\t' + eb70MLST + '\n')
    # ST_annots.write(i + '\t' + str(ST[i]["count"]) + '\t' + common + '\t' + str(common_count_perc) + '\t' + mindate + '\t' + maxdate + '\t' + source + '\n')

ST_annots.close()