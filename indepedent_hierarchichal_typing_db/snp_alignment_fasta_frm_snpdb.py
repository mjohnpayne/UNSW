from time import sleep as sl
slp = sl(0.3)

import sys


def merge_two_dicts(x, y):
    """Given two dicts, merge them into a new dict as a shallow copy."""
    z = x.copy()
    z.update(y)
    return z


def main():
    reference = True
    outdict = {}
    prefix = sys.argv[1]
    strain_ls = open(sys.argv[2],"r").read().strip('\n').split('\n')
    print strain_ls
    assignments = open(prefix + "_assignments.txt","r").read().strip('\n').split('\n')
    assignments = {x.split('\t')[0]:x.split('\t')[1] for x in assignments}
    profiles = open(prefix + "_profiles.txt","r").read().strip('\n').split('\n')
    profiles = [x.split('\t') for x in profiles]
    profdict = {}
    names = profiles[0]
    snplsd = {}

    # dict of ST and list of alleles


    for i in profiles[1:]:
        name = ""
        for j in range(len(i)-1):
            if j == 0:
                name = i[j]
                profdict[name] = []
            else:
                if i[j] != "1":
                    profdict[name].append(names[j]+"-"+i[j])

    # for i in profdict:
    #     print i,profdict[i][:10]
    #     sl(0.3)
    snps = open(prefix + "_snpdb.txt","r").read().splitlines()
    snpdict = {}

    # dict of alleles with dict of position and mutation nucleotide
    # also dict of position and reference for all loci with a snp
    for i in snps:
        col = i.split('\t')
        if col[1] != "1":
            d = {}
            for t in col[2].split(","):
                d[int(t[1:-1])] = t[-1]
            snpdict[col[0]+ "-" + col[1]] = d

            for i in col[2].split(","):
                p = int(i[1:-1])
                ref = i[0]
                new = i[-1]
                if p not in snplsd:
                    snplsd[p] = ref
    # for i in snpdict:
    #     print i,snpdict[i]
    #     sl(0.3)

    #dict of strain name to dict of position to snp change

    convdict = {}
    for i in profdict:
        # print i
        # sl(0.3)
        for t in assignments:
            if assignments[t] == i:
                # print i
                for j in profdict[i]:
                    # print t,assignments[t]
                    # print j
                    # print snpdict[j]
                    # sl(0.3)
                    if t not in convdict:
                        convdict[t] = snpdict[j]
                    else:
                        convdict[t] = merge_two_dicts(convdict[t],snpdict[j])
    # for i in convdict:
    #     print i,sorted(convdict[i].keys())[:5]
    #     sl(0.3)
    # for i in sorted(snplsd.keys()):
    #     print i
    #     sl(0.2)
    # out = prefix + "_concat_snps.fasta"
    out = sys.argv[2].replace(".txt","")+"_concat_snps.fasta"
    outf = open(out,"w")
    newout = {}

    # dict of strain with nucleotide string of all loci with snps , if no snp then inserts reference nucleotide

    for i in strain_ls:
        if i in convdict:
            outdict[i] = ""
            newout[i] = ""
            # outf.write(">"+i+"\n")
            for j in sorted(snplsd.keys()):
                if j in convdict[i].keys():
                    outdict[i] += convdict[i][j]
                else:
                    outdict[i] += snplsd[j]

    # if wanted makes string for reference strain

    if reference == True:
        outdict["LT2"] = ""
        newout["LT2"] = ""
        for i in sorted(snplsd.keys()):
            outdict["LT2"] += snplsd[i]


    # remove positions where there are only reference snps (comes about because snp profiles built for all snps in db

    k = outdict.keys()[0]
    print outdict.keys()
    for i in range(len(outdict[k])):
        count = 0
        cur = ""
        for j in outdict:
            if outdict[j][i] != cur:
                count+=1
                cur = outdict[j][i]
        if count > 1:
            for j in outdict:
                newout[j] += outdict[j][i]

    #output concatenated snp fasta file for tree building

    for i in newout:
        outf.write(">" + i +"\n"+ newout[i]+"\n")
    outf.close()


if __name__ == '__main__':
    main()