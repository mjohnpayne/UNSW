from time import sleep as sl
import pandas as pd
import sys

from Bio import SeqIO

prefix = sys.argv[1]

profiles = prefix + "_profiles.txt"

SNPdb = prefix + "_snpdb.txt"

allele_location = prefix + "_allele_locations.txt"

insnps = sys.argv[2]

strain = "_".join(insnps.split("/")[-1].split("_")[:-1])

# gen = SeqIO.parse(sys.argv[3],"fasta")

gen = SeqIO.parse('/Users/michaelpayne/Documents/UNSW/Salmonella/new_typing_scheme/LT2_info/LT2_genome.fa',"fasta")

alleles = prefix + "_alleles.fasta"

allele_folder = prefix + "_alleles/"

annots = prefix + "_assignments.txt"

hier_folder = prefix + "_hierarchical_ST/"

hier_assign = prefix + "_hierarchical_assignments.txt"

loci_lists = prefix+"_loci_accession_lists/"

LT2=""
for i in gen:
    LT2=i.seq


def check_strain(strain,strainls):
    inf = open(strainls).read().strip('\n').split('\n')
    ids = [x.split('\t')[0] for x in inf]
    if strain in ids:
        return True
    else:
        return False

def get_allele_pos(allele_location):
    inf = [x.split('\t') for x in open(allele_location,"r").read().split('\n')]
    outd = {x[0]:x[1:] for x in inf}
    return outd


def parse_snps(snps):
    inf = [x.split('\t') for x in open(snps,"r").read().strip('\n').split('\n')]
    outd = {x[1] : [x[3], x[4]] for x in inf if x[2]=="snp" and len(x) > 2 and len(x[3]) ==1}
    return outd

def replace_pos(pos,nuc,seq):
    if len(pos) == 0:
        return seq
    else:
        seq = seq[:pos[0]]+nuc[0]+seq[pos[0]+1:]
        pos = pos[1:]
        nuc = nuc[1:]
        return replace_pos(pos,nuc,seq)

def get_snp_db(indb):
    inf = [x.split('\t') for x in open(indb,"r").read().split('\n')][:-1]
    # sl(0.3)
    outdb = {}
    for i in inf:
        if i[0] not in outdb:
            outdb[i[0]] = {i[2]:i[1]}
        else:
            outdb[i[0]][i[2]] = i[1]
    return outdb

def test_unique(profsdict,newprof):
    uniq = True
    ident = ""
    for i in profsdict:
        # print profsdict[i][:10],newprof[:10]
        # print len(profsdict[i]),len(newprof)
        if profsdict[i] == newprof:
            uniq = False
            ident = i
    return uniq,ident


def determine_or_assign_st(new,st,profile,inaccs):
    profiles = pd.read_csv(prefix + "_profiles.txt", sep="\t", header=0)
    type = ""
    if new:
        inf = open(prefix + "_gene_accessions.txt", "r").read().strip("\n").splitlines()
        scheme_profiles = pd.read_csv(prefix + "_profiles.txt", sep="\t", header=0)
        # print scheme_profiles
        scheme_profiles = [list(scheme_profiles.columns.values)] + scheme_profiles.values.tolist()
        # print scheme_profiles
        newprofile = pd.DataFrame(data=[profile],columns=profiles.columns.values.tolist())
        # print newprofile
        nprofile = newprofile[["ST"]+inf]
        # print nprofile
        newprofilels = map(int,nprofile.values.tolist()[0])
        # print newprofilels
        new2 = True
        for q in scheme_profiles[1:]:
            # if i == "7":
            #     print q,newprofilels
            #     sl(0.4)
            if q[1:] == newprofilels[1:]:
                type = q[0]
                new2=False
        if new2:
            no = len(scheme_profiles)
            type = no
            outf = open(prefix +"_gene_profiles.txt","a+")
            outf.write(str(no)+"\t" + "\t".join(map(str,newprofilels[1:]))+"\n")
            outf.close()
        # else:
            # full_profs = open(hier_assign,"r").read().strip('\n').splitlines()
            # notnewprof = ''
            # for j in full_profs:
            #     j=j.split("\t")
            #     if j[-1] == st:
            #         notnewprof = strain + "\t".join(i[1:])+'\n'
            # write_profs = open(hier_assign, "a+")
            # write_profs.write(notnewprof)
            # write_profs.close()
    print type
    # outstring = open(prefix + "_hierarchical_assignments.txt", "a+")
    # outstring.write(strain)
    # string = ""
    #
    # for i in schemels:
    #     stype = str(outdict)
    #     if len(string) == 0:
    #         outstring.write('\t'+stype+"\t"+stype)
    #         string = stype
    #     else:
    #         string = string + "-" + stype
    #         outstring.write("\t"+string + "\t" + stype)
    #     outf = open(hier_folder + i + "_gene_annotation.txt","a+")
    #     outf.write(str(outdict["cgMLST"]) + "\t" + stype + "\n")
    #     outf.close()
    # outstring.write('\n')
    # outstring.close()




def main(profile,al_loc,snpdb,snp,genome,annot,STM=False):
    if check_strain(strain,annot)==True:
        sys.exit("\n\nA strain with the name " + strain + " is already in the database.\n")
    else:
        locs = get_allele_pos(al_loc)
        snp = parse_snps(snp)
        snpd = get_snp_db(snpdb)
        # for i in map(str,sorted(map(int,snp.keys()))):
        #     print i,snp[i]
        #     sl(0.3)
        genome = genome
        strain_prof = {}
        newst = False
        # outf = open("/Users/michaelpayne/Documents/UNSW/Salmonella/new_typing_scheme/cgMLST/cgMLST_snpdb.txt","w")
        for i in locs:
            # outf.write(i + "\t1\t\n")
            if len(locs[i]) >0:
                st = int(locs[i][0])
                end = int(locs[i][1])
                sq = genome[st-1:end]
                plis = []
                nuclis = []
                snpls = []
                for j in snp:
                    if st <= int(j) <= end:
                        plis.append(int(j)-st)
                        snpls.append(snp[j][0]+j+snp[j][1])
                        nuclis.append(snp[j][1])
                if len(plis) > 0:
                    snpls = ",".join(snpls)
                    if snpls in snpd[i]:
                        strain_prof[i] = snpd[i][snpls]
                        # print snpls, snpd[i],snpd[i][snpls]
                        # sl(0.5)
                    else:
                        new_allele = str(len(snpd[i])+1)
                        strain_prof[i] = new_allele
                        otf = open(alleles,"a+")
                        otfa = open(allele_folder+i+".tfa","a+")
                        sq = replace_pos(plis,nuclis,sq)
                        otf.write(">"+str(i)+"-"+str(new_allele)+"\n"+str(sq)+"\n")
                        otfa.write(">" + str(i) + "-" + str(new_allele) + "\n" + str(sq) + "\n")
                        otf.close()
                        otfa.close()
                        newst = True
                        otf2 = open(snpdb,"a+")
                        otf2.write(i+'\t'+new_allele+"\t"+snpls+"\n")
                        otf2.close()
                else:
                    strain_prof[i] = "1"
        # if newst == True:
        profs = open(profile,"r").read().split('\n')[:-1]
        pdict = {x[0]:x[1:] for x in [y.split('\t') for y in profs[1:]]}
        names = profs[0].split('\t')[1:]
        sts = {}
        testls = []
        annotations = open(annot,"a+")
        num = len(profs[1:])
        for i in names:
            testls.append(strain_prof[i])
        val,ids = test_unique(pdict,testls)
        # print val
        if val == True:
            profls = []
            new = str(num+1)
            profs = open(profile,"a")
            profs.write(new)
            profls.append(new)
            annotations.write(strain+'\t'+new+'\n')
            annotations.close()
            for i in names:
                profs.write('\t'+strain_prof[i])
                profls.append(strain_prof[i])
            profs.write('\n')
            profs.close()
            determine_or_assign_st(True,new,profls,loci_lists)
        else:
            annotations.write(strain+'\t'+ids+'\n')
            annotations.close()


    # outf.close()


main(profiles,allele_location,SNPdb,insnps,LT2,annots,STM=True)