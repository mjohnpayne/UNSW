from time import sleep as sl
import glob
import pylab
from matplotlib.patches import Circle, Ellipse
from itertools import chain
from collections import Iterable
# from matplotlib_venn import venn3,venn3_unweighted
from matplotlib import pyplot as plt
from matplotlib.pylab import savefig
import numpy as np
from stackedBarGraph import StackedBarGrapher
from venn_func import venn3
from venn_func import get_labels
import numpy as np

alignment = {'horizontalalignment':'center', 'verticalalignment':'baseline'}


nullarbor_in = "/Users/michaelpayne/Documents/UNSW/Salmonella/Multiple_SNP_calls_testing/nullarbor_random_20_snp_files/*snpscg*"

nullarbor_in = glob.glob(nullarbor_in)

enterobase_in = "/Users/michaelpayne/Documents/UNSW/Salmonella/Multiple_SNP_calls_testing/enterobase_snps/cgMLST_snps_*.txt"

enterobase_in = glob.glob(enterobase_in)

SNPfilt_in = "/Users/michaelpayne/Documents/UNSW/Salmonella/Multiple_SNP_calls_testing/snp_filt_outputs/*.txt"

SNPfilt_in = glob.glob(SNPfilt_in)

lanlab_pipeline_in = "/Users/michaelpayne/Documents/UNSW/Salmonella/Multiple_SNP_calls_testing/LanPipe/final_snps/*_cgMLST_only.txt"

lanlab_in = glob.glob(lanlab_pipeline_in)

snappers = "/Users/michaelpayne/Documents/UNSW/fastq_files/snpdb/*_cgMLST_only.vcf"

snap_in = glob.glob(snappers)


def parse_null(inf):
    lst = []
    for i in inf.readlines()[1:]:
        col = i.split('\t')
        if col[2] == "snp":
            snp = col[3]+col[1]+col[4]
            lst.append(snp)
    return lst

def parse_ent(inf):
    lst = []
    for i in inf.readlines()[1:]:
        col = i.strip('\n').split('\t')
        snp = col[1]+col[0]+col[2]
        lst.append(snp)
    return lst

def parse_snpf(inf):
    lst = []
    for i in inf.readlines():
        if i[0] != "#" and i[0] != "f":
            col = i.split(',')
            if col[5] == "X" and col[0] in ["0","1"]:
                snp = col[3]+col[2]+col[4]
                lst.append(snp)
    return lst

def parse_ll(inf):
    lst = []
    for i in inf.readlines()[1:]:
        col = i.strip('\n').split('\t')
        snp = col[2]+col[1]+col[-2]
        lst.append(snp)
    return lst

def parse_snap(inf):
    lst = []
    for i in inf.readlines():
        if i[0] != "#":
            col = i.strip('\n').split('\t')
            if col[4] != "." and col[6] == "PASS":
                snp = col[3]+col[1]+col[4]
                lst.append(snp)
    return lst

# print "strain\ttotal\tnull_only\tentero_only\tsnpf_only\tnull and snp\t null and ent\t ent and snp\t all3"
print "strain\ttotal\tnull\tentero\tsnpf\tlanpipe\tsnapper\tin all\tnull_only\tentero_only\tsnpf_only\tlanpipe_only\tsnapper_only"
plotdict = {}

def venn4_frm_labels(data=None, names=None, labels=None, show_names=True, show_plot=True, out=True, path=None, **kwds):

    if (data is None) or len(data) != 4:
        raise Exception("length of data should be 4!")
    if (names is None) or (len(names) != 4):
        names = ("set 1", "set 2", "set 3", "set 4")

    labels = labels#get_labels(data, fill=fill)

    # set figure size
    if 'figsize' in kwds and len(kwds['figsize']) == 2:
        # if 'figsize' is in kwds, and it is a list or tuple with length of 2
        figsize = kwds['figsize']
    else: # default figure size
        figsize = (10, 10)

    # set colors for different Circles or ellipses
    if 'colors' in kwds and isinstance(kwds['colors'], Iterable) and len(kwds['colors']) >= 4:
        colors = kwds['colors']
    else:
        colors = ['r', 'g', 'b', 'c']

    # draw ellipse, the coordinates are hard coded in the rest of the function
    fig = pylab.figure(figsize=figsize)   # set figure size
    ax = fig.gca()
    patches = []
    width, height = 170, 110  # width and height of the ellipses
    patches.append(Ellipse((170, 170), width, height, -45, color=colors[0], alpha=0.5))
    patches.append(Ellipse((200, 200), width, height, -45, color=colors[1], alpha=0.5))
    patches.append(Ellipse((200, 200), width, height, -135, color=colors[2], alpha=0.5))
    patches.append(Ellipse((230, 170), width, height, -135, color=colors[3], alpha=0.5))
    for e in patches:
        ax.add_patch(e)
    ax.set_xlim(80, 320); ax.set_ylim(80, 320)
    ax.set_xticks([]); ax.set_yticks([]);
    ax.set_aspect("equal")

    ### draw text
    # 1
    pylab.text(120, 200, labels['1000'], **alignment)
    pylab.text(280, 200, labels['0100'], **alignment)
    pylab.text(155, 250, labels['0010'], **alignment)
    pylab.text(245, 250, labels['0001'], **alignment)
    # 2
    pylab.text(200, 115, labels['1100'], **alignment)
    pylab.text(140, 225, labels['1010'], **alignment)
    pylab.text(145, 155, labels['1001'], **alignment)
    pylab.text(255, 155, labels['0110'], **alignment)
    pylab.text(260, 225, labels['0101'], **alignment)
    pylab.text(200, 240, labels['0011'], **alignment)
    # 3
    pylab.text(235, 205, labels['0111'], **alignment)
    pylab.text(165, 205, labels['1011'], **alignment)
    pylab.text(225, 135, labels['1101'], **alignment)
    pylab.text(175, 135, labels['1110'], **alignment)
    # 4
    pylab.text(200, 175, labels['1111'], **alignment)
    # names of different groups
    if show_names:
        pylab.text(110, 110, names[0], fontsize=16, **alignment)
        pylab.text(290, 110, names[1], fontsize=16, **alignment)
        pylab.text(130, 275, names[2], fontsize=16, **alignment)
        pylab.text(270, 275, names[3], fontsize=16, **alignment)

    leg = ax.legend(names, loc='best', fancybox=True)
    leg.get_frame().set_alpha(0.5)

    if show_plot:
        pylab.show()

    if out:
        pylab.savefig(path)


def find_common(nullin,enterin,snpfilt,lanpipe_in,snapper_in):


    inter = []
    ne = []
    ns = []
    es = []
    n = []
    e = []
    s = []
    strains = []
    vennnums = {}
    sizes = {}

    for i in nullin:
        # plt.figure()
        strain = i.split('/')[-1][:-16]
        loc = "/".join(i.split('/')[:-2]) + "/venns/"
        # print strain
        ent_file = ""
        snpf_file = ""
        lanlab_file = ""
        snapper_file = ""

        for j in enterin:
            if strain in j:
                ent_file = j

        for j in snpfilt:
            if strain in j:
                snpf_file = j

        for k in lanpipe_in:
            if strain in k:
                lanlab_file = k

        for t in snapper_in:
            if strain in t:
                snapper_file = t

        null = open(i,"r")
        ent = open(ent_file,"r")
        snpf = open(snpf_file)
        lnlb = open(lanlab_file,"r")
        snap = open(snapper_file,"r")


        null = set(parse_null(null))
        ent = set(parse_ent(ent))
        snpf = set(parse_snpf(snpf))
        lnlb = set(parse_ll(lnlb))
        snap = set(parse_snap(snap))

        # venplot = venn3([null,ent,snpf],("Nullarbor","Enterobase","SNPFilt"))
        all = set()
        all.update(null,ent,snpf,lnlb,snap)
        size = len(all)
        sizes[strain] = size
        # # all.update(snpf)
        null_only = null.difference(ent,snpf,lnlb,snap)
        ent_only = ent.difference(null,snpf,lnlb,snap)
        snpf_only = snpf.difference(null,ent,lnlb,snap)
        lnlb_only = lnlb.difference(null,ent,snpf,snap)
        snap_only = snap.difference(null,ent,snpf,lnlb)
        # nullent = null.intersection(ent)
        # nullent = nullent.difference(snpf)
        # nullsnp = null.intersection(snpf)
        # nullsnp = nullsnp.difference(ent)
        # entsnp = ent.intersection(snpf)
        # entsnp = entsnp.difference(null)
        intersection = ent.intersection(null).intersection(snpf).intersection(lnlb).intersection(snap)
        # strains.append(strain)
        # inter.append(len(intersection))
        # ne.append(len(nullent))
        # ns.append(len(nullsnp))
        # es.append(len(entsnp))
        # n.append(len(null_only))
        # s.append(len(snpf_only))
        # e.append(len(ent_only))
        # plotdict[strain] = [len(intersection),len(nullent),len(nullsnp),len(entsnp),len(null_only),len(ent_only),len(snpf_only)]
        # print strain+'\t'+str(len(all))+'\t'+str(len(null_only))+'\t'+str(len(ent_only))+'\t'+str(len(snpf_only))+'\t'+str(len(nullent))+'\t'+str(len(nullsnp))+'\t'+str(len(entsnp))+'\t'+str(len(intersection))

        print strain + '\t' + str(len(all)) + '\t' + str(len(null)) + '\t' + str(len(ent)) + '\t' + str(len(snpf)) + '\t' + str(len(lnlb)) + '\t' + str(len(snap)) +'\t' + str(len(intersection)) + '\t' + str(len(null_only)) + '\t' + str(len(ent_only)) + '\t' + str(len(snpf_only)) + '\t' + str(len(lnlb_only)) + '\t' + str(len(snap_only))

        labels = get_labels([list(null),list(ent),list(snpf),list(lnlb)])
        vennnums[strain] = labels

        # venn4_frm_labels([list(null),list(ent),list(snpf),list(lnlb)],names=["nullarbor","enterobase","SNP_filt","Lanpipe"],labels=vennnums[strain], figsize=(12, 12),out=False,path="/Users/michaelpayne/Documents/UNSW/Salmonella/Multiple_SNP_calls_testing/4venns/" + strain + "_4way_venn.pdf",show_plot=True)

        venn3([list(null),list(ent),list(snap)],names=["nullarbor","enterobase","SnapperDB"], figsize=(12, 12),out=True,path="/Users/michaelpayne/Documents/UNSW/Salmonella/Multiple_SNP_calls_testing/3venns/" + strain + "_null_ent_snapper_venn.pdf",show_plot=False)

        # sl(0.5)
        # venplot = venn3_unweighted([null, ent, snpf], ("Nullarbor", "Enterobase", "SNPFilt"))
        # plt.title(strain)
        # savefig(loc + strain + "-01-venn.png")


        # shared = set(null).intersection(ent).intersection(snpf)

        # print strain,len(shared),len(null)-len(shared),len(ent)-len(shared)
        # print strain, len(shared), len(null),len(ent),len(snpf)
        # sl(0.2)
    return vennnums,sizes





venndata,datasize = find_common(nullarbor_in,enterobase_in,SNPfilt_in,lanlab_in,snap_in)

# percvenndata = {}
#
#
#
# for i in sorted(venndata.keys()):
#     percvenndata[i] = {x:(float(venndata[i][x])/datasize[i])*100 for x in venndata[i]}
#     # print i, venndata[i]['1111'], percvenndata[i]['1111']
#     sl(0.2)
#
# combined_venn = {}
#
# for i in percvenndata["ERR1040444"]:
#     numlis = []
#     for j in percvenndata:
#         numlis.append(percvenndata[j][i])
#     combined_venn[i] = str("%.2f" % np.mean(numlis))+"%"
#
# venn4_frm_labels([1,2,3,4],
#                  names=["Nullarbor", "Enterobase", "SNP_filt", "Lanpipe"], labels=combined_venn,
#                  figsize=(12, 12), out=True,
#                  path="/Users/michaelpayne/Documents/UNSW/Salmonella/Multiple_SNP_calls_testing/4venns/Overall_percentages_4way_venn.pdf",
#                  show_plot=False)



    ######### Average percentage for one venn