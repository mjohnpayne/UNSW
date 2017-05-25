__author__ = 'mjohnpayne'


from Bio import SeqIO
from Bio import Phylo
from Bio.Seq import Seq
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
from Bio.Blast.Applications import NcbiblastpCommandline as blastp
import sys
import os
import subprocess
from time import sleep as sl
import glob
import math
import shutil
from ete2 import Tree, faces, AttrFace, TreeStyle, NodeStyle, PhyloTree, PieChartFace

indb = "/Volumes/MP_HD/pm1_2161_CI/All_strains_blast_dbs/CI_talaro_db_dedup.fasta"
in_id = "/Users/mjohnpayne/Documents/PhD/pop paper/tm_PrtT_id.txt"#"/Volumes/MP_HD/pm1_2161_CI/pm1_MADS_genes.txt"
in_eval = "1e-5"

align_mao = "/Volumes/MP_HD/pm1_2161_CI/top_hits_pm1_madss/clustal_align_protein.mao"
tree_mao = "/Volumes/MP_HD/pm1_2161_CI/top_hits_pm1_madss/infer_NJ_protein.mao"

def blast_gene(ids,eval,database,of):
    fasta_sequences = SeqIO.parse(open(database),"fasta")
    for seq in fasta_sequences:
        if seq.id == ids:
            SeqIO.write(seq,"temp.fasta", "fasta")
    run = blastp(query='temp.fasta',db=database,num_threads=6,outfmt=5,word_size=4,evalue=eval,out='temp.xml')
    run()
    result_handle = open('temp.xml')
    result = NCBIXML.read(result_handle)
    rets = []
    for i in result.descriptions:
        ttl = i.title
        e = i.e
        species = ttl.split(' ')[0]
        rets.append(species)
        rets.append(str(e))
    # for i in result.alignments:
    #     for j in i.hsps:
    #         rets.append(str(j.frame[1]))
    #         rets.append(str(j.query))
    #         rets.append(str(j.match))
    #         rets.append(str(j.sbjct_start))
    os.remove('temp.fasta')
    os.remove('temp.xml')
    genlis = ["PMAA_018770"]
    for i in range(0,len(rets),2):
        genlis.append(rets[i])
        print rets[i]
    fasta_sequences = SeqIO.parse(open(database),"fasta")
    fs_len = 0
    seqs = []
    for seq in fasta_sequences:
        if seq.id in genlis:
            if seq not in seqs:
                seqs.append(seq)
                fs_len += 1
    SeqIO.write(seqs,of, "fasta")
    return fs_len

align_dir = "/".join(in_id.split('/')[:-1])+"/alignments"
tree_dir = "/".join(in_id.split('/')[:-1])+"/tree"
if os.path.exists(align_dir):
    shutil.rmtree(align_dir)
if os.path.exists(tree_dir):
    shutil.rmtree(tree_dir)


os.mkdir(align_dir)
os.mkdir(tree_dir)

ts = TreeStyle()
# ts.mode = "c"

for i in open(in_id,"r").readlines():
    i=i.strip('\n')
    print i
    outf = "/".join(in_id.split('/')[:-1])+"/top_hits_pm1_madss/"+i+"_blastp_hits_"+in_eval+".fasta"
    no_hits = blast_gene(i,in_eval,indb,outf)
    print no_hits
    align_args = "/usr/local/bin/megacc -a "+ align_mao +" -o "+align_dir+" -s -d " + outf
    subprocess.Popen(align_args, shell=True).wait()
    sl(2)
    align_lis = glob.glob(align_dir + "/*.meg")
    alignpath = ''
    for j in align_lis:
        if i in j:
            tree_args = "/usr/local/bin/megacc -a "+ tree_mao +" -o "+tree_dir+" -d " + j
            subprocess.Popen(tree_args, shell=True).wait()
    tree_ls = glob.glob(tree_dir + "/*.nwk")
    for j in tree_ls:
        if i in j and "consensus" not in j:
            t = PhyloTree(j, format=1)
            #t.show()
            # t = Phylo.read(j,"newick")
            # #t.ladderize()
            # #Phylo.draw(t)
            # Phylo.write(t,j.replace(".nwk",".xml"),"phyloxml")
            # Phylo.draw_graphviz(t,prog="neato")
            t.render(tree_dir+"/"+i+"_blastp_hits_"+in_eval+".pdf",tree_style=ts,dpi=200)
