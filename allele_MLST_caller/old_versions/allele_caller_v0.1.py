from time import sleep as sl
import glob
from os.path import basename
from os import remove
import os
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Blast import NCBIXML
import glob
from multiprocessing.dummy import Pool as ThreadPool
import time
import subprocess as sub
import itertools
import shutil
import sys
from Bio import SeqIO
import inspect


pool = ThreadPool(4)

start_time = time.time()

"""
1 - Pick scheme
2 - Blast all loci in scheme
3 - Perfect matches store match and fact of match
4 - any loci with no perfect match store location in query genome from matches to intact alleles
5 - blast newly defined allele region against complete allele set
6 - call alleles
"""



from multiprocessing.dummy import Pool as ThreadPool


def makedb(fasta,dbfolder):
    out = fasta.split("/")[-1]
    path = dbfolder+"/"+out
    p = sub.Popen(['makeblastdb', '-dbtype','nucl','-in',fasta,"-out",path],stdout=sub.PIPE,stderr=sub.PIPE)
    p.communicate()


def mean(numbers):
    return float(sum(numbers)) / max(len(numbers), 1)

def run_blast(query_genome,locus_db,wordsize):
    gene = basename(locus_db).replace(".fasta","")
    tmp_out = "tmp/"+gene+"_tmp_blast.xml"
    cline = NcbiblastnCommandline(
        query=query_genome,
        db=locus_db,
        evalue=10,
        perc_identity=90,
        out=tmp_out,
        outfmt=5,
        max_target_seqs=10000,
        max_hsps=3,
        word_size=wordsize)

    stdout, stderr = cline()

    # r_handle = open(
    #     tmp_out+"gene")

    # blast_records = list(NCBIXML.parse(r_handle))

    # remove(tmp_out)

    # return blast_records

def return_scheme_loci_ls(schemefolder):
    """
    :param schemefolder: accession list folder regex
    :return: dictionary of each scheme as key with list of loci as value
    """
    inschemes = glob.glob(schemefolder)
    schemes = {}

    for scheme in inschemes:
        sc = scheme.split("/")[-1].replace("_gene_accessions.txt", "")
        s = open(scheme, "r").read().splitlines()
        schemes[sc] = s

    return schemes


def get_allele_fastas(allelefolder):
    allele_list = glob.glob(allelefolder)
    allele_fastas = {}
    for allele in allele_list:
        name = basename(allele).replace(".fasta","")
        allele_fastas[name] = allele
    return allele_fastas


def get_word_size(allele_fasta):
    alleles = SeqIO.parse(allele_fasta,"fasta")
    sizes = []
    sizesdict = {}
    for i in alleles:
        sizes.append(len(i.seq))
    avsize = mean(sizes)
    wordsize = int(avsize*0.2)
    # print(wordsize)
    return wordsize


def main(schemelists, alleles, query,scheme):

    if os.path.exists("tmp"):
        shutil.rmtree('tmp')
        os.mkdir("tmp")
    else:
        os.mkdir("tmp")
    inschemes = schemelists + "/*.txt"

    alleles_folder = alleles + "/*.fasta"

    query_genome = query

    scheme_gene_lists = return_scheme_loci_ls(inschemes)

    allele_fastas = get_allele_fastas(alleles_folder)


    target_scheme = scheme

    schemelist = scheme_gene_lists[target_scheme]
    partial_hit = schemelist
    exacthits = 0
    totalhits = 0

    locusls = [allele_fastas[x] for x in schemelist]

    wordsizes = []
    scheme_fastas = []

    for i in schemelist:
        allele_file = allele_fastas[i]
        scheme_fastas.append(allele_file)
        ws = get_word_size(allele_file)
        wordsizes.append(ws)


    tmpdb = "tmp/"+target_scheme+"_alleles.fasta"

    with open(tmpdb, 'w') as outfile:
        for fname in scheme_fastas:
            with open(fname) as infile:
                outfile.write(infile.read())

    makedb(tmpdb, "tmp/")

    allele_sizes = {}

    scheme_alleles = SeqIO.parse(tmpdb,"fasta")


    for allele in scheme_alleles:
        allele_sizes[allele.id] = len(allele.seq)

    #pool.starmap(run_blast, zip(itertools.repeat(query_genome),locusls,wordsizes))

    run_blast(query_genome, tmpdb, 20)

    r_handle = open("tmp/" + target_scheme + "_alleles_tmp_blast.xml")

    blast_hits = list(NCBIXML.parse(r_handle))

    print(dir(blast_hits))
    # remove("tmp/"+gene+"_tmp_blast.xml")

    # print(inspect.getmembers(blast_hits))

    '''
    
    result attributes: 'alignments', 'application', 'blast_cutoff', 'database', 'database_length', 'database_letters', 'database_name', 'database_sequences', 'date', 'descriptions', 'dropoff_1st_pass', 'effective_database_length', 'effective_hsp_length', 'effective_query_length', 'effective_search_space', 'effective_search_space_used', 'expect', 'filter', 'frameshift', 'gap_penalties', 'gap_trigger', 'gap_x_dropoff', 'gap_x_dropoff_final', 'gapped', 'hsps_gapped', 'hsps_no_gap', 'hsps_prelim_gapped', 'hsps_prelim_gapped_attemped', 'ka_params', 'ka_params_gap', 'matrix', 'multiple_alignment', 'num_good_extends', 'num_hits', 'num_letters_in_database', 'num_seqs_better_e', 'num_sequences', 'num_sequences_in_database', 'posted_date', 'query', 'query_id', 'query_length', 'query_letters', 'reference', 'sc_match', 'sc_mismatch', 'threshold', 'version', 'window_size']
    alignment attributes: 'accession', 'hit_def', 'hit_id', 'hsps', 'length', 'title']
    hsp attributes: 'align_length', 'bits', 'expect', 'frame', 'gaps', 'identities', 'match', 'num_alignments', 'positives', 'query', 'query_end', 'query_start', 'sbjct', 'sbjct_end', 'sbjct_start', 'score', 'strand'
    
    '''

    for result in blast_hits:
        # print(inspect.getmembers(result))
        # print(dir(result))
        #print(result.query)
        for alignment in result.alignments:
            # print(dir(alignment))
            # print(alignment.title) #gnl|BL_ORD_ID|43 STMMW_16561--5_1
            # print(alignment.accession) #43
            # print(alignment.hit_def) # STMMW_16561--5_1
            # print(alignment.hit_id) #gnl|BL_ORD_ID|43
            # print(alignment.length)
            for hsp in alignment.hsps:
                # print(dir(hsp))
                if hsp.identities == alignment.length:
                    perfect_hit_allele = alignment.title.split(" ")[-1]
                    perfect_hit_locus = perfect_hit_allele.split("-")[0]
                    # print(result.query)
                    # print(result.query_length)
                    # print(perfect_hit_allele)
                    # try:
                    partial_hit.remove(perfect_hit_locus)
                    # if "STM3752-" in perfect_hit_allele:
                    #     print("ERROR: " + perfect_hit_locus)
                    print("\n\n"+alignment.hit_def)
                    print('allele length:',allele_sizes[alignment.hit_def])
                    print('alignment length:', alignment.length)
                    print('identities:', hsp.identities)
                    print('gaps:', hsp.gaps)
                    #     print('query start', hsp.query_start)
                    #     print('query end: ', hsp.query_start + alignment.length)
                    #     print(hsp.query)
                    #     print(hsp.match)
                    #     print(hsp.sbjct)
                    # except:
                    #     print("ERROR: "+ perfect_hit_locus)
                    #     print(alignment.title)
                    #     print('length:', alignment.length)
                    #     print('identities:', hsp.identities)
                    #     print('gaps:', hsp.gaps)
                    #     print('query start: ', hsp.query_start)
                    #     print('query end: ',hsp.query_start+alignment.length)
                    #     print(hsp.query)
                    #     print(hsp.match)
                    #     print(hsp.sbjct)
                    # print(perfect_hit_allele)
                    exacthits +=1
                    totalhits += 1
                elif hsp.identities > (alignment.length*0.5) and hsp.identities < alignment.length:
                    # print(gene)
                    # print('length:', alignment.length)
                    # print('identities:', hsp.identities)
                    # print(hsp.query)
                    # print(hsp.match)
                    # print(hsp.sbjct)
                    totalhits +=1
    print(",".join(partial_hit))

    # for gene in schemelist:
    #     r_handle = open("tmp/"+gene+"_tmp_blast.xml")
    #
    #     blast_hits = list(NCBIXML.parse(r_handle))
    #     remove("tmp/"+gene+"_tmp_blast.xml")
    #     for result in blast_hits:
    #         for alignment in result.alignments:
    #             for hsp in alignment.hsps:
    #                 if hsp.identities > (alignment.length*0.5) and hsp.identities < alignment.length:
    #                     # print(gene)
    #                     # print('length:', alignment.length)
    #                     # print('identities:', hsp.identities)
    #                     # print(hsp.query)
    #                     # print(hsp.match)
    #                     # print(hsp.sbjct)
    #                     totalhits +=1
    #                 elif hsp.identities == alignment.length:
    #                     exacthits +=1
    #                     totalhits += 1

    # shutil.rmtree('tmp')

    print("exact hits: "+str(exacthits)+", total hits: "+str(totalhits)+", remaining loci to check: "+ str(len(partial_hit)))
    print(("--- %s seconds ---" % (time.time() - start_time)))





# run_blast()
#20 ~6s
#30 ~10s
#100 ~30s

if __name__ == "__main__":
    main(sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4])