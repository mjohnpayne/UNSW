from time import sleep as sl
import glob
from os.path import basename
from os import remove
import os
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Blast import NCBIXML
import glob
from multiprocessing.dummy import Pool as ThreadPool
import multiprocessing
import time
import subprocess as sub
import itertools
import shutil
import sys
from Bio import SeqIO
from Bio import SeqRecord
from Bio.Seq import Seq
import inspect
import re


start_time = time.time()

"""
1 - Pick scheme
2 - Blast all loci in scheme
3 - Perfect matches store match and fact of match
4 - any loci with no perfect match store location in query genome from matches to intact alleles
5 - blast newly defined allele region against complete allele set
6 - call alleles
"""

def get_allele_pos(allele_location):
    """
    :param allele_location: allele location file
    :return: dictionary of {allele name: [allele_st,allele_end,allele_direction]}
    """
    inf = [x.split('\t') for x in open(allele_location,"r").readlines()]
    outd = {y[0]:y[1:] for y in inf}
    return outd


def makedb(fasta,dbfolder):
    """

    :param fasta: fasta file to make database with
    :param dbfolder: folder to write db to
    :return: nothing BUT runds makeblastdb to generate the blast database required
    """
    out = fasta.split("/")[-1]
    path = dbfolder+"/"+out
    p = sub.Popen(['makeblastdb', '-dbtype','nucl','-in',fasta,"-out",path],stdout=sub.PIPE,stderr=sub.PIPE)
    p.communicate()


def mean(numbers):
    return float(sum(numbers)) / max(len(numbers), 1)


def run_blast(query_seq,locus_db,wordsize):
    """

    :param query_seq: query sequence - can be multiple fasta seqs
    :param locus_db: blastdb path
    :param wordsize: blast word size
    :return: returns list of blast results
    """
    cpus = multiprocessing.cpu_count()
    gene = basename(locus_db).replace(".fasta","")
    tmp_out = "tmp/"+gene+"_tmp_blast.xml"
    cline = NcbiblastnCommandline(
        query=query_seq,
        db=locus_db,
        evalue=10,
        perc_identity=90,
        out=tmp_out,
        outfmt=5,
        max_target_seqs=10000,
        max_hsps=3,
        word_size=wordsize,
        num_threads=cpus,
        task="dc-megablast")

    stdout, stderr = cline()

    r_handle = open(tmp_out)

    blast_records = list(NCBIXML.parse(r_handle))

    remove(tmp_out)
    """
    blast_records structure: 
    list of results (if multifasta input, one result per fasta seq) 
    
    result attributes: >>>'alignments'<<<, 'application', 'blast_cutoff', 'database', 'database_length', 'database_letters', 'database_name', 'database_sequences', 'date', 'descriptions', 'dropoff_1st_pass', 'effective_database_length', 'effective_hsp_length', 'effective_query_length', 'effective_search_space', 'effective_search_space_used', 'expect', 'filter', 'frameshift', 'gap_penalties', 'gap_trigger', 'gap_x_dropoff', 'gap_x_dropoff_final', 'gapped', 'hsps_gapped', 'hsps_no_gap', 'hsps_prelim_gapped', 'hsps_prelim_gapped_attemped', 'ka_params', 'ka_params_gap', 'matrix', 'multiple_alignment', 'num_good_extends', 'num_hits', 'num_letters_in_database', 'num_seqs_better_e', 'num_sequences', 'num_sequences_in_database', 'posted_date', 'query', 'query_id', 'query_length', 'query_letters', 'reference', 'sc_match', 'sc_mismatch', 'threshold', 'version', 'window_size']
        alignment attributes: 'accession', 'hit_def', 'hit_id', >>>'hsps'<<<, 'length', 'title']
            hsp attributes: 'align_length', 'bits', 'expect', 'frame', 'gaps', 'identities', 'match', 'num_alignments', 'positives', 'query', 'query_end', 'query_start', 'sbjct', 'sbjct_end', 'sbjct_start', 'score', 'strand'
    """

    return blast_records


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
    """

    :param allele_fasta: fasta path for alleles of a locus
    :return: word size for use in blast that is 20% of the average locus length
    """
    alleles = SeqIO.parse(allele_fasta,"fasta")
    sizes = []
    sizesdict = {}
    for i in alleles:
        sizes.append(len(i.seq))
    avsize = mean(sizes)
    wordsize = int(avsize*0.2)
    # print(wordsize)
    return wordsize

def print_blasthsp(result,alignment,hsp,alignment_locus):
    print("\n\n")
    print(alignment.hit_def)
    # print('allele size:', sizes[alignment_locus + "-1"])
    print('alignment length:', hsp.align_length)
    print('identities:', hsp.identities)
    print('gaps:', hsp.gaps)
    print('query contig', result.query)
    print('query start', hsp.query_start)
    print('query end: ', hsp.query_end)
    print('subject start', hsp.sbjct_start)
    print('subject end: ', hsp.sbjct_end)
    print(hsp.query)
    print(hsp.match)
    print(hsp.sbjct)

def get_exactmatches(blasthits,partial_hit,allele_sizes):
    """

    :param blasthits: parsed blast output as list of qury seq results
    :param partial_hit:list of loci, when an exact hit is found it is removed from list - resulting in a list of loci still to be examined
    :param allele_sizes: allele lengths of reference alleles
    :return: alleles that match exactly for length with no gaps and 100% identity - one for each locus - alleles with Ns are ignored
    """
    exacthits = 0
    totalhits = 0
    perfect_hit = {}
    # print(partial_hit)
    for result in blasthits:
        for alignment in result.alignments:
            for hsp in alignment.hsps:
                if hsp.identities == alignment.length and hsp.gaps == 0 and alignment.length == allele_sizes[alignment.hit_def]:
                    perfect_hit_allele = alignment.title.split(" ")[-1]
                    perfect_hit_locus = perfect_hit_allele.split("-")[0]
                    # print(perfect_hit_locus)
                    perfect_hit[perfect_hit_locus] = perfect_hit_allele
                    partial_hit.remove(perfect_hit_locus)
                    exacthits +=1
                    totalhits += 1
                elif hsp.identities > (alignment.length*0.5) and hsp.identities < alignment.length:
                    totalhits +=1

    print("exact hits: " + str(exacthits) + ", total hits: " + str(totalhits) + ", remaining loci to check: " + str(len(partial_hit)))

    return partial_hit, perfect_hit


def get_partial_match_query_region(blast_results,partial_matches,sizes):
    """

    :param blast_results: 1st round blast results
    :param partial_matches: loci without exact matches
    :return: partials dictionary {locus:list of tuples} tuple -> (hsp matching reference allele,query contig where match was found)
    Also writes fasta file of these hits for each locus to be used to blast all alleles of the locus
    """

    # print(partial_matches)
    partials = {}
    for result in blast_results:
        for alignment in result.alignments:
            alignment_locus = alignment.hit_def.rsplit("-")[0]
            if alignment_locus in partial_matches and alignment.hit_def == alignment_locus+"-1":
                for hsp in alignment.hsps:
                    if 'STMMW_00001' in alignment_locus:
                        print("\n\n")
                        print(alignment.hit_def)
                        # print('allele size:',sizes[alignment_locus+"-1"])
                        print('alignment length:', hsp.align_length)
                        print('identities:', hsp.identities)
                        print('gaps:', hsp.gaps)
                        print('query contig', result.query)
                        print('query start', hsp.query_start)
                        print('query end: ', hsp.query_end)
                        print('subject start', hsp.sbjct_start)
                        print('subject end: ', hsp.sbjct_end)
                        print(hsp.query)
                        print(hsp.match)
                        print(hsp.sbjct)
                    if alignment_locus not in partials:
                        partials[alignment_locus] = [(hsp,result.query)]
                    else:
                        partials[alignment_locus].append((hsp,result.query))

                    ### for queries with non intact hit to ref allele check other hits to see if remainder of locus is present before calling deletion - could be caused by assembly error or repetitive insertion
                    ### need to deal with split matches that overlap - STMMW_45221-1 in current example 100 gene scheme
    for locus in partials:
        hsplis = []
        c = 1
        for hsp in partials[locus]:
            hitseq = Seq(hsp[0].query.replace("-",""))
            # hitseq = only_nucs(hitseq)
            s = SeqRecord.SeqRecord(hitseq,locus+"_hit_no_"+str(c),description="")
            hsplis.append(s)
            c+=1
        SeqIO.write(hsplis,"tmp/"+locus+"_hits.fasta","fasta")
    return partials


def only_nucs(instring):
    outstring = ''
    for n in instring:
        if n in ["A","T","G","C","N","a","t","c","g","n"]:
            outstring+=n
    return outstring


def run_secondary_N_blast(partialhsps,allelefastas):
    """
    Make blast db from alleles for locus of interest
    blast hit regions in query genome against db
    :param partialhsps: all hsps matching to loci without exact alleles
    :param allele_fastas: dictionary of {locus:alleles path}
    :return: blast results in same format as run_blast() function
    """
    secresults = {}
    for locus in partialhsps:
        for i in allelefastas:
            if i == locus:
                makedb(allelefastas[locus],"tmp")
        sec_blast_results = run_blast("tmp/"+locus+"_hits.fasta","tmp/"+locus+".fasta",11)
        secresults[locus] = sec_blast_results
    return secresults


def return_n_ignored_hit_identities(hsp):
    # TODO:may need to check for locations where N nucleotides align with gaps - not sure what to do there...
    """
    :param hsp: hsp from a blast output
    :return: number of matches in a HSP including positions where either query or subject have an N
    """
    ident = int(hsp.identities)
    for pos in range(len(hsp.query)):
        if hsp.query[pos] == "N" or hsp.query[pos] == "n" or hsp.sbjct[pos] == "N" or hsp.sbjct[pos] == "n":
            ident +=1
    return ident


def reverse_complement(dna):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A',"-":"-","N":"N"}
    return ''.join([complement[base] for base in dna[::-1]])


def find_muts_frm_sec_results(s_results,allele_sizes,ref_pos):
    # TODO: loci that have two intact hits with the reference - may need to mark as 0 unless we implement synteny identification
    """

    :param s_results: all secondary blast results
    :param allele_sizes: length of all alleles
    :param ref_pos: positions of loci in the reference (allows genome position naming of SNPs)
    :return: list of snps in (ref,pos,change) format relative to the allele
    also return deletions and insertions IN THE UNKNOWN ALLELE i.e. treat reference as unchanged
    del format (ref,pos,"del")
    ins format ("ins",pos,nuc) - insertions will be a series of changes at the same reference position
    remaining loci - loci that do not have a single hit that alignes to the reference from start to end
    multiple_hits - loci that have two intact hits with the reference - may need to mark as 0 unless we implement synteny identification somewhere
    newalleles - for a locus where there is one hit that spans the reference, the locus and sequence(with "-" removed) returned in dictionary of {locus:sequence}
    """

    snps = {}
    deletions = {}
    insertions = {}
    remaining_loci = []
    multiple_hits = []
    newalleles = {}
    anymuts = {}
    for locus in s_results:
        remaining_loci.append(locus)
        snps[locus] = []
        deletions[locus] = []
        insertions[locus] = []
        anymuts[locus] = []
        genome_pos = int(ref_pos[locus][0])
        for result in s_results[locus]:
            for alignment in result.alignments:
                if alignment.hit_def == locus + "-1":
                    alignment_locus = alignment.hit_def.rsplit("-")[0]
                    for hsp in alignment.hsps:
                        orient = "+"
                        unknown_allele = hsp.query
                        ref_allele = hsp.sbjct
                        matches = hsp.match
                        refstart = hsp.sbjct_start
                        refend = hsp.sbjct_end
                        ## re-orient matches so positions are always in direction of reference allele
                        if int(hsp.sbjct_start) > int(hsp.sbjct_end):
                            orient = "-"
                            unknown_allele = reverse_complement(unknown_allele)
                            ref_allele = reverse_complement(ref_allele)
                            matches = matches[::-1]
                            refstart,refend = refend,refstart

                        intact_nucs = ["A", "T", "G", "C", "a", "t", "c", "g"]
                        for i in range(len(unknown_allele)):
                            if matches[i] == " ":
                                if unknown_allele[i] in intact_nucs and ref_allele[i] in intact_nucs:
                                    pos = i - ref_allele[:i].count(
                                        "-")  # correct for reference allele length if it is changed by upstream insertions in the query
                                    nuc = (ref_allele[i], genome_pos + pos, unknown_allele[i])
                                    snps[locus].append(nuc)
                                    anymuts[locus].append(nuc)
                                elif unknown_allele[i] == "-" and ref_allele[i] in intact_nucs:
                                    pos = i - ref_allele[:i].count("-")
                                    d = (ref_allele[i], genome_pos + pos, "del")
                                    deletions[locus].append(d)
                                    anymuts[locus].append(d)
                                elif unknown_allele[i] in intact_nucs and ref_allele[i] == "-":
                                    pos = i - ref_allele[:i].count("-")
                                    ins=("ins", genome_pos + pos, unknown_allele[i])
                                    insertions[locus].append(ins)
                                    anymuts[locus].append(ins)

                        if refstart == 1 and refend == allele_sizes[locus+"-1"]: # this step gets hit that covers reference allele
                            if locus in remaining_loci and "N" not in unknown_allele:
                                remaining_loci.remove(locus)
                            else:
                                multiple_hits.append(locus)
                            if "N" not in unknown_allele:
                                newalleles[locus] = unknown_allele.replace("-","")

    return anymuts,snps,deletions,insertions,remaining_loci,multiple_hits,newalleles


def make_intact_allele_fasta(scheme_fastas,dbname):
    """

    :param scheme_fastas: fasta file paths for scheme of interest
    :param dbname: output fasta path
    :return:fasta file with all intact alleles derived from fasta files of loci making up a scheme
    """
    outseqs = []
    for fname in scheme_fastas:
        infile = SeqIO.parse(fname,"fasta")
        for allele in infile:
            if "--" not in allele.id:
                outseqs.append(allele)
    SeqIO.write(outseqs,dbname,"fasta")

def remove_Ns(inseq):
    outname = inseq.replace(".fasta","_no_n.fasta")
    seqs = SeqIO.parse(inseq,"fasta")
    no_n_ls = []
    for s in seqs:
        s2 = s
        s2.seq = Seq(str(s.seq).replace("N",""))
        no_n_ls.append(s2)
    SeqIO.write(no_n_ls,outname,"fasta")
    return outname

def check_ends(contig,qstart,qend,sstart,send,reflen):
    ## if hit in reverse orientation have to check opposite end of hit in query genome

    new_qstart = 0
    new_qend = 0
    added_start = ''
    added_end = ''



    if int(sstart) > int(send):
        new_qstart = int(qstart) - (reflen - int(sstart))
        new_qend = int(qend) + int(send) - 1


    elif int(sstart) < int(send):
        new_qstart = int(qstart) - int(sstart)
        new_qend = int(qend) + (reflen - int(send) - 1)

    added_start = contig[new_qstart:int(qstart) - 1]
    added_end = contig[int(qend):new_qend - 1]

    ## if missing flanking region is not all Ns then use original query hit edge - allows for rearrangements/deletions where flanking regions are really gone
    if added_start.count("N") != len(added_start):
        new_qstart = qstart
        added_start = ["del",len(added_start)]
    if added_end.count("N") != len(added_end):
        new_qend = qend
        added_end = ["del", len(added_end)]
    full_allele = contig[new_qstart:new_qend]

    return full_allele,added_start,added_end

def check_all_orient(hsp_list):
    orient = []
    for tup in hsp_list:
        hsp = tup[0]
        if hsp.sbjct_end < hsp.sbjct_start:
            orient.append("negative")
        else:
            orient.append("positive")
    if len(list(set(orient))) > 1:
        return "mixed"
    else:
        return orient[0]

def largest_nonn_strings(string):
    matches = re.findall(r"[ATGC]*",string)
    mx = max([len(x) for x in matches])
    return mx

def check_mid(contig,hspls,q_genome,wordsize,reflen):
    """
    use word length ratio as cutoff for allowing non-N letters in restored gap sequence(i.e. if word length is 12. non-N sequence chunk can be a max of 18 before mid regions is called something else)
    :param contig:
    :param hspls:
    :param q_genome:
    :return:
    """


    orient = check_all_orient(hspls)
    if orient == "mixed":
        # TODO work out what to do with mixed orientation
        return "mixed_orientation"
    else:
        ordered_by_query = []
        order = {}
        for tup in hspls:
            hsp = tup[0]
            order[hsp.query_start] = hsp

        sorted_hsps = sorted(map(int, order.keys()))
        full_allele = order[sorted_hsps[0]].query
        sstart = order[sorted_hsps[0]].sbjct_start
        send = order[sorted_hsps[-1]].sbjct_end
        qstart = order[sorted_hsps[0]].query_start
        qend = order[sorted_hsps[-1]].query_start

        for start in range(len(sorted_hsps)-1):
            hsp = order[sorted_hsps[start]]
            hspnext = order[sorted_hsps[start+1]]
            mid_section_st = int(hsp.query_end)
            mid_section_en = int(hspnext.query_start)
            mid_section_seq = q_genome[contig][mid_section_st:mid_section_en-1]
            full_allele += mid_section_seq
            full_allele += hspnext.query

            nonn_size = largest_nonn_strings(mid_section_seq)
            if nonn_size > (wordsize*2):
                return "possible_insertion"

        full_allele, added_start, added_end = check_ends(contig, qstart, qend, sstart, send, reflen)

        return full_allele











    #     nmin = 0
    #     nmax = 0
    #     if hsp.sbjct_end < hsp.sbjct_start:
    #         nmin = hsp.sbjct_end
    #         nmax = hsp.sbjct_start
    #     else:
    #         nmin = hsp.sbjct_start
    #         nmax = hsp.sbjct_end
    #     if nmin < min:
    #         min = nmin
    #     if nmax > max:
    #         max = nmax
    #
    #     hsp = tup[0]
    #     queryid = tup[1]
    #     contig = q_genome[queryid]
    #     print("\n")
    #     # print('allele length:', reflen)
    #     print('alignment length:', hsp.align_length)
    #     print('identities:', hsp.identities)
    #     print('gaps:', hsp.gaps)
    #     print('query contig', queryid)
    #     print('query start', hsp.query_start)
    #     print('query end: ', hsp.query_end)
    #     print('subject start', hsp.sbjct_start)
    #     print('subject end: ', hsp.sbjct_end)
    #     print(hsp.query)
    #     print(hsp.match)
    #     print(hsp.sbjct)
    # print(min,max)

def reconstruct_n_interupted_loci(locus,hspls,query_genome,reflen,missing_perc_cutoff,wordsize,s_results):
    query_genome = SeqIO.parse(query_genome,"fasta")
    q_genome = {}

    full_allele = ""

    for s in query_genome:
        q_genome[s.id] = str(s.seq)

    ##check that number of identities in blast hits is at least X fraction of normal reference allele length
    ident_no = 0
    for tup in hspls:
        matches = int(tup[0].identities)
        ident_no += matches
    fraction_ident = float(ident_no)/float(reflen)
    if fraction_ident < float(missing_perc_cutoff):
        # print(locus,"unscorable: number of Ns exceeds cutoff")
        return "unscorable_too_much_missing",""

    if len(hspls) > 1:
        contigs = {}
        hsp_to_investigate = []
        for tup in hspls:
            hsp = tup[0]
            contig = tup[1].split(" ")[0]
            if contig not in contigs:
                contigs[contig] = [tup]
            else:
                contigs[contig].append(tup)

        # if hits come from > 1 contig
        if len(contigs.keys()) > 1:
            return "split over contigs",""
        else:
            for c in contigs:
                fixed_mid = check_mid(c, contigs[c],q_genome,wordsize,reflen)
                if fixed_mid == "possible_insertion":
                    return "possible_insertion",""
                elif fixed_mid == "mixed_orientation":
                    return "mixed_orientation",""
                else:
                    full_allele = fixed_mid

    else:
        # print(locus+": "+hspls[0][1])
        hsp = hspls[0][0]
        queryid = hspls[0][1].split(" ")[0]
        contig = q_genome[queryid]
        # print("\n")
        # print('allele length:', reflen)
        # print('alignment length:', hsp.align_length)
        # print('identities:', hsp.identities)
        # print('gaps:', hsp.gaps)
        # print('query contig', queryid)
        # print('query start', hsp.query_start)
        # print('query end: ', hsp.query_end)
        # print('subject start', hsp.sbjct_start)
        # print('subject end: ', hsp.sbjct_end)
        # print(hsp.query)
        # print(hsp.match)
        # print(hsp.sbjct)

        full_allele,added_start,added_end = check_ends(contig,hsp.query_start,hsp.query_end,hsp.sbjct_start,hsp.sbjct_end,reflen)

    """
    Get partial matches to loci:
    when match has end missing check input genome for Ns back to start/end of locus
    - if Ns store this in allele
    - if not Ns store as truncated allele

    when multiple hits check if on same contig:
    - if yes then check if gap between is Ns or something else
        - if Ns reconstruct allele
        - if not Ns - potential large indel or other rearrangement
    - if no split allele - return pieces for other analysis

     - if no snps or indels of another type - need to generate -1_X
     - if snps or indels store
    :param partial_hits_to_ref_allele: for each locus existing hits to the ref (-1) allele
    :param query_genome:
    :return:
    """
    if float(len(full_allele)) > 1.5*float(reflen):
        # print(locus,"unscorable: N-locus too long!!")
        return "unscorable_too_long",""

    # need to get closest hit of non N region/s and record muts

    # if single region then get list of alleles with 0 indels or snps, if either is present > new allele
    exact_ns = ""
    idents = []
    allele = []
    mutno = []
    sbjct_n_no = []
    #count how many mismatches there are to each alignment as well as how many identities (i.e. across multiple hsps if present), also need to take into account large dels from
    for result in s_results[locus]:
        for alignment in result.alignments:
            muts = 0
            ident = 0
            sbjct_n = 0
            alignment_locus = alignment.hit_def.rsplit("-")[0]
            for hsp in alignment.hsps:
                ident += int(hsp.identities)
                for pos in range(len(hsp.match)):
                    if " " in hsp.match[pos] and hsp.query[pos] != "N" and hsp.sbjct[pos] != "N":
                        muts+=1
                    if hsp.sbjct[pos] == "N":
                        sbjct_n += 1
            allele.append(alignment.hit_def)
            idents.append(ident)
            mutno.append(muts)
            sbjct_n_no.append(sbjct_n)
    #get the highest identity score in all alignments
    maxident = max(idents)
    # get all allele matches that have the highest identity
    tophits = []
    muts = {}
    for hit in range(len(idents)):
        if idents[hit] == maxident:
            if mutno[hit] == 0 and sbjct_n_no[hit] == 0:
                # print(allele[hit],idents[hit])
                tophits.append(allele[hit])
                muts[allele[hit]]="no"
            elif sbjct_n_no[hit] > 0 and mutno[hit] == 0:
                tophits.append(allele[hit])
                print("Ns",allele[hit], idents[hit])
                muts[allele[hit]] = "no"
            elif sbjct_n_no[hit] > 0 and mutno[hit] > 0:
                print("MUTS + Ns", allele[hit], idents[hit])
                tophits.append(allele[hit])
                muts[allele[hit]] = "yes"
            else:
                print("MUTS", allele[hit], idents[hit])
                tophits.append(allele[hit])
                muts[allele[hit]] = "yes"

    # if >1 +ve hit mark as unsure
    # if 1 +ve mark as -ve version of allele and check if it already exists
    # if 1 or more -ve versions of the same allele with no others mark as -ve version of allele and check if it already exists
    intact = [x for x in tophits if "--" not in x]
    if len(intact) > 1:
        return "mutiple_tophits",full_allele
    elif len(intact) == 1 and muts[intact[0]] == "no":
        return "new neg allele",[intact[0].split("-")[-1],full_allele]
    elif len(intact) == 1 and muts[intact[0]] == "yes":
        if "N" in full_allele:
            return "novel neg allele", full_allele
        else:
            # TODO : need to record that these cases represent large deletions relative to the reference allele
            return "novel pos allele", full_allele
    elif len(intact) == 0:
        alleles = list(set([x.split("-")[-1].split("_")[0] for x in tophits]))

        if len(alleles) > 1:
            return "mutiple_neg_tophits",full_allele
        elif len(alleles) == 1 and muts[alleles[0]] == "no":
            return "new neg allele",["-"+alleles[0],full_allele]
        elif len(intact) == 1 and muts[alleles[0]] == "yes":
            if "N" in full_allele:
                return "novel neg allele", full_allele
            else:
                #TODO : need to record that these cases represent large deletions relative to the reference allele
                return "novel pos allele", full_allele
        else:
            return "nohits",full_allele


    return "other_investigate",full_allele

    # if >1 region get alleles that are exact mathes to all regions if any snps or indels to top hit > new allele


    # print("\nSequence added at start: \n",added_start,"\n","Sequence added at end: \n",added_end,"\n",full_allele)

    # print(locus,hspls)
    # return full_allele


def main(schemelists, alleles, query_genome,locus_refrence_locs,scheme,wordsize):

    cpus = multiprocessing.cpu_count()

    if os.path.exists("tmp"):
        shutil.rmtree('tmp')
        os.mkdir("tmp")
    else:
        os.mkdir("tmp")
    inschemes = schemelists + "/*.txt"

    alleles_folder = alleles + "/*.fasta"

    # gets dictionary with lists of loci for each scheme
    scheme_gene_lists = return_scheme_loci_ls(inschemes)

    # gets list of allele fasta file paths
    allele_fastas = get_allele_fastas(alleles_folder)

    # gets list of loci in scheme of interest
    schemelist = scheme_gene_lists[scheme]

    # gets allele reference contig, positions and orientation dict
    locus_positions = get_allele_pos(locus_refrence_locs)

    wordsizes = []
    scheme_fastas = []

    ## makes a list of all of the fasta files in the scheme of interest
    ## also makes list of word size to use depending on ref allele size
    for i in schemelist:
        allele_file = allele_fastas[i]
        scheme_fastas.append(allele_file)
        ws = get_word_size(allele_file)
        wordsizes.append(ws)


    ##writes all fasta files in scheme to tmpdb location -- effectively a python "cat"
    tmpdb = "tmp/"+scheme+"_alleles.fasta"

    # with open(tmpdb, 'w') as outfile:
    #     for fname in scheme_fastas:
    #         with open(fname) as infile:
    #             outfile.write(infile.read())
    make_intact_allele_fasta(scheme_fastas,tmpdb)

    ##generates blast db from above concatenated fasta file
    makedb(tmpdb, "tmp/")


    ##load scheme allele fastas as Seq object
    allele_sizes = {}

    scheme_alleles = SeqIO.parse(tmpdb,"fasta")

    ## make dict of allele sizes
    for allele in scheme_alleles:
        allele_sizes[allele.id] = len(allele.seq)


    no_n_query_genome = remove_Ns(query_genome)

    ##return blast results
    # blast_hits = run_blast(no_n_query_genome, tmpdb, 12)
    blast_hits = run_blast(query_genome, tmpdb, wordsize)

    ##complete loci are exact matches to an allele already in the db
    ##partial loci are the rest
    partial_loci,complete_loci = get_exactmatches(blast_hits,schemelist,allele_sizes)

    ##returns hsps that give the positions of matches in the query genomes that match the "1"/reference allele
    ## also writes these regions to one fasta file for each locus different hsps are called "_hit_no_1", "hit_no_2" etc
    partial_hsps = get_partial_match_query_region(blast_hits,partial_loci,allele_sizes)

    # for i in complete_loci:
    #     print(i,complete_loci[i])
    #     sl(0.2)

    # for locus in partial_hsps:
    #     print("\n")
    #     print(locus)
    #     print("ref allele length:",allele_sizes[locus+"-1"])
    #     for tup in partial_hsps[locus]:
    #         hsp = tup[0]
    #         queryid = tup[1]
    #         print("\n")
    #         # print('allele length:', allele_sizes[alignment.hit_def])
    #         print('alignment length:', hsp.align_length)
    #         print('identities:', hsp.identities)
    #         print('gaps:', hsp.gaps)
    #         print('query contig', queryid)
    #         print('query start', hsp.query_start)
    #         print('query end: ', hsp.query_end)
    #         print('subject start', hsp.sbjct_start)
    #         print('subject end: ', hsp.sbjct_end)
    #         print(hsp.query)
    #         print(hsp.match)
    #         print(hsp.sbjct)

    ## run blast for each locus with 'query genome region that corresponds to "reference" allele hit regions' against all alleles for that locus

    sec_results = run_secondary_N_blast(partial_hsps,allele_fastas)

# -----------------
    # different process for split hits
    # check for Ns
    # check for missing ends
    # check for indels - always call relative to existing allele
    #       if gap in subject = insertion in query
    #       if gap in query = deletion in query
    # check for SNPs
# ------------------

    ### if match against reference runs from reference start to end then snps/indesl can be called without checking for missing ends
    ### therefore can extract snps from alignment against "-1" reference allele using locus blast results - find_snps function

    all_muts,snp_dict,deldict,insdict,remaining_loci,mult_hits,new_alleles = find_muts_frm_sec_results(sec_results,allele_sizes,locus_positions)

    # for i in new_alleles:
    #     print(i,snp_dict[i],deldict[i],insdict[i],new_alleles[i])
    #     sl(0.2)

    has_muts = []

    inss = 0
    dels = 0
    snps = 0
    bothinsdel = 0

    for locus in all_muts:
        if len(snp_dict[locus])==0 and len(deldict[locus])>0:
            dels+=1
        if len(snp_dict[locus])==0 and (len(deldict[locus])>0 or len(insdict[locus])>0):
            bothinsdel+=1
        if len(snp_dict[locus])==0 and len(insdict[locus])>0:
            inss+=1
        if len(snp_dict[locus])>0:
            snps+=1
        if len(snp_dict[locus])>0 or len(deldict[locus]) > 0 or len(insdict[locus]) > 0:
            has_muts.append(locus)
        print("\n",locus)
        print(all_muts[locus])
        print(snp_dict[locus])
        print(deldict[locus])
        print(insdict[locus])
        # print(new_alleles[locus])

    muts_only = has_muts
    unscorable = []
    split_over_contigs = []
    multiple_hits_mixed_orientation = []
    possible_mid_locus_insertion = []
    split_over_contigs_with_muts = []
    n_loci = {}
    nloci_ls = []
    no_hits = []
    unaccounted = schemelist
    # print("exact match")
    # for locus in complete_loci:
    #     print(locus,complete_loci[locus])

    print("Non-exact match no Ns")

    for locus in has_muts:
        print(locus,"new pos allele")
        if locus in new_alleles:
            print(new_alleles[locus])

    print("non-exact match with Ns")
    for locus in remaining_loci:
        ## return query loci with N positions restored (large regions or regions at loci ends will be missing from hsps)
        out1,out2 = reconstruct_n_interupted_loci(locus,partial_hsps[locus],query_genome,allele_sizes[locus + "-1"],0.50,wordsize,sec_results)
        n_loci[locus] = out1
        print(locus,out1,out2)
        if out1 == "other_investigate":
            print(out2)
        if out1 == "unscorable":
            unscorable.append(locus)
            unaccounted.remove(locus)
            if locus in has_muts:
                muts_only.remove(locus)
        elif out1 == "split over contigs":
            if locus not in has_muts:
                split_over_contigs.append(locus)
                unaccounted.remove(locus)
            else:
                split_over_contigs_with_muts.append(locus)
                muts_only.remove(locus)
                unaccounted.remove(locus)
        elif out1 == "mixed_orientation":
            multiple_hits_mixed_orientation.append(locus)
            unaccounted.remove(locus)
        elif out1 == "possible_insertion":
            possible_mid_locus_insertion.append(locus)
            unaccounted.remove(locus)
        elif out1 == "":
            no_hits.append(locus)
            unaccounted.remove(locus)
        elif out1 == "new neg allele":
            # print(locus,out1)
            nloci_ls.append(locus)
            unaccounted.remove(locus)
        elif out1 == "nohits":
            continue

    for i in muts_only:
        unaccounted.remove(i)

    ###TODO need to get categories properly organised so that totals add up - probably overlap of alleles with mutations and alleles with N regions
    # print("Perfect_hits: ",len(complete_loci.keys()))
    # print("intact with mutations: ",len(muts_only),"\n","\n".join(muts_only))
    # print("reconstructed 'N' containing:",len(nloci_ls),"\n","\n".join(nloci_ls))
    # print("hit split over different contigs:", len(split_over_contigs),"\n","\n".join(split_over_contigs))
    # print("hit split over different contigs with muts: \n", len(split_over_contigs_with_muts),"\n", "\n".join(split_over_contigs_with_muts))
    # print("mixed orientation:",len(multiple_hits_mixed_orientation),"\n","\n".join(multiple_hits_mixed_orientation))
    # print("possible mid locus insertion:", len(possible_mid_locus_insertion),"\n","\n".join(possible_mid_locus_insertion))
    # print("unscorable:",len(unscorable),"\n".join(unscorable))
    # print("loci with no hits in query:", len(no_hits),"\n","\n".join(no_hits))
    # print("loci remaining:", len(unaccounted),"\n","\n".join(unaccounted))


    # print(unscorable)
    # print(split_over_contigs)

        # print("\n")
        # print(locus)
        # print("ref allele length:", allele_sizes[locus + "-1"])
        # for tup in partial_hsps[locus]:
        #     hsp = tup[0]
        #     queryid = tup[1]
        #     print("\n")
        #     # print('allele length:', allele_sizes[alignment.hit_def])
        #     print('alignment length:', hsp.align_length)
        #     print('identities:', hsp.identities)
        #     print('gaps:', hsp.gaps)
        #     print('query contig', queryid)
        #     print('query start', hsp.query_start)
        #     print('query end: ', hsp.query_end)
        #     print('subject start', hsp.sbjct_start)
        #     print('subject end: ', hsp.sbjct_end)
        #     print(hsp.query)
        #     print(hsp.match)
        #     print(hsp.sbjct)



    #
    # print("loci with snps: ",snps)
    # print("loci with deletions no snps: ",dels)
    # print("loci with insertions no snps:",inss)
    # print("loci with insertions or deletions, no snps:",bothinsdel)
    # print("Remaining with 1 or less intact hits:",remaining_loci)
    # print("Remaining with multiple intact hits:",mult_hits)

    # print("query_hit\tallele\tfull allele length\talignment length\tn_corr_ident\tgaps")
    # for locus in sec_results:
    #     if locus == "STMMW_25991":
    #         print(locus)
    #         for result in sec_results[locus]:
    #             for alignment in result.alignments:
    #                 alignment_locus = alignment.hit_def.rsplit("-")[0]
    #                 for hsp in alignment.hsps:
    #                     if hsp.align_length >= allele_sizes["STMMW_25991-1"]:
    #                         ncor_ident = return_n_ignored_hit_identities(hsp)
    #                         # if hsp.align_length == ncor_ident:
    #                         outp = [result.query,alignment.hit_def,allele_sizes["STMMW_25991-1"],hsp.align_length,ncor_ident,hsp.gaps]
    #                         outp = map(str,outp)
    #                         print("\t".join(outp))
    #                         # print("\n\n")
    #                         # print(alignment.hit_def)
    #                         # print('alignment length:', hsp.align_length)
    #                         # print('identities:', hsp.identities)
    #                         # print('n_corr identities',ncor_ident)
    #                         # print('gaps:', hsp.gaps)
    #                         # print('query contig', result.query)
    #                         # print('query start', hsp.query_start)
    #                         # print('query end: ', hsp.query_end)
    #                         # print('subject start', hsp.sbjct_start)
    #                         # print('subject end: ', hsp.sbjct_end)
    #                         print(hsp.query)
    #                         print(hsp.match)
    #                         print(hsp.sbjct)

    '''
    blast hits structure: 
    list of results
    
    result attributes: >>>'alignments'<<<, 'application', 'blast_cutoff', 'database', 'database_length', 'database_letters', 'database_name', 'database_sequences', 'date', 'descriptions', 'dropoff_1st_pass', 'effective_database_length', 'effective_hsp_length', 'effective_query_length', 'effective_search_space', 'effective_search_space_used', 'expect', 'filter', 'frameshift', 'gap_penalties', 'gap_trigger', 'gap_x_dropoff', 'gap_x_dropoff_final', 'gapped', 'hsps_gapped', 'hsps_no_gap', 'hsps_prelim_gapped', 'hsps_prelim_gapped_attemped', 'ka_params', 'ka_params_gap', 'matrix', 'multiple_alignment', 'num_good_extends', 'num_hits', 'num_letters_in_database', 'num_seqs_better_e', 'num_sequences', 'num_sequences_in_database', 'posted_date', 'query', 'query_id', 'query_length', 'query_letters', 'reference', 'sc_match', 'sc_mismatch', 'threshold', 'version', 'window_size']
        alignment attributes: 'accession', 'hit_def', 'hit_id', >>>'hsps'<<<, 'length', 'title']
            hsp attributes: 'align_length', 'bits', 'expect', 'frame', 'gaps', 'identities', 'match', 'num_alignments', 'positives', 'query', 'query_end', 'query_start', 'sbjct', 'sbjct_end', 'sbjct_start', 'score', 'strand'
    
    '''


    # print(",".join(partial_hit))
    # shutil.rmtree('tmp')


    print(("--- %s seconds ---" % (time.time() - start_time)))





# run_blast()
#20 ~6s
#30 ~10s
#100 ~30s


# main("/Users/michaelpayne/Documents/UNSW/Salmonella/new_typing_scheme/scheme_development/betas/cgMLST_b4/hierMLST_stm_loci_accession_lists",
# "/Users/michaelpayne/Documents/UNSW/Salmonella/new_typing_scheme/allele_scheme_dev/blast_testing/hierMLST_stm_alleles",
# "/Users/michaelpayne/Documents/UNSW/Salmonella/new_typing_scheme/allele_scheme_dev/dev_genomes/SRR2532791_nullarbor_megahit.fa",
# "/Users/michaelpayne/Documents/UNSW/Salmonella/new_typing_scheme/allele_scheme_dev/blast_testing/hierMLST_stm_allele_locations.txt",sys.argv[1],12)


# main("/Users/michaelpayne/Documents/UNSW/Salmonella/new_typing_scheme/scheme_development/betas/cgMLST_b4/hierMLST_stm_loci_accession_lists",
# "/Users/michaelpayne/Documents/UNSW/Salmonella/new_typing_scheme/scheme_development/betas/cgMLST_b4/hierMLST_stm_alleles",
# "/Users/michaelpayne/Documents/UNSW/Salmonella/new_typing_scheme/allele_scheme_dev/dev_genomes/SRR2532791_skesa_assembly.fasta",
# "/Users/michaelpayne/Documents/UNSW/Salmonella/new_typing_scheme/allele_scheme_dev/blast_testing/hierMLST_stm_allele_locations.txt",sys.argv[1],12)

# main("/Users/michaelpayne/Documents/UNSW/Salmonella/new_typing_scheme/scheme_development/betas/cgMLST_b4/hierMLST_stm_loci_accession_lists",
# "/Users/michaelpayne/Documents/UNSW/Salmonella/new_typing_scheme/scheme_development/betas/cgMLST_b4/hierMLST_stm_alleles",
# "/Users/michaelpayne/Documents/UNSW/Salmonella/new_typing_scheme/allele_scheme_dev/dev_genomes/SRR2532791.fasta",'100')

# main("/Users/michaelpayne/Documents/UNSW/Salmonella/new_typing_scheme/scheme_development/betas/cgMLST_b4/hierMLST_stm_loci_accession_lists",
# "/Users/michaelpayne/Documents/UNSW/Salmonella/new_typing_scheme/scheme_development/betas/cgMLST_b4/hierMLST_stm_alleles",
# "/Users/michaelpayne/Documents/UNSW/Salmonella/new_typing_scheme/allele_scheme_dev/dev_genomes/SRR2532791_snpfilt_masked_contigs.fasta",
# "/Users/michaelpayne/Documents/UNSW/Salmonella/new_typing_scheme/allele_scheme_dev/blast_testing/hierMLST_stm_allele_locations.txt",sys.argv[1],12)

main("/Users/michaelpayne/Documents/UNSW/Salmonella/new_typing_scheme/scheme_development/betas/cgMLST_b4/hierMLST_stm_loci_accession_lists",
"/Users/michaelpayne/Documents/UNSW/Salmonella/new_typing_scheme/scheme_development/betas/cgMLST_b4/hierMLST_stm_alleles",
"/Users/michaelpayne/Documents/UNSW/Salmonella/new_typing_scheme/allele_scheme_dev/dev_genomes/SRR2532791_nullarbor_megahit.fa",
"/Users/michaelpayne/Documents/UNSW/Salmonella/new_typing_scheme/allele_scheme_dev/blast_testing/hierMLST_stm_allele_locations.txt",sys.argv[1],12)

# main(
#     "/Users/michaelpayne/Documents/UNSW/Salmonella/new_typing_scheme/scheme_development/betas/cgMLST_b4/hierMLST_stm_loci_accession_lists",
#     "/Users/michaelpayne/Documents/UNSW/Salmonella/new_typing_scheme/allele_scheme_dev/blast_testing/hierMLST_stm_alleles",
#     "/Users/michaelpayne/Documents/UNSW/Salmonella/new_typing_scheme/allele_scheme_dev/dev_genomes/SRR2532791_snpfilt_scaffolds.fasta",
#     '560')




# if __name__ == "__main__":
#     main(sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4])