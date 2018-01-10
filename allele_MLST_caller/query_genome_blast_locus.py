from time import sleep as sl

from Bio import SeqIO
from Bio.Blast.Applications import NcbiblastnCommandline

from Bio.Blast import NCBIXML
import time


def run_blast(query_genome,locus_db,tmp_out):

    cline = NcbiblastnCommandline(
        query=query_genome,
        db=locus_db,
        evalue=10,
        perc_identity=100,
        out=tmp_out,
        outfmt=5,
        max_target_seqs=10000000,
        max_hsps=3,
        word_size=20)

    stdout, stderr = cline()

    r_handle = open(
        "/Users/michaelpayne/Documents/UNSW/Salmonella/new_typing_scheme/allele_scheme_dev/blast_testing/LT2_blast_all_alleles.xml")

    blast_records = list(NCBIXML.parse(r_handle))

    return blast_records


if __name__ == "__main__":
    run_blast()