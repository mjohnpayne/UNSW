from Bio import Phylo
from Bio.Phylo.Consensus import *
from Bio import AlignIO
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor
from Bio.Phylo.TreeConstruction import DistanceCalculator
import sys


calculator = DistanceCalculator('blosum62')
#constructor = DistanceTreeConstructor(calculator)



#msa = AlignIO.read('/Users/mjohnpayne/Documents/UNSW/Salmonella/data_aquisition/enterobase_data/concat_MLST_seqs.fasta', 'fasta')

msa = AlignIO.read(sys.argv[1], 'fasta')

constructor = DistanceTreeConstructor(calculator, 'nj')

#tree = constructor.build_tree(msa)

consensus_tree = bootstrap_consensus(msa, sys.argv[3], constructor, majority_consensus)

Phylo.write(consensus_tree,[sys.argv[2]],"newick")