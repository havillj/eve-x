from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio import SeqIO, Entrez
from dna_features_viewer import GraphicFeature, GraphicRecord, BiopythonTranslator
import matplotlib.pyplot as pyplot
import matplotlib
from pathlib import Path
import re
from lxml import etree as ET
import os
import numpy as np
import math
from Bio.Phylo.TreeConstruction import DistanceCalculator
from Bio import AlignIO
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor
from Bio import Phylo
from Bio.Cluster import treecluster 


def drawTree(fileName):
    aln = AlignIO.read(fileName, 'fasta') # Bio.Align.MultipleSeqAlignment object
    aln = aln[1:]
    calculator = DistanceCalculator('identity')
#    dm = calculator.get_distance(aln)
    constructor = DistanceTreeConstructor(calculator, 'nj')
    tree = constructor.build_tree(aln)
    tree.ladderize()
    
    # scorer = Phylo.TreeConstruction.ParsimonyScorer()
#     searcher = Phylo.TreeConstruction.NNITreeSearcher(scorer)
#     constructor = Phylo.TreeConstruction.ParsimonyTreeConstructor(searcher, tree)
#     pars_tree = constructor.build_tree(aln)


    fig = pyplot.figure(figsize=(10, 20), dpi=100)
    axes = fig.add_subplot(1, 1, 1)
    matplotlib.rc('font', size=6)
    Phylo.draw(tree, axes=axes)
    pyplot.savefig('tree.pdf')
    Phylo.write(tree, 'tree.newick', 'newick')
    
def cluster(fileName):
    aln = AlignIO.read(fileName, 'fasta') # Bio.Align.MultipleSeqAlignment object
    aln = aln[1:]
    calculator = DistanceCalculator('identity')
    dm = calculator.get_distance(aln)
    cluster = treecluster(distancematrix=dm)
    print(cluster)

    
#drawTree('/home/havill/data2/results4/viruses/sequences/Flaviviridae/sequences_NC_001564.2.fasta')
drawTree('/home/havill/data2/ortho_flanks.fasta')