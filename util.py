import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as pyplot
from Bio import SeqIO, Entrez, AlignIO
from Bio.Align import MultipleSeqAlignment
from Bio.SeqFeature import SeqFeature, FeatureLocation
from dna_features_viewer import GraphicFeature, GraphicRecord, BiopythonTranslator
from intervaltree import Interval, IntervalTree
from lxml import etree as ET
from multiprocessing import Process
from pathlib import Path
import copy
import math
import numpy as np
import os
import pickle
import pysam
import re
import sys
import time
import urllib.request as web

SPADES_EXEC = '/Volumes/Data/bin/SPAdes-3.14.1-Linux/bin/spades.py'
VIRUS_DB = 'virusdb'

ROOT_DIR = '/Volumes/Data2/'  # '/home/havill/data/aegypti/analyzed/'
SPECIMENS_DIR = ROOT_DIR + 'specimens/'

RESULTS_DIR = ROOT_DIR + 'results/'
VIRUSES_DIR = RESULTS_DIR + 'viruses/'
SEQUENCES_DIR = VIRUSES_DIR + 'sequences/'

GB_DIR = ROOT_DIR + 'gb/'
FASTA_DIR = ROOT_DIR + 'fasta/'
GFF3_DIR = ROOT_DIR + 'gff3/'
BASE_FEATURES_GFF3_FILENAME = GFF3_DIR + 'Aedes-aegypti-LVP_AGWG_BASEFEATURES_AaegL5.2.gff3'
REPEAT_FEATURES_GFF3_FILENAME = GFF3_DIR + 'Aedes-aegypti-LVP_AGWG_REPEATFEATURES_AaegL5.gff3'
FEATURE_TREES_PICKLED_FILENAME = ROOT_DIR + 'featuretrees.pickle'

FAMILY_CSV = ROOT_DIR + 'families.csv'
LOGFILE_NAME = 'eve.log'
LOGFILE_PATH = Path(ROOT_DIR) / LOGFILE_NAME

MIN_PIDENT = 80
EVALUE_V = 1e-10
EVALUE_A = 1e-10
ALLOWED_OVERLAP = 30           # bp aedes hits may overlap ends of virus hits or each other
OVERLAP_FRACTION = 0.95        # if this fraction of a hit overlaps with an Aedes hit, it is considered to be present in the reference genome
MAX_FLANK_DISTANCE = 20000
COMPLEXITY_K = 3
COMPLEXITY_CUTOFF = 0.75
MIN_VIRUS_HIT_LENGTH = 100
#MIN_FLANKING_QSTART = 10        # min distance a flanking hit must be from either end of the contig  # REMOVE PROBABLY
MIN_FLANKING_HIT_LENGTH = 100    # min length of a flanking hit if it passes above test
MAX_FLANKING_DISTANCE = 200    
#MAX_FLANKING_HITS = 10
FLANK_DRAW_LIMIT = 5
FIND_FEATURES = True
FEATURE_SEARCH_DIST = 10000
OMIT_EVES_IN_REFERENCE = False
BEST_HITS_ONLY = True
DO_CLUSTERING = False
MIN_K = 2
MAX_K = 12
KMEANS_TRIALS = 10
MAX_KMEANS_ITERATIONS = 100
PREFERRED_ACCS = {'Aedes anphevirus':               ['gb|MH037149.1|'], 
                  'Australian Anopheles totivirus': ['ref|NC_035674.1|'], 
                  'Culex ohlsrhavirus':             ['ref|NC_035132.1|'], 
                  'Phasi Charoen-like phasivirus':  ['ref|NC_038263.1|'], 
                  'Wuhan Mosquito Virus 6':         ['gb|MF176251.1|']}
MIN_HITS_TO_SHOW_POSITION = 1  # min aedes hits found at a position to show position in insertpositions text file
MIN_HITS_TO_SHOW_VIRUS = 2     # min specimens found with a particular virus hit to show virus in results.tsv

POPULATIONS = {'Angola': ['Angola'], 
               'Argentina': ['Argentina', 'US_U'], 
               'Australia': ['Australia'], 
               'Brazil': ['Brazil'], 
               'French-Polynesia': ['FrenchPolynesia'],
               'Gabon': ['Gabon'], 
               'Mexico': ['Mexico'], 
               'Philippines': ['Philippines'], 
               'South-Africa': ['South_Africa'], 
               'Thailand': ['Thailand'], 
               'USA': ['USA', 'AZ'], 
               'Vietnam': ['Vietnam']}
               
logFile = open(LOGFILE_PATH, 'a')

###############################################################################

def reverseComplement(dna):
	'''Return the reverse complement of dna.'''
	
	dna = dna.upper()
	basecomplement = {'A': 'T', 'C': 'G', 'T': 'A', 'G': 'C'} 
	complement = ''
	for base in reversed(dna):
		complement = complement + basecomplement.get(base, base)
	return complement

###############################################################################

def getSpecimenLabel(specimen):     
    region = 'Unknown'
    for popName in POPULATIONS:
        for pattern in POPULATIONS[popName]:
            if pattern in specimen:
                region = popName
                break
                
    p = re.compile(r'[0-9]+')
    try:
        num = int(p.findall(specimen)[0])
    except:
        num = 0
    
    return region, num

###############################################################################

def cleanACC(acc):
    acc = acc.rstrip().rstrip('|')  # remove trailing whitespace then '|'
    if '|' in acc:
        acc = acc.split('|')[1]     # remove leading text plus '|'
    return acc.strip()              # remove any remaining whitespace, just in case
    
###############################################################################

def complexity(dna, W = 0):
    """https://doi.org/10.1016/S0097-8485(99)00007-8"""
    
    dna = dna.upper()
    dna = dna.replace('-', '')
    
    if W == 0:
        W = int(math.log(len(dna), 4)) + 2
    
    C = 1
    for k in range(1, W + 1):
        kmers = set()
        for i in range(0, len(dna) - k + 1):
            kmers.add(dna[i:i+k])
        C *= (len(kmers) / (4 ** k))
    return C

###############################################################################

def getFamily(accessionID):
    """
    Query NCBI to get family of virus corresponding to an accession number.
    """
    
    Entrez.email = 'havill@denison.edu'
    try:
        handle = Entrez.efetch(db = 'nucleotide', id = accessionID, rettype = 'native', retmode = 'xml')
        result = handle.read()
    except:
        writelog('Exception:', sys.exc_info()[1])
        return ''
        
    root = ET.fromstring(result)
    lineage = root.find('.//OrgName_lineage')

    if lineage is not None:
        lineageParts = lineage.text.split('; ')
        if len(lineageParts) <= 6:    # 1-6
            return lineageParts[-1]
        elif len(lineageParts) <= 8:  # 7-8
            return lineageParts[6]
        else:
            return lineageParts[7]    # 9-
    else:
        return ''
    
def readFamFile():
    try:
        famFile = open(FAMILY_CSV, 'r')
    except FileNotFoundError:
        return {}
    fams = {}
    for line in famFile:
        cols = line.rstrip().split(',')
        fams[cols[0]] = cols[1]
    famFile.close()
    return fams
    
def addFamily(acc, name):
    famFile = open(FAMILY_CSV, 'a')
    famFile.write(acc + ',' + name + '\n')
    famFile.close()
    
###############################################################################

def writelog(message, alsoPrint = False):
    logFile.write(message + '\n')
    if alsoPrint:
        print(message)
        
###############################################################################
        
def getContig(dirName, contigName, start, end = None):
    scaffoldsFilename = str(Path(dirName) / 'results/scaffolds/scaffolds.fasta')
    scaffoldRecords = SeqIO.index(scaffoldsFilename, 'fasta')
    contigRecord = scaffoldRecords[contigName]
    if end is None:
        return str(contigRecord.seq)[start - 1:]
    return str(contigRecord.seq)[start - 1:end]  # account for 1-based indexing
        
def main():
    pass
    
if __name__ == '__main__':
    main()
    
