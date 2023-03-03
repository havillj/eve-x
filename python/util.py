import matplotlib
matplotlib.use('agg')
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
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as pyplot
import numpy as np
import os
import pickle
import pysam
import random
import re
import sys
import time
import urllib.request as web

from config import *

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
        num = int(p.findall(specimen)[-1])  # number assigned to be last number in name
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
    """https://doi.org/10.1016/S0097-8485(99)00007-8
       We assume that dna is long enough to contain all possible k-mers for 
       k=1 to W."""
    
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

def percentOverlap(iv1, iv2):
    """Return the percent of overlap between intervals iv1 and iv2."""
    
    return max(0, min(iv1[1], iv2[1]) - max(iv1[0], iv2[0])) / (iv1[1] - iv1[0])

###############################################################################

def getSeqRecord(fileNameFASTA, seqid):
    try:
        refFile = open(fileNameFASTA, 'r')
    except:
        prefix1 = 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi'
        prefix2 = '?db=nuccore&id='
        suffix = '&rettype=fasta&retmode=text'
        url = prefix1 + prefix2 + seqid + suffix
        try:
            fastaFile = web.urlopen(url)
        except:
            writelog('Error getting reference sequence ' + seqid + ' from NCBI.  Aligned FASTA file not created.')
            return None
        fastaText = fastaFile.read().decode('utf-8')
        fastaFile.close()
        refFile = open(fileNameFASTA, 'w')
        refFile.write(fastaText)
        refFile.close()
        refFile = open(fileNameFASTA, 'r')
        
    refRecord = SeqIO.read(refFile, 'fasta')
    refFile.close()
    
    return refRecord
    
def getGenBank(accessionID):
	"""
	Query NCBI to get a gb file.
	"""
	
	Entrez.email = EMAIL
	try:
		handle = Entrez.efetch(db = 'nucleotide', id = accessionID, rettype = 'gb', retmode = 'text')
	except:
		writelog('Error retrieving GenBank record for ' + accessionID + ': ' + sys.exc_info()[1] + '\nDiagram not created.')
		return False
		
	gbFile = open(GB_DIR + accessionID + '.gb', 'w')
	gbFile.write(handle.read())
	gbFile.close()
	
	return True
    
###############################################################################

def getFamily(accessionID):
    """
    Query NCBI to get family of virus corresponding to an accession number.
    """
    
    Entrez.email = EMAIL
    try:
        handle = Entrez.efetch(db = 'nucleotide', id = accessionID, rettype = 'native', retmode = 'xml')
        result = handle.read()
    except:
        writelog('Exception:', sys.exc_info()[1])
        return ''
        
    result = bytes(result, 'utf-8')
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
        
def writeConfig():
    writelog('\nConfiguration parameters:')
    writelog('ROOT_DIR = ' + ROOT_DIR)
    writelog('RESULTS_DIR = ' + RESULTS_DIR)
    writelog('SPECIMENS_DIR = ' + SPECIMENS_DIR)
    writelog('SPECIMEN_RESULTS_DIR = ' + SPECIMEN_RESULTS_DIR)
    writelog('HOST_DB = ' + HOST_DB)
    writelog('VIRUS_DB = ' + VIRUS_DB)
    for key in config:
        writelog(key + ' = ' + str(config[key]))
    writelog('')
        
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
    