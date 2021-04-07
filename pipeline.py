from Bio import SeqIO, Entrez
from Bio.SeqFeature import SeqFeature, FeatureLocation
from dna_features_viewer import GraphicFeature, GraphicRecord, BiopythonTranslator
from intervaltree import Interval, IntervalTree
from lxml import etree as ET
from multiprocessing import Process
from pathlib import Path
import copy
import math
import matplotlib
import matplotlib.pyplot as pyplot
import numpy as np
import os
import pickle
import pysam
import re
import sys
import time
import urllib.request as web

SPADES_EXEC = '/home/havill/bin/SPAdes-3.14.1-Linux/bin/spades.py'
ROOT_DIR = '/home/havill/data/aegypti/analyzed/'
RESULTS_DIR = '/home/havill/data2/results4/'
FAMILY_CSV = '/home/havill/data2/results/families.csv'
LOGFILE_NAME = 'log.log'
MIN_PIDENT = 80
EVALUE_V = 1e-10
EVALUE_A = 1e-10
ALLOWED_OVERLAP = 20  # bp aedes hits may overlap ends of virus hits or each other
OVERLAP_FRACTION = 0.9  # if this fraction of a hit overlaps with an Aedes hit, it is considered to be present in the reference genome
MAX_FLANK_DISTANCE = 10000
COMPLEXITY_K = 3
COMPLEXITY_CUTOFF = 0.75
#MAX_FLANKING_HITS = 10
FLANK_DRAW_LIMIT = 5

iTrees = None

def cleanACC(acc):
    acc = acc.rstrip().rstrip('|')  # remove trailing whitespace then '|'
    if '|' in acc:
        acc = acc.split('|')[1]     # remove leading text plus '|'
    return acc.strip()              # remove any remaining whitespace, just in case

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

def getUnmappedReads(filenameBAM):
    """Get unmapped reads AND MATES from a BAM file and write them to a FASTA file.
    
       Parameter:
           filenameBAM: absolute path of input BAM file
           
       Return value: absolute path of the output FASTA file
    """
    
    writelog('\nAnalyzing ' + filenameBAM, True)
    
    try:
        index = filenameBAM.index('.sorted')
    except:
        index = filenameBAM.index('.bam')
    newFilenameBAM = filenameBAM[:index] + '_unmapped_with_mates.bam'  # absolute path

    writelog('   Getting unmapped reads...', True)

    if not Path(newFilenameBAM).exists():
        bamfile = pysam.AlignmentFile(filenameBAM, 'rb', threads = 8)
        newBAM = pysam.AlignmentFile(newFilenameBAM, 'wb', template = bamfile, threads = 8)

        count = 0
        count_total = 0
        for read in bamfile:
            count_total += 1
            if read.is_unmapped or read.mate_is_unmapped:
                seq = read.query_sequence
                if seq.count('N') / len(seq) <= 0.5:
                    count += 1
                    newBAM.write(read)
        bamfile.close()
        newBAM.close()

        writelog('   Total reads =', count_total)
        writelog('   Unmapped reads =', count, '({:.2f}%)'.format(100*count/count_total))
    else:
        writelog('   Skipping - ' + newFilenameBAM + ' exists.', True)

    # Convert unmapped reads to FASTA.

    writelog('   Writing unmapped reads to a FASTA file...', True)
    newFilenameFASTA = newFilenameBAM[:-4] + '.fasta'  # absolute path
    if not Path(newFilenameFASTA).exists():
        os.system('samtools fasta ' + newFilenameBAM + ' > ' + newFilenameFASTA)
    else:
        writelog('   Skipping - ' + newFilenameFASTA + ' exists.', True)
        
    return newFilenameFASTA  # absolute path
    
###############################################################################

def writeReadsFASTQ(viralFilenameBAM, virusName):
    """Write paired-end reads from BAM file to 3 FASTQ files.
    
       Parameters:
           viralFilenameBAM: absolute path of BAM file containing viral reads to assemble
           virusName:        string containing the name of the virus
           
       Return value: None
    """
    
    writelog('Writing paired-end reads to FASTQ files...', True)
    
    spadesPath = Path(viralFilenameBAM).parent / ('spades_' + virusName)
    
    if not spadesPath.exists():
        os.system('mkdir ' + str(spadesPath))
        
    if (spadesPath / 'viral1.fastq').exists():
        writelog('   Skipping ... fastq files exist.', True)
        return
    
    os.system('samtools fastq -1 ' + str(spadesPath / 'viral1.fastq') + ' -2 ' + str(spadesPath / 'viral2.fastq') + ' -s ' + str(spadesPath / 'viral_single.fastq') + ' ' + viralFilenameBAM)
    
###############################################################################

def assembleReads(dirName, virusName):
    """Assemble reads in dirName / spades_virusName into scaffolds.
    
       Parameter: 
           dirName:   absolute path of parent of directory containing fastq files to assemble 
           virusName: string containing the name of the virus
            
       Return value: Boolean indicating whether assembly was successful
    """

    writelog('Assembling reads with SPAdes...', True)
    
    spadesPath = Path(dirName) / ('spades_' + virusName)
    if (spadesPath / 'scaffolds.fasta').exists():
        writelog('   Skipping ... scaffolds.fasta exists.', True)
        return True
        
    viral1Name = str(spadesPath / 'viral1.fastq')
    viral2Name = str(spadesPath / 'viral2.fastq')
    viralSName = str(spadesPath / 'viral_single.fastq')
    
    if os.path.getsize(viralSName) > 0:
        result = os.system(SPADES_EXEC + ' --pe1-1 ' + viral1Name + ' --pe1-2 ' + viral2Name + ' --s1 ' + viralSName + ' --careful -o ' + str(spadesPath))
    else:
        os.remove(viralSName)
        result = os.system(SPADES_EXEC + ' --pe1-1 ' + viral1Name + ' --pe1-2 ' + viral2Name + ' --careful -o ' + str(spadesPath))
                
    return (spadesPath / 'scaffolds.fasta').exists()  # assembly succeeded
    
###############################################################################
    
def blastScaffolds(dirName, virusName, virusDB, evalueV, force = False):
    """BLAST scaffolds against viral and aegypti genomes and combine to identify putative viral insertions."""
    
    spadesPath = Path(dirName) / ('spades_' + virusName)
    scaffoldsName = str(spadesPath / 'scaffolds.fasta')
    outVirusCSV = spadesPath / ('blast_scaffolds_' + virusName + '.csv')
#    outAACSV = str(spadesPath / 'blast_scaffolds_aa.csv')
    
    writelog('BLASTing scaffolds against viral database...', True)
    
    if not force and outVirusCSV.exists():
        writelog('   Skipping - ' + str(outVirusCSV) + ' exists.', True)
    else:
        os.system('blastn -query ' + scaffoldsName + ' -db ' + virusDB + ' -num_threads 16 -task blastn -evalue ' + str(evalueV) +
                  ' -max_target_seqs 5000' +
                  ' -outfmt "10 qseqid qstart qend qseq sstart send sseq evalue bitscore sseqid stitle pident"' + 
                  ' -out ' + str(outVirusCSV))

###############################################################################

def readCSV(csvFilename, minLength):
    """Read a BLAST CSV file and return a dictionary of hits.
     
       Parameter:
           csvFilename: absolute path of CSV file containing BLAST results
           
       Return value: a dictionary with contig_name keys and values that are 
                     lists of tuples containing hits in that contig
    """
    
    try:
        csv = open(csvFilename, 'r')
    except FileNotFoundError:
        return None
    hits = {}
    accsFound = {}
    for line in csv:
        line = line.rstrip()
        row = line.split(',')
        contig = row[0]            # qseqid = contig name
#         qstart = row[1]          # start pos in contig
#         qend = row[2]            # end pos in contig
#         qseq = row[3]            # aligned sequence in contig
#        length = len(row[3])       # length of hit in contig
        sstart = int(row[4])             # start pos in viral genome
        send = int(row[5])               # end pos in viral genome
        length = abs(send - sstart) + 1  # length of viral hit
#         sseq = row[6]            # aligned sequence in viral genome
#         evalue = float(row[7])   # e-value of alignment
#         bitscore = float(row[8]) # bit score of alignment
        sseqid = row[9]          # viral genome accession number
        stitle = row[10]         # viral genome title
#         pident = row[11]         # % identity of alignment

#         if length >= minLength:
#             if contig in hits:
#                 if stitle not in accsUsed[contig]:
#                     accsUsed[contig][stitle] = sseqid
#                     hits[contig].append(tuple(row[1:]))
#                 elif accsUsed[contig][stitle] == sseqid:
#                     hits[contig].append(tuple(row[1:]))
#                 # else skip this line
#             else:
#                 hits[contig] = [tuple(row[1:])]
#                 accsUsed[contig] = {stitle: sseqid}

        if length >= minLength:
            if contig in hits:
                if stitle not in accsFound[contig]:
                    accsFound[contig][stitle] = [(sseqid, length)]
                else: 
                    accsFound[contig][stitle].append((sseqid, length))
                hits[contig].append(tuple(row[1:]))
            else:
                hits[contig] = [tuple(row[1:])]
                accsFound[contig] = {stitle: [(sseqid, length)]}
    csv.close()
    
    accsToUse = {}
    for contig in accsFound:
        accsToUse[contig] = {}
        for stitle in accsFound[contig]:
            refAccs = []
            maxLength = 0
            maxLengthAccs = []
            for acc, length in accsFound[contig][stitle]:
                if acc[:3] == 'ref':
                    refAccs.append(acc)
                elif length > maxLength:
                    maxLength = length
                    maxLengthAccs = [acc]
                elif length == maxLength:
                    maxLengthAccs.append(acc)
            if refAccs != []:
                accsToUse[contig][stitle] = refAccs
            else:
                maxLengthAccs.sort()
                accsToUse[contig][stitle] = [maxLengthAccs[0]]  # use acc with longest hit; if tie use first alphabetically for consistency
                
    filteredHits = {}
    for contig in hits:
        filteredHits[contig] = []
        for hit in hits[contig]:
            sseqid = hit[8]
            stitle = hit[9]
            if sseqid in accsToUse[contig][stitle]:
                filteredHits[contig].append(hit)
            else:
                writelog(contig + ': discarded hit "' + stitle.lstrip('|') + '" with duplicate species')

    return filteredHits
    
def readCSVAA(csvFilename, minLength):
    """Read a BLAST CSV file and return a dictionary of hits.
     
       Parameter:
           csvFilename: absolute path of CSV file containing BLAST results
           
       Return value: a dictionary with contig_name keys and values that are 
                     lists of tuples containing hits in that contig
    """
    
    # qseqid qstart qend qseq sseqid sstart send sseq evalue bitscore
    
    try:
        csv = open(csvFilename, 'r')
    except FileNotFoundError:
        return {}
        
    hits = {}
    for line in csv:
        line = line.rstrip()
        row = line.split(',')
        contig = row[0]            # qseqid = contig name
#         qstart = row[1]          # start pos in contig
#         qend = row[2]            # end pos in contig
#         qseq = row[3]            # aligned sequence in contig
        length = len(row[3])       # length of hit in contig
#         sseqid = row[4]          # viral genome accession number
#         sstart = row[5]          # start pos in viral genome
#         send = row[6]            # end pos in viral genome
#         sseq = row[7]            # aligned sequence in viral genome
#         evalue = float(row[8])   # e-value of alignment
#         bitscore = float(row[9]) # bit score of alignment

        if length >= minLength:
            if contig in hits:
                hits[contig].append(tuple(row[1:]))
            else:
                hits[contig] = [tuple(row[1:])]

    return hits

def getHits(dirName, virusName, minLength, evalueA): #, maxFlankingHits):
    """Combine viral and vector BLAST results to find putative insertions.
    
       Parameters:
           dirName:   absolute path of parent of directory containing SPAdes results 
           virusName: name of virus under consideration
           
       Return value: absolute path of XML file containing results
    """
    
    writelog('Combining BLAST results to locate putative insertions...', True)
    
    spadesPath = Path(dirName) / ('spades_' + virusName)
    csvVirusFilename = str(spadesPath / ('blast_scaffolds_' + virusName + '.csv'))
    csvAedesFilename = str(spadesPath / 'blast_scaffolds_aa.csv')
    xmlFilename = str(spadesPath / (Path(dirName).name + '_' + virusName + '_hits.xml'))
    scaffoldsFilename = str(spadesPath / 'scaffolds.fasta')

    virusHits = readCSV(csvVirusFilename, 0)
    if virusHits is None:
        return None
        
    # virusHits[contig] = [hit0, hit1, ...]
    # where each hit = [qstart, qend, qseq, sstart, send, sseq, evalue, bitscore, sseqid, stitle, pident]
        
    QSTART = 0   # start pos in contig
    QEND = 1     # end pos in contig
    QSEQ = 2     # sequence in contig aligned to sequence in viral genome
    SSTART = 3   # start pos in virus genome
    SEND = 4     # end pos in virus genome
    SSEQ = 5     # sequence in virus genome aligned to sequence in contig
    EVALUE = 6
    BITSCORE = 7
    SEQID = 8    # accession number of viral genome hit
    STITLE = 9   # species name of viral genome hit
    PIDENT = 10  # % identities in alignment
    
    A_QSTART = 0   # start pos in contig
    A_QEND = 1     # end pos in contig
    A_QSEQ = 2     # sequence in contig aligned to sequence in aedes genome
    A_SEQID = 3    # accession number of aedes genome hit
    A_SSTART = 4   # start pos in aedes genome
    A_SEND = 5     # end pos in aedes genome
    A_SSEQ = 6     # sequence in aedes genome aligned to sequence in contig
    A_EVALUE = 7
    A_BITSCORE = 8

##### Combine hits in each contig from the same virus
    
    newVirusHits = {}
    maxLenPIdent = {}
    for contig in virusHits:
        newVirusHits[contig] = []
        leftIndex = 0
        rightIndex = 0
        virusHits[contig].sort(key = lambda hit: (hit[SEQID], int(hit[QSTART]))) # (sseqid, qstart)
        virusHit = virusHits[contig][0]   # first hit in this contig
        sseqid = virusHit[SEQID]          # virus accession number
        hitLength = abs(int(virusHit[QEND]) - int(virusHit[QSTART])) + 1
        hitSeq = virusHit[QSEQ]
#        hitLength = len(virusHit[QSEQ])   # length of hit in contig
        pident = float(virusHit[PIDENT])  # max pident among all hits for current virus
        maxHitLength = hitLength
        maxHitLengthIndex = 0
        maxPIdent = pident
        maxPIdentIndex = 0
        
        for index in range(1, len(virusHits[contig])):  # iterate over remaining hits in contig
            virusHit = virusHits[contig][index]
            nextSeqID = virusHit[SEQID]
            if (nextSeqID == sseqid): # and (nextLeft <= right + 1):
                hitLength += (abs(int(virusHit[QEND]) - int(virusHit[QSTART])) + 1)
                hitSeq += virusHit[QSEQ]
#                hitLength += len(virusHit[2])
                pident = max(pident, float(virusHit[PIDENT]))
            else:  # finish this combined hit and start a new one
                if hitLength >= minLength:
                    hitComplexity = complexity(hitSeq, COMPLEXITY_K)  # COMPLEXITY TEST
                    if hitComplexity < COMPLEXITY_CUTOFF:
                        writelog(contig + ': discarded hit "' + virusHit[STITLE].lstrip('|') + '" with complexity ' + str(hitComplexity) + ' < ' + str(COMPLEXITY_CUTOFF))
                    else:
                        newVirusHits[contig].append((leftIndex, rightIndex))
                        if hitLength > maxHitLength:
                            maxHitLength = hitLength
                            maxHitLengthIndex = leftIndex
                        if pident > maxPIdent:
                            maxPIdent = pident
                            maxPIdentIndex = leftIndex
                else:
                    writelog(contig + ': discarded hit "' + virusHit[STITLE].lstrip('|') + '" with length ' + str(hitLength) + ' < ' + str(minLength))

                # start new combined hit
                leftIndex = index
                hitLength = abs(int(virusHit[QEND]) - int(virusHit[QSTART])) + 1
                hitSeq = virusHit[QSEQ]
#                hitLength = len(virusHit[2])
                pident = float(virusHit[PIDENT])
            rightIndex = index
            sseqid = nextSeqID
            
        # finish last combined hit
        if hitLength >= minLength:
            hitComplexity = complexity(hitSeq, COMPLEXITY_K)  # COMPLEXITY TEST
            if hitComplexity < COMPLEXITY_CUTOFF:
                writelog(contig + ': discarded hit "' + virusHit[STITLE].lstrip('|') + '" with complexity ' + str(hitComplexity) + ' < ' + str(COMPLEXITY_CUTOFF))
            else:
                newVirusHits[contig].append((leftIndex, rightIndex))
                if hitLength > maxHitLength:
                    maxHitLength = hitLength
                    maxHitLengthIndex = leftIndex
                if pident > maxPIdent:
                    maxPIdent = pident
                    maxPIdentIndex = leftIndex
        else:
            writelog(contig + ': discarded hit "' + virusHit[STITLE].lstrip('|') + '" with length ' + str(hitLength) + ' < ' + str(minLength))

        maxLenPIdent[contig] = (maxHitLength, maxHitLengthIndex, maxPIdent, maxPIdentIndex)
        
##### Search contigs for aedes hits to locate putative insertions.
    
    # create a new XML tree
    root = ET.Element('root')
    tree = ET.ElementTree(root)

    all_viral_hits = []
    
    scaffoldRecords = SeqIO.index(scaffoldsFilename, 'fasta')
    
    for contig in virusHits:
        # search for aedes hits
        SeqIO.write(scaffoldRecords[contig], str(spadesPath / 'contig_temp.fasta'), 'fasta')
        os.system('blastn -query ' + str(spadesPath / 'contig_temp.fasta') 
                                   + ' -db aegyptidb -num_threads 8 -evalue ' + str(evalueA)
                                   + ' -max_target_seqs 500'
                                   + ' -outfmt "10 qseqid qstart qend qseq sseqid sstart send sseq evalue bitscore"'
                                   + ' -out ' + str(spadesPath / 'blast_contig_temp.csv'))
        aedesHits = readCSVAA(str(spadesPath / 'blast_contig_temp.csv'), 0)

#         if contig in aedesHits:
#             writelog('  ' + str(len(aedesHits[contig])) + ' aa hits')
#         else:
#             writelog('  no aa hits')
#             continue

#         writelog('  ' + str([v[SEQID] for v in virusHits[contig]]))

        all_viral_hits.extend([(v[SEQID], v[STITLE]) for v in virusHits[contig]])  
        
        # get min evalue, max bitscore, max length, max pident for this contig
        
        minEvalueIndex = 0
        for index in range(1, len(virusHits[contig])):
            if float(virusHits[contig][index][EVALUE]) < float(virusHits[contig][minEvalueIndex][EVALUE]):
                minEvalueIndex = index
                
        maxBitscoreIndex = 0
        for index in range(1, len(virusHits[contig])):
            if float(virusHits[contig][index][BITSCORE]) > float(virusHits[contig][maxBitscoreIndex][BITSCORE]):
                maxBitscoreIndex = index
                
        maxHitLength, maxHitLengthIndex, maxPIdent, maxPIdentIndex = maxLenPIdent[contig]
        
        # Create a separate xml node for each virus found in the contig
        
        thisContigCount = 0
        for leftIndex, rightIndex in newVirusHits[contig]:
            thisContigCount += 1
            virus_qstart = int(virusHits[contig][leftIndex][QSTART])
            virus_qend = int(virusHits[contig][rightIndex][QEND])
            
            # tree of intervals (qstart, qend, (sstart, send))
            virusIntervals = IntervalTree([Interval(int(virusHits[contig][i][QSTART]), int(virusHits[contig][i][QEND]), (int(virusHits[contig][i][SSTART]), int(virusHits[contig][i][SEND]))) for i in range(leftIndex, rightIndex + 1)])
            virusIntervalsLength = sum([iv.length() for iv in virusIntervals])
                
#             writelog(contig + ' ' + str(leftIndex) + ', ' + str(rightIndex))   
            
            bestHit = False
            if (maxPIdent >= MIN_PIDENT):
                if maxPIdentIndex == leftIndex:
                    bestHit = True
            else:
                if maxHitLengthIndex == leftIndex:
                    bestHit = True

            if len(newVirusHits[contig]) > 1:
                contigName = contig + '__' + str(thisContigCount)
            else:
                contigName = contig
                
            contigTree = ET.SubElement(root, 'contig', {'name': contigName, 'besthit': str(bestHit)})
            for index in range(leftIndex, rightIndex + 1):                
                virusHit = virusHits[contig][index]
                vElement = ET.SubElement(contigTree, 'virushit', {'seqid': virusHit[SEQID], 'stitle': virusHit[STITLE], 
                                                                  'minevalue': str(index == minEvalueIndex), 
                                                                  'maxbitscore': str(index == maxBitscoreIndex)})
                ET.SubElement(vElement, 'qstart').text = virusHit[QSTART]
                ET.SubElement(vElement, 'qend').text = virusHit[QEND]
                ET.SubElement(vElement, 'qseq').text = virusHit[QSEQ]
                ET.SubElement(vElement, 'sstart').text = virusHit[SSTART]
                ET.SubElement(vElement, 'send').text = virusHit[SEND]
                ET.SubElement(vElement, 'sseq').text = virusHit[SSEQ]
                ET.SubElement(vElement, 'evalue').text = virusHit[EVALUE]
                ET.SubElement(vElement, 'bitscore').text = virusHit[BITSCORE]
                ET.SubElement(vElement, 'pident').text = virusHit[PIDENT]
                
            # partition aedes hits into before and after viral hit
            
            before = []
            after = []
            overlap = []
            virusIntervalsNotInRef = copy.deepcopy(virusIntervals)
            doesOverlap = False
            if contig in aedesHits:
                for index in range(len(aedesHits[contig])):
                    aedesHit = aedesHits[contig][index]
                    aedes_qstart = int(aedesHit[A_QSTART])
                    aedes_qend = int(aedesHit[A_QEND])
                    virusIntervalsNotInRef.chop(aedes_qstart, aedes_qend)  # detect overlap with virus hit
#                    writelog('chopped', aedes_qstart, aedes_qend)
#                     aedes_chr = aedesHit[A_SEQID]
#                     aedes_sstart = int(aedesHit[A_SSTART])
#                     aedes_send = int(aedesHit[A_SEND])
                    if aedes_qstart < virus_qstart and aedes_qend < virus_qend:
                        before.append(index)
                    elif aedes_qstart > virus_qstart and aedes_qend > virus_qend:
                        after.append(index)
                    else:  # as >= vs and ae <= ve or as <= vs and ae >= ve
                        overlap.append(index)  # aedes hit completely overlaps or is contained in virus hit

#                     if aedes_qend < virus_qstart + ALLOWED_OVERLAP:
#                         before.append(index)
#                     elif aedes_qstart > virus_qend - ALLOWED_OVERLAP:
#                         after.append(index)
#                     else:
#                         overlap.append(index)

#                     if aedes_qend < virus_qstart + ALLOWED_OVERLAP:
#                         before.append((index, aedes_sstart, aedes_send, aedes_chr))
#                     elif aedes_qstart > virus_qend - ALLOWED_OVERLAP:
#                         after.append((index, aedes_sstart, aedes_send, aedes_chr))
#                     else:
#                         overlap.append((index, aedes_sstart, aedes_send, aedes_chr))

                doesOverlap = sum([iv.length() for iv in virusIntervalsNotInRef]) <= (1 - OVERLAP_FRACTION) * virusIntervalsLength
                                        
            contigTree.attrib['inreference'] = str(doesOverlap)
            
            # get intervals that overlap with reference
            choppedIntervals = virusIntervalsNotInRef - virusIntervals  # remove intervals not affected by chopping
            originalChoppedIntervals = IntervalTree()                   # get original versions of chopped intervals
            for iv in choppedIntervals:
                for iv2 in virusIntervals:
                    if iv2.contains_interval(iv):
                        originalChoppedIntervals.add(iv2)
            virusIntervalsInRef = originalChoppedIntervals | choppedIntervals  # combine chopped with originals
            virusIntervalsInRef.split_overlaps()                               # and split at the overlaps
            virusIntervalsInRef -= choppedIntervals                            # remove the chopped intervals
            
            # adjust subject (viral genome) positions based on what was chopped from query positions
            virusIntervalsInRef2 = IntervalTree()
            for iv in virusIntervalsInRef:
                for iv2 in originalChoppedIntervals:
                    if iv.data == iv2.data:
                        if iv.data[0] < iv.data[1]:
                            virusIntervalsInRef2.addi(iv[0], iv[1], (iv.data[0] + (iv[0] - iv2[0]), iv.data[1] - (iv2[1] - iv[1])))
                        else:
                            virusIntervalsInRef2.addi(iv[0], iv[1], (iv.data[0] - (iv[0] - iv2[0]), iv.data[1] + (iv2[1] - iv[1])))
                            
#             writelog(contigName, virusHit[SEQID])
#             writelog(virusIntervals)
#             writelog(choppedIntervals)
#             writelog(virusIntervalsInRef)
#             writelog(virusIntervalsInRef2)

            contigTree.attrib['referenceoverlaps'] = str([iv.data for iv in virusIntervalsInRef2])
            
            matchingFlanks = ET.SubElement(contigTree, 'flanks')
            # find pairs of before/after aedes hits that are on the same strand and within 10Kbp of each other
            for b in before:
                beforeHit = aedesHits[contig][b]
                for a in after:
                    afterHit = aedesHits[contig][a]
                    if beforeHit[A_SEQID] == afterHit[A_SEQID]:  # same chr
                        if ((int(beforeHit[A_SSTART]) < int(beforeHit[A_SEND]) < int(afterHit[A_SSTART]) + ALLOWED_OVERLAP < int(afterHit[A_SEND]) + ALLOWED_OVERLAP) or 
                            (int(beforeHit[A_SSTART]) > int(beforeHit[A_SEND]) > int(afterHit[A_SSTART]) - ALLOWED_OVERLAP > int(afterHit[A_SEND]) - ALLOWED_OVERLAP)) \
                        and (abs(int(beforeHit[A_SEND]) - int(afterHit[A_SSTART])) < MAX_FLANK_DISTANCE):
                            ET.SubElement(matchingFlanks, 'match', {'leftid': str(b), 'rightid': str(a)})
                    elif ((int(beforeHit[A_SSTART]) < int(beforeHit[A_SEND]) < int(afterHit[A_SEND]) + ALLOWED_OVERLAP < int(afterHit[A_SSTART]) + ALLOWED_OVERLAP) or 
                            (int(beforeHit[A_SSTART]) > int(beforeHit[A_SEND]) > int(afterHit[A_SEND]) - ALLOWED_OVERLAP > int(afterHit[A_SSTART]) - ALLOWED_OVERLAP)) \
                        and (abs(int(beforeHit[A_SEND]) - int(afterHit[A_SSTART])) < MAX_FLANK_DISTANCE):
                            ET.SubElement(matchingFlanks, 'inversion', {'leftid': str(b), 'rightid': str(a)})

#             for b in before:
#                 for a in after:
#                     if b[3] == a[3]:  # same chr
#                         if ((b[1] < b[2] < a[1] + 20 < a[2] + 20) or (b[1] > b[2] > a[1] + 20 > a[2] + 20)) and (abs(b[2] - a[1]) < 10000):
#                             ET.SubElement(matchingFlanks, 'match', {'leftid': str(b[0]), 'rightid': str(a[0])})
                            
            for indices, name in [(before, 'vectorhitleft'), (after, 'vectorhitright'), (overlap, 'vectorhitoverlap')]:
                for i in indices:
                    aedesHit = aedesHits[contig][i]
                    aedesElement = ET.SubElement(contigTree, name, {'id': str(i), 'seqid': aedesHit[A_SEQID]})
                    ET.SubElement(aedesElement, 'qstart').text = aedesHit[A_QSTART]
                    ET.SubElement(aedesElement, 'qend').text = aedesHit[A_QEND]
                    ET.SubElement(aedesElement, 'qseq').text = aedesHit[A_QSEQ]
                    ET.SubElement(aedesElement, 'sstart').text = aedesHit[A_SSTART]
                    ET.SubElement(aedesElement, 'send').text = aedesHit[A_SEND]
                    ET.SubElement(aedesElement, 'sseq').text = aedesHit[A_SSEQ]
                    ET.SubElement(aedesElement, 'evalue').text = aedesHit[A_EVALUE]
                    ET.SubElement(aedesElement, 'bitscore').text = aedesHit[A_BITSCORE]
                    
#             for b in before:
#                 beforeHit = aedesHits[contig][b[0]]
#                 aedesElement = ET.SubElement(contigTree, 'vectorhitleft', {'id': str(b[0]), 'seqid': beforeHit[3]})
#                 ET.SubElement(aedesElement, 'qstart').text = beforeHit[0]
#                 ET.SubElement(aedesElement, 'qend').text = beforeHit[1]
#                 ET.SubElement(aedesElement, 'qseq').text = beforeHit[2]
#                 ET.SubElement(aedesElement, 'sstart').text = beforeHit[4]
#                 ET.SubElement(aedesElement, 'send').text = beforeHit[5]
#                 ET.SubElement(aedesElement, 'sseq').text = beforeHit[6]
#                 ET.SubElement(aedesElement, 'evalue').text = beforeHit[7]
#                 ET.SubElement(aedesElement, 'bitscore').text = beforeHit[8]
#             if len(after) > 0:
#                 for a in after:
#                     afterHit = aedesHits[contig][a[0]]
#                     aedesElement = ET.SubElement(contigTree, 'vectorhitright', {'id': str(a[0]), 'seqid': afterHit[3]})
#                     ET.SubElement(aedesElement, 'qstart').text = afterHit[0]
#                     ET.SubElement(aedesElement, 'qend').text = afterHit[1]
#                     ET.SubElement(aedesElement, 'qseq').text = afterHit[2]
#                     ET.SubElement(aedesElement, 'sstart').text = afterHit[4]
#                     ET.SubElement(aedesElement, 'send').text = afterHit[5]
#                     ET.SubElement(aedesElement, 'sseq').text = afterHit[6]
#                     ET.SubElement(aedesElement, 'evalue').text = afterHit[7]
#                     ET.SubElement(aedesElement, 'bitscore').text = afterHit[8]
     
    tree.write(xmlFilename, xml_declaration=True, pretty_print=True)
    
#     writelog(set(all_viral_hits))
#     for a in set(all_viral_hits):
#         writelog(a[0], a[1])
    
    return xmlFilename

###############################################################################

def getGFFFeatures(fileName, iTrees):
    """Read repeat features from a GFF3 file and insert into dictionary of interval trees."""
    
    chromNumbers = {'1': 'NC_035107.1', '2': 'NC_035108.1', '3': 'NC_035109.1', 'MT': 'NC_035159.1'}
    strands = {'+': +1, '-': -1, '.': 0}
   
    gff = open(fileName, 'r')
    
    for line in gff:
        if line[0] == '#':
            continue
        
        cols = line.strip().split('\t')
        
        seqid = cols[0]
#        if seqid == 'MT':
#            continue
        if seqid in chromNumbers:
            seqid = chromNumbers[seqid]
        else:
            seqid += '.1'
            
        if seqid not in iTrees:
            iTrees[seqid] = IntervalTree()
            
        type = cols[2]
#         if type in ['chromosome']:
#             continue
            
        attrib = cols[8]
        
        if type == 'supercontig':
            altID = re.findall(r'Alias=' + seqid + r',(.*)', attrib)  # in supercontig record
            if len(altID) > 0:
                altID = altID[0]
                iTrees[altID] = iTrees[seqid]
        elif type in ('gene', 'ncRNA_gene', 'exon', 'repeat_region'):
            start = int(cols[3])
            end = int(cols[4])
            strand = strands[cols[6]]
            f = SeqFeature()
            f.location = FeatureLocation(start, end + 1, strand = strand)
            f.type = type
            findName = re.findall(r'Name=(.*?);.*', attrib)
            if len(findName) > 0:
                f.id = findName[0]
            else:
                findName = re.findall(r'ID=(.*?);.*', attrib)  # try again
                if len(findName) > 0:
                    f.id = findName[0]
                else:
                    f.id = 'Unknown'
            iTrees[seqid][start:end+1] = f  # interval start <= x < end+1
    
    gff.close()
                
def makeIntervalTrees(dirName):
    pickleName = str(Path(dirName) / 'featuretrees.pickle')
    if Path(pickleName).exists():
        writelog('Pickled feature trees found.  Reading...', True)
        picklediTreesDict = open(pickleName, 'rb')
        iTrees = pickle.load(picklediTreesDict)
        picklediTreesDict.close()
    else:
        iTrees = {}
        writelog('Reading base features...', True)
        # https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/002/204/515/GCF_002204515.2_AaegL5.0/
        getGFFFeatures('/home/havill/data/aegypti/gb/Aedes-aegypti-LVP_AGWG_BASEFEATURES_AaegL5.2.gff3', iTrees)
    
        writelog('Reading repeat features...', True)
        getGFFFeatures('/home/havill/data/aegypti/gb/Aedes-aegypti-LVP_AGWG_REPEATFEATURES_AaegL5.gff3', iTrees)
    
        writelog('Writing feature trees to disk...', True)
        picklediTreesDict = open(pickleName, 'wb')
        pickle.dump(iTrees, picklediTreesDict)
        picklediTreesDict.close()
        
    return iTrees
    
###############################################################################

def searchAA(iTrees, seqid, start, end, hitElement, counts):
    """Get features for a given interval and write to xml tree."""
    
#    searchInterval = Interval(start, end+1)
    
    overlapIntervals = iTrees[seqid][start:end]
    allFeatureTypes = []
    for interval in overlapIntervals:
        f = interval.data                       # feature associated with interval
#         if f.location_operator == 'join':       # mRNA with a CompoundLocation
#             introns = getIntrons(f.location)
#             featureType = None
#             for intron in introns:
#                 if intron.contains_interval(searchInterval):
#                     featureType = 'in_intron'                 # wholly contained in an intron
#                     break
#                 elif intron.overlaps(searchInterval):
#                     featureType = 'overlaps_intron'           # crosses an intron/exon boundary
#                     break
#                     
#             if featureType is None:
#                 if start in f.location and end in f.location:
#                     featureType = 'in_exon'                   # wholly contained in an exon
#                 else:
#                     featureType = 'overlaps_end_exon'         # overlaps first or last exon
#         else:
        featureType = f.type
        allFeatureTypes.append(featureType)
                 
        feature = ET.SubElement(hitElement, 'feature', {'type': featureType})
        ET.SubElement(feature, 'location').text = str(f.location)
        ET.SubElement(feature, 'type').text = f.type
        ET.SubElement(feature, 'id').text = f.id
        
    if 'gene' in allFeatureTypes:
        if 'exon' in allFeatureTypes:
            counts['exon'] = counts.get('exon', 0) + 1
        else:
            counts['intron'] = counts.get('intron', 0) + 1
    elif 'ncRNA_gene' in allFeatureTypes:
        if 'exon' in allFeatureTypes:
            counts['exon_ncRNA'] = counts.get('exon_ncRNA', 0) + 1
        else:
            counts['intron_ncRNA'] = counts.get('intron_ncRNA', 0) + 1
    elif 'pseudogene' in allFeatureTypes:
        if 'exon' in allFeatureTypes:
            counts['exon_pseudogene'] = counts.get('exon_pseudogene', 0) + 1
        else:
            counts['intron_pseudogene'] = counts.get('intron_pseudogene', 0) + 1
        
    for featureType in set(allFeatureTypes):    # only count one repeat_region per sequence
        if featureType not in ('gene', 'ncRNA_gene', 'pseudogene', 'exon'):
            counts[featureType] = counts.get(featureType, 0) + 1
    
    if len(overlapIntervals) == 0:
        counts['None'] += 1
            
def findFeatures(fileName, iTrees):
    """Find all features in an xml file and insert them.
    
       Parameters:
           fileName: absolute path of XML file containing hit results
           iTrees:   dict of interval trees containing AA features
           
       Return value: absolute path of XML file augmented with features
    """
    
    writelog('Annotating features...', True)
    
    parser = ET.XMLParser(remove_blank_text=True)  # lxml
    tree = ET.parse(fileName, parser)
    root = tree.getroot()
    for contig in root:
        counts = {'None': 0}
        hitsLeft = contig.findall('vectorhitleft')
        for v in hitsLeft:
            start = int(v.find('sstart').text)
            end = int(v.find('send').text)
            seqid = v.attrib['seqid']
            searchAA(iTrees, seqid, start, end, v, counts)
            
        hitsRight = contig.findall('vectorhitright')
        for v in hitsRight:
            start = int(v.find('sstart').text)
            end = int(v.find('send').text)
            seqid = v.attrib['seqid']
            searchAA(iTrees, seqid, start, end, v, counts)
            
        featureSummary = ET.SubElement(contig, 'featuresummary')  # number of flanking hits with each feature
        for ftype in counts:
            ET.SubElement(featureSummary, ftype).text = str(counts[ftype])
        
    xmlFeaturesFilename = fileName[:-4] + '_features.xml'
    tree.write(xmlFeaturesFilename, xml_declaration=True, pretty_print=True)
    
    return xmlFeaturesFilename
    
###############################################################################
    
def displayXML(fileName, outFileName, displaySeq = True):
    """Write XML results to human-readable text file.
    
       Parameters:
           fileName:    absolute path of XML results file
           outFileName: absolute path of output text file
           displaySeq:  whether to write long version with all sequence alignments
           
       Return value: None
    """
    
    outFile = open(outFileName, 'w')
    
    totalCounts = {}
    totalFeatures = 0
    totalContigs = 0
    
    tree = ET.parse(fileName)
    root = tree.getroot()
    for contig in root:
        totalContigs += 1
        outFile.write('Contig: ' + contig.attrib['name'] + '\n')
        outFile.write('   Viral hit: \n')
        for v in contig.findall('virushit'):
#            if displaySeq:
            outFile.write('      contig         {0:>5} - {1:<5} {2}\n'.format(v.find('qstart').text, v.find('qend').text, v.find('qseq').text))
            outFile.write('      {3:<14} {0:>5} - {1:<5} {2}\n\n'.format(v.find('sstart').text, v.find('send').text, v.find('sseq').text, v.attrib['seqid']))
#             else:
#                 outFile.write('      contig         {0:>5} - {1:<5}\n'.format(v.find('qstart').text, v.find('qend').text, v.find('qseq').text))
#                 outFile.write('      {3:<14} {0:>5} - {1:<5}\n\n'.format(v.find('sstart').text, v.find('send').text, v.find('sseq').text, v.attrib['seqid']))

        if displaySeq:
            hitsLeft = contig.findall('vectorhitleft')
            hitsRight = contig.findall('vectorhitright')
            hitsOverlap = contig.findall('vectorhitoverlap')
            for hits, locString in ((hitsLeft, 'left'), (hitsRight, 'right'), (hitsOverlap, 'overlapping')):
                if len(hits) > 0:
                    hitsDict = {}  # consolidate query hit: [subject hits]
                    for v in hits:
                        pos = (v.find('qstart').text, v.find('qend').text, v.find('qseq').text)
                        features = v.findall('feature')
                        featureString = '| '
                        for f in features:
                            featureString += f.attrib['type'] + '; ' + f.find('id').text + '; ' + f.find('location').text + ' | '
                        if pos not in hitsDict:
                            hitsDict[pos] = [(v.find('sstart').text, v.find('send').text, v.find('sseq').text, v.attrib['seqid'], featureString)]
                        else:
                            hitsDict[pos].append((v.find('sstart').text, v.find('send').text, v.find('sseq').text, v.attrib['seqid'], featureString))
                    outFile.write('   Vector hits ' + locString + ': \n')
                    posSort = list(hitsDict.keys())
                    posSort.sort(key = lambda pos: (int(pos[0]), int(pos[1])))
                    for q in posSort:
                        outFile.write('      contig         {0:>9} - {1:<9} {2}\n'.format(*q))
                        for s in hitsDict[q]:
                            outFile.write('      {3:<14} {0:>9} - {1:<9} {2} {4}\n'.format(*s))
                        outFile.write('\n')

#             hitsLeft = contig.findall('vectorhitleft')
#             if len(hitsLeft) > 0:
#                 hitsLeftDict = {}  # consolidate query hit: [subject hits]
#                 for v in hitsLeft:
#                     pos = (v.find('qstart').text, v.find('qend').text, v.find('qseq').text)
#                     features = v.findall('feature')
#                     featureString = '| '
#                     for f in features:
# #                         if displaySeq:
#                         featureString += f.attrib['type'] + '; ' + f.find('id').text + '; ' + f.find('location').text + ' | '
# #                         else:
# #                             featureString += f.attrib['type'] + '; ' + f.find('id').text + ' | '
#                     if pos not in hitsLeftDict:
#                         hitsLeftDict[pos] = [(v.find('sstart').text, v.find('send').text, v.find('sseq').text, v.attrib['seqid'], featureString)]
#                     else:
#                         hitsLeftDict[pos].append((v.find('sstart').text, v.find('send').text, v.find('sseq').text, v.attrib['seqid'], featureString))
# #                 if displaySeq:
#                 outFile.write('   Vector hits left: \n')
#                 posSort = list(hitsLeftDict.keys())
#                 posSort.sort(key = lambda pos: (int(pos[0]), int(pos[1])))
#                 for q in posSort:
# #                     if displaySeq:
#                     outFile.write('      contig         {0:>9} - {1:<9} {2}\n'.format(*q))
#                     for s in hitsLeftDict[q]:
#                         outFile.write('      {3:<14} {0:>9} - {1:<9} {2} {4}\n'.format(*s))
#                     outFile.write('\n')
# #                     else:
# #                         outFile.write('      contig         {0:>9} - {1:<9}\n'.format(*q))
# #                         for s in hitsLeftDict[q]:
# #                             outFile.write('      {3:<14} {0:>9} - {1:<9} {4}\n'.format(*s))
# #                         outFile.write('\n')
#         
#             hitsRight = contig.findall('vectorhitright')
#             if len(hitsRight) > 0:
#                 hitsRightDict = {}  # consolidate query hit: [subject hits]
#                 for v in hitsRight:
#                     pos = (v.find('qstart').text, v.find('qend').text, v.find('qseq').text)
#                     features = v.findall('feature')
#                     featureString = '| '
#                     for f in features:
# #                         if displaySeq:
#                         featureString += f.attrib['type'] + '; ' + f.find('id').text + '; ' + f.find('location').text + ' | '
# #                         else:
# #                             featureString += f.attrib['type'] + '; ' + f.find('id').text + ' | '
#                     if pos not in hitsRightDict:
#                         hitsRightDict[pos] = [(v.find('sstart').text, v.find('send').text, v.find('sseq').text, v.attrib['seqid'], featureString)]
#                     else:
#                         hitsRightDict[pos].append((v.find('sstart').text, v.find('send').text, v.find('sseq').text, v.attrib['seqid'], featureString))
# #                if displaySeq:
#                 outFile.write('   Vector hits right: \n')
#                 posSort = list(hitsRightDict.keys())
#                 posSort.sort(key = lambda pos: (int(pos[0]), int(pos[1])))
#                 for q in posSort:
# #                     if displaySeq:
#                     outFile.write('      contig         {0:>9} - {1:<9} {2}\n'.format(*q))
#                     for s in hitsRightDict[q]:
#                         outFile.write('      {3:<14} {0:>9} - {1:<9} {2} {4}\n'.format(*s))
#                     outFile.write('\n')
# #                     else:
# #                         outFile.write('      contig         {0:>9} - {1:<9}\n'.format(*q))
# #                         for s in hitsRightDict[q]:
# #                             outFile.write('      {3:<14} {0:>9} - {1:<9} {4}\n'.format(*s))
# #                         outFile.write('\n')

            flanks = contig.find('flanks')
            matches = flanks.findall('match')
            if len(matches) > 0:
                outFile.write('   Matching flanking regions: \n')
                counter = 0
                for match in matches:
                    outFile.write('      *** HIT ' + str(counter) + ' ***\n')
                    for v in hitsLeft:
                        if v.attrib['id'] == match.attrib['leftid']:
                            features = v.findall('feature')
                            featureString = '| '
                            for f in features:
#                                 if displaySeq:
                                featureString += f.attrib['type'] + '; ' + f.find('id').text + '; ' + f.find('location').text + ' | '
#                                 else:
#                                     featureString += f.attrib['type'] + '; ' + f.find('id').text + ' | '
#                             if displaySeq:
                            outFile.write('         contig         {0:>9} - {1:<9} {2}\n'.format(v.find('qstart').text, v.find('qend').text, v.find('qseq').text))
                            outFile.write('         {3:<14} {0:>9} - {1:<9} {2} {4}\n\n'.format(v.find('sstart').text, v.find('send').text, v.find('sseq').text, v.attrib['seqid'], featureString))
#                             else:
#                                 outFile.write('         contig         {0:>9} - {1:<9}\n'.format(v.find('qstart').text, v.find('qend').text, v.find('qseq').text))
#                                 outFile.write('         {3:<14} {0:>9} - {1:<9} {4}\n\n'.format(v.find('sstart').text, v.find('send').text, v.find('sseq').text, v.attrib['seqid'], featureString))
                            break
                        
                    for v in contig.findall('virushit'):
#                         if displaySeq:
                        outFile.write('         contig         {0:>9} - {1:<9} {2}\n'.format(v.find('qstart').text, v.find('qend').text, v.find('qseq').text))
                        outFile.write('         {3:<14} {0:>9} - {1:<9} {2}\n\n'.format(v.find('sstart').text, v.find('send').text, v.find('sseq').text, v.attrib['seqid']))
#                         else:
#                             outFile.write('         contig         {0:>9} - {1:<9}\n'.format(v.find('qstart').text, v.find('qend').text, v.find('qseq').text))
#                             outFile.write('         {3:<14} {0:>9} - {1:<9}\n\n'.format(v.find('sstart').text, v.find('send').text, v.find('sseq').text, v.attrib['seqid']))

                    for v in hitsRight:
                        if v.attrib['id'] == match.attrib['rightid']:
                            features = v.findall('feature')
                            featureString = '| '
                            for f in features:
#                                 if displaySeq:
                                featureString += f.attrib['type'] + '; ' + f.find('id').text + '; ' + f.find('location').text + ' | '
#                                 else:
#                                     featureString += f.attrib['type'] + '; ' + f.find('id').text + ' | '

#                             if displaySeq:
                            outFile.write('         contig         {0:>9} - {1:<9} {2}\n'.format(v.find('qstart').text, v.find('qend').text, v.find('qseq').text))
                            outFile.write('         {3:<14} {0:>9} - {1:<9} {2} {4}\n\n'.format(v.find('sstart').text, v.find('send').text, v.find('sseq').text, v.attrib['seqid'], featureString))
#                             else:
#                                 outFile.write('         contig         {0:>9} - {1:<9}\n'.format(v.find('qstart').text, v.find('qend').text, v.find('qseq').text))
#                                 outFile.write('         {3:<14} {0:>9} - {1:<9} {4}\n\n'.format(v.find('sstart').text, v.find('send').text, v.find('sseq').text, v.attrib['seqid'], featureString))
                            break
                    counter += 1
        
        featureSummary = contig.find('featuresummary')
        if featureSummary is not None:
            outFile.write('   Feature summary:\n')
            if len(featureSummary) == 1:  # only 'None'
                outFile.write('      None found\n\n')
            else:
                total = 0
                for f in featureSummary:
                    if f.tag != 'None':
                        total += int(f.text)
                noneCount = 0
                for f in featureSummary:
                    if f.tag != 'None':
                        outFile.write('      {0:>18}: {1:>5}  {2:>7.2%}\n'.format(f.tag, f.text, int(f.text) / total))
                    else:
                        noneCount = int(f.text)
                    totalCounts[f.tag] = totalCounts.get(f.tag, 0) + int(f.text)
                outFile.write('      Sequences with no features found: ' + str(noneCount) + '\n')
                totalFeatures += total
                outFile.write('\n')
                
    outFile.write('Overall summary:\n')
    outFile.write('   Total viral hits: ' + str(totalContigs) + '\n\n')
    outFile.write('   Feature counts:\n')
    for f in totalCounts:
        if f == 'None':
            outFile.write('      {0:>18}: {1:>7}\n'.format(f, totalCounts[f]))
        else:
            outFile.write('      {0:>18}: {1:>7}  {2:>7.2%}\n'.format(f, totalCounts[f], totalCounts[f] / totalFeatures))
    outFile.close()
    
def drawXML(fileName, outputDir, flankLimit, bestHitsOnly = True):
    """Write XML results to human-readable text file.
    
       Parameters:
           fileName:    absolute path of XML results file
           outFileName: absolute path of output text file
           displaySeq:  whether to write long version with all sequence alignments
           
       Return value: None
    """
    
    chromNames = {'NC_035107.1': 'Chr1', 'NC_035108.1': 'Chr2', 'NC_035109.1': 'Chr3'}
    
    if not Path(outputDir).exists():
        os.mkdir(outputDir)
    else:
        os.system('rm ' + outputDir + '*/*.png')
        
    fams = readFamFile()
    
    totalCounts = {}
    totalFeatures = 0
    totalContigs = 0
    
    contigFeaturesV = {}
    contigFeaturesA = {}
    contigLengths = {}
    families = {}
    
    tree = ET.parse(fileName)
    root = tree.getroot()
    for contig in root:
        family = None
        if bestHitsOnly and (contig.attrib['besthit'] != 'True'):
            continue
            
        p = re.compile(r'(.*)_length')
        m = p.search(contig.attrib['name'])
        c = m.group(1)                       # short contig name (NODE_xxx)
        if c not in contigFeaturesV:
            contigFeaturesV[c] = []
            contigFeaturesA[c] = []
            
#        writelog(contig.attrib['name'])
        p = re.compile(r'length_(.*)_cov')
        m = p.search(contig.attrib['name'])
        contigLength = int(m.group(1))
        contigLengths[c] = contigLength
        features = []
        totalContigs += 1
        for v in contig.findall('virushit'):
            qstart = int(v.find('qstart').text)
            qend = int(v.find('qend').text)
            sstart = int(v.find('sstart').text)
            send = int(v.find('send').text)
            length = abs(sstart - send) + 1
            if send > sstart:
                strandStr = '+'
            else:
                strandStr = '-'
            evalue = float(v.find('evalue').text)
            bitscore = float(v.find('bitscore').text)
            pident = float(v.find('pident').text)
            stitle = v.attrib['stitle']
            seqid = cleanACC(v.attrib['seqid'])
            
            if family is None:               # set family with first viral hit; overwrite later if find best hit
                if seqid not in fams:
                    fam = getFamily(seqid)
                    fams[seqid] = fam
                    addFamily(seqid, fam)
                family = fams[seqid]
                
#            if v.attrib['minevalue'] == 'True':
#            if v.attrib['maxbitscore'] == 'True':
            stitle += ' (' + str(seqid) + ')'
            if not bestHitsOnly and (contig.attrib['besthit'] == 'True'):
                stitle = '**' + stitle.lstrip('|') + '**'
                
                if seqid not in fams:
                    fam = getFamily(seqid)
                    fams[seqid] = fam
                    addFamily(seqid, fam)
                family = fams[seqid]
                
            if qend > qstart:
                strand = 1
            else:
                strand = -1
            features.append(GraphicFeature(start=qstart, end=qend, strand=strand, 
                                           color='#ff0000', fontdict = {'size': 9},
                                           label = stitle.lstrip('|') + ' ' + str(sstart) + '-' + str(send) + ' (' + strandStr + '; ' + str(length) + ' bp; ' + str(pident) + '%)'))  #; ' + str(evalue) + ')'))
            contigFeaturesV[c].append(copy.deepcopy(features[-1]))
        families[c] = family
        hitsLeft = contig.findall('vectorhitleft')
        hitsRight = contig.findall('vectorhitright')
        hitsOverlap = contig.findall('vectorhitoverlap')
        flanks = contig.find('flanks')
        matches = flanks.findall('match')
        hitsDone = []
        if len(matches) > 0:
            matchCount = 1
            for match in matches[:flankLimit]:
                for hits, attribName in [(hitsLeft, 'leftid'), (hitsRight, 'rightid')]:
                    for v in hits:
                        if v.attrib['id'] == match.attrib[attribName]:
                            hitsDone.append(v.attrib['id'])
                            qstart = int(v.find('qstart').text)
                            qend = int(v.find('qend').text)
                            if qend > qstart:
                                strand = 1
                            else:
                                strand = -1
                            sstart = int(v.find('sstart').text)
                            send = int(v.find('send').text)
                            if send > sstart:
                                strandStr = '+'
                            else:
                                strandStr = '-'
                            seqid = v.attrib['seqid']
                            if seqid in chromNames:
                                seqid = chromNames[seqid]
                            elif seqid[:8] == 'NW_01873':
                                seqid = 'Scaffold ' + seqid[8:12]
                            features.append(GraphicFeature(start=qstart, end=qend, strand=strand, 
                                                       color='#777777', fontdict = {'size': 8},
                                                       label = 'M' + str(matchCount) + ': ' + seqid + ' ' + str(sstart) + '-' + str(send) + ' (' + strandStr + ')'))
                            break
                matchCount += 1
                
        inversions = flanks.findall('inversion')
        if len(inversions) > 0:
            matchCount = 1
            for match in inversions[:flankLimit]:
                for hits, attribName in [(hitsLeft, 'leftid'), (hitsRight, 'rightid')]:
                    for v in hits:
                        if v.attrib['id'] == match.attrib[attribName]:
                            hitsDone.append(v.attrib['id'])
                            qstart = int(v.find('qstart').text)
                            qend = int(v.find('qend').text)
                            if qend > qstart:
                                strand = 1
                            else:
                                strand = -1
                            sstart = int(v.find('sstart').text)
                            send = int(v.find('send').text)
                            if send > sstart:
                                strandStr = '+'
                            else:
                                strandStr = '-'
                            seqid = v.attrib['seqid']
                            if seqid in chromNames:
                                seqid = chromNames[seqid]
                            elif seqid[:8] == 'NW_01873':
                                seqid = 'Scaffold ' + seqid[8:12]
                            features.append(GraphicFeature(start=qstart, end=qend, strand=strand, 
                                                       color='#333333', fontdict = {'size': 8},
                                                       label = 'I' + str(matchCount) + ': ' + seqid + ' ' + str(sstart) + '-' + str(send) + ' (' + strandStr + ')'))
                            break
                matchCount += 1

        aaCoverage = [0] * contigLength
        for hits in (hitsLeft, hitsRight, hitsOverlap):
            if len(hits) > 0:
                hitsDict = {}  # consolidate query hit: [subject hits]
                for v in hits:
                    pos = (v.find('qstart').text, v.find('qend').text, v.find('qseq').text)
                    for x in range(int(pos[0]), int(pos[1])):
                        aaCoverage[x] += 1
                    if v.attrib['id'] in hitsDone:
                        continue
#                    featureString = '| '
#                     features = v.findall('feature')
#                     for f in features:
#                         featLocation = f.find('location').text
#                         p = re.compile(r'\[(\d*):(\d*)\]')
#                         m = p.search(featLocation)
#                         featStart, featEnd = int(m.group(1)), int(m.group(2))
#                         featID = f.find('id').text
#                         featType = f.attrib['type']
                        
                    if pos not in hitsDict:
                        hitsDict[pos] = [(v.find('sstart').text, v.find('send').text, v.find('sseq').text, v.attrib['seqid'])] #, featureString)]
                    else:
                        hitsDict[pos].append((v.find('sstart').text, v.find('send').text, v.find('sseq').text, v.attrib['seqid'])) #, featureString))
                posSort = list(hitsDict.keys())
                if hits != hitsRight:
                    posSort.sort(key = lambda pos: (int(pos[1]), int(pos[0])), reverse = True)  # closest first
                else:
                    posSort.sort(key = lambda pos: (int(pos[0]), int(pos[1])))  # closest first
                for q in posSort[:max(0, flankLimit-len(matches))]:
                    qstart = int(q[0])
                    qend = int(q[1])
                    if qend > qstart:
                        strand = 1
                    else:
                        strand = -1
    #                for s in hitsDict[q]:
                    s = hitsDict[q][0]
                    sstart = s[0]
                    send = s[1]
                    if send > sstart:
                        strandStr = '+'
                    else:
                        strandStr = '-'
                    seqid = s[3]
#                     if len(hitsDict[q]) > 1:
#                         send += '*'
                    if seqid in chromNames:
                        seqid = chromNames[seqid]
                    elif seqid[:8] == 'NW_01873':
                        seqid = 'Scaffold ' + seqid[8:12]
                    if hits == hitsOverlap:
                        color = '#00ff00'
                    else:
                        color = '#0000ff'
                    features.append(GraphicFeature(start=qstart, end=qend, strand=strand, 
                                               color=color, fontdict = {'size': 8},
                                               label = seqid + ' ' + sstart + '-' + send + ' (' + strandStr + ')'))
                                
        record = GraphicRecord(sequence_length = contigLength, features = features)
        
        fig, (ax1, ax2) = pyplot.subplots(2, 1, sharex = True, figsize = (10, 6), gridspec_kw = {'height_ratios': [5, 1]})
        
        record.plot(max_label_length = 80, ax = ax1, with_ruler = False)
        ax2.fill_between(range(contigLength), aaCoverage, step = 'mid', alpha = 0.3)
#        ax2.scatter(range(contigLength), aaCoverage)
        ax2.set_ylim(bottom = 0)
        ax2.set_ylabel('Aa hits')

        familyDir = str(Path(outputDir) / family)
        if not Path(familyDir).exists():
            os.mkdir(familyDir)
        if contig.attrib['besthit'] == 'True':
            fig.savefig(Path(familyDir) / (contig.attrib['name'] + '_BEST-HIT.png'))
        else:
            fig.savefig(Path(familyDir) / (contig.attrib['name'] + '.png'))
        pyplot.close(fig)
        
    if not bestHitsOnly:
        for c in contigFeaturesV:
            record = GraphicRecord(sequence_length = contigLengths[c], features = contigFeaturesV[c])
        
            ax, _ = record.plot(max_label_length = 80, figure_width = 10)
    #        ax2.fill_between(range(contigLength), aaCoverage, step = 'mid', alpha = 0.3)
    #        ax2.set_ylim(bottom = 0)
    #        ax2.set_ylabel('Aa hits')

            familyDir = str(Path(outputDir) / families[c])
            ax.figure.savefig(Path(familyDir) / (c + '_all-viral-hits.png'))
            pyplot.close(ax.figure)
                    
###############################################################################

def copyResults():
    if Path(RESULTS_DIR).exists():
        answer = ' '
        while answer not in 'yYnN':
            answer = input('Results directory exists.  Overwrite (y/n)? ')
        if answer[0] not in 'yY':
            writelog('OK, quitting.', True)
            return False
    else:
        os.system('mkdir ' + RESULTS_DIR)
        
    SPECIMENS_DIR = RESULTS_DIR + 'specimens/'
    if not Path(SPECIMENS_DIR).exists():
        os.system('mkdir ' + SPECIMENS_DIR)
        
    dirList = [dir for dir in Path(ROOT_DIR).iterdir() if dir.is_dir() and 'combined' not in dir.name]
    for dir in dirList:
        writelog('Copying results from ' + dir.name, , True)
        if not Path(SPECIMENS_DIR + dir.name).exists():
            os.system('mkdir ' + SPECIMENS_DIR + dir.name)
            os.system('mkdir ' + SPECIMENS_DIR + dir.name + '/xml')
            os.system('mkdir ' + SPECIMENS_DIR + dir.name + '/txt')
#            os.system('mkdir ' + SPECIMENS_DIR + dir.name + '/sequences')
            os.system('mkdir ' + SPECIMENS_DIR + dir.name + '/scaffolds')
            os.system('mkdir ' + SPECIMENS_DIR + dir.name + '/diagrams')
        
        os.system('cp ' + str(dir / 'spades_all') + '/*_all_hits*.xml ' + SPECIMENS_DIR + dir.name + '/xml')
        os.system('cp ' + str(dir / 'spades_all') + '/*_all_hits*.txt ' + SPECIMENS_DIR + dir.name + '/txt')
        os.system('cp ' + str(dir / 'spades_all') + '/scaffolds.fasta ' + SPECIMENS_DIR + dir.name + '/scaffolds')
        os.system('cp ' + str(dir / 'spades_all') + '/blast_scaffolds_all.csv ' + SPECIMENS_DIR + dir.name + '/scaffolds')
        os.system('cp -r ' + str(dir / 'spades_all') + '/png/* ' + SPECIMENS_DIR + dir.name + '/diagrams')
#        os.system('cp ' + ROOT + '/../results/HITS_all/' + dir.name + '_all_hits.fasta ' + SPECIMENS_DIR + dir.name + '/sequences')
        
    return True

###############################################################################

def reverseComplement(dna):
	'''Return the reverse complement of dna.'''
	
	dna = dna.upper()
	basecomplement = {'A': 'T', 'C': 'G', 'T': 'A', 'G': 'C'} 
	complement = ''
	for base in reversed(dna):
		complement = complement + basecomplement.get(base, base)
	return complement

def consolidateAll(virusName, bestHitsOnly = True, omitInReference = True):
    """Consolidate all results for a particular virus.
       Write a tab-delimited table summarizing hit regions for all specimens.
       Write a FASTA file containing specimen EV regions aligned to viral reference genome
    
       Parameters:
           virusName:   name of the virus
           virusLength: length of the virus reference genome
           
       Return value: None
    """
    
#    cp */spades_all/*all_hits.xml ../results/HITS_all
    
    writelog('Consolidating results in ' + RESULTS_DIR + '...', True)
    
#     dirList = [dir for dir in Path(ROOT).iterdir() if dir.is_dir() and 'combined' not in dir.name]
#     for dir in dirList:
#         spadesPath = dir / ('spades_' + virusName)
#         xmlFilename = spadesPath / (Path(dir).name + '_' + virusName + '_hits.xml')
#         if xmlFilename.exists():
    
    viralHits1 = {}
    viralHits2 = {}
    viralHits0 = {}   # with AA overlap
    viralHits = {}  # all together
    allHits = []
    viralSeqs = {}
    refSeqs = {}
    
    xmlFiles = []
    specimensPath = Path(RESULTS_DIR + 'specimens/').resolve()
    for specimenDir in specimensPath.iterdir():
        for file in (specimenDir / 'xml').iterdir():
            if ('_hits_features.xml' in file.name) or (('_hits.xml' in file.name) and not Path(str(file)[:-4] + '_features.xml').exists()):
                xmlFiles.append(file.resolve())
    
    for file in xmlFiles:
        specimen = file.name.split('_' + virusName)[0]
        specimenSequencesDir = RESULTS_DIR + 'specimens/' + specimen + '/sequences/'
        if not Path(specimenSequencesDir).exists():
            os.system('mkdir ' + specimenSequencesDir)
        virus_fasta = open(specimenSequencesDir + specimen + '_' + virusName + '_hits.fasta', 'w')
        virus_fasta_percontig = open(specimenSequencesDir + specimen + '_' + virusName + '_hits2.fasta', 'w')
#             viralSeqs[specimen] = {}
#             refSeqs[specimen] = {}

        # Parse XML tree.
        
        tree = ET.parse(str(file))
        root = tree.getroot()
        for contig in root:   # one per virus per contig
            if bestHitsOnly and (contig.attrib['besthit'] != 'True'):
                continue
                
            hit = {}
            v = contig.find('virushit')
            seqid = v.attrib['seqid']            # should be same for all v
            if '|' in seqid:
                seqid = seqid.split('|')[1]
            stitle = v.attrib['stitle']
            if '|' in stitle:
                stitle = stitle.lstrip('|')
            if (seqid, stitle) not in viralSeqs:
                viralSeqs[(seqid, stitle)] = {}
                refSeqs[(seqid, stitle)] = {}
            if specimen not in viralSeqs[(seqid, stitle)]:
                viralSeqs[(seqid, stitle)][specimen] = {}
                refSeqs[(seqid, stitle)][specimen] = {}
            for v in contig.findall('virushit'):
                hit[(int(v.find('qstart').text), int(v.find('qend').text))] = (int(v.find('sstart').text), int(v.find('send').text))
                viralSeqs[(seqid, stitle)][specimen][(int(v.find('sstart').text), int(v.find('send').text))] = v.find('qseq').text
                refSeqs[(seqid, stitle)][specimen][(int(v.find('sstart').text), int(v.find('send').text))] = v.find('sseq').text
                virus_fasta.write('>' + seqid + '_' + stitle.replace(' ', '-') + '_' + v.find('sstart').text + '-' + v.find('send').text + '_(host)\n')
                virus_fasta.write(v.find('qseq').text + '\n')
                virus_fasta.write('>' + seqid + '_' + stitle.replace(' ', '-') + '_' + v.find('sstart').text + '-' + v.find('send').text + '_(reference)\n')
                virus_fasta.write(v.find('sseq').text + '\n')
            
            qpos = list(hit.keys())
            qpos.sort()
            combinedSeq = ''
            sposString = ''
            revSeq = hit[qpos[0]][0] > hit[qpos[0]][1]  # reverse complement if first hit is reversed
            for qse in qpos:
                combinedSeq += viralSeqs[(seqid, stitle)][specimen][hit[qse]].replace('-', '')  # remove gaps
                if revSeq:
                    sposString = str(hit[qse][1]) + '-' + str(hit[qse][0]) + ',' + sposString
                else:
                    sposString += str(hit[qse][0]) + '-' + str(hit[qse][1]) + ','
            virus_fasta_percontig.write('>' + seqid + '_' + stitle.replace(' ', '-') + '_' + sposString[:-1] + '\n')
            if revSeq:
                combinedSeq = reverseComplement(combinedSeq)
            virus_fasta_percontig.write(combinedSeq + '\n') 

#                 if hit[0] > hit[-1]:
#                     hit = hit[::-1]
#                hit = tuple(hit)

            numVectorHits = len(contig.findall('vectorhitleft') + contig.findall('vectorhitright'))
#            numOverlapVectorHits = len(contig.findall('vectorhitoverlap'))
            features = []
            featureSummary = contig.find('featuresummary')
            if featureSummary is not None:
#                     total = 0
#                     for f in featureSummary:
#                         if f.tag != 'None':
#                             total += int(f.text)
                for f in featureSummary:
                    if f.tag != 'None':
                        features.append((f.tag, int(f.text) / numVectorHits))
                        
                    
            if specimen not in viralHits:
                viralHits[specimen] = []
            viralHits[specimen].append((seqid, stitle, hit, features))
            
            flanks = contig.find('flanks')
            matches = flanks.findall('match')
#                numVectorHits = len(contig.findall('vectorhitleft') + contig.findall('vectorhitright'))
            if len(matches) > 0:
                if specimen not in viralHits2:
                    viralHits2[specimen] = []
                viralHits2[specimen].append((seqid, stitle, hit, features))
            elif numVectorHits > 0:
                if specimen not in viralHits1:
                    viralHits1[specimen] = []
                viralHits1[specimen].append((seqid, stitle, hit, features))
                
#            if numOverlapVectorHits > 0:
            if contig.attrib['inreference'] == 'True':
                if specimen not in viralHits0:
                    viralHits0[specimen] = []
                viralHits0[specimen].append((seqid, stitle, hit, features))
                
            if (seqid, stitle, hit) not in allHits:
                allHits.append((seqid, stitle, hit, features))
        virus_fasta.close()
        virus_fasta_percontig.close()

    VIRUSES_DIR = RESULTS_DIR + 'viruses/'
    if not Path(VIRUSES_DIR).exists():
        os.system('mkdir ' + VIRUSES_DIR)
    SEQUENCES_DIR = VIRUSES_DIR + 'sequences/'
    if not Path(SEQUENCES_DIR).exists():
        os.system('mkdir ' + SEQUENCES_DIR)
        
    fams = readFamFile()
    for (seqid, stitle) in viralSeqs:
        if seqid not in fams:
            fam = getFamily(seqid)
            fams[seqid] = fam
            addFamily(seqid, fam)
        for specimen in viralSeqs[(seqid, stitle)]:
            newDict = {}
            newDictRef = {}
            for (start, end) in viralSeqs[(seqid, stitle)][specimen]:
                if start > end:
                    newDict[(end, start)] = reverseComplement(viralSeqs[(seqid, stitle)][specimen][(start, end)])
                    newDictRef[(end, start)] = reverseComplement(refSeqs[(seqid, stitle)][specimen][(start, end)])
                else:
                    newDict[(start, end)] = viralSeqs[(seqid, stitle)][specimen][(start, end)]
                    newDictRef[(start, end)] = refSeqs[(seqid, stitle)][specimen][(start, end)]
            viralSeqs[(seqid, stitle)][specimen] = newDict
            refSeqs[(seqid, stitle)][specimen] = newDictRef
        
        refIndices = {}    # specimen seq at each reference index
        for specimen in viralSeqs[(seqid, stitle)]:
            refIndices[specimen] = {}
            for (start, end) in viralSeqs[(seqid, stitle)][specimen]:
                refIndex = start - 2  # -1 for 0-indexing, -1 for initial increment
                seq = viralSeqs[(seqid, stitle)][specimen][(start, end)]
                refSeq = refSeqs[(seqid, stitle)][specimen][(start, end)]
                for index in range(len(seq)):
                    if refSeq[index] != '-':
                        refIndex += 1
                        refIndices[specimen][refIndex] = seq[index]
                    else:
                        refIndices[specimen][refIndex] += seq[index]
            viralSeqs[(seqid, stitle)][specimen] = ''
                    
        referenceFilename = '/home/havill/data/aegypti/genomes/' + seqid + '.fasta'
        try:
            refFile = open(referenceFilename, 'r')
        except:
            prefix1 = 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi'
            prefix2 = '?db=nuccore&id='
            suffix = '&rettype=fasta&retmode=text'
            url = prefix1 + prefix2 + seqid + suffix
            try:
                fastaFile = web.urlopen(url)
            except:
                writelog('Error getting reference sequence ' + seqid + ' from NCBI.  Aligned FASTA file not created.')
                continue
            fastaText = fastaFile.read().decode('utf-8')
            fastaFile.close()
            refFile = open(referenceFilename, 'w')
            refFile.write(fastaText)
            refFile.close()
            refFile = open(referenceFilename, 'r')
            
        refRecord = SeqIO.read(refFile, 'fasta')
        refFile.close()
        referenceSeq = str(refRecord.seq)
        virusLength = len(referenceSeq)
        referenceID = refRecord.id
        assert referenceID == seqid
        viralSeqs[(seqid, stitle)][referenceID] = ''
        for refIndex in range(virusLength):
            maxLength = max([len(refIndices[specimen].get(refIndex, '-')) for specimen in refIndices])
            viralSeqs[(seqid, stitle)][referenceID] += '{0:-<{1}}'.format(referenceSeq[refIndex], maxLength)
            for specimen in refIndices:
                viralSeqs[(seqid, stitle)][specimen] += '{0:-<{1}}'.format(refIndices[specimen].get(refIndex, ''), maxLength)

        FAMILY_DIR = SEQUENCES_DIR + fams[seqid]
        if not Path(FAMILY_DIR).exists():
            os.system('mkdir ' + FAMILY_DIR)
        
        seqOut = open(str(Path(FAMILY_DIR) / ('sequences_' + seqid + '.fasta')), 'w')
        sortedSpecimens = list(viralSeqs[(seqid, stitle)].keys())
        sortedSpecimens.remove(referenceID)
        sortedSpecimens.sort()
        seqOut.write('>' + referenceID + ' ' + stitle + '\n')
        seqOut.write(viralSeqs[(seqid, stitle)][referenceID] + '\n')
        for specimen in sortedSpecimens:
            seqOut.write('>' + specimen + '\n')
            seqOut.write(viralSeqs[(seqid, stitle)][specimen] + '\n')
        seqOut.close()
#
    out = open(str(Path(RESULTS_DIR) / ('results_' + virusName + '.tsv')), 'w')
    
    allViruses = []
    for entry in allHits:
        allViruses.append((entry[1], entry[0]))  # stitle, seqid
    allViruses = list(set(allViruses))
    allViruses.sort()
#    allHits.sort()                        # write (start, end) header
#     hitCounts1 = {}
#     hitCounts2 = {}
    out.write('\t')
    for stitle, seqid in allViruses:
        out.write(seqid + '\t')
    out.write('\n')
    out.write('\t')
    for stitle, seqid in allViruses:
        out.write(stitle + '\t')
    out.write('\n')
    
    fams = readFamFile()
    allFamilies = []
    for stitle, seqid in allViruses:
        if seqid not in fams:
            fam = getFamily(seqid)
            fams[seqid] = fam
            addFamily(seqid, fam)
        allFamilies.append(fams[seqid])
    out.write('\t')
    for fam in allFamilies:
        out.write(fam + '\t')
    out.write('\n')
    
#     for hit in allHits:
#         out.write(str(hit) + '\t')
#         hitCounts1[hit] = 0
#         hitCounts2[hit] = 0
#     out.write('\n')
#     out.write('\t')                       # write length header
#     for hit in allHits:
#         length = sum([abs(hit[i+1] - hit[i]) + 1 for i in range(0, len(hit), 2)])
#         out.write(str(length) + '\t')
#     out.write('\n')
    specimens = list(viralHits.keys())
    specimenLabels = {}
    for specimen in specimens:
        specimenLabels[getSpecimenLabel(specimen)] = specimen
    specimenLabelsList = list(specimenLabels.keys())
    specimenLabelsList.sort()
    for region, num in specimenLabelsList:
        specimen = specimenLabels[(region, num)]
        out.write(region + ' ' + num + '\t')
        for stitle, seqid in allViruses:
            length = 0
#            rangeString = ''
            rangeList = []
            for hit in viralHits[specimen]:  # hit = (seqid, stitle, hitDict) in one contig
                if hit[0] == seqid:
                    hitDict = hit[2]
                    qpos = list(hitDict.keys())
                    qpos.sort()                  # within hit, sort by qstart
                    s = []
                    for p in qpos:
                        length += abs(hitDict[p][1] - hitDict[p][0]) + 1
#                        rangeString += str(hitDict[p][0]) + '-' + str(hitDict[p][1]) + ', '
                        s.extend([hitDict[p][0], hitDict[p][1]])
                    if s[0] > s[-1]:
                        s = s[::-1]
                    if (specimen in viralHits0) and (hit in viralHits0[specimen]):  # overlapping vector hit
                        if omitInReference:
                            continue
                        else:
                            s.append(3)
                    elif (specimen in viralHits2) and (hit in viralHits2[specimen]):
                        s.append(2)
                    elif (specimen in viralHits1) and (hit in viralHits1[specimen]):
                        s.append(1)
                    else:
                        s.append(0)
                        
                    if len(hit) > 3:
                        features = hit[3]
                    else:
                        features = []
                    
                    featureString = ''
                    for f in features:
                        featureString += '{0}: {1:4.2f}, '.format(f[0], f[1])
                    if len(featureString) > 0:
                        featureString = featureString[:-2]
                    s.append(featureString)

                    rangeList.append(s)  # s = (sstart, send, sstart, send, ..., 0/1/2/3, features)
                        
                    
            rangeList.sort()    # among all hits, sort by first sstart
            rangeString = ''
            for s in rangeList:
                for i in range(0, len(s) - 2, 2):
                    rangeString += str(s[i]) + '-' + str(s[i+1]) + ', '
                rangeString = rangeString[:-2] + ('*' * s[-2]) + ' (' + s[-1] + ') | '
            if length > 0:
                out.write(rangeString[:-1] + '| ' + str(length))
            out.write('\t')
                    
            # if (specimen in viralHits2) and (hit in viralHits2[specimen]):
#                 out.write('2\t')
#                 hitCounts2[hit] += 1
#             elif (specimen in viralHits1) and (hit in viralHits1[specimen]):
#                 out.write('1\t')
#                 hitCounts1[hit] += 1
#             else:
#                 out.write('0\t')
        out.write(specimen + '\n')
#     out.write('Total 2\t')
#     for hit in allHits:
#         out.write(str(hitCounts2[hit]) + '\t')
#     out.write('\n')
#     out.write('Total 1\t')
#     for hit in allHits:
#         out.write(str(hitCounts1[hit]) + '\t')
#     out.write('\n')
#     out.write('Total all\t')
#     for hit in allHits:
#         out.write(str(hitCounts1[hit] + hitCounts2[hit]) + '\t')
#     out.write('\n')
    out.close()   

###############################################################################

populations = {'Angola': ['Angola'], 
               'Argentina': ['Argentina', 'US_U'], 
               'Australia': ['Australia'], 
               'Brazil': ['Brazil'], 
               'French Polynesia': ['FrenchPolynesia'],
               'Gabon': ['Gabon'], 
               'Mexico': ['Mexico'], 
               'Philippines': ['Philippines'], 
               'South Africa': ['South_Africa'], 
               'Thailand': ['Thailand'], 
               'USA': ['USA', 'AZ'], 
               'Vietnam': ['Vietnam']}
               
def getSpecimenLabel(specimen):
    region = 'Unknown'
    for popName in populations:
        for pattern in populations[popName]:
            if pattern in specimen:
                region = popName
                break
                
    p = re.compile(r'[0-9]+')
    try:
        num = int(p.findall(specimen)[0])
    except:
        num = 0
    
    return region, num
               
AaegL5_hits = {'NC_001564.2': [9330, 8303],   # Cell fusing agent virus
               'NC_034017.1': [1677, 839, 2269, 3801, 5188, 6416, 7537, 8296, 10083, 8327],  # Xishuangbanna aedes flavivirus
               'NC_038512.1': [2579, 3877],   # Trichoplusia ni TED virus
               'NC_035674.1': [446, 3299],    # Australian Anopheles totivirus
               'NC_035132.1': [6568, 7294],   # Culex ohlsrhavirus
               'NC_028484.1': [6555, 7282],   # Tongilchon ohlsrhavirus
               'NC_040669.1': [6773, 7497],   # Riverside virus 1
               'NC_025384.1': [6789, 6267]}   # Culex tritaeniorhynchus rhabdovirus
               
class Hit:
    def __init__(self, pos = [], overlap = False, primary = True, data = None):
        self._pos = pos
        self._overlap = overlap
        self._primary = primary
        self._data = data
        
    def doesOverlap(self):
        return self._overlap
        
    def isPrimary(self):
        return self._primary
        
    def getData(self):
        return self._data
        
    def getPos(self):
        return self._pos
        
    def __len__(self):
        return len(self._pos)
        
    def __getitem__(self, index):
        return self._pos[index]
        
    def __contains__(self, interval):
        try:
            iterator = iter(interval)
        except TypeError:
            start = end = interval
        else:
            start, end = interval
            if start > end:
                start, end = end, start
                
        for index in range(0, len(self._pos), 2):
            if (self._pos[index] <= self._pos[index + 1]) and ((start >= self._pos[index]) and (end <= self._pos[index + 1])):
                return True
            elif (self._pos[index] > self._pos[index + 1]) and ((start <= self._pos[index]) and (end >= self._pos[index + 1])):
                return True
        return False
        
class Hits:
    def __init__(self):
        self._hits = []   # list of hits
        
    def addHit(self, pos, overlap = False, primary = True, data = None):
        self._hits.append(Hit(pos, overlap, primary, data))
        
    def addHits(self, otherHitsList):
        self._hits.extend(otherHitsList)
        
    def getHits(self):
        return self._hits
        
    def getPrimaryHitsOnly(self):
        primaryHits = []
        for hit in self._hits:
            if hit.isPrimary():
                primaryHits.append(hit)
        return primaryHits
        
    def adjust(self, adjustment):
        for index in range(len(self._hits)):
            for index2 in range(len(self._hits[index]._pos)):
                self._hits[index]._pos[index2] += adjustment
                
    def __contains__(self, interval):
        for hit in self._hits:
            if interval in hit:
                return True
        return False
                
    def __len__(self):
        return len(self._hits)

def getGenBank(accessionID):
	"""
	Query NCBI to get a gb file.
	"""
	
	Entrez.email = 'havill@denison.edu'
	try:
		handle = Entrez.efetch(db = 'nucleotide', id = accessionID, rettype = 'gb', retmode = 'text')
	except:
		writelog('Error retrieving GenBank record for ' + accessionID + ': ' + sys.exc_info()[1] + '\nDiagram not created.')
		return False
		
	gbFile = open(GB_DIR + accessionID + '.gb', 'w')
	gbFile.write(handle.read())
	gbFile.close()
	
	return True

def drawVirus(acc, family, hits, allSpecimens, separatePops, isFamily, showAaegL5_hits, showPlot, small, omitInReference):    

#    hits = getHits(acc)

    VIRUSES_DIR = RESULTS_DIR + 'viruses/'
    if not Path(VIRUSES_DIR).exists():
        os.system('mkdir ' + VIRUSES_DIR)
    DIAGRAMS_DIR = VIRUSES_DIR + 'diagrams/'
    if not Path(DIAGRAMS_DIR).exists():
        os.system('mkdir ' + DIAGRAMS_DIR)
    FAMILY_DIR = DIAGRAMS_DIR + family
    if not Path(FAMILY_DIR).exists():
        os.system('mkdir ' + FAMILY_DIR)
    
#     if not Path(outputDir).exists():
#         os.mkdir(outputDir)

    if small:
        width = 3
        thickness = 1
        ref_thickness = 2
        fontsize = 1
        plot_linewidth = 0.25 # 0.75
        draw_line = False
        with_ruler = False
    #   ax.xaxis.set_tick_params(width=5)
    else:
        width = 10
        thickness = 10
        fontsize = 8
        ref_thickness = 10
        plot_linewidth = 0.75
        draw_line = True
        with_ruler = True

    
    if not Path(GB_DIR + acc + '.gb').exists():
        if not getGenBank(acc):
            return
    virusRecord = MyCustomTranslator().translate_record(GB_DIR + acc + '.gb')  # GraphicRecord
#    virusRecord.plots_indexing = 'genbank'  # start at 1
    codingStart = 1
    codingEnd = virusRecord.sequence_length
    for f in virusRecord.features:
        f.thickness = ref_thickness
        f.linewidth = 0
        f.fontdict = {'size': fontsize}
        if 'UTR' in f.label:
            if '5' in f.label:
                codingStart = f.end + 1
            elif '3' in f.label:
                codingEnd = f.start   # starts are actual start - 1 for some reason
                
    virusRecord2 = SeqIO.read('/home/havill/data/aegypti/gb/' + acc + '.gb', 'gb')  # SeqRecord
    if isFamily:
        virusName = 'Family ' + family + ' (' + virusRecord2.id + ': ' + virusRecord2.description[:80]
        if len(virusRecord2.description) > 80:
            virusName += '...'
        virusName += ')'
        writelog('Creating diagram' + 's' * separatePops + ' for family ' + family + ' (using representative ' + virusRecord2.id + ')...', True)
    else:
        virusName = virusRecord2.id + ': ' + virusRecord2.description[:80]
        if len(virusRecord2.description) > 80:
            virusName += '...'
        virusName += ' (' + family + ')'
        writelog('Creating diagram' + 's' * separatePops + ' for ' + virusName + '...', True)
                
    features = {}
    for specimen in allSpecimens:
        features[specimen] = []
        if specimen in hits:
            g = 1.0
            for hit in hits[specimen].getHits():
                if hit.isPrimary():
                    transparency = 0.75  # so overlap can be seen
                else:
                    transparency = 0.25
                if len(hit) == 2:
                    color = (0, 0, 1) + (transparency,) # blue
                else:
                    color = (g, 0, 0) + (transparency,) # red
                    g = max(g - 0.25, 0)
#                 if hit.doesOverlap():
#                     lc = (1, 1, 1) + (transparency,)    # white
#                 else:
                lc = (0, 0, 0, 0)
                for index in range(0, len(hit), 2):
                    start = hit[index]
                    end = hit[index + 1]
                    features[specimen].append(GraphicFeature(start=start, end=end, strand=1, color=color, label = None, thickness = thickness, linewidth = 0, linecolor = lc))
                if hit.getData() != []:  # add reference overlap regions
                    lc = (1, 1, 1) + (transparency,)    # white
                    for start, end in hit.getData():
                        features[specimen].append(GraphicFeature(start=start, end=end, strand=1, color=color, label = None, thickness = thickness, linewidth = 1, linecolor = lc))

            
#     for specimen in hits2:
#         g = 1.0
#         for hit in hits2[specimen].getHits():
#             if len(hit) == 2:
#                 color = (0, 0, 1, 0.25)   # translucent blue
#             else:
#                 color = (g, 0, 0, 0.25)   # translucent red
#                 g = g - 0.25
#             if hit.doesOverlap():
#                 lc = (1, 1, 1, 0.25)      # white
#             else:
#                 lc = (0, 0, 0, 0)
#             for index in range(0, len(hit), 2):
#                 start = hit[index]
#                 end = hit[index + 1]
#                 features[specimen].append(GraphicFeature(start=start, end=end, strand=1, color=color, label = None, thickness = thickness, linewidth = 0, linecolor = lc))
        
    newFeatures = {}
    for specimen in features:
        region, num = getSpecimenLabel(specimen)
        newFeatures[(region, num)] = features[specimen]
        
    pops = ['all']
    if separatePops:
        pops.extend(list(populations.keys()))
        
    for pop in pops:
        if pop == 'all':
            specimens = list(newFeatures.keys())
        else:
            specimens = [specimen for specimen in newFeatures if specimen[0] == pop]
        specimens.sort()
        
        height = (6 * thickness / 225) * (len(specimens) + 3)
        
        numAxes = len(specimens) + 2
        if showPlot:
            numAxes +=1
        if showAaegL5_hits and (acc in AaegL5_hits):
            numAxes += 1
            
        height = (6 * thickness / 225) * numAxes
        
        fig, ax = pyplot.subplots(numAxes, 1, sharex = False, sharey = False, figsize = (width, height), gridspec_kw = {'height_ratios': [4] + [1] * (numAxes - 1)})
        fig.subplots_adjust(top = 0.99, bottom = 0.01, left = 0.18, right = 0.95)
        fig.suptitle(virusName, fontsize = fontsize, fontweight = 'bold', y = 0.995)
        
        virusRecord.box_linewidth = 0
        virusRecord.plot(ax = ax[0], figure_width=width, draw_line = draw_line, with_ruler = with_ruler, max_line_length = 30, max_label_length = 30) # prevent multiline labels
        x0, y0, w, h = ax[0].get_position().bounds
        ax[0].set_position([x0, y0+h*0.2, w, h])
        ax[0].tick_params(labelsize = fontsize)
        xlims = ax[0].get_xlim()
        
        AaegL5_features = []
        if acc in AaegL5_hits:
            for index in range(0, len(AaegL5_hits[acc]), 2):
                start = AaegL5_hits[acc][index]
                end = AaegL5_hits[acc][index + 1]
                AaegL5_features.append(GraphicFeature(start=start, end=end, strand=0, color='lightgray', label = None, thickness = thickness, linewidth = 0))
    
        x = range(virusRecord.sequence_length)
        y = [0] * virusRecord.sequence_length
        
        for specimen in specimens:
            for feature in newFeatures[specimen]:
                start = min(max(feature.start, 1), virusRecord.sequence_length - 1)
                end = min(feature.end, virusRecord.sequence_length - 1)
                for i in range(min(start, end), max(start, end) + 1):
                y[i-1] += 1
                    
        for feature in AaegL5_features:
            start = feature.start
            end = feature.end
            for i in range(min(start, end), max(start, end) + 1):
                y[i-1] += 1
                
        axisIndex = 1
        y = [y[i] / len(specimens) for i in range(len(y))]
        if showPlot:
            ax[1].plot(x, y, color = 'blue', linewidth = plot_linewidth)
            ax[1].set_xlim(*xlims)
            ax[1].set_ylim(0, 2)
            ax[1].axis('off')
            ax[1].set_frame_on(False)
            axisIndex += 1
    
        y = [y[i] > 0 for i in range(len(y))]
        ax[axisIndex].bar(x, y, color = 'blue')
        ax[axisIndex].set_xlim(*xlims)
        ax[axisIndex].set_ylim(0, 2)
        ax[axisIndex].axis('off')
        ax[axisIndex].set_frame_on(False)
        coverage = (sum(y[codingStart-1:codingEnd]) / (codingEnd - codingStart + 1)) * 100
        ax[axisIndex].text(virusRecord.sequence_length, -0.1, '{0:4.1f}%'.format(coverage) + ' ', horizontalalignment = 'left', fontdict = {'size': fontsize})
    
        axisIndex += 1
        if showAaegL5_hits and (acc in AaegL5_hits):
            record = GraphicRecord(features = AaegL5_features, sequence_length = virusRecord.sequence_length, feature_level_height = 0)
            record.plot(ax = ax[axisIndex], figure_width=width, with_ruler = False, draw_line = draw_line)
            
            start, end = record.span
            plot_start, plot_end = start - 0.8, end - 0.2
            ax[axisIndex].plot([plot_start, plot_end], [0, 0], linewidth = 0.25, zorder=-1000, c="k")
            
            ax[axisIndex].text(0, -0.1, 'AaegL5 ', horizontalalignment = 'right', fontdict = {'size': fontsize})
            axisIndex += 1
       
        for specimen in specimens:
            specimenName = specimen[0] + ' ' + str(specimen[1])
#            writelog(specimenName)
            record = GraphicRecord(features = newFeatures[specimen], sequence_length = virusRecord.sequence_length, feature_level_height = 0)
            record.plot(ax = ax[axisIndex], figure_width=width, with_ruler = False, draw_line = draw_line)
            
            start, end = record.span
            plot_start, plot_end = start - 0.8, end - 0.2
            ax[axisIndex].plot([plot_start, plot_end], [0, 0], linewidth = 0.25, zorder=-1000, c="k")
            
            ax[axisIndex].text(0, -0.1, specimenName + ' ', horizontalalignment = 'right', fontdict = {'size': fontsize})
            
            for p in ax[axisIndex].patches:
                if p.get_edgecolor()[:3] == (1, 1, 1):  # white
                    p.set_hatch('////')
                    
            
            axisIndex += 1
        if isFamily:
            figName = family + '_' + pop + '.pdf'
        else:
            figName = acc + '_' + pop + '.pdf'
        fig.savefig(Path(FAMILY_DIR) / figName)
        pyplot.close(fig)

def getHits(omitInReference = True):
    dir = Path(RESULTS_DIR + 'specimens/')
    subdirs = [d for d in dir.iterdir()]
    files = [d / ('xml/' + d.name + '_all_hits_features.xml') for d in subdirs ]
    files.sort()

    allVirusHits = {}
    allSpecimens = []
    for file in files:
        specimen = file.name.rstrip('_all_hits_features.xml')
        writelog('Getting hits for ' + specimen + '...', True)
        allSpecimens.append(specimen)
        tree = ET.parse(str(file))
        root = tree.getroot()
        for contig in root:
            if 'inreference' in contig.attrib:
                overlap = contig.attrib['inreference'] == 'True'
            else:
                overlap = False
            if 'referenceoverlaps' in contig.attrib:
                referenceOverlaps = eval(contig.attrib['referenceoverlaps'])
            else:
                referenceOverlaps = []
            if omitInReference and overlap:
                continue
            virushits = contig.findall('virushit')
            v = virushits[0]  # first virus hit
            seqid = v.attrib['seqid']
            if '|' in seqid:
                seqid = seqid.split('|')[1]
            if seqid not in allVirusHits:
                allVirusHits[seqid] = {}
            if specimen not in allVirusHits[seqid]:
                allVirusHits[seqid][specimen] = Hits()

            hit = []
            for virushit in virushits:
                start = int(virushit.find('sstart').text)
                end = int(virushit.find('send').text)
                hit.extend([start, end])
            primary = contig.attrib['besthit'] == 'True'  # and (len(virushits) < 5):
            allVirusHits[seqid][specimen].addHit(hit, overlap, primary, referenceOverlaps)
    return allVirusHits, allSpecimens

def drawFamily(families, famACCs, separatePops, showAaegL5_hits, showPlot, small, omitInReference):
    allVirusHits, allSpecimens = getHits(omitInReference)
        
    if famACCs is None:
        fams = readFamFile()
    
    for family, repACC in families:
        if famACCs is None:  # if False, must be only one family in families
            famACCs = []
            for acc in allVirusHits:
                if acc not in fams:
                    fam = getFamily(acc)
                    fams[acc] = fam
                    addFamily(acc, fam)
                if fams[acc] == family:
                    famACCs.append(acc)
        
        startPositions = {}
        for acc in famACCs:
            if not Path(GB_DIR + acc + '.gb').exists():
                if not getGenBank(acc):
                    continue
            virusRecord = SeqIO.read('/home/havill/data/aegypti/gb/' + acc + '.gb', 'gb')  # SeqRecord
            minCDSPos = math.inf
            for feature in virusRecord.features:
                if feature.type == 'CDS':
                    minCDSPos = min(minCDSPos, feature.location.start)
            startPositions[acc] = minCDSPos
    
        famHits = {}
        for acc in famACCs:
            for specimen in allVirusHits[acc]:
                if specimen not in famHits:
                    famHits[specimen] = Hits()
                allVirusHits[acc][specimen].adjust(startPositions[repACC] - startPositions[acc])
                famHits[specimen].addHits(allVirusHits[acc][specimen].getPrimaryHitsOnly())  # only primary hits so one hit per contig
                
        drawVirus(repACC, family, famHits, allSpecimens, separatePops, True, showAaegL5_hits, showPlot, small, omitInReference)
        famACCs = None
    
def drawAll(separatePops, omitInReference):
    allVirusHits, allSpecimens = getHits(omitInReference)
    
    fams = readFamFile()

    for acc in allVirusHits:
        doVirus = False  # only draw virus if it has at least one primary hit
        for specimen in allVirusHits[acc]:
            for hit in allVirusHits[acc][specimen].getHits():
                if hit.isPrimary():
                    doVirus = True
                    break
            if doVirus:
                break
        if doVirus:
            if acc not in fams:
                fam = getFamily(acc)
                fams[acc] = fam
                addFamily(acc, fam)
            drawVirus(acc, fams[acc], allVirusHits[acc], allSpecimens, separatePops, False, True, True, False, omitInReference)

###############################################################################

def doItAllViruses(dirName, filenameBAM, virusName, virusDB, alsoFindFeatures = True):
    
    filenameBAM = str(Path(dirName) / filenameBAM)
    newFilenameFASTA = getUnmappedReads(filenameBAM)
    newFilenameBAM = newFilenameFASTA[:-6] + '.bam'  # all unmapped plus mates
    
# Alternatively, blast reads first and then assemble only those that hit => shorter contigs
#    filenameCSV = blastViral(newFilenameFASTA, virusName, virusDB)
#    newFilenameBAM = getReads(filenameCSV, newFilenameBAM, True, virusName)  # just viral hits
    
    if virusName is not None:
        writeReadsFASTQ(newFilenameBAM, virusName)

        success = assembleReads(dirName, virusName)  #1 = assemble
        if success:
            blastScaffolds(dirName, virusName, virusDB, EVALUE_V, False)  #2 = blast scaffolds
            
#            hitsDirName = str(Path(dirName).parent.parent / ('results/HITS_' + virusName))
#            hitsCombinedDirName = str(Path(dirName).parent.parent / ('results/HITS_' + virusName + '_combined'))
        
            xmlFilename = getHits(dirName, virusName, 100, EVALUE_A)  #3 = get hits
            if xmlFilename is not None:
                drawXML(xmlFilename, str(Path(xmlFilename).parent / 'png'), FLANK_DRAW_LIMIT, False)
#                 if not Path(hitsDirName).exists():
#                     os.mkdir(hitsDirName)
                
                if alsoFindFeatures:
                    global iTrees
                    if iTrees is None:
                        iTrees = makeIntervalTrees(str(Path(dirName).parent))
                    xmlFeaturesFilename = findFeatures(xmlFilename, iTrees)  #4 = get features
    
                    displayXML(xmlFeaturesFilename, xmlFeaturesFilename[:-4] + '.txt', True)
                    displayXML(xmlFeaturesFilename, xmlFeaturesFilename[:-4] + '_short.txt', False)
                    
                    # remove any old txt files
                    if Path(xmlFilename[:-4] + '.txt').exists():
                        os.system('rm ' + xmlFilename[:-4] + '.txt')
                    if Path(xmlFilename[:-4] + '_short.txt').exists():
                        os.system('rm ' + xmlFilename[:-4] + '_short.txt')
        
#                     writelog('Copying hits/features to ' + hitsDirName + '...')
#                     os.system('cp -v ' + xmlFeaturesFilename + ' ' + hitsDirName)
                else:
                    displayXML(xmlFilename, xmlFilename[:-4] + '.txt', True)
                    displayXML(xmlFilename, xmlFilename[:-4] + '_short.txt', False)
#                     writelog('Copying hits to ' + hitsDirName + '...')
#                     os.system('cp -v ' + xmlFilename + ' ' + hitsDirName)
            else:
                writelog('*** ' + str(Path(dirName).name) + ': no results found. ***', True)
        else:
            writelog('*** ' + str(Path(dirName).name) + ': assembly failed; no results found. ***', True)
        
    writelog('*** ' + str(Path(dirName).name) + ' done. ***', True)

def doAll():
    dirList = [dir for dir in Path(ROOT_DIR).iterdir() if dir.is_dir() and 'combined' not in dir.name]
#    dirList = ['Amacuzac-Mexico-10.LIN210A1625', 'Bangkok_Thailand_10.LIN210A1688']
#    dirList = ['La_Lope-Gabon-10.LIN210A1646']
#    dirList = ['Tahiti_FrenchPolynesia_10.LIN210A1708']
#    dirList = [Path(ROOT_DIR + 'El_Dorado_Argentina_U_8.LIN210A2746')]
#    dirList = [Path(ROOT_DIR + 'Bangkok_Thailand_03.LIN210A1681')]

# import datetime
#     dirListNew = []
#     for dir in dirList:
#         dt = datetime.datetime.fromtimestamp(os.path.getmtime(str(dir / ('spades_all/' + dir.name + '_all_hits_features.xml'))))
#         if dt.day >= 12 and dt.hour >= 8:
#             continue
#         else:
#             dirListNew.append(dir)
#     dirList = dirListNew
    
    count = 0
    for dir in dirList:
        for file in dir.iterdir():
            if not file.is_dir() and (file.name[-26:] == '.sorted.deduped.merged.bam'):
                count += 1
                writelog(str(count) + '/' + str(len(dirList)) + ': ' + dir.name, True)
                doItAllViruses(str(dir.resolve()), file.name, 'all', 'virusdb5', True)  # **virusdb5, features**
                break
                    
def doAllParallel2():
    dirList = [dir for dir in Path(ROOT_DIR).iterdir() if dir.is_dir() and ('combined' not in dir.name)]
    workers = []
    count = 0
    for dir in dirList:
        for file in dir.iterdir():
            if not file.is_dir() and (file.name[-26:] == '.sorted.deduped.merged.bam'):
                p = Process(target = doItAllViruses, args = (str(dir.resolve()), file.name, 'all', 'virusdb', True))
                workers.append(p)
                p.start()
                count += 1
                
                while count == 2:
                    for index in range(len(workers)):
                        p = workers[index]
                        p.join(3)
                        if not p.is_alive():
                            count -= 1
                            workers.pop(index)
                            break
                        
                break
        
    for p in workers:
        p.join()
        
    os.system('cp ' + ROOT_DIR + '/*/spades_all/*all_hits.xml ' + str(Path(ROOT).parent) + '/results/HITS_all')
    consolidateAll(str(Path(ROOT_DIR).parent) + '/results/HITS_all', 'all')

logFile = ''

def writelog(message, alsoPrint = False):
    logFile.write(message + '\n')
    if alsoPrint:
        writelog(message)

def setup():
    global logFile
    
    logFilePath = Path(RESULTS_DIR) / LOGFILE_NAME
    logFile = open(logFilePath, 'a')
    writelog('\n************************************************')
    writelog('Starting pipeline at ' + time.strftime('%c'))

def main():
    setup()
    
    doAll()

    if copyResults():
        consolidateAll('all', True, True)
        drawAll(False, False)
        drawFamily([('Flaviviridae', 'NC_001564.2'), ('Orthomyxoviridae', 'MF176251.1'), ('Phenuiviridae', 'NC_038263.1'), ('Rhabdoviridae', 'NC_035132.1'), ('Totiviridae', 'NC_035674.1'), ('Xinmoviridae', 'MH037149.1')],None, False, False, True, False, False)

    
main()
#doAll()
#consolidateAll('all', True, True)
#xmlFilename = getHits(str(Path(ROOT_DIR + 'Cuanda_Angola_02.LIN210A1720')), 'all', 100, EVALUE_A)  #3 = get hits
