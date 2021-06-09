import matplotlib
matplotlib.use('agg')
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as pyplot
from Bio import SeqIO, Entrez, AlignIO
from Bio.Align import MultipleSeqAlignment
from Bio.SeqFeature import SeqFeature, FeatureLocation
from dna_features_viewer import GraphicFeature, GraphicRecord, BiopythonTranslator
#from intervaltree import Interval, IntervalTree
from lxml import etree as ET
from multiprocessing import Process
from pathlib import Path
import copy
import math
import numpy as np
import os
#import pickle
import pysam
import re
import sys
import time
import urllib.request as web

from util import *
from hits import Hit, Hits
from kmeans import findClusters
from features import *

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

    writelog('Getting unmapped reads...', True)

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
        
    return newFilenameBAM  # absolute path
    
###############################################################################

def writeReadsFASTQ(viralFilenameBAM):
    """Write paired-end reads from BAM file to 3 FASTQ files.
    
       Parameters:
           viralFilenameBAM: absolute path of BAM file containing viral reads to assemble
           
       Return value: None
    """
    
    writelog('Writing paired-end reads to FASTQ files...', True)
    
    spadesPath = Path(viralFilenameBAM).parent / 'spades'
    
    if not spadesPath.exists():
        os.system('mkdir ' + str(spadesPath))
        
    if (spadesPath / 'viral1.fastq').exists():
        writelog('   Skipping - FASTQ files exist.', True)
        return
    
    os.system('samtools fastq -1 ' + str(spadesPath / 'viral1.fastq') + ' -2 ' + str(spadesPath / 'viral2.fastq') + ' -s ' + str(spadesPath / 'viral_single.fastq') + ' ' + viralFilenameBAM)
    
###############################################################################

def assembleReads(dirName):
    """Assemble reads in dirName / spades into scaffolds.
    
       Parameter: 
           dirName:   absolute path of parent of directory containing fastq files to assemble 
            
       Return value: Boolean indicating whether assembly was successful
    """

    writelog('Assembling reads with SPAdes...', True)
    
    spadesPath = Path(dirName) / 'spades'
    if (spadesPath / 'scaffolds.fasta').exists():
        writelog('   Skipping - scaffolds.fasta exists.', True)
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
    
def blastScaffolds(dirName, force = False):
    """BLAST scaffolds against viral and aegypti genomes and combine to identify putative viral insertions."""
    
    spadesPath = Path(dirName) / 'spades'
    scaffoldsName = str(spadesPath / 'scaffolds.fasta')
    outVirusCSV = spadesPath / 'blast_scaffolds.csv'
    
    writelog('BLASTing scaffolds against viral database...', True)
    
    if not force and outVirusCSV.exists():
        writelog('   Skipping - ' + str(outVirusCSV) + ' exists.', True)
    else:
        os.system('blastn -query ' + scaffoldsName + ' -db ' + VIRUS_DB + ' -num_threads 16 -task blastn -evalue ' + str(EVALUE_V) +
                  ' -max_target_seqs 5000' +
                  ' -outfmt "10 qseqid qstart qend qseq sstart send sseq evalue bitscore sseqid stitle pident"' + 
                  ' -out ' + str(outVirusCSV))

###############################################################################
                  
def readCSV(csvFilename):
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
        contig = row[0]                   # qseqid = contig name
        sstart = int(row[4])              # start pos in viral genome
        send = int(row[5])                # end pos in viral genome
        length = abs(send - sstart) + 1   # length of viral hit
        sseqid = row[9]                   # viral genome accession number
        stitle = row[10].strip('|')       # viral genome title

        if contig in hits:
            if stitle not in accsFound[contig]:
                accsFound[contig][stitle] = [(sseqid, length)]
            else: 
                accsFound[contig][stitle].append((sseqid, length))
            hits[contig].append(tuple(row[1:10] + [row[10].strip('|')] + row[11:]))
        else:
            hits[contig] = [tuple(row[1:10] + [row[10].strip('|')] + row[11:])]
            accsFound[contig] = {}
            accsFound[contig][stitle] = [(sseqid, length)]
    csv.close()
    
    accsToUse = {}
    for contig in accsFound:
        accsToUse[contig] = {}
        for stitle in accsFound[contig]:
            refAccs = []
            maxLength = 0
            maxLengthAccs = []
            prefACCFound = False
            for acc, length in accsFound[contig][stitle]:
                if (stitle in PREFERRED_ACCS) and (acc in PREFERRED_ACCS[stitle]):
                    prefACCFound = True
                    break
                elif acc[:3] == 'ref':
                    refAccs.append(acc)
                elif length > maxLength:
                    maxLength = length
                    maxLengthAccs = [acc]
                elif length == maxLength:
                    maxLengthAccs.append(acc)
            if prefACCFound:
                accsToUse[contig][stitle] = PREFERRED_ACCS[stitle]
            elif refAccs != []:
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
    
def readCSVAA(csvFilename):
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
        length = len(row[3])       # length of hit in contig

        if contig in hits:
            hits[contig].append(tuple(row[1:]))
        else:
            hits[contig] = [tuple(row[1:])]

    return hits

def getHits(dirName): #, maxFlankingHits):
    """Combine viral and vector BLAST results to find putative insertions.
    
       Parameters:
           dirName:   absolute path of parent of directory containing SPAdes results 
           
       Return value: absolute path of XML file containing results
    """
    
    writelog('Combining BLAST results to locate putative insertions...', True)
    
    spadesPath = Path(dirName) / 'spades'
    resultsPath = Path(dirName) / 'results'
    if not resultsPath.exists():
        os.system('mkdir ' + str(resultsPath))
    scaffoldsPath = resultsPath / 'scaffolds'
    if not scaffoldsPath.exists():
        os.system('mkdir ' + str(scaffoldsPath))
    os.system('cp ' + str(spadesPath / 'scaffolds.fasta') + ' ' + str(scaffoldsPath))
    os.system('cp ' + str(spadesPath / 'blast_scaffolds.csv') + ' ' + str(scaffoldsPath))
    scaffoldsFilename = str(scaffoldsPath / 'scaffolds.fasta')
    csvVirusFilename = str(scaffoldsPath / 'blast_scaffolds.csv')
    
    xmlPath = resultsPath / 'xml'
    if not xmlPath.exists():
        os.system('mkdir ' + str(xmlPath))
    xmlFilename = str(xmlPath / (Path(dirName).name + '_hits.xml'))
   
    virusHits = readCSV(csvVirusFilename)
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
                pident = max(pident, float(virusHit[PIDENT]))
            else:    # finish this combined hit and start a new one
                if hitLength >= MIN_VIRUS_HIT_LENGTH:
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
                    writelog(contig + ': discarded hit "' + virusHit[STITLE].lstrip('|') + '" with length ' + str(hitLength) + ' < ' + str(MIN_VIRUS_HIT_LENGTH))

                # start new combined hit
                leftIndex = index
                hitLength = abs(int(virusHit[QEND]) - int(virusHit[QSTART])) + 1
                hitSeq = virusHit[QSEQ]
                pident = float(virusHit[PIDENT])
            rightIndex = index
            sseqid = nextSeqID
            
        # finish last combined hit
        if hitLength >= MIN_VIRUS_HIT_LENGTH:
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
            writelog(contig + ': discarded hit "' + virusHit[STITLE].lstrip('|') + '" with length ' + str(hitLength) + ' < ' + str(MIN_VIRUS_HIT_LENGTH))

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
                                   + ' -db aegyptidb -num_threads 8 -evalue ' + str(EVALUE_A)
                                   + ' -max_target_seqs 500'
                                   + ' -outfmt "10 qseqid qstart qend qseq sseqid sstart send sseq evalue bitscore"'
                                   + ' -out ' + str(spadesPath / 'blast_contig_temp.csv'))
        aedesHits = readCSVAA(str(spadesPath / 'blast_contig_temp.csv'))

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
            if contig in aedesHits:  # if there were aedes hits
                p = re.compile(r'.*length_(\d+)_.*')
                m = p.match(contig)
                contigLength = int(m.group(1))
                for index in range(len(aedesHits[contig])):
                    aedesHit = aedesHits[contig][index]
                    aedes_qseq = aedesHit[A_QSEQ]
                    aedes_qstart = int(aedesHit[A_QSTART])
                    aedes_qend = int(aedesHit[A_QEND])
                    if ((virus_qstart - aedes_qend > MAX_FLANKING_DISTANCE) or (aedes_qstart - virus_qend > MAX_FLANKING_DISTANCE)) and (len(aedes_qseq) < MIN_FLANKING_HIT_LENGTH):
                        writelog(contig + ': discarded distant Aedes hit with length ' + str(len(aedes_qseq)) + ' < ' + str(MIN_FLANKING_HIT_LENGTH))
                        continue
                    hitComplexity = complexity(aedes_qseq, COMPLEXITY_K)  # COMPLEXITY TEST
                    if hitComplexity < COMPLEXITY_CUTOFF:
                        writelog(contig + ': discarded Aedes hit with complexity ' + str(hitComplexity) + ' < ' + str(COMPLEXITY_CUTOFF))
                        continue
    
                    virusIntervalsNotInRef.chop(aedes_qstart, aedes_qend)  # detect overlap with virus hit
                    if aedes_qstart < virus_qstart and aedes_qend < virus_qend:
                        before.append(index)
                    elif aedes_qstart > virus_qstart and aedes_qend > virus_qend:
                        after.append(index)
                    else:  # as >= vs and ae <= ve or as <= vs and ae >= ve
                        overlap.append(index)  # aedes hit completely overlaps or is contained in virus hit

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
                            
            contigTree.attrib['referenceoverlaps'] = str([iv.data for iv in virusIntervalsInRef2])
                                        
            for indices, name in [(before, 'vectorhitleft'), (after, 'vectorhitright'), (overlap, 'vectorhitoverlap')]:
                for i in indices:
                    aedesHit = aedesHits[contig][i]
                    aedesElement = ET.SubElement(contigTree, name, {'id': str(i), 'seqid': aedesHit[A_SEQID]})
                    ET.SubElement(aedesElement, 'qstart').text = aedesHit[A_QSTART]
                    ET.SubElement(aedesElement, 'qend').text = aedesHit[A_QEND]
#                    ET.SubElement(aedesElement, 'qseq').text = aedesHit[A_QSEQ]
                    ET.SubElement(aedesElement, 'sstart').text = aedesHit[A_SSTART]
                    ET.SubElement(aedesElement, 'send').text = aedesHit[A_SEND]
#                    ET.SubElement(aedesElement, 'sseq').text = aedesHit[A_SSEQ]
                    ET.SubElement(aedesElement, 'evalue').text = aedesHit[A_EVALUE]
                    ET.SubElement(aedesElement, 'bitscore').text = aedesHit[A_BITSCORE]
                    
            matchingFlanks = ET.SubElement(contigTree, 'flanks')
            # find pairs of before/after aedes hits that are on the same strand and within 10Kbp of each other
            for b in before:
                beforeHit = aedesHits[contig][b]
                for a in after:
                    afterHit = aedesHits[contig][a]
                    if beforeHit[A_SEQID] == afterHit[A_SEQID]:  # same chr
                        beforeStart = int(beforeHit[A_SSTART])
                        beforeEnd = int(beforeHit[A_SEND])
                        afterStart = int(afterHit[A_SSTART])
                        afterEnd = int(afterHit[A_SEND])
                        flankDistance = max(beforeStart, beforeEnd, afterStart, afterEnd) - min(beforeStart, beforeEnd, afterStart, afterEnd)
                        if flankDistance <= MAX_FLANK_DISTANCE:
                            if (beforeStart < beforeEnd < afterStart + ALLOWED_OVERLAP < afterEnd + ALLOWED_OVERLAP) or \
                               (beforeStart > beforeEnd > afterStart - ALLOWED_OVERLAP > afterEnd - ALLOWED_OVERLAP):
                                ET.SubElement(matchingFlanks, 'match', {'leftid': str(b), 'rightid': str(a)})
                            else:
                                ET.SubElement(matchingFlanks, 'inversion', {'leftid': str(b), 'rightid': str(a)})
     
    tree.write(xmlFilename, xml_declaration=True, pretty_print=True)
        
    return xmlFilename

###############################################################################
    
def displayXML(fileName, displaySeq = True):
    """Write XML results to human-readable text file.
    
       Parameters:
           fileName:    absolute path of XML results file
           displaySeq:  whether to write long version with all sequence alignments
           
       Return value: None
    """
    
    txtPath = Path(fileName).parent.parent / 'txt'
    if not txtPath.exists():
        os.system('mkdir ' + str(txtPath))
        
    if displaySeq:
        outFileName = str(txtPath / (Path(fileName).name[:-4] + '.txt'))
    else:
        outFileName = str(txtPath / (Path(fileName).name[:-4] + '_short.txt'))
        
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
            outFile.write('      contig         {0:>5} - {1:<5} {2}\n'.format(v.find('qstart').text, v.find('qend').text, v.find('qseq').text))
            outFile.write('      {3:<14} {0:>5} - {1:<5} {2}\n\n'.format(v.find('sstart').text, v.find('send').text, v.find('sseq').text, v.attrib['seqid']))

        if displaySeq:
            hitsLeft = contig.findall('vectorhitleft')
            hitsRight = contig.findall('vectorhitright')
            hitsOverlap = contig.findall('vectorhitoverlap')
            for hits, locString in ((hitsLeft, 'left'), (hitsRight, 'right'), (hitsOverlap, 'overlapping')):
                if len(hits) > 0:
                    hitsDict = {}  # consolidate query hit: [subject hits]
                    for v in hits:
                        pos = (v.find('qstart').text, v.find('qend').text)
                        features = v.findall('feature')
                        featureString = '| '
                        for f in features:
                            featureString += f.attrib['type'] + '; ' + f.find('id').text + '; ' + f.find('location').text + ' | '

                        if pos not in hitsDict:
                            hitsDict[pos] = [(v.find('sstart').text, v.find('send').text, v.attrib['seqid'], featureString)]
                        else:
                            hitsDict[pos].append((v.find('sstart').text, v.find('send').text, v.attrib['seqid'], featureString))
                    outFile.write('   Vector hits ' + locString + ': \n')
                    posSort = list(hitsDict.keys())
                    posSort.sort(key = lambda pos: (int(pos[0]), int(pos[1])))
                    for q in posSort:
                        outFile.write('      contig         {0:>9} - {1:<9}\n'.format(*q))
                        for s in hitsDict[q]:
                            outFile.write('      {2:<14} {0:>9} - {1:<9} {3}\n'.format(*s))
                        outFile.write('\n')

            flanks = contig.find('flanks')
            for matches, message in [(flanks.findall('match'), 'Matching flanking regions:'), (flanks.findall('inversion'), 'Inverted matching flanking regions:')]:
                if len(matches) > 0:
                    outFile.write('   ' + message + '\n')
                    counter = 0
                    for match in matches:
                        outFile.write('      *** HIT ' + str(counter) + ' ***\n')
                        for v in hitsLeft:
                            if v.attrib['id'] == match.attrib['leftid']:
                                features = v.findall('feature')
                                featureString = '| '
                                for f in features:
                                    featureString += f.attrib['type'] + '; ' + f.find('id').text + '; ' + f.find('location').text + ' | '
                                outFile.write('         contig         {0:>9} - {1:<9}\n'.format(v.find('qstart').text, v.find('qend').text))
                                outFile.write('         {2:<14} {0:>9} - {1:<9} {3}\n\n'.format(v.find('sstart').text, v.find('send').text, v.attrib['seqid'], featureString))
                                break
                        
                        for v in contig.findall('virushit'):
                            outFile.write('         contig         {0:>9} - {1:<9} {2}\n'.format(v.find('qstart').text, v.find('qend').text, v.find('qseq').text))
                            outFile.write('         {3:<14} {0:>9} - {1:<9} {2}\n\n'.format(v.find('sstart').text, v.find('send').text, v.find('sseq').text, v.attrib['seqid']))

                        for v in hitsRight:
                            if v.attrib['id'] == match.attrib['rightid']:
                                features = v.findall('feature')
                                featureString = '| '
                                for f in features:
                                    featureString += f.attrib['type'] + '; ' + f.find('id').text + '; ' + f.find('location').text + ' | '
                                outFile.write('         contig         {0:>9} - {1:<9}\n'.format(v.find('qstart').text, v.find('qend').text))
                                outFile.write('         {2:<14} {0:>9} - {1:<9} {3}\n\n'.format(v.find('sstart').text, v.find('send').text, v.attrib['seqid'], featureString))
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
    
def drawXML(fileName):
    """Draw diagrams of contigs from XML file.
    
       Parameters:
           fileName:      absolute path of XML results file
           
       Return value: None
    """
    
    chromNames = {'NC_035107.1': 'Chr1', 'NC_035108.1': 'Chr2', 'NC_035109.1': 'Chr3'}
    
    diagramsPath = Path(fileName).parent.parent / 'diagrams'
    
    if not diagramsPath.exists():
        os.mkdir(str(diagramsPath))
    else:
        os.system('rm -r ' + str(diagramsPath) + '/*')
        
    fams = readFamFile()
    
    totalCounts = {}
    totalFeatures = 0
    totalContigs = 0
    
    contigFeaturesV = {}
    contigFeaturesA = {}
    contigLengths = {}
    families = {}
    
    labelFontSize = 6
    axesFontSize = 8
    thickness = 12
    
    tree = ET.parse(fileName)
    root = tree.getroot()
    for contig in root:
        family = None
        if BEST_HITS_ONLY and (contig.attrib['besthit'] != 'True'):
            continue
            
        p = re.compile(r'(.*)_length')
        m = p.search(contig.attrib['name'])
        c = m.group(1)                       # short contig name (NODE_xxx)
        if c not in contigFeaturesV:
            contigFeaturesV[c] = []
            contigFeaturesA[c] = []
            
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
                
            stitle += ' (' + str(seqid) + ')'
            if not BEST_HITS_ONLY and (contig.attrib['besthit'] == 'True'):                
                if seqid not in fams:
                    fam = getFamily(seqid)
                    fams[seqid] = fam
                    addFamily(seqid, fam)
                family = fams[seqid]
                
            if qend > qstart:
                strand = 1
            else:
                strand = -1
            gf = GraphicFeature(start=qstart, end=qend, strand=strand, thickness=thickness, linewidth=0,
                                           color='#ff0000', fontdict = {'size': labelFontSize},
                                           label = stitle.lstrip('|') + ' ' + str(sstart) + '-' + str(send) + ' (' + str(length) + ' bp; ' + str(pident) + '%)')  #; ' + str(evalue) + ')')
            features.append(gf)
            gfCopy = copy.deepcopy(gf)
            if not BEST_HITS_ONLY and (contig.attrib['besthit'] == 'True'):
                gfCopy.label = '**' + stitle.lstrip('|') + '** ' + str(sstart) + '-' + str(send) + ' (' + str(length) + ' bp; ' + str(pident) + '%)'
            contigFeaturesV[c].append(gfCopy)
        families[c] = family
        hitsLeft = contig.findall('vectorhitleft')
        hitsRight = contig.findall('vectorhitright')
        hitsOverlap = contig.findall('vectorhitoverlap')
        flanks = contig.find('flanks')
        matches = flanks.findall('match')
        hitsDone = []
        if len(matches) > 0:
            matchCount = 1
            for match in matches[:FLANK_DRAW_LIMIT]:
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
                            features.append(GraphicFeature(start=qstart, end=qend, strand=strand, thickness=thickness, linewidth=0,
                                                       color='#009051', fontdict = {'size': labelFontSize},
                                                       label = 'M' + str(matchCount) + ': ' + seqid + ' {0:,}-{1:,} '.format(sstart, send) + '(' + strandStr + ')'))
                            break
                matchCount += 1
                
        inversions = flanks.findall('inversion')
        if len(inversions) > 0:
            matchCount = 1
            for match in inversions[:max(0, FLANK_DRAW_LIMIT - len(matches))]:
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
                            features.append(GraphicFeature(start=qstart, end=qend, strand=strand, thickness=thickness, linewidth=0,
                                                       color='indigo', fontdict = {'size': labelFontSize},
                                                       label = 'I' + str(matchCount) + ': ' + seqid + ' {0:,}-{1:,} '.format(sstart, send) + '(' + strandStr + ')'))
                            break
                matchCount += 1

        aaCoverage = [0] * contigLength
        for hits in (hitsLeft, hitsRight, hitsOverlap):
            if len(hits) > 0:
                hitsDict = {}  # consolidate query hit: [subject hits]
                for v in hits:
                    pos = (v.find('qstart').text, v.find('qend').text)
                    for x in range(int(pos[0]), int(pos[1])):
                        aaCoverage[x] += 1
                    if v.attrib['id'] in hitsDone:
                        continue
                        
                    if pos not in hitsDict:
                        hitsDict[pos] = [(v.find('sstart').text, v.find('send').text, v.attrib['seqid'])] #, featureString)]
                    else:
                        hitsDict[pos].append((v.find('sstart').text, v.find('send').text, v.attrib['seqid'])) #, featureString))
                posSort = list(hitsDict.keys())
                if hits != hitsRight:
                    posSort.sort(key = lambda pos: int(pos[1]), reverse = True)  # closest first
                else:
                    posSort.sort(key = lambda pos: int(pos[0]))  # closest first
                for q in posSort[:max(0, FLANK_DRAW_LIMIT - len(matches) - len(inversions))]:
                    qstart = int(q[0])
                    qend = int(q[1])
                    if qend > qstart:
                        strand = 1
                    else:
                        strand = -1
                    s = hitsDict[q][0]
                    sstart = int(s[0])
                    send = int(s[1])
                    if send > sstart:
                        strandStr = '+'
                    else:
                        strandStr = '-'
                    seqid = s[2]
                    if seqid in chromNames:
                        seqid = chromNames[seqid]
                    elif seqid[:8] == 'NW_01873':
                        seqid = 'Scaffold ' + seqid[8:12]
                    if hits == hitsOverlap:
                        color = '#00ff00'
                    else:
                        color = '#0000ff'
                    features.append(GraphicFeature(start=qstart, end=qend, strand=strand, thickness=thickness, linewidth=0,
                                               color=color, fontdict = {'size': labelFontSize},
                                               label = seqid + ' {0:,}-{1:,} '.format(sstart, send) + '(' + strandStr + ')'))
                                
        record = GraphicRecord(sequence_length = contigLength, features = features)
        
        fig, (ax1, ax2) = pyplot.subplots(2, 1, sharex = True, figsize = (10, 6), gridspec_kw = {'height_ratios': [5, 1]})
        
        record.plot(max_label_length = 80, ax = ax1, with_ruler = False)
        ax2.fill_between(range(contigLength), aaCoverage, step = 'mid', alpha = 0.3)
        ax2.tick_params(axis='both', which='major', labelsize=axesFontSize)
        ax2.set_ylim(bottom = 0, top = max(aaCoverage + [1]))
        ax2.set_yticks([0, max(aaCoverage + [1]) // 2, max(aaCoverage + [1])])
        ax2.set_ylabel('Aa hits', fontsize = axesFontSize)

        familyDir = str(diagramsPath / family)
        if not Path(familyDir).exists():
            os.mkdir(familyDir)
        if contig.attrib['besthit'] == 'True':
            pp = PdfPages(str(Path(familyDir) / (contig.attrib['name'] + '_BEST-HIT.pdf')))
            pp.savefig(fig)
            pp.close()
        else:
            pp = PdfPages(str(Path(familyDir) / (contig.attrib['name'] + '.pdf')))
            pp.savefig(fig)
            pp.close()
        pyplot.close(fig)
        
    if not BEST_HITS_ONLY:
        for c in contigFeaturesV:
            record = GraphicRecord(sequence_length = contigLengths[c], features = contigFeaturesV[c])
        
            ax, _ = record.plot(max_label_length = 80, figure_width = 10)
            familyDir = str(diagramsPath / families[c])
            pp = PdfPages(str(Path(familyDir) / (c + '_all-viral-hits.pdf')))
            pp.savefig(ax.figure)
            pp.close()
            pyplot.close(ax.figure)
                    
###############################################################################

def getSeqRecord(fileNameFASTA):
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
        
def writeSummaryTable(fileName, viralHits, viralHits0, viralHits1, viralHits2, allFamilies, allViruses):
    out = open(fileName, 'w')

    famVir = list(zip(allFamilies, allViruses))
        
    specimens = list(viralHits.keys())
    specimenLabels = {}
    for specimen in specimens:
        specimenLabels[getSpecimenLabel(specimen)] = specimen
    specimenLabelsList = list(specimenLabels.keys())
    specimenLabelsList.sort()
    
    counts = {}
    for fam, (stitle, seqid) in famVir:
        counts[(fam, (stitle, seqid))] = 0
        
    for region, num in specimenLabelsList:
        specimen = specimenLabels[(region, num)]
        for fam, (stitle, seqid) in famVir:
            foundInSpecimen = False
            for hit in viralHits[specimen]:  # hit = (seqid, stitle, hitDict) in one contig
                if hit[0] == seqid:
                    if OMIT_EVES_IN_REFERENCE and (specimen in viralHits0) and (hit in viralHits0[specimen]):  # overlapping vector hit
                        continue
                    foundInSpecimen = True
                    break
            if foundInSpecimen:
                counts[(fam, (stitle, seqid))] += 1
                
    newFamVir = []
    for fam, (stitle, seqid) in famVir:
        if counts[(fam, (stitle, seqid))] >= MIN_HITS_TO_SHOW_VIRUS:
            newFamVir.append((fam, (stitle, seqid)))
    famVir = newFamVir
    
    famVir.sort(key=lambda x: x[1][0])                  # sort by stitle
    famVir.sort(key=lambda x: counts[x], reverse=True)  # then by counts
    famVir.sort(key=lambda x: x[0])                     # then by family
    
    out.write('\t')
    for fam, (stitle, seqid) in famVir:
        out.write(seqid + '\t')
    out.write('\n')
    out.write('\t')
    for fam, (stitle, seqid) in famVir:
        out.write(stitle + '\t')
    out.write('\n')
    
    out.write('\t')
    for fam, (stitle, seqid) in famVir:
        out.write(fam + '\t')
    out.write('\n')
    
    # Write total virus hit counts
    out.write('\t')
    for fam, (stitle, seqid) in famVir:
        out.write(str(counts[(fam, (stitle, seqid))]) + '\t')
    out.write('\n')
    
    for region, num in specimenLabelsList:
        specimen = specimenLabels[(region, num)]
        out.write(region + ' ' + str(num) + '\t')
        for fam, (stitle, seqid) in famVir:
            length = 0
            rangeList = []
            for hit in viralHits[specimen]:  # hit = (seqid, stitle, hitDict) in one contig
                if hit[0] == seqid:
                    if OMIT_EVES_IN_REFERENCE and (specimen in viralHits0) and (hit in viralHits0[specimen]):  # overlapping vector hit
                        continue
                    hitDict = hit[2]
                    qpos = list(hitDict.keys())
                    qpos.sort()                  # within hit, sort by qstart
                    s = []
                    for p in qpos:
                        length += abs(hitDict[p][1] - hitDict[p][0]) + 1
                        s.extend([hitDict[p][0], hitDict[p][1]])
                    if s[0] > s[-1]:
                        s = s[::-1]
                    if (specimen in viralHits0) and (hit in viralHits0[specimen]):  # overlapping vector hit
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
                rangeString = rangeString[:-2] + ('*' * s[-2]) 
                if len(s[-1]) > 0:     # if features to show
                    rangeString += ' (' + s[-1] + ')'
                rangeString += ' | '
            if length > 0:
                out.write(rangeString[:-1] + '| ' + str(length))
                counts[(fam, (stitle, seqid))] += 1
            out.write('\t')
        out.write(specimen + '\n')
        
    out.close()   

def consolidateAll():
    """Consolidate all results.
       Write a tab-delimited table summarizing hit regions for all specimens.
       Write a FASTA file containing specimen EV regions aligned to viral reference genome
    
       Parameters:
           virusLength: length of the virus reference genome
           
       Return value: None
    """
    
    writelog('Consolidating results in ' + RESULTS_DIR + '...', True)
        
    viralHits1 = {}
    viralHits2 = {}
    viralHits0 = {}   # with AA overlap
    viralHits = {}    # all together
    allHits = []    
    viralSeqs = {}
    refSeqs = {}
    pIdents = {}
    insertSites = {}
    
    xmlFiles = []
    specimensPath = Path(SPECIMENS_DIR)
    for specimenDir in specimensPath.iterdir():
        if specimenDir.name[0] != '.':
            for file in (specimenDir / 'results/xml').iterdir():
                if (file.name[0] != '.') and (('_hits_features.xml' in file.name) or (('_hits.xml' in file.name) and not Path(str(file)[:-4] + '_features.xml').exists())):
                    xmlFiles.append(file.resolve())
    xmlFiles.sort()
    
    specimen2Label = {}
    label2Specimen = {}
    specimenCount = 1
    for file in xmlFiles:
        specimen = file.name.split('_hits')[0]
        region, num = getSpecimenLabel(specimen)
        region = region.replace(' ', '-')
        specimen2Label[specimen] = (region, num)
        label2Specimen[(region, num)] = specimen
        writelog('  Processing ' + str(specimenCount) + '/' + str(len(xmlFiles)) + ': ' + specimen, True)
        specimenCount += 1
        seqPath = file.parent.parent / 'sequences'
        if not seqPath.exists():
            os.system('mkdir ' + str(seqPath))
        virus_fasta = open(str(seqPath / (specimen + '_hits_aligned.fasta')), 'w')
        virus_fasta_percontig = open(str(seqPath / (specimen + '_hits_unaligned.fasta')), 'w')

        # Parse XML tree.
        
        tree = ET.parse(str(file))
        root = tree.getroot()
        for contig in root.iterchildren():   # one per virus per contig
            contigName = contig.attrib['name']
            hit = {}
            v = contig.find('virushit')
            seqid = v.attrib['seqid']            # should be same for all v
            if '|' in seqid:
                seqid = seqid.split('|')[1]
            stitle = v.attrib['stitle']
            if '|' in stitle:
                stitle = stitle.lstrip('|')
            
            if BEST_HITS_ONLY and (contig.attrib['besthit'] != 'True'):  # and (stitle not in PREFERRED_ACCS):  # last condition may result in >1 hit per original contig
                continue
                
            if OMIT_EVES_IN_REFERENCE and (contig.attrib['inreference'] == 'True'):
                continue
            
            if (seqid, stitle) not in viralSeqs:
                viralSeqs[(seqid, stitle)] = {}
                refSeqs[(seqid, stitle)] = {}
                pIdents[(seqid, stitle)] = {}
                insertSites[(seqid, stitle)] = {}
            if specimen not in viralSeqs[(seqid, stitle)]:
                viralSeqs[(seqid, stitle)][specimen] = {}
                refSeqs[(seqid, stitle)][specimen] = {}
                pIdents[(seqid, stitle)][specimen] = {}
                insertSites[(seqid, stitle)][specimen] = {}
            viralSeqs[(seqid, stitle)][specimen][contigName] = {}
            refSeqs[(seqid, stitle)][specimen][contigName] = {}
            pIdents[(seqid, stitle)][specimen][contigName] = {}
            insertSites[(seqid, stitle)][specimen][contigName] = {}
            for v in contig.findall('virushit'):
                hit[(int(v.find('qstart').text), int(v.find('qend').text))] = (int(v.find('sstart').text), int(v.find('send').text))
                viralSeqs[(seqid, stitle)][specimen][contigName][(int(v.find('sstart').text), int(v.find('send').text))] = v.find('qseq').text
                refSeqs[(seqid, stitle)][specimen][contigName][(int(v.find('sstart').text), int(v.find('send').text))] = v.find('sseq').text
                pIdents[(seqid, stitle)][specimen][contigName][(int(v.find('sstart').text), int(v.find('send').text))] = float(v.find('pident').text)
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
                combinedSeq += viralSeqs[(seqid, stitle)][specimen][contigName][hit[qse]].replace('-', '')  # remove gaps
                if revSeq:
                    sposString = str(hit[qse][1]) + '-' + str(hit[qse][0]) + ',' + sposString
                else:
                    sposString += str(hit[qse][0]) + '-' + str(hit[qse][1]) + ','
            virus_fasta_percontig.write('>' + seqid + '_' + stitle.replace(' ', '-') + '_' + sposString[:-1] + '\n')
            if revSeq:
                combinedSeq = reverseComplement(combinedSeq)
            virus_fasta_percontig.write(combinedSeq + '\n') 

            recordedIDs = []
            insertSites[(seqid, stitle)][specimen][contigName]['match'] = []
            flanks = contig.find('flanks')
            vectorHitsLeft = contig.findall('vectorhitleft')
            vectorHitsRight = contig.findall('vectorhitright')
            for match in flanks.findall('match'):
                recordedIDs.append(match.attrib['leftid'])
                recordedIDs.append(match.attrib['rightid'])
                for v in vectorHitsLeft:
                    if v.attrib['id'] == match.attrib['leftid']:
                        matchLeft = int(v.find('send').text)
                        break
                for v in vectorHitsRight:
                    if v.attrib['id'] == match.attrib['rightid']:
                        insertSites[(seqid, stitle)][specimen][contigName]['match'].append((v.attrib['seqid'], (matchLeft, int(v.find('sstart').text))))
                        break
                        
            insertSites[(seqid, stitle)][specimen][contigName]['inversion'] = []
            for match in flanks.findall('inversion'):
                recordedIDs.append(match.attrib['leftid'])
                recordedIDs.append(match.attrib['rightid'])
                for v in vectorHitsLeft:
                    if v.attrib['id'] == match.attrib['leftid']:
                        matchLeft = int(v.find('send').text)
                        break
                for v in vectorHitsRight:
                    if v.attrib['id'] == match.attrib['rightid']:
                        insertSites[(seqid, stitle)][specimen][contigName]['inversion'].append((v.attrib['seqid'], (matchLeft, int(v.find('sstart').text))))
                        break
        
            insertSites[(seqid, stitle)][specimen][contigName]['left'] = []
            for v in vectorHitsLeft:
                if v.attrib['id'] not in recordedIDs:
                    aedes_qstart = int(v.find('qstart').text)
                    aedes_qend = int(v.find('qend').text)
                    virus_qstart = qpos[0][0]
                    virus_qend = qpos[-1][1]
                    if (-ALLOWED_OVERLAP <= virus_qstart - aedes_qend <= MAX_FLANKING_DISTANCE) or (-ALLOWED_OVERLAP <= aedes_qstart - virus_qend <= MAX_FLANKING_DISTANCE): # or (len(aedes_qseq) >= MIN_FLANKING_HIT_LENGTH):
                        insertSites[(seqid, stitle)][specimen][contigName]['left'].append((v.attrib['seqid'], int(v.find('send').text)))
                
            insertSites[(seqid, stitle)][specimen][contigName]['right'] = []
            for v in vectorHitsRight:
                if v.attrib['id'] not in recordedIDs:
                    aedes_qstart = int(v.find('qstart').text)
                    aedes_qend = int(v.find('qend').text)
                    virus_qstart = qpos[0][0]
                    virus_qend = qpos[-1][1]
                    if (-ALLOWED_OVERLAP <= virus_qstart - aedes_qend <= MAX_FLANKING_DISTANCE) or (-ALLOWED_OVERLAP <= aedes_qstart - virus_qend <= MAX_FLANKING_DISTANCE): # or (len(aedes_qseq) >= MIN_FLANKING_HIT_LENGTH):
                        insertSites[(seqid, stitle)][specimen][contigName]['right'].append((v.attrib['seqid'], int(v.find('sstart').text)))
        
            numVectorHits = len(vectorHitsLeft + vectorHitsRight)
        
            features = []
            featureSummary = contig.find('featuresummary')
            if featureSummary is not None:
                for f in featureSummary:
                    if f.tag != 'None':
                        features.append((f.tag, int(f.text) / numVectorHits))
                    
            if specimen not in viralHits:
                viralHits[specimen] = []
            viralHits[specimen].append((seqid, stitle, hit, features))
        
            flanks = contig.find('flanks')
            matches = flanks.findall('match')
            if len(matches) > 0:
                if specimen not in viralHits2:
                    viralHits2[specimen] = []
                viralHits2[specimen].append((seqid, stitle, hit, features))
            elif numVectorHits > 0:
                if specimen not in viralHits1:
                    viralHits1[specimen] = []
                viralHits1[specimen].append((seqid, stitle, hit, features))
            
            if contig.attrib['inreference'] == 'True':
                if specimen not in viralHits0:
                    viralHits0[specimen] = []
                viralHits0[specimen].append((seqid, stitle, hit, features))
            
            if (seqid, stitle, hit) not in allHits:
                allHits.append((seqid, stitle, hit, features))
            
        virus_fasta.close()
        virus_fasta_percontig.close()
        
    if not Path(VIRUSES_DIR).exists():
        os.system('mkdir ' + VIRUSES_DIR)
    if not Path(SEQUENCES_DIR).exists():
        os.system('mkdir ' + SEQUENCES_DIR)
        
    fams = readFamFile()
    positions = {}
    viralSeqsBySpecimen = {}
    for (seqid, stitle) in viralSeqs:
        if seqid not in fams:
            fam = getFamily(seqid)
            fams[seqid] = fam
            addFamily(seqid, fam)
        for specimen in viralSeqs[(seqid, stitle)]:
            for contigName in viralSeqs[(seqid, stitle)][specimen]:
                newDict = {}
                newDictRef = {}
                for (start, end) in viralSeqs[(seqid, stitle)][specimen][contigName]:
                    if start > end:
                        newDict[(end, start)] = reverseComplement(viralSeqs[(seqid, stitle)][specimen][contigName][(start, end)])
                        newDictRef[(end, start)] = reverseComplement(refSeqs[(seqid, stitle)][specimen][contigName][(start, end)])
                    else:
                        newDict[(start, end)] = viralSeqs[(seqid, stitle)][specimen][contigName][(start, end)]
                        newDictRef[(start, end)] = refSeqs[(seqid, stitle)][specimen][contigName][(start, end)]
                viralSeqs[(seqid, stitle)][specimen][contigName] = newDict
                refSeqs[(seqid, stitle)][specimen][contigName] = newDictRef
        
        refIndices = {}    # specimen seq at each reference index
        positions[(seqid, stitle)] = {}     # virus positions for fasta headers
        for specimen in viralSeqs[(seqid, stitle)]:
            refIndices[specimen] = {}
            positions[(seqid, stitle)][specimen] = {}
            for contigName in viralSeqs[(seqid, stitle)][specimen]:
                refIndices[specimen][contigName] = {}
                positions[(seqid, stitle)][specimen][contigName] = []
                for (start, end) in viralSeqs[(seqid, stitle)][specimen][contigName]:
                    refIndex = start - 2  # -1 for 0-indexing, -1 for initial increment
                    seq = viralSeqs[(seqid, stitle)][specimen][contigName][(start, end)]
                    refSeq = refSeqs[(seqid, stitle)][specimen][contigName][(start, end)]
                    positions[(seqid, stitle)][specimen][contigName].append((start, end))
                    for index in range(len(seq)):
                        if refSeq[index] != '-':
                            refIndex += 1
                            refIndices[specimen][contigName][refIndex] = seq[index]
                        else:
                            refIndices[specimen][contigName][refIndex] += seq[index]
                viralSeqs[(seqid, stitle)][specimen][contigName] = ''
                    
        referenceFilename = FASTA_DIR + seqid + '.fasta'
        refRecord = getSeqRecord(referenceFilename)
        referenceSeq = str(refRecord.seq)
        virusLength = len(referenceSeq)
        referenceID = refRecord.id
        assert referenceID == seqid
        viralSeqs[(seqid, stitle)][referenceID] = ''
        refIndicesBySpecimen = {}
        viralSeqsBySpecimen[(seqid, stitle)] = {}
        for specimen in refIndices:
            refIndicesBySpecimen[specimen] = {}
            viralSeqsBySpecimen[(seqid, stitle)][specimen] = ''
            for contigName in refIndices[specimen]:
                for refIndex in refIndices[specimen][contigName]:
                    if (refIndex not in refIndicesBySpecimen[specimen]) or (len(refIndices[specimen][contigName][refIndex]) > len(refIndicesBySpecimen[specimen][refIndex])):
                        refIndicesBySpecimen[specimen][refIndex] = refIndices[specimen][contigName][refIndex]
        
        for refIndex in range(virusLength):
            maxLength = max([len(refIndicesBySpecimen[specimen].get(refIndex, '-')) for specimen in refIndices])
            viralSeqs[(seqid, stitle)][referenceID] += '{0:-<{1}}'.format(referenceSeq[refIndex], maxLength)
            for specimen in refIndices:
                viralSeqsBySpecimen[(seqid, stitle)][specimen] += '{0:-<{1}}'.format(refIndicesBySpecimen[specimen].get(refIndex, ''), maxLength)
                for contigName in refIndices[specimen]:
                    viralSeqs[(seqid, stitle)][specimen][contigName] += '{0:-<{1}}'.format(refIndices[specimen][contigName].get(refIndex, ''), maxLength)

        FAMILY_DIR = SEQUENCES_DIR + fams[seqid]
        if not Path(FAMILY_DIR).exists():
            os.system('mkdir ' + FAMILY_DIR)
        
        seqOutFileName = str(Path(FAMILY_DIR) / ('sequences_' + seqid + '_per_specimen_aligned.fasta'))
        seqOutPerContigFileName = str(Path(FAMILY_DIR) / ('sequences_' + seqid + '_per_contig_aligned.fasta'))
        seqOutPerContigUnalignedFileName = str(Path(FAMILY_DIR) / ('sequences_' + seqid + '_per_contig_unaligned.fasta'))
        seqOut = open(seqOutFileName, 'w')
        seqOutPerContig = open(seqOutPerContigFileName, 'w')
        seqOutPerContigUnaligned = open(seqOutPerContigUnalignedFileName, 'w')
        sortedSpecimens = list(viralSeqs[(seqid, stitle)].keys())
        sortedSpecimens.remove(referenceID)
        sortedSpecimens.sort()
        
        theLength = len(viralSeqs[(seqid, stitle)][referenceID])
        
        seqOut.write('>' + referenceID + ' ' + stitle + '\n')
        seqOut.write(viralSeqs[(seqid, stitle)][referenceID] + '\n')
        seqOutPerContig.write('>' + referenceID + ' ' + stitle + '\n')
        seqOutPerContig.write(viralSeqs[(seqid, stitle)][referenceID] + '\n')
        for specimen in sortedSpecimens:
            region, num = specimen2Label[specimen]
            posSpecimen = []
            for contigName in positions[(seqid, stitle)][specimen]:
                posSpecimen.extend(positions[(seqid, stitle)][specimen][contigName])
            posSpecimen.sort()
            posString = ''
            for start, end in posSpecimen:
                posString += str(start) + '-' + str(end) + ','
                
            seqOut.write('>' + region + '_' + str(num) + '_' + seqid + '_' + posString[:-1] + '\n')
            seqOut.write(viralSeqsBySpecimen[(seqid, stitle)][specimen] + '\n')
            contigCount = 1
            for contigName in viralSeqs[(seqid, stitle)][specimen]:
            
                if len(viralSeqs[(seqid, stitle)][specimen][contigName]) != theLength:
                    writelog('Sequence length error: ' + seqid + ', ' + specimen + ', ' + contigName + ' (length = ' + str(len(viralSeqs[(seqid, stitle)][specimen][contigName])) + ', should be ' + str(theLength) + ') - FIXED', True)
                    viralSeqs[(seqid, stitle)][specimen][contigName] = viralSeqs[(seqid, stitle)][specimen][contigName][:theLength]
                    
                positions[(seqid, stitle)][specimen][contigName].sort()
                posString = ''
                for start, end in positions[(seqid, stitle)][specimen][contigName]:
                    posString += str(start) + '-' + str(end) + ','
                    
                # compute weighted pident and add to headers
                top = bottom = 0
                for start, end in pIdents[(seqid, stitle)][specimen][contigName]:
                    top += pIdents[(seqid, stitle)][specimen][contigName][(start, end)] * (abs(end - start) + 1)
                    bottom += abs(end - start) + 1
                pident = top / bottom
                
                seqOutPerContig.write('>' + specimen + '_|_' + contigName + '_|_' + seqid + '_' + posString[:-1] + '_pident_' + '{:5.3f}'.format(pident) + '\n')
                seqOutPerContig.write(viralSeqs[(seqid, stitle)][specimen][contigName] + '\n')
                
                seqOutPerContigUnaligned.write('>' + specimen + '_|_' + contigName + '_|_' + seqid + '_' + posString[:-1] + '_pident_' + '{:5.3f}'.format(pident) + '\n')
                seqOutPerContigUnaligned.write(viralSeqs[(seqid, stitle)][specimen][contigName].replace('-', '') + '\n')

                contigCount += 1
        seqOut.close()
        seqOutPerContig.close()
        seqOutPerContigUnaligned.close()
        
        # Cluster.
        if DO_CLUSTERING:
            clusteredFileName = findClusters(seqOutPerContigFileName)
    
    allViruses = []
    for entry in allHits:
        allViruses.append((entry[1], entry[0]))  # stitle, seqid
    allViruses = list(set(allViruses))
    allViruses.sort()
    
    allFamilies = []
    for stitle, seqid in allViruses:
        if seqid not in fams:
            fam = getFamily(seqid)
            fams[seqid] = fam
            addFamily(seqid, fam)
        allFamilies.append(fams[seqid])

    # Write family sequence files
        
    for family in set(allFamilies):
        FAMILY_DIR = SEQUENCES_DIR + family
        if not Path(FAMILY_DIR).exists():
            os.system('mkdir ' + FAMILY_DIR)
            
        # get seqids in family
        famACCs = []
        for stitle, seqid in allViruses:
            if seqid not in fams:
                fam = getFamily(seqid)
                fams[seqid] = fam
                addFamily(seqid, fam)
            if fams[seqid] == family:
                famACCs.append((seqid, stitle))
                
        # get seqid for most hits in family to use as reference/representative
        repACC = famACCs[0]
        for seqid, stitle in famACCs[1:]:
            if len(viralSeqs[(seqid, stitle)]) > len(viralSeqs[repACC]):
                repACC = (seqid, stitle)
            
        # find start positions for leftmost CDS
        startPositions = {}
        for seqid, stitle in famACCs:
            if not Path(GB_DIR + seqid + '.gb').exists():
                if not getGenBank(seqid):
                    continue
            virusRecord = SeqIO.read('/home/havill/data/aegypti/gb/' + seqid + '.gb', 'gb')  # SeqRecord
            minCDSPos = math.inf
            for feature in virusRecord.features:
                if feature.type == 'CDS':
                    minCDSPos = min(minCDSPos, feature.location.start)
            startPositions[(seqid, stitle)] = minCDSPos
                    
        seqOut = open(str(Path(FAMILY_DIR) / ('sequences_' + family + '_per_specimen_aligned.fasta')), 'w')
        seqOutPerContig = open(str(Path(FAMILY_DIR) / ('sequences_' + family + '_per_contig_aligned.fasta')), 'w')
        seqOutPerContigUnaligned = open(str(Path(FAMILY_DIR) / ('sequences_' + family + '_per_contig_unaligned.fasta')), 'w')
        
        seqOut.write('>' + repACC[0] + ' ' + repACC[1] + '\n')
        seqOut.write(viralSeqs[repACC][repACC[0]] + '\n')
        seqOutPerContig.write('>' + repACC[0] + ' ' + repACC[1] + '\n')
        seqOutPerContig.write(viralSeqs[repACC][repACC[0]] + '\n')
        repLength = len(viralSeqs[repACC][repACC[0]])
        
        for seqid, stitle in famACCs:
            adjust = startPositions[repACC] - startPositions[(seqid, stitle)]
            addGaps = abs(adjust)
            sortedSpecimens = list(viralSeqs[(seqid, stitle)].keys())
            try:
                sortedSpecimens.remove(seqid)
            except ValueError:
                pass
            sortedSpecimens.sort()
            for specimen in sortedSpecimens:
                region, num = getSpecimenLabel(specimen)
                region = region.replace(' ', '_')
                posSpecimen = []
                for contigName in positions[(seqid, stitle)][specimen]:
                    posSpecimen.extend(positions[(seqid, stitle)][specimen][contigName])
                posSpecimen.sort()
                posString = ''
                for start, end in posSpecimen:
                    posString += str(start) + '-' + str(end) + ','
                
                # compute weighted pident and add to headers
                
                seqOut.write('>' + region + '_' + str(num) + '_' + seqid + '_' + posString[:-1] + '\n')
                if adjust > 0:
                    seqOut.write(('{:-<' + str(repLength) + '}').format('-' * addGaps + viralSeqsBySpecimen[(seqid, stitle)][specimen][:repLength-addGaps]) + '\n')
                else:
                    seqOut.write(('{:-<' + str(repLength) + '}').format(viralSeqsBySpecimen[(seqid, stitle)][specimen][addGaps:repLength]) + '\n')
                contigCount = 1
                for contigName in viralSeqs[(seqid, stitle)][specimen]:
                    positions[(seqid, stitle)][specimen][contigName].sort()
                    posString = ''
                    for start, end in positions[(seqid, stitle)][specimen][contigName]:
                        posString += str(start) + '-' + str(end) + ','
                    seqOutPerContig.write('>' + region + '_' + str(num) + '__contig' + str(contigCount) + '_' + seqid + '_' + posString[:-1] + '\n')
                    seqOutPerContigUnaligned.write('>' + region + '_' + str(num) + '__contig' + str(contigCount) + '_' + seqid + '_' + posString[:-1] + '\n')
                    seqOutPerContigUnaligned.write(viralSeqs[(seqid, stitle)][specimen][contigName].replace('-', '') + '\n')
                    if adjust > 0:
                        seqOutPerContig.write(('{:-<' + str(repLength) + '}').format('-' * addGaps + viralSeqs[(seqid, stitle)][specimen][contigName][:repLength-addGaps]) + '\n')
                    else:
                        seqOutPerContig.write(('{:-<' + str(repLength) + '}').format(viralSeqs[(seqid, stitle)][specimen][contigName][addGaps:repLength]) + '\n')
                    contigCount += 1
        seqOut.close()
        seqOutPerContig.close()
        seqOutPerContigUnaligned.close()

    writeSummaryTable(str(Path(RESULTS_DIR) / 'results.tsv'), viralHits, viralHits0, viralHits1, viralHits2, allFamilies, allViruses)

###############################################################################
               
AaegL5_hits = {'NC_001564.2': [9330, 8303],   # Cell fusing agent virus
               'NC_034017.1': [1677, 839, 2269, 3801, 5188, 6416, 7537, 8296, 10083, 8327],  # Xishuangbanna aedes flavivirus
               'NC_038512.1': [2579, 3877],   # Trichoplusia ni TED virus
               'NC_035674.1': [446, 3299],    # Australian Anopheles totivirus
               'NC_035132.1': [6568, 7294],   # Culex ohlsrhavirus
               'NC_028484.1': [6555, 7282],   # Tongilchon ohlsrhavirus
               'NC_040669.1': [6773, 7497],   # Riverside virus 1
               'NC_025384.1': [6789, 6267]}   # Culex tritaeniorhynchus rhabdovirus
               
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
	
class MyCustomTranslator(BiopythonTranslator):
    def compute_feature_color(self, feature):
        if 'UTR' in feature.type:
            return "red"
        elif feature.type in ["mat_peptide", "CDS"]:
            return "gold"
        else:
            return "blue"

    def compute_feature_label(self, feature):
        if feature.type == "mat_peptide":
            return feature.qualifiers['product'][0]
        elif feature.type == "CDS":
            return feature.qualifiers['product'][0]
        else:
            return BiopythonTranslator.compute_feature_label(self, feature)

    def compute_filtered_features(self, features):
        filtered_features1 = [feature for feature in features if feature.type not in ('source', 'gene')]
        
        if len(filtered_features1) <= 4:
            return filtered_features1
        else:
            return [feature for feature in filtered_features1 if 'product' not in feature.qualifiers or feature.qualifiers['product'][0] != 'polyprotein']

def drawVirus(acc, family, hits, allSpecimens, separatePops, isFamily, showAaegL5_hits, showPlot, small):    

    if not Path(VIRUSES_DIR).exists():
        os.system('mkdir ' + VIRUSES_DIR)
    DIAGRAMS_DIR = VIRUSES_DIR + 'diagrams/'
    if not Path(DIAGRAMS_DIR).exists():
        os.system('mkdir ' + DIAGRAMS_DIR)
    FAMILY_DIR = DIAGRAMS_DIR + family
    if not Path(FAMILY_DIR).exists():
        os.system('mkdir ' + FAMILY_DIR)
    
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
                
    virusRecord2 = SeqIO.read(GB_DIR + acc + '.gb', 'gb')  # SeqRecord
    if isFamily:
        virusName = 'Family ' + family + ' (' + virusRecord2.id + ': ' + virusRecord2.description[:80]
        if len(virusRecord2.description) > 80:
            virusName += '...'
        virusName += ')'
        writelog('  Creating diagram' + 's' * separatePops + ' for family ' + family + ' (using representative ' + virusRecord2.id + ')...', True)
    else:
        virusName = virusRecord2.id + ': ' + virusRecord2.description[:80]
        if len(virusRecord2.description) > 80:
            virusName += '...'
        virusName += ' (' + family + ')'
        writelog('  Creating diagram' + 's' * separatePops + ' for ' + virusName + '...', True)
                
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
                lc = (0, 0, 0, 0)
                for index in range(0, len(hit), 2):
                    start = hit[index]
                    end = hit[index + 1]
                    features[specimen].append(GraphicFeature(start=start, end=end, strand=1, color=color, label = None, thickness = thickness, linewidth = 0, linecolor = lc))
                if hit.getData() != []:  # add reference overlap regions
                    lc = (1, 1, 1) + (transparency,)    # white
                    for start, end in hit.getData():
                        features[specimen].append(GraphicFeature(start=start, end=end, strand=1, color=color, label = None, thickness = thickness, linewidth = 1, linecolor = lc))
        
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
        
        pp = PdfPages(str(Path(FAMILY_DIR) / figName))
        pp.savefig(fig)
        pp.close()
        pyplot.close(fig)

def getHitsForDiagram():
    writelog('  Reading hits for diagrams...', True)
    dir = Path(SPECIMENS_DIR)
    subdirs = [d for d in dir.iterdir()]
    files = [d / ('results/xml/' + d.name + '_hits.xml') for d in subdirs ]
    files.sort()

    allVirusHits = {}
    allSpecimens = []
    for file in files:
        specimen = file.name.split('_hits')[0]
#        writelog('Getting hits for ' + specimen + '...', True)
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
            if OMIT_EVES_IN_REFERENCE and overlap:
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

def drawFamily(families, famACCs, separatePops, showAaegL5_hits, showPlot, small):
    allVirusHits, allSpecimens = getHitsForDiagram()
        
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
            virusRecord = SeqIO.read(GB_DIR + acc + '.gb', 'gb')  # SeqRecord
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
                
        drawVirus(repACC, family, famHits, allSpecimens, separatePops, True, showAaegL5_hits, showPlot, small)
        famACCs = None
    
def drawAll(separatePops = False):
    writelog('Creating diagrams...', True)
    allVirusHits, allSpecimens = getHitsForDiagram()
    
    fams = readFamFile()

    for acc in allVirusHits:
        hitCount = sum(hit.isPrimary() for specimen in allVirusHits[acc] for hit in allVirusHits[acc][specimen].getHits())
        if hitCount >= MIN_HITS_TO_SHOW_VIRUS:
            if acc not in fams:
                fam = getFamily(acc)
                fams[acc] = fam
                addFamily(acc, fam)
            drawVirus(acc, fams[acc], allVirusHits[acc], allSpecimens, separatePops, False, True, True, False)

###############################################################################

def doIt(dirName, filenameBAM):
    
    filenameBAM = str(Path(dirName) / filenameBAM)
    newFilenameBAM = getUnmappedReads(filenameBAM)
    writeReadsFASTQ(newFilenameBAM)

    success = assembleReads(dirName)  #1 = assemble
    if success:
        blastScaffolds(dirName, False)  #2 = blast scaffolds
    
        xmlFilename = getHits(dirName)  #3 = get hits
        if xmlFilename is not None:
            drawXML(xmlFilename)
            
            if FIND_FEATURES:
                makeIntervalTrees()
                xmlFeaturesFilename = findFeaturesXML(xmlFilename)  #4 = get features
                displayXML(xmlFeaturesFilename, True)
                displayXML(xmlFeaturesFilename, False)
                
                # remove any old txt files
                if Path(xmlFilename[:-4] + '.txt').exists():
                    os.system('rm ' + xmlFilename[:-4] + '.txt')
                if Path(xmlFilename[:-4] + '_short.txt').exists():
                    os.system('rm ' + xmlFilename[:-4] + '_short.txt')
            else:
                displayXML(xmlFilename, True)
                displayXML(xmlFilename, False)
        else:
            writelog('*** ' + str(Path(dirName).name) + ': no results found. ***', True)
    else:
        writelog('*** ' + str(Path(dirName).name) + ': assembly failed; no results found. ***', True)
        
    writelog('*** ' + str(Path(dirName).name) + ' done. ***', True)

def doAll():
    dirList = [dir for dir in Path(SPECIMENS_DIR).iterdir() if dir.is_dir() and 'combined' not in dir.name]
    
    count = 0
    dirList.sort()
    for dir in dirList:
        for file in dir.iterdir():
            if not file.is_dir() and (file.name[-26:] == '.sorted.deduped.merged.bam'):
                count += 1
                writelog('\n' + str(count) + '/' + str(len(dirList)) + ': ' + dir.name, True)
                doIt(str(dir.resolve()), file.name)
                break
                    
# def doAllParallel2():
#     dirList = [dir for dir in Path(SPECIMENS_DIR).iterdir() if dir.is_dir() and ('combined' not in dir.name)]
#     workers = []
#     count = 0
#     for dir in dirList:
#         for file in dir.iterdir():
#             if not file.is_dir() and (file.name[-26:] == '.sorted.deduped.merged.bam'):
#                 p = Process(target = doIt, args = (str(dir.resolve()), file.name, 'virusdb', True))
#                 workers.append(p)
#                 p.start()
#                 count += 1
#                 
#                 while count == 2:
#                     for index in range(len(workers)):
#                         p = workers[index]
#                         p.join(3)
#                         if not p.is_alive():
#                             count -= 1
#                             workers.pop(index)
#                             break
#                         
#                 break
#         
#     for p in workers:
#         p.join()
#         
#     os.system('cp ' + ROOT_DIR + '/*/spades_all/*hits.xml ' + str(Path(ROOT).parent) + '/results/HITS_all')
#     consolidateAll(str(Path(ROOT_DIR).parent) + '/results/HITS_all', 'all')

def main():
    if Path(RESULTS_DIR).exists():
        answer = ' '
        while answer not in 'yYnN':
            answer = input('Results directory ' + RESULTS_DIR + ' exists.  OK to overwrite (y/n)? ')
        if answer[0] not in 'yY':
            writelog('OK, quitting.', True)
            return 1
    else:
        os.system('mkdir ' + RESULTS_DIR)
        
    writelog('\n************************************************')
    writelog('Starting pipeline at ' + time.strftime('%c'))
    
    doAll()

    consolidateAll()
    drawAll(False)
#    drawFamily([('Flaviviridae', 'NC_001564.2'), ('Orthomyxoviridae', 'MF176251.1'), ('Phenuiviridae', 'NC_038263.1'), ('Rhabdoviridae', 'NC_035132.1'), ('Totiviridae', 'NC_035674.1'), ('Xinmoviridae', 'MH037149.1')],None, False, False, True, False, False)
    
main()
#consolidateAll()
