import sys
import os
import pysam
from pathlib import Path
from multiprocessing import Process
#import xml.etree.ElementTree as ET
from lxml import etree as ET
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio import SeqIO
from intervaltree import Interval, IntervalTree
import re
import pickle
from dna_features_viewer import GraphicFeature, GraphicRecord
import matplotlib.pyplot as pyplot
import copy
import urllib.request as web

SPADES_EXEC = '/home/havill/bin/SPAdes-3.14.1-Linux/bin/spades.py'
ROOT = '/home/havill/data/aegypti/analyzed'
MIN_PIDENT = 80

###############################################################################

def getUnmappedReads(filenameBAM):
    """Get unmapped reads AND MATES from a BAM file and write them to a FASTA file.
    
       Parameter:
           filenameBAM: absolute path of input BAM file
           
       Return value: absolute path of the output FASTA file
    """
    
    print('\nAnalyzing ' + filenameBAM)
    
    try:
        index = filenameBAM.index('.sorted')
    except:
        index = filenameBAM.index('.bam')
    newFilenameBAM = filenameBAM[:index] + '_unmapped_with_mates.bam'  # absolute path

    print('Getting unmapped reads...')

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

        print('Total reads =', count_total)
        print('Unmapped reads =', count, '({:.2f}%)'.format(100*count/count_total))
    else:
        print('Skipping - ' + newFilenameBAM + ' exists.')

    # Convert unmapped reads to FASTA.

    print('     Writing unmapped reads to a FASTA file...')
    newFilenameFASTA = newFilenameBAM[:-4] + '.fasta'  # absolute path
    if not Path(newFilenameFASTA).exists():
        os.system('samtools fasta ' + newFilenameBAM + ' > ' + newFilenameFASTA)
    else:
        print('Skipping - ' + newFilenameFASTA + ' exists.')
        
    return newFilenameFASTA  # absolute path
    
###############################################################################

# def blastViral(filenameFASTA, virusName, virusDB):
#     """BLAST unmapped reads against viral genome.
#     
#        Parameter:
#            filenameFASTA: absolute path of FASTA file containing unmapped reads
#            virusName:     string containing name of virus
#            virusDB:       string containing name of BLAST database for virus
#            
#        Return value: absolute path of CSV file containing BLAST results
#     """
# 
#     print('Searching for viral reads with BLAST...')
#     
#     blastOutputFileNameCSV = Path(filenameFASTA).parent / ('blast_unmapped_' + virusName + '.csv')
#     if blastOutputFileNameCSV.exists():
#         print('Skipping - ' + str(blastOutputFileNameCSV) + ' exists.')
#     else:
#         os.system('blastn -query ' + filenameFASTA + ' -db ' + virusDB + ' -num_threads 8 -task blastn -evalue 0.0001' + 
#                   ' -outfmt "10 qaccver saccver stitle pident length mismatch gapopen qstart qend sstart send evalue bitscore"' + 
#                   ' -out ' + str(blastOutputFileNameCSV))
# 
#     return str(blastOutputFileNameCSV)  # absolute path
    
###############################################################################
    
def getReads(filenameCSV, filenameBAM, getMates, virusName):
    """Get reads identified in BLAST CSV file into a BAM file.
    
       Parameters:
           filenameCSV: absolute path of BLAST CSV file
           filenameBAM: absolute path of original BAM file
           getMates:    Boolean indicating whether to include paired-end mates of hits
           virusName:   name of the virus of interest
           
       Return value: absolute path of the BAM file containing the viral reads
    """
    
    print('Getting viral sequences (plus paired-end mates)...')
    
    csv = open(filenameCSV, 'r')
    queries = {}
    for line in csv:
        cols = line.split(',')
        query = cols[0]
        queries[query] = ''
    csv.close()

    bamfile = pysam.AlignmentFile(filenameBAM, 'rb', threads = 8)
    try:
        index = filenameBAM.index('.sorted')
    except ValueError:
        index = filenameBAM.index('.bam')
    viralFilenameBAM = filenameBAM[:index] + '_' + virusName + '.bam'
    bamfile_out = pysam.AlignmentFile(viralFilenameBAM, 'wb', template = bamfile, threads = 8)

    written = {}
    for read in bamfile:
        if read.is_read1:
            qn = read.query_name + '/1'
            qn_mate = read.query_name + '/2'
        elif read.is_read2:
            qn = read.query_name + '/2'
            qn_mate = read.query_name + '/1'
        else:
            qn = read.query_name
        read.query_name += '_VIRAL'
        if getMates:
            if (qn in queries or qn_mate in queries) and qn not in written:
                written[qn] = True
                bamfile_out.write(read)
        else:
            if qn in queries and qn not in written:
                bamfile_out.write(read)
                written[qn] = True

    bamfile.close()
    bamfile_out.close()
    
    return viralFilenameBAM
    
###############################################################################

def writeReadsFASTQ(viralFilenameBAM, virusName):
    """Write paired-end reads from BAM file to 3 FASTQ files.
    
       Parameters:
           viralFilenameBAM: absolute path of BAM file containing viral reads to assemble
           virusName:        string containing the name of the virus
           
       Return value: None
    """
    
    print('Writing paired-end reads to FASTQ files...')
    
    spadesPath = Path(viralFilenameBAM).parent / ('spades_' + virusName)
    
    if not spadesPath.exists():
        os.system('mkdir ' + str(spadesPath))
        
    if (spadesPath / 'viral1.fastq').exists():
        print('  Skipping ... fastq files exist.')
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

    print('Assembling reads with SPAdes...')
    
    spadesPath = Path(dirName) / ('spades_' + virusName)
    if (spadesPath / 'scaffolds.fasta').exists():
        print('  Skipping ... scaffolds.fasta exists.')
        return True
        
#     if (spadesPath / 'params.txt').exists():
#         print('  Skipping ... params.txt exists.')
#         return True
        
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
    
def blastScaffolds(dirName, virusName, virusDB):
    """BLAST scaffolds against viral and aegypti genomes and combine to identify putative viral insertions."""
    
    spadesPath = Path(dirName) / ('spades_' + virusName)
    scaffoldsName = str(spadesPath / 'scaffolds.fasta')
    outVirusCSV = spadesPath / ('blast_scaffolds_' + virusName + '.csv')
#    outAACSV = str(spadesPath / 'blast_scaffolds_aa.csv')
    
    print('BLASTing scaffolds against viral database...')
    
    if outVirusCSV.exists():
        print('Skipping - ' + str(outVirusCSV) + ' exists.')
    else:
        os.system('blastn -query ' + scaffoldsName + ' -db ' + virusDB + ' -num_threads 16 -task blastn -evalue 1e-10' +
                  ' -outfmt "10 qseqid qstart qend qseq sstart send sseq evalue bitscore sseqid stitle pident"' + 
                  ' -out ' + str(outVirusCSV))

#     print('BLASTing scaffolds against Aedes aegypti genome...')
#     os.system('blastn -query ' + scaffoldsName + ' -db aegyptidb -num_threads 8 -max_target_seqs 10' +
#               ' -outfmt "10 qseqid qstart qend qseq sseqid sstart send sseq evalue bitscore" -out ' + outAACSV)

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
    for line in csv:
        line = line.rstrip()
        row = line.split(',')
        contig = row[0]            # qseqid = contig name
#         qstart = row[1]          # start pos in contig
#         qend = row[2]            # end pos in contig
#         qseq = row[3]            # aligned sequence in contig
        length = len(row[3])       # length of hit in contig
#         sstart = row[4]          # start pos in viral genome
#         send = row[5]            # end pos in viral genome
#         sseq = row[6]            # aligned sequence in viral genome
#         evalue = float(row[7])   # e-value of alignment
#         bitscore = float(row[8]) # bit score of alignment
#         sseqid = row[9]          # viral genome accession number
#         stitle = row[10]         # viral genome title
#         pident = row[11]         # % identity of alignment

        if length >= minLength:
            if contig in hits:
                hits[contig].append(tuple(row[1:]))
            else:
                hits[contig] = [tuple(row[1:])]
    csv.close()

    return hits
    
def getName(accessionID):
	"""
	Query NCBI to get name corresponding to an accession number.
	"""
	
	Entrez.email = 'havill@denison.edu'
	try:
		handle = Entrez.efetch(db = 'nucleotide', id = accessionID, rettype = 'gb', retmode = 'text')
		result = handle.read().split('\n')
		for line in result:
			if 'ORGANISM' in line:
				return ' '.join(line.split()[1:])
	except:
		print('Exception:', sys.exc_info()[1])
		return ''


def getHits(dirName, virusName, minLength):
    """Combine viral and vector BLAST results to find putative insertions.
    
       Parameters:
           dirName:   absolute path of parent of directory containing SPAdes results 
           virusName: name of virus under consideration
           
       Return value: absolute path of XML file containing results
    """
    
    print('Combining BLAST results to locate putative insertions...')
    
    spadesPath = Path(dirName) / ('spades_' + virusName)
    
    csvVirusFilename = str(spadesPath / ('blast_scaffolds_' + virusName + '.csv'))
    csvAedesFilename = str(spadesPath / 'blast_scaffolds_aa.csv')
    
    xmlFilename = str(spadesPath / (Path(dirName).name + '_' + virusName + '_hits.xml'))

    virusHits = readCSV(csvVirusFilename, 0)
    if virusHits is None:
        return None
#    aedesHits = readCSV(csvAedesFilename, 0)
    
    scaffoldsFilename = str(spadesPath / 'scaffolds.fasta')
    scaffoldRecords = SeqIO.index(scaffoldsFilename, 'fasta')

    ### Combine overlapping viral hits in each contig into long hits
    # Combine hits from the same virus
    newVirusHits = {}
    maxLenPIdent = {}
    for contig in virusHits:
        newVirusHits[contig] = []
        leftIndex = 0
        rightIndex = 0
        virusHits[contig].sort(key = lambda hit: (hit[8], int(hit[0]))) # (sseqid, qstart)
        virusHit = virusHits[contig][0]
        hitLength = len(virusHit[2])
        pident = float(virusHit[10])  # max pident among all hits for current virus
        maxHitLength = hitLength
        maxHitLengthIndex = 0
        maxPIdent = pident
        maxPIdentIndex = 0
        sseqid = virusHit[8]  # virus accession number
        
        for index in range(1, len(virusHits[contig])):
            virusHit = virusHits[contig][index]
            nextSeqID = virusHit[8]
            if (nextSeqID == sseqid): # and (nextLeft <= right + 1):
                hitLength += len(virusHit[2])
                pident = max(pident, float(virusHit[10]))
            else:
                if hitLength >= minLength:
                    newVirusHits[contig].append((leftIndex, rightIndex))
                    if hitLength > maxHitLength:
                        maxHitLength = hitLength
                        maxHitLengthIndex = leftIndex
                    if pident > maxPIdent:
                        maxPIdent = pident
                        maxPIdentIndex = leftIndex
                else:
                    print(contig + ': hit "' + virusHit[9][:20].lstrip('|') + '" with length ' + str(hitLength) + ' less than ' + str(minLength))
                leftIndex = index
                hitLength = len(virusHit[2])
                pident = float(virusHit[10])
            rightIndex = index
            sseqid = nextSeqID
            
        if hitLength >= minLength:
            newVirusHits[contig].append((leftIndex, rightIndex))
            if hitLength > maxHitLength:
                maxHitLength = hitLength
                maxHitLengthIndex = leftIndex
            if pident > maxPIdent:
                maxPIdent = pident
                maxPIdentIndex = leftIndex
        else:
            print(contig + ': hit "' + virusHit[9][:20].lstrip('|') + '" with length ' + str(hitLength) + ' less than ' + str(minLength))

        maxLenPIdent[contig] = (maxHitLength, maxHitLengthIndex, maxPIdent, maxPIdentIndex)
        
    # Combine virus and aedes hit results to locate putative insertions.
    
    root = ET.Element('root')
    tree = ET.ElementTree(root)

#
    all_viral_hits = []
#
    
    for contig in virusHits:
        print(contig)
        SeqIO.write(scaffoldRecords[contig], str(spadesPath / 'contig_temp.fasta'), 'fasta')
        os.system('blastn -query ' + str(spadesPath / 'contig_temp.fasta') + ' -db aegyptidb -num_threads 8' +
                  ' -outfmt "10 qseqid qstart qend qseq sseqid sstart send sseq evalue bitscore" -out ' + str(spadesPath / 'blast_contig_temp.csv'))
        aedesHits = readCSV(str(spadesPath / 'blast_contig_temp.csv'), 0)

        if contig not in aedesHits:   # only include contig if also an aa hit
            print('  no aa hits')
#            continue
#       
        if contig in aedesHits:
            print('  ' + str(len(aedesHits[contig])) + ' aa hits')
        print('  ' + str([v[8] for v in virusHits[contig]]))
        all_viral_hits.extend([(v[8], v[9]) for v in virusHits[contig]])
#        
        
        minEvalueIndex = 0
        for index in range(1, len(virusHits[contig])):
            if float(virusHits[contig][index][6]) < float(virusHits[contig][minEvalueIndex][6]):
                minEvalueIndex = index
        maxBitscoreIndex = 0
        for index in range(1, len(virusHits[contig])):
            if float(virusHits[contig][index][7]) > float(virusHits[contig][maxBitscoreIndex][7]):
                maxBitscoreIndex = index
                
        maxHitLength, maxHitLengthIndex, maxPIdent, maxPIdentIndex = maxLenPIdent[contig]
        
        thisContigCount = 0
        for leftIndex, rightIndex in newVirusHits[contig]:       # partition aedes hits into before and after viral hit
            thisContigCount += 1
            virus_qstart = int(virusHits[contig][leftIndex][0])
            virus_qend = int(virusHits[contig][rightIndex][1])
            
            print(contig + ' ' + str(leftIndex) + ', ' + str(rightIndex))
            
            bestHit = False
            if (maxPIdent >= MIN_PIDENT):
                if maxPIdentIndex == leftIndex:
                    bestHit = True
            else:
                if maxHitLengthIndex == leftIndex:
                    bestHit = True
            
            before = []
            after = []            
            if contig in aedesHits:
                for index in range(len(aedesHits[contig])):
                    aedesHit = aedesHits[contig][index]
                    aedes_qstart = int(aedesHit[0])
                    aedes_qend = int(aedesHit[1])
                    aedes_chr = aedesHit[3]
                    aedes_sstart = int(aedesHit[4])
                    aedes_send = int(aedesHit[5])
                    if aedes_qstart <= virus_qstart: # aedes_qend < virus_qstart + 15:
                        before.append((index, aedes_sstart, aedes_send, aedes_chr))
                    else: # if aedes_qstart > virus_qend - 15:
                        after.append((index, aedes_sstart, aedes_send, aedes_chr))

            if len(newVirusHits[contig]) > 1:                     # one "contig" in XML per virus
                contigName = contig + '_' + str(thisContigCount)
            else:
                contigName = contig
            contigTree = ET.SubElement(root, 'contig', {'name': contigName, 'besthit': str(bestHit)})
            #outFileBothEnds.write('Contig: ' + contig + '\n')
            #outFileBothEnds.write('\n   Viral hit:' + '\n')
            for index in range(leftIndex, rightIndex + 1):                
                virusHit = virusHits[contig][index]
                vElement = ET.SubElement(contigTree, 'virushit', {'seqid': virusHit[8], 'stitle': virusHit[9], 
                                                                  'minevalue': str(index == minEvalueIndex), 
                                                                  'maxbitscore': str(index == maxBitscoreIndex)})
                ET.SubElement(vElement, 'qstart').text = virusHit[0]
                ET.SubElement(vElement, 'qend').text = virusHit[1]
                ET.SubElement(vElement, 'qseq').text = virusHit[2]
                ET.SubElement(vElement, 'sstart').text = virusHit[3]
                ET.SubElement(vElement, 'send').text = virusHit[4]
                ET.SubElement(vElement, 'sseq').text = virusHit[5]
                ET.SubElement(vElement, 'evalue').text = virusHit[6]
                ET.SubElement(vElement, 'bitscore').text = virusHit[7]
                ET.SubElement(vElement, 'pident').text = virusHit[10]
                                        
            matchingFlanks = ET.SubElement(contigTree, 'flanks')
            # find pairs of before/after aedes hits that are on the same strand and within 10Kbp of each other
            for b in before:
                for a in after:
                    if b[3] == a[3]:
                        if ((b[1] < b[2] < a[1] + 20 < a[2] + 20) or (b[1] > b[2] > a[1] + 20 > a[2] + 20)) and (abs(b[2] - a[1]) < 10000):
                            ET.SubElement(matchingFlanks, 'match', {'leftid': str(b[0]), 'rightid': str(a[0])})
            for b in before:
                beforeHit = aedesHits[contig][b[0]]
                aedesElement = ET.SubElement(contigTree, 'vectorhitleft', {'id': str(b[0]), 'seqid': beforeHit[3]})
                ET.SubElement(aedesElement, 'qstart').text = beforeHit[0]
                ET.SubElement(aedesElement, 'qend').text = beforeHit[1]
                ET.SubElement(aedesElement, 'qseq').text = beforeHit[2]
                ET.SubElement(aedesElement, 'sstart').text = beforeHit[4]
                ET.SubElement(aedesElement, 'send').text = beforeHit[5]
                ET.SubElement(aedesElement, 'sseq').text = beforeHit[6]
                ET.SubElement(aedesElement, 'evalue').text = beforeHit[7]
                ET.SubElement(aedesElement, 'bitscore').text = beforeHit[8]
            if len(after) > 0:
                for a in after:
                    afterHit = aedesHits[contig][a[0]]
                    aedesElement = ET.SubElement(contigTree, 'vectorhitright', {'id': str(a[0]), 'seqid': afterHit[3]})
                    ET.SubElement(aedesElement, 'qstart').text = afterHit[0]
                    ET.SubElement(aedesElement, 'qend').text = afterHit[1]
                    ET.SubElement(aedesElement, 'qseq').text = afterHit[2]
                    ET.SubElement(aedesElement, 'sstart').text = afterHit[4]
                    ET.SubElement(aedesElement, 'send').text = afterHit[5]
                    ET.SubElement(aedesElement, 'sseq').text = afterHit[6]
                    ET.SubElement(aedesElement, 'evalue').text = afterHit[7]
                    ET.SubElement(aedesElement, 'bitscore').text = afterHit[8]
     
    tree.write(xmlFilename, xml_declaration=True, pretty_print=True)
    
    print(set(all_viral_hits))
    for a in set(all_viral_hits):
        print(a[0], a[1])
    
    return xmlFilename

###############################################################################

def getGFFFeatures(fileName, iTrees):
    """Read repeat features from a GFF3 file and insert into dictionary of interval trees."""
    
    chromNumbers = {'1': 'NC_035107.1', '2': 'NC_035108.1', '3': 'NC_035109.1'}
    strands = {'+': +1, '-': -1, '.': 0}
   
    gff = open(fileName, 'r')
    
    for line in gff:
        if line[0] == '#':
            continue
        
        cols = line.strip().split('\t')
        
        seqid = cols[0]
        if seqid == 'MT':
            continue
        elif seqid in chromNumbers:
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

# def getGFFRepeatFeatures(fileName, iTrees):
#     """Read repeat features from a GFF3 file and insert into dictionary of interval trees."""
#     
#     chromNumbers = {'1': 'NC_035107.1', '2': 'NC_035108.1', '3': 'NC_035109.1'}
#     strands = {'+': +1, '-': -1, '.': 0}
#    
#     gff = open(fileName, 'r')
#     
#     for line in gff:
#         if line[0] == '#':
#             continue
#         
#         cols = line.strip().split('\t')
#         seqid = cols[0]
#         if seqid == 'MT':
#             continue
#         if seqid in chromNumbers:
#             seqid = chromNumbers[seqid]
#         else:
#             seqid += '.1'
#         type = cols[2]
#         if type in ['chromosome']:
#             continue
#         attrib = cols[8]
#         
#         if seqid not in iTrees:
#             iTrees[seqid] = IntervalTree()
#         if type == 'supercontig':
#             altID = re.findall(r'Alias=' + seqid + r',(.*)', attrib)  # in supercontig record
#             if len(altID) > 0:
#                 altID = altID[0]
#                 iTrees[altID] = iTrees[seqid]
#         else:
#             start = int(cols[3])
#             end = int(cols[4])
#             strand = strands[cols[6]]
#             f = SeqFeature()
#             f.location = FeatureLocation(start, end + 1, strand = strand)
#             f.type = type
#             findName = re.findall(r'Name=(.*?);.*', attrib)
#             if len(findName) > 0:
#                 f.id = findName[0]
#             else:
#                 findName = re.findall(r'ID=(.*?);.*', attrib)  # try again
#                 if len(findName) > 0:
#                     f.id = findName[0]
#                 else:
#                     f.id = 'Unknown'
#             iTrees[seqid][start:end+1] = f  # interval start <= x < end+1
#     
#     gff.close()

# def getGenBankFeatures(fileName, iTrees):
#     """Read features from a multi-record flat GenBank file and insert into dictionary of interval trees."""
# 
#     for seqrecord in SeqIO.parse(fileName, 'genbank'):
#         iTrees[seqrecord.id] = IntervalTree()
#         altID = re.findall(r'NIGP01.*\.1', seqrecord.annotations['comment'])
#         if len(altID) > 0:
#             altID = altID[0]
#             iTrees[altID] = iTrees[seqrecord.id]
#         for f in seqrecord.features:
#             if f.type != 'source':
#                 iTrees[seqrecord.id][f.location.start:f.location.end - 1] = f
                
def makeIntervalTrees(dirName):
    pickleName = str(Path(dirName) / 'featuretrees.pickle')
    if Path(pickleName).exists():
        print('Pickled feature trees found.  Reading...')
        picklediTreesDict = open(pickleName, 'rb')
        iTrees = pickle.load(picklediTreesDict)
        picklediTreesDict.close()
    else:
        iTrees = {}
        print('Reading base features...')
        # https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/002/204/515/GCF_002204515.2_AaegL5.0/
        getGFFFeatures('/home/havill/data/aegypti/gb/Aedes-aegypti-LVP_AGWG_BASEFEATURES_AaegL5.2.gff3', iTrees)
    
        print('Reading repeat features...')
        getGFFFeatures('/home/havill/data/aegypti/gb/Aedes-aegypti-LVP_AGWG_REPEATFEATURES_AaegL5.gff3', iTrees)
    
        print('Writing feature trees to disk...')
        picklediTreesDict = open(pickleName, 'wb')
        pickle.dump(iTrees, picklediTreesDict)
        picklediTreesDict.close()
        
    return iTrees
    
###############################################################################

def populateFeature(feature, f):
    """Populate an XML feature element from a SeqFeature object.
       f has only location, type, id attributes."""
    
    
    ET.SubElement(feature, 'location').text = str(f.location)
    ET.SubElement(feature, 'type').text = f.type
    #ET.SubElement(feature, 'qualifiers').text = str(f.qualifiers)
    
    ET.SubElement(feature, 'id').text = f.id
    
#     if 'gene' in f.qualifiers:
#         ET.SubElement(feature, 'gene').text = f.qualifiers['gene'][0]
#     else:
#         ET.SubElement(feature, 'gene').text = 'None'
#         
#     if 'product' in f.qualifiers:
#         ET.SubElement(feature, 'product').text = f.qualifiers['product'][0]
#     else:
#         ET.SubElement(feature, 'product').text = 'None'
#         
#     if 'unknown' not in f.id:
#         ET.SubElement(feature, 'id').text = f.id
#     elif 'protein_id' in f.qualifiers:
#         ET.SubElement(feature, 'id').text = f.qualifiers['protein_id'][0]
#     elif 'transcript_id' in f.qualifiers:
#         ET.SubElement(feature, 'id').text = f.qualifiers['transcript_id'][0]
#     elif 'gene' in f.qualifiers:
#         ET.SubElement(feature, 'id').text = f.qualifiers['gene'][0]
#     else:
#         ET.SubElement(feature, 'id').text = 'Unknown'

# def getIntrons(location):
#     """Get a list of intron intervals from a CompoundLocation object representing exons."""
#     
#     prev = location.start - 1
#     intervals = []
#     for pos in sorted(location):
#         if pos > prev + 1:
#             intervals.append(Interval(prev+1, pos-1))
#         prev = pos
#     return intervals

def searchAA(iTrees, seqid, start, end, hitElement, counts):
    """Get features for a given interval and write to xml tree."""
    
    searchInterval = Interval(start, end+1)
    
    overlapIntervals = iTrees[seqid][start:end]
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
                 
        feature = ET.SubElement(hitElement, 'feature', {'type': featureType})
        populateFeature(feature, f)
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
    
    print('Annotating features...')
    
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
            
        featureSummary = ET.SubElement(contig, 'featuresummary')
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
            if len(hitsLeft) > 0:
                hitsLeftDict = {}  # consolidate query hit: [subject hits]
                for v in hitsLeft:
                    pos = (v.find('qstart').text, v.find('qend').text, v.find('qseq').text)
                    features = v.findall('feature')
                    featureString = '| '
                    for f in features:
#                         if displaySeq:
                        featureString += f.attrib['type'] + '; ' + f.find('id').text + '; ' + f.find('location').text + ' | '
#                         else:
#                             featureString += f.attrib['type'] + '; ' + f.find('id').text + ' | '
                    if pos not in hitsLeftDict:
                        hitsLeftDict[pos] = [(v.find('sstart').text, v.find('send').text, v.find('sseq').text, v.attrib['seqid'], featureString)]
                    else:
                        hitsLeftDict[pos].append((v.find('sstart').text, v.find('send').text, v.find('sseq').text, v.attrib['seqid'], featureString))
#                 if displaySeq:
                outFile.write('   Vector hits left: \n')
                posSort = list(hitsLeftDict.keys())
                posSort.sort(key = lambda pos: (int(pos[0]), int(pos[1])))
                for q in posSort:
#                     if displaySeq:
                    outFile.write('      contig         {0:>9} - {1:<9} {2}\n'.format(*q))
                    for s in hitsLeftDict[q]:
                        outFile.write('      {3:<14} {0:>9} - {1:<9} {2} {4}\n'.format(*s))
                    outFile.write('\n')
#                     else:
#                         outFile.write('      contig         {0:>9} - {1:<9}\n'.format(*q))
#                         for s in hitsLeftDict[q]:
#                             outFile.write('      {3:<14} {0:>9} - {1:<9} {4}\n'.format(*s))
#                         outFile.write('\n')
        
            hitsRight = contig.findall('vectorhitright')
            if len(hitsRight) > 0:
                hitsRightDict = {}  # consolidate query hit: [subject hits]
                for v in hitsRight:
                    pos = (v.find('qstart').text, v.find('qend').text, v.find('qseq').text)
                    features = v.findall('feature')
                    featureString = '| '
                    for f in features:
#                         if displaySeq:
                        featureString += f.attrib['type'] + '; ' + f.find('id').text + '; ' + f.find('location').text + ' | '
#                         else:
#                             featureString += f.attrib['type'] + '; ' + f.find('id').text + ' | '
                    if pos not in hitsRightDict:
                        hitsRightDict[pos] = [(v.find('sstart').text, v.find('send').text, v.find('sseq').text, v.attrib['seqid'], featureString)]
                    else:
                        hitsRightDict[pos].append((v.find('sstart').text, v.find('send').text, v.find('sseq').text, v.attrib['seqid'], featureString))
#                if displaySeq:
                outFile.write('   Vector hits right: \n')
                posSort = list(hitsRightDict.keys())
                posSort.sort(key = lambda pos: (int(pos[0]), int(pos[1])))
                for q in posSort:
#                     if displaySeq:
                    outFile.write('      contig         {0:>9} - {1:<9} {2}\n'.format(*q))
                    for s in hitsRightDict[q]:
                        outFile.write('      {3:<14} {0:>9} - {1:<9} {2} {4}\n'.format(*s))
                    outFile.write('\n')
#                     else:
#                         outFile.write('      contig         {0:>9} - {1:<9}\n'.format(*q))
#                         for s in hitsRightDict[q]:
#                             outFile.write('      {3:<14} {0:>9} - {1:<9} {4}\n'.format(*s))
#                         outFile.write('\n')

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
        outFile.write('      {0:>18}: {1:>7}  {2:>7.2%}\n'.format(f, totalCounts[f], totalCounts[f] / totalFeatures))
    outFile.close()
    
def drawXML(fileName, outputDir, bestHitsOnly = True):
    """Write XML results to human-readable text file.
    
       Parameters:
           fileName:    absolute path of XML results file
           outFileName: absolute path of output text file
           displaySeq:  whether to write long version with all sequence alignments
           
       Return value: None
    """
    
    FLANK_LIMIT = 3
    chromNames = {'NC_035107.1': 'Chr1', 'NC_035108.1': 'Chr2', 'NC_035109.1': 'Chr3'}
    
    if not Path(outputDir).exists():
        os.mkdir(outputDir)
    else:
        os.system('rm ' + outputDir + '/*.png')
    
    totalCounts = {}
    totalFeatures = 0
    totalContigs = 0
    
    contigFeaturesV = {}
    contigFeaturesA = {}
    contigLengths = {}
    
    tree = ET.parse(fileName)
    root = tree.getroot()
    for contig in root:
    
        if bestHitsOnly and (contig.attrib['besthit'] != 'True'):
            continue
            
        p = re.compile(r'(.*)_length')
        m = p.search(contig.attrib['name'])
        c = m.group(1)
        if c not in contigFeaturesV:
            contigFeaturesV[c] = []
            contigFeaturesA[c] = []
            
        print(contig.attrib['name'])
        p = re.compile(r'length_(.*)_cov')
        m = p.search(contig.attrib['name'])
        contigLength = int(m.group(1))
        contigLengths[c] = contigLength
        features = []
        totalContigs += 1
        for v in contig.findall('virushit'):
            qstart = int(v.find('qstart').text)
            qend = int(v.find('qend').text)
            sstart = v.find('sstart').text
            send = v.find('send').text
            length = abs(int(sstart) - int(send) + 1)
            if send > sstart:
                strandStr = '+'
            else:
                strandStr = '-'
            evalue = float(v.find('evalue').text)
            bitscore = float(v.find('bitscore').text)
            pident = float(v.find('pident').text)
            stitle = v.attrib['stitle']
#            if v.attrib['minevalue'] == 'True':
#            if v.attrib['maxbitscore'] == 'True':
            if not bestHitsOnly and (contig.attrib['besthit'] == 'True'):
                stitle = '**' + stitle.lstrip('|') + '**'
            if qend > qstart:
                strand = 1
            else:
                strand = -1
            features.append(GraphicFeature(start=qstart, end=qend, strand=strand, 
                                           color='#ff0000', fontdict = {'size': 9},
                                           label = stitle.lstrip('|') + ' ' + sstart + '-' + send + ' (' + strandStr + '; ' + str(length) + ' bp; ' + str(pident) + '%; ' + str(evalue) + ')'))
            contigFeaturesV[c].append(copy.deepcopy(features[-1]))
        hitsLeft = contig.findall('vectorhitleft')
        hitsRight = contig.findall('vectorhitright')
        flanks = contig.find('flanks')
        matches = flanks.findall('match')
        hitsDone = []
        if len(matches) > 0:
            matchCount = 1
            for match in matches[:FLANK_LIMIT]:
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
                            sstart = v.find('sstart').text
                            send = v.find('send').text
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
                                                       label = 'M' + str(matchCount) + ': ' + seqid + ' ' + sstart + '-' + send + ' (' + strandStr + ')'))
                            break
                matchCount += 1

        aaCoverage = [0] * contigLength
        for hits in (hitsLeft, hitsRight):
            if len(hits) > 0:
                hitsDict = {}  # consolidate query hit: [subject hits]
                for v in hits:
                    pos = (v.find('qstart').text, v.find('qend').text, v.find('qseq').text)
                    for x in range(int(pos[0]), int(pos[1])):
                        aaCoverage[x] += 1
                    if v.attrib['id'] in hitsDone:
                        continue
                    featureString = '| '
    #                     features = v.findall('feature')
    #                     for f in features:
    # #                         if displaySeq:
    #                         featureString += f.attrib['type'] + '; ' + f.find('id').text + '; ' + f.find('location').text + ' | '
    # #                         else:
    # #                             featureString += f.attrib['type'] + '; ' + f.find('id').text + ' | '
                    if pos not in hitsDict:
                        hitsDict[pos] = [(v.find('sstart').text, v.find('send').text, v.find('sseq').text, v.attrib['seqid'], featureString)]
                    else:
                        hitsDict[pos].append((v.find('sstart').text, v.find('send').text, v.find('sseq').text, v.attrib['seqid'], featureString))
                posSort = list(hitsDict.keys())
                if hits == hitsLeft:
                    posSort.sort(key = lambda pos: (int(pos[1]), int(pos[0])), reverse = True)  # closest first
                else:
                    posSort.sort(key = lambda pos: (int(pos[0]), int(pos[1])))  # closest first
                for q in posSort[:max(0, FLANK_LIMIT-len(matches))]:
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
                    features.append(GraphicFeature(start=qstart, end=qend, strand=strand, 
                                               color='#0000ff', fontdict = {'size': 8},
                                               label = seqid + ' ' + sstart + '-' + send + ' (' + strandStr + ')'))
                
#         
#         featureSummary = contig.find('featuresummary')
#         if featureSummary is not None:
#             outFile.write('   Feature summary:\n')
#             if len(featureSummary) == 1:  # only 'None'
#                 outFile.write('      None found\n\n')
#             else:
#                 total = 0
#                 for f in featureSummary:
#                     if f.tag != 'None':
#                         total += int(f.text)
#                 noneCount = 0
#                 for f in featureSummary:
#                     if f.tag != 'None':
#                         outFile.write('      {0:>18}: {1:>5}  {2:>7.2%}\n'.format(f.tag, f.text, int(f.text) / total))
#                     else:
#                         noneCount = int(f.text)
#                     totalCounts[f.tag] = totalCounts.get(f.tag, 0) + int(f.text)
#                 outFile.write('      Sequences with no features found: ' + str(noneCount) + '\n')
#                 totalFeatures += total
#                 outFile.write('\n')
                
        record = GraphicRecord(sequence_length = contigLength, features = features)
        
        fig, (ax1, ax2) = pyplot.subplots(2, 1, sharex = True, figsize = (10, 6), gridspec_kw = {'height_ratios': [5, 1]})
        
        record.plot(max_label_length = 80, ax = ax1, with_ruler = False)
        ax2.fill_between(range(contigLength), aaCoverage, step = 'mid', alpha = 0.3)
#        ax2.scatter(range(contigLength), aaCoverage)
        ax2.set_ylim(bottom = 0)
        ax2.set_ylabel('Aa hits')

        fig.savefig(Path(outputDir) / (contig.attrib['name'] + '.png'))
        pyplot.close(fig)
        
    if not bestHitsOnly:
        for c in contigFeaturesV:
            record = GraphicRecord(sequence_length = contigLengths[c], features = contigFeaturesV[c])
        
            ax, _ = record.plot(max_label_length = 80, figure_width = 10)
    #        ax2.fill_between(range(contigLength), aaCoverage, step = 'mid', alpha = 0.3)
    #        ax2.set_ylim(bottom = 0)
    #        ax2.set_ylabel('Aa hits')

            ax.figure.savefig(Path(outputDir) / (c + '.png'))
            pyplot.close(ax.figure)
                
#     outFile.write('Overall summary:\n')
#     outFile.write('   Total viral hits: ' + str(totalContigs) + '\n\n')
#     outFile.write('   Feature counts:\n')
#     for f in totalCounts:
#         outFile.write('      {0:>18}: {1:>7}  {2:>7.2%}\n'.format(f, totalCounts[f], totalCounts[f] / totalFeatures))
#     outFile.close()


    
###############################################################################

def reverseComplement(dna):
	'''Return the reverse complement of dna.'''
	
	dna = dna.upper()
	basecomplement = {'A': 'T', 'C': 'G', 'T': 'A', 'G': 'C'} 
	complement = ''
	for base in reversed(dna):
		complement = complement + basecomplement.get(base, base)
	return complement

def consolidate(dirName, virusName, virusLength):
    """Consolidate all results for a particular virus.
       Write a tab-delimited table summarizing hit regions for all specimens.
       Write a FASTA file containing specimen EV regions aligned to viral reference genome
    
       Parameters:
           dirName:     absolute path of directory containing results XML files
           virusName:   name of the virus
           virusLength: length of the virus reference genome
           
       Return value: None
    """
    
    print('Consolidating results in ' + dirName + '...')
    
    dir = Path(dirName).resolve()
    viralHits1 = {}
    viralHits2 = {}
    allHits = []
    viralSeqs = {}
    refSeqs = {}
    for file in dir.iterdir():
        if '_hits.xml' in file.name:
            index = file.name.index('_' + virusName)
            specimen = file.name[:index]
            viralSeqs[specimen] = {}
            refSeqs[specimen] = {}
            tree = ET.parse(str(file))
            root = tree.getroot()
            for contig in root:
                hit = []
                for v in contig.findall('virushit'):
                    hit.extend([int(v.find('sstart').text), int(v.find('send').text)])
                    viralSeqs[specimen][(int(v.find('sstart').text), int(v.find('send').text))] = v.find('qseq').text
                    refSeqs[specimen][(int(v.find('sstart').text), int(v.find('send').text))] = v.find('sseq').text
                if hit[0] > hit[-1]:
                    hit = hit[::-1]
                hit = tuple(hit)
                
                flanks = contig.find('flanks')
                matches = flanks.findall('match')
                if len(matches) > 0:
                    if specimen not in viralHits2:
                        viralHits2[specimen] = []
                    viralHits2[specimen].append(hit)
                else:
                    if specimen not in viralHits1:
                        viralHits1[specimen] = []
                    viralHits1[specimen].append(hit)
                    
                if hit not in allHits:
                    allHits.append(hit)
                        
    for specimen in viralSeqs:
        newDict = {}
        newDictRef = {}
        for (start, end) in viralSeqs[specimen]:
            if start > end:
                newDict[(end, start)] = reverseComplement(viralSeqs[specimen][(start, end)])
                newDictRef[(end, start)] = reverseComplement(refSeqs[specimen][(start, end)])
            else:
                newDict[(start, end)] = viralSeqs[specimen][(start, end)]
                newDictRef[(start, end)] = refSeqs[specimen][(start, end)]
        viralSeqs[specimen] = newDict
        refSeqs[specimen] = newDictRef

#     for specimen in viralSeqs:
#         hits = viralHits1[specimen] + viralHits2[specimen]
#         newDict = {}
#         for hit in hits:
#             if len(hit) == 2:
#                 if hit not in viralSeqs[specimen]:  # reversed
#                     newDict[hit] = reverseComplement(viralSeqs[specimen][(hit[1], hit[0])])
#                 else:
#                     newDict[hit] = viralSeqs[specimen][hit]
#             else:
#                 newDict[hit] = ''
#                 for index in range(0, len(hit), 2):
#                     hit2 = (hit[index], hit[index + 1])
#                     if hit2 not in viralSeqs[specimen]:  # reversed
#                         newDict[hit] += reverseComplement(viralSeqs[specimen][(hit2[1], hit2[0])])
#                     else:
#                         newDict[hit] += viralSeqs[specimen][hit2]
#         viralSeqs[specimen] = newDict
    
#     for specimen in viralSeqs:
#         items = list(viralSeqs[specimen].items())
#         items.sort()
#         seq = ''
#         for item in items:
#             seq = seq + item[1]
#         viralSeqs[specimen] = seq.replace('-', '')   # remove gaps

#     for specimen in viralSeqs:
#         seq = ['-'] * virusLength
#         for (start, end) in viralSeqs[specimen]:
#             seq[start:end+1] = list(viralSeqs[specimen][(start, end)])
#         viralSeqs[specimen] = ''.join(seq)

#     for specimen in viralSeqs:
#         seq = '-' * virusLength
#         for (start, end) in viralSeqs[specimen]:
#             seq = seq[:start - 1] + viralSeqs[specimen][(start, end)] + seq[end:]
#         #assert len(seq) == virusLength
#         viralSeqs[specimen] = seq
        
    refIndices = {}    # specimen seq at each reference index
    for specimen in viralSeqs:
        refIndices[specimen] = {}
        for (start, end) in viralSeqs[specimen]:
            refIndex = start - 2  # -1 for 0-indexing, -1 for initial increment
            seq = viralSeqs[specimen][(start, end)]
            refSeq = refSeqs[specimen][(start, end)]
            for index in range(len(seq)):
                if refSeq[index] != '-':
                    refIndex += 1
                    refIndices[specimen][refIndex] = seq[index]
                else:
                    refIndices[specimen][refIndex] += seq[index]
        viralSeqs[specimen] = ''
                    
    referenceFilename = '/home/havill/data/aegypti/genomes/' + virusName + '.fasta'
    refFile = open(referenceFilename, 'r')
    refRecord = SeqIO.read(refFile, 'fasta')
    refFile.close()
    referenceSeq = str(refRecord.seq)
    referenceID = refRecord.id
    viralSeqs[referenceID] = ''
    for refIndex in range(virusLength):
        maxLength = max([len(refIndices[specimen].get(refIndex, '-')) for specimen in refIndices])
        viralSeqs[referenceID] += '{0:-<{1}}'.format(referenceSeq[refIndex], maxLength)
        for specimen in refIndices:
            viralSeqs[specimen] += '{0:-<{1}}'.format(refIndices[specimen].get(refIndex, ''), maxLength)

    seqOut = open(str(dir / ('sequences_' + virusName + '.fasta')), 'w')
#     referenceFilename = '/home/havill/data/aegypti/genomes/' + virusName + '.fasta'
#     if Path(referenceFilename).exists():
#         refFile = open(referenceFilename, 'r')
#         for line in refFile:
#             seqOut.write(line)
#         refFile.close()
    sortedSpecimens = list(viralSeqs.keys())
    sortedSpecimens.remove(referenceID)
    sortedSpecimens.sort()
    seqOut.write('>' + referenceID + '\n')
    seqOut.write(viralSeqs[referenceID] + '\n')
    for specimen in sortedSpecimens:
        seqOut.write('>' + specimen + '\n')
        seqOut.write(viralSeqs[specimen] + '\n')
    seqOut.close()

    out = open(str(dir / ('results_' + virusName + '.txt')), 'w')
    allHits.sort()                        # write (start, end) header
    hitCounts1 = {}
    hitCounts2 = {}
    out.write('\t')
    for hit in allHits:
        out.write(str(hit) + '\t')
        hitCounts1[hit] = 0
        hitCounts2[hit] = 0
    out.write('\n')
    out.write('\t')                       # write length header
    for hit in allHits:
        length = sum([abs(hit[i+1] - hit[i]) + 1 for i in range(0, len(hit), 2)])
        out.write(str(length) + '\t')
    out.write('\n')
    specimens = list(viralHits1.keys())
    specimens.sort()
    for specimen in specimens:
        out.write(specimen + '\t')
        for hit in allHits:
            if (specimen in viralHits2) and (hit in viralHits2[specimen]):
                out.write('2\t')
                hitCounts2[hit] += 1
            elif (specimen in viralHits1) and (hit in viralHits1[specimen]):
                out.write('1\t')
                hitCounts1[hit] += 1
            else:
                out.write('0\t')
        out.write('\n')
    out.write('Total 2\t')
    for hit in allHits:
        out.write(str(hitCounts2[hit]) + '\t')
    out.write('\n')
    out.write('Total 1\t')
    for hit in allHits:
        out.write(str(hitCounts1[hit]) + '\t')
    out.write('\n')
    out.write('Total all\t')
    for hit in allHits:
        out.write(str(hitCounts1[hit] + hitCounts2[hit]) + '\t')
    out.write('\n')
    out.close()
    

# def consolidateAll(dirName, virusName, bestHitsOnly = True):
#     """Consolidate all results for a particular virus.
#        Write a tab-delimited table summarizing hit regions for all specimens.
#        Write a FASTA file containing specimen EV regions aligned to viral reference genome
#     
#        Parameters:
#            dirName:     absolute path of directory containing results XML files
#            virusName:   name of the virus
#            virusLength: length of the virus reference genome
#            
#        Return value: None
#     """
#     
# #    cp */spades_all/*all_hits.xml HITS_all
#     
#     print('Consolidating results in ' + dirName + '...')
#     
#     dir = Path(dirName).resolve()
#     viralHits1 = {}
#     viralHits2 = {}
#     viralHits = {}  # all together
#     allHits = []
#     viralSeqs = {}
#     refSeqs = {}
#     for file in dir.iterdir():
#         if '_hits.xml' in file.name:
#             index = file.name.index('_' + virusName)
#             specimen = file.name[:index]
#             viralSeqs[specimen] = {}
#             refSeqs[specimen] = {}
#             tree = ET.parse(str(file))
#             root = tree.getroot()
#             for contig in root:   # one per virus
#                 if bestHitsOnly and (contig.attrib['besthit'] != 'True'):
#                     continue
#                 hit = []
#                 for v in contig.findall('virushit'):
#                     hit.extend([int(v.find('sstart').text), int(v.find('send').text)])
#                     viralSeqs[specimen][(int(v.find('sstart').text), int(v.find('send').text))] = v.find('qseq').text
#                     refSeqs[specimen][(int(v.find('sstart').text), int(v.find('send').text))] = v.find('sseq').text
#                 seqid, stitle = v.attrib['seqid'].strip('ref|'), v.attrib['stitle'].lstrip('|')  # should be same for all v
#                 if hit[0] > hit[-1]:
#                     hit = hit[::-1]
# #                hit = tuple(hit)
#                 
# # ALSO APPEND CONTIG ID TO IDENTIFY HITS IN DIFFERENT CONTIGS
# # FOR NOW, MAYBE INSERT 367-3; 1223-1503, 291-426; 367-141; ETC.
# # ALSO, LIAO NING MAYBE 75%+
# 
#                 if specimen not in viralHits:
#                     viralHits[specimen] = []
#                 viralHits[specimen].append([seqid, stitle] + hit)
#                 
#                 flanks = contig.find('flanks')
#                 matches = flanks.findall('match')
#                 if len(matches) > 0:
#                     if specimen not in viralHits2:
#                         viralHits2[specimen] = []
#                     viralHits2[specimen].append([seqid, stitle] + hit)
#                 else:
#                     if specimen not in viralHits1:
#                         viralHits1[specimen] = []
#                     viralHits1[specimen].append([seqid, stitle] + hit)
#                     
#                 if [seqid, stitle] + hit not in allHits:
#                     allHits.append([seqid, stitle] + hit)
#                         
# #     for specimen in viralSeqs:
# #         newDict = {}
# #         newDictRef = {}
# #         for (start, end) in viralSeqs[specimen]:
# #             if start > end:
# #                 newDict[(end, start)] = reverseComplement(viralSeqs[specimen][(start, end)])
# #                 newDictRef[(end, start)] = reverseComplement(refSeqs[specimen][(start, end)])
# #             else:
# #                 newDict[(start, end)] = viralSeqs[specimen][(start, end)]
# #                 newDictRef[(start, end)] = refSeqs[specimen][(start, end)]
# #         viralSeqs[specimen] = newDict
# #         refSeqs[specimen] = newDictRef
# #         
# #     refIndices = {}    # specimen seq at each reference index
# #     for specimen in viralSeqs:
# #         refIndices[specimen] = {}
# #         for (start, end) in viralSeqs[specimen]:
# #             refIndex = start - 2  # -1 for 0-indexing, -1 for initial increment
# #             seq = viralSeqs[specimen][(start, end)]
# #             refSeq = refSeqs[specimen][(start, end)]
# #             for index in range(len(seq)):
# #                 if refSeq[index] != '-':
# #                     refIndex += 1
# #                     refIndices[specimen][refIndex] = seq[index]
# #                 else:
# #                     refIndices[specimen][refIndex] += seq[index]
# #         viralSeqs[specimen] = ''
# #                     
# #     referenceFilename = '/home/havill/data/aegypti/genomes/' + virusName + '.fasta'
# #     refFile = open(referenceFilename, 'r')
# #     refRecord = SeqIO.read(refFile, 'fasta')
# #     refFile.close()
# #     referenceSeq = str(refRecord.seq)
# #     referenceID = refRecord.id
# #     viralSeqs[referenceID] = ''
# #     for refIndex in range(virusLength):
# #         maxLength = max([len(refIndices[specimen].get(refIndex, '-')) for specimen in refIndices])
# #         viralSeqs[referenceID] += '{0:-<{1}}'.format(referenceSeq[refIndex], maxLength)
# #         for specimen in refIndices:
# #             viralSeqs[specimen] += '{0:-<{1}}'.format(refIndices[specimen].get(refIndex, ''), maxLength)
# # 
# #     seqOut = open(str(dir / ('sequences_' + virusName + '.fasta')), 'w')
# #     sortedSpecimens = list(viralSeqs.keys())
# #     sortedSpecimens.remove(referenceID)
# #     sortedSpecimens.sort()
# #     seqOut.write('>' + referenceID + '\n')
# #     seqOut.write(viralSeqs[referenceID] + '\n')
# #     for specimen in sortedSpecimens:
# #         seqOut.write('>' + specimen + '\n')
# #         seqOut.write(viralSeqs[specimen] + '\n')
# #     seqOut.close()
# 
#     out = open(str(dir / ('results_' + virusName + '.txt')), 'w')
#     
#     allViruses = []
#     for entry in allHits:
#         allViruses.append((entry[1], entry[0]))  # stitle, seqid
#     allViruses = list(set(allViruses))
#     allViruses.sort()
# #    allHits.sort()                        # write (start, end) header
# #     hitCounts1 = {}
# #     hitCounts2 = {}
#     out.write('\t')
#     for stitle, seqid in allViruses:
#         out.write(seqid + '\t')
#     out.write('\n')
#     out.write('\t')
#     for stitle, seqid in allViruses:
#         out.write(stitle + '\t')
#     out.write('\n')
#     
# #     for hit in allHits:
# #         out.write(str(hit) + '\t')
# #         hitCounts1[hit] = 0
# #         hitCounts2[hit] = 0
# #     out.write('\n')
# #     out.write('\t')                       # write length header
# #     for hit in allHits:
# #         length = sum([abs(hit[i+1] - hit[i]) + 1 for i in range(0, len(hit), 2)])
# #         out.write(str(length) + '\t')
# #     out.write('\n')
#     specimens = list(viralHits1.keys())
#     specimens.sort()
#     for specimen in specimens:
#         out.write(specimen + '\t')
#         for stitle, seqid in allViruses:
#             length = 0
#             for hit in viralHits[specimen]:
#                 if hit[0] == seqid:
#                     length += sum([abs(hit[i+1] - hit[i]) + 1 for i in range(2, len(hit), 2)])
#             out.write(str(length) + '\t')
#                     
#             # if (specimen in viralHits2) and (hit in viralHits2[specimen]):
# #                 out.write('2\t')
# #                 hitCounts2[hit] += 1
# #             elif (specimen in viralHits1) and (hit in viralHits1[specimen]):
# #                 out.write('1\t')
# #                 hitCounts1[hit] += 1
# #             else:
# #                 out.write('0\t')
#         out.write('\n')
# #     out.write('Total 2\t')
# #     for hit in allHits:
# #         out.write(str(hitCounts2[hit]) + '\t')
# #     out.write('\n')
# #     out.write('Total 1\t')
# #     for hit in allHits:
# #         out.write(str(hitCounts1[hit]) + '\t')
# #     out.write('\n')
# #     out.write('Total all\t')
# #     for hit in allHits:
# #         out.write(str(hitCounts1[hit] + hitCounts2[hit]) + '\t')
# #     out.write('\n')
#     out.close()

def consolidateAll(dirName, virusName, bestHitsOnly = True):
    """Consolidate all results for a particular virus.
       Write a tab-delimited table summarizing hit regions for all specimens.
       Write a FASTA file containing specimen EV regions aligned to viral reference genome
    
       Parameters:
           dirName:     absolute path of directory containing results XML files
           virusName:   name of the virus
           virusLength: length of the virus reference genome
           
       Return value: None
    """
    
#    cp */spades_all/*all_hits.xml ../results/HITS_all
    
    print('Consolidating results in ' + dirName + '...')
    
#     dirList = [dir for dir in Path(ROOT).iterdir() if dir.is_dir() and 'combined' not in dir.name]
#     for dir in dirList:
#         spadesPath = dir / ('spades_' + virusName)
#         xmlFilename = spadesPath / (Path(dir).name + '_' + virusName + '_hits.xml')
#         if xmlFilename.exists():
    
    dir = Path(dirName).resolve()
    viralHits1 = {}
    viralHits2 = {}
    viralHits = {}  # all together
    allHits = []
    viralSeqs = {}
    refSeqs = {}
    for file in dir.iterdir():
        if '_hits.xml' in file.name:
            index = file.name.index('_' + virusName)
            specimen = file.name[:index]
            virus_fasta = open(str(file)[:-4] + '.fasta', 'w')
#             viralSeqs[specimen] = {}
#             refSeqs[specimen] = {}
            tree = ET.parse(str(file))
            root = tree.getroot()
            for contig in root:   # one per virus
                if bestHitsOnly and (contig.attrib['besthit'] != 'True'):
                    continue
                hit = {}
                v = contig.find('virushit')
                seqid, stitle = v.attrib['seqid'].strip('ref|'), v.attrib['stitle'].lstrip('|')  # should be same for all v
                if (seqid, stitle) not in viralSeqs:
                    viralSeqs[(seqid, stitle)] = {}
                    refSeqs[(seqid, stitle)] = {}
                viralSeqs[(seqid, stitle)][specimen] = {}
                refSeqs[(seqid, stitle)][specimen] = {}
                for v in contig.findall('virushit'):
                    hit[(int(v.find('qstart').text), int(v.find('qend').text))] = (int(v.find('sstart').text), int(v.find('send').text))
                    viralSeqs[(seqid, stitle)][specimen][(int(v.find('sstart').text), int(v.find('send').text))] = v.find('qseq').text
                    refSeqs[(seqid, stitle)][specimen][(int(v.find('sstart').text), int(v.find('send').text))] = v.find('sseq').text
                    virus_fasta.write('>' + seqid + ' ' + stitle + ' ' + v.find('sstart').text + '-' + v.find('send').text + ' (host)\n')
                    virus_fasta.write(v.find('qseq').text + '\n')
                    virus_fasta.write('>' + seqid + ' ' + stitle + ' ' + v.find('sstart').text + '-' + v.find('send').text + ' (reference)\n')
                    virus_fasta.write(v.find('sseq').text + '\n')
#                seqid, stitle = v.attrib['seqid'].strip('ref|'), v.attrib['stitle'].lstrip('|')  # should be same for all v
#                 if hit[0] > hit[-1]:
#                     hit = hit[::-1]
#                hit = tuple(hit)

                if specimen not in viralHits:
                    viralHits[specimen] = []
                viralHits[specimen].append((seqid, stitle, hit))
                
                flanks = contig.find('flanks')
                matches = flanks.findall('match')
                numVectorHits = len(contig.findall('vectorhitleft') + contig.findall('vectorhitright'))
                if len(matches) > 0:
                    if specimen not in viralHits2:
                        viralHits2[specimen] = []
                    viralHits2[specimen].append((seqid, stitle, hit))
                elif numVectorHits > 0:
                    if specimen not in viralHits1:
                        viralHits1[specimen] = []
                    viralHits1[specimen].append((seqid, stitle, hit))
                    
                if (seqid, stitle, hit) not in allHits:
                    allHits.append((seqid, stitle, hit))
            virus_fasta.close()
#    
    for (seqid, stitle) in viralSeqs:
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
                print('Error getting reference sequence ' + seqid + ' from NCBI.')
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

        seqOut = open(str(dir / ('sequences_' + seqid + '.fasta')), 'w')
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
    out = open(str(dir / ('results_' + virusName + '.txt')), 'w')
    
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
    specimens.sort()
    for specimen in specimens:
        out.write(specimen + '\t')
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
                    if (specimen in viralHits2) and (hit in viralHits2[specimen]):
                        s.append(2)
                    elif (specimen in viralHits1) and (hit in viralHits1[specimen]):
                        s.append(1)
                    else:
                        s.append(0)
                    rangeList.append(s)  # s = (sstart, send, sstart, send, ..., 0/1/2)
                    
            rangeList.sort()    # among all hits, sort by first sstart
            rangeString = ''
            for s in rangeList:
                for i in range(0, len(s) - 1, 2):
                    rangeString += str(s[i]) + '-' + str(s[i+1]) + ', '
                rangeString = rangeString[:-2] + ('*' * s[-1]) + ' | '
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
        out.write('\n')
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

def assembleAnalyze(dirName, virusName, virusDB, minLength, startStep, alsoFindFeatures):
    """Assemble reads and identify putative viral insertions and features.
    
       Parameters:
           dirName:      absolute path of parent of directory containing fastq files to assemble 
           virusName:    string containing the name of the virus
           virusDB:      string containing name of BLAST virus database
           startStep:    step at which to start
           findFeatures: whether to locate features
           
       Return value: None
    """
    
    if startStep <= 1:
        success = assembleReads(dirName, virusName)  #1 = assemble
    else:
        success = True

    if success:
        if startStep <= 2:
            blastScaffolds(dirName, virusName, virusDB)  #2 = blast scaffolds
            
        hitsDirName = str(Path(dirName).parent.parent / ('results/HITS_' + virusName))
        hitsCombinedDirName = str(Path(dirName).parent.parent / ('results/HITS_' + virusName + '_combined'))
        
        if startStep <= 3:
            xmlFilename = getHits(dirName, virusName, minLength)  #3 = get hits
            if xmlFilename is None:
                print('Assembly failed.  No results found.')
                print('Done.')
                return
            
            displayXML(xmlFilename, xmlFilename[:-4] + '.txt', True)
            displayXML(xmlFilename, xmlFilename[:-4] + '_short.txt', False)
            
            print('Copying hits to ' + hitsDirName + '...')
            if not Path(hitsDirName).exists():
                os.mkdir(hitsDirName)
            os.system('cp -v ' + xmlFilename + ' ' + hitsDirName)

            if 'combined' in dirName:
                if not Path(hitsCombinedDirName).exists():
                    os.mkdir(hitsCombinedDirName)
                os.system('cp -v ' + xmlFilename + ' ' + hitsCombinedDirName)
        else:
            spadesPath = Path(dirName) / ('spades_' + virusName)
            xmlFilename = str(spadesPath / (Path(dirName).name + '_' + virusName + '_hits.xml'))
    
        if (startStep <= 4) and alsoFindFeatures and Path(xmlFilename).exists():
            global iTrees
            if iTrees is None:
                iTrees = makeIntervalTrees(str(Path(dirName).parent))
            xmlFeaturesFilename = findFeatures(xmlFilename, iTrees)  #4 = get features
    
            displayXML(xmlFeaturesFilename, xmlFeaturesFilename[:-4] + '.txt', True)
            displayXML(xmlFeaturesFilename, xmlFeaturesFilename[:-4] + '_short.txt', False)
        
            print('Copying hits/features to ' + hitsDirName + '...')
            if 'combined' in dirName:
                if not Path(hitsCombinedDirName).exists():
                    os.mkdir(hitsCombinedDirName)
                os.system('cp -v ' + xmlFeaturesFilename + ' ' + hitsCombinedDirName)
            else:
                if not Path(hitsDirName).exists():
                    os.mkdir(hitsDirName)
                os.system('cp -v ' + xmlFeaturesFilename + ' ' + hitsDirName)

    else:
        print('Assembly failed.  No results found.')
        
    print('Done.')
    

def doIt(dirName, filenameBAM, virusName, virusDB, minLength, alsoFindFeatures):
    
    filenameBAM = str(Path(dirName) / filenameBAM)
    newFilenameFASTA = getUnmappedReads(filenameBAM)
    
    if virusName is not None:
        blastOutputFileNameCSV = blastViral(newFilenameFASTA, virusName, virusDB)

        viralFilenameBAM = getReads(blastOutputFileNameCSV, filenameBAM, True, virusName)
    
        writeReadsFASTQ(viralFilenameBAM, virusName)
    
        assembleAnalyze(dirName, virusName, virusDB, minLength, 1, alsoFindFeatures)
        

# def process_cmdline():
#     if (sys.argv[1] == '-a'):
#         if len(sys.argv) >= 3:
#             pattern = sys.argv[2]
#         else:
#             pattern = ''
#         dir = Path('.').resolve()   # ~/data/aegypti/analyzed/
#         # workers = []
#         # for subdir in dir.iterdir():
#         #   if subdir.is_dir() and subdir.name[:len(pattern)] == pattern and \
#         #     not Path(subdir / ('spades_' + VIRUS_NAME)).exists():
#         #       for file in subdir.iterdir():
#         #           if not file.is_dir() and (file.name[-26:] == '.sorted.deduped.merged.bam'):
#         #               p = Process(target = doit, args = (str(subdir.resolve()), file.name, VIRUS_NAME, VIRUS_DB))
#         #               workers.append(p)
#         #               p.start()
#         #               break
#         #
#         # for p in workers:
#         #   p.join()
# 
#         print('Consolidating results...')
#         os.chdir(dir)
#         if not Path('HITS_' + VIRUS_NAME).exists():
#             os.mkdir('HITS_' + VIRUS_NAME)
#         #os.system('cp -v */spades_' + VIRUS_NAME + '/*viral_hits* ' + 'HITS_' + VIRUS_NAME)
#         consolidate(VIRUS_NAME)
#     elif len(sys.argv) == 2:
#         filenameBAM = sys.argv[1]
#         doit(Path('.').resolve(), filenameBAM)
        
# def doAllGetUnmappedReads():
#     dirList = [dir for dir in Path(ROOT).iterdir() if 'combined' not in dir.name]
#     workers = []
#     for dir in dirList:
#         if dir.is_dir():
#             for file in dir.iterdir():
#                 if not file.is_dir() and (file.name[-26:] == '.sorted.deduped.merged.bam'):
#                     p = Process(target = getUnmappedReads, args = (str(file.resolve()),))
#                     workers.append(p)
#                     p.start()
#                     break
#         
#     for p in workers:
#         p.join()
#         
# doAllGetUnmappedReads()
# exit()

def doItAllViruses(dirName, filenameBAM, virusName, virusDB):
    
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
            blastScaffolds(dirName, virusName, virusDB)  #2 = blast scaffolds
            
#             hitsDirName = str(Path(dirName).parent.parent / ('results/HITS_' + virusName))
#             hitsCombinedDirName = str(Path(dirName).parent.parent / ('results/HITS_' + virusName + '_combined'))
        
            xmlFilename = getHits(dirName, virusName, 100)  #3 = get hits
            if xmlFilename is None:
                print('*** ' + str(Path(dirName).name) + ': assembly failed; no results found. ***')
                return
            
            drawXML(xmlFilename, str(Path(xmlFilename).parent / 'png'), False)
#            displayXML(xmlFilename, xmlFilename[:-4] + '.txt', True)
#            displayXML(xmlFilename, xmlFilename[:-4] + '_short.txt', False)
            
#             print('Copying hits to ' + hitsDirName + '...')
#             if not Path(hitsDirName).exists():
#                 os.mkdir(hitsDirName)
#             os.system('cp -v ' + xmlFilename + ' ' + hitsDirName)
                
#             if alsoFindFeatures and Path(xmlFilename).exists():
#                 global iTrees
#                 if iTrees is None:
#                     iTrees = makeIntervalTrees(str(Path(dirName).parent))
#                 xmlFeaturesFilename = findFeatures(xmlFilename, iTrees)  #4 = get features
#     
#                 displayXML(xmlFeaturesFilename, xmlFeaturesFilename[:-4] + '.txt', True)
#                 displayXML(xmlFeaturesFilename, xmlFeaturesFilename[:-4] + '_short.txt', False)
#         
#                 print('Copying hits/features to ' + hitsDirName + '...')
#                 if not Path(hitsDirName).exists():
#                     os.mkdir(hitsDirName)
#                 os.system('cp -v ' + xmlFeaturesFilename + ' ' + hitsDirName)

        else:
            print('*** ' + str(Path(dirName).name) + ': assembly failed; no results found. ***')
        
    print('*** ' + str(Path(dirName).name) + ' done. ***')

#dir = '/home/havill/data/aegypti/analyzed/El_Dorado_Argentina_U_11.LIN210A2770'
#bamFile = 'Debug019_aegypti_El_Dorado_Argentina_U_11.LIN210A2770.sorted.deduped.merged.bam'

# dirs = []
# dirs.append(('/home/havill/data/aegypti/analyzed/La_Lope-Gabon-10.LIN210A1646',
# 'La_Lope-Gabon-10.LIN210A1646.sorted.deduped.merged.bam'))
# 
# dirs.append(('/home/havill/data/aegypti/analyzed/Paqueta-Brazil-1.LIN210A2127',
# 'Debug05-10162017-aegy-Paqueta-Brazil-1.LIN210A2127.sorted.deduped.merged.bam'))
# 
# dirs.append(('/home/havill/data/aegypti/analyzed/Cairns-Australia-1.LIN210A2151',
# 'Debug05-10162017-aegy-Cairns-Australia-1.LIN210A2151.sorted.deduped.merged.bam'))
# 
# for dir, bamFile in dirs:
#     doItAllViruses(dir, bamFile, 'all', 'virusdb')

### Redo /Volumes/Data/aegypti/analyzed/Maricopa_County_Arizona-USA-20.LIN210A1615, /Volumes/Data/aegypti/analyzed/La_Lope-Gabon-11.LIN210A1647

# dir = '/home/havill/data/aegypti/analyzed/Bangkok_Thailand_12.LIN210A1690'
# blastScaffolds(dir, 'all', 'virusdb')
# xmlFilename = getHits(dir, 'all', 100)
# xmlFilename = '/home/havill/data/aegypti/analyzed/La_Lope-Gabon-10.LIN210A1646/spades_all/La_Lope-Gabon-10.LIN210A1646_all_hits.xml'
# xmlFilename = '/home/havill/data/aegypti/analyzed/Paqueta-Brazil-1.LIN210A2127/spades_all/Paqueta-Brazil-1.LIN210A2127_all_hits.xml'
# xmlFilename = '/home/havill/data/aegypti/analyzed/Bangkok_Thailand_12.LIN210A1690/spades_all/Bangkok_Thailand_12.LIN210A1690_all_hits.xml'
# drawXML(xmlFilename, str(Path(xmlFilename).parent / 'png'), False)
# exit()
# displayXML(xmlFilename, xmlFilename[:-4] + '.txt', True)
# displayXML(xmlFilename, xmlFilename[:-4] + '_short.txt', False)
consolidateAll('/home/havill/data/aegypti/results/HITS_all', 'all')
exit()

def checkSuccess():
    dirList = [dir for dir in Path(ROOT).iterdir() if 'combined' not in dir.name]
    dirList =  ['La_Lope_Gabon-_22.2-116710',
                'Skukusa_South_Africa_09.LIN210A1747',
                'Amacuzac_Mexico-_2.2-116616',
                'Amacuzac_Mexico-_21.2-116635',
                'Cairns-Australia-16.LIN210A2166',
                'Cairns-Australia-9.LIN210A2159',
                'El_Dorado_US_U_31',
                'Maricopa_County_AZ-_11.2-116601',
                'Maricopa_County_AZ-_4.2-116544',
                'Amacuzac_Mexico-_12.2-116626',
                'Maricopa_County_AZ-_14.2-116604',
                'Cebu_City_Philippines-_8.2-116720',
                'Cebu_City_Philippines-_17.2-116729',
                'Cebu_City_Philippines-_5.2-116717',
                'El_Dorado_US_U_28',
                'Maricopa_County_AZ-_1.2-116541',
                'Paqueta-Brazil-5.LIN210A2131',
                'La_Lope_Gabon-_21.2-116709',
                'Paqueta-Brazil-9.LIN210A2135',
                'Amacuzac_Mexico-_9.2-116623']
    for dir in dirList:
        dir = Path(ROOT) / dir
        if dir.is_dir():
            spadesDir = dir / 'spades_all'
            scaffoldsFilename = dir / 'spades_all/scaffolds.fasta'
            print(dir.name, end = '\t')
            if spadesDir.exists():
                if scaffoldsFilename.exists():
                    print('assembly OK')
                else:
                    print('assembly failed')
            else:
                print('no spades dir')

#checkSuccess()
#exit()

def rerun():
    dirList = [dir for dir in Path(ROOT).iterdir() if 'combined' not in dir.name]
#    dirList = [Path('/home/havill/data/aegypti/analyzed/Bangkok_Thailand_09.LIN210A1687')]
    for dir in dirList:
        if dir.is_dir():
            blastResultsFilename = dir / 'spades_all/blast_scaffolds_all.csv'
            scaffoldsFilename = dir / 'spades_all/scaffolds.fasta'
            print(dir.name, end = ': ')
            if blastResultsFilename.exists():
                print('re-blasting contigs')
                scaffoldRecords = SeqIO.index(str(scaffoldsFilename), 'fasta')
                hits = readCSV(str(blastResultsFilename), 0)
                for contig in hits:
                    SeqIO.write(scaffoldRecords[contig], '/tmp/query_temp.fasta', 'fasta')
                    os.system('blastn -query /tmp/query_temp.fasta -db virusdb -task blastn -evalue 1e-10' +
                              ' -outfmt "10 qseqid qstart qend qseq sstart send sseq evalue bitscore sseqid stitle pident" >> /tmp/blast_scaffolds_all_temp.csv')
                os.system('mv ' + str(blastResultsFilename) + ' ' + str(blastResultsFilename)[:-4] + '_backup.csv')
                os.system('mv /tmp/blast_scaffolds_all_temp.csv ' + str(blastResultsFilename))
                
#                 xmlFilename = getHits(dir, 'all', 100)
#                 drawXML(xmlFilename, str(Path(xmlFilename).parent / 'png'), False)
            else:
                print('skipping')
                
# rerun()
# exit()

def redoHits():
#    dirList = [dir for dir in Path(ROOT).iterdir() if 'combined' not in dir.name]
#    dirList = [Path(ROOT) / 'El_Dorado_US_U_28']
#     dirNames = ['Cairns-Australia-9.LIN210A2159',
#                 'Maricopa_County_Arizona-USA-23.LIN210A1617',
#                 'Cuanda_Angola_10.LIN210A1728',
#                 'La_Lope-Gabon-24.LIN210A1658',
#                 'Paqueta-Brazil-8.LIN210A2134',
#                 'El_Dorado_US_U_31',
#                 'Cairns-Australia-12.LIN210A2162',
#                 'Tahiti_FrenchPolynesia_02.LIN210A1700',
#                 'Cuanda_Angola_08.LIN210A1726',
#                 'featuretrees.pickle.bak',
#                 'Paqueta-Brazil-11.LIN210A2137',
#                 'Amacuzac_Mexico-_2.2-116616',
#                 'Amacuzac-Mexico-16.LIN210A1629',
#                 'La_Lope-Gabon-9.LIN210A1645',
#                 'HoChiMin_Vietnam_12.LIN210A1670',
#                 'La_Lope-Gabon-12.LIN210A1648',
#                 'Paqueta-Brazil-13.LIN210A2139',
#                 'Maricopa_County_AZ-_11.2-116601',
#                 'Cairns-Australia-21.LIN210A2171',
#                 'Amacuzac-Mexico-20.LIN210A1633',
#                 'Cebu_City_Philippines-_16.2-116728',
#                 'La_Lope-Gabon-5.LIN210A1641',
#                 'Tahiti_FrenchPolynesia_08.LIN210A1706',
#                 'Tahiti_FrenchPolynesia_09.LIN210A1707',
#                 'Bangkok_Thailand_20.LIN210A1698',
#                 'Cuanda_Angola_20.LIN210A1738',
#                 'Maricopa_County_AZ-_4.2-116544',
#                 'Cebu_City_Philippines-_10.2-116722',
#                 'Cairns-Australia-19.LIN210A2169',
#                 'Amacuzac_Mexico-_21.2-116635',
#                 'HoChiMin_Vietnam_03.LIN210A1661',
#                 'Maricopa_County_AZ-_1.2-116541',
#                 'Bangkok_Thailand_17.LIN210A1695',
#                 'Cebu_City_Philippines-_14.2-116726',
#                 'Bangkok_Thailand_07.LIN210A1685',
#                 'Cebu_City_Philippines-_9.2-116721',
#                 'Bangkok_Thailand_03.LIN210A1681',
#                 'La_Lope-Gabon-14.LIN210A1650',
#                 'Amacuzac-Mexico-13.LIN210A1627',
#                 'Cuanda_Angola_03.LIN210A1721',
#                 'Skukusa_South_Africa_08.LIN210A1746',
#                 'Maricopa_County_AZ-_18.2-116608',
#                 'HoChiMin_Vietnam_01.LIN210A1659',
#                 'Cairns-Australia-15.LIN210A2165',
#                 'La_Lope-Gabon-2.LIN210A1638',
#                 'Skukusa_South_Africa_12.LIN210A1750',
#                 'Tahiti_FrenchPolynesia_18.LIN210A1716',
#                 'Cebu_City_Philippines-_11.2-116723',
#                 'Maricopa_County_AZ-_12.2-116602',
#                 'Tahiti_FrenchPolynesia_04.LIN210A1702',
#                 'El_Dorado_Argentina_U_3.LIN210A2706',
#                 'Cairns-Australia-1.LIN210A2151',
#                 'HoChiMin_Vietnam_06.LIN210A1664',
#                 'Amacuzac_Mexico-_12.2-116626',
#                 'La_Lope-Gabon-18.LIN210A1654',
#                 'Cairns-Australia-22.LIN210A2172',
#                 'La_Lope-Gabon-3.LIN210A1639',
#                 'Cairns-Australia-7.LIN210A2157',
#                 'Cebu_City_Philippines-_20.2-116732',
#                 'Paqueta-Brazil-19.LIN210A2145',
#                 'Skukusa_South_Africa_13.LIN210A1751',
#                 'Maricopa_County_AZ-_8.2-116548',
#                 'Skukusa_South_Africa_14.LIN210A1752',
#                 'Tahiti_FrenchPolynesia_14.LIN210A1712',
#                 'HoChiMin_Vietnam_09.LIN210A1667',
#                 'Amacuzac-Mexico-6.LIN210A1622',
#                 'Cebu_City_Philippines-_4.2-116716',
#                 'Maricopa_County_Arizona-USA-10.LIN210A1614',
#                 'Cuanda_Angola_14.LIN210A1732',
#                 'Amacuzac-Mexico-24.LIN210A1636',
#                 'Amacuzac-Mexico-3.LIN210A1619',
#                 'Cairns-Australia-17.LIN210A2167',
#                 'Cebu_City_Philippines-_21.2-116733',
#                 'Amacuzac-Mexico-1.LIN210A1618',
#                 'Cairns-Australia-24.LIN210A2174',
#                 'Maricopa_County_Arizona-USA-5.LIN210A1612',
#                 'Tahiti_FrenchPolynesia_05.LIN210A1703',
#                 'El_Dorado_Argentina_U_11.LIN210A2770',
#                 'Maricopa_County_Arizona-USA-7.LIN210A1613',
#                 'Paqueta-Brazil-10.LIN210A2136',
#                 'Cairns-Australia-18.LIN210A2168',
#                 'Tahiti_FrenchPolynesia_10.LIN210A1708',
#                 'Bangkok_Thailand_04.LIN210A1682',
#                 'Maricopa_County_AZ-_14.2-116604',
#                 'La_Lope-Gabon-20.LIN210A1655',
#                 'Cairns-Australia-20.LIN210A2170',
#                 'La_Lope_Gabon-_19.2-116707',
#                 'La_Lope-Gabon-10.LIN210A1646',
#                 'Maricopa_County_Arizona-USA-2.LIN210A1611',
#                 'Maricopa_County_AZ-_13.2-116603',
#                 'Tahiti_FrenchPolynesia_03.LIN210A1701',
#                 'Paqueta-Brazil-15.LIN210A2141',
#                 'Amacuzac_Mexico-_15.2-116629',
#                 'Cebu_City_Philippines-_19.2-116731',
#                 'Paqueta-Brazil-5.LIN210A2131',
#                 'Amacuzac-Mexico-17.LIN210A1630',
#                 'Cebu_City_Philippines-_8.2-116720',
#                 'Skukusa_South_Africa_15.LIN210A1753',
#                 'Maricopa_County_AZ-_16.2-116606',
#                 'Cairns-Australia-10.LIN210A2160',
#                 'HoChiMin_Vietnam_20.LIN210A1678',
#                 'HoChiMin_Vietnam_13.LIN210A1671',
#                 'Maricopa_County_Arizona-USA-22.LIN210A1616',
#                 'Cuanda_Angola_18.LIN210A1736',
#                 'Amacuzac-Mexico-11.LIN210A1626',
#                 'La_Lope_Gabon-_21.2-116709',
#                 'La_Lope-Gabon-16.LIN210A1652',
#                 'Paqueta-Brazil-9.LIN210A2135',
#                 'Bangkok_Thailand_18.LIN210A1696',
#                 'Cebu_City_Philippines-_17.2-116729',
#                 'Skukusa_South_Africa_03.LIN210A1741',
#                 'Cuanda_Angola_11.LIN210A1729',
#                 'Cairns-Australia-13.LIN210A2163',
#                 'La_Lope-Gabon-7.LIN210A1643',
#                 'Cairns-Australia-2.LIN210A2152',
#                 'Cebu_City_Philippines-_23.2-116735',
#                 'Cebu_City_Philippines-_6.2-116718',
#                 'Paqueta-Brazil-3.LIN210A2129',
#                 'Amacuzac-Mexico-7.LIN210A1623',
#                 'Maricopa_County_AZ-_24.2-116614',
#                 'HoChiMin_Vietnam_08.LIN210A1666',
#                 'Paqueta-Brazil-12.LIN210A2138',
#                 'HoChiMin_Vietnam_05.LIN210A1663',
#                 'La_Lope-Gabon-17.LIN210A1653',
#                 'Maricopa_County_AZ-_6.2-116546',
#                 'Amacuzac-Mexico-10.LIN210A1625',
#                 'Bangkok_Thailand_14.LIN210A1692',
#                 'Skukusa_South_Africa_07.LIN210A1745',
#                 'HoChiMin_Vietnam_02.LIN210A1660',
#                 'Tahiti_FrenchPolynesia_16.LIN210A1714',
#                 'Amacuzac_Mexico-_9.2-116623',
#                 'Bangkok_Thailand_02.LIN210A1680',
#                 'Bangkok_Thailand_13.LIN210A1691',
#                 'Bangkok_Thailand_11.LIN210A1689',
#                 'featuretrees.pickle',
#                 'Bangkok_Thailand_16.LIN210A1694',
#                 'HoChiMin_Vietnam_11.LIN210A1669',
#                 'Cebu_City_Philippines-_24.2-116736',
#                 'Maricopa_County_AZ-_3.2-116543',
#                 'Cebu_City_Philippines-_5.2-116717',
#                 'Cebu_City_Philippines-_3.2-116715',
#                 'Cebu_City_Philippines-_22.2-116734',
#                 'La_Lope-Gabon-13.LIN210A1649',
#                 'Cebu_City_Philippines-_13.2-116725',
#                 'Cairns-Australia-5.LIN210A2155',
#                 'Bangkok_Thailand_10.LIN210A1688',
#                 'Amacuzac-Mexico-14.LIN210A1628',
#                 'Tahiti_FrenchPolynesia_12.LIN210A1710',
#                 'HoChiMin_Vietnam_15.LIN210A1673',
#                 'Bangkok_Thailand_08.LIN210A1686',
#                 'HoChiMin_Vietnam_14.LIN210A1672',
#                 'Cairns-Australia-23.LIN210A2173',
#                 'Tahiti_FrenchPolynesia_06.LIN210A1704',
#                 'Cuanda_Angola_12.LIN210A1730',
#                 'Cuanda_Angola_16.LIN210A1734',
#                 'Paqueta-Brazil-4.LIN210A2130',
#                 'Amacuzac-Mexico-23.LIN210A1635',
#                 'Paqueta-Brazil-1.LIN210A2127',
#                 'Amacuzac-Mexico-19.LIN210A1632',
#                 'Paqueta-Brazil-18.LIN210A2144',
#                 'Tahiti_FrenchPolynesia_15.LIN210A1713',
#                 'Cuanda_Angola_19.LIN210A1737',
#                 'Bangkok_Thailand_01.LIN210A1679',
#                 'Cuanda_Angola_04.LIN210A1722',
#                 'Tahiti_FrenchPolynesia_20.LIN210A1718',
#                 'Skukusa_South_Africa_05.LIN210A1743',
#                 'HoChiMin_Vietnam_04.LIN210A1662',
#                 'La_Lope-Gabon-8.LIN210A1644',
#                 'Bangkok_Thailand_19.LIN210A1697',
#                 'Cebu_City_Philippines-_7.2-116719',
#                 'Cebu_City_Philippines-_12.2-116724',
#                 'Amacuzac-Mexico-22.LIN210A1634',
#                 'Paqueta-Brazil-24.LIN210A2150',
#                 'La_Lope-Gabon-6.LIN210A1642',
#                 'Paqueta-Brazil-7.LIN210A2133',
#                 'Bangkok_Thailand_06.LIN210A1684',
#                 'Maricopa_County_AZ-_15.2-116605',
#                 'Skukusa_South_Africa_11.LIN210A1749',
#                 'Cairns-Australia-11.LIN210A2161',
#                 'Cuanda_Angola_02.LIN210A1720',
#                 'Paqueta-Brazil-14.LIN210A2140',
#                 'Maricopa_County_AZ-_9.2-116549',
#                 'Bangkok_Thailand_05.LIN210A1683',
#                 'La_Lope-Gabon-23.LIN210A1656',
#                 'HoChiMin_Vietnam_07.LIN210A1665']

    dirNames = ['La_Lope-Gabon-11.LIN210A1647',
                'Maricopa_County_Arizona-USA-20.LIN210A1615',
                'Maricopa_County_AZ-_19.2-116609',
                'El_Dorado_Argentina_U_2.LIN210A2698',
                'Tahiti_FrenchPolynesia_13.LIN210A1711',
                'Tahiti_FrenchPolynesia_19.LIN210A1717',
                'Cairns-Australia-8.LIN210A2158',
                'La_Lope_Gabon-_22.2-116710',
                'Cebu_City_Philippines-_2.2-116714',
                'Amacuzac-Mexico-8.LIN210A1624',
                'Cuanda_Angola_01.LIN210A1719',
                'HoChiMin_Vietnam_17.LIN210A1675',
                'La_Lope-Gabon-1.LIN210A1637',
                'Cuanda_Angola_06.LIN210A1724',
                'Cuanda_Angola_05.LIN210A1723',
                'Paqueta-Brazil-23.LIN210A2149',
                'Bangkok_Thailand_09.LIN210A1687',
                'Amacuzac-Mexico-4.LIN210A1620',
                'La_Lope-Gabon-15.LIN210A1651',
                'Paqueta-Brazil-2.LIN210A2128',
                'Amacuzac-Mexico-5.LIN210A1621',
                'Cuanda_Angola_09.LIN210A1727',
                'HoChiMin_Vietnam_18.LIN210A1676',
                'Cebu_City_Philippines-_18.2-116730',
                'Cairns-Australia-16.LIN210A2166',
                'El_Dorado_Argentina_U_8.LIN210A2746',
                'Cuanda_Angola_17.LIN210A1735',
                'HoChiMin_Vietnam_16.LIN210A1674',
                'Skukusa_South_Africa_02.LIN210A1740',
                'Amacuzac-Mexico-18.LIN210A1631',
                'Paqueta-Brazil-6.LIN210A2132',
                'Bangkok_Thailand_12.LIN210A1690',
                'Skukusa_South_Africa_09.LIN210A1747',
                'Skukusa_South_Africa_04.LIN210A1742',
                'HoChiMin_Vietnam_10.LIN210A1668',
                'Maricopa_County_AZ-_21.2-116611',
                'Cebu_City_Philippines-_1.2-116713',
                'Tahiti_FrenchPolynesia_01.LIN210A1699',
                'HoChiMin_Vietnam_19.LIN210A1677',
                'Skukusa_South_Africa_06.LIN210A1744',
                'Tahiti_FrenchPolynesia_17.LIN210A1715',
                'Skukusa_South_Africa_10.LIN210A1748',
                'Tahiti_FrenchPolynesia_11.LIN210A1709',
                'Cebu_City_Philippines-_15.2-116727',
                'Tahiti_FrenchPolynesia_07.LIN210A1705',
                'El_Dorado_US_U_28',
                'Bangkok_Thailand_15.LIN210A1693',
                'La_Lope-Gabon-4.LIN210A1640']

    for dir in dirNames:
        dir = Path(ROOT) / dir
        if dir.is_dir():
            print(dir.name)
            xmlFilename = getHits(dir, 'all', 100)  #3 = get hits
            if xmlFilename is None:
                print('*** ' + str(Path(dir).name) + ': assembly failed; no results found. ***')
            else:
                drawXML(xmlFilename, str(Path(xmlFilename).parent / 'png'), False)
            
# redoHits()
# exit()

# La_Lope-Gabon-11.LIN210A1647
# Maricopa_County_Arizona-USA-20.LIN210A1615
# Maricopa_County_AZ-_19.2-116609
# El_Dorado_Argentina_U_2.LIN210A2698
# Tahiti_FrenchPolynesia_13.LIN210A1711
# Tahiti_FrenchPolynesia_19.LIN210A1717
# Cairns-Australia-8.LIN210A2158
# La_Lope_Gabon-_22.2-116710
# Cebu_City_Philippines-_2.2-116714
# Amacuzac-Mexico-8.LIN210A1624
# Cuanda_Angola_01.LIN210A1719
# HoChiMin_Vietnam_17.LIN210A1675
# La_Lope-Gabon-1.LIN210A1637
# Cuanda_Angola_06.LIN210A1724
# Cuanda_Angola_05.LIN210A1723
# Paqueta-Brazil-23.LIN210A2149
# Bangkok_Thailand_09.LIN210A1687
# Amacuzac-Mexico-4.LIN210A1620
# La_Lope-Gabon-15.LIN210A1651
# Paqueta-Brazil-2.LIN210A2128
# Amacuzac-Mexico-5.LIN210A1621
# Cuanda_Angola_09.LIN210A1727
# HoChiMin_Vietnam_18.LIN210A1676
# Cebu_City_Philippines-_18.2-116730
# Cairns-Australia-16.LIN210A2166
# El_Dorado_Argentina_U_8.LIN210A2746
# Cuanda_Angola_17.LIN210A1735
# HoChiMin_Vietnam_16.LIN210A1674
# Skukusa_South_Africa_02.LIN210A1740
# Amacuzac-Mexico-18.LIN210A1631
# Paqueta-Brazil-6.LIN210A2132
# Bangkok_Thailand_12.LIN210A1690
# Skukusa_South_Africa_09.LIN210A1747
# Skukusa_South_Africa_04.LIN210A1742
# HoChiMin_Vietnam_10.LIN210A1668
# Maricopa_County_AZ-_21.2-116611
# Cebu_City_Philippines-_1.2-116713
# Tahiti_FrenchPolynesia_01.LIN210A1699
# HoChiMin_Vietnam_19.LIN210A1677
# Skukusa_South_Africa_06.LIN210A1744
# Tahiti_FrenchPolynesia_17.LIN210A1715
# Skukusa_South_Africa_10.LIN210A1748
# Tahiti_FrenchPolynesia_11.LIN210A1709
# Cebu_City_Philippines-_15.2-116727
# Tahiti_FrenchPolynesia_07.LIN210A1705
# El_Dorado_US_U_28
# Bangkok_Thailand_15.LIN210A1693
# La_Lope-Gabon-4.LIN210A1640

# Cairns-Australia-9.LIN210A2159
# Maricopa_County_Arizona-USA-23.LIN210A1617
# Cuanda_Angola_10.LIN210A1728
# La_Lope-Gabon-24.LIN210A1658
# Paqueta-Brazil-8.LIN210A2134
# El_Dorado_US_U_31
# Cairns-Australia-12.LIN210A2162
# Tahiti_FrenchPolynesia_02.LIN210A1700
# Cuanda_Angola_08.LIN210A1726
# featuretrees.pickle.bak
# Paqueta-Brazil-11.LIN210A2137
# Amacuzac_Mexico-_2.2-116616
# Amacuzac-Mexico-16.LIN210A1629
# La_Lope-Gabon-9.LIN210A1645
# HoChiMin_Vietnam_12.LIN210A1670
# La_Lope-Gabon-12.LIN210A1648
# Paqueta-Brazil-13.LIN210A2139
# Maricopa_County_AZ-_11.2-116601
# Cairns-Australia-21.LIN210A2171
# Amacuzac-Mexico-20.LIN210A1633
# Cebu_City_Philippines-_16.2-116728
# La_Lope-Gabon-5.LIN210A1641
# Tahiti_FrenchPolynesia_08.LIN210A1706
# Tahiti_FrenchPolynesia_09.LIN210A1707
# Bangkok_Thailand_20.LIN210A1698
# Cuanda_Angola_20.LIN210A1738
# Maricopa_County_AZ-_4.2-116544
# Cebu_City_Philippines-_10.2-116722
# Cairns-Australia-19.LIN210A2169
# Amacuzac_Mexico-_21.2-116635
# HoChiMin_Vietnam_03.LIN210A1661
# Maricopa_County_AZ-_1.2-116541
# Bangkok_Thailand_17.LIN210A1695
# Cebu_City_Philippines-_14.2-116726
# Bangkok_Thailand_07.LIN210A1685
# Cebu_City_Philippines-_9.2-116721
# Bangkok_Thailand_03.LIN210A1681
# La_Lope-Gabon-14.LIN210A1650
# Amacuzac-Mexico-13.LIN210A1627
# Cuanda_Angola_03.LIN210A1721
# Skukusa_South_Africa_08.LIN210A1746
# Maricopa_County_AZ-_18.2-116608
# HoChiMin_Vietnam_01.LIN210A1659
# Cairns-Australia-15.LIN210A2165
# La_Lope-Gabon-2.LIN210A1638
# Skukusa_South_Africa_12.LIN210A1750
# Tahiti_FrenchPolynesia_18.LIN210A1716
# Cebu_City_Philippines-_11.2-116723
# Maricopa_County_AZ-_12.2-116602
# Tahiti_FrenchPolynesia_04.LIN210A1702
# El_Dorado_Argentina_U_3.LIN210A2706
# Cairns-Australia-1.LIN210A2151
# HoChiMin_Vietnam_06.LIN210A1664
# Amacuzac_Mexico-_12.2-116626
# La_Lope-Gabon-18.LIN210A1654
# Cairns-Australia-22.LIN210A2172
# La_Lope-Gabon-3.LIN210A1639
# Cairns-Australia-7.LIN210A2157
# Cebu_City_Philippines-_20.2-116732
# Paqueta-Brazil-19.LIN210A2145
# Skukusa_South_Africa_13.LIN210A1751
# Maricopa_County_AZ-_8.2-116548
# Skukusa_South_Africa_14.LIN210A1752
# Tahiti_FrenchPolynesia_14.LIN210A1712
# HoChiMin_Vietnam_09.LIN210A1667
# Amacuzac-Mexico-6.LIN210A1622
# Cebu_City_Philippines-_4.2-116716
# Maricopa_County_Arizona-USA-10.LIN210A1614
# Cuanda_Angola_14.LIN210A1732
# Amacuzac-Mexico-24.LIN210A1636
# Amacuzac-Mexico-3.LIN210A1619
# Cairns-Australia-17.LIN210A2167
# Cebu_City_Philippines-_21.2-116733
# Amacuzac-Mexico-1.LIN210A1618
# Cairns-Australia-24.LIN210A2174
# Maricopa_County_Arizona-USA-5.LIN210A1612
# Tahiti_FrenchPolynesia_05.LIN210A1703
# El_Dorado_Argentina_U_11.LIN210A2770
# Maricopa_County_Arizona-USA-7.LIN210A1613
# Paqueta-Brazil-10.LIN210A2136
# Cairns-Australia-18.LIN210A2168
# Tahiti_FrenchPolynesia_10.LIN210A1708
# Bangkok_Thailand_04.LIN210A1682
# Maricopa_County_AZ-_14.2-116604
# La_Lope-Gabon-20.LIN210A1655
# Cairns-Australia-20.LIN210A2170
# La_Lope_Gabon-_19.2-116707
# La_Lope-Gabon-10.LIN210A1646
# Maricopa_County_Arizona-USA-2.LIN210A1611
# Maricopa_County_AZ-_13.2-116603
# Tahiti_FrenchPolynesia_03.LIN210A1701
# Paqueta-Brazil-15.LIN210A2141
# Amacuzac_Mexico-_15.2-116629
# Cebu_City_Philippines-_19.2-116731
# Paqueta-Brazil-5.LIN210A2131
# Amacuzac-Mexico-17.LIN210A1630
# Cebu_City_Philippines-_8.2-116720
# Skukusa_South_Africa_15.LIN210A1753
# Maricopa_County_AZ-_16.2-116606
# Cairns-Australia-10.LIN210A2160
# HoChiMin_Vietnam_20.LIN210A1678
# HoChiMin_Vietnam_13.LIN210A1671
# Maricopa_County_Arizona-USA-22.LIN210A1616
# Cuanda_Angola_18.LIN210A1736
# Amacuzac-Mexico-11.LIN210A1626
# La_Lope_Gabon-_21.2-116709
# La_Lope-Gabon-16.LIN210A1652
# Paqueta-Brazil-9.LIN210A2135
# Bangkok_Thailand_18.LIN210A1696
# Cebu_City_Philippines-_17.2-116729
# Skukusa_South_Africa_03.LIN210A1741
# Cuanda_Angola_11.LIN210A1729
# Cairns-Australia-13.LIN210A2163
# La_Lope-Gabon-7.LIN210A1643
# Cairns-Australia-2.LIN210A2152
# Cebu_City_Philippines-_23.2-116735
# Cebu_City_Philippines-_6.2-116718
# Paqueta-Brazil-3.LIN210A2129
# Amacuzac-Mexico-7.LIN210A1623
# Maricopa_County_AZ-_24.2-116614
# HoChiMin_Vietnam_08.LIN210A1666
# Paqueta-Brazil-12.LIN210A2138
# HoChiMin_Vietnam_05.LIN210A1663
# La_Lope-Gabon-17.LIN210A1653
# Maricopa_County_AZ-_6.2-116546
# Amacuzac-Mexico-10.LIN210A1625
# Bangkok_Thailand_14.LIN210A1692
# Skukusa_South_Africa_07.LIN210A1745
# HoChiMin_Vietnam_02.LIN210A1660
# Tahiti_FrenchPolynesia_16.LIN210A1714
# Amacuzac_Mexico-_9.2-116623
# Bangkok_Thailand_02.LIN210A1680
# Bangkok_Thailand_13.LIN210A1691
# Bangkok_Thailand_11.LIN210A1689
# featuretrees.pickle
# Bangkok_Thailand_16.LIN210A1694
# HoChiMin_Vietnam_11.LIN210A1669
# Cebu_City_Philippines-_24.2-116736
# Maricopa_County_AZ-_3.2-116543
# Cebu_City_Philippines-_5.2-116717
# Cebu_City_Philippines-_3.2-116715
# Cebu_City_Philippines-_22.2-116734
# La_Lope-Gabon-13.LIN210A1649
# Cebu_City_Philippines-_13.2-116725
# Cairns-Australia-5.LIN210A2155
# Bangkok_Thailand_10.LIN210A1688
# Amacuzac-Mexico-14.LIN210A1628
# Tahiti_FrenchPolynesia_12.LIN210A1710
# HoChiMin_Vietnam_15.LIN210A1673
# Bangkok_Thailand_08.LIN210A1686
# HoChiMin_Vietnam_14.LIN210A1672
# Cairns-Australia-23.LIN210A2173
# Tahiti_FrenchPolynesia_06.LIN210A1704
# Cuanda_Angola_12.LIN210A1730
# Cuanda_Angola_16.LIN210A1734
# Paqueta-Brazil-4.LIN210A2130
# Amacuzac-Mexico-23.LIN210A1635
# Paqueta-Brazil-1.LIN210A2127
# Amacuzac-Mexico-19.LIN210A1632
# Paqueta-Brazil-18.LIN210A2144
# Tahiti_FrenchPolynesia_15.LIN210A1713
# Cuanda_Angola_19.LIN210A1737
# Bangkok_Thailand_01.LIN210A1679
# Cuanda_Angola_04.LIN210A1722
# Tahiti_FrenchPolynesia_20.LIN210A1718
# Skukusa_South_Africa_05.LIN210A1743
# HoChiMin_Vietnam_04.LIN210A1662
# La_Lope-Gabon-8.LIN210A1644
# Bangkok_Thailand_19.LIN210A1697
# Cebu_City_Philippines-_7.2-116719
# Cebu_City_Philippines-_12.2-116724
# Amacuzac-Mexico-22.LIN210A1634
# Paqueta-Brazil-24.LIN210A2150
# La_Lope-Gabon-6.LIN210A1642
# Paqueta-Brazil-7.LIN210A2133
# Bangkok_Thailand_06.LIN210A1684
# Maricopa_County_AZ-_15.2-116605
# Skukusa_South_Africa_11.LIN210A1749
# Cairns-Australia-11.LIN210A2161
# Cuanda_Angola_02.LIN210A1720
# Paqueta-Brazil-14.LIN210A2140
# Maricopa_County_AZ-_9.2-116549
# Bangkok_Thailand_05.LIN210A1683
# La_Lope-Gabon-23.LIN210A1656
# HoChiMin_Vietnam_07.LIN210A1665

def doAll():
    dirList = [dir for dir in Path(ROOT).iterdir() if 'combined' not in dir.name]

    for dir in dirList:
        dir = Path(ROOT) / dir
        if dir.is_dir():
            for file in dir.iterdir():
                if not file.is_dir() and (file.name[-26:] == '.sorted.deduped.merged.bam'):
                    doItAllViruses(str(dir.resolve()), file.name, 'all', 'virusdb')
                    break
                    
def doAllParallel2():
    dirList = [dir for dir in Path(ROOT).iterdir() if 'combined' not in dir.name]
    workers = []
    count = 0
    for dir in dirList:
        if dir.is_dir():
            for file in dir.iterdir():
                if not file.is_dir() and (file.name[-26:] == '.sorted.deduped.merged.bam'):
                    p = Process(target = doItAllViruses, args = (str(dir.resolve()), file.name, 'all', 'virusdb'))
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
        
    os.system('cp ' + ROOT + '/*/spades_all/*all_hits.xml ' + str(Path(ROOT).parent) + '/results/HITS_all')
    consolidateAll(str(Path(ROOT).parent) + '/results/HITS_all', 'all')


# def doAllParallel():
#     dirList = [dir for dir in Path(ROOT).iterdir() if 'combined' not in dir.name]
#     workers = []
#     for dir in dirList:q
#         if dir.is_dir():
#             for file in dir.iterdir():
#                 if not file.is_dir() and (file.name[-26:] == '.sorted.deduped.merged.bam'):
#                     p = Process(target = doItAllViruses, args = (dir.resolve(), file.name, 'all', 'virusdb'))
#                     workers.append(p)
#                     p.start()
#                     break
#         
#     for p in workers:
#         p.join()
        
doAll()
exit()

iTrees = None

def main():
    VIRUS_NAMES = ['cfav', 'xincheng', 'anphevirus', 'horseradish', 'aus_anopheles'] # 'liao_ning'
    VIRUS_DBS = ['cfavdb', 'xinchengdb', 'anphevirusdb', 'horseradish_db', 'aus_anophelesdb'] # 'liao_ning_db'
    VIRUS_LENGTHS = [10682, 12774, 12916, 7954, 6203]
#    LIAO_NING_SEGMENT_LENGTHS = [3740, 3055, 2404, 2062, 1870, 1750, 1208, 1147, 943, 903, 897, 760]
    combinedDirs = ['Amacuzac_Mexico_combined', 'Bangkok_Thailand_combined', 'Cairns_Australia_combined', 'Cebu_City_Philippines_combined', 'Cuanda_Angola_combined', 'El_Dorado_Argentina_combined', 'HoChiMin_Vietnam_combined', 'La_Lope-Gabon_combined', 'Maricopa_County_Arizona-USA_combined', 'Paqueta-Brazil_combined', 'Skukusa_South_Africa_combined', 'Tahiti_FrenchPolynesia_combined']
#                   'ALL_combined'

    MIN_LENGTH = 50
    
    for virusName, virusDB, virusLength in zip(VIRUS_NAMES, VIRUS_DBS, VIRUS_LENGTHS):
        print('Starting ' + virusName + '...')
#        dirList = combinedDirs
        dirList = [str(dir) for dir in Path(ROOT).iterdir() if ('combined' not in dir.name) and (dir / 'spades_cfav').exists()]
        for dirName in dirList:
            print('Analyzing ' + dirName + '...')
            dirName = str(Path(ROOT) / dirName)
#            assembleAnalyze(dirName, virusName, virusDB, MIN_LENGTH, 3, False)  # get hits, no features
            assembleAnalyze(dirName, virusName, virusDB, MIN_LENGTH, 4, True)   # get features
            
#        hitsDirName = str(Path(ROOT).parent) + '/results/HITS_' + virusName + '_combined'
        hitsDirName = str(Path(ROOT).parent) + '/results/HITS_' + virusName
        consolidate(hitsDirName, virusName, virusLength)
            
#     for index in range(len(VIRUS_NAMES)):
#         VIRUS_NAME = VIRUS_NAMES[index]
#         VIRUS_DB = VIRUS_DBS[index]
#         process_cmdline()

# MIN_LENGTH = 50
# virusName = 'cfav'
# virusDB = 'cfavdb'
# virusLength = 10682
# 
# dirName = 'Skukusa_South_Africa_04.LIN210A1742' # 'Skukusa_South_Africa_15.LIN210A1753'
# print('Analyzing ' + dirName + '...')
# dirName = str(Path(ROOT) / dirName)
# assembleAnalyze(dirName, virusName, virusDB, MIN_LENGTH, 2, False)  # get hits, no features
#             
# hitsDirName = str(Path(ROOT).parent) + '/results/HITS_' + virusName
# consolidate(hitsDirName, virusName, virusLength)
# exit()

if __name__ == '__main__':
    main()

# import matplotlib.pyplot as pyplot
#     
# def plotLengths():
#     VIRUS_NAMES = ['cfav', 'xincheng', 'anphevirus', 'horseradish', 'aus_anopheles'] # 'liao_ning'
# 
#     lengths = []
#     evalues = []
#     for dir in Path(ROOT).iterdir():
#         for virus in VIRUS_NAMES:
#             fn = dir / ('spades_cfav/blast_scaffolds_' + virus + '.csv')
#             if fn.exists():
#                 hits = readCSV(str(fn), 0)
#                 for contig in hits:
#                     for hit in hits[contig]:
#                         if (len(hit[2]) < 100): # and (float(hit[6]) > 1e-10):
#                             lengths.append(len(hit[2]))
#                             evalues.append(float(hit[6]))
#     pyplot.scatter(lengths, evalues, s = 3)
#     pyplot.xlabel('Length')
#     pyplot.ylabel('e-value')
#     pyplot.show()
    
#     pyplot.hist(evalues)
#     pyplot.show()
    
# plotLengths()

# dirName = '/home/havill/data/aegypti/analyzed/Bangkok_Thailand_combined'
# virusName = 'cfav'
# virusDB = 'cfavdb'
# # xmlResultFile = dir + '/spades_' + virusName + '/' + dir + '_' + virusName + '_hits.xml'
# # displayXML(xmlResultFile[:-4] + '_features.xml', xmlResultFile[:-4] + '_features_short.txt', False)
# assembleAnalyze(dirName, virusName, virusDB, 3, False)
# 
# virusName = 'cfav'
# virusLength = 10682
# hitsDirName = '/home/havill/data/aegypti/results/HITS_' + virusName + '_combined'
# consolidate(hitsDirName, virusName, virusLength)
