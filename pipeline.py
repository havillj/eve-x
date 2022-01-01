#!/usr/bin/env python3

from util import *
from config import *
from virushit import VirusHit
from hosthit import HostHit
from diagramhit import DiagramHits
from kmeans import findClusters
from features import *

###############################################################################

def getUnmappedReads(filenameBAM):
    """Get unmapped reads and their mates from a BAM file and write them to 
       a new BAM file.
    
       Parameter:
           filenameBAM: absolute path of input BAM file as a Path object
           
       Return value: absolute path of the output BAM file as a Path object
    """
    
    shortName = filenameBAM.name.split('.')[0]
    specimenResultsPath = Path(SPECIMEN_RESULTS_DIR) / shortName
    if not specimenResultsPath.exists():
        os.system('mkdir ' + str(specimenResultsPath))
    
    newFilenameBAM = specimenResultsPath / (shortName + '_unmapped_with_mates.bam')  # absolute path

    writelog('Getting unmapped reads ...', True)

    if not newFilenameBAM.exists():
        bamfile = pysam.AlignmentFile(str(filenameBAM), 'rb', threads=8)
        newBAM = pysam.AlignmentFile(str(newFilenameBAM), 'wb', template=bamfile, threads=8)

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

        writelog('   Total reads = ' + str(count_total), True)
        writelog('   Unmapped reads = ' + str(count) + ' ({:.2f}%)'.format(100*count/count_total), True)
    else:
        writelog('   Skipping - ' + str(newFilenameBAM) + ' exists.', True)
        
    return newFilenameBAM  # absolute path
    
###############################################################################

def writeReadsFASTQ(filenameBAM):
    """Write paired-end reads from BAM file to 3 FASTQ files.
    
       Parameters:
           viralFilenameBAM: absolute path of BAM file containing potential 
                             viral reads as a Path object
           
       Return value: absolute path of SPAdes working directory as a Path object
    """
    
    writelog('Writing paired-end reads to FASTQ files ...', True)
    
    spadesPath = filenameBAM.parent / SPADES_DIR
    
    if not spadesPath.exists():
        os.system('mkdir ' + str(spadesPath))
        
    if (spadesPath / 'viral1.fastq').exists():
        writelog('   Skipping - FASTQ files exist.', True)
    else:
        os.system(SAMTOOLS_EXEC + ' fastq -1 ' + str(spadesPath / 'viral1.fastq') 
                  + ' -2 ' + str(spadesPath / 'viral2.fastq') 
                  + ' -s ' + str(spadesPath / 'viral_single.fastq') + ' ' 
                  + str(filenameBAM))
              
    return spadesPath
    
###############################################################################

def assembleReads(dirName):
    """Assemble reads in dirName/spades into scaffolds.
    
       Parameter: 
           dirName: absolute path of parent of directory containing fastq files 
            
       Return value: Boolean indicating whether assembly was successful
    """

    writelog('Assembling reads with SPAdes ...', True)
    
    spadesPath = Path(dirName)
    if (spadesPath / 'scaffolds.fasta').exists():
        writelog('   Skipping - scaffolds.fasta exists.', True)
        return True
        
    viral1Name = str(spadesPath / 'viral1.fastq')
    viral2Name = str(spadesPath / 'viral2.fastq')
    viralSName = str(spadesPath / 'viral_single.fastq')
    
    if os.path.getsize(viralSName) > 0:
        result = os.system(SPADES_EXEC + ' --pe1-1 ' + viral1Name 
                                       + ' --pe1-2 ' + viral2Name 
                                       + ' --s1 ' + viralSName 
                                       + ' --careful -o ' + str(spadesPath))
    else:
        os.remove(viralSName)
        result = os.system(SPADES_EXEC + ' --pe1-1 ' + viral1Name 
                                       + ' --pe1-2 ' + viral2Name 
                                       + ' --careful -o ' + str(spadesPath))
                
    return (spadesPath / 'scaffolds.fasta').exists()  # assembly succeeded
    
###############################################################################
    
def blastScaffolds(dirName, force = False):
    """BLAST scaffolds against viral database and write results to a CSV file.
    
       Parameters: 
           dirName: absolute path of parent of directory containing scaffold
                    FASTA files
           force:   whether to overwrite existing blast_scaffolds.csv
            
       Return value: None
    """
    
    spadesPath = Path(dirName) / SPADES_DIR
    scaffoldsName = str(spadesPath / 'scaffolds.fasta')
    outVirusCSV = spadesPath / 'blast_scaffolds.csv'
    
    writelog('BLASTing scaffolds against viral database ...', True)
    
    if not force and outVirusCSV.exists():
        writelog('   Skipping - ' + str(outVirusCSV) + ' exists.', True)
    else:
        os.system(BLAST_EXEC + ' -query ' + scaffoldsName 
                      + ' -db ' + VIRUS_DB 
                      + ' -num_threads 16'
                      + ' -task blastn' 
                      + ' -evalue ' + str(config['EVALUE_VIRUS'])
                      + ' -max_target_seqs 5000'
                      + ' -gapopen 2'
                      + ' -outfmt "10 qseqid qstart qend qseq sstart send sseq evalue bitscore sseqid stitle pident"'
                      + ' -out ' + str(outVirusCSV))

#                 ' -gapopen 5 -gapextend 2 -penalty -3 -reward 2 -word_size 11'

###############################################################################
                  
def readCSV_viral(csvFilename):
    """Read a BLAST CSV file containing viral hits.
     
       Parameter:
           csvFilename: absolute path of CSV file containing BLAST results
           
       Return value: dictionary with contig keys and values of lists of VirusHits
    """
    
    try:
        csv = open(csvFilename, 'r')
    except FileNotFoundError:
        return None
        
    hits = {}           # hits in this file
    accsFound = {}      # accession ids found for each species
    for line in csv:
        row = line.rstrip().split(',')
        row[10] = row[10].strip('|')   # viral genome title                  
        hit = VirusHit(row)
        
        contig = hit.contig
        if contig not in hits:
            hits[contig] = []
            accsFound[contig] = {}
        hits[contig].append(hit)
        if hit.stitle not in accsFound[contig]:
             accsFound[contig][hit.stitle] = []
        accsFound[contig][hit.stitle].append((hit.sseqid, hit.slength))
    csv.close()
    
    # Find preferred accession id(s) for each species.
    
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
        
    # Filter hits by preferred accession id(s).
    
    filteredHits = {}
    for contig in hits:
        filteredHits[contig] = []
        for hit in hits[contig]:
            if hit.sseqid in accsToUse[contig][hit.stitle]:
                filteredHits[contig].append(hit)
            else:
                writelog(contig + ': discarded hit "' + hit.stitle + '" with duplicate species')

    return filteredHits
    
def readCSV_host(csvFilename):
    """Read a BLAST CSV file containing host hists.
     
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
        row = line.rstrip().split(',')
        contig = row[0]  # qseqid = contig name

        if contig not in hits:
            hits[contig] = []
        hits[contig].append(HostHit(row))

    return hits
    
def getHits(dirName):
    """Combine viral and vector BLAST results to find putative insertions.
    
       Parameters:
           dirName: absolute path of parent of directory containing SPAdes results 
           
       Return value: absolute path of XML file containing results
    """
    
    writelog('Combining BLAST results to locate putative EVEs...', True)
    
    spadesPath = Path(dirName) / SPADES_DIR              # specimen assembly dir
    
    resultsPath = Path(dirName) / SPECIMEN_RESULTS_DIR   # specimen results dir
    if not resultsPath.exists():
        os.system('mkdir ' + str(resultsPath))
        
    scaffoldsPath = resultsPath / SCAFFOLDS_DIR
    if not scaffoldsPath.exists():
        os.system('mkdir ' + str(scaffoldsPath))
        
    os.system('cp ' + str(spadesPath / 'scaffolds.fasta') + ' ' + str(scaffoldsPath))
    os.system('cp ' + str(spadesPath / 'blast_scaffolds.csv') + ' ' + str(scaffoldsPath))
    scaffoldsFilename = str(scaffoldsPath / 'scaffolds.fasta')
    csvVirusFilename = str(scaffoldsPath / 'blast_scaffolds.csv')
    
    xmlPath = resultsPath / XML_DIR
    if not xmlPath.exists():
        os.system('mkdir ' + str(xmlPath))
    xmlFilename = str(xmlPath / (Path(dirName).name + '_hits.xml'))
   
    virusHits = readCSV_viral(csvVirusFilename)
    if virusHits is None:
        return None
        
    # virusHits[contig] = [hit0, hit1, ...]
    # where each hit = [qstart, qend, qseq, sstart, send, sseq, evalue, bitscore, sseqid, stitle, pident]

    ##### Combine hits in each contig into mosaic hits #####
    
    hitIntervals = {}
    maxLenPIdent = {}
    for contig in virusHits:
        hitIntervals[contig] = []  # list of (left, right) intervals delineating hits
        leftIndex = 0
        rightIndex = 0
        
        # sort hits by virus seq id, then position in contig
        virusHits[contig].sort(key = lambda hit: (hit.sseqid, hit.qstart))
        
        virusHit = virusHits[contig][0]   # first hit in this contig
        sseqid = virusHit.sseqid           # virus accession number
        hitLength = virusHit.qlength
        hitSeq = virusHit.qseq
        pident = virusHit.pident          # max pident among all hits for virus
        maxHitLength = hitLength
        maxHitLengthIndex = 0
        maxPIdent = pident
        maxPIdentIndex = 0
        
        # iterate over remaining hits in the contig
        # (one extra to close last interval)
        for index in range(1, len(virusHits[contig]) + 1): 
            if index < len(virusHits[contig]):
                virusHit = virusHits[contig][index]
                nextSeqID = virusHit.sseqid
            else:                              # this is the extra iteration
                virusHit = None
                nextSeqID = ''
                
            if (nextSeqID == sseqid):             # continue same mosaic hit
                hitLength += virusHit.qlength
                hitSeq += virusHit.qseq
                pident = max(pident, virusHit.pident)
            else:                                 # finish mosaic hit and start a new one
                if hitLength >= config['MIN_VIRUS_HIT_LENGTH']:
                    hitComplexity = complexity(hitSeq, config['COMPLEXITY_K'])  # COMPLEXITY TEST
                    if hitComplexity < config['COMPLEXITY_CUTOFF']:
                        writelog(contig + ': discarded hit "' + sseqid \
                                        + '" with complexity ' + str(hitComplexity) \
                                        + ' < ' + str(config['COMPLEXITY_CUTOFF']))
                    else:
                        hitIntervals[contig].append((leftIndex, rightIndex))
                        if hitLength > maxHitLength:
                            maxHitLength, maxHitLengthIndex = hitLength, leftIndex
                        if pident > maxPIdent:
                            maxPIdent, maxPIdentIndex = pident, leftIndex
                else:
                    writelog(contig + ': discarded hit "' + sseqid \
                                    + '" with length ' + str(hitLength) \
                                    + ' < ' + str(config['MIN_VIRUS_HIT_LENGTH']))

                # start new mosaic hit
                if index < len(virusHits[contig]):
                    leftIndex = index
                    hitLength = virusHit.qlength
                    hitSeq = virusHit.qseq
                    pident = virusHit.pident
            rightIndex = index
            sseqid = nextSeqID

        maxLenPIdent[contig] = (maxHitLength, maxHitLengthIndex, maxPIdent, maxPIdentIndex)
        
##### Search contigs for host hits to locate putative insertions.
    
    # create a new XML tree
    root = ET.Element('root')
    tree = ET.ElementTree(root)

    all_viral_hits = []
    
    scaffoldRecords = SeqIO.index(scaffoldsFilename, 'fasta')
    
    for contig in virusHits:
        # search for host hits using megablast
        SeqIO.write(scaffoldRecords[contig], str(spadesPath / 'contig_temp.fasta'), 'fasta')
        os.system(BLAST_EXEC + ' -query ' + str(spadesPath / 'contig_temp.fasta')
                      + ' -db ' + HOST_DB 
                      + ' -num_threads 8'
                      + ' -evalue ' + str(config['EVALUE_HOST'])
                      + ' -max_target_seqs 100'
                      + ' -outfmt "10 qseqid qstart qend qseq sseqid sstart send sseq evalue bitscore"'
                      + ' -out ' + str(spadesPath / 'blast_contig_temp.csv'))
        hostHits = readCSV_host(str(spadesPath / 'blast_contig_temp.csv'))

        all_viral_hits.extend([(v.sseqid, v.stitle) for v in virusHits[contig]])  
        
        # get min evalue, max bitscore, max length, max pident for this contig
        
        minEvalueIndex = 0
        for index in range(1, len(virusHits[contig])):
            if virusHits[contig][index].evalue < virusHits[contig][minEvalueIndex].evalue:
                minEvalueIndex = index
                
        maxBitscoreIndex = 0
        for index in range(1, len(virusHits[contig])):
            if virusHits[contig][index].bitscore > virusHits[contig][maxBitscoreIndex].bitscore:
                maxBitscoreIndex = index
                
        maxHitLength, maxHitLengthIndex, maxPIdent, maxPIdentIndex = maxLenPIdent[contig]
        
        # Create a separate xml node for each virus found in the contig
        
        thisContigCount = 0
        for leftIndex, rightIndex in hitIntervals[contig]:
            thisContigCount += 1
            virus_qstart = virusHits[contig][leftIndex].qstart
            virus_qend = virusHits[contig][rightIndex].qend
            
            # tree of intervals (qstart, qend, (sstart, send))
            ivs = [Interval(virusHits[contig][i].qstart, 
                            virusHits[contig][i].qend, 
                            (virusHits[contig][i].sstart, virusHits[contig][i].send)) 
                   for i in range(leftIndex, rightIndex + 1)]
            virusIntervals = IntervalTree(ivs)
            virusIntervalsLength = sum([iv.length() for iv in virusIntervals])
            
            bestHit = False
            if (maxPIdent >= config['MIN_PIDENT']):
                if maxPIdentIndex == leftIndex:
                    bestHit = True
            else:
                if maxHitLengthIndex == leftIndex:
                    bestHit = True

            if len(hitIntervals[contig]) > 1:
                contigName = contig + '__' + str(thisContigCount)
            else:
                contigName = contig
                
            contigTree = ET.SubElement(root, 'contig', {'name': contigName, 'besthit': str(bestHit)})
            for index in range(leftIndex, rightIndex + 1):                
                virusHit = virusHits[contig][index]
                vElement = ET.SubElement(contigTree, 'virushit', 
                                                     {'seqid': virusHit.sseqid, 
                                                      'stitle': virusHit.stitle, 
                                                      'minevalue': str(index == minEvalueIndex), 
                                                      'maxbitscore': str(index == maxBitscoreIndex)})
                ET.SubElement(vElement, 'qstart').text = str(virusHit.qstart)
                ET.SubElement(vElement, 'qend').text = str(virusHit.qend)
                ET.SubElement(vElement, 'qseq').text = virusHit.qseq
                ET.SubElement(vElement, 'sstart').text = str(virusHit.sstart)
                ET.SubElement(vElement, 'send').text = str(virusHit.send)
                ET.SubElement(vElement, 'sseq').text = virusHit.sseq
                ET.SubElement(vElement, 'evalue').text = str(virusHit.evalue)
                ET.SubElement(vElement, 'bitscore').text = str(virusHit.bitscore)
                ET.SubElement(vElement, 'pident').text = str(virusHit.pident)
                
            # partition host hits into before and after viral hit
            
            before = []    # indices of host hits before, 
            after = []     # after, and
            overlap = []   # overlapping viral hit in contig string
            virusIntervalsNotInRef = copy.deepcopy(virusIntervals)
            doesOverlap = False
            if contig in hostHits:  # if there were host hits, examine them
                p = re.compile(r'.*length_(\d+)_.*')
                contigLength = int(p.match(contig).group(1))
                
                for index in range(len(hostHits[contig])):
                    hostHit = hostHits[contig][index]
                    host_qseq = hostHit.qseq
                    host_qstart = hostHit.qstart
                    host_qend = hostHit.qend
                    if ((virus_qstart - host_qend > config['MAX_HOST_VIRUS_DISTANCE']) or (host_qstart - virus_qend > config['MAX_HOST_VIRUS_DISTANCE'])) and (len(host_qseq) < config['MIN_FLANKING_HIT_LENGTH']):
                        writelog(contig + ': discarded distant host hit with length ' + str(len(host_qseq)) + ' < ' + str(config['MIN_FLANKING_HIT_LENGTH']))
                        continue
                    virusIntervalsNotInRef.chop(host_qstart, host_qend)  # detect overlap with virus hit
                    if percentOverlap((virus_qstart, virus_qend), (host_qstart, host_qend)) >= config['HOST_OVERLAP_FRACTION']:
                        overlap.append(index)
                        if host_qstart < virus_qstart and host_qend >= virus_qend:
                            doesOverlap = True
                    elif host_qstart < virus_qstart and host_qend < virus_qend:
                        before.append(index)
                    elif host_qstart > virus_qstart and host_qend > virus_qend:
                        after.append(index)
                    else:  # as >= vs and ae <= ve or as <= vs and ae >= ve
                        overlap.append(index)  # host hit completely contained in virus hit

                doesOverlap = doesOverlap or (sum([iv.length() for iv in virusIntervalsNotInRef]) <= (1 - config['OVERLAP_FRACTION']) * virusIntervalsLength)
                                        
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
                    hostHit = hostHits[contig][i]
                    hostElement = ET.SubElement(contigTree, name, {'id': str(i), 'seqid': hostHit.sseqid})
                    ET.SubElement(hostElement, 'qstart').text = str(hostHit.qstart)
                    ET.SubElement(hostElement, 'qend').text = str(hostHit.qend)
#                    ET.SubElement(hostElement, 'qseq').text = hostHit.qseq
                    ET.SubElement(hostElement, 'sstart').text = str(hostHit.sstart)
                    ET.SubElement(hostElement, 'send').text = str(hostHit.send)
#                    ET.SubElement(hostElement, 'sseq').text = hostHit.sseq
                    ET.SubElement(hostElement, 'evalue').text = str(hostHit.evalue)
                    ET.SubElement(hostElement, 'bitscore').text = str(hostHit.bitscore)
                    
            matchingFlanks = ET.SubElement(contigTree, 'flanks')
            # find pairs of before/after host hits that are on the same strand and within config['MAX_FLANK_DISTANCE'] of each other
            for b in before:
                beforeHit = hostHits[contig][b]
                for a in after:
                    afterHit = hostHits[contig][a]
                    if beforeHit.sseqid == afterHit.sseqid:  # same chromosome
                        beforeStart = beforeHit.sstart
                        beforeEnd = beforeHit.send
                        afterStart = afterHit.sstart
                        afterEnd = afterHit.send
                        flankDistance = max(beforeStart, beforeEnd, afterStart, afterEnd) - min(beforeStart, beforeEnd, afterStart, afterEnd)
                        if flankDistance <= config['MAX_FLANK_DISTANCE']:
                            if (beforeStart < beforeEnd < afterStart + config['ALLOWED_OVERLAP'] < afterEnd + config['ALLOWED_OVERLAP']) or \
                               (beforeStart > beforeEnd > afterStart - config['ALLOWED_OVERLAP'] > afterEnd - config['ALLOWED_OVERLAP']):
                                ET.SubElement(matchingFlanks, 'match', {'leftid': str(b), 'rightid': str(a)})
                            else:
                                ET.SubElement(matchingFlanks, 'inversion', {'leftid': str(b), 'rightid': str(a)})
     
    tree.write(xmlFilename, xml_declaration=True, pretty_print=True)
        
    return xmlFilename

###############################################################################
    
def drawContigs(fileName):
    """Draw diagrams of contigs from XML file.
    
       Parameters:
           fileName: absolute path of XML results file
           
       Return value: None
    """
    
    diagramsPath = Path(fileName).parent.parent / DIAGRAMS_DIR
    
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
        if config['BEST_HITS_ONLY'] and (contig.attrib['besthit'] != 'True'):
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
            virus_stitle = v.attrib['stitle']
            virus_seqid = cleanACC(v.attrib['seqid'])
            
            if family is None:               # set family with first viral hit; overwrite later if find best hit
                if virus_seqid not in fams:
                    fam = getFamily(virus_seqid)
                    fams[virus_seqid] = fam
                    addFamily(virus_seqid, fam)
                family = fams[virus_seqid]
                
            stitle_seqid = virus_stitle + ' (' + str(virus_seqid) + ')'
            if not config['BEST_HITS_ONLY'] and (contig.attrib['besthit'] == 'True'):                
                if virus_seqid not in fams:
                    fam = getFamily(virus_seqid)
                    fams[virus_seqid] = fam
                    addFamily(virus_seqid, fam)
                family = fams[virus_seqid]
                
            if qend > qstart:
                strand = 1
            else:
                strand = -1
            gf = GraphicFeature(start=qstart, end=qend, strand=strand, thickness=thickness, linewidth=0,
                                color='#ff0000', fontdict = {'size': labelFontSize},
                                label = stitle_seqid.lstrip('|') + ' ' + str(sstart) + '-' + str(send) + ' (' + str(length) + ' bp; ' + str(pident) + '%)')  #; ' + str(evalue) + ')')
            features.append(gf)
            gfCopy = copy.deepcopy(gf)
            if not config['BEST_HITS_ONLY'] and (contig.attrib['besthit'] == 'True'):
                gfCopy.label = '**' + stitle_seqid.lstrip('|') + '** ' + str(sstart) + '-' + str(send) + ' (' + str(length) + ' bp; ' + str(pident) + '%)'
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
            for match in matches[:config['FLANK_DRAW_LIMIT']]:
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
                            if seqid in CHR_NAMES:
                                seqid = CHR_NAMES[seqid]
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
            for match in inversions[:max(0, config['FLANK_DRAW_LIMIT'] - len(matches))]:
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
                            if seqid in CHR_NAMES:
                                seqid = CHR_NAMES[seqid]
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
                for q in posSort[:max(0, config['FLANK_DRAW_LIMIT'] - len(matches) - len(inversions))]:
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
                    if seqid in CHR_NAMES:
                        seqid = CHR_NAMES[seqid]
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
        specimen = Path(fileName).parent.parent.parent.name
        if contig.attrib['besthit'] == 'True':
            pp = PdfPages(str(Path(familyDir) / (specimen + '_' + contig.attrib['name'].replace('__', '-') + '_' + virus_stitle.replace(' ', '_') + '_' + virus_seqid + '_BEST-HIT.pdf')))
            pp.savefig(fig)
            pp.close()
        else:
            pp = PdfPages(str(Path(familyDir) / (specimen + '_' + contig.attrib['name'].replace('__', '-') + '_' + virus_stitle.replace(' ', '_') + '_' + virus_seqid + '.pdf')))
            pp.savefig(fig)
            pp.close()
        pyplot.close(fig)
        
    if not config['BEST_HITS_ONLY']:
        for c in contigFeaturesV:
            record = GraphicRecord(sequence_length = contigLengths[c], features = contigFeaturesV[c])
        
            ax, _ = record.plot(max_label_length = 80, figure_width = 10)
            familyDir = str(diagramsPath / families[c])
            pp = PdfPages(str(Path(familyDir) / (specimen + '_' + c + '_' + virus_stitle.replace(' ', '_') + '_' + virus_seqid + '_all-viral-hits.pdf')))
            pp.savefig(ax.figure)
            pp.close()
            pyplot.close(ax.figure)
                    
###############################################################################
        
def writeSummaryTable(fileName, viralHits, viralHits0, viralHits1, viralHits2, allFamilies, allViruses):
    """Create summary TSV of all viral hits.
    
       Parameters:
           fileName: absolute path of output file
           viralHits:
           viralHits0:
           viralHits1:
           viralHits2:
           allFamilies:
           allViruses:
           
       Return value: None
    """
    
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
                    if config['OMIT_EVES_IN_REFERENCE'] and (specimen in viralHits0) and (hit in viralHits0[specimen]):  # overlapping vector hit
                        continue
                    foundInSpecimen = True
                    break
            if foundInSpecimen:
                counts[(fam, (stitle, seqid))] += 1
                
    newFamVir = []
    for fam, (stitle, seqid) in famVir:
        if counts[(fam, (stitle, seqid))] >= config['MIN_HITS_TO_SHOW_VIRUS']:
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
                    if config['OMIT_EVES_IN_REFERENCE'] and (specimen in viralHits0) and (hit in viralHits0[specimen]):  # overlapping vector hit
                        continue
                    hitDict = hit[2]
                    qpos = list(hitDict.keys())
                    qpos.sort()                  # within hit, sort by qstart
                    s = []
                    for p in qpos:
                        length += abs(hitDict[p][1] - hitDict[p][0]) + 1
                        s.extend([hitDict[p][0], hitDict[p][1]])
#                    if s[0] > s[-1]:
#                        s = s[::-1]
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
    """
    
    writelog('Consolidating results in ' + RESULTS_DIR + '...', True)
        
    viralHits1 = {}
    viralHits2 = {}
    viralHits0 = {}   # with host overlap
    viralHits = {}    # all together
    allHits = []    
    viralSeqs = {}
    viralSeqStrings = {}
    viralSeqStringsWithFlanks = {}
    viralSeqPositions = {}
    refSeqs = {}
    pIdents = {}
    
    xmlFiles = []
    specimensPath = Path(SPECIMENS_RESULTS_DIR)
    for specimenDir in specimensPath.iterdir():
        if specimenDir.name[0] != '.':
            for file in (specimenDir / XML_DIR).iterdir():
                if (file.name[0] != '.') and ('_hits.xml' in file.name):
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
        seqPath = file.parent.parent / SEQUENCES_DIR
        if not seqPath.exists():
            os.system('mkdir ' + str(seqPath))
        virus_fasta = open(str(seqPath / (specimen + '_hits_aligned.fasta')), 'w')
        virus_fasta_percontig = open(str(seqPath / (specimen + '_hits_unaligned.fasta')), 'w')

        contigs = SeqIO.index(str(file.parent.parent / SCAFFOLDS_DIR / 'scaffolds.fasta'), 'fasta')

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
            
            if config['BEST_HITS_ONLY'] and (contig.attrib['besthit'] != 'True'):  # and (stitle not in PREFERRED_ACCS):  # last condition may result in >1 hit per original contig
                continue
                
            if config['OMIT_EVES_IN_REFERENCE'] and (contig.attrib['inreference'] == 'True'):
                continue
            
            if (seqid, stitle) not in viralSeqs:
                viralSeqs[(seqid, stitle)] = {}
                viralSeqStrings[(seqid, stitle)] = {}
                viralSeqStringsWithFlanks[(seqid, stitle)] = {}
                viralSeqPositions[(seqid, stitle)] = {}
                refSeqs[(seqid, stitle)] = {}
                pIdents[(seqid, stitle)] = {}
#                insertSites[(seqid, stitle)] = {}
            if specimen not in viralSeqs[(seqid, stitle)]:
                viralSeqs[(seqid, stitle)][specimen] = {}
                viralSeqStrings[(seqid, stitle)][specimen] = {}
                viralSeqStringsWithFlanks[(seqid, stitle)][specimen] = {}
                viralSeqPositions[(seqid, stitle)][specimen] = {}
                refSeqs[(seqid, stitle)][specimen] = {}
                pIdents[(seqid, stitle)][specimen] = {}
#                insertSites[(seqid, stitle)][specimen] = {}
            viralSeqs[(seqid, stitle)][specimen][contigName] = {}
            refSeqs[(seqid, stitle)][specimen][contigName] = {}
            pIdents[(seqid, stitle)][specimen][contigName] = {}
#            insertSites[(seqid, stitle)][specimen][contigName] = {}
            for v in contig.findall('virushit'):
                hit[(int(v.find('qstart').text), int(v.find('qend').text))] = (int(v.find('sstart').text), int(v.find('send').text))
                viralSeqs[(seqid, stitle)][specimen][contigName][(int(v.find('sstart').text), int(v.find('send').text))] = v.find('qseq').text
                refSeqs[(seqid, stitle)][specimen][contigName][(int(v.find('sstart').text), int(v.find('send').text))] = v.find('sseq').text
                pIdents[(seqid, stitle)][specimen][contigName][(int(v.find('sstart').text), int(v.find('send').text))] = float(v.find('pident').text)
                virus_fasta.write('>' + seqid + '_' + stitle.replace(' ', '-') + '_' + v.find('sstart').text + '-' + v.find('send').text + '_(host)\n')
                virus_fasta.write(v.find('qseq').text + '\n')
                virus_fasta.write('>' + seqid + '_' + stitle.replace(' ', '-') + '_' + v.find('sstart').text + '-' + v.find('send').text + '_(reference)\n')
                virus_fasta.write(v.find('sseq').text + '\n')
        
            qpos = list(hit.keys())  # list of (qstart, qend) tuples
            qpos.sort()
            viralSeqStrings[(seqid, stitle)][specimen][contigName] = str(contigs[contigName.split('__')[0]].seq)[qpos[0][0]-1:qpos[-1][1]] # getContig(str(file.parent.parent.parent), contigName.split('__')[0], qpos[0][0], qpos[-1][1])
            viralSeqStringsWithFlanks[(seqid, stitle)][specimen][contigName] = str(contigs[contigName.split('__')[0]].seq)
            sstart1, send1 = hit[qpos[0]]
            if (len(hit) == 1) and (sstart1 > send1):   # if not a mosaic hit, take rev comp as necessary
                viralSeqStrings[(seqid, stitle)][specimen][contigName] = reverseComplement(viralSeqStrings[(seqid, stitle)][specimen][contigName])
                viralSeqStringsWithFlanks[(seqid, stitle)][specimen][contigName] = reverseComplement(viralSeqStringsWithFlanks[(seqid, stitle)][specimen][contigName])
                viralSeqPositions[(seqid, stitle)][specimen][contigName] = str(send1) + '-' + str(sstart1)
            else:
                sposString = ''
                for (qstart, qend) in qpos:
                    sposString += str(hit[(qstart, qend)][0]) + '-' + str(hit[(qstart, qend)][1]) + ','
                viralSeqPositions[(seqid, stitle)][specimen][contigName] = sposString[:-1]
            virus_fasta_percontig.write('>' + seqid + '_' + stitle.replace(' ', '-') + '_' + viralSeqPositions[(seqid, stitle)][specimen][contigName] + '\n')
            virus_fasta_percontig.write(viralSeqStrings[(seqid, stitle)][specimen][contigName] + '\n')

#             recordedIDs = []
#             insertSites[(seqid, stitle)][specimen][contigName]['match'] = []
#             flanks = contig.find('flanks')
            vectorHitsLeft = contig.findall('vectorhitleft')
            vectorHitsRight = contig.findall('vectorhitright')
#             for match in flanks.findall('match'):
#                 recordedIDs.append(match.attrib['leftid'])
#                 recordedIDs.append(match.attrib['rightid'])
#                 for v in vectorHitsLeft:
#                     if v.attrib['id'] == match.attrib['leftid']:
#                         matchLeft = int(v.find('send').text)
#                         break
#                 for v in vectorHitsRight:
#                     if v.attrib['id'] == match.attrib['rightid']:
#                         insertSites[(seqid, stitle)][specimen][contigName]['match'].append((v.attrib['seqid'], (matchLeft, int(v.find('sstart').text))))
#                         break
#                         
#             insertSites[(seqid, stitle)][specimen][contigName]['inversion'] = []
#             for match in flanks.findall('inversion'):
#                 recordedIDs.append(match.attrib['leftid'])
#                 recordedIDs.append(match.attrib['rightid'])
#                 for v in vectorHitsLeft:
#                     if v.attrib['id'] == match.attrib['leftid']:
#                         matchLeft = int(v.find('send').text)
#                         break
#                 for v in vectorHitsRight:
#                     if v.attrib['id'] == match.attrib['rightid']:
#                         insertSites[(seqid, stitle)][specimen][contigName]['inversion'].append((v.attrib['seqid'], (matchLeft, int(v.find('sstart').text))))
#                         break
#         
#             insertSites[(seqid, stitle)][specimen][contigName]['left'] = []
#             for v in vectorHitsLeft:
#                 if v.attrib['id'] not in recordedIDs:
#                     host_qstart = int(v.find('qstart').text)
#                     host_qend = int(v.find('qend').text)
#                     virus_qstart = qpos[0][0]
#                     virus_qend = qpos[-1][1]
#                     if (-ALLOWED_OVERLAP <= virus_qstart - host_qend <= MAX_FLANKING_DISTANCE) or (-ALLOWED_OVERLAP <= host_qstart - virus_qend <= MAX_FLANKING_DISTANCE): # or (len(host_qseq) >= MIN_FLANKING_HIT_LENGTH):
#                         insertSites[(seqid, stitle)][specimen][contigName]['left'].append((v.attrib['seqid'], int(v.find('send').text)))
#                 
#             insertSites[(seqid, stitle)][specimen][contigName]['right'] = []
#             for v in vectorHitsRight:
#                 if v.attrib['id'] not in recordedIDs:
#                     host_qstart = int(v.find('qstart').text)
#                     host_qend = int(v.find('qend').text)
#                     virus_qstart = qpos[0][0]
#                     virus_qend = qpos[-1][1]
#                     if (-ALLOWED_OVERLAP <= virus_qstart - host_qend <= MAX_FLANKING_DISTANCE) or (-ALLOWED_OVERLAP <= host_qstart - virus_qend <= MAX_FLANKING_DISTANCE): # or (len(host_qseq) >= MIN_FLANKING_HIT_LENGTH):
#                         insertSites[(seqid, stitle)][specimen][contigName]['right'].append((v.attrib['seqid'], int(v.find('sstart').text)))
        
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
        
    if not Path(VIRUS_RESULTS_DIR).exists():
        os.system('mkdir ' + VIRUS_RESULTS_DIR)
    if not (Path(VIRUS_RESULTS_DIR) / SEQUENCES_DIR).exists():
        os.system('mkdir ' + str((Path(VIRUS_RESULTS_DIR) / SEQUENCES_DIR)))
        
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
        refRecord = getSeqRecord(referenceFilename, seqid)
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

        FAMILY_DIR = Path(VIRUS_RESULTS_DIR) / SEQUENCES_DIR / fams[seqid]
        if not Path(FAMILY_DIR).exists():
            os.system('mkdir ' + str(FAMILY_DIR))
        
        seqOutFileName = str(Path(FAMILY_DIR) / (seqid + '_per_specimen_aligned.fasta'))
        seqOutPerContigFileName = str(Path(FAMILY_DIR) / (seqid + '_per_contig_aligned.fasta'))
        seqOutPerContigUnalignedFileName = str(Path(FAMILY_DIR) / (seqid + '_per_contig_unaligned.fasta'))
        seqOutPerContigWithFlanksFileName = str(Path(FAMILY_DIR) / (seqid + '_per_contig_with_flanks.fasta'))
        seqOut = open(seqOutFileName, 'w')
        seqOutPerContig = open(seqOutPerContigFileName, 'w')
        seqOutPerContigUnaligned = open(seqOutPerContigUnalignedFileName, 'w')
        seqOutPerContigWithFlanks = open(seqOutPerContigWithFlanksFileName, 'w')
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
#            posSpecimen.sort()
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
                    
#                positions[(seqid, stitle)][specimen][contigName].sort()
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
                
                seqOutPerContigUnaligned.write('>' + specimen + '_|_' + contigName + '_|_' + seqid + '_' + viralSeqPositions[(seqid, stitle)][specimen][contigName] + '_pident_' + '{:5.3f}'.format(pident) + '\n')
                seqOutPerContigUnaligned.write(viralSeqStrings[(seqid, stitle)][specimen][contigName] + '\n')
                
                seqOutPerContigWithFlanks.write('>' + specimen + '_|_' + contigName + '_|_' + seqid + '_' + viralSeqPositions[(seqid, stitle)][specimen][contigName] + '_pident_' + '{:5.3f}'.format(pident) + '_with_flanks\n')
                seqOutPerContigWithFlanks.write(viralSeqStringsWithFlanks[(seqid, stitle)][specimen][contigName] + '\n')

                contigCount += 1
        seqOut.close()
        seqOutPerContig.close()
        seqOutPerContigUnaligned.close()
        seqOutPerContigWithFlanks.close()
        
        # Cluster.
        if config['DO_CLUSTERING']:
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
    
    # blastn -task blastn -query results20k/viruses/sequences/Flaviviridae/sequences_NC_027819.1_per_contig_unaligned.fasta 
    #        -subject fasta/NC_001564.2.fasta -evalue 1e-10 -no_greedy
    # -outfmt "10 qseqid qstart qend qseq sstart send sseq evalue bitscore sseqid stitle pident"
        
    for family in set(allFamilies):
        FAMILY_DIR = Path(VIRUS_RESULTS_DIR) / SEQUENCES_DIR / family
        if not Path(FAMILY_DIR).exists():
            os.system('mkdir ' + str(FAMILY_DIR))
            
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
            virusRecord = SeqIO.read(GB_DIR + seqid + '.gb', 'gb')  # SeqRecord
            minCDSPos = math.inf
            for feature in virusRecord.features:
                if feature.type == 'CDS':
                    minCDSPos = min(minCDSPos, feature.location.start)
            startPositions[(seqid, stitle)] = minCDSPos
                    
        seqOut = open(str(Path(FAMILY_DIR) / (family + '_per_specimen_aligned.fasta')), 'w')
        seqOutPerContig = open(str(Path(FAMILY_DIR) / (family + '_per_contig_aligned.fasta')), 'w')
        seqOutPerContigUnaligned = open(str(Path(FAMILY_DIR) / (family + '_per_contig_unaligned.fasta')), 'w')
        
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
#                posSpecimen.sort()
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
#                    positions[(seqid, stitle)][specimen][contigName].sort()
                    posString = ''
                    for start, end in positions[(seqid, stitle)][specimen][contigName]:
                        posString += str(start) + '-' + str(end) + ','
                    seqOutPerContig.write('>' + region + '_' + str(num) + '__contig' + str(contigCount) + '_' + seqid + '_' + posString[:-1] + '\n')
                    seqOutPerContigUnaligned.write('>' + region + '_' + str(num) + '__contig' + str(contigCount) + '_' + seqid + '_' + viralSeqPositions[(seqid, stitle)][specimen][contigName] + '\n')
                    seqOutPerContigUnaligned.write(viralSeqStrings[(seqid, stitle)][specimen][contigName] + '\n')
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

def drawVirus(acc, family, hits, allSpecimens, separatePops, isFamily, showPlot, small):    
    
    if not Path(VIRUS_RESULTS_DIR).exists():
        os.system('mkdir ' + VIRUS_RESULTS_DIR)
    if not (Path(VIRUS_RESULTS_DIR) / DIAGRAMS_DIR).exists():
        os.system('mkdir ' + str((Path(VIRUS_RESULTS_DIR) / DIAGRAMS_DIR)))
    FAMILY_DIR = Path(VIRUS_RESULTS_DIR) / DIAGRAMS_DIR / family
    if not Path(FAMILY_DIR).exists():
        os.system('mkdir ' + str(FAMILY_DIR))
    
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
        pops.extend(list(POPULATIONS.keys()))
        
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
            
        x = range(virusRecord.sequence_length)
        y = [0] * virusRecord.sequence_length
        
        for specimen in specimens:
            for feature in newFeatures[specimen]:
                start = min(max(feature.start, 1), virusRecord.sequence_length - 1)
                end = min(feature.end, virusRecord.sequence_length - 1)
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
    writelog('  Reading hits for diagrams ...', True)
    dir = Path(SPECIMENS_RESULTS_DIR)
    subdirs = [d for d in dir.iterdir()]
    files = [d / XML_DIR / (d.name + '_hits.xml') for d in subdirs]
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
            if config['OMIT_EVES_IN_REFERENCE'] and overlap:
                continue
            virushits = contig.findall('virushit')
            v = virushits[0]  # first virus hit
            seqid = v.attrib['seqid']
            if '|' in seqid:
                seqid = seqid.split('|')[1]
            if seqid not in allVirusHits:
                allVirusHits[seqid] = {}
            if specimen not in allVirusHits[seqid]:
                allVirusHits[seqid][specimen] = DiagramHits()

            hit = []
            for virushit in virushits:
                start = int(virushit.find('sstart').text)
                end = int(virushit.find('send').text)
                hit.extend([start, end])
            primary = contig.attrib['besthit'] == 'True'  # and (len(virushits) < 5):
            allVirusHits[seqid][specimen].addHit(hit, overlap, primary, referenceOverlaps)
    return allVirusHits, allSpecimens

def drawFamily(families, famACCs, separatePops, showPlot, small):
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
                    famHits[specimen] = DiagramHits()
                allVirusHits[acc][specimen].adjust(startPositions[repACC] - startPositions[acc])
                famHits[specimen].addHits(allVirusHits[acc][specimen].getPrimaryHitsOnly())  # only primary hits so one hit per contig
                
        drawVirus(repACC, family, famHits, allSpecimens, separatePops, True, showPlot, small)
        famACCs = None
    
def drawAll(separatePops = False):
    writelog('Creating diagrams ...', True)
    allVirusHits, allSpecimens = getHitsForDiagram()
    
    fams = readFamFile()

    for acc in allVirusHits:
        hitCount = sum(hit.isPrimary() for specimen in allVirusHits[acc] for hit in allVirusHits[acc][specimen].getHits())
        if hitCount >= config['MIN_HITS_TO_SHOW_VIRUS']:
            if acc not in fams:
                fam = getFamily(acc)
                fams[acc] = fam
                addFamily(acc, fam)
            drawVirus(acc, fams[acc], allVirusHits[acc], allSpecimens, separatePops, False, True, False)

###############################################################################

def doIt(filenameBAM):
    if not Path(SPECIMEN_RESULTS_DIR).exists():
        os.system('mkdir ' + str(SPECIMEN_RESULTS_DIR))
        
    newFilenameBAM = getUnmappedReads(filenameBAM)
    spadesPath = writeReadsFASTQ(newFilenameBAM)
    success = assembleReads(spadesPath) # 1. Assemble
    if success:
        dirName = spadesPath.parent
        blastScaffolds(dirName)         # 2. Blast scaffolds
        xmlFilename = getHits(dirName)  # 3. Get hits
        if xmlFilename is not None:
            drawContigs(xmlFilename)    # 4. Draw contigs
        else:
            writelog('*** ' + str(Path(dirName).name) + ': no results found. ***', True)
    else:
        writelog('*** ' + str(Path(dirName).name) + ': assembly failed; no results found. ***', True)
        
    writelog('*** ' + str(Path(dirName).name) + ' done. ***', True)

def doAll():
#    fileList = sorted(Path(SPECIMENS_DIR).rglob('*.bam'))
    fileList = [Path('/Volumes/Data2/specimen_copies/Angola/Debug010_aegypti_Cuanda_Angola_01.LIN210A1719.sorted.deduped.merged.bam')]

    count = 0
    for file in fileList:
        count += 1
        writelog('\n' + str(count) + '/' + str(len(fileList)) + ': ' + file.name, True)
        doIt(file.resolve())
                    
# def doAllParallel2():
#     dirList = [dir for dir in Path(SPECIMENS_DIR).iterdir() if dir.is_dir() and ('combined' not in dir.name)]
#     workers = []
#     count = 0
#     for dir in dirList:
#         for file in dir.iterdir():
#             if not file.is_dir() and (file.name[-26:] == '.sorted.deduped.merged.bam'):
#                 p = Process(target = doIt, args = (str(dir.resolve()), file.name))
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
            answer = input('Results directory ' + RESULTS_DIR + ' exists.  OK to write (y/n)? ')
        if answer[0] not in 'yY':
            writelog('OK, quitting.', True)
            return 1
    else:
        os.system('mkdir ' + RESULTS_DIR)
        
    writelog('\n************************************************')
    writelog('Starting pipeline on ' + time.strftime('%c'))
    writeConfig()
    
    doAll()
    
    consolidateAll()
    drawAll(False)
#    drawFamily([('Flaviviridae', 'NC_001564.2'), ('Orthomyxoviridae', 'MF176251.1'), ('Phenuiviridae', 'NC_038263.1'), ('Rhabdoviridae', 'NC_035132.1'), ('Totiviridae', 'NC_035674.1'), ('Xinmoviridae', 'MH037149.1')],None, False, False, False, False)
    
main()
#drawAll(False)
