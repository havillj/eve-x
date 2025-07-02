#!/usr/bin/env python3

from util import *
from features import *

def removeSite(sites, posPair):
    for index in range(len(sites)):
        if sites[index][:2] == posPair:
            features = sites[index][2]
#            sites.pop(index)
            return features
    return []
    
def removeRuns(values):
    newValues = []
    prevValue = ''
    for value in values:
        if value != prevValue:
            newValues.append(value)
        prevValue = value
    return newValues
    
def getFlanks(specimen, contigName, leftStart, leftEnd, rightStart, rightEnd):
    scaffoldsFilename = str(Path(SPECIMEN_RESULTS_DIR) / specimen / SCAFFOLDS_DIR / 'scaffolds.fasta')
    scaffoldRecords = SeqIO.index(scaffoldsFilename, 'fasta')
    contigRecord = scaffoldRecords[contigName]
    seq = str(contigRecord.seq)
    return seq[leftStart - 1:leftEnd], seq[rightStart - 1:rightEnd]  # account for 1-based indexing

def getInsertSites(findFeatures, desiredSeqIDs = None, bestHitsOnly = True):
    '''bestHitsOnly must match value given to getHits that created these xml files.'''
    
    xmlFiles = sorted(Path(SPECIMEN_RESULTS_DIR).rglob('*.xml'))
    
    insertSites = {}
    
    specimen2Label = {}
    specimenCount = 1
    for file in xmlFiles:
        specimen = file.name.split('_hits')[0]
        region, num = getSpecimenLabel(specimen)
        region = region.replace(' ', '-')
        specimen2Label[specimen] = (region, num)
        
        writelog('Processing ' + str(specimenCount) + '/' + str(len(xmlFiles)) + ': ' + specimen, VERBOSE)
        
        specimenCount += 1
        
        # Parse XML tree.
        
        tree = ET.parse(str(file))
        root = tree.getroot()
        skip = False
        for element in root.iter():   # one per virus per contig
            if element.tag == 'contig':
                contig = element
                contigName = element.attrib['name']
                skip = (bestHitsOnly and (element.attrib['besthit'] != 'True'))
                if not skip:
                    hit = {}
                    vectorHitsLeft = []
                    vectorHitsRight = []
#                    referenceOverlaps = eval(element.attrib['referenceoverlaps'])
            elif not skip and element.tag == 'virushit':
                seqid = element.attrib['seqid']            # should be same for all v
                if '|' in seqid:
                    seqid = seqid.split('|')[1]
                skip = (desiredSeqIDs is not None) and (seqid not in desiredSeqIDs)
                if not skip:
                    if seqid not in insertSites:
                        insertSites[seqid] = {}
                    if specimen not in insertSites[seqid]:
                        insertSites[seqid][specimen] = {}
                    insertSites[seqid][specimen][contigName] = {'match': [], 'inversion': [], 'left': [], 'right': [], 'referenceoverlaps': []}
                    v = element
                    hit[(int(v.find('qstart').text), int(v.find('qend').text))] = (int(v.find('sstart').text), int(v.find('send').text))
                    qpos = list(hit.keys())
                    qpos.sort()
                    virus_qstart = qpos[0][0]
                    virus_qend = qpos[-1][1]
                    virus_sstart = hit[qpos[0]][0]
                    virus_send = hit[qpos[-1]][1]
            elif not skip and element.tag == 'vectorhitleft':
                vectorHitsLeft.append(element)
                v = element
                host_qstart = int(v.find('qstart').text)
                host_qend = int(v.find('qend').text)
                if (virus_qstart - host_qend <= config['MAX_HOST_VIRUS_DISTANCE']) or (host_qstart - virus_qend <= config['MAX_HOST_VIRUS_DISTANCE']):
                    chrSeqID = v.attrib['seqid']
                    sstart = int(v.find('sstart').text)
                    send = int(v.find('send').text)
                    if sstart < send:
                        if findFeatures:
                            site = (chrSeqID, send, searchiTrees(chrSeqID, send - config['FEATURE_SEARCH_DIST'], send), virus_qstart - host_qend, (virus_sstart, virus_send))
                        else:
                            site = (chrSeqID, send, [], virus_qstart - host_qend, (virus_sstart, virus_send))
                        insertSites[seqid][specimen][contigName]['left'].append(site + ((host_qstart, host_qend),))
                    else:  # this is actually a right flanking region
                        if findFeatures:
                            site = (chrSeqID, send, searchiTrees(chrSeqID, send, send + config['FEATURE_SEARCH_DIST']), virus_qstart - host_qend, (virus_send, virus_sstart))
                        else:
                            site = (chrSeqID, send, [], virus_qstart - host_qend, (virus_send, virus_sstart))
                        insertSites[seqid][specimen][contigName]['right'].append(site + ((host_qstart, host_qend),))
            elif not skip and element.tag == 'vectorhitright':
                vectorHitsRight.append(element)
                v = element
                host_qstart = int(v.find('qstart').text)
                host_qend = int(v.find('qend').text)
                if (virus_qstart - host_qend <= config['MAX_HOST_VIRUS_DISTANCE']) or (host_qstart - virus_qend <= config['MAX_HOST_VIRUS_DISTANCE']):
                    chrSeqID = v.attrib['seqid']
                    sstart = int(v.find('sstart').text)
                    send = int(v.find('send').text)
                    if sstart < send:
                        if findFeatures:
                            site = (chrSeqID, sstart, searchiTrees(chrSeqID, sstart, sstart + config['FEATURE_SEARCH_DIST']), host_qstart - virus_qend, (virus_sstart, virus_send))
                        else:
                            site = (chrSeqID, sstart, [], host_qstart - virus_qend, (virus_sstart, virus_send))
                        insertSites[seqid][specimen][contigName]['right'].append(site + ((host_qstart, host_qend),))
                    else:  # this is actually a left flanking region
                        if findFeatures:
                            site = (chrSeqID, sstart, searchiTrees(chrSeqID, sstart - config['FEATURE_SEARCH_DIST'], sstart), host_qstart - virus_qend, (virus_send, virus_sstart))
                        else:
                            site = (chrSeqID, sstart, [], host_qstart - virus_qend, (virus_send, virus_sstart))
                        insertSites[seqid][specimen][contigName]['left'].append(site + ((host_qstart, host_qend),))
            elif not skip and element.tag == 'vectorhitoverlap':
                v = element
                host_qstart = int(v.find('qstart').text)
                host_qend = int(v.find('qend').text)
                chrSeqID = v.attrib['seqid']
                sstart = int(v.find('sstart').text)
                send = int(v.find('send').text)
                overlap = None
                flanks = ['', '', virus_qstart - host_qstart, host_qend - virus_qend]  # dists beyond ends of EVE; negative means doesn't cover that end of EVE
                if sstart > send:
                    sstart, send = send, sstart
                    flanks[2], flanks[3] = flanks[3], flanks[2]
                    flanks = tuple(flanks)
                    vse = (virus_send, virus_sstart)
                else:
                    vse = (virus_sstart, virus_send)
                if findFeatures:
                    site = (chrSeqID, (sstart, send), searchiTrees(chrSeqID, sstart, send), overlap, flanks, vse)
                else:
                    site = (chrSeqID, (sstart, send), [], overlap, flanks, vse)
                insertSites[seqid][specimen][contigName]['referenceoverlaps'].append(site + ((host_qstart, host_qend),))
            elif not skip and (element.tag == 'match'):
                match = element
                for v in vectorHitsLeft:
                    if v.attrib['id'] == match.attrib['leftid']:
                        chrSeqID = v.attrib['seqid']
                        leftStart = int(v.find('sstart').text)
                        leftEnd = int(v.find('send').text)
#                        leftQStart = int(v.find('qstart').text)
                        leftQEnd = int(v.find('qend').text)
                        break
                for v in vectorHitsRight:
                    if v.attrib['id'] == match.attrib['rightid']:
                        rightStart = int(v.find('sstart').text)
                        rightEnd = int(v.find('send').text)
                        rightQStart = int(v.find('qstart').text)
#                        rightQEnd = int(v.find('qend').text)
                        break
                
                if abs(rightStart - leftEnd) > config['MAX_FLANK_DISTANCE']:
                    continue
                 
                if (leftStart < leftEnd) and (rightStart < rightEnd):  # host hits on forward strand
                    if rightStart < leftEnd:  # overlapping flanking region
                        overlap = (rightStart, leftEnd)
                        flankLength = leftEnd - rightStart + 1
                        flanks = getFlanks(specimen, contigName.split('__')[0], leftQEnd - flankLength + 1, leftQEnd, rightQStart, rightQStart + flankLength - 1)
                        flanks += (virus_qstart - leftQEnd, rightQStart - virus_qend)
                    else:
                        overlap = None
                        flanks = ('', '', virus_qstart - leftQEnd, rightQStart - virus_qend)
#                     featuresLeft = removeSite(insertSites[seqid][specimen][contigName]['left'], (v.attrib['seqid'], leftEnd))
#                     featuresRight = removeSite(insertSites[seqid][specimen][contigName]['right'], (v.attrib['seqid'], rightStart))
#                     features = (featuresLeft, featuresRight, searchiTrees(chrSeqID, leftEnd, rightStart))
                    vse = (virus_sstart, virus_send)  # virus start end
                elif (leftStart > leftEnd) and (rightStart > rightEnd):  # host hits on reverse strand
                    if rightStart > leftEnd:  # overlapping flanking region
                        overlap = (leftEnd, rightStart)  # convert everything to forward strand
                        flankLength = rightStart - leftEnd + 1
                        flanks = getFlanks(specimen, contigName.split('__')[0], leftQEnd - flankLength + 1, leftQEnd, rightQStart, rightQStart + flankLength - 1)
                        flanks = (reverseComplement(flanks[1]), reverseComplement(flanks[0]), rightQStart - virus_qend, virus_qstart - leftQEnd)
                    else:
                        overlap = None
                        flanks = ('', '', rightQStart - virus_qend, virus_qstart - leftQEnd)
#                     featuresLeft = removeSite(insertSites[seqid][specimen][contigName]['left'], (v.attrib['seqid'], rightStart))
#                     featuresRight = removeSite(insertSites[seqid][specimen][contigName]['right'], (v.attrib['seqid'], leftEnd))
#                     features = (featuresLeft, featuresRight, searchiTrees(chrSeqID, leftEnd, rightStart))
                    leftEnd, rightStart = rightStart, leftEnd
                    vse = (virus_send, virus_sstart)
                else:
                    print('Error!')
                    exit()
                featuresLeft = removeSite(insertSites[seqid][specimen][contigName]['left'], (v.attrib['seqid'], leftEnd))
                featuresRight = removeSite(insertSites[seqid][specimen][contigName]['right'], (v.attrib['seqid'], rightStart))
                features = (featuresLeft, featuresRight, searchiTrees(chrSeqID, leftEnd, rightStart))
                insertSites[seqid][specimen][contigName][element.tag].append((v.attrib['seqid'], (leftEnd, rightStart), features, overlap, flanks, vse))
            elif not skip and (element.tag == 'inversion'):
                match = element
                for v in vectorHitsLeft:
                    if v.attrib['id'] == match.attrib['leftid']:
                        chrSeqID = v.attrib['seqid']
                        leftStart = int(v.find('sstart').text)
                        leftEnd = int(v.find('send').text)
#                        leftQStart = int(v.find('qstart').text)
                        leftQEnd = int(v.find('qend').text)
                        featuresLeft = removeSite(insertSites[seqid][specimen][contigName]['left'], (v.attrib['seqid'], leftEnd)) + \
                                       removeSite(insertSites[seqid][specimen][contigName]['right'], (v.attrib['seqid'], leftEnd))
                        break
                for v in vectorHitsRight:
                    if v.attrib['id'] == match.attrib['rightid']:
                        rightStart = int(v.find('sstart').text)
                        rightEnd = int(v.find('send').text)
                        rightQStart = int(v.find('qstart').text)
#                        rightQEnd = int(v.find('qend').text)
                        featuresRight = removeSite(insertSites[seqid][specimen][contigName]['left'], (v.attrib['seqid'], rightStart)) + \
                                        removeSite(insertSites[seqid][specimen][contigName]['right'], (v.attrib['seqid'], rightStart))
                        break
                                
                if abs(rightStart - leftEnd) > config['MAX_FLANK_DISTANCE']:
                    continue
                 
                overlap = None
                flanks = ('', '', virus_qstart - leftQEnd - 1, rightQStart - virus_qend - 1)
                if leftEnd > rightStart:
                    leftEnd, rightStart = rightStart, leftEnd
                    features = (featuresRight, featuresLeft, searchiTrees(chrSeqID, leftEnd, rightStart))
                    vse = (virus_send, virus_sstart)
                else:
                    features = (featuresLeft, featuresRight, searchiTrees(chrSeqID, leftEnd, rightStart))
                    vse = (virus_sstart, virus_send)
                
                insertSites[seqid][specimen][contigName][element.tag].append((v.attrib['seqid'], (leftEnd, rightStart), features, overlap, flanks, vse))

    return insertSites, specimen2Label
       
def drawInsertSites(clusteredFileName, findFeatures, bestHitsOnly = True):
    if not Path(clusteredFileName).exists():
        writelog('File does not exist: ' + clusteredFileName, VERBOSE)
        return
        
    unalignedFastaName = str(Path(clusteredFileName).resolve().parent.parent / (Path(clusteredFileName).name.split('_aligned_clustered')[0] + '_unaligned.fasta'))
    withFlanksFastaName = str(Path(clusteredFileName).resolve().parent.parent / (Path(clusteredFileName).name.split('_aligned_clustered')[0] + '_with_flanks.fasta'))
        
    foundFirst = False
    seqids = set()
    for record in SeqIO.parse(clusteredFileName, 'fasta'):
        if not foundFirst:
            foundFirst = True
            continue
        parts = record.description.split('_|_')
        if len(parts) < 4:  # no cluster ids
            seqid = parts[2].split('.')[0] 
            seqid += parts[2].split(seqid)[1][:2]
        else:
            seqid = parts[3].split('.')[0]
            seqid += parts[3].split(seqid)[1][:2]
        seqids.add(seqid)
    
    insertSites, specimen2Label = getInsertSites(findFeatures, seqids, bestHitsOnly)
    
    writelog('Consolodating insertion site data...', VERBOSE)
    
    y = 1           # leave a blank y tick
    yvalues = {}
    xvalues = {}
    ylabels  = []
    colors = {}
    positions = {}
    features = {}
    pidents = {}
    totalClusterPidents = {}
    totalClusterSize = {}
    clusterCounts = {}
    clusterSeqs = {}
        
    numRecords = 0
    foundFirst = False
    for record in SeqIO.parse(clusteredFileName, 'fasta'):
        if not foundFirst:
            foundFirst = True
            continue
        desc = record.description
        parts = desc.split('_|_')
        if len(parts) < 4:  # no cluster ids
            parts.insert(0, 'C1')
        cluster = parts[0]
        specimen = parts[1]
        contig = parts[2]
        seqid = parts[3].split('.')[0]
        seqid += parts[3].split(seqid)[1][:2]
        pident = float(parts[3].split('_pident_')[1])
        
        numRecords += 1
        region, num = specimen2Label[specimen]
        ylabels.append(cluster + ' ' + region + ' ' + str(num) + ' ' + contig.split('_length')[0])
        
        if cluster not in yvalues:
            yvalues[cluster] = {chr: [] for chr in CHR_NAMES}
            xvalues[cluster] = {chr: [] for chr in CHR_NAMES}
            colors[cluster] = {chr: [] for chr in CHR_NAMES}
            positions[cluster] = {chr: {} for chr in CHR_NAMES}
            features[cluster] = {chr: {} for chr in CHR_NAMES}
            pidents[cluster] = {chr: {} for chr in CHR_NAMES}
            totalClusterPidents[cluster] = 0
            totalClusterSize[cluster] = 0
            clusterCounts[cluster] = {}
            clusterSeqs[cluster] = []
#            ylabels[cluster] = []
#            y[cluster] = 0

        clusterSeqs[cluster].append('_|_'.join(parts[1:]))
        
        if region not in clusterCounts[cluster]:
            clusterCounts[cluster][region] = set()
        clusterCounts[cluster][region].add(num)
        
        totalClusterPidents[cluster] += pident
        totalClusterSize[cluster] += 1
        
        hostContigIntervals = IntervalTree()
        hostContigIntervalCounts = {}
        for chromID, position, feat, dist, vse, hostqpos in insertSites[seqid][specimen][contig]['left']:
            if chromID in CHR_NAMES:
                if hostqpos not in hostContigIntervalCounts:
                    hostContigIntervalCounts[hostqpos] = 1
                    hostContigIntervals.addi(*hostqpos)
                else:
                    hostContigIntervalCounts[hostqpos] += 1
                
        for chromID, position, feat, dist, vse, hostqpos in insertSites[seqid][specimen][contig]['left']:
            if chromID in CHR_NAMES:
                if (hostContigIntervalCounts[hostqpos] > config['MAX_REPEATS_TO_SHOW_POSITION']) or \
                   (sum([hostContigIntervalCounts[(iv.begin, iv.end)] for iv in hostContigIntervals[hostqpos[0]:hostqpos[1]]]) > config['MAX_REPEATS_TO_SHOW_POSITION']):
                    continue
                yvalues[cluster][chromID].append(y)
                xvalues[cluster][chromID].append(position)
                colors[cluster][chromID].append('#5e5e5e')    # dark gray (was green/asparagus '#87a96b')
                if position not in positions[cluster][chromID]:
                    positions[cluster][chromID][position] = []
                    features[cluster][chromID][position] = []
                    pidents[cluster][chromID][position] = 0
                positions[cluster][chromID][position].append((region + '_' + str(num), contig, 'left', dist, None, vse))
                features[cluster][chromID][position].append((feat, [], []))
                pidents[cluster][chromID][position] += pident
                
        hostContigIntervals = IntervalTree()
        hostContigIntervalCounts = {}
        for chromID, position, feat, dist, vse, hostqpos in insertSites[seqid][specimen][contig]['right']:
            if chromID in CHR_NAMES:
                if hostqpos not in hostContigIntervalCounts:
                    hostContigIntervalCounts[hostqpos] = 1
                    hostContigIntervals.addi(*hostqpos)
                else:
                    hostContigIntervalCounts[hostqpos] += 1
                    
        for chromID, position, feat, dist, vse, hostqpos in insertSites[seqid][specimen][contig]['right']:
            if chromID in CHR_NAMES:
                if (hostContigIntervalCounts[hostqpos] >  config['MAX_REPEATS_TO_SHOW_POSITION']) or \
                   (sum([hostContigIntervalCounts[(iv.begin, iv.end)] for iv in hostContigIntervals[hostqpos[0]:hostqpos[1]]]) > config['MAX_REPEATS_TO_SHOW_POSITION']):
                    continue
                yvalues[cluster][chromID].append(y)
                xvalues[cluster][chromID].append(position)
                colors[cluster][chromID].append('#5e5e5e')       # dark gray (was red '#ff2600')
                if position not in positions[cluster][chromID]:
                    positions[cluster][chromID][position] = []
                    features[cluster][chromID][position] = []
                    pidents[cluster][chromID][position] = 0
                positions[cluster][chromID][position].append((region + '_' + str(num), contig, 'right', dist, None, vse))
                features[cluster][chromID][position].append(([], feat, []))
                pidents[cluster][chromID][position] += pident
        for chromID, position, feat, overlap, flanks, vse in insertSites[seqid][specimen][contig]['match']:
            if chromID in CHR_NAMES:
                yvalues[cluster][chromID].append(y)
                xvalues[cluster][chromID].append(position[0])
                colors[cluster][chromID].append('#0000ff')       # blue (was magenta '#ff40ff')
                if position not in positions[cluster][chromID]:
                    positions[cluster][chromID][position] = []
                    features[cluster][chromID][position] = []
                    pidents[cluster][chromID][position] = 0
                positions[cluster][chromID][position].append((region + '_' + str(num), contig, 'match', overlap, flanks, vse))
                features[cluster][chromID][position].append(feat)
                pidents[cluster][chromID][position] += pident
        for chromID, position, feat, overlap, flanks, vse in insertSites[seqid][specimen][contig]['inversion']:
            if chromID in CHR_NAMES:
                yvalues[cluster][chromID].append(y)
                xvalues[cluster][chromID].append(position[0])
                colors[cluster][chromID].append('#ff2600')       # red (was dark gray '#5e5e5e')
                if position not in positions[cluster][chromID]:
                    positions[cluster][chromID][position] = []
                    features[cluster][chromID][position] = []
                    pidents[cluster][chromID][position] = 0
                positions[cluster][chromID][position].append((region + '_' + str(num), contig, 'inversion', overlap, flanks, vse))
                features[cluster][chromID][position].append(feat)
                pidents[cluster][chromID][position] += pident
        for chromID, position, feat, overlap, flanks, vse in insertSites[seqid][specimen][contig]['referenceoverlaps']:
            if chromID in CHR_NAMES:
                yvalues[cluster][chromID].append(y)
                xvalues[cluster][chromID].append(position[0])
                colors[cluster][chromID].append('#0000ff')   # same as match
                if position not in positions[cluster][chromID]:
                    positions[cluster][chromID][position] = []
                    features[cluster][chromID][position] = []
                    pidents[cluster][chromID][position] = 0
                positions[cluster][chromID][position].append((region + '_' + str(num), contig, 'referenceoverlap', overlap, flanks, vse))
                features[cluster][chromID][position].append(([], [], feat))
                pidents[cluster][chromID][position] += pident
        y += 1
        
    minChromLength = min(CHR_LENGTHS.values())
    widthRatios = [CHR_LENGTHS[chr] / minChromLength for chr in CHR_LENGTHS]
    fig, ax = pyplot.subplots(1, 3, sharey = True, figsize = (7, 0.045 * (numRecords + 1)), gridspec_kw = {'width_ratios': widthRatios})
    fig.subplots_adjust(top = 0.97, bottom = 0.03, left = 0.25, right = 0.99, wspace = 0.0)
    
    ax[0].set_yticks(range(numRecords + 1))   # leave position 0 empty for readability
    ax[0].set_yticklabels([''] + ylabels)
    ax[0].set_ylim(0, numRecords + 1)
    ax[0].invert_yaxis()
    
    chrCount = 0
    xtickPositions = [list(range(0, CHR_LENGTHS[chr], int(1e8))) for chr in CHR_LENGTHS]
    xtickLabels = [[str(n//int(1e6))+'M' for n in range(0, CHR_LENGTHS[chr], int(1e8))] for chr in CHR_LENGTHS]
    for chromID in CHR_NAMES:
        for cluster in yvalues:
            ax[chrCount].scatter(xvalues[cluster][chromID], yvalues[cluster][chromID], color=colors[cluster][chromID], s = 1)
        ax[chrCount].tick_params(labelsize = 3)
        if chrCount > 0:
            ax[chrCount].tick_params(left=False)
        ax[chrCount].set_xticks(xtickPositions[chrCount])
        ax[chrCount].set_xticklabels(xtickLabels[chrCount])
        ax[chrCount].set_xlim(0, CHR_LENGTHS[chromID])
        ax[chrCount].set_xlabel('Chromosome ' + str(chrCount + 1), fontsize=4)
        chrCount += 1

    if not (Path(clusteredFileName).parent / seqid).exists():
        os.system('mkdir ' + str(Path(clusteredFileName).parent / seqid))
    saveFileName = str(Path(clusteredFileName).parent / seqid / (seqid + '_insertpositions'))
    if len(seqids) > 1:
        saveFileName += '_plus'
    pp = PdfPages(saveFileName + '.pdf')
    pp.savefig(fig)
    pp.close()
    pyplot.close(fig)
    
    groupLabels = {}
    nextGroupLabel = {}
    groups = {}
    clusterPidents = {}
    for chromID in CHR_NAMES:
        tsvPositionFile = open(saveFileName + '_' + CHR_NAMES[chromID] + '.tsv', 'w')
        tsvPositionFile.write('Cluster\tStart\tEnd\tRepeats\tIntron\tExon\t5\' Repeats\t5\' Intron\t5\' Exon\t3\' Repeats\t3\' Intron\t3\' Exon\tGroup\t% Identity\t' + '\t'.join(POPULATIONS.keys()) + '\tSpecimens\n')
        for cluster in positions:
            if cluster not in groupLabels:
                groupLabels[cluster] = {}
                nextGroupLabel[cluster] = 1
                groups[cluster] = {}
            sortedPositions = list(positions[cluster][chromID].keys())
            sortedPositions.sort(key=lambda x: int(x[0]) if isinstance(x, tuple) else int(x))
            for position in sortedPositions:
                if len(positions[cluster][chromID][position]) >= config['MIN_HITS_TO_SHOW_POSITION']:
                    groupNames = tuple([t[0] for t in positions[cluster][chromID][position]])
                    if groupNames not in groupLabels[cluster]:
                        groupLabels[cluster][groupNames] = nextGroupLabel[cluster]
                        groups[cluster][nextGroupLabel[cluster]] = [groupNames, []]
                        nextGroupLabel[cluster] += 1
                    
                    featDict = {'left': {'repeat_region': [], 'intron': [], 'exon': []}, 'right': {'repeat_region': [], 'intron': [], 'exon': []}, 'inregion': {'repeat_region': [], 'intron': [], 'exon': []}}
                    for featLeft, featRight, featIn in features[cluster][chromID][position]:
                        for featType, featID, featPos in featLeft:
                            if featType not in featDict['left']:
                                featDict['left'][featType] = []
                            featDict['left'][featType].append((featType, featID, featPos))
                        for featType, featID, featPos in featRight:
                            if featType not in featDict['right']:
                                featDict['right'][featType] = []
                            featDict['right'][featType].append((featType, featID, featPos))
                        for featType, featID, featPos in featIn:
                            if featType not in featDict['inregion']:
                                featDict['inregion'][featType] = []
                            featDict['inregion'][featType].append((featType, featID, featPos))
                            
                        for featType in featDict['left']:
                            featDict['left'][featType].sort(key=lambda t: t[2].end, reverse=True)  # closest first
                            featDict['left'][featType] = removeRuns(featDict['left'][featType])
                        for featType in featDict['right']:
                            featDict['right'][featType].sort(key=lambda t: t[2].start)  # closest first
                            featDict['right'][featType] = removeRuns(featDict['right'][featType])
                         
                    featStrings = {'left': {}, 'right': {}, 'inregion': {}}
                    for featType in featDict['left']:
                        if featType in ['gene', 'ncRNA_gene', 'pseudogene']:
                            if len(featDict['left']['exon']) > 0:  # omit this; just show exon
                                continue
                            else:                      # position is in an intron
                                featStrings['left']['intron'] = ', '.join([t[1] for t in featDict['left'][featType]])  # feature ids
                        else:
                            featStrings['left'][featType] = ', '.join([t[1] for t in featDict['left'][featType]])
                    for featType in featDict['right']:
                        if featType in ['gene', 'ncRNA_gene', 'pseudogene']:
                            if len(featDict['right']['exon']) > 0:  # omit this; just show exon
                                continue
                            else:                      # position is in an intron
                                featStrings['right']['intron'] = ', '.join([t[1] for t in featDict['right'][featType]])
                        else:
                            featStrings['right'][featType] = ', '.join([t[1] for t in featDict['right'][featType]])
                    for featType in featDict['inregion']:
                        if featType in ['gene', 'ncRNA_gene', 'pseudogene']:
                            if len(featDict['inregion']['exon']) > 0:  # omit this; just show exon
                                continue
                            else:                      # position is in an intron
                                featStrings['inregion']['intron'] = ', '.join(set([t[1] for t in featDict['inregion'][featType]])) # order doesn't matter here
                        else:
                            featStrings['inregion'][featType] = ', '.join(set([t[1] for t in featDict['inregion'][featType]]))
                            
                    if len(featDict['left']['repeat_region']) > 0:
                        upstreamRepeat = featDict['left']['repeat_region'][0][1]
                    else:
                        upstreamRepeat = ''
                    if len(featDict['right']['repeat_region']) > 0:
                        downstreamRepeat = featDict['right']['repeat_region'][0][1]
                    else:
                        downstreamRepeat = ''
                    if len(featDict['inregion']['repeat_region']) > 0:
                        inregionRepeats = featStrings['inregion']['repeat_region']
                    else:
                        inregionRepeats = ''
                    groups[cluster][groupLabels[cluster][groupNames]][1].append((CHR_NAMES[chromID], position, (upstreamRepeat, downstreamRepeat, inregionRepeats, featStrings['inregion']['intron']), positions[cluster][chromID][position]))
                    
                    if isinstance(position, tuple):
                        start = str(position[0])
                        end = str(position[1])
                    else:
                        start = str(position)
                        end = ''
                    tsvPositionFile.write(cluster[1:] + '\t' + start + '\t' + end + '\t' + featStrings['inregion']['repeat_region'] + '\t' + featStrings['inregion']['intron'] + '\t' + featStrings['inregion']['exon'] + '\t' + featStrings['left']['repeat_region'] + '\t' + featStrings['left']['intron'] + '\t' + featStrings['left']['exon'] + '\t' + featStrings['right']['repeat_region'] + '\t' + featStrings['right']['intron'] + '\t' + featStrings['right']['exon'] + '\t' + str(groupLabels[cluster][groupNames]) + '\t' + '{:5.2f}'.format(pidents[cluster][chromID][position] / len(positions[cluster][chromID][position])) + '\t')
                    
                    countsByRegion = {region: 0 for region in POPULATIONS}
                    line = ''
                    for specimen, contig, type, overlap, flanks, vse in positions[cluster][chromID][position]:
                        region = specimen.split('_')[0]
                        countsByRegion[region] += 1
                        if isinstance(position, tuple):
                            if overlap is None:
                                line += specimen + ' (' + type + ')\t'   # include match/inversion if tuple
                            else:
                                line += specimen + ' (' + type + ', ' + str(overlap) + ', ' + str(flanks) + ')\t'
                        else:
                            line += specimen + ' (' + type + ', dist = ' + str(overlap) + ')\t'  # overlap == distance
                    for region in countsByRegion:
                        tsvPositionFile.write(str(countsByRegion[region]) + '\t')
                    tsvPositionFile.write(line[:-1] + '\n')
        tsvPositionFile.close()
        
    seqs = SeqIO.index(unalignedFastaName, 'fasta')
    seqsFlanks = SeqIO.index(withFlanksFastaName, 'fasta')
    textFile = open(saveFileName + '.txt', 'w')
        
    textFile.write('Cluster statistics\n\n')
    for cluster in clusterCounts:
        clusterRecords = []
        clusterRecordsFlanks = []
        eveLengthFreq = {}
        for desc in clusterSeqs[cluster]:
            try:
                clusterRecords.append(seqs[desc])
                clusterRecordsFlanks.append(seqsFlanks[desc + '_with_flanks'])
                eveLength = len(seqs[desc].seq)
                if eveLength not in eveLengthFreq:
                    eveLengthFreq[eveLength] = 0
                eveLengthFreq[eveLength] += 1
            except:  # description may have been adjusted
                for desc2 in seqs:
                    if desc.split('_|_')[:2] == desc2.split('_|_')[:2]:
                        clusterRecords.append(seqs[desc2])
                        clusterRecordsFlanks.append(seqsFlanks[desc2 + '_with_flanks'])
                        eveLength = len(seqs[desc2].seq)
                        if eveLength not in eveLengthFreq:
                            eveLengthFreq[eveLength] = 0
                        eveLengthFreq[eveLength] += 1
                        break
        SeqIO.write(clusterRecords, str(Path(clusteredFileName).parent / seqid / (seqid + '_' + cluster + '.fasta')), 'fasta')
        SeqIO.write(clusterRecordsFlanks, str(Path(clusteredFileName).parent / seqid / (seqid + '_' + cluster + '_with_flanks.fasta')), 'fasta')
        
        textFile.write('  Cluster ' + cluster + '\n')
        textFile.write('    EVE length distribution:\n')
        eveLengths = list(eveLengthFreq.keys())
        eveLengths.sort()
        for eveLength in eveLengths:
            textFile.write('      {0:>4}: {1:>4}\n'.format(eveLength, eveLengthFreq[eveLength]))
#         textFile.write('    Average EVE length: {:<7.2f}\n'.format(sum(eveLengths) / len(eveLengths)))
#         textFile.write('    Minimum EVE length: {:}\n'.format(min(eveLengths)))
#         textFile.write('    Maximum EVE length: {:}\n\n'.format(max(eveLengths)))
        textFile.write('    Percent identity: {:5.2f}%\n\n'.format(totalClusterPidents[cluster] / totalClusterSize[cluster]))
        regionsSorted = list(clusterCounts[cluster].keys())
        regionsSorted.sort()
        for region in regionsSorted:
            textFile.write('    {0:<20}: {1:>3}\n'.format(region, len(clusterCounts[cluster][region])))
        textFile.write('\n')
            
    textFile.write('\nInsertion site groups:\n')
    for cluster in groups:
        textFile.write('\n  Cluster ' + cluster + '\n')
        for groupLabel in groups[cluster]:
            textFile.write('    Group ' + str(groupLabel) + ' (size ' + str(len(groups[cluster][groupLabel][0])) + '): ' + ', '.join(groups[cluster][groupLabel][0]) + '\n')
    
    textFile.write('\nInsertion positions by group (sorted by group size)\n')
    for cluster in positions:
        textFile.write('\n  Cluster ' + str(cluster) + '\n')
        groupLabelsSorted = list(groups[cluster].keys())
        groupLabelsSorted.sort(key = lambda label: len(groups[cluster][label][0]), reverse = True)
        for groupLabel in groupLabelsSorted:
            textFile.write('    Group ' + str(groupLabel) + ' (size ' + str(len(groups[cluster][groupLabel][0])) + ')\n      ') # + ' (' + ', '.join(groups[cluster][groupLabel][0]) + ')\n      ')
            for chr, pos, feat, hits in groups[cluster][groupLabel][1]:
                if isinstance(pos, tuple):
                    textFile.write('{0}: {1:,}--{2:,} ({3})'.format(chr, pos[0], pos[1], hits[0][2]) + '; ')
                else:
                    textFile.write('{0}: {1:,}'.format(chr, pos) + '; ')
            textFile.write('\n')
            
    textFile.write('\nMatching flanks only by group (sorted by group size)\n')
    for cluster in positions:
        textFile.write('\n  Cluster ' + str(cluster) + '\n')
        groupLabelsSorted = list(groups[cluster].keys())
        groupLabelsSorted.sort(key = lambda label: len(groups[cluster][label][0]), reverse = True)
        for groupLabel in groupLabelsSorted:
            found2 = False
            for chr, pos, feat, hits in groups[cluster][groupLabel][1]:
#                if isinstance(pos[1], tuple):
                found1 = False
                for hit in hits:
                    if (hit[2] in ['match', 'inversion', 'referenceoverlap']) and (abs(pos[0] - pos[1]) <= config['MAX_FLANK_DISTANCE']):
                        if not found2:
                            textFile.write('    Group ' + str(groupLabel) + ' (size ' + str(len(groups[cluster][groupLabel][0])) + ')\n')
                            found2 = True
                        if not found1:
                            textFile.write('      {0}: {1:,}--{2:,} {3}\n'.format(chr, pos[0], pos[1], ', '.join(feat[2:])))  # inregion features only
                            found1 = True
                        if hit[3] is not None:
#                            textFile.write('        {0:<20} {1:<10} {2:,} - {3:,} {4} ({5}, {6})\n'.format(hit[0], hit[2], hit[3][0], hit[3][1], ', '.join(hit[4][:2]), hit[4][2], hit[4][3]))
                            textFile.write('        {0:<20} {1:<10} {2} ({3}, {4})  virus {5}-{6}\n'.format(hit[0], hit[2], ', '.join(hit[4][:2]), hit[4][2], hit[4][3], hit[5][0], hit[5][1]))
                        else:
                            textFile.write('        {0:<20} {1:<10} ({2}, {3})  virus {4}-{5}\n'.format(hit[0], hit[2], hit[4][2], hit[4][3], hit[5][0], hit[5][1]))
            if found2:
                textFile.write('\n\n')
            
#     textFile.write('\nPercent identities\n')
#     for cluster in positions:
#         textFile.write('  Cluster ' + str(cluster) + ': {:5.2f}%\n'.format(totalClusterPidents[cluster] / totalClusterSize[cluster]))
                
def usage():
    print('Usage: ' + sys.argv[0] + ' [--find_features] /path/to/clustered_filename.fasta')

def main():
    if len(sys.argv) < 2:
        usage()
        return
        
    findFeatures = False
    for arg in sys.argv[1:-1]:
        if arg == '--find_features':
            findFeatures = True
        else:
            usage()
            return
    
    if sys.argv[-1][0] == '-':
        usage()
        return
       
    if findFeatures: 
        makeIntervalTrees()
    drawInsertSites(sys.argv[-1], findFeatures)
    
if __name__ == '__main__':
    main()
