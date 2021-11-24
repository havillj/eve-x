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
    scaffoldsFilename = str(Path(SPECIMENS_DIR) / (specimen + '/' + SPECIMEN_RESULTS_DIR + '/scaffolds/scaffolds.fasta'))
    scaffoldRecords = SeqIO.index(scaffoldsFilename, 'fasta')
    contigRecord = scaffoldRecords[contigName]
    seq = str(contigRecord.seq)
    return seq[leftStart - 1:leftEnd], seq[rightStart - 1:rightEnd]  # account for 1-based indexing

def getInsertSites(desiredSeqIDs = None, bestHitsOnly = True):
    '''bestHitsOnly must match value given to getHits that created these xml files.'''
    
    xmlFiles = []
    specimensPath = Path(SPECIMENS_DIR)
    for specimenDir in specimensPath.iterdir():
        if specimenDir.name[0] != '.':
            for file in (specimenDir / (SPECIMEN_RESULTS_DIR + '/xml')).iterdir():
#                if (file.name[0] != '.') and (('_hits_features.xml' in file.name) or (('_hits.xml' in file.name) and not Path(str(file)[:-4] + '_features.xml').exists())):
                if (file.name[0] != '.') and file.name.endswith('_hits.xml'):
                    xmlFiles.append(file.resolve())
    xmlFiles.sort()
    
    insertSites = {}
    
    specimen2Label = {}
    specimenCount = 1
    for file in xmlFiles:
        specimen = file.name.split('_hits')[0]
        region, num = getSpecimenLabel(specimen)
        region = region.replace(' ', '-')
        specimen2Label[specimen] = (region, num)
        
        writelog('  Processing ' + str(specimenCount) + '/' + str(len(xmlFiles)) + ': ' + specimen, True)
        
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
#                if (-ALLOWED_OVERLAP <= virus_qstart - host_qend <= config['MAX_HOST_VIRUS_DISTANCE'])) or (-ALLOWED_OVERLAP <= host_qstart - virus_qend <= config['MAX_HOST_VIRUS_DISTANCE'])):
                if (virus_qstart - host_qend <= config['MAX_HOST_VIRUS_DISTANCE'])) or (host_qstart - virus_qend <= config['MAX_HOST_VIRUS_DISTANCE'])):
                    chrSeqID = v.attrib['seqid']
                    sstart = int(v.find('sstart').text)
                    send = int(v.find('send').text)
                    if sstart < send:
                        if FIND_FEATURES:
                            site = (chrSeqID, send, searchiTrees(chrSeqID, send - config['FEATURE_SEARCH_DIST'], send), virus_qstart - host_qend, (virus_sstart, virus_send))
                        else:
                            site = (chrSeqID, send, [], virus_qstart - host_qend, (virus_sstart, virus_send))
                        insertSites[seqid][specimen][contigName]['left'].append(site)
                    else:  # this is actually a right flanking region
                        if FIND_FEATURES:
                            site = (chrSeqID, send, searchiTrees(chrSeqID, send, send + config['FEATURE_SEARCH_DIST']), virus_qstart - host_qend, (virus_send, virus_sstart))
                        else:
                            site = (chrSeqID, send, [], virus_qstart - host_qend, (virus_send, virus_sstart))
                        insertSites[seqid][specimen][contigName]['right'].append(site)
            elif not skip and element.tag == 'vectorhitright':
                vectorHitsRight.append(element)
                v = element
                host_qstart = int(v.find('qstart').text)
                host_qend = int(v.find('qend').text)
#                if (-ALLOWED_OVERLAP <= virus_qstart - host_qend <= config['MAX_HOST_VIRUS_DISTANCE'])) or (-ALLOWED_OVERLAP <= host_qstart - virus_qend <= config['MAX_HOST_VIRUS_DISTANCE'])):
                if (virus_qstart - host_qend <= config['MAX_HOST_VIRUS_DISTANCE'])) or (host_qstart - virus_qend <= config['MAX_HOST_VIRUS_DISTANCE'])):
                    chrSeqID = v.attrib['seqid']
                    sstart = int(v.find('sstart').text)
                    send = int(v.find('send').text)
                    if sstart < send:
                        if FIND_FEATURES:
                            site = (chrSeqID, sstart, searchiTrees(chrSeqID, sstart, sstart + config['FEATURE_SEARCH_DIST']), host_qstart - virus_qend, (virus_sstart, virus_send))
                        else:
                            site = (chrSeqID, sstart, [], host_qstart - virus_qend, (virus_sstart, virus_send))
                        insertSites[seqid][specimen][contigName]['right'].append(site)
                    else:  # this is actually a left flanking region
                        if FIND_FEATURES:
                            site = (chrSeqID, sstart, searchiTrees(chrSeqID, sstart - config['FEATURE_SEARCH_DIST'], sstart), host_qstart - virus_qend, (virus_send, virus_sstart))
                        else:
                            site = (chrSeqID, sstart, [], host_qstart - virus_qend, (virus_send, virus_sstart))
                        insertSites[seqid][specimen][contigName]['left'].append(site)
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
                if FIND_FEATURES:
                    site = (chrSeqID, (sstart, send), searchiTrees(chrSeqID, sstart, send), overlap, flanks, vse)
                else:
                    site = (chrSeqID, (sstart, send), [], overlap, flanks, vse)
                insertSites[seqid][specimen][contigName]['referenceoverlaps'].append(site)
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
                
#                 if (leftStart < leftEnd) and (rightStart < rightEnd):
#                     if rightStart < leftEnd:
#                         overlap = (rightStart, leftEnd)
#                         flankLength = leftEnd - rightStart + 1
#                         flanks = getFlanks(specimen, contigName.split('__')[0], leftQEnd - flankLength + 1, leftQEnd, rightQStart, rightQStart + flankLength - 1)
#                         flanks += (virus_qstart - leftQEnd - 1, rightQStart - virus_qend - 1)
#                     else:
#                         overlap = None
#                         flanks = ('', '', virus_qstart - leftQEnd - 1, rightQStart - virus_qend - 1)
#                 else:  # if (leftStart > leftEnd) and (rightStart > rightEnd):
#                     if rightStart > leftEnd:
#                         overlap = (leftEnd, rightStart)
#                         flankLength = rightStart - leftEnd + 1
#                         flanks = getFlanks(specimen, contigName.split('__')[0], leftQEnd - flankLength + 1, leftQEnd, rightQStart, rightQStart + flankLength - 1)
#                         flanks = (reverseComplement(flanks[1]), reverseComplement(flanks[0]), rightQStart - virus_qend - 1, virus_qstart - leftQEnd - 1)
#                     else:
#                         overlap = None
#                         flanks = ('', '', rightQStart - virus_qend - 1, virus_qstart - leftQEnd - 1)
                
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


            # elif not skip and element.tag == 'inversion':
#                 match = element
#                 for v in vectorHitsLeft:
#                     if v.attrib['id'] == match.attrib['leftid']:
#                         matchLeft = int(v.find('send').text)
#                         try:
#                             insertSites[seqid][specimen][contigName]['left'].remove((v.attrib['seqid'], matchLeft))
#                         except ValueError:
#                             pass
#                         try:
#                             insertSites[seqid][specimen][contigName]['right'].remove((v.attrib['seqid'], matchLeft))
#                         except ValueError:
#                             pass
#                         break
#                 for v in vectorHitsRight:
#                     if v.attrib['id'] == match.attrib['rightid']:
#                         matchRight = int(v.find('sstart').text)
#                         try:
#                             insertSites[seqid][specimen][contigName]['left'].remove((v.attrib['seqid'], matchRight))
#                         except ValueError:
#                             pass
#                         try:
#                             insertSites[seqid][specimen][contigName]['right'].remove((v.attrib['seqid'], matchRight))
#                         except ValueError:
#                             pass
#                         insertSites[seqid][specimen][contigName]['inversion'].append((v.attrib['seqid'], (matchLeft, matchRight)))
#                         break

    return insertSites, specimen2Label
       
def drawInsertSites(clusteredFileName, seqids, omitFirst, maxDist, unalignedFastaName, bestHitsOnly = True):
    if not Path(clusteredFileName).exists():
        writelog('File does not exist: ' + clusteredFileName, True)
        return
        
    paramOmitFirst = omitFirst
    if seqids is None:
        seqids = []
        for record in SeqIO.parse(clusteredFileName, 'fasta'):
            if omitFirst:
                omitFirst = False
                continue
            parts = record.description.split('_|_')
            if len(parts) < 4:  # no cluster ids
                seqid = parts[2].split('.')[0] 
                seqid += parts[2].split(seqid)[1][:2]
            else:
                seqid = parts[3].split('.')[0]
                seqid += parts[3].split(seqid)[1][:2]
            seqids.append(seqid)
        
        omitFirst = paramOmitFirst
        
    lengths = {'NC_035107.1': 310827022, 'NC_035108.1': 474425716, 'NC_035109.1': 409777670}
    chromNumber = {'NC_035107.1': 1, 'NC_035108.1': 2, 'NC_035109.1': 3}
    
    insertSites, specimen2Label = getInsertSites(list(set(seqids)), bestHitsOnly)
    
    y = 0
    yvalues = {} # {'NC_035107.1': [], 'NC_035108.1': [], 'NC_035109.1': []}
    xvalues = {} # {'NC_035107.1': [], 'NC_035108.1': [], 'NC_035109.1': []}
    ylabels  = []
    colors = {} # {'NC_035107.1': [], 'NC_035108.1': [], 'NC_035109.1': []}
    positions = {} # {'NC_035107.1': {}, 'NC_035108.1': {}, 'NC_035109.1': {}}
    features = {}
    pidents = {}
    totalClusterPidents = {}
    totalClusterSize = {}
    clusterCounts = {}
    clusterSeqs = {}
        
    numRecords = 0
    for record in SeqIO.parse(clusteredFileName, 'fasta'):
        if omitFirst:
            omitFirst = False
            continue
        numRecords += 1
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
        
        region, num = specimen2Label[specimen]
        ylabels.append(cluster + ' ' + region + ' ' + str(num) + ' ' + contig)
        
        if cluster not in yvalues:
            yvalues[cluster] = {'NC_035107.1': [], 'NC_035108.1': [], 'NC_035109.1': []}
            xvalues[cluster] = {'NC_035107.1': [], 'NC_035108.1': [], 'NC_035109.1': []}
            colors[cluster] = {'NC_035107.1': [], 'NC_035108.1': [], 'NC_035109.1': []}
            positions[cluster] = {'NC_035107.1': {}, 'NC_035108.1': {}, 'NC_035109.1': {}}
            features[cluster] = {'NC_035107.1': {}, 'NC_035108.1': {}, 'NC_035109.1': {}}
            pidents[cluster] = {'NC_035107.1': {}, 'NC_035108.1': {}, 'NC_035109.1': {}}
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
        
        for chromID, position, feat, dist, vse in insertSites[seqid][specimen][contig]['left']:
            if chromID in lengths:
                yvalues[cluster][chromID].append(y)
                xvalues[cluster][chromID].append(position)
                colors[cluster][chromID].append('#87a96b')    # green/asparagus
                if position not in positions[cluster][chromID]:
                    positions[cluster][chromID][position] = []
                    features[cluster][chromID][position] = []
                    pidents[cluster][chromID][position] = 0
                positions[cluster][chromID][position].append((region + '_' + str(num), contig, 'left', dist, None, vse))
                features[cluster][chromID][position].append((feat, [], []))
                pidents[cluster][chromID][position] += pident
        for chromID, position, feat, dist, vse in insertSites[seqid][specimen][contig]['right']:
            if chromID in lengths:
                yvalues[cluster][chromID].append(y)
                xvalues[cluster][chromID].append(position)
                colors[cluster][chromID].append('#ff2600')       # red
                if position not in positions[cluster][chromID]:
                    positions[cluster][chromID][position] = []
                    features[cluster][chromID][position] = []
                    pidents[cluster][chromID][position] = 0
                positions[cluster][chromID][position].append((region + '_' + str(num), contig, 'right', dist, None, vse))
                features[cluster][chromID][position].append(([], feat, []))
                pidents[cluster][chromID][position] += pident
        for chromID, position, feat, overlap, flanks, vse in insertSites[seqid][specimen][contig]['match']:
            if chromID in lengths:
#                 if position[0] > position[1]:
#                     position = position[::-1]
#                     flanks[0], flanks[1] = reverseComplement(flanks[1]), reverseComplement(flanks[0])
#                     flanks[2], flanks[3] = flanks[3], flanks[2]
                yvalues[cluster][chromID].append(y)
                xvalues[cluster][chromID].append(position[0])
                colors[cluster][chromID].append('#ff40ff')   # magenta
                if position not in positions[cluster][chromID]:
                    positions[cluster][chromID][position] = []
                    features[cluster][chromID][position] = []
                    pidents[cluster][chromID][position] = 0
                positions[cluster][chromID][position].append((region + '_' + str(num), contig, 'match', overlap, flanks, vse))
                features[cluster][chromID][position].append(feat)
                pidents[cluster][chromID][position] += pident
        for chromID, position, feat, overlap, flanks, vse in insertSites[seqid][specimen][contig]['inversion']:
            if chromID in lengths:
#                 if position[0] > position[1]:
#                     position = position[::-1]
#                     flanks[0], flanks[1] = reverseComplement(flanks[1]), reverseComplement(flanks[0])
#                     flanks[2], flanks[3] = flanks[3], flanks[2]
                yvalues[cluster][chromID].append(y)
                xvalues[cluster][chromID].append(position[0])
                colors[cluster][chromID].append('#5e5e5e')       # dark gray
                if position not in positions[cluster][chromID]:
                    positions[cluster][chromID][position] = []
                    features[cluster][chromID][position] = []
                    pidents[cluster][chromID][position] = 0
                positions[cluster][chromID][position].append((region + '_' + str(num), contig, 'inversion', overlap, flanks, vse))
                features[cluster][chromID][position].append(feat)
                pidents[cluster][chromID][position] += pident
        for chromID, position, feat, overlap, flanks, vse in insertSites[seqid][specimen][contig]['referenceoverlaps']:
            if chromID in lengths:
#                 if position[0] > position[1]:
#                     position = position[::-1]
#                     flanks[2], flanks[3] = flanks[3], flanks[2]
                yvalues[cluster][chromID].append(y)
                xvalues[cluster][chromID].append(position[0])
                colors[cluster][chromID].append('#ff40ff')   # same as match
                if position not in positions[cluster][chromID]:
                    positions[cluster][chromID][position] = []
                    features[cluster][chromID][position] = []
                    pidents[cluster][chromID][position] = 0
                positions[cluster][chromID][position].append((region + '_' + str(num), contig, 'referenceoverlap', overlap, flanks, vse))
                features[cluster][chromID][position].append(([], [], feat))
                pidents[cluster][chromID][position] += pident
        y += 1
        
    fig, ax = pyplot.subplots(1, 3, sharey = True, figsize = (7, 10), gridspec_kw = {'width_ratios': [1, lengths['NC_035108.1']/lengths['NC_035107.1'], lengths['NC_035109.1']/lengths['NC_035107.1']]})
    fig.subplots_adjust(top = 0.97, bottom = 0.03, left = 0.25, right = 0.99, wspace = 0.0)
    ax[0].set_yticks(range(numRecords))
    ax[0].set_yticklabels(ylabels)
    ax[0].set_ylim(0, numRecords)
    ax[0].invert_yaxis()
    chrCount = 0
    xtickPositions = [[0, int(1e8), int(2e8), int(3e8)], [0, int(1e8), int(2e8), int(3e8), int(4e8)], [0, int(1e8), int(2e8), int(3e8), int(4e8)]]
    xtickLabels = [['', '100M', '200M', '300M'], ['', '100M', '200M', '300M', '400M'], ['', '100M', '200M', '300M', '400M']]
    for chromID in lengths:
        for cluster in yvalues:
            ax[chrCount].scatter(xvalues[cluster][chromID], yvalues[cluster][chromID], color=colors[cluster][chromID], s = 1)
        ax[chrCount].tick_params(labelsize = 3)
        if chrCount > 0:
            ax[chrCount].tick_params(left=False)
        #pyplot.axvline(x = lengths['NC_035107.1'], linewidth=0.5, color='black')
        #pyplot.axvline(x = lengths['NC_035107.1'] + lengths['NC_035108.1'], linewidth=0.5, color='black')
        #pyplot.xticks(range(int(1e8), sum(lengths.values()), int(1e8)), [str(x)+'M' for x in range(100, lengths['NC_035107.1'] // int(1e6), 100)])
        ax[chrCount].set_xticks(xtickPositions[chrCount])
        ax[chrCount].set_xticklabels(xtickLabels[chrCount])
        ax[chrCount].set_xlim(0, lengths[chromID])
        ax[chrCount].set_xlabel('Chromosome ' + str(chrCount + 1), fontsize=4)
        chrCount += 1

#     fig = pyplot.figure(figsize=(7, 10))
#     xoffset = 0
#     for chromID in lengths:
#         for cluster in yvalues:
#             pyplot.scatter([x + xoffset for x in xvalues[cluster][chromID]], yvalues[cluster][chromID], color=colors[cluster][chromID], s = 1)
#         xoffset += lengths[chromID]
#             
#     pyplot.tick_params(labelsize = 3)
#     pyplot.subplots_adjust(top = 0.97, bottom = 0.03, left = 0.25, right = 0.99)
#     pyplot.yticks(range(len(aln) - 1), ylabels)
#     pyplot.gca().invert_yaxis()
#     pyplot.axvline(x = lengths['NC_035107.1'], linewidth=0.5, color='black')
#     pyplot.axvline(x = lengths['NC_035107.1'] + lengths['NC_035108.1'], linewidth=0.5, color='black')
#     #pyplot.xticks(range(int(1e8), sum(lengths.values()), int(1e8)), [str(x)+'M' for x in range(100, lengths['NC_035107.1'] // int(1e6), 100)])
#     xtickPositions = [0, int(1e8), int(2e8), int(3e8), lengths['NC_035107.1']+int(1e8), lengths['NC_035107.1']+int(2e8), lengths['NC_035107.1']+int(3e8), lengths['NC_035107.1']+int(4e8), lengths['NC_035107.1']+lengths['NC_035108.1']+int(1e8), lengths['NC_035107.1']+lengths['NC_035108.1']+int(2e8), lengths['NC_035107.1']+lengths['NC_035108.1']+int(3e8), lengths['NC_035107.1']+lengths['NC_035108.1']+int(4e8)]
#     xtickLabels = ['0', '100M', '200M', '300M', '100M', '200M', '300M', '400M', '100M', '200M', '300M', '400M']
#     pyplot.xticks(xtickPositions, xtickLabels)
#     pyplot.xlim(0, sum(lengths.values()))
#    fig.suptitle('EVE insertion positions for ' + seqid, fontsize = 6, fontweight = 'bold', y = 0.98)
    
    if not Path(VIRUSES_DIR).exists():
        os.system('mkdir ' + VIRUSES_DIR)
    if not Path(SEQUENCES_DIR).exists():
        os.system('mkdir ' + SEQUENCES_DIR)
    fams = readFamFile()
    familyDir = SEQUENCES_DIR + fams[seqid]
    if not Path(familyDir).exists():
        os.system('mkdir ' + familyDir)

    saveFileName = familyDir + '/insertpositions_' + seqid
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
    for chromID in lengths:
        tsvPositionFile = open(saveFileName + '_chr' + str(chromNumber[chromID]) + '.tsv', 'w')
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
#                    groups[cluster][groupLabels[cluster][groupNames]][1].append((chromNumber[chromID], position, positions[cluster][chromID][position]))
                    
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
                    groups[cluster][groupLabels[cluster][groupNames]][1].append((chromNumber[chromID], position, (upstreamRepeat, downstreamRepeat, inregionRepeats, featStrings['inregion']['intron']), positions[cluster][chromID][position]))
                    
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
    textFile = open(saveFileName + '.txt', 'w')
    
#     for cluster in clusterSeqs:
#         clusterRecords = []
#         totalLength = 0
#         for desc in clusterSeqs[cluster]:
#             clusterRecords.append(seqs[desc])
#             totalLength += len(seqs[desc].seq)
#         SeqIO.write(clusterRecords, unalignedFastaName.split('.')[0] + '_' + cluster + '.fasta', 'fasta')
    
    textFile.write('Cluster statistics\n\n')
    for cluster in clusterCounts:
        clusterRecords = []
        eveLengthFreq = {}
        for desc in clusterSeqs[cluster]:
            try:
                clusterRecords.append(seqs[desc])
                eveLength = len(seqs[desc].seq)
                if eveLength not in eveLengthFreq:
                    eveLengthFreq[eveLength] = 0
                eveLengthFreq[eveLength] += 1
            except:  # description may have been adjusted
                for desc2 in seqs:
                    if desc.split('_|_')[:2] == desc2.split('_|_')[:2]:
                        clusterRecords.append(seqs[desc2])
                        eveLength = len(seqs[desc2].seq)
                        if eveLength not in eveLengthFreq:
                            eveLengthFreq[eveLength] = 0
                        eveLengthFreq[eveLength] += 1
                        break
        SeqIO.write(clusterRecords, unalignedFastaName.split('.')[0] + '_' + cluster + '.fasta', 'fasta')
        
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
                    if (hit[2] in ['match', 'inversion', 'referenceoverlap']) and (abs(pos[0] - pos[1]) <= maxDist):
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
    print('Usage: ' + sys.argv[0] + ' [--omit_first] [--accs=XXX[,YYY,...]] [--max_dist=N] clustered_filename.fasta unaligned_sequences.fasta')

def main():
    if len(sys.argv) < 2:
        usage()
        return
        
    accs = None
    omitFirst = False
    maxDist = MAX_FLANK_DISTANCE
    for arg in sys.argv[1:-2]:
        if arg.startswith('--accs'):
            accs = arg.split('=')[1]
            accs = accs.split(',')
            if accs == ['']:
                usage()
                return
        elif arg == '--omit_first':
            omitFirst = True
        elif arg.startswith('--max_dist'):
            maxDist = arg.split('=')[1]
            try:
                maxDist = int(maxDist)
            except:
                usage()
                return
        else:
            usage()
            return
    
    if (sys.argv[-2][0] == '-') or (sys.argv[-1][0] == '-'):
        usage()
        return
       
    if FIND_FEATURES: 
        makeIntervalTrees()
    drawInsertSites(sys.argv[-2], accs, omitFirst, maxDist, sys.argv[-1])
    
if __name__ == '__main__':
    main()
    
#python3 insertsites.py '/Volumes/Data2/sequences_NC_001564.2_per_contig_aligned_clustered_k8_relabeled.fasta' 'NC_001564.2'
