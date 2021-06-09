from pathlib import Path
import matplotlib
matplotlib.use('agg')
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as pyplot
from lxml import etree as ET

from util import *
from features import *

def removeSite(sites, posPair):
    for index in range(len(sites)):
        if sites[index][:2] == posPair:
            features = sites[index][2]
            sites.pop(index)
            return features
    return []
    
def getFlanks(specimen, contigName, leftStart, leftEnd, rightStart, rightEnd):
    scaffoldsFilename = str(Path(SPECIMENS_DIR) / (specimen + '/results/scaffolds/scaffolds.fasta'))
    scaffoldRecords = SeqIO.index(scaffoldsFilename, 'fasta')
    contigRecord = scaffoldRecords[contigName]
    seq = str(contigRecord.seq)
    return seq[leftStart - 1:leftEnd], seq[rightStart - 1:rightEnd]  # account for 1-based indexing

def getInsertSites(desiredSeqID = None, bestHitsOnly = True, omitInReference = True):
    '''bestHitsOnly must match value given to getHits that created these xml files.'''
    
    xmlFiles = []
    specimensPath = Path(SPECIMENS_DIR)
    for specimenDir in specimensPath.iterdir():
        if specimenDir.name[0] != '.':
            for file in (specimenDir / 'results/xml').iterdir():
                if (file.name[0] != '.') and (('_hits_features.xml' in file.name) or (('_hits.xml' in file.name) and not Path(str(file)[:-4] + '_features.xml').exists())):
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
                skip = (desiredSeqID is not None) and (seqid != desiredSeqID)
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
            elif not skip and element.tag == 'vectorhitleft':
                vectorHitsLeft.append(element)
                v = element
                aedes_qstart = int(v.find('qstart').text)
                aedes_qend = int(v.find('qend').text)
                if (-ALLOWED_OVERLAP <= virus_qstart - aedes_qend <= MAX_FLANKING_DISTANCE) or (-ALLOWED_OVERLAP <= aedes_qstart - virus_qend <= MAX_FLANKING_DISTANCE): # or (len(aedes_qseq) >= MIN_FLANKING_HIT_LENGTH):
                    chrSeqID = v.attrib['seqid']
                    sstart = int(v.find('sstart').text)
                    send = int(v.find('send').text)
                    insertSites[seqid][specimen][contigName]['left'].append((chrSeqID, send, searchiTrees(chrSeqID, sstart, send)))
            elif not skip and element.tag == 'vectorhitright':
                vectorHitsRight.append(element)
                v = element
                aedes_qstart = int(v.find('qstart').text)
                aedes_qend = int(v.find('qend').text)
                if (-ALLOWED_OVERLAP <= virus_qstart - aedes_qend <= MAX_FLANKING_DISTANCE) or (-ALLOWED_OVERLAP <= aedes_qstart - virus_qend <= MAX_FLANKING_DISTANCE): # or (len(aedes_qseq) >= MIN_FLANKING_HIT_LENGTH):
                    chrSeqID = v.attrib['seqid']
                    sstart = int(v.find('sstart').text)
                    send = int(v.find('send').text)
                    insertSites[seqid][specimen][contigName]['right'].append((chrSeqID, sstart, searchiTrees(chrSeqID, sstart, send)))
            elif not skip and element.tag == 'vectorhitoverlap':
                v = element
                aedes_qstart = int(v.find('qstart').text)
                aedes_qend = int(v.find('qend').text)
                chrSeqID = v.attrib['seqid']
                sstart = int(v.find('sstart').text)
                send = int(v.find('send').text)
                overlap = None
                flanks = ('', '', virus_qstart - aedes_qstart - 1, aedes_qend - virus_qend - 1)  # dists beyond ends of EVE; negative means doesn't cover that end of EVE
                insertSites[seqid][specimen][contigName]['referenceoverlaps'].append((chrSeqID, (sstart, send), searchiTrees(chrSeqID, sstart, send), overlap, flanks))
            elif not skip and (element.tag == 'match' or element.tag == 'inversion'):
                match = element
                for v in vectorHitsLeft:
                    if v.attrib['id'] == match.attrib['leftid']:
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
                
                if (leftStart < leftEnd) and (rightStart < rightEnd):
                    if rightStart < leftEnd:
                        overlap = (rightStart, leftEnd)
                        flankLength = leftEnd - rightStart + 1
                        flanks = getFlanks(specimen, contigName.split('__')[0], leftQEnd - flankLength + 1, leftQEnd, rightQStart, rightQStart + flankLength - 1)
                        flanks += (virus_qstart - leftQEnd - 1, rightQStart - virus_qend - 1)
                    else:
                        overlap = None
                        flanks = ('', '', virus_qstart - leftQEnd - 1, rightQStart - virus_qend - 1)
                else:  # if (leftStart > leftEnd) and (rightStart > rightEnd):
                    if rightStart > leftEnd:
                        overlap = (leftEnd, rightStart)
                        flankLength = rightStart - leftEnd + 1
                        flanks = getFlanks(specimen, contigName.split('__')[0], leftQEnd - flankLength + 1, leftQEnd, rightQStart, rightQStart + flankLength - 1)
                        flanks = (reverseComplement(flanks[1]), reverseComplement(flanks[0]), rightQStart - virus_qend - 1, virus_qstart - leftQEnd - 1)
                    else:
                        overlap = None
                        flanks = ('', '', rightQStart - virus_qend - 1, virus_qstart - leftQEnd - 1)
                
                insertSites[seqid][specimen][contigName][element.tag].append((v.attrib['seqid'], (leftEnd, rightStart), featuresLeft + featuresRight, overlap, flanks))
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
       
def drawInsertSites(clusteredFileName, seqid, bestHitsOnly = True, omitInReference = True):
    if not Path(clusteredFileName).exists():
        writelog('File does not exist: ' + clusteredFileName, True)
        return
        
    paramSeqID = seqid
        
    lengths = {'NC_035107.1': 310827022, 'NC_035108.1': 474425716, 'NC_035109.1': 409777670}
    chromNumber = {'NC_035107.1': 1, 'NC_035108.1': 2, 'NC_035109.1': 3}
    
    insertSites, specimen2Label = getInsertSites(seqid, bestHitsOnly, omitInReference)
    
    aln = AlignIO.read(clusteredFileName, 'fasta') # Bio.Align.MultipleSeqAlignment object
    
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
        
    for i in range(1, len(aln)):
        desc = aln[i].description
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
#            ylabels[cluster] = []
#            y[cluster] = 0

        totalClusterPidents[cluster] += pident
        totalClusterSize[cluster] += 1
        
        for chromID, position, feat in insertSites[seqid][specimen][contig]['left']:
            if chromID in lengths:
                yvalues[cluster][chromID].append(y)
                xvalues[cluster][chromID].append(position)
                colors[cluster][chromID].append('#929000')
                if position not in positions[cluster][chromID]:
                    positions[cluster][chromID][position] = []
                    features[cluster][chromID][position] = []
                    pidents[cluster][chromID][position] = 0
                positions[cluster][chromID][position].append((region + '_' + str(num), contig, 'left', None, None))
                features[cluster][chromID][position].extend(feat)
                pidents[cluster][chromID][position] += pident
        for chromID, position, feat in insertSites[seqid][specimen][contig]['right']:
            if chromID in lengths:
                yvalues[cluster][chromID].append(y)
                xvalues[cluster][chromID].append(position)
                colors[cluster][chromID].append('#ff2600')
                if position not in positions[cluster][chromID]:
                    positions[cluster][chromID][position] = []
                    features[cluster][chromID][position] = []
                    pidents[cluster][chromID][position] = 0
                positions[cluster][chromID][position].append((region + '_' + str(num), contig, 'right', None, None))
                features[cluster][chromID][position].extend(feat)
                pidents[cluster][chromID][position] += pident
        for chromID, position, feat, overlap, flanks in insertSites[seqid][specimen][contig]['match']:
            if chromID in lengths:
                if position[0] > position[1]:
                    position = position[::-1]
                yvalues[cluster][chromID].append(y)
                xvalues[cluster][chromID].append(position[0])
                colors[cluster][chromID].append('#ff40ff')
                if position not in positions[cluster][chromID]:
                    positions[cluster][chromID][position] = []
                    features[cluster][chromID][position] = []
                    pidents[cluster][chromID][position] = 0
                positions[cluster][chromID][position].append((region + '_' + str(num), contig, 'match', overlap, flanks))
                features[cluster][chromID][position].extend(feat)
                pidents[cluster][chromID][position] += pident
        for chromID, position, feat, overlap, flanks in insertSites[seqid][specimen][contig]['inversion']:
            if chromID in lengths:
                if position[0] > position[1]:
                    position = position[::-1]
                yvalues[cluster][chromID].append(y)
                xvalues[cluster][chromID].append(position[0])
                colors[cluster][chromID].append('#5e5e5e')
                if position not in positions[cluster][chromID]:
                    positions[cluster][chromID][position] = []
                    features[cluster][chromID][position] = []
                    pidents[cluster][chromID][position] = 0
                positions[cluster][chromID][position].append((region + '_' + str(num), contig, 'inversion', overlap, flanks))
                features[cluster][chromID][position].extend(feat)
                pidents[cluster][chromID][position] += pident
        for chromID, position, feat, overlap, flanks in insertSites[seqid][specimen][contig]['referenceoverlaps']:
            if chromID in lengths:
                if position[0] > position[1]:
                    position = position[::-1]
                yvalues[cluster][chromID].append(y)
                xvalues[cluster][chromID].append(position[0])
                colors[cluster][chromID].append('#ff40ff')   # same as match
                if position not in positions[cluster][chromID]:
                    positions[cluster][chromID][position] = []
                    features[cluster][chromID][position] = []
                    pidents[cluster][chromID][position] = 0
                positions[cluster][chromID][position].append((region + '_' + str(num), contig, 'referenceoverlap', overlap, flanks))
                features[cluster][chromID][position].extend(feat)
                pidents[cluster][chromID][position] += pident
        y += 1
        
    fig, ax = pyplot.subplots(1, 3, sharey = True, figsize = (7, 10), gridspec_kw = {'width_ratios': [1, lengths['NC_035108.1']/lengths['NC_035107.1'], lengths['NC_035109.1']/lengths['NC_035107.1']]})
    fig.subplots_adjust(top = 0.97, bottom = 0.03, left = 0.25, right = 0.99, wspace = 0.0)
    ax[0].set_yticks(range(len(aln) - 1))
    ax[0].set_yticklabels(ylabels)
    ax[0].set_ylim(0, len(aln) - 1)
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
    if paramSeqID is None:
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
        tsvPositionFile.write('Cluster\tPosition\tRepeat\tIntron\tExon\tGroup\t% Identity\t' + '\t'.join(POPULATIONS.keys()) + '\tSpecimens\n')
        for cluster in positions:
            if cluster not in groupLabels:
                groupLabels[cluster] = {}
                nextGroupLabel[cluster] = 1
                groups[cluster] = {}
            sortedPositions = list(positions[cluster][chromID].keys())
            sortedPositions.sort(key=lambda x: int(x[0]) if isinstance(x, tuple) else int(x))
            for position in sortedPositions:
                if len(positions[cluster][chromID][position]) >= MIN_HITS_TO_SHOW_POSITION:
                    groupNames = tuple([t[0] for t in positions[cluster][chromID][position]])
                    if groupNames not in groupLabels[cluster]:
                        groupLabels[cluster][groupNames] = nextGroupLabel[cluster]
                        groups[cluster][nextGroupLabel[cluster]] = [groupNames, []]
                        nextGroupLabel[cluster] += 1
                    groups[cluster][groupLabels[cluster][groupNames]][1].append((chromNumber[chromID], position, positions[cluster][chromID][position]))
                    
                    featDict = {'repeat_region': [], 'intron': [], 'exon': []}
                    for featType, featID, featPos in features[cluster][chromID][position]:
                        if featType not in featDict:
                            featDict[featType] = []
                        featDict[featType].append(featID)
                    featStrings = {}
                    for featType in featDict:
                        if featType in ['gene', 'ncRNA_gene', 'pseudogene']:
                            if len(featDict['exon']) > 0:  # omit; just show exon
                                continue
                            else:                      # position is in an intron
                                featStrings['intron'] = ', '.join(set(featDict[featType]))
                        else:
                            featStrings[featType] = ', '.join(set(featDict[featType]))
                        
                    tsvPositionFile.write(cluster[1:] + '\t' + str(position) + '\t' + featStrings['repeat_region'] + '\t' + featStrings['intron'] + '\t' + featStrings['exon'] + '\t' + str(groupLabels[cluster][groupNames]) + '\t' + '{:5.2f}'.format(pidents[cluster][chromID][position] / len(positions[cluster][chromID][position])) + '\t')
                    
                    countsByRegion = {region: 0 for region in POPULATIONS}
                    line = ''
                    for specimen, contig, type, overlap, flanks in positions[cluster][chromID][position]:
                        region = specimen.split('_')[0]
                        countsByRegion[region] += 1
                        if isinstance(position, tuple):
                            if overlap is None:
                                line += specimen + ' (' + type + ')\t'   # include match/inversion if tuple
                            else:
                                line += specimen + ' (' + type + ', ' + str(overlap) + ', ' + str(flanks) + ')\t'
                        else:
                            line += specimen + '\t'
                    for region in countsByRegion:
                        tsvPositionFile.write(str(countsByRegion[region]) + '\t')
                    tsvPositionFile.write(line[:-1] + '\n')
        tsvPositionFile.close()
        
    textFile = open(saveFileName + '.txt', 'w')
    
    textFile.write('Insertion positions by group\n')
    for cluster in positions:
        textFile.write('\n  Cluster ' + str(cluster) + '\n')
        for groupLabel in groups[cluster]:
            textFile.write('    Group ' + str(groupLabel) + ' (' + ', '.join(groups[cluster][groupLabel][0]) + ')\n      ')
            for chr, pos, hits in groups[cluster][groupLabel][1]:
                if isinstance(pos, tuple):
                    textFile.write('{0}: {1:,}--{2:,} ({3})'.format(chr, pos[0], pos[1], hits[0][2]) + '; ')
                else:
                    textFile.write('{0}: {1:,}'.format(chr, pos) + '; ')
            textFile.write('\n')
            
    textFile.write('\nMatching flanks only by group\n')
    for cluster in positions:
        textFile.write('\n  Cluster ' + str(cluster) + '\n')
        for groupLabel in groups[cluster]:
            found2 = False
            for chr, pos, hits in groups[cluster][groupLabel][1]:
#                if isinstance(pos[1], tuple):
                found1 = False
                for hit in hits:
                    if hit[2] == 'match' or hit[2] == 'referenceoverlap':
                        if not found2:
                            textFile.write('    Group ' + str(groupLabel) + ' (' + ', '.join(groups[cluster][groupLabel][0]) + ')\n')
                            found2 = True
                        if not found1:
                            textFile.write('      {0}: {1:,}--{2:,}\n'.format(chr, pos[0], pos[1]))
                            found1 = True
                        if hit[3] is not None:
                            textFile.write('        {0:<20} {1:<10} {2:,} - {3:,} {4} ({5}, {6})\n'.format(hit[0], hit[2], hit[3][0], hit[3][1], ', '.join(hit[4][:2]), hit[4][2], hit[4][3]))
                        else:
                            textFile.write('        {0:<20} {1:<10} ({2}, {3})\n'.format(hit[0], hit[2], hit[4][2], hit[4][3]))
            if found2:
                textFile.write('\n\n')
            
    textFile.write('\nPercent identities\n')
    for cluster in positions:
        textFile.write('  Cluster ' + str(cluster) + ': {:5.2f}%\n'.format(totalClusterPidents[cluster] / totalClusterSize[cluster]))
        
    

def main():
    if len(sys.argv) < 2:
        print('Usage: ' + sys.argv[0] + 'clustered_filename.fasta [accession_number]')
        return
    if len(sys.argv) == 2:
        acc = None
    else:
        acc = sys.argv[2]
        
    makeIntervalTrees()
    drawInsertSites(sys.argv[1], acc)
    
if __name__ == '__main__':
    main()
    
#python3 insertsites.py '/Volumes/Data2/sequences_NC_001564.2_per_contig_aligned_clustered_k8_relabeled.fasta' 'NC_001564.2'
