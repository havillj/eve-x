import sys
import os
from pathlib import Path
import xml.etree.ElementTree as ET
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio import SeqIO
from intervaltree import Interval, IntervalTree
import re
import pickle

def readCSV(csvFilename):
    """ Return a dictionary from BLAST hits. 
        key = contig name
        value = list of tuples, one for each hit in that contig
    """
    
    csv = open(csvFilename, 'r')
    hits = {}
    for line in csv:
        line = line.rstrip()
        row = line.split(',')
        contig = row[0]
        if contig in hits:
            hits[contig].append(tuple(row[1:]))
        else:
            hits[contig] = [tuple(row[1:])]
    csv.close()

    return hits

def getHits(dir, virusName):
    csvVirusFilename = str(dir / ('spades_' + virusName + '/blast_scaffolds_' + virusName + '.csv'))
    csvAedesFilename = str(dir / ('spades_' + virusName + '/blast_scaffolds_aa.csv'))

#   outFilenameBothEnds = str(dir / ('spades_' + virusName) / (dir.name + '_viral_hits2.txt'))
#   outFilenameOneEnd = str(dir / ('spades_' + virusName) / (dir.name + '_viral_hits1.txt'))
    
    xmlFilename = str(dir / ('spades_' + virusName) / (dir.name + '_' + virusName + '_hits.xml'))

    # os.system('blastn -query ' + str(dir) + \
    #   '/spadesCFAV/scaffolds.fasta -db aegyptidb -num_threads 8 -outfmt "10 qseqid qstart qend qseq sseqid sstart send sseq evalue bitscore" -out ' + \
    #   csvAedesFilename)
    #
    # os.system('blastn -query ' + str(dir) + \
    #   '/spadesCFAV/scaffolds.fasta -db cfavdb -num_threads 8 -task blastn -evalue 0.0001 -outfmt "10 qseqid qstart qend qseq sstart send sseq evalue bitscore" -out ' + \
    #   csvVirusFilename)

    virusHits = readCSV(csvVirusFilename)
    aedesHits = readCSV(csvAedesFilename)

    # Combine overlapping viral hits in each contig into long hits
    newVirusHits = {}
    for contig in virusHits:
        newVirusHits[contig] = []
        leftIndex = 0
        rightIndex = 0
        if len(virusHits[contig]) > 1:
            virusHits[contig].sort(key = lambda hit: int(hit[0]))
            virusHit = virusHits[contig][0]
            left = int(virusHit[0])
            right = int(virusHit[1])
            for index in range(1, len(virusHits[contig])):
                virusHit = virusHits[contig][index]
                nextLeft = int(virusHit[0])
                if nextLeft <= right + 1:
                    right =int(virusHit[1])
                    rightIndex = index
                else:
                    newVirusHits[contig].append((leftIndex, rightIndex))
                    leftIndex = index
                    rightIndex = index
                    left = int(virusHit[0])
                    right = int(virusHit[1])
        newVirusHits[contig].append((leftIndex, rightIndex))

    # Combine virus and aedes hit results to locate putative insertions.

#     outFileBothEnds = open(outFilenameBothEnds, 'w')
#     outFileOneEnd = open(outFilenameOneEnd, 'w')
#     summaryBothEnds = 'Summary of hits:\n   '
#     summaryOneEnd = 'Summary of hits:\n   '
    
    root = ET.Element('root')
    tree = ET.ElementTree(root)

    for contig in virusHits:
        if contig in aedesHits:   # only include contig if also an aa hit
            for leftIndex, rightIndex in newVirusHits[contig]:       # partition aedes hits into before and after viral hit
                virus_qstart = int(virusHits[contig][leftIndex][0])
                virus_qend = int(virusHits[contig][rightIndex][1])
                before = []
                after = []
                for index in range(len(aedesHits[contig])):
                    aedesHit = aedesHits[contig][index]
                    aedes_qstart = int(aedesHit[0])
                    aedes_qend = int(aedesHit[1])
                    aedes_chr = aedesHit[3]
                    aedes_sstart = int(aedesHit[4])
                    aedes_send = int(aedesHit[5])
                    if aedes_qend < virus_qstart + 15:
                        before.append((index, aedes_sstart, aedes_send, aedes_chr))
                    elif aedes_qstart > virus_qend - 15:
                        after.append((index, aedes_sstart, aedes_send, aedes_chr))

                # beforeHits = []      # find pairs of before/after aedes hits that are within 10Kbp of each other
#                 afterHits = []
#                 for b in before:
#                     for a in after:
#                         if (b[3] == a[3]) and (abs(b[2] - a[1]) < 10000):
#                             beforeHits.append(aedesHits[contig][b[0]])
#                             afterHits.append(aedesHits[contig][a[0]])

                contigTree = ET.SubElement(root, 'contig', {'name': contig})
                #outFileBothEnds.write('Contig: ' + contig + '\n')
                #outFileBothEnds.write('\n   Viral hit:' + '\n')
                for index in range(leftIndex, rightIndex + 1):
                    virusHit = virusHits[contig][index]
                    vElement = ET.SubElement(contigTree, 'virushit', {'seqid': virusHit[8]})
                    ET.SubElement(vElement, 'qstart').text = virusHit[0]
                    ET.SubElement(vElement, 'qend').text = virusHit[1]
                    ET.SubElement(vElement, 'qseq').text = virusHit[2]
                    ET.SubElement(vElement, 'sstart').text = virusHit[3]
                    ET.SubElement(vElement, 'send').text = virusHit[4]
                    ET.SubElement(vElement, 'sseq').text = virusHit[5]
                    ET.SubElement(vElement, 'evalue').text = virusHit[6]
                    ET.SubElement(vElement, 'bitscore').text = virusHit[7]
                    
                    #outFileBothEnds.write('      ' + ', '.join(virusHit[:2]+virusHit[3:5]+virusHit[6:]) + '\n')
#                     summaryBothEnds += virusHit[3] + '-' + virusHit[4]
#                     if index != rightIndex:
#                         summaryBothEnds += ':'
#                     else:
#                         summaryBothEnds += '\n   '
                        
                matchingFlanks = ET.SubElement(contigTree, 'flanks')
                # find pairs of before/after aedes hits that are within 10Kbp of each other
                for b in before:
                    for a in after:
                        if (b[3] == a[3]) and (abs(b[2] - a[1]) < 10000):
                            ET.SubElement(matchingFlanks, 'match', {'leftid': str(b[0]), 'rightid': str(a[0])})
                #beforeHits = list(set(beforeHits))
                #beforeHits.sort(key = lambda hit: int(hit[0]))
                #afterHits = list(set(afterHits))
                #afterHits.sort(key = lambda hit: int(hit[0]))
                #outFileBothEnds.write('\n   Aedes aegypti hits:' + '\n')

                #outFileBothEnds.write('      Left:  ')
                #spaces = 0
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
                    #outFileBothEnds.write(' ' * spaces + ', '.join(beforeHit[:2]+beforeHit[3:6]+beforeHit[7:]) + '\n')
                    #spaces = 13
                if len(after) > 0:
                    #outFileBothEnds.write('\n      Right: ')
                    #spaces = 0
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
                        #outFileBothEnds.write(' ' * spaces + ', '.join(afterHit[:2]+afterHit[3:6]+afterHit[7:]) + '\n')
                        #spaces = 13
                #outFileBothEnds.write('\n')
                
#     outFileBothEnds.write(summaryBothEnds + '\n')
#     outFileOneEnd.write(summaryOneEnd + '\n')
#     outFileBothEnds.close()
#     outFileOneEnd.close()
    
    tree.write(xmlFilename, xml_declaration=True)
    
def displayXML(fileName, outFileName, displaySeq = True):
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
            if displaySeq:
                outFile.write('      contig         {0:>5} - {1:<5} {2}\n'.format(v.find('qstart').text, v.find('qend').text, v.find('qseq').text))
                outFile.write('      {3:<14} {0:>5} - {1:<5} {2}\n\n'.format(v.find('sstart').text, v.find('send').text, v.find('sseq').text, v.attrib['seqid']))
            else:
                outFile.write('      contig         {0:>5} - {1:<5}\n'.format(v.find('qstart').text, v.find('qend').text, v.find('qseq').text))
                outFile.write('      {3:<14} {0:>5} - {1:<5}\n\n'.format(v.find('sstart').text, v.find('send').text, v.find('sseq').text, v.attrib['seqid']))

        hitsLeft = contig.findall('vectorhitleft')
        if len(hitsLeft) > 0:
            hitsLeftDict = {}  # consolidate query hit: [subject hits]
            for v in hitsLeft:
                pos = (v.find('qstart').text, v.find('qend').text, v.find('qseq').text)
                features = v.findall('feature')
                featureString = '| '
                for f in features:
                    if displaySeq:
                        featureString += f.attrib['type'] + '; ' + f.find('id').text + '; ' + f.find('location').text + ' | '
                    else:
                        featureString += f.attrib['type'] + '; ' + f.find('id').text + ' | '
                if pos not in hitsLeftDict:
                    hitsLeftDict[pos] = [(v.find('sstart').text, v.find('send').text, v.find('sseq').text, v.attrib['seqid'], featureString)]
                else:
                    hitsLeftDict[pos].append((v.find('sstart').text, v.find('send').text, v.find('sseq').text, v.attrib['seqid'], featureString))
        
            outFile.write('   Vector hits left: \n')
            posSort = list(hitsLeftDict.keys())
            posSort.sort(key = lambda pos: (int(pos[0]), int(pos[1])))
            for q in posSort:
                if displaySeq:
                    outFile.write('      contig         {0:>9} - {1:<9} {2}\n'.format(*q))
                    for s in hitsLeftDict[q]:
                        outFile.write('      {3:<14} {0:>9} - {1:<9} {2} {4}\n'.format(*s))
                else:
                    outFile.write('      contig         {0:>9} - {1:<9}\n'.format(*q))
                    for s in hitsLeftDict[q]:
                        outFile.write('      {3:<14} {0:>9} - {1:<9} {4}\n'.format(*s))
                outFile.write('\n')
        
        hitsRight = contig.findall('vectorhitright')
        if len(hitsRight) > 0:
            hitsRightDict = {}  # consolidate query hit: [subject hits]
            for v in hitsRight:
                pos = (v.find('qstart').text, v.find('qend').text, v.find('qseq').text)
                features = v.findall('feature')
                featureString = '| '
                for f in features:
                    if displaySeq:
                        featureString += f.attrib['type'] + '; ' + f.find('id').text + '; ' + f.find('location').text + ' | '
                    else:
                        featureString += f.attrib['type'] + '; ' + f.find('id').text + ' | '
                if pos not in hitsRightDict:
                    hitsRightDict[pos] = [(v.find('sstart').text, v.find('send').text, v.find('sseq').text, v.attrib['seqid'], featureString)]
                else:
                    hitsRightDict[pos].append((v.find('sstart').text, v.find('send').text, v.find('sseq').text, v.attrib['seqid'], featureString))
        
            outFile.write('   Vector hits right: \n')
            posSort = list(hitsRightDict.keys())
            posSort.sort(key = lambda pos: (int(pos[0]), int(pos[1])))
            for q in posSort:
                if displaySeq:
                    outFile.write('      contig         {0:>9} - {1:<9} {2}\n'.format(*q))
                    for s in hitsRightDict[q]:
                        outFile.write('      {3:<14} {0:>9} - {1:<9} {2} {4}\n'.format(*s))
                else:
                    outFile.write('      contig         {0:>9} - {1:<9}\n'.format(*q))
                    for s in hitsRightDict[q]:
                        outFile.write('      {3:<14} {0:>9} - {1:<9} {4}\n'.format(*s))
                outFile.write('\n')

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
                            if displaySeq:
                                featureString += f.attrib['type'] + '; ' + f.find('id').text + '; ' + f.find('location').text + ' | '
                            else:
                                featureString += f.attrib['type'] + '; ' + f.find('id').text + ' | '
                        if displaySeq:
                            outFile.write('         contig         {0:>9} - {1:<9} {2}\n'.format(v.find('qstart').text, v.find('qend').text, v.find('qseq').text))
                            outFile.write('         {3:<14} {0:>9} - {1:<9} {2} {4}\n\n'.format(v.find('sstart').text, v.find('send').text, v.find('sseq').text, v.attrib['seqid'], featureString))
                        else:
                            outFile.write('         contig         {0:>9} - {1:<9}\n'.format(v.find('qstart').text, v.find('qend').text, v.find('qseq').text))
                            outFile.write('         {3:<14} {0:>9} - {1:<9} {4}\n\n'.format(v.find('sstart').text, v.find('send').text, v.find('sseq').text, v.attrib['seqid'], featureString))
                        break
                        
                for v in contig.findall('virushit'):
                    if displaySeq:
                        outFile.write('         contig         {0:>9} - {1:<9} {2}\n'.format(v.find('qstart').text, v.find('qend').text, v.find('qseq').text))
                        outFile.write('         {3:<14} {0:>9} - {1:<9} {2}\n\n'.format(v.find('sstart').text, v.find('send').text, v.find('sseq').text, v.attrib['seqid']))
                    else:
                        outFile.write('         contig         {0:>9} - {1:<9}\n'.format(v.find('qstart').text, v.find('qend').text, v.find('qseq').text))
                        outFile.write('         {3:<14} {0:>9} - {1:<9}\n\n'.format(v.find('sstart').text, v.find('send').text, v.find('sseq').text, v.attrib['seqid']))

                for v in hitsRight:
                    if v.attrib['id'] == match.attrib['rightid']:
                        features = v.findall('feature')
                        featureString = '| '
                        for f in features:
                            if displaySeq:
                                featureString += f.attrib['type'] + '; ' + f.find('id').text + '; ' + f.find('location').text + ' | '
                            else:
                                featureString += f.attrib['type'] + '; ' + f.find('id').text + ' | '

                        if displaySeq:
                            outFile.write('         contig         {0:>9} - {1:<9} {2}\n'.format(v.find('qstart').text, v.find('qend').text, v.find('qseq').text))
                            outFile.write('         {3:<14} {0:>9} - {1:<9} {2} {4}\n\n'.format(v.find('sstart').text, v.find('send').text, v.find('sseq').text, v.attrib['seqid'], featureString))
                        else:
                            outFile.write('         contig         {0:>9} - {1:<9}\n'.format(v.find('qstart').text, v.find('qend').text, v.find('qseq').text))
                            outFile.write('         {3:<14} {0:>9} - {1:<9} {4}\n\n'.format(v.find('sstart').text, v.find('send').text, v.find('sseq').text, v.attrib['seqid'], featureString))
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

# REWRITE for XML

# def consolidate(virusName):
#     dir = Path('HITS_' + virusName).resolve()
#     viralHits1 = {}
#     viralHits2 = {}
#     allHits = []
#     for file in dir.iterdir():
#         if '_viral_hits' in file.name:
#             if 'hits1' in file.name:
#                 number = 1
#             else:
#                 number = 2
#             index = file.name.index('_viral_hits')
#             specimen = file.name[:index]
#             f = open(str(file), 'r')
#             line = f.readline()
#             while line:
#                 if line == 'Summary of hits:\n':
#                     # print(specimen, end = ', ')
#                     if number == 1:
#                         viralHits1[specimen] = []
#                     else:
#                         viralHits2[specimen] = []
#                     line = f.readline()
#                     line = line.strip()
#                     while line:
#                         if '-' in line and ':' not in line:
#                             l, r = line.split('-')
#                             if int(l) > int(r):
#                                 hit = (int(r), int(l))
#                             else:
#                                 hit = (int(l), int(r))
#                         elif ':' in line:
#                             parts = line.split('-')
#                             if int(parts[0]) > int(parts[-1]):
#                                 parts = parts[::-1]
#                                 hit = [int(parts[0])]
#                                 for index in range(1, len(parts) - 1):
#                                     r, l = parts[index].split(':')
#                                     hit.append(int(l))
#                                     hit.append(int(r))
#                                 hit.append(int(parts[-1]))
#                                 hit = tuple(hit)
#                             else:
#                                 hit = [int(parts[0])]
#                                 for index in range(1, len(parts) - 1):
#                                     r, l = parts[index].split(':')
#                                     hit.append(int(r))
#                                     hit.append(int(l))
#                                 hit.append(int(parts[-1]))
#                                 hit = tuple(hit)
#                         # print(hit, end = ', ')
#                         if number == 1:
#                             viralHits1[specimen].append(hit)
#                         else:
#                             viralHits2[specimen].append(hit)
#                         if hit not in allHits:
#                             allHits.append(hit)
#                         line = f.readline()
#                         line = line.strip()
#                     # print()
#                     break
#                 line = f.readline()
#             f.close()
# 
#     out = open(str(dir / ('results_' + virusName + '.txt')), 'w')
#     allHits.sort()                        # write (start, end) header
#     hitCounts1 = {}
#     hitCounts2 = {}
#     out.write('\t')
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
#     specimens = list(viralHits1.keys())
#     specimens.sort()
#     for specimen in specimens:
#         out.write(specimen + '\t')
#         for hit in allHits:
#             if hit in viralHits2[specimen]:
#                 out.write('2\t')
#                 hitCounts2[hit] += 1
#             elif hit in viralHits1[specimen]:
#                 out.write('1\t')
#                 hitCounts1[hit] += 1
#             else:
#                 out.write('0\t')
#         out.write('\n')
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
#     out.close()

def printFeature(f):
    print('Location: ', f.location)
    print('Operator: ', f.location_operator)
#    print('Strand: ', f.strand)
    print('Type: ' + f.type)
    print('Id: ' + f.id)
    print('Qualifiers: ' + str(f.qualifiers))
    
def populateFeature(feature, f):
    """Populate an XML feature element from a SeqFeature object."""
    
    ET.SubElement(feature, 'location').text = str(f.location)
    ET.SubElement(feature, 'type').text = f.type
    #ET.SubElement(feature, 'qualifiers').text = str(f.qualifiers)
    
    if 'gene' in f.qualifiers:
        ET.SubElement(feature, 'gene').text = f.qualifiers['gene'][0]
    else:
        ET.SubElement(feature, 'gene').text = 'None'
        
    if 'product' in f.qualifiers:
        ET.SubElement(feature, 'product').text = f.qualifiers['product'][0]
    else:
        ET.SubElement(feature, 'product').text = 'None'
        
    if 'unknown' not in f.id:
        ET.SubElement(feature, 'id').text = f.id
    elif 'protein_id' in f.qualifiers:
        ET.SubElement(feature, 'id').text = f.qualifiers['protein_id'][0]
    elif 'transcript_id' in f.qualifiers:
        ET.SubElement(feature, 'id').text = f.qualifiers['transcript_id'][0]
    elif 'gene' in f.qualifiers:
        ET.SubElement(feature, 'id').text = f.qualifiers['gene'][0]
    else:
        ET.SubElement(feature, 'id').text = 'Unknown'

def getIntrons(location):
    """Get a list of intron intervals from a CompoundLocation object representing exons."""
    
    prev = location.start - 1
    intervals = []
    for pos in sorted(location):
        if pos > prev + 1:
            intervals.append(Interval(prev+1, pos-1))
        prev = pos
    return intervals

def searchAA(iTrees, seqid, start, end, hitElement, counts):
    """Get features for a given interval and write to xml tree."""
    
    searchInterval = Interval(start, end)
    
    overlapIntervals = iTrees[seqid][start:end]
    for interval in overlapIntervals:
        f = interval.data                       # feature associated with interval
        if f.location_operator == 'join':       # mRNA with a CompoundLocation
            introns = getIntrons(f.location)
            featureType = None
            for intron in introns:
                if intron.contains_interval(searchInterval):
                    featureType = 'in_intron'                 # wholly contained in an intron
                    break
                elif intron.overlaps(searchInterval):
                    featureType = 'overlaps_intron'           # crosses an intron/exon boundary
                    break
                    
            if featureType is None:
                if start in f.location and end in f.location:
                    featureType = 'in_exon'                   # wholly contained in an exon
                else:
                    featureType = 'overlaps_end_exon'         # overlaps first or last exon
        else:
            featureType = f.type
                 
        feature = ET.SubElement(hitElement, 'feature', {'type': featureType})
        populateFeature(feature, f)
        counts[featureType] = counts.get(featureType, 0) + 1
    
    if len(overlapIntervals) == 0:
        counts['None'] += 1
            
def findFeatures(fileName, iTrees):
    """Find all features in an xml file and insert them."""
    
    tree = ET.parse(fileName)
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
        
    tree.write(fileName[:-4] + '_features.xml', xml_declaration=True)
    
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
        if seqid in chromNumbers:
            seqid = chromNumbers[seqid]
        else:
            seqid += '.1'
        type = cols[2]
        if type in ['chromosome', 'supercontig']:
            continue
        start = int(cols[3])
        end = int(cols[4])
        strand = strands[cols[6]]
        attrib = cols[8]
        
        f = SeqFeature()
        f.location = FeatureLocation(start, end + 1, strand = strand)
        f.type = type
        findName = re.findall(r'Name=(.*);.*', 'Name=Gypsy-58_AA-I;class=LTR/Gypsy')
        if len(findName) > 0:
            f.id = findName[0]
        else:
            f.id = 'Unknown'
        if seqid in iTrees:
            iTrees[seqid][start:end] = f
    
    gff.close()

def getGenBankFeatures(fileName, iTrees):
    """Read features from a multi-record flat GenBank file and insert into dictionary of interval trees."""

    for seqrecord in SeqIO.parse(fileName, 'genbank'):
        iTrees[seqrecord.id] = IntervalTree()
        altID = re.findall(r'NIGP01.*\.1', seqrecord.annotations['comment'])
        if len(altID) > 0:
            altID = altID[0]
            iTrees[altID] = iTrees[seqrecord.id]
        for f in seqrecord.features:
            if f.type != 'source':
                iTrees[seqrecord.id][f.location.start:f.location.end - 1] = f
                
def makeIntervalTrees():
    if Path('featuretrees.pickle').exists():
        print('Pickled feature trees found.  Reading...')
        picklediTreesDict = open('featuretrees.pickle', 'rb')
        iTrees = pickle.load(picklediTreesDict)
        picklediTreesDict.close()
    else:
        iTrees = {}
        print('Reading GenBank features...')
        # https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/002/204/515/GCF_002204515.2_AaegL5.0/
        getGenBankFeatures('/home/havill/data/aegypti/gb/GCF_002204515.2_AaegL5.0_genomic.gbff', iTrees)
    
        print('Reading repeat features...')
        getGFFFeatures('/home/havill/data/aegypti/gb/Aedes-aegypti-LVP_AGWG_REPEATFEATURES_AaegL5.gff3', iTrees)
    
        print('Writing feature trees to disk...')
        picklediTreesDict = open('featuretrees.pickle', 'wb')
        pickle.dump(iTrees, picklediTreesDict)
        picklediTreesDict.close()
        
    return iTrees

def process(dir, virusName, iTrees, displaySeqs):
    print('Processing ' + dir + '...')
    print('Getting viral hits...')
    getHits(Path(dir), virusName)
    xmlResultFile = dir + '/spades_' + virusName + '/' + dir + '_' + virusName + '_hits.xml'
#    displayXML(xmlResultFile, xmlResultFile[:-4] + '.txt')
    
    print('Annotating features...')
    findFeatures(xmlResultFile, iTrees)
    
    if displaySeqs:
        displayXML(xmlResultFile[:-4] + '_features.xml', xmlResultFile[:-4] + '_features.txt', True)
    else:
        displayXML(xmlResultFile[:-4] + '_features.xml', xmlResultFile[:-4] + '_features_short.txt', False)

def main():
    combinedDirs = ['ALL_combined', 'Amacuzac_Mexico_combined', 'Bangkok_Thailand_combined', 'Cairns_Australia_combined', 'Cebu_City_Philippines_combined', 'Cuanda_Angola_combined', 'El_Dorado_Argentina_combined', 'HoChiMin_Vietnam_combined', 'La_Lope-Gabon_combined', 'Maricopa_County_Arizona-USA_combined', 'Paqueta-Brazil_combined', 'Skukusa_South_Africa_combined', 'Tahiti_FrenchPolynesia_combined']
    virusName = 'cfav'
    
    iTrees = makeIntervalTrees()
    displaySeqs = False
    
    for dir in combinedDirs:
        process(dir, virusName, iTrees, displaySeqs)

# generating tree: 8 min; + write tree: 10 min
# reading pickle: 4:45 min


# def main():
#     if (len(sys.argv) == 3) and (sys.argv[1] == '-c'):
#         virusName = sys.argv[2]
#         if not Path('HITS_' + virusName).exists():
#             os.mkdir('HITS_' + virusName)
#         os.system('cp -v */spades_' + virusName + '/*viral_hits* ' + 'HITS_' + virusName)
#         consolidate(virusName)
#     elif len(sys.argv) == 3:
#         dir = Path(sys.argv[1])
#         dir = dir.resolve()
#         virusName = sys.argv[2]
#         getHits(dir, virusName)
#     elif len(sys.argv) == 2:
#         dir = Path('.').resolve()
#         virusName = sys.argv[1]
#         for subdir in dir.iterdir():
#             if subdir.is_dir() and ('test' not in str(subdir)) and (subdir / 'spades_' + virusName / 'scaffolds.fasta').exists():
#                 print(subdir)
#                 getHits(subdir, virusName)

#if __name__ == '__main__':
#    main()

dir = 'Bangkok_Thailand_combined'
virusName = 'cfav'
xmlResultFile = dir + '/spades_' + virusName + '/' + dir + '_' + virusName + '_hits.xml'
displayXML(xmlResultFile[:-4] + '_features.xml', xmlResultFile[:-4] + '_features_short.txt', False)
