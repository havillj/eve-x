#!/Users/havill/opt/anaconda3/bin/python3

from Bio.Graphics import GenomeDiagram
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio import GenBank
from reportlab.lib import colors
from reportlab.lib.units import cm
import sys
from pathlib import Path

INTERVAL = 1000

def addFeature(gd_diagram, index, read, arrow_color, strand):
    track = gd_diagram.get_tracks()[index]
    feature_set = track.get_sets()[0]

    feature = SeqFeature(FeatureLocation(read[0], read[1]), strand = strand)
#   arrow_color = colors.blue
    
    feature_set.add_feature(feature, color = arrow_color, name = read[2], label = True, label_position = 'middle', label_size = 6, label_angle = 0, label_strand = +1, sigil = 'ARROW', arrowshaft_height = 0.5)

def drawChimera():
    insertion = GenomeDiagram.Diagram('Insertion 2')
    virusLength = 739
    
    tick_track = insertion.new_track(2, name = 'labels', greytrack = 0, greytrack_labels = 0, scale_format = 'SInt', scale_fontangle = 45, scale_largetick_interval = 100, scale_color = colors.green, scale_fontsize = 8)
    ts = tick_track.new_set()
    
    virus_track = insertion.new_track(0, name = 'insertion', scale = False, greytrack = 0, greytrack_labels = 0)
    fs = virus_track.new_set()
    
    feature = SeqFeature(FeatureLocation(0, 738), strand = +1)
    fs.add_feature(feature, sigil = 'BIGARROW', arrowshaft_height = 1.0, color = 'green', name = 'Insertion 2', label = True, label_position = 'middle', label_size = 10, label_angle = 0, label_strand = +1)
    
    virus_tracks = [insertion.new_track(i+3, name = 'insertion', scale = False, greytrack = 0, greytrack_labels = 0) for i in range(4)]
    fss = [virus_tracks[i].new_set() for i in range(4)]
    
    feature = SeqFeature(FeatureLocation(0, 341), strand = -1)
    fss[0].add_feature(feature, sigil = 'BIGARROW', arrowshaft_height = 1.0, color = 'orange', name = '4254-3914', label = True, label_position = 'middle', label_size = 10, label_angle = 0, label_strand = +1)
    
    feature = SeqFeature(FeatureLocation(340, 396), strand = -1)
    fss[1].add_feature(feature, sigil = 'BIGARROW', arrowshaft_height = 1.0, color = 'orange', name = '  6369-6313', label = True, label_position = 'start', label_size = 10, label_angle = 0, label_strand = +1)

    feature = SeqFeature(FeatureLocation(387, 564), strand = +1)
    fss[2].add_feature(feature, sigil = 'BIGARROW', arrowshaft_height = 1.0, color = 'orange', name = '7311-7486', label = True, label_position = 'middle', label_size = 10, label_angle = 0, label_strand = +1)

    feature = SeqFeature(FeatureLocation(550, 738), strand = +1)
    fss[3].add_feature(feature, sigil = 'BIGARROW', arrowshaft_height = 1.0, color = 'orange', name = '8549-8737', label = True, label_position = 'middle', label_size = 10, label_angle = 0, label_strand = +1)

    orf_track = insertion.new_track(1, name = 'ORFs', scale = True, greytrack = 1, greytrack_labels = 5, scale_format = 'SInt', scale_fontangle = 45, scale_largetick_interval = 1000, scale_color = colors.green, scale_fontsize = 8)
    fs3 = orf_track.new_set()
    
    orfs = [(84, 269, -1), (369, 437, -1), (409, 558, 1), (465, 563, -1), (515, 700, -1)]

    for o in orfs:
        feature = SeqFeature(FeatureLocation(o[0], o[1]), strand = o[2])
        fs3.add_feature(feature, sigil = 'ARROW', arrowshaft_height = 0.25, color = 'lightgreen')

    insertion.draw(format = 'linear', fragments = 1, fragment_size = 0.4, track_size = 0.35, tracklines = False, orientation = 'landscape', pagesize = 'LETTER', start = 0, end = virusLength)
    insertion.write('insertion2.pdf', 'PDF')
    
def drawCFAV():
    virus = GenomeDiagram.Diagram('CFAV')
    virusLength = 10682
    
    tick_track = virus.new_track(0, name = 'labels', greytrack = 0, greytrack_labels = 0, scale_format = 'SInt', scale_fontangle = 45, scale_largetick_interval = 1000, scale_color = colors.green, scale_fontsize = 8)
    tick_track.new_set()
    
    virus_track = virus.new_track(1, name = 'cfav', scale = False, greytrack = 0, greytrack_labels = 0)
    fs = virus_track.new_set()
    
    parser = GenBank.FeatureParser()
    fhandle = open('cfav.gb', 'r')
    genbank_entry = parser.parse(fhandle)
    fhandle.close()
    
    for feature in genbank_entry.features:
        if feature.type == 'mat_peptide': 
            if feature.qualifiers['product'][0] == '2K protein':
                feature.label_strand = -1
            fs.add_feature(feature, name = feature.qualifiers['product'][0], label = True, label_position = 'middle', label_size = 10, label_angle = 60, sigil = 'BIGARROW', arrowshaft_height = 1.0, color = 'blue')
            
    insert_track = virus.new_track(2, name = 'Insertions', scale = False, greytrack = 0, greytrack_labels = 0)
    fs2 = insert_track.new_set()
    addFeature(virus, 2, (323, 3780, ' 1'), 'red', 1)
    
    addFeature(virus, 2, (3914, 4254, '2A'), 'green', -1)
    addFeature(virus, 2, (6313, 6369, '2B'), 'green', -1)
    addFeature(virus, 2, (7311, 7486, '2C'), 'green', +1)
    addFeature(virus, 2, (8549, 8737, '2D'), 'green', +1)
    
    orf_track = virus.new_track(3, name = 'ORFs', scale = True, greytrack = 1, greytrack_labels = 5, scale_format = 'SInt', scale_fontangle = 45, scale_largetick_interval = 1000, scale_color = colors.green, scale_fontsize = 8)
    fs3 = orf_track.new_set()
    
    orfs = [(29, 91, -1),
            (51, 203, -1),
            (53, 145, 1),
            (239, 337, 1),
            (261, 674, -1),
            (752, 1012, 1),
            (775, 960, -1),
            (804, 1067, -1),
            (873, 995, 1),
            (1194, 1364, -1),
            (1377, 1850, -1),
            (1520, 1597, -1),
            (1686, 1841, 1),
            (1867, 1995, -1),
            (1884, 1982, -1),
            (1886, 2038, 1),
            (2013, 2225, -1),
            (2048, 2182, -1),
            (2162, 2233, 1),
            (2255, 2512, 1),
            (2411, 2473, -1),
            (2448, 2543, -1),
            (2516, 2716, 1),
            (2604, 2780, -1),
            (2674, 2829, -1),
            (2819, 2890, -1),
            (2931, 3374, -1),
            (3158, 3295, -1)]

    for o in orfs:
        feature = SeqFeature(FeatureLocation(o[0]+323, o[1]+323), strand = o[2])
        fs3.add_feature(feature, sigil = 'ARROW', arrowshaft_height = 0.25, color = 'red')

#   orfs2 = [(84, 269, -1), (369, 437, -1), (409, 558, 1), (465, 563, -1), (515, 700, -1)]
    orfs2 = [(84+3914, 269+3914, -1), (409-387+7311, 558-387+7311, 1), (465-387+7311, 563-387+7311, -1)]

    for o in orfs2:
        feature = SeqFeature(FeatureLocation(o[0], o[1]), strand = o[2])
        fs3.add_feature(feature, sigil = 'ARROW', arrowshaft_height = 0.25, color = 'green')
        
    virus.draw(format = 'linear', fragments = 1, fragment_size = 0.4, track_size = 0.25, tracklines = False, orientation = 'landscape', pagesize = 'LETTER', start = 0, end = virusLength)
    virus.write('cfav.pdf', 'PDF')
    
###############################################################################
# one diagram per host
    
def drawVirus(dir, virusName, virusLength, gbFilename, hostName, features):

    virus = GenomeDiagram.Diagram(virusName)
#   virusLength = 10682
    
    tick_track = virus.new_track(0, name = 'labels', greytrack = 0, greytrack_labels = 0, scale_format = 'SInt', scale_fontangle = 45, scale_largetick_interval = 1000, scale_color = colors.green, scale_fontsize = 8)
    tick_track.new_set()
    
    virus_track = virus.new_track(1, name = virusName, scale = False, greytrack = 1, greytrack_labels = 1, greytrack_fontsize = 12, greytrack_font_color = 'white')
    fs = virus_track.new_set()
    
    parser = GenBank.FeatureParser()
    fhandle = open(gbFilename, 'r')
    genbank_entry = parser.parse(fhandle)
    fhandle.close()
    
    for feature in genbank_entry.features:
        if ((virusName == 'cfav') and (feature.type == 'mat_peptide')) or ((virusName != 'cfav') and (feature.type == 'CDS')): 
            fs.add_feature(feature, name = feature.qualifiers['product'][0], label = True, label_position = 'middle', label_size = 10, label_angle = 60, sigil = 'BIGARROW', arrowshaft_height = 1.0, color = 'darkblue')

    if features is not None:
        insert_track = virus.new_track(2, name = hostName, scale = False, greytrack = 1, greytrack_labels = 1, greytrack_fontsize = 12, greytrack_font_color = 'white')
        fs2 = insert_track.new_set()
        for feature in features:
            addFeature(virus, 2, feature[:-1], 'lightgreen', feature[-1])
    
    virus.draw(format = 'linear', fragments = 1, fragment_size = 0.4, track_size = 0.25, tracklines = False, orientation = 'landscape', pagesize = 'LETTER', start = 0, end = virusLength)
    virus.write(dir + '/' + hostName + '_' + virusName + '.pdf', 'PDF')
    
def drawResults(virusName, virusLength, gbName, resultsFileName, diagramDir, whichHits):    
    resultsFile = open(resultsFileName, 'r')
    header = resultsFile.readline()
    header = header.rstrip().split('\t')
    lengths = resultsFile.readline()
    
    for line in resultsFile:
        features = []
        label = 1
        cols = line.rstrip().split('\t')
        if cols[0][:5] == 'Total':
            break
        hostName = cols[0]
        for index in range(1, len(cols)):
            if cols[index] in whichHits:
                feature = header[index]
                feature = feature.strip('()').split(',')
                for j in range(0, len(feature), 2):
                    if int(feature[j]) < int(feature[j+1]):
                        features.append((int(feature[j]), int(feature[j+1]), str(label), 1))
                    else:
                        features.append((int(feature[j+1]), int(feature[j]), str(label), -1))
                label += 1
#       print(features)
        drawVirus(diagramDir, virusName, virusLength, gbName, hostName, features)
    resultsFile.close()
    
###############################################################################
# all hosts on one diagram

def addFeature1(gd_diagram, index, read, arrow_color, strand):
    track = gd_diagram.get_tracks()[index]
    feature_set = track.get_sets()[0]

    feature = SeqFeature(FeatureLocation(read[0], read[1]), strand = strand)
#   arrow_color = colors.blue
    
    feature_set.add_feature(feature, color = arrow_color, name = '', label = True, label_position = 'middle', label_size = 6, label_angle = 0, label_strand = +1, sigil = 'ARROW', arrowshaft_height = 0.5)

    
def drawViruses(dir, virusName, virusLength, gbFilename, hostNames, featureList):
    virus = GenomeDiagram.Diagram(virusName)
#   virusLength = 10682
    
    count = 0
    for hostName, features in zip(hostNames, featureList):
        hostName = hostName[:-9]
        if features is not None:
            insert_track = virus.new_track(count, name = hostName, scale = False, greytrack = 1, greytrack_labels = 1, greytrack_fontsize = 8, greytrack_font_color = 'white')
            fs2 = insert_track.new_set()
            for feature in features:
                if feature[2] == 'single':
                    addFeature1(virus, count, feature[:-1], 'lightgreen', feature[-1])
                else:
                    addFeature1(virus, count, feature[:-1], 'lightblue', feature[-1])
            count += 1
    
    virus_track = virus.new_track(count, name = '', scale = False, greytrack = 1, greytrack_labels = 1, greytrack_fontsize = 10, greytrack_font_color = 'white')
    fs = virus_track.new_set()
    
    parser = GenBank.FeatureParser()
    fhandle = open(gbFilename, 'r')
    genbank_entry = parser.parse(fhandle)
    fhandle.close()
    
    for feature in genbank_entry.features:
        if ((virusName == 'cfav') and (feature.type == 'mat_peptide')) or ((virusName != 'cfav') and (feature.type == 'CDS')): 
            fs.add_feature(feature, name = feature.qualifiers['product'][0].replace('protein', ''), label = True, label_position = 'start', label_size = 8, label_angle = 0, sigil = 'BIGARROW', arrowshaft_height = 0.75, color = 'cyan')

    tick_track = virus.new_track(count + 2, name = 'labels', greytrack = 0, greytrack_labels = 0, scale_format = 'SInt', scale_fontangle = 45, scale_largetick_interval = 1000, scale_color = colors.green, scale_fontsize = 8)
    tick_track.new_set()
    
    virus.draw(format = 'linear', fragments = 1, fragment_size = 1, track_size = 0.5, y = 0.1, tracklines = False, orientation = 'landscape', pagesize = 'LETTER', start = -1000, end = virusLength)
    virus.write(dir + '/' + virusName + '.pdf', 'PDF')
    
def drawResults1(virusName, virusLength, gbName, resultsFileName, diagramDir, whichHits):    
    resultsFile = open(resultsFileName, 'r')
    header = resultsFile.readline()
    header = header.rstrip().split('\t')
    lengths = resultsFile.readline()
    
    hostNames = []
    featureList = []
    
    for line in resultsFile:
        features = []
        label = 1
        cols = line.rstrip().split('\t')
        if cols[0][:5] == 'Total':
            break
        hostName = cols[0]
        for index in range(1, len(cols)):
            if cols[index] in whichHits:
                feature = header[index]
                feature = feature.strip('()').split(',')
                if len(feature) == 2:
                    type = 'single'
                else:
                    type = 'multi'
                for j in range(0, len(feature), 2):
                    if int(feature[j]) < int(feature[j+1]):
                        features.append((int(feature[j]), int(feature[j+1]), type, 1))  # no labels
                    else:
                        features.append((int(feature[j+1]), int(feature[j]), type, -1))
                label += 1
#       print(features)
        hostNames.append(hostName)
        featureList.append(features)
        
    resultsFile.close()
    drawViruses(diagramDir, virusName, virusLength, gbName, hostNames, featureList)
    
###############################################################################

def main():
#   drawCFAV()
#   drawChimera()

    VIRUS_NAMES = ['cfav', 'xincheng', 'anphevirus', 'horseradish',
                   'aus_anopheles'] #, 'liao_ning']
    VIRUS_LENGTHS = [10682, 12774, 12916, 7954, 6203]
    LIAO_NING_SEGMENT_LENGTHS = [3740, 3055, 2404, 2062, 1870, 1750, 1208, 1147, 943, 903, 897, 760]
    
    for index in range(0, len(VIRUS_NAMES)):
        virusName = VIRUS_NAMES[index]
        virusLength = VIRUS_LENGTHS[index]
        gbName = '../gb/' + virusName + '.gb'
        resultsFileName = 'HITS_' + virusName + '/results_' + virusName + '.txt'
        diagramDir = 'HITS_' + virusName

#         resultsFileName = 'HITS_' + virusName + '_combined/results_' + virusName + '.txt'
#         diagramDir = 'HITS_' + virusName + '_combined'
        whichHits = '12'  # '12' or '2'
        print(virusName)
        
        drawResults1(virusName, virusLength, gbName, resultsFileName, diagramDir, whichHits)
    
main()
