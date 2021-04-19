from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio import SeqIO, Entrez
from dna_features_viewer import GraphicFeature, GraphicRecord, BiopythonTranslator
import matplotlib.pyplot as pyplot
import matplotlib
from pathlib import Path
import re
from lxml import etree as ET
import os
import numpy as np
import math
import sys

RESULTS_DIR = '/home/havill/data2/results4/'
GB_DIR = '/home/havill/data/aegypti/gb/'

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

def drawFromFile(fileName, virusName):
    f = open(fileName, 'r')
    features = []
    for line in f:
        cols = line.rstrip().split('\t')
        
        features.append(GraphicFeature(start=int(cols[0]), end=int(cols[1]), strand=0, color='gray', label = None, thickness = 10, linewidth = 0))
    f.close()
    
    width = 10
    height = 2
    thickness = 10
    fontsize = 8
    
    virusRecord = MyCustomTranslator().translate_record('/home/havill/data/aegypti/gb/' + virusName + '.gb')
    
    fig, ax = pyplot.subplots(2, 1, sharex = True, sharey = False, figsize = (width, height), gridspec_kw = {'height_ratios': [4, 1]})
    virusRecord.plot(ax = ax[0], figure_width=width, strand_in_label_threshold=7, with_ruler = False)
#    fig.subplots_adjust(top = 0.99, bottom = 0.01, left = 0.18, right = 0.95)
    
    
    for f in virusRecord.features:
        f.thickness = thickness
        f.linewidth = 0.5
        f.fontdict = {'size': fontsize}
                
    record = GraphicRecord(features = features, sequence_length = virusRecord.sequence_length, feature_level_height = 0)
    record.plot(ax = ax[1], figure_width=width, with_ruler = False)
    fig.savefig(virusName + '.pdf')
    
#drawFromFile('flav.tsv', 'xishuangbanna')

###############################################################################

def drawEVEs(EVEs, openEnds, virusName, title):

    width = 10
    height = 3
    thickness = 10
    fontsize = 8
    
    virusRecord = MyCustomTranslator().translate_record('/home/havill/data/aegypti/gb/' + virusName + '.gb')
    for f in virusRecord.features:
        f.thickness = thickness
        f.linewidth = 0.5
        f.fontdict = {'size': fontsize}
    
    fig, ax = pyplot.subplots(len(EVEs) + 1, 1, sharex = False, sharey = False, figsize = (width, height), gridspec_kw = {'height_ratios': [5] + [1] * len(EVEs)})
    virusRecord.plot(ax = ax[0], figure_width=width, with_ruler = True, max_line_length = 80, max_label_length = 80)
    x0, y0, w, h = ax[0].get_position().bounds
    ax[0].set_position([x0, y0+h*0.1, w, h])
    ax[0].tick_params(labelsize = fontsize)
#    fig.suptitle(title, fontsize = fontsize, fontweight = 'bold', y = 0.9)
#    fig.subplots_adjust(top = 0.99, bottom = 0.01, left = 0.18, right = 0.95)
        
    eveNum = 1
    for eve, open in zip(EVEs, openEnds):
        eveName = eve[0]
        features = []
        for index in range(1, len(eve), 2):
            features.append(GraphicFeature(start=eve[index], end=eve[index+1], strand=1, color='blue', label = None, thickness = 10, linewidth = 0, open_left = open[index - 1], open_right = open[index]))
           
        record = GraphicRecord(features = features, sequence_length = virusRecord.sequence_length, feature_level_height = 0)
        record.plot(ax = ax[eveNum], figure_width=width, with_ruler = False)
        ax[eveNum].text(0, -0.3, eveName + ' ', horizontalalignment = 'right', fontdict = {'size': fontsize})
        eveNum += 1
        
    fig.savefig(title + '.pdf')
    
# CFAV_EVEs = [('CFAV EVE 1', 1816, 4663, 739, 1094, 4953, 6564),
#              ('CFAV EVE 2', 1925, 2498, 4559, 5270),
#              ('CFAV EVE 3', 4254, 3914, 6369, 6313, 7311, 7486, 8549, 8737),
#              ('CFAV EVE 4', 3585, 3915),
#              ('CFAV EVE 5', 323, 3780)]
#              
# drawEVEs(CFAV_EVEs, 'NC_001564.2', 'CFAV EVEs')

# XIN_EVEs = [('XIN EVE 1', 10, 1401),
#              ('XIN EVE 2', 3775, 5551),
#              ('XIN EVE 3', 4286, 5551)]
# openEnds = [(False, False), (False, True), (True, True)]
#              
# drawEVEs(XIN_EVEs, openEnds, 'MH037149.1', 'XIN EVEs')

# ORTH_EVEs = [('ORTH EVE 1', 148, 824),
#              ('ORTH EVE 2', 455, 966)]
# openEnds = [(False, False), (False, False)]
#              
# drawEVEs(ORTH_EVEs, openEnds, 'MF176251.1', 'ORTH EVEs')
# exit()

###############################################################################

# code from https://matplotlib.org/stable/gallery/images_contours_and_fields/image_annotated_heatmap.html
    
def heatmap(data, row_labels, col_labels, ax=None, cbar_kw={}, cbarlabel="", **kwargs):
    """
    Create a heatmap from a numpy array and two lists of labels.

    Parameters
    ----------
    data
        A 2D numpy array of shape (N, M).
    row_labels
        A list or array of length N with the labels for the rows.
    col_labels
        A list or array of length M with the labels for the columns.
    ax
        A `matplotlib.axes.Axes` instance to which the heatmap is plotted.  If
        not provided, use current axes or create a new one.  Optional.
    cbar_kw
        A dictionary with arguments to `matplotlib.Figure.colorbar`.  Optional.
    cbarlabel
        The label for the colorbar.  Optional.
    **kwargs
        All other arguments are forwarded to `imshow`.
    """

    if not ax:
        ax = pyplot.gca()

    # Plot the heatmap
    im = ax.imshow(data, **kwargs)

    # Create colorbar
    cbar = ax.figure.colorbar(im, ax=ax, **cbar_kw)
    cbar.ax.set_ylabel(cbarlabel, rotation=-90, va="bottom")

    # We want to show all ticks...
    ax.set_xticks(np.arange(data.shape[1]))
    ax.set_yticks(np.arange(data.shape[0]))
    # ... and label them with the respective list entries.
    ax.set_xticklabels(col_labels, fontfamily='serif', fontweight='bold')
    ax.set_yticklabels(row_labels, fontfamily='serif', fontweight='bold')

    # Let the horizontal axes labeling appear on top.
    ax.tick_params(top=False, bottom=False, left=False,
                   labeltop=True, labelbottom=False)

    # Rotate the tick labels and set their alignment.
    pyplot.setp(ax.get_xticklabels(), rotation=0, ha="center",
             rotation_mode="anchor")

    # Turn spines off and create white grid.
    for edge, spine in ax.spines.items():
        spine.set_visible(False)

    ax.set_xticks(np.arange(data.shape[1]+1)-.5, minor=True)
    ax.set_yticks(np.arange(data.shape[0]+1)-.5, minor=True)
    ax.grid(which="minor", color="w", linestyle='-', linewidth=3)
    ax.tick_params(which="minor", bottom=False, left=False)

    return im, cbar

def annotate_heatmap(im, data=None, valfmt="{x:.2f}", textcolors=("black", "white"), threshold=None, **textkw):
    """
    A function to annotate a heatmap.

    Parameters
    ----------
    im
        The AxesImage to be labeled.
    data
        Data used to annotate.  If None, the image's data is used.  Optional.
    valfmt
        The format of the annotations inside the heatmap.  This should either
        use the string format method, e.g. "$ {x:.2f}", or be a
        `matplotlib.ticker.Formatter`.  Optional.
    textcolors
        A pair of colors.  The first is used for values below a threshold,
        the second for those above.  Optional.
    threshold
        Value in data units according to which the colors from textcolors are
        applied.  If None (the default) uses the middle of the colormap as
        separation.  Optional.
    **kwargs
        All other arguments are forwarded to each call to `text` used to create
        the text labels.
    """

    if not isinstance(data, (list, np.ndarray)):
        data = im.get_array()

    # Normalize the threshold to the images color range.
    if threshold is not None:
        threshold = im.norm(threshold)
    else:
        threshold = im.norm(data.max())/2.

    # Set default alignment to center, but allow it to be
    # overwritten by textkw.
    kw = dict(horizontalalignment="center",
              verticalalignment="center")
    kw.update(textkw)

    # Get the formatter in case a string is supplied
    if isinstance(valfmt, str):
        valfmt = matplotlib.ticker.StrMethodFormatter(valfmt)

    # Loop over the data and create a `Text` for each "pixel".
    # Change the text's color depending on the data.
    texts = []
    for i in range(data.shape[0]):
        for j in range(data.shape[1]):
            kw.update(color=textcolors[int(im.norm(data[i, j]) > threshold)])
            if data[i,j] > 0.0:
                text = im.axes.text(j, i, valfmt(data[i, j], None), **kw)
                texts.append(text)

    return texts
    
def makeHeatMap(values, xlabels, ylabels):
    fig, ax = pyplot.subplots()

    im, cbar = heatmap(values, ylabels, xlabels, ax=ax, cmap="Blues", cbarlabel="")
    texts = annotate_heatmap(im, valfmt="{x:.2f}", fontsize=8)
    
    fig.savefig('heatmap.pdf')
    
def readPopCSV(fileName):
    csvFile = open(fileName, 'r')
    header = csvFile.readline()
    columns = header.rstrip().split(',')
    dataCols = []
    for index in range(len(columns)):
        if 'EVE' in columns[index] and 'frequency' in columns[index]:
            dataCols.append(index)
            
    regions = []
    data = []
    for line in csvFile:
        columns = line.rstrip().split(',')
        if columns[0] != '':
            regions.append(columns[0])
            row = []
            for index in dataCols:
                row.append(float(columns[index]))
            data.append(row)
            
    csvFile.close()
    
    return np.array(data), regions
    
# data, regions = readPopCSV('/home/havill/data2/ortho freq.csv')
# makeHeatMap(data, ['1', '2'], regions)
# exit()

# eves_cfav = np.array([[0.00, 0.35, 0.12, 0.00, 0.00],
# [0.00, 0.52, 0.04, 0.00, 0.00],
# [0.00, 0.50, 0.00, 0.00, 0.00],
# [0.00, 0.00, 0.59, 0.00, 0.00],
# [0.00, 0.00, 0.09, 0.18, 0.00],
# [0.00, 0.00, 0.95, 0.00, 0.00],
# [0.55, 0.00, 0.45, 0.05, 0.00],
# [0.91, 0.00, 0.00, 0.00, 0.00],
# [0.00, 0.00, 0.20, 0.00, 0.00],
# [0.00, 0.00, 0.94, 0.00, 0.00],
# [0.00, 0.00, 1.00, 0.00, 0.65],
# [0.00, 0.00, 0.80, 0.00, 0.90]])

#ylabels = ['Angola', 'Gabon', 'South Africa', 'Mexico', 'USA', 'Australia', 'French Polynesia', 'Philippines', 'Argentina', 'Brazil', 'Thailand', 'Vietnam']
#ylabels = [r'\textbf{' + s + r'}' for s in ylabels]
#xlabels = ['EVE 1', 'EVE 2', 'EVE 3', 'EVE 4', 'EVE 5']
#xlabels = ['1', '2', '3', '4', '5']
#xlabels = [r'\textbf{' + s + r'}' for s in xlabels]

# eves_xin = np.array([[0.12, 0.24, 0.24, 0.24, 0.00, 0.71],
# [0.13, 0.00, 0.13, 0.04, 0.13, 0.70],
# [0.50, 0.21, 0.00, 0.00, 0.00, 0.21],
# [0.55, 0.00, 0.00, 0.05, 0.00, 0.95],
# [0.55, 0.09, 0.05, 0.05, 0.05, 0.91],
# [0.00, 0.00, 0.21, 0.00, 0.00, 1.00],
# [0.00, 0.00, 0.20, 0.25, 0.10, 0.65],
# [0.00, 0.00, 0.00, 0.00, 0.00, 0.96],
# [0.20, 0.00, 0.00, 0.20, 0.20, 0.60],
# [0.06, 0.00, 0.29, 0.00, 0.29, 0.71],
# [0.00, 0.00, 0.00, 0.00, 0.00, 1.00],
# [0.00, 0.00, 0.00, 0.00, 0.00, 0.95]])
# 
# xlabels = ['1a', '1b', '2a', '2b', '2c', '2d']
# xlabels = [r'\textbf{' + s + r'}' for s in xlabels]
# 
# eves_xin = np.array([[0.35, 1.00],
# [0.13, 1.00],
# [0.71, 0.21],
# [0.55, 1.00],
# [0.64, 1.00],
# [0.00, 1.00],
# [0.00, 1.00],
# [0.00, 0.96],
# [0.20, 1.00],
# [0.06, 1.00],
# [0.00, 1.00],
# [0.00, 0.95]])

# eves_xin = np.array([[0.35, 0.24, 0.94],
# [0.13, 0.13, 0.87],
# [0.71, 0.00, 0.21],
# [0.55, 0.00, 0.95],
# [0.59, 0.05, 1.00],
# [0.00, 0.21, 1.00],
# [0.00, 0.20, 1.00],
# [0.00, 0.00, 0.96],
# [0.00, 0.80, 1.00],
# [0.06, 0.29, 1.00],
# [0.00, 0.00, 1.00],
# [0.00, 0.00, 0.90]])
# # 
# xlabels = ['1', '2', '3']
# xlabels = [r'\textbf{' + s + r'}' for s in xlabels]
# # 
# pyplot.rcParams.update({
#     "text.usetex": True,
#     "font.family": "serif",
#     "font.serif": ["Computer Modern Roman"],
#     "font.sans-serif": ["Computer Modern Sans serif"]})
# #     
# # ##Needed to install additional latex packages under ubuntu:
# # ##sudo apt-get install dvipng texlive-latex-extra texlive-fonts-recommended cm-super
# # 
# makeHeatMap(eves_xin, xlabels, ylabels)
# exit()

###############################################################################

def getHits1(virusACC):
    dir = Path(RESULTS_DIR + 'specimens/')
    subdirs = [d for d in dir.iterdir()]
    files = [d / ('xml/' + d.name + '_all_hits_features.xml') for d in subdirs ]
    files.sort()

    virusHits = {}
    allSpecimens = []
    for file in files:
        specimen = file.name.rstrip('_all_hits_features.xml')
        print(specimen)
        allSpecimens.append(specimen)
        tree = ET.parse(str(file))
        root = tree.getroot()
        for contig in root:
            virushits = contig.findall('virushit')
            v = virushits[0]  # first virus hit
            seqid = v.attrib['seqid']
            if '|' in seqid:
                seqid = seqid.split('|')[1]
            if seqid == virusACC:
                if specimen not in virusHits:
                    virusHits[specimen] = Hits()

                hit = []
                for virushit in virushits:
                    start = int(virushit.find('sstart').text)
                    end = int(virushit.find('send').text)
                    hit.extend([start, end])
                if 'inreference' in contig.attrib:
                    overlap = contig.attrib['inreference'] == 'True'
                else:
                    overlap = False
                primary = contig.attrib['besthit'] == 'True'  # and (len(virushits) < 5):
                virusHits[specimen].addHit(hit, overlap, primary)
    return virusHits, allSpecimens
    
def includesInterval(virusACC, interval):
    virusHits, allSpecimens = getHits1(virusACC)
    
    intervalHits = {}
    for specimen in allSpecimens:
        intervalHits[specimen] = (specimen in virusHits) and (interval in virusHits[specimen])
    
    return intervalHits
    
# ih = includesInterval('NC_001564.2', 4100)
# keys = list(ih.keys())
# keys.sort()
# for key in keys:
#     print(key, ih[key])
# exit()

from pipeline_utilities import readFamFile, getFamily, addFamily

# Write all virus sequences for ACCs for family
# Remove gaps
def getVirusSequences(family, virusACCs, outFileName):
    if family is not None:
        fams = readFamFile()
            
    dir = Path(RESULTS_DIR + 'specimens/')
    subdirs = [d for d in dir.iterdir()]
    subdirs.sort()
    
    outFile = open(outFileName, 'w')
    for subdir in subdirs:
        fileName = str(subdir / ('sequences/' + subdir.name + '_all_hits2.fasta'))  # combined in contigs
        f = open(fileName, 'r')
        headerHost = f.readline()
        while headerHost != '':
            seqHost = f.readline()
#            headerReference = f.readline()
#            seqReference = f.readline()
            parts = headerHost[1:].split('_')
            if len(parts[0]) < 4:
                virusACC = parts[0] + '_' + parts[1]
            else:
                virusACC = parts[0]
            doWrite = False
            if family is not None:
                if virusACC not in fams:
                    fam = getFamily(virusACC)
                    fams[virusACC] = fam
                    addFamily(virusACC, fam)
                if fams[virusACC] == family:
                    doWrite = True
            elif virusACC in virusACCs:
                doWrite = True
            if doWrite:
#                outFile.write('>' + subdir.name + '_' + headerHost[1:].rstrip(' (host)\n').replace(' ', '_') + '\n')
                outFile.write('>' + subdir.name + '_' + headerHost[1:].replace(' ', '-'))
                outFile.write(seqHost.replace('-', ''))  # remove gaps
            headerHost = f.readline()
        f.close()
    outFile.close()
            
# getVirusSequences(None, ['NC_001564.2'], 'cfav.fasta')
# getVirusSequences('Orthomyxoviridae', None, 'orthomyxoviridae.fasta')
# exit()

# from pipeline, move to utilities
def reverseComplement(dna):
	'''Return the reverse complement of dna.'''
	
	dna = dna.upper()
	basecomplement = {'A': 'T', 'C': 'G', 'T': 'A', 'G': 'C'} 
	complement = ''
	for base in reversed(dna):
		complement = complement + basecomplement.get(base, base)
	return complement

# only works for mosaics if parts appear in same order in both contig and virus
# for everything to align properly, can also do this for one virus
def getVirusSequences2(virusACC, vStart, vEnd, flankLength, outFileName, sort = 'left'):
    dir = Path(RESULTS_DIR + 'specimens/')
    subdirs = [d for d in dir.iterdir()]
    subdirs.sort()
    
    refIndices = {}    # specimen seq at each reference index
    info = {}

    outFile = open(outFileName + '.fasta', 'w')
    outFile2 = open(outFileName + '.csv', 'w')
    for subdir in subdirs:
        fileName = str(subdir / ('xml/' + subdir.name + '_all_hits_features.xml'))
        specimen = str(subdir.name)
        scaffoldsFileName = str(subdir / ('scaffolds/scaffolds.fasta'))
        scaffoldRecords = SeqIO.index(scaffoldsFileName, 'fasta')
        print(specimen)
        tree = ET.parse(str(fileName))
        root = tree.getroot()
        for contig in root:
            virushits = contig.findall('virushit')
            v = virushits[0]  # first virus hit
            seqid = v.attrib['seqid']
            if '|' in seqid:
                seqid = seqid.split('|')[1]
            if seqid == virusACC:
                contigName = contig.attrib['name'].split('__')[0]
                contigSeq = str(scaffoldRecords[contigName].seq)
                refIndices[(specimen, contigName)] = {}
                
                # minSS = int(v.find('sstart').text)
#                 maxSE = int(v.find('send').text)
#                 minQS = int(v.find('qstart').text)
#                 maxQE = int(v.find('qend').text)
#                 if minSS > maxSE:
#                     minSS, maxSE = maxSE, minSS
#                     minQS, maxQE = len(contigSeq) - maxQE + 1, len(contigSeq) - minQS + 1
#                     contigSeq = reverseComplement(contigSeq)
                minSS = math.inf
                maxSE = 0
                minQS = math.inf
                maxQE = 0
                reverse = False
                for virushit in virushits:
                    qstart = int(virushit.find('qstart').text)
                    qend = int(virushit.find('qend').text)
                    sstart = int(virushit.find('sstart').text)
                    send = int(virushit.find('send').text)
                    sseq = virushit.find('sseq').text
                    qseq = virushit.find('qseq').text
                    if sstart > send:
                        sstart, send = send, sstart
                        qstart, qend = len(contigSeq) - qend + 1, len(contigSeq) - qstart + 1
                        sseq = reverseComplement(sseq)
                        qseq = reverseComplement(qseq)
                        reverse = True
                    minSS = min(minSS, sstart)
                    maxSE = max(maxSE, send)
                    minQS = min(minQS, qstart)
                    maxQE = max(maxQE, qend)
                    
                    refIndex = sstart - 2  # -1 for 0-indexing, -1 for initial increment
                    for index in range(len(sseq)):
                        if sseq[index] != '-':
                            refIndex += 1
                            refIndices[(specimen, contigName)][refIndex] = qseq[index]
                        else:
                            refIndices[(specimen, contigName)][refIndex] += qseq[index]
                if reverse:
                    contigSeq = reverseComplement(contigSeq)
                info[(specimen, contigName)] = (minSS, maxSE, minQS, maxQE, contigSeq)
#                 print(contigName, minSS, maxSE, minQS, maxQE, contigSeq)
#                 print(refIndices[(specimen, contigName)])

    referenceFilename = '/home/havill/data/aegypti/genomes/' + virusACC + '.fasta'
#     try:
#        refFile = open(referenceFilename, 'r')
#     except:
#         prefix1 = 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi'
#         prefix2 = '?db=nuccore&id='
#         suffix = '&rettype=fasta&retmode=text'
#         url = prefix1 + prefix2 + seqid + suffix
#         try:
#             fastaFile = web.urlopen(url)
#         except:
#             print('Error getting reference sequence ' + seqid + ' from NCBI.  Aligned FASTA file not created.')
#             continue
#         fastaText = fastaFile.read().decode('utf-8')
#         fastaFile.close()
#         refFile = open(referenceFilename, 'w')
#         refFile.write(fastaText)
#         refFile.close()
#         refFile = open(referenceFilename, 'r')
        
    refRecord = SeqIO.read(referenceFilename, 'fasta')
#    refFile.close()
    referenceSeq = str(refRecord.seq)
    virusLength = len(referenceSeq)
    referenceID = refRecord.id
    assert referenceID == virusACC
    referenceSeqAligned = ''
    
    if vEnd is None:
        vEnd = virusLength - 1
    
    for refIndex in range(virusLength):
        maxLength = max([len(refIndices[(specimen, contigName)].get(refIndex, '-')) for (specimen, contigName) in refIndices])
        referenceSeqAligned += '{0:-<{1}}'.format(referenceSeq[refIndex], maxLength)
                
    viralSeqs = {}
    for (specimen, contigName) in refIndices:
        viralSeqs[(specimen, contigName)] = ''
    alignedRefLength = 0
    for refIndex in range(vStart - 1, vEnd):
        maxLength = max([len(refIndices[(specimen, contigName)].get(refIndex, '-')) for (specimen, contigName) in refIndices])
        alignedRefLength += maxLength
        for (specimen, contigName) in refIndices:
            viralSeqs[(specimen, contigName)] += '{0:-<{1}}'.format(refIndices[(specimen, contigName)].get(refIndex, ''), maxLength)
        
    refStart = vStart - flankLength
    refEnd = vStart + alignedRefLength + flankLength - 1
    addGapsLeft = addGapsRight = 0
    if refStart < 1:
        refStart = 1
        addGapsLeft = 1 - refStart
    if refEnd > len(referenceSeqAligned):
        refEnd = len(referenceSeqAligned)
        addGapsRight = refEnd - len(referenceSeqAligned)
    outFile.write('>' + virusACC + '_' + str(max(1, vStart - flankLength)) + '-' + str(min(virusLength, vEnd + flankLength)) + '\n')
    outFile.write('-' * addGapsLeft + referenceSeqAligned[refStart-1:refEnd] + '-' * addGapsRight + '\n')
        
    fastaDict = {}
    csvDict = {}
    for (specimen, contigName) in viralSeqs:
        minSS, maxSE, minQS, maxQE, contigSeq = info[(specimen, contigName)]
        viralSeq = viralSeqs[(specimen, contigName)]
        
        index = 0
        while viralSeq[index] == '-':
            index += 1
        countGapsLeft = index
        
        index = len(viralSeq) - 1
        while viralSeq[index] == '-':
            index -= 1
        countGapsRight = len(viralSeq) - 1 - index
        
        viralSeq = viralSeq.strip('-')
            
        sliceStart = minQS - (countGapsLeft + flankLength) # additional virus + left flank
        sliceEnd = maxQE + (countGapsRight + flankLength)
        addGapsLeft = addGapsRight = 0
        if sliceStart < 1:
            addGapsLeft = 1 - sliceStart
            sliceStart = 1
        if sliceEnd > len(contigSeq):
            addGapsRight = sliceEnd - len(contigSeq)
            sliceEnd = len(contigSeq)
        
        header = specimen + '|' + contigName + '_(' + str(sliceStart) + '-' + str(sliceEnd) + ')'
        fastaDict[header] = '-' * addGapsLeft + contigSeq[sliceStart-1:minQS-1] + viralSeq + contigSeq[maxQE:sliceEnd] + '-' * addGapsRight
        csvDict[header] = header + ',' + str(addGapsLeft + (minQS - 2) - (sliceStart - 1) + 1) + ',' + str(len(viralSeq))
        
    items = list(fastaDict.items())
    if sort.lower() == 'left':
        items.sort(key = lambda t: t[1], reverse = True)  # reverse causes gaps to be at the end
    elif sort.lower() == 'right':
        items.sort(key = lambda t: t[1][::-1], reverse = True)  # reverse causes gaps to be at the end
    
    for pair in items:
        header = pair[0]
        outFile.write('>' + header + '\n')
        outFile.write(fastaDict[header] + '\n')
        outFile2.write(csvDict[header] + '\n')
        
    outFile.close()
    outFile2.close()
    
#getVirusSequences2('MF176251.1', 148, 966, 100, 'ortho_flanks', 'right')
getVirusSequences2('NC_001564.2', 0, None, 0, 'cfav', 'left')
exit()

###############################################################################

# from Bio.Phylo.TreeConstruction import DistanceCalculator
# from Bio import AlignIO
# from Bio.Phylo.TreeConstruction import DistanceTreeConstructor
# from Bio import Phylo
# 
# 
# def drawTree(fileName):
#     aln = AlignIO.read(fileName, 'fasta') # Bio.Align.MultipleSeqAlignment object
#     aln = aln[1:]
#     calculator = DistanceCalculator('identity')
# #    dm = calculator.get_distance(aln)
#     constructor = DistanceTreeConstructor(calculator, 'nj')
#     tree = constructor.build_tree(aln)
#     fig = pyplot.figure(figsize=(10, 20), dpi=100)
#     axes = fig.add_subplot(1, 1, 1)
#     matplotlib.rc('font', size=6)
#     Phylo.draw(tree, axes=axes)
#     pyplot.savefig('tree.pdf')
    
#drawTree('/home/havill/data2/results4/viruses/sequences/Flaviviridae/sequences_NC_001564.2.fasta')
#exit()

# populations = {'Angola': ['Angola'], 
#                'Argentina': ['Argentina', 'US_U'], 
#                'Australia': ['Australia'], 
#                'Brazil': ['Brazil'], 
#                'French Polynesia': ['FrenchPolynesia'],
#                'Gabon': ['Gabon'], 
#                'Mexico': ['Mexico'], 
#                'Philippines': ['Philippines'], 
#                'South Africa': ['South_Africa'], 
#                'Thailand': ['Thailand'], 
#                'USA': ['USA', 'AZ'], 
#                'Vietnam': ['Vietnam']}
#                
# AaegL5_hits = {'NC_001564.2': [9330, 8303],   # Cell fusing agent virus
#                'NC_034017.1': [1677, 839, 2269, 3801, 5188, 6416, 7537, 8296, 10083, 8327],  # Xishuangbanna aedes flavivirus
#                'NC_038512.1': [2579, 3877],   # Trichoplusia ni TED virus
#                'NC_035674.1': [446, 3299],    # Australian Anopheles totivirus
#                'NC_035132.1': [6568, 7294],   # Culex ohlsrhavirus
#                'NC_028484.1': [6555, 7282],   # Tongilchon ohlsrhavirus
#                'NC_040669.1': [6773, 7497],   # Riverside virus 1
#                'NC_025384.1': [6789, 6267]}   # Culex tritaeniorhynchus rhabdovirus
#                
# class Hit:
#     def __init__(self, pos = [], overlap = False, primary = True, data = None):
#         self._pos = pos
#         self._overlap = overlap
#         self._primary = primary
#         self._data = data
#         
#     def doesOverlap(self):
#         return self._overlap
#         
#     def isPrimary(self):
#         return self._primary
#         
#     def getData(self):
#         return self._data
#         
#     def getPos(self):
#         return self._pos
#         
#     def __len__(self):
#         return len(self._pos)
#         
#     def __getitem__(self, index):
#         return self._pos[index]
#         
#     def __contains__(self, interval):
#         try:
#             iterator = iter(interval)
#         except TypeError:
#             start = end = interval
#         else:
#             start, end = interval
#             if start > end:
#                 start, end = end, start
#                 
#         for index in range(0, len(self._pos), 2):
#             if (self._pos[index] <= self._pos[index + 1]) and ((start >= self._pos[index]) and (end <= self._pos[index + 1])):
#                 return True
#             elif (self._pos[index] > self._pos[index + 1]) and ((start <= self._pos[index]) and (end >= self._pos[index + 1])):
#                 return True
#         return False
#         
# #    def __setitem__(self, index, value)
#         
# #    def __iadd__(self, other)
#         
# class Hits:
#     def __init__(self):
#         self._hits = []   # list of hits
#         
#     def addHit(self, pos, overlap = False, primary = True, data = None):
#         self._hits.append(Hit(pos, overlap, primary, data))
#         
#     def addHits(self, otherHitsList):
#         self._hits.extend(otherHitsList)
#         
#     def getHits(self):
#         return self._hits
#         
#     def getPrimaryHitsOnly(self):
#         primaryHits = []
#         for hit in self._hits:
#             if hit.isPrimary():
#                 primaryHits.append(hit)
#         return primaryHits
#         
#     def adjust(self, adjustment):
#         for index in range(len(self._hits)):
#             for index2 in range(len(self._hits[index]._pos)):
#                 self._hits[index]._pos[index2] += adjustment
#                 
#     def __contains__(self, interval):
#         for hit in self._hits:
#             if interval in hit:
#                 return True
#         return False
#                 
#     def __len__(self):
#         return len(self._hits)
#             
# #    def __iadd__(self, other)
#         
# # class Hits:
# #     def __init__(self):
# #         self._hits = []   # list of lists
# #         self._overlaps = []
# #         
# #     def addHit(self, hit, overlap = False):
# #         self._hits.append(hit)
# #         self._overlaps.append(overlap)
# #         
# #     def getHits(self):
# #         return self._hits, self._overlaps
# 
# def getGenBank(accessionID):
# 	"""
# 	Query NCBI to get a gb file.
# 	"""
# 	
# 	Entrez.email = 'havill@denison.edu'
# 	try:
# 		handle = Entrez.efetch(db = 'nucleotide', id = accessionID, rettype = 'gb', retmode = 'text')
# 	except:
# 		print('Error retrieving GenBank record for ' + accessionID + ': ' + sys.exc_info()[1] + '\nDiagram not created.')
# 		return False
# 		
# 	gbFile = open(GB_DIR + accessionID + '.gb', 'w')
# 	gbFile.write(handle.read())
# 	gbFile.close()
# 	
# 	return True
# 
# def drawVirus(acc, family, hits, allSpecimens, separatePops, isFamily, showAaegL5_hits, showPlot, small, omitInReference):    
# 
# #    hits = getHits(acc)
# 
#     VIRUSES_DIR = RESULTS_DIR + 'viruses/'
#     if not Path(VIRUSES_DIR).exists():
#         os.system('mkdir ' + VIRUSES_DIR)
#     DIAGRAMS_DIR = VIRUSES_DIR + 'diagrams/'
#     if not Path(DIAGRAMS_DIR).exists():
#         os.system('mkdir ' + DIAGRAMS_DIR)
#     FAMILY_DIR = DIAGRAMS_DIR + family
#     if not Path(FAMILY_DIR).exists():
#         os.system('mkdir ' + FAMILY_DIR)
#     
# #     if not Path(outputDir).exists():
# #         os.mkdir(outputDir)
# 
#     if small:
#         width = 3
#         thickness = 1
#         ref_thickness = 2
#         fontsize = 1
#         plot_linewidth = 0.25 # 0.75
#         draw_line = False
#         with_ruler = False
#     #   ax.xaxis.set_tick_params(width=5)
#     else:
#         width = 10
#         thickness = 10
#         fontsize = 8
#         ref_thickness = 10
#         plot_linewidth = 0.75
#         draw_line = True
#         with_ruler = True
# 
#     
#     if not Path(GB_DIR + acc + '.gb').exists():
#         if not getGenBank(acc):
#             return
#     virusRecord = MyCustomTranslator().translate_record(GB_DIR + acc + '.gb')  # GraphicRecord
# #    virusRecord.plots_indexing = 'genbank'  # start at 1
#     codingStart = 1
#     codingEnd = virusRecord.sequence_length
#     for f in virusRecord.features:
#         f.thickness = ref_thickness
#         f.linewidth = 0
#         f.fontdict = {'size': fontsize}
#         if 'UTR' in f.label:
#             if '5' in f.label:
#                 codingStart = f.end + 1
#             elif '3' in f.label:
#                 codingEnd = f.start   # starts are actual start - 1 for some reason
#                 
#     virusRecord2 = SeqIO.read('/home/havill/data/aegypti/gb/' + acc + '.gb', 'gb')  # SeqRecord
#     if isFamily:
#         virusName = 'Family ' + family + ' (' + virusRecord2.id + ': ' + virusRecord2.description[:80]
#         if len(virusRecord2.description) > 80:
#             virusName += '...'
#         virusName += ')'
#         print('Creating diagram' + 's' * separatePops + ' for family ' + family + ' (using representative ' + virusRecord2.id + ')...')
#     else:
#         virusName = virusRecord2.id + ': ' + virusRecord2.description[:80]
#         if len(virusRecord2.description) > 80:
#             virusName += '...'
#         virusName += ' (' + family + ')'
#         print('Creating diagram' + 's' * separatePops + ' for ' + virusName + '...')
#                 
#     features = {}
#     for specimen in allSpecimens:
#         features[specimen] = []
#         if specimen in hits:
#             g = 1.0
#             for hit in hits[specimen].getHits():
#                 if hit.isPrimary():
#                     transparency = 0.75  # so overlap can be seen
#                 else:
#                     transparency = 0.25
#                 if len(hit) == 2:
#                     color = (0, 0, 1) + (transparency,) # blue
#                 else:
#                     color = (g, 0, 0) + (transparency,) # red
#                     g = max(g - 0.25, 0)
# #                 if hit.doesOverlap():
# #                     lc = (1, 1, 1) + (transparency,)    # white
# #                 else:
#                 lc = (0, 0, 0, 0)
#                 for index in range(0, len(hit), 2):
#                     start = hit[index]
#                     end = hit[index + 1]
#                     features[specimen].append(GraphicFeature(start=start, end=end, strand=1, color=color, label = None, thickness = thickness, linewidth = 0, linecolor = lc))
#                 if hit.getData() != []:  # add reference overlap regions
#                     lc = (1, 1, 1) + (transparency,)    # white
#                     for start, end in hit.getData():
#                         features[specimen].append(GraphicFeature(start=start, end=end, strand=1, color=color, label = None, thickness = thickness, linewidth = 1, linecolor = lc))
# 
#             
# #     for specimen in hits2:
# #         g = 1.0
# #         for hit in hits2[specimen].getHits():
# #             if len(hit) == 2:
# #                 color = (0, 0, 1, 0.25)   # translucent blue
# #             else:
# #                 color = (g, 0, 0, 0.25)   # translucent red
# #                 g = g - 0.25
# #             if hit.doesOverlap():
# #                 lc = (1, 1, 1, 0.25)      # white
# #             else:
# #                 lc = (0, 0, 0, 0)
# #             for index in range(0, len(hit), 2):
# #                 start = hit[index]
# #                 end = hit[index + 1]
# #                 features[specimen].append(GraphicFeature(start=start, end=end, strand=1, color=color, label = None, thickness = thickness, linewidth = 0, linecolor = lc))
#         
#     newFeatures = {}
#     for specimen in features:
#         for popName in populations:
#             for pattern in populations[popName]:
#                 if pattern in specimen:
#                     region = popName
#                     break
#         p = re.compile(r'[0-9]+')
#         num = int(p.findall(specimen)[0])
#         newFeatures[(region, num)] = features[specimen]
#         
#     pops = ['all']
#     if separatePops:
#         pops.extend(list(populations.keys()))
#         
#     for pop in pops:
#         if pop == 'all':
#             specimens = list(newFeatures.keys())
#         else:
#             specimens = [specimen for specimen in newFeatures if specimen[0] == pop]
#         specimens.sort()
#         
#         height = (6 * thickness / 225) * (len(specimens) + 3)
#         
#         numAxes = len(specimens) + 2
#         if showPlot:
#             numAxes +=1
#         if showAaegL5_hits and (acc in AaegL5_hits):
#             numAxes += 1
#             
#         height = (6 * thickness / 225) * numAxes
#         
#         fig, ax = pyplot.subplots(numAxes, 1, sharex = False, sharey = False, figsize = (width, height), gridspec_kw = {'height_ratios': [4] + [1] * (numAxes - 1)})
#         fig.subplots_adjust(top = 0.99, bottom = 0.01, left = 0.18, right = 0.95)
#         fig.suptitle(virusName, fontsize = fontsize, fontweight = 'bold', y = 0.995)
#         
#         virusRecord.box_linewidth = 0
#         virusRecord.plot(ax = ax[0], figure_width=width, draw_line = draw_line, with_ruler = with_ruler, max_line_length = 30, max_label_length = 30) # prevent multiline labels
#         x0, y0, w, h = ax[0].get_position().bounds
#         ax[0].set_position([x0, y0+h*0.2, w, h])
#         ax[0].tick_params(labelsize = fontsize)
#         xlims = ax[0].get_xlim()
#         
#         AaegL5_features = []
#         if acc in AaegL5_hits:
#             for index in range(0, len(AaegL5_hits[acc]), 2):
#                 start = AaegL5_hits[acc][index]
#                 end = AaegL5_hits[acc][index + 1]
#                 AaegL5_features.append(GraphicFeature(start=start, end=end, strand=0, color='lightgray', label = None, thickness = thickness, linewidth = 0))
#     
#         x = range(virusRecord.sequence_length)
#         y = [0] * virusRecord.sequence_length
#         
#         for specimen in specimens:
#             for feature in newFeatures[specimen]:
#                 start = min(max(feature.start, 1), virusRecord.sequence_length - 1)
#                 end = min(feature.end, virusRecord.sequence_length - 1)
#                 for i in range(min(start, end), max(start, end) + 1):
#                     try:
#                         y[i-1] += 1
#                     except IndexError:
#                         print(i, virusRecord.sequence_length - 1)
#                         return
#                     
#         for feature in AaegL5_features:
#             start = feature.start
#             end = feature.end
#             for i in range(min(start, end), max(start, end) + 1):
#                 y[i-1] += 1
#                 
#         axisIndex = 1
#         y = [y[i] / len(specimens) for i in range(len(y))]
#         if showPlot:
#             ax[1].plot(x, y, color = 'blue', linewidth = plot_linewidth)
#             ax[1].set_xlim(*xlims)
#             ax[1].set_ylim(0, 2)
#             ax[1].axis('off')
#             ax[1].set_frame_on(False)
#             axisIndex += 1
#     
#         y = [y[i] > 0 for i in range(len(y))]
#         ax[axisIndex].bar(x, y, color = 'blue')
#         ax[axisIndex].set_xlim(*xlims)
#         ax[axisIndex].set_ylim(0, 2)
#         ax[axisIndex].axis('off')
#         ax[axisIndex].set_frame_on(False)
#         coverage = (sum(y[codingStart-1:codingEnd]) / (codingEnd - codingStart + 1)) * 100
#         ax[axisIndex].text(virusRecord.sequence_length, -0.1, '{0:4.1f}%'.format(coverage) + ' ', horizontalalignment = 'left', fontdict = {'size': fontsize})
#     
#         axisIndex += 1
#         if showAaegL5_hits and (acc in AaegL5_hits):
#             record = GraphicRecord(features = AaegL5_features, sequence_length = virusRecord.sequence_length, feature_level_height = 0)
#             record.plot(ax = ax[axisIndex], figure_width=width, with_ruler = False, draw_line = draw_line)
#             
#             start, end = record.span
#             plot_start, plot_end = start - 0.8, end - 0.2
#             ax[axisIndex].plot([plot_start, plot_end], [0, 0], linewidth = 0.25, zorder=-1000, c="k")
#             
#             ax[axisIndex].text(0, -0.1, 'AaegL5 ', horizontalalignment = 'right', fontdict = {'size': fontsize})
#             axisIndex += 1
#        
#         for specimen in specimens:
#             specimenName = specimen[0] + ' ' + str(specimen[1])
# #            print(specimenName)
#             record = GraphicRecord(features = newFeatures[specimen], sequence_length = virusRecord.sequence_length, feature_level_height = 0)
#             record.plot(ax = ax[axisIndex], figure_width=width, with_ruler = False, draw_line = draw_line)
#             
#             start, end = record.span
#             plot_start, plot_end = start - 0.8, end - 0.2
#             ax[axisIndex].plot([plot_start, plot_end], [0, 0], linewidth = 0.25, zorder=-1000, c="k")
#             
#             ax[axisIndex].text(0, -0.1, specimenName + ' ', horizontalalignment = 'right', fontdict = {'size': fontsize})
#             
#             for p in ax[axisIndex].patches:
#                 if p.get_edgecolor()[:3] == (1, 1, 1):  # white
#                     p.set_hatch('////')
#                     
#             
#             axisIndex += 1
#         if isFamily:
#             figName = family + '_' + pop + '.pdf'
#         else:
#             figName = acc + '_' + pop + '.pdf'
#         fig.savefig(Path(FAMILY_DIR) / figName)
#         pyplot.close(fig)

# def getHits_old(acc):
#     """Return all hits for given accession number as a dictionary
#        with key = specimen and value = [start, end, start, end, ...]"""
#        
#     resultsFile = open(RESULTS_DIR + 'results_all.tsv', 'r')
#     header1 = resultsFile.readline()
#     header2 = resultsFile.readline()
#     accs = header1.split('\t')
#     index = accs.index(acc)
#     p = re.compile(r'[0-9]+-[0-9]+')
#     allHits = {}
#     for line in resultsFile:
#         cols = line.rstrip().split('\t')
#         specimen = cols[0]
#         hits = cols[index]
#         hits = p.findall(hits)
#         allHits[specimen] = []
#         for h in hits:
#             pos = h.split('-')
#             start = int(pos[0])
#             end = int(pos[1])
#             allHits[specimen].append([start, end])
#     resultsFile.close()
#     
#     return allHits
    
# def getHits(acc):
#     dir = Path(RESULTS_DIR + 'specimens/')
#     subdirs = [d for d in dir.iterdir()]
#     files = [d / ('xml/' + d.name + '_all_hits_features.xml') for d in subdirs ]
#     files.sort()
# 
#     hits = {}
#     hits2 = {}
#     for file in files:
#         specimen = file.name.rstrip('_all_hits_features.xml')
#         print(specimen)
#         hits[specimen] = Hits()
#         hits2[specimen] = Hits()
#         tree = ET.parse(str(file))
#         root = tree.getroot()
#         for contig in root:
#             v = contig.find('virushit')  # first virus hit
#             seqid = v.attrib['seqid']
#             if '|' in seqid:
#                 seqid = seqid.split('|')[1]
#             if seqid == acc:
#                 virushits = contig.findall('virushit')
#                 hit = []
#                 for virushit in virushits:
#                     start = int(virushit.find('sstart').text)
#                     end = int(virushit.find('send').text)
#                     hit.extend([start, end])
#                 overlap = len(contig.findall('vectorhitoverlap')) > 0
#                 if (contig.attrib['besthit'] == 'True') and (len(virushits) < 5):
#                     hits[specimen].addHit(hit, overlap)
#                 else:
#                     hits2[specimen].addHit(hit, overlap)
#                     print('*', end='')
#                 print(hit)
#     return hits, hits2

# def getHits(acc):
#     dir = Path(RESULTS_DIR + 'specimens/')
#     subdirs = [d for d in dir.iterdir()]
#     files = [d / ('xml/' + d.name + '_all_hits_features.xml') for d in subdirs ]
#     files.sort()
# 
#     hits = {}
#     for file in files:
#         specimen = file.name.rstrip('_all_hits_features.xml')
#         print(specimen)
#         hits[specimen] = Hits()
#         tree = ET.parse(str(file))
#         root = tree.getroot()
#         for contig in root:
#             v = contig.find('virushit')  # first virus hit
#             seqid = v.attrib['seqid']
#             if '|' in seqid:
#                 seqid = seqid.split('|')[1]
#             if seqid == acc:
#                 virushits = contig.findall('virushit')
#                 hit = []
#                 for virushit in virushits:
#                     start = int(virushit.find('sstart').text)
#                     end = int(virushit.find('send').text)
#                     hit.extend([start, end])
#                 overlap = len(contig.findall('vectorhitoverlap')) > 0
#                 if (contig.attrib['besthit'] == 'True'): # and (len(virushits) < 5):
#                     hits[specimen].addHit(hit, overlap, True)
#                 else:
#                     hits[specimen].addHit(hit, overlap, False)
#                     print('*', end='')
#                 print(hit)
#     return hits

# def getHits(omitInReference = True):
#     dir = Path(RESULTS_DIR + 'specimens/')
#     subdirs = [d for d in dir.iterdir()]
#     files = [d / ('xml/' + d.name + '_all_hits_features.xml') for d in subdirs ]
#     files.sort()
# 
#     allVirusHits = {}
#     allSpecimens = []
#     for file in files:
#         specimen = file.name.rstrip('_all_hits_features.xml')
#         print(specimen)
#         allSpecimens.append(specimen)
#         tree = ET.parse(str(file))
#         root = tree.getroot()
#         for contig in root:
#             if 'inreference' in contig.attrib:
#                 overlap = contig.attrib['inreference'] == 'True'
#             else:
#                 overlap = False
#             if 'referenceoverlaps' in contig.attrib:
#                 referenceOverlaps = eval(contig.attrib['referenceoverlaps'])
#             else:
#                 referenceOverlaps = []
#             if omitInReference and overlap:
#                 continue
#             virushits = contig.findall('virushit')
#             v = virushits[0]  # first virus hit
#             seqid = v.attrib['seqid']
#             if '|' in seqid:
#                 seqid = seqid.split('|')[1]
#             if seqid not in allVirusHits:
#                 allVirusHits[seqid] = {}
#             if specimen not in allVirusHits[seqid]:
#                 allVirusHits[seqid][specimen] = Hits()
# 
#             hit = []
#             for virushit in virushits:
#                 start = int(virushit.find('sstart').text)
#                 end = int(virushit.find('send').text)
#                 hit.extend([start, end])
#             primary = contig.attrib['besthit'] == 'True'  # and (len(virushits) < 5):
#             allVirusHits[seqid][specimen].addHit(hit, overlap, primary, referenceOverlaps)
# #             if not primary:
# #                 print('*', end='')
# #             print(seqid, hit)
#     return allVirusHits, allSpecimens

# def drawFamily(family, repACC, separatePops = False):
#     allVirusHits, allSpecimens = getHits()
#         
#     fams = readFamFile()
#     
#     famHits = {}
#     for acc in allVirusHits:
#         if acc not in fams:
#             fam = getFamily(acc)
#             fams[acc] = fam
#             addFamily(acc, fam)
#         if fams[acc] == family:
#             for specimen in allVirusHits[acc]:
#                 if specimen not in famHits:
#                     famHits[specimen] = Hits()
#                 famHits[specimen].addHits(allVirusHits[acc][specimen])
#                 
#     drawVirus(repACC, family, famHits, allSpecimens, separatePops, True)

# def drawFamily(families, famACCs, separatePops, showAaegL5_hits, showPlot, small, omitInReference):
# #     if isinstance(families, str):
# #         if not isinstance(repACCs, str):
# #             print('The family and repACC must both be strings.')
# #             return
# #         families = [families]
# #         repACCs = [repACCs]
# #     else:
# #         if famACCs is not None:
# #             print('family and famACCs cannot both be lists.')
# #             return
# #         if not isinstance(repACCs, list):
# #             print('The family and repACC must both be iterable.')
# #             return
# 
#     allVirusHits, allSpecimens = getHits(omitInReference)
#         
#     if famACCs is None:
#         fams = readFamFile()
#     
#     for family, repACC in families:
#         if famACCs is None:  # if False, must be only one family in families
#             famACCs = []
#             for acc in allVirusHits:
#                 if acc not in fams:
#                     fam = getFamily(acc)
#                     fams[acc] = fam
#                     addFamily(acc, fam)
#                 if fams[acc] == family:
#                     famACCs.append(acc)
#         
#         startPositions = {}
#         for acc in famACCs:
#             if not Path(GB_DIR + acc + '.gb').exists():
#                 if not getGenBank(acc):
#                     continue
#             virusRecord = SeqIO.read('/home/havill/data/aegypti/gb/' + acc + '.gb', 'gb')  # SeqRecord
#             minCDSPos = math.inf
#             for feature in virusRecord.features:
#                 if feature.type == 'CDS':
#                     minCDSPos = min(minCDSPos, feature.location.start)
#             startPositions[acc] = minCDSPos
#     
#         famHits = {}
#         for acc in famACCs:
#             for specimen in allVirusHits[acc]:
#                 if specimen not in famHits:
#                     famHits[specimen] = Hits()
#                 allVirusHits[acc][specimen].adjust(startPositions[repACC] - startPositions[acc])
#                 famHits[specimen].addHits(allVirusHits[acc][specimen].getPrimaryHitsOnly())  # only primary hits so one hit per contig
#                 
#         drawVirus(repACC, family, famHits, allSpecimens, separatePops, True, showAaegL5_hits, showPlot, small, omitInReference)
#         famACCs = None
#     
# def drawAll(separatePops, omitInReference):
#     allVirusHits, allSpecimens = getHits(omitInReference)
#     
#     fams = readFamFile()
# 
#     for acc in allVirusHits:
#         doVirus = False  # only draw virus if it has at least one primary hit
#         for specimen in allVirusHits[acc]:
#             for hit in allVirusHits[acc][specimen].getHits():
#                 if hit.isPrimary():
#                     doVirus = True
#                     break
#             if doVirus:
#                 break
#         if doVirus:
#             if acc not in fams:
#                 fam = getFamily(acc)
#                 fams[acc] = fam
#                 addFamily(acc, fam)
#             drawVirus(acc, fams[acc], allVirusHits[acc], allSpecimens, separatePops, False, True, True, False, omitInReference)

# def main():
# #     viruses = ['NC_001564.2',   # cfav
# #                'NC_007740.1',   # liao ning segment 5
# #                'NC_035132.1',   # culex ohlsrhavirus
# #                'NC_038512.1',   # ted
# #                'NC_035674.1',   # australian-anopheles 
# #                'NC_034017.1',   # xishuangbanna
# #                'NC_028484.1',   # tongilchon-ohlsrhavirus
# #                'NC_025384.1',   # culex-tritaeniorhynchus
# #                'NC_040669.1',   # riverside-virus-1
# #                'MH430659.1', 
# #                'MH037149.1',
# #                'KY768856.1',
# #                'NC_038263.1',
# #                'MF176337.1',
# #                'MF176381.1',
# #                'MH188035.1',
# #                'MT822181.1']
# #            
# #     viruses = ['NC_001564.2'] # cfav
#     
#     allVirusHits, allSpecimens = getHits()
#     
#     # HACK
# #     for specimen in allVirusHits['MF176381.1']:
# #         if specimen not in allVirusHits['MF176337.1']:
# #             allVirusHits['MF176337.1'][specimen] = Hits()
# #         hits = allVirusHits['MF176381.1'][specimen].getHits()
# #         for hit in  hits:
# #             for index in range(len(hit._pos)):
# #                  hit._pos[index] -= 3
# #             allVirusHits['MF176337.1'][specimen].addHit(hit)
#     
#     fams = readFamFile()
#     
#     family = 'Flaviviridae'
#     reprACC = 'NC_001564.2'
#     famHits = {}
#     for acc in allVirusHits:
#         if acc not in fams:
#             fam = getFamily(acc)
#             fams[acc] = fam
#             addFamily(acc, fam)
#         if fams[acc] == family:
#             for specimen in allVirusHits[acc]:
#                 if specimen not in famHits:
#                     famHits[specimen] = Hits()
#                 famHits[specimen].addHits(allVirusHits[acc][specimen])
#     drawVirus(reprACC, family, famHits, allSpecimens, False)
#     return
#     
# #     acc = 'MF176337.1'
# #     drawVirus(acc, fams[acc], allVirusHits[acc], allSpecimens, False)
#         
#     for acc in allVirusHits:
#         doVirus = False  # only draw virus if it has at least one primary hit
#         for specimen in allVirusHits[acc]:
#             for hit in allVirusHits[acc][specimen].getHits():
#                 if hit.isPrimary():
#                     doVirus = True
#                     break
#             if doVirus:
#                 break
#         if doVirus:
#             if acc not in fams:
#                 fam = getFamily(acc)
#                 fams[acc] = fam
#                 addFamily(acc, fam)
#             drawVirus(acc, fams[acc], allVirusHits[acc], allSpecimens, True)

#if __name__ == '__main__':
#    drawFamily(family, repACC, famACCs, separatePops, showAaegL5_hits, showPlot, small, omitInReference)
#    drawFamily('Flaviviridae', 'NC_001564.2', None, False, False, True, False, False)
#    drawFamily('Xinmoviridae', 'MH037149.1', None, False, False, True, False, False)
#    drawFamily('Rhabdoviridae', 'NC_035132.1', None, False, False, True, False, False)
#    drawFamily('Rhabdoviridae', 'NC_035132.1', ['NC_035132.1', 'KY768856.1', 'MH188003.1', 'MF176358.1', 'NC_028484.1'], False)
#    drawFamily('Orthomyxoviridae', 'MF176251.1', None, False, False, True, False, False)
#    drawFamily('Phenuiviridae', 'NC_038263.1', None, False, False, True, False, False)
#    drawFamily('Totiviridae', 'NC_035674.1', None, False, False, True, False, False)
    
#    drawFamily(['Flaviviridae', 'Orthomyxoviridae', 'Phenuiviridae', 'Rhabdoviridae', 'Totiviridae', 'Xinmoviridae'], ['NC_001564.2', 'MF176251.1', 'NC_038263.1', 'NC_035132.1', 'NC_035674.1', 'MH037149.1'], None, False, False, True, False, False)
#    drawFamily([('Flaviviridae', 'NC_001564.2'), ('Orthomyxoviridae', 'MF176251.1'), ('Phenuiviridae', 'NC_038263.1'), ('Rhabdoviridae', 'NC_035132.1'), ('Totiviridae', 'NC_035674.1'), ('Xinmoviridae', 'MH037149.1')],None, False, False, True, False, False)
#    drawAll(False, False)
