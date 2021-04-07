from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio import SeqIO, Entrez
from dna_features_viewer import GraphicFeature, GraphicRecord, BiopythonTranslator
import matplotlib.pyplot as pyplot
from pathlib import Path
import re
from lxml import etree as ET

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


def tryThis():
    acc = 'MF176337.1'
    pop = 'all'
    family = 'foo'
    
    width = 10
    thickness = 10
    fontsize = 8
    
    virusRecord = MyCustomTranslator().translate_record(GB_DIR + acc + '.gb')  # GraphicRecord
    virusRecord.plots_indexing = 'genbank'  # start at 1
    codingStart = 1
    codingEnd = virusRecord.sequence_length
    for f in virusRecord.features:
        f.thickness = thickness
        f.linewidth = 0.5
        f.fontdict = {'size': fontsize}
        if 'UTR' in f.label:
            if '5' in f.label:
                codingStart = f.end + 1
            elif '3' in f.label:
                codingEnd = f.start   # starts are actual start - 1 for some reason
                
    virusRecord2 = SeqIO.read('/home/havill/data/aegypti/gb/' + acc + '.gb', 'gb')  # SeqRecord
    virusName = virusRecord2.id + ': ' + virusRecord2.description[:40]
    if len(virusRecord2.description) > 40:
        virusName += '...'
    virusName += ' (' + family + ')'
                
    fig, ax = pyplot.subplots(2, 1, sharex = False, sharey = False)
#     fig.subplots_adjust(top = 0.99, bottom = 0.01, left = 0.18, right = 0.95)
#     fig.suptitle(virusName, fontsize = fontsize + 1, fontweight = 'bold', y = 0.995)
    
    virusRecord.plot(ax = ax[0], figure_width=width, strand_in_label_threshold=7, max_line_length = 30, max_label_length = 30) # prevent multiline labels
#     x0, y0, w, h = ax[0].get_position().bounds
#     print(x0, y0, w, h)
#     ax[0].set_position([x0, y0+h*0.5, w, h*1.5])
#     print(ax[0].get_xticks())
        
    fig.savefig(acc + '_' + pop + '.pdf')
        
tryThis()
