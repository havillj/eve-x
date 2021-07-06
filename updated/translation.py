import matplotlib
matplotlib.use('agg')
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as pyplot
from Bio import SeqIO, Entrez, AlignIO
from Bio.Seq import Seq
from Bio.Align import MultipleSeqAlignment
from Bio.SeqFeature import SeqFeature, FeatureLocation
from dna_features_viewer import GraphicFeature, GraphicRecord, BiopythonTranslator
import math


path = "/Volumes/Data2/specimens/Amacuzac-Mexico-13.LIN210A1627/results/sequences/Amacuzac-Mexico-13.LIN210A1627_hits_aligned.fasta"

# with open(path, 'r') as file:
#     count = 1
#     for line in file:
#         print(count)
#         if line[0] != ">":
#             print("New contig:\n")
#             seq = Seq(line.replace('-', 'N'))
#             for index in range(3):
#                 end = math.floor(len(seq[index:])/3)*3
#                 print(end)
#                 print(seq[index:end].translate(to_stop=True))
#         count += 1

#
# with open(path, 'r') as f:
#     for line in f:
#         if line[0] != '>':
#             print("New contig:")
#             seq = Seq(line.replace('-', 'N').strip())
#             for index in range(3):
#                 to_translate = seq[index:len(seq[index:]) - (len(seq[index:])%3)]
#                 print(to_translate.translate(to_stop = True))


sequence = "GTCAACAGGACGATCATGAAGCATTTTATCCGGTTCTTTAAGAATCGCAACTTACGACCGCGTATTTTATCAAAAGAAGAATTCATCGCCAACGTTCGAAATGATGCGGCTATAGGATCGTGGAGCAGGGATGTGCCATGGCGAGATGTGCAAGAAGCCATTCAGGACCAGTGCTTCTGGGATCCCATC"
#seq2 = "ATGCT--TGTGACTGGAGTTCAGACGTGTGCTCTTCCGAT"
seq = Seq(sequence.strip().replace('-', 'N'))
startLoc = [seq.find('ATG')]
for x in startLoc:
    loc = seq.find('ATG',x+1)
    #print("running")
    if loc != -1:
        startLoc.append(loc)
print(startLoc)
for loc in startLoc:
    to_translate = seq[loc:len(seq[loc:]) - (len(seq[loc:])%3)]
    print(seq[loc:].translate(to_stop = True,stop_symbol='*'))
