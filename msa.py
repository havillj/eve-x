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

RESULTS_DIR = '/home/havill/data2/results4/'
GB_DIR = '/home/havill/data/aegypti/gb/'
    
def heatmap(data, row_labels, col_labels, ax=None,
            cbar_kw={}, cbarlabel="", **kwargs):
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

    return im


def annotate_heatmap(im, data=None, valfmt="{x:.2f}",
                     textcolors=("black", "white"),
                     threshold=None, **textkw):
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
    
def makeHeatMap(values, xlabels, ylabels, ids, boundsFile):
    fig, ax = pyplot.subplots(2, 1, figsize = (100, 20), dpi=200, sharex = True, gridspec_kw = {'height_ratios': [3, 222], 'hspace': 0.01})
    fig.subplots_adjust(top = 0.95, bottom = 0.05, left = 0, right = 1.0)
    im = ax[1].imshow(values, interpolation='none')
    # Turn spines off and create white grid.
    for edge, spine in ax[0].spines.items():
        spine.set_visible(False)
    ax[1].set_yticks(np.arange(values.shape[0]))
    ax[1].set_yticklabels(ylabels, fontfamily='serif', fontsize=5)
    xlims = ax[1].get_xlim()
    ax[0].bar(range(len(ids)), ids, color = 'blue')
    ax[0].set_xlim(*xlims)
    x1, y1, w1, h1 = ax[1].get_position().bounds
    x0, y0, w0, h0 = ax[0].get_position().bounds
    ax[0].set_position([x1, y0, w1, h0])
    ax[0].axis('off')
    ax[0].set_frame_on(False)
    
    kw = dict(horizontalalignment="center", verticalalignment="center", fontsize=6, color='black')
    try:
        bfile = open(boundsFile, 'r')
    except:
        pass
    else:
        i = 1
        for line in bfile:
            cols = line.rstrip().split(',')
            left = int(cols[1])
            right = int(cols[1]) + int(cols[2]) - 1
            for j in range(left, right + 1):
                if not np.array_equal(values[i,j], [1., 1., 1.]):
                    im.axes.text(j, i, '+', **kw)
            i += 1
        bfile.close()
        
#    im = heatmap(values, ylabels, xlabels, ax=ax, cmap="Blues", cbarlabel="")
#    texts = annotate_heatmap(im, valfmt="{x:.2f}", fontsize=8)
    
#    fig.tight_layout()
    return fig


d = {'A': (0., 0., 1.), 'C': (1., 0., 0.), 'G': (0., 1., 0.), 'T': (1., 1., 0.)}

def drawMSA(fileName, refID, boundsFile):
    table1 = []
    table2 = []
    table3 = []
    ylabels = []
    records = SeqIO.index(fileName, 'fasta')
    for id in records:
        if refID in id:
            refseq = str(records[id].seq)
            break
    ids = np.zeros(len(refseq))
    totals = np.zeros(len(refseq))
    for seqid in records:
        ylabels.append(seqid)
        row1 = []
        row2 = []
        row3 = []
        index = 0
        for char in str(records[seqid].seq).upper():
            row1.append(d.get(char, (1., 1., 1.))[0])
            row2.append(d.get(char, (1., 1., 1.))[1])
            row3.append(d.get(char, (1., 1., 1.))[2])
            if seqid != refID:
                if char != '-':
                    if char == refseq[index]:
                        ids[index] += 1
                    totals[index] += 1
            index += 1
        table1.append(row1)
        table2.append(row2)
        table3.append(row3)
    
    ids /= totals

    # pyplot.rcParams.update({
    #     "text.usetex": True,
    #     "font.family": "serif",
    #     "font.serif": ["Computer Modern Roman"],
    #     "font.sans-serif": ["Computer Modern Sans serif"]})
    # #     
    # # ##Needed to install additional latex packages under ubuntu:
    # # ##sudo apt-get install dvipng texlive-latex-extra texlive-fonts-recommended cm-super
    # # 
    values = np.dstack([table1, table2, table3])
    #print(values)
    fig = makeHeatMap(values, [], ylabels, ids, boundsFile)
    fig.savefig(fileName.split('.fasta')[0] + '_msa.png', dpi=200)

if __name__ == '__main__':
    drawMSA('/home/havill/data2/ortho_flanks.fasta', 'MF176251.1', '/home/havill/data2/ortho_flanks.csv')
