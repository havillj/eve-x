from lxml import etree
import os
from pathlib import Path
import sys
from Bio import SeqIO, Entrez, AlignIO, Seq
from util import *
import csv

###############################################################################

###############################################################################

def getOverlaps(dir):
    """
    Finds the overlapping flanking regions from the cell fusing agent virus in a specified document.

    Parameters:
        Path: A path to the specified xml file

    Return value: a dictionary of dictionaries
    """

    path = f"/Volumes/Data2/specimens/{dir}/results/xml/{dir}_hits.xml"

    parser = etree.XMLParser(remove_blank_text = True)
    mosc = etree.parse(path, parser).getroot()

    xstring = f"""/root/contig[@besthit="True"]/@name"""
    contigs = mosc.xpath(xstring)

    overlaps = {}
    for i in contigs:
        contigname = mosc.xpath(f"""/root/contig[@besthit="True"][@name="{i}"]/virushit/@stitle""")[0]

        contigOverlaps = {}
        left  = mosc.xpath(f"""/root/contig[@name="{i}"]/flanks/match/@leftid""")
        right = mosc.xpath(f"""/root/contig[@name="{i}"]/flanks/match/@rightid""")

        for index in range(len(left)):
            acc = mosc.xpath(f"""/root/contig[@name="{i}"]/vectorhitleft[@id="{left[index]}"]/@seqid""")[0] #accession id for node

            lstart = int(mosc.xpath(f"""/root/contig[@name="{i}"]/vectorhitleft[@id="{left[index]}"]/sstart/text()""")[0])  #start location in aa genome
            lend = int(mosc.xpath(f"""/root/contig[@name="{i}"]/vectorhitleft[@id="{left[index]}"]/send/text()""")[0])  #end location in aa genome
            qend = int(mosc.xpath(f"""/root/contig[@name="{i}"]/vectorhitleft[@id="{left[index]}"]/qend/text()""")[0]) #end location in contig for node
            vqstart = int(mosc.xpath(f"""/root/contig[@name="{i}"]/virushit/qstart/text()""")[0]) #start location in contig for viral hit

            if lstart > lend:   #check if sequence is backwards
                lastLeft = lstart
            else:
                lastLeft = lend

            rstart = int(mosc.xpath(f"""/root/contig[@name="{i}"]/vectorhitright[@id="{right[index]}"]/sstart/text()""")[0])
            rend = int(mosc.xpath(f"""/root/contig[@name="{i}"]/vectorhitright[@id="{right[index]}"]/send/text()""")[0])
            qstart = int(mosc.xpath(f"""/root/contig[@name="{i}"]/vectorhitright[@id="{right[index]}"]/qstart/text()""")[0]) #start location in contig for node
            vqend = int(mosc.xpath(f"""/root/contig[@name="{i}"]/virushit/qend/text()""")[0])   #end location in contig for viral hit

            if rstart > rend:   #check if sequence is backwards
                firstRight = rend
            else:
                firstRight = rstart

            querySeq = mosc.xpath(f"""/root/contig[@name="{i}"]/virushit/qseq/text()""")[0]     #qseq for viralhit
            leftOverlap = ""
            rightOverlap = ""
            overlapSeq = ""

            scaffoldsPath = f"/Volumes/Data2/specimens/{dir}"

            if qend > vqstart:
                leftOverlap = querySeq[: qend - vqstart]
            if qstart < vqend:
                rightOverlap = querySeq[- (vqend - qstart):]
            if lastLeft > firstRight:
                overlapSeq = getContig(scaffoldsPath, i.split("__")[0], firstRight, lastLeft)

            overlapInfo = [dir, contigname,lastLeft,firstRight, overlapSeq, qend - qstart, firstRight - lastLeft, leftOverlap, rightOverlap]    #distance in the contig, distance in the aa genome

            """if overlapInfo[3] >= 0 and overlapInfo[4] >= 0:  #checking for insertions in contig
                overlapInfo.append(overlapInfo[3] - overlapInfo[4])#inserted bases in contig
            else:
                overlapInfo.append("NaN")"""

            contigOverlaps[acc] = overlapInfo
            overlaps[i] = contigOverlaps


    return overlaps

###############################################################################

def getGen(accessionID, start, stop):
    """Retrieves Aedes Aegypti reference genome.

       Parameters:
            accessionID: an accession ID
            start:       starting position of sequence
            stop:        end position of sequence

    Return value: nucleotide sequence of the specific genome region
    """
    file = f"""{accessionID}.gbk"""
    Entrez.email = 'havill@denison.edu'
    with Entrez.efetch(db = "nucleotide", rettype = "gb", retmode = "text", id = accessionID, seq_start = start , seq_stop = stop ) as handle:
        seq_record = SeqIO.read(handle, "gb")
        myseq = seq_record.seq
        return myseq

###############################################################################
def consolidate():
    """Finds and writes overlaps to a TSV file

       Parameters:
            None:

    Return value: None
    """
    xmlFiles = []
    paths = os.listdir(SPECIMENS_DIR)

    for path in paths:
        xmlFiles.append(path)

    with open("allOverlaps.tsv", "w") as outfile:

        tsv_writer = csv.writer(outfile, delimiter = "\t")
        tsv_writer.writerow(["Specimen Name", "Node", "Accession ID","Virus Name", "End Position", "Start Position","Sequence Overlap", "Contig Flanks Distance", "AA Genome Flanks Distance", "Left Overlap", "Right Overlap"])

        count = 0
        for file in xmlFiles:
        #for file in range(2):
            count += 1
            overlaps = getOverlaps(file)
            print(count)
            for i in overlaps:
                for contig in overlaps[i]:
                    sub = overlaps[i][contig]
                    overlapL = [sub[0]] + [i, contig] + sub[1:]
                    tsv_writer.writerow(overlapL)


###############################################################################


def main():
    consolidate()


main()
