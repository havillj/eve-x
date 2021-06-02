from lxml import etree
import os
from pathlib import Path
import sys
from Bio import SeqIO, Entrez, AlignIO, Seq
from util import *
import csv

###############################################################################

###############################################################################

def getOverlaps(path):
    """
    Finds the overlapping flanking regions from the cell fusing agent virus in a specified document.

    Parameters:
        Path: A path to the specified xml file

    Return value: a dictionary of dictionaries
    """

    parser = etree.XMLParser(remove_blank_text = True)
    mosc = etree.parse(path, parser).getroot()

    xstring = f"""/root/contig[@besthit="True"]/@name"""
    contigs = mosc.xpath(xstring)

    overlaps = {}
    for i in contigs:
        contigname = mosc.xpath(f"""/root/contig[@besthit="True"][@name="{i}"]/virushit/@stitle""")[0]

        contigOverlaps = {}
        l = f"""/root/contig[@name="{i}"]/flanks/match/@leftid"""
        r = f"""/root/contig[@name="{i}"]/flanks/match/@rightid"""
        left  = mosc.xpath(l)
        right = mosc.xpath(r)

        for index in range(len(left)):
            lstart = f"""/root/contig[@name="{i}"]/vectorhitleft[@id="{left[index]}"]/sstart/text()"""
            lend = f"""/root/contig[@name="{i}"]/vectorhitleft[@id="{left[index]}"]/send/text()"""
            qend = int(mosc.xpath(f"""/root/contig[@name="{i}"]/vectorhitleft[@id="{left[index]}"]/qend/text()""")[0]) #qend for aa left hit
            vqstart = int(mosc.xpath(f"""/root/contig[@name="{i}"]/virushit/qstart/text()""")[0]) #qstart for viralhit
            acc = mosc.xpath(f"""/root/contig[@name="{i}"]/vectorhitleft[@id="{left[index]}"]/@seqid""")[0]
            if mosc.xpath(lstart)[0] > mosc.xpath(lend)[0]:
                Lrev = True
                lastLeft = int(mosc.xpath(lstart)[0])
            else:
                Lrev = False
                lastLeft = int(mosc.xpath(lend)[0])

            rstart = f"""/root/contig[@name="{i}"]/vectorhitright[@id="{right[index]}"]/sstart/text()"""
            rend = f"""/root/contig[@name="{i}"]/vectorhitright[@id="{right[index]}"]/send/text()"""
            qstart = int(mosc.xpath(f"""/root/contig[@name="{i}"]/vectorhitright[@id="{right[index]}"]/qstart/text()""")[0]) #qstart for aa right hit
            vqend = int(mosc.xpath(f"""/root/contig[@name="{i}"]/virushit/qend/text()""")[0])#qend for viralhit

            if mosc.xpath(rstart)[0] > mosc.xpath(rend)[0]:
                Rrev = True
                firstRight = int(mosc.xpath(rend)[0])
            else:
                Rrev = False
                firstRight = int(mosc.xpath(rstart)[0])
            if firstRight < lastLeft:
                querySeq = mosc.xpath(f"""/root/contig[@name="{i}"]/virushit/qseq/text()""")[0]#qseq for viralhit
                #overlapSeq = str(getGen(acc, firstRight, lastLeft))
                overlapInfo = {"virus name": contigname, "lastLeftPosition": lastLeft, "firstRightPosition": firstRight}#, "overlap": overlapSeq}
                contigOverlaps[acc] = overlapInfo
                if qend > vqstart:
                    contigOverlaps[acc]["queryOverlapLeft"] = querySeq[: qend - vqstart]
                if qstart < vqend:
                    contigOverlaps[acc]["queryOverlapRight"] = querySeq[- (vqend - qstart):]
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
    paths = os.listdir("/Volumes/Data2/specimens")
    for path in paths:
        xmlFiles.append(f"/Volumes/Data2/specimens/{path}/results/xml/{path}_hits.xml")

    with open("allOverlaps.tsv", "w") as outfile:

        tsv_writer = csv.writer(outfile, delimiter = "\t")
        tsv_writer.writerow(["Specimen Name", "Node", "Accession ID","Virus Name", "End Position", "Start Position", "Sequence Overlap", "Left Overlap", "Right Overlap"])


        for file in xmlFiles:
            specimenName = file.split("/")[-1][:-9]
            overlaps = getOverlaps(file)
            for i in overlaps:
                for contig in overlaps[i]:
                    sub = list(overlaps[i][contig].values())
                    overlapL = [specimenName, i, contig] + sub
                    tsv_writer.writerow(overlapL)


###############################################################################


def main():

    consolidate()

main()
