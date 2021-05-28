from lxml import etree
import os
from pathlib import Path
import sys
from Bio import SeqIO, Entrez, AlignIO, Seq


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

    xstring = f"""/root/contig[@besthit="True"][virushit/@stitle="Cell fusing agent virus"]/@name"""
    contigs = mosc.xpath(xstring)
    overlaps = {}
    for i in contigs:
        contigOverlaps = {}
        l = f"""/root/contig[@name="{i}"]/flanks/match/@leftid"""
        r = f"""/root/contig[@name="{i}"]/flanks/match/@rightid"""
        left  = mosc.xpath(l)
        right = mosc.xpath(r)

        for index in range(len(left)):
            lstart = f"""/root/contig[@name="{i}"]/vectorhitleft[@id="{left[index]}"]/sstart/text()"""
            lend = f"""/root/contig[@name="{i}"]/vectorhitleft[@id="{left[index]}"]/send/text()"""
            acc = mosc.xpath(f"""/root/contig[@name="{i}"]/vectorhitleft[@id="{left[index]}"]/@seqid""")[0]
            if mosc.xpath(lstart)[0] > mosc.xpath(lend)[0]:
                Lrev = True
                lastLeft = int(mosc.xpath(lstart)[0])
            else:
                Lrev = False
                lastLeft = int(mosc.xpath(lend)[0])

            rstart = f"""/root/contig[@name="{i}"]/vectorhitright[@id="{right[index]}"]/sstart/text()"""
            rend = f"""/root/contig[@name="{i}"]/vectorhitright[@id="{right[index]}"]/send/text()"""

            if mosc.xpath(rstart)[0] > mosc.xpath(rend)[0]:
                Rrev = True
                firstRight = int(mosc.xpath(rend)[0])
            else:
                Rrev = False
                firstRight = int(mosc.xpath(rstart)[0])

            if firstRight < lastLeft:
                overlapSeq = str(getGen(acc, firstRight, lastLeft))
                overlapInfo = {"lastLeftPosition": lastLeft, "firstRightPosition": firstRight, "overlap": overlapSeq}
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

def main():
    path = os.path.join(".", "../../../../../Data2/specimens/Bangkok_Thailand_15.LIN210A1693/results/xml/Bangkok_Thailand_15.LIN210A1693_hits.xml")
    print(getOverlaps(path))


main()
