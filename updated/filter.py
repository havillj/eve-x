from lxml import etree
import os
from pathlib import Path
import sys
#from Bio import SeqIO, Entrez, AlignIO

#import util
parser = etree.XMLParser(remove_blank_text = True)
path = os.path.join(".", "../../../../../Data2/specimens/Bangkok_Thailand_15.LIN210A1693/results/xml/Bangkok_Thailand_15.LIN210A1693_hits.xml")
#path = os.path.join(".", "../../../Data2/specimens/Amacuzac-Mexico-13.LIN210A1627/results/xml/Amacuzac-Mexico-13.LIN210A1627_hits.xml")
#print(path)

mosc = etree.parse(path, parser).getroot()

xstring = f"""/root/contig[@besthit="True"][virushit/@stitle="Cell fusing agent virus"]/@name"""
contigs = mosc.xpath(xstring)
overlaps = {}
for i in contigs:
    contigOverlaps = []
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
            overlapInfo = {"accessionID": acc, "lastLeftPosition": lastLeft, "firstRightPosition": firstRight, "leftFlankReversed": Lrev, "rightFlankReversed": Rrev, "overlap": lastLeft - firstRight}
            contigOverlaps.append(overlapInfo)
    overlaps[i] = contigOverlaps

print(overlaps)
