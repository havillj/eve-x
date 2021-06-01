from lxml import etree
import os
from pathlib import Path
import sys
from Bio import SeqIO, Entrez, AlignIO, Seq
import glob
from util import *
import csv

###############################################################################

def specimenFiles():
    """
    Retuns a dictionary mapping a specimen region and num to the specimen name

    Parameters:
        None

    Return value: a dictionary mapping a specimen region and num to the specimen name
    """
    xmlFiles = []
    specimensPath = Path("/Volumes/Data2/specimens")
    for specimenDir in specimensPath.iterdir():
        for file in (specimenDir / 'results/xml').iterdir():
            if ('_hits_features.xml' in file.name) or (('_hits.xml' in file.name) and not Path(str(file)[:-4] + '_features.xml').exists()):
                xmlFiles.append(file.resolve())
    xmlFiles.sort()
    specimenList = []
    label2Specimen = {}
    for file in xmlFiles:
        specimen = file.name.split('_hits')[0]
        region, num = getSpecimenLabel(specimen)
        region = region.replace(' ', '-')
        label2Specimen[(region, num)] = specimen

    return label2Specimen

###############################################################################

def getPath(chromosome):
    """
    Retuns a list of file paths of specimens with possible overlapping regions

    Parameters:
        chromosome: the number of the chromosome to use

    Return value: a list of file paths
    """

    tsvPath = f"/Volumes/Data2/results/viruses/sequences/Flaviviridae/insertpositions_NC_001564.2_chr{chromosome}.tsv"
    filepaths = []
    file = open(tsvPath, 'r')
    next(file)
    for line in file:
        row = line.strip().split("\t")
        if line[0] not in ["1","2","3","4","5","6"]:    #break if line is empty
            break

        if row[4] == "" and  "(" in row[1]:             #exclude exons and only use flanking regions
            if "(" in row[19]:
                spec = row[19].split(" (")[0]

            if "French-Polynesia" or "South-Africa" in spec:
                name = spec.split("_")[0]
                num = int(spec.split("_")[1])

            else:
                spec = spec.replace("_","-")
                name = spec.split("-")[0]
                num = int(spec.split("-")[1])

            specimensL =  specimenFiles()
            specP = specimensL[(name,num)]

            filepaths.append(f"/Volumes/Data2/specimens/{specP}/results/xml/{specP}_hits.xml")

    file.close()

    filepaths = list(dict.fromkeys(filepaths))

    return filepaths

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
                overlapInfo = {"lastLeftPosition": lastLeft, "firstRightPosition": firstRight}#, "overlap": overlapSeq}
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

def main():
    paths = getPath(2)

    with open("overlaps.tsv", "w") as outfile:
        tsv_writer = csv.writer(outfile, delimiter = "\t")
        tsv_writer.writerow(["Node", "Accession ID", "End Position", "Start Postion", "Left Overlap", "Right Overlap"])
        for file in paths:
            overlapData = getOverlaps(file)
            overlapL = []
            for i in overlapData:
                #print(i)
                for contig in overlapData[i]:
                    sub = list(overlapData[i][contig].values())
                    overlapL = [i, contig] + sub
                    #print(overlapL)
                    tsv_writer.writerow(overlapL)

main()
