
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio import SeqIO
from pathlib import Path
import xml.etree.ElementTree as ET
from intervaltree import Interval, IntervalTree
import re
import pickle

def printFeature(f):
    print('Location: ', f.location)
    print('Operator: ', f.location_operator)
#    print('Strand: ', f.strand)
    print('Type: ' + f.type)
    print('Id: ' + f.id)
    print('Qualifiers: ' + str(f.qualifiers))
    
def populateFeature(feature, f):
    """Populate an XML feature element from a SeqFeature object."""
    
    ET.SubElement(feature, 'location').text = str(f.location)
    ET.SubElement(feature, 'type').text = f.type
    #ET.SubElement(feature, 'qualifiers').text = str(f.qualifiers)
    
    if 'gene' in f.qualifiers:
        ET.SubElement(feature, 'gene').text = f.qualifiers['gene'][0]
    else:
        ET.SubElement(feature, 'gene').text = 'None'
        
    if 'product' in f.qualifiers:
        ET.SubElement(feature, 'product').text = f.qualifiers['product'][0]
    else:
        ET.SubElement(feature, 'product').text = 'None'
        
    if 'unknown' not in f.id:
        ET.SubElement(feature, 'id').text = f.id
    elif 'protein_id' in f.qualifiers:
        ET.SubElement(feature, 'id').text = f.qualifiers['protein_id'][0]
    elif 'transcript_id' in f.qualifiers:
        ET.SubElement(feature, 'id').text = f.qualifiers['transcript_id'][0]
    elif 'gene' in f.qualifiers:
        ET.SubElement(feature, 'id').text = f.qualifiers['gene'][0]
    else:
        ET.SubElement(feature, 'id').text = 'Unknown'

def getIntrons(location):
    """Get a list of intron intervals from a CompoundLocation object representing exons."""
    
    prev = location.start - 1
    intervals = []
    for pos in sorted(location):
        if pos > prev + 1:
            intervals.append(Interval(prev+1, pos-1))
        prev = pos
    return intervals

def searchAA(iTrees, seqid, start, end, hitElement, counts):
    """Get features for a given interval and write to xml tree."""
    
    searchInterval = Interval(start, end)
    
    overlapIntervals = iTrees[seqid][start:end]
    for interval in overlapIntervals:
        f = interval.data                       # feature associated with interval
        if f.location_operator == 'join':       # mRNA with a CompoundLocation
            introns = getIntrons(f.location)
            featureType = None
            for intron in introns:
                if intron.contains_interval(searchInterval):
                    featureType = 'in_intron'                 # wholly contained in an intron
                    break
                elif intron.overlaps(searchInterval):
                    featureType = 'overlaps_intron'           # crosses an intron/exon boundary
                    break
                    
            if featureType is None:
                if start in f.location and end in f.location:
                    featureType = 'in_exon'                   # wholly contained in an exon
                else:
                    featureType = 'overlaps_end_exon'         # overlaps first or last exon
        else:
            featureType = f.type
                 
        feature = ET.SubElement(hitElement, 'feature', {'type': featureType})
        populateFeature(feature, f)
        counts[featureType] = counts.get(featureType, 0) + 1
    
            
def findFeatures(fileName, iTrees):
    """Find all features in an xml file and insert them."""
    
    tree = ET.parse(fileName)
    root = tree.getroot()
    for contig in root:
        counts = {}
        hitsLeft = contig.findall('vectorhitleft')
        for v in hitsLeft:
            start = int(v.find('sstart').text)
            end = int(v.find('send').text)
            seqid = v.attrib['seqid']
            searchAA(iTrees, seqid, start, end, v, counts)
            
        hitsRight = contig.findall('vectorhitright')
        for v in hitsRight:
            start = int(v.find('sstart').text)
            end = int(v.find('send').text)
            seqid = v.attrib['seqid']
            searchAA(iTrees, seqid, start, end, v, counts)
            
        featureSummary = ET.SubElement(contig, 'featuresummary')
        for ftype in counts:
            ET.SubElement(featureSummary, ftype).text = str(counts[ftype])
        
    tree.write(fileName[:-4] + '_features.xml')
    
def getGFFFeatures(fileName, iTrees):
    """Read repeat features from a GFF3 file and insert into dictionary of interval trees."""
    
    chromNumbers = {'1': 'NC_035107.1', '2': 'NC_035108.1', '3': 'NC_035109.1'}
    strands = {'+': +1, '-': -1, '.': 0}
   
    gff = open(fileName, 'r')
    
    for line in gff:
        if line[0] == '#':
            continue
        
        cols = line.strip().split('\t')
        seqid = cols[0]
        if seqid == 'MT':
            continue
        if seqid in chromNumbers:
            seqid = chromNumbers[seqid]
        else:
            seqid += '.1'
        type = cols[2]
        if type in ['chromosome', 'supercontig']:
            continue
        start = int(cols[3])
        end = int(cols[4])
        strand = strands[cols[6]]
        attrib = cols[8]
        
        f = SeqFeature()
        f.location = FeatureLocation(start, end + 1, strand = strand)
        f.type = type
        findName = re.findall(r'Name=(.*);.*', 'Name=Gypsy-58_AA-I;class=LTR/Gypsy')
        if len(findName) > 0:
            f.id = findName[0]
        else:
            f.id = 'Unknown'
        if seqid in iTrees:
            iTrees[seqid][start:end] = f
    
    gff.close()

def getGenBankFeatures(fileName, iTrees):
    """Read features from a multi-record flat GenBank file and insert into dictionary of interval trees."""

    for seqrecord in SeqIO.parse(fileName, 'genbank'):
        iTrees[seqrecord.id] = IntervalTree()
        altID = re.findall(r'NIGP01.*\.1', seqrecord.annotations['comment'])
        if len(altID) > 0:
            altID = altID[0]
            iTrees[altID] = iTrees[seqrecord.id]
        for f in seqrecord.features:
            if f.type != 'source':
                iTrees[seqrecord.id][f.location.start:f.location.end - 1] = f

def main():
    if Path('featuretrees.pickle').exists():
        print('Pickled feature trees found.  Reading...')
        picklediTreesDict = open('featuretrees.pickle', 'rb')
        iTrees = pickle.load(picklediTreesDict)
        picklediTreesDict.close()
    else:
        iTrees = {}
        print('Reading GenBank features...')
        # https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/002/204/515/GCF_002204515.2_AaegL5.0/
        getGenBankFeatures('/home/havill/data/aegypti/gb/GCF_002204515.2_AaegL5.0_genomic.gbff', iTrees)
    
        print('Reading repeat features...')
        getGFFFeatures('/home/havill/data/aegypti/gb/Aedes-aegypti-LVP_AGWG_REPEATFEATURES_AaegL5.gff3', iTrees)
    
        print('Writing feature trees to disk...')
        picklediTreesDict = open('featuretrees.pickle', 'wb')
        pickle.dump(iTrees, picklediTreesDict)
        picklediTreesDict.close()
    
    print('Annotating features...')
    findFeatures('Amacuzac_Mexico_combined/spades_cfav/Amacuzac_Mexico_combined_cfav_hits.xml', iTrees)
#    findFeatures('La_Lope-Gabon-10.LIN210A1646/spades_cfav/La_Lope-Gabon-10.LIN210A1646_viral_hits.xml', iTrees)
    
main()

# generating tree: 8 min; + write tree: 10 min
# reading pickle: 4:45 min
