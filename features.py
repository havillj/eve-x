import re
from Bio.SeqFeature import SeqFeature, FeatureLocation
from intervaltree import Interval, IntervalTree
from lxml import etree as ET
from pathlib import Path
import pickle

from util import *

iTrees = None

def getGFFFeatures(fileName):
    """Read repeat features from a GFF3 file and insert into dictionary of interval trees."""
    
    global iTrees
    
    chromNumbers = {'1': 'NC_035107.1', '2': 'NC_035108.1', '3': 'NC_035109.1', 'MT': 'NC_035159.1'}
    strands = {'+': +1, '-': -1, '.': 0}
   
    gff = open(fileName, 'r')
    
    for line in gff:
        if line[0] == '#':
            continue
        
        cols = line.strip().split('\t')
        
        seqid = cols[0]
        if seqid in chromNumbers:
            seqid = chromNumbers[seqid]
        else:
            seqid += '.1'
            
        if seqid not in iTrees:
            iTrees[seqid] = IntervalTree()
            
        type = cols[2]            
        attrib = cols[8]
        
        if type == 'supercontig':
            altID = re.findall(r'Alias=' + seqid + r',(.*)', attrib)  # in supercontig record
            if len(altID) > 0:
                altID = altID[0]
                iTrees[altID] = iTrees[seqid]
        elif type in ('gene', 'ncRNA_gene', 'pseudogene', 'exon', 'repeat_region'):
            start = int(cols[3])
            end = int(cols[4])
            strand = strands[cols[6]]
            f = SeqFeature()
            f.location = FeatureLocation(start, end + 1, strand = strand)
            f.type = type
            findName = re.findall(r'Name=(.*?);.*', attrib)
            if len(findName) > 0:
                f.id = findName[0]
            else:
                findName = re.findall(r'ID=(.*?);.*', attrib)  # try again
                if len(findName) > 0:
                    f.id = findName[0]
                else:
                    f.id = 'Unknown'
            iTrees[seqid][start:end+1] = f  # interval start <= x < end+1
    
    gff.close()
                
def makeIntervalTrees():
    global iTrees
    if iTrees is not None:
        return
        
    pickleName = str(FEATURE_TREES_PICKLED_FILENAME)
    if Path(pickleName).exists():
        writelog('Pickled feature trees found.  Reading...', True)
        picklediTreesDict = open(pickleName, 'rb')
        iTrees = pickle.load(picklediTreesDict)
        picklediTreesDict.close()
    else:
        iTrees = {}
        writelog('Reading base features...', True)
        # https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/002/204/515/GCF_002204515.2_AaegL5.0/
        getGFFFeatures(BASE_FEATURES_GFF3_FILENAME)
    
        writelog('Reading repeat features...', True)
        getGFFFeatures(REPEAT_FEATURES_GFF3_FILENAME)
    
        writelog('Writing feature trees to disk...', True)
        picklediTreesDict = open(pickleName, 'wb')
        pickle.dump(iTrees, picklediTreesDict)
        picklediTreesDict.close()
    
###############################################################################

def searchiTreesXML(seqid, start, end, hitElement, counts):
    """Get features for a given interval and write to xml tree."""
    
#    searchInterval = Interval(start, end+1)
    
    overlapIntervals = iTrees[seqid][start:end]
    allFeatureTypes = []
    for interval in overlapIntervals:
        f = interval.data                       # feature associated with interval
        featureType = f.type
        allFeatureTypes.append(featureType)
                 
        feature = ET.SubElement(hitElement, 'feature', {'type': featureType})
        ET.SubElement(feature, 'location').text = str(f.location)
        ET.SubElement(feature, 'type').text = f.type
        ET.SubElement(feature, 'id').text = f.id
        
    if 'gene' in allFeatureTypes:
        if 'exon' in allFeatureTypes:
            counts['exon'] = counts.get('exon', 0) + 1
        else:
            counts['intron'] = counts.get('intron', 0) + 1
    elif 'ncRNA_gene' in allFeatureTypes:
        if 'exon' in allFeatureTypes:
            counts['exon_ncRNA'] = counts.get('exon_ncRNA', 0) + 1
        else:
            counts['intron_ncRNA'] = counts.get('intron_ncRNA', 0) + 1
    elif 'pseudogene' in allFeatureTypes:
        if 'exon' in allFeatureTypes:
            counts['exon_pseudogene'] = counts.get('exon_pseudogene', 0) + 1
        else:
            counts['intron_pseudogene'] = counts.get('intron_pseudogene', 0) + 1
        
    for featureType in set(allFeatureTypes):    # only count one repeat_region per sequence
        if featureType not in ('gene', 'ncRNA_gene', 'pseudogene', 'exon'):
            counts[featureType] = counts.get(featureType, 0) + 1
    
    if len(overlapIntervals) == 0:
        counts['None'] += 1
            
def findFeaturesXML(fileName):
    """Find all features in an xml file and insert them.
    
       Parameters:
           fileName: absolute path of XML file containing hit results
           iTrees:   dict of interval trees containing AA features
           
       Return value: absolute path of XML file augmented with features
    """
    
    writelog('Annotating features...', True)
    
    parser = ET.XMLParser(remove_blank_text=True)  # lxml
    tree = ET.parse(fileName, parser)
    root = tree.getroot()
    for contig in root:
        counts = {'None': 0}
        hitsLeft = contig.findall('vectorhitleft')
        for v in hitsLeft:
            start = int(v.find('sstart').text)
            end = int(v.find('send').text)
            seqid = v.attrib['seqid']
            searchiTreesXML(seqid, min(start, end), max(start, end), v, counts)
            
        hitsRight = contig.findall('vectorhitright')
        for v in hitsRight:
            start = int(v.find('sstart').text)
            end = int(v.find('send').text)
            seqid = v.attrib['seqid']
            searchiTreesXML(seqid, min(start, end), max(start, end), v, counts)
            
        featureSummary = ET.SubElement(contig, 'featuresummary')  # number of flanking hits with each feature
        for ftype in counts:
            ET.SubElement(featureSummary, ftype).text = str(counts[ftype])
        
    xmlFeaturesFilename = fileName[:-4] + '_features.xml'
    tree.write(xmlFeaturesFilename, xml_declaration=True, pretty_print=True)
    
    return xmlFeaturesFilename
    
###############################################################################
    
def searchiTrees(seqid, start, end):
    """Get features for a given interval."""
    
    if start > end:
        start, end = end, start
        
    start -= FEATURE_SEARCH_DIST
    end += FEATURE_SEARCH_DIST
    
    if iTrees[seqid][start:end] is None:
        return []
    else:
        return [(iv.data.type, iv.data.id, str(iv.data.location)) for iv in iTrees[seqid][start:end]]
    