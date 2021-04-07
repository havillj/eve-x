import matplotlib.pyplot as pyplot
from pathlib import Path
from Bio import Entrez
from lxml import etree as ET
import os
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio import SeqIO
from dna_features_viewer import GraphicFeature, GraphicRecord, BiopythonTranslator
import math

ROOT = '/home/havill/data/aegypti/analyzed'
RESULTS_ROOT = '/home/havill/data2/results4/'
    
def plotLengths():
    VIRUS_NAMES = ['cfav', 'xincheng', 'anphevirus', 'horseradish', 'aus_anopheles'] # 'liao_ning'

    lengths = []
    evalues = []
    for dir in Path(ROOT).iterdir():
        for virus in VIRUS_NAMES:
            fn = dir / ('spades_cfav/blast_scaffolds_' + virus + '.csv')
            if fn.exists():
                hits = readCSV(str(fn), 0)
                for contig in hits:
                    for hit in hits[contig]:
                        if (len(hit[2]) < 100): # and (float(hit[6]) > 1e-10):
                            lengths.append(len(hit[2]))
                            evalues.append(float(hit[6]))
    pyplot.scatter(lengths, evalues, s = 3)
    pyplot.xlabel('Length')
    pyplot.ylabel('e-value')
    pyplot.show()
    
    pyplot.hist(evalues)
    pyplot.show()
    
##############################################################################

def checkSuccess():
    dirList = [dir for dir in Path(ROOT).iterdir() if 'combined' not in dir.name]
    dirList =  ['La_Lope_Gabon-_22.2-116710',
                'Skukusa_South_Africa_09.LIN210A1747',
                'Amacuzac_Mexico-_2.2-116616',
                'Amacuzac_Mexico-_21.2-116635',
                'Cairns-Australia-16.LIN210A2166',
                'Cairns-Australia-9.LIN210A2159',
                'El_Dorado_US_U_31',
                'Maricopa_County_AZ-_11.2-116601',
                'Maricopa_County_AZ-_4.2-116544',
                'Amacuzac_Mexico-_12.2-116626',
                'Maricopa_County_AZ-_14.2-116604',
                'Cebu_City_Philippines-_8.2-116720',
                'Cebu_City_Philippines-_17.2-116729',
                'Cebu_City_Philippines-_5.2-116717',
                'El_Dorado_US_U_28',
                'Maricopa_County_AZ-_1.2-116541',
                'Paqueta-Brazil-5.LIN210A2131',
                'La_Lope_Gabon-_21.2-116709',
                'Paqueta-Brazil-9.LIN210A2135',
                'Amacuzac_Mexico-_9.2-116623']
    for dir in dirList:
        dir = Path(ROOT) / dir
        if dir.is_dir():
            spadesDir = dir / 'spades_all'
            scaffoldsFilename = dir / 'spades_all/scaffolds.fasta'
            print(dir.name, end = '\t')
            if spadesDir.exists():
                if scaffoldsFilename.exists():
                    print('assembly OK')
                else:
                    print('assembly failed')
            else:
                print('no spades dir')

#checkSuccess()
#exit()

def rerun():
    dirList = [dir for dir in Path(ROOT).iterdir() if 'combined' not in dir.name]
#    dirList = [Path('/home/havill/data/aegypti/analyzed/Bangkok_Thailand_09.LIN210A1687')]
    for dir in dirList:
        if dir.is_dir():
            blastResultsFilename = dir / 'spades_all/blast_scaffolds_all.csv'
            scaffoldsFilename = dir / 'spades_all/scaffolds.fasta'
            print(dir.name, end = ': ')
            if blastResultsFilename.exists():
                print('re-blasting contigs')
                scaffoldRecords = SeqIO.index(str(scaffoldsFilename), 'fasta')
                hits = readCSV(str(blastResultsFilename), 0)
                for contig in hits:
                    SeqIO.write(scaffoldRecords[contig], '/tmp/query_temp.fasta', 'fasta')
                    os.system('blastn -query /tmp/query_temp.fasta -db virusdb -task blastn -evalue 1e-10' +
                              ' -outfmt "10 qseqid qstart qend qseq sstart send sseq evalue bitscore sseqid stitle pident" >> /tmp/blast_scaffolds_all_temp.csv')
                os.system('mv ' + str(blastResultsFilename) + ' ' + str(blastResultsFilename)[:-4] + '_backup.csv')
                os.system('mv /tmp/blast_scaffolds_all_temp.csv ' + str(blastResultsFilename))
                
#                 xmlFilename = getHits(dir, 'all', 100)
#                 drawXML(xmlFilename, str(Path(xmlFilename).parent / 'png'), False)
            else:
                print('skipping')
                
# rerun()
# exit()

def redoHits():
#    dirList = [dir for dir in Path(ROOT).iterdir() if 'combined' not in dir.name]
#    dirList = [Path(ROOT) / 'El_Dorado_US_U_28']
#     dirNames = ['Cairns-Australia-9.LIN210A2159',
#                 'Maricopa_County_Arizona-USA-23.LIN210A1617',
#                 'Cuanda_Angola_10.LIN210A1728',
#                 'La_Lope-Gabon-24.LIN210A1658',
#                 'Paqueta-Brazil-8.LIN210A2134',
#                 'El_Dorado_US_U_31',
#                 'Cairns-Australia-12.LIN210A2162',
#                 'Tahiti_FrenchPolynesia_02.LIN210A1700',
#                 'Cuanda_Angola_08.LIN210A1726',
#                 'featuretrees.pickle.bak',
#                 'Paqueta-Brazil-11.LIN210A2137',
#                 'Amacuzac_Mexico-_2.2-116616',
#                 'Amacuzac-Mexico-16.LIN210A1629',
#                 'La_Lope-Gabon-9.LIN210A1645',
#                 'HoChiMin_Vietnam_12.LIN210A1670',
#                 'La_Lope-Gabon-12.LIN210A1648',
#                 'Paqueta-Brazil-13.LIN210A2139',
#                 'Maricopa_County_AZ-_11.2-116601',
#                 'Cairns-Australia-21.LIN210A2171',
#                 'Amacuzac-Mexico-20.LIN210A1633',
#                 'Cebu_City_Philippines-_16.2-116728',
#                 'La_Lope-Gabon-5.LIN210A1641',
#                 'Tahiti_FrenchPolynesia_08.LIN210A1706',
#                 'Tahiti_FrenchPolynesia_09.LIN210A1707',
#                 'Bangkok_Thailand_20.LIN210A1698',
#                 'Cuanda_Angola_20.LIN210A1738',
#                 'Maricopa_County_AZ-_4.2-116544',
#                 'Cebu_City_Philippines-_10.2-116722',
#                 'Cairns-Australia-19.LIN210A2169',
#                 'Amacuzac_Mexico-_21.2-116635',
#                 'HoChiMin_Vietnam_03.LIN210A1661',
#                 'Maricopa_County_AZ-_1.2-116541',
#                 'Bangkok_Thailand_17.LIN210A1695',
#                 'Cebu_City_Philippines-_14.2-116726',
#                 'Bangkok_Thailand_07.LIN210A1685',
#                 'Cebu_City_Philippines-_9.2-116721',
#                 'Bangkok_Thailand_03.LIN210A1681',
#                 'La_Lope-Gabon-14.LIN210A1650',
#                 'Amacuzac-Mexico-13.LIN210A1627',
#                 'Cuanda_Angola_03.LIN210A1721',
#                 'Skukusa_South_Africa_08.LIN210A1746',
#                 'Maricopa_County_AZ-_18.2-116608',
#                 'HoChiMin_Vietnam_01.LIN210A1659',
#                 'Cairns-Australia-15.LIN210A2165',
#                 'La_Lope-Gabon-2.LIN210A1638',
#                 'Skukusa_South_Africa_12.LIN210A1750',
#                 'Tahiti_FrenchPolynesia_18.LIN210A1716',
#                 'Cebu_City_Philippines-_11.2-116723',
#                 'Maricopa_County_AZ-_12.2-116602',
#                 'Tahiti_FrenchPolynesia_04.LIN210A1702',
#                 'El_Dorado_Argentina_U_3.LIN210A2706',
#                 'Cairns-Australia-1.LIN210A2151',
#                 'HoChiMin_Vietnam_06.LIN210A1664',
#                 'Amacuzac_Mexico-_12.2-116626',
#                 'La_Lope-Gabon-18.LIN210A1654',
#                 'Cairns-Australia-22.LIN210A2172',
#                 'La_Lope-Gabon-3.LIN210A1639',
#                 'Cairns-Australia-7.LIN210A2157',
#                 'Cebu_City_Philippines-_20.2-116732',
#                 'Paqueta-Brazil-19.LIN210A2145',
#                 'Skukusa_South_Africa_13.LIN210A1751',
#                 'Maricopa_County_AZ-_8.2-116548',
#                 'Skukusa_South_Africa_14.LIN210A1752',
#                 'Tahiti_FrenchPolynesia_14.LIN210A1712',
#                 'HoChiMin_Vietnam_09.LIN210A1667',
#                 'Amacuzac-Mexico-6.LIN210A1622',
#                 'Cebu_City_Philippines-_4.2-116716',
#                 'Maricopa_County_Arizona-USA-10.LIN210A1614',
#                 'Cuanda_Angola_14.LIN210A1732',
#                 'Amacuzac-Mexico-24.LIN210A1636',
#                 'Amacuzac-Mexico-3.LIN210A1619',
#                 'Cairns-Australia-17.LIN210A2167',
#                 'Cebu_City_Philippines-_21.2-116733',
#                 'Amacuzac-Mexico-1.LIN210A1618',
#                 'Cairns-Australia-24.LIN210A2174',
#                 'Maricopa_County_Arizona-USA-5.LIN210A1612',
#                 'Tahiti_FrenchPolynesia_05.LIN210A1703',
#                 'El_Dorado_Argentina_U_11.LIN210A2770',
#                 'Maricopa_County_Arizona-USA-7.LIN210A1613',
#                 'Paqueta-Brazil-10.LIN210A2136',
#                 'Cairns-Australia-18.LIN210A2168',
#                 'Tahiti_FrenchPolynesia_10.LIN210A1708',
#                 'Bangkok_Thailand_04.LIN210A1682',
#                 'Maricopa_County_AZ-_14.2-116604',
#                 'La_Lope-Gabon-20.LIN210A1655',
#                 'Cairns-Australia-20.LIN210A2170',
#                 'La_Lope_Gabon-_19.2-116707',
#                 'La_Lope-Gabon-10.LIN210A1646',
#                 'Maricopa_County_Arizona-USA-2.LIN210A1611',
#                 'Maricopa_County_AZ-_13.2-116603',
#                 'Tahiti_FrenchPolynesia_03.LIN210A1701',
#                 'Paqueta-Brazil-15.LIN210A2141',
#                 'Amacuzac_Mexico-_15.2-116629',
#                 'Cebu_City_Philippines-_19.2-116731',
#                 'Paqueta-Brazil-5.LIN210A2131',
#                 'Amacuzac-Mexico-17.LIN210A1630',
#                 'Cebu_City_Philippines-_8.2-116720',
#                 'Skukusa_South_Africa_15.LIN210A1753',
#                 'Maricopa_County_AZ-_16.2-116606',
#                 'Cairns-Australia-10.LIN210A2160',
#                 'HoChiMin_Vietnam_20.LIN210A1678',
#                 'HoChiMin_Vietnam_13.LIN210A1671',
#                 'Maricopa_County_Arizona-USA-22.LIN210A1616',
#                 'Cuanda_Angola_18.LIN210A1736',
#                 'Amacuzac-Mexico-11.LIN210A1626',
#                 'La_Lope_Gabon-_21.2-116709',
#                 'La_Lope-Gabon-16.LIN210A1652',
#                 'Paqueta-Brazil-9.LIN210A2135',
#                 'Bangkok_Thailand_18.LIN210A1696',
#                 'Cebu_City_Philippines-_17.2-116729',
#                 'Skukusa_South_Africa_03.LIN210A1741',
#                 'Cuanda_Angola_11.LIN210A1729',
#                 'Cairns-Australia-13.LIN210A2163',
#                 'La_Lope-Gabon-7.LIN210A1643',
#                 'Cairns-Australia-2.LIN210A2152',
#                 'Cebu_City_Philippines-_23.2-116735',
#                 'Cebu_City_Philippines-_6.2-116718',
#                 'Paqueta-Brazil-3.LIN210A2129',
#                 'Amacuzac-Mexico-7.LIN210A1623',
#                 'Maricopa_County_AZ-_24.2-116614',
#                 'HoChiMin_Vietnam_08.LIN210A1666',
#                 'Paqueta-Brazil-12.LIN210A2138',
#                 'HoChiMin_Vietnam_05.LIN210A1663',
#                 'La_Lope-Gabon-17.LIN210A1653',
#                 'Maricopa_County_AZ-_6.2-116546',
#                 'Amacuzac-Mexico-10.LIN210A1625',
#                 'Bangkok_Thailand_14.LIN210A1692',
#                 'Skukusa_South_Africa_07.LIN210A1745',
#                 'HoChiMin_Vietnam_02.LIN210A1660',
#                 'Tahiti_FrenchPolynesia_16.LIN210A1714',
#                 'Amacuzac_Mexico-_9.2-116623',
#                 'Bangkok_Thailand_02.LIN210A1680',
#                 'Bangkok_Thailand_13.LIN210A1691',
#                 'Bangkok_Thailand_11.LIN210A1689',
#                 'featuretrees.pickle',
#                 'Bangkok_Thailand_16.LIN210A1694',
#                 'HoChiMin_Vietnam_11.LIN210A1669',
#                 'Cebu_City_Philippines-_24.2-116736',
#                 'Maricopa_County_AZ-_3.2-116543',
#                 'Cebu_City_Philippines-_5.2-116717',
#                 'Cebu_City_Philippines-_3.2-116715',
#                 'Cebu_City_Philippines-_22.2-116734',
#                 'La_Lope-Gabon-13.LIN210A1649',
#                 'Cebu_City_Philippines-_13.2-116725',
#                 'Cairns-Australia-5.LIN210A2155',
#                 'Bangkok_Thailand_10.LIN210A1688',
#                 'Amacuzac-Mexico-14.LIN210A1628',
#                 'Tahiti_FrenchPolynesia_12.LIN210A1710',
#                 'HoChiMin_Vietnam_15.LIN210A1673',
#                 'Bangkok_Thailand_08.LIN210A1686',
#                 'HoChiMin_Vietnam_14.LIN210A1672',
#                 'Cairns-Australia-23.LIN210A2173',
#                 'Tahiti_FrenchPolynesia_06.LIN210A1704',
#                 'Cuanda_Angola_12.LIN210A1730',
#                 'Cuanda_Angola_16.LIN210A1734',
#                 'Paqueta-Brazil-4.LIN210A2130',
#                 'Amacuzac-Mexico-23.LIN210A1635',
#                 'Paqueta-Brazil-1.LIN210A2127',
#                 'Amacuzac-Mexico-19.LIN210A1632',
#                 'Paqueta-Brazil-18.LIN210A2144',
#                 'Tahiti_FrenchPolynesia_15.LIN210A1713',
#                 'Cuanda_Angola_19.LIN210A1737',
#                 'Bangkok_Thailand_01.LIN210A1679',
#                 'Cuanda_Angola_04.LIN210A1722',
#                 'Tahiti_FrenchPolynesia_20.LIN210A1718',
#                 'Skukusa_South_Africa_05.LIN210A1743',
#                 'HoChiMin_Vietnam_04.LIN210A1662',
#                 'La_Lope-Gabon-8.LIN210A1644',
#                 'Bangkok_Thailand_19.LIN210A1697',
#                 'Cebu_City_Philippines-_7.2-116719',
#                 'Cebu_City_Philippines-_12.2-116724',
#                 'Amacuzac-Mexico-22.LIN210A1634',
#                 'Paqueta-Brazil-24.LIN210A2150',
#                 'La_Lope-Gabon-6.LIN210A1642',
#                 'Paqueta-Brazil-7.LIN210A2133',
#                 'Bangkok_Thailand_06.LIN210A1684',
#                 'Maricopa_County_AZ-_15.2-116605',
#                 'Skukusa_South_Africa_11.LIN210A1749',
#                 'Cairns-Australia-11.LIN210A2161',
#                 'Cuanda_Angola_02.LIN210A1720',
#                 'Paqueta-Brazil-14.LIN210A2140',
#                 'Maricopa_County_AZ-_9.2-116549',
#                 'Bangkok_Thailand_05.LIN210A1683',
#                 'La_Lope-Gabon-23.LIN210A1656',
#                 'HoChiMin_Vietnam_07.LIN210A1665']

    dirNames = ['La_Lope-Gabon-11.LIN210A1647',
                'Maricopa_County_Arizona-USA-20.LIN210A1615',
                'Maricopa_County_AZ-_19.2-116609',
                'El_Dorado_Argentina_U_2.LIN210A2698',
                'Tahiti_FrenchPolynesia_13.LIN210A1711',
                'Tahiti_FrenchPolynesia_19.LIN210A1717',
                'Cairns-Australia-8.LIN210A2158',
                'La_Lope_Gabon-_22.2-116710',
                'Cebu_City_Philippines-_2.2-116714',
                'Amacuzac-Mexico-8.LIN210A1624',
                'Cuanda_Angola_01.LIN210A1719',
                'HoChiMin_Vietnam_17.LIN210A1675',
                'La_Lope-Gabon-1.LIN210A1637',
                'Cuanda_Angola_06.LIN210A1724',
                'Cuanda_Angola_05.LIN210A1723',
                'Paqueta-Brazil-23.LIN210A2149',
                'Bangkok_Thailand_09.LIN210A1687',
                'Amacuzac-Mexico-4.LIN210A1620',
                'La_Lope-Gabon-15.LIN210A1651',
                'Paqueta-Brazil-2.LIN210A2128',
                'Amacuzac-Mexico-5.LIN210A1621',
                'Cuanda_Angola_09.LIN210A1727',
                'HoChiMin_Vietnam_18.LIN210A1676',
                'Cebu_City_Philippines-_18.2-116730',
                'Cairns-Australia-16.LIN210A2166',
                'El_Dorado_Argentina_U_8.LIN210A2746',
                'Cuanda_Angola_17.LIN210A1735',
                'HoChiMin_Vietnam_16.LIN210A1674',
                'Skukusa_South_Africa_02.LIN210A1740',
                'Amacuzac-Mexico-18.LIN210A1631',
                'Paqueta-Brazil-6.LIN210A2132',
                'Bangkok_Thailand_12.LIN210A1690',
                'Skukusa_South_Africa_09.LIN210A1747',
                'Skukusa_South_Africa_04.LIN210A1742',
                'HoChiMin_Vietnam_10.LIN210A1668',
                'Maricopa_County_AZ-_21.2-116611',
                'Cebu_City_Philippines-_1.2-116713',
                'Tahiti_FrenchPolynesia_01.LIN210A1699',
                'HoChiMin_Vietnam_19.LIN210A1677',
                'Skukusa_South_Africa_06.LIN210A1744',
                'Tahiti_FrenchPolynesia_17.LIN210A1715',
                'Skukusa_South_Africa_10.LIN210A1748',
                'Tahiti_FrenchPolynesia_11.LIN210A1709',
                'Cebu_City_Philippines-_15.2-116727',
                'Tahiti_FrenchPolynesia_07.LIN210A1705',
                'El_Dorado_US_U_28',
                'Bangkok_Thailand_15.LIN210A1693',
                'La_Lope-Gabon-4.LIN210A1640']

    for dir in dirNames:
        dir = Path(ROOT) / dir
        if dir.is_dir():
            print(dir.name)
            xmlFilename = getHits(dir, 'all', 100)  #3 = get hits
            if xmlFilename is None:
                print('*** ' + str(Path(dir).name) + ': assembly failed; no results found. ***')
            else:
                drawXML(xmlFilename, str(Path(xmlFilename).parent / 'png'), False)

def getName(accessionID):
	"""
	Query NCBI to get name corresponding to an accession number.
	"""
	
	Entrez.email = 'havill@denison.edu'
	try:
		handle = Entrez.efetch(db = 'nucleotide', id = accessionID, rettype = 'gb', retmode = 'text')
		result = handle.read().split('\n')
		for line in result:
			if 'ORGANISM' in line:
				return ' '.join(line.split()[1:])
	except:
		print('Exception:', sys.exc_info()[1])
		return ''
		
def getFamily(accessionID):
    """
    Query NCBI to get family of virus corresponding to an accession number.
    """
    
    Entrez.email = 'havill@denison.edu'
    try:
        handle = Entrez.efetch(db = 'nucleotide', id = accessionID, rettype = 'native', retmode = 'xml')
        result = handle.read()
    except:
        print('Exception:', sys.exc_info()[1])
        return ''
        
    root = ET.fromstring(result)
    lineage = root.find('.//OrgName_lineage')

    if lineage is not None:
        lineageParts = lineage.text.split('; ')
        if len(lineageParts) <= 6:    # 1-6
            return lineageParts[-1]
        elif len(lineageParts) <= 8:  # 7-8
            return lineageParts[6]
        else:
            return lineageParts[7]    # 9-
    else:
        return ''
        
def getFamsForResults():
    f = open('/home/havill/data/aegypti/results/HITS_all/results_all.tsv', 'r')
    header = f.readline()
    f.close()
    accs = header.strip().split('\t')
    print(len(accs))
    fout = open('/home/havill/fams.txt', 'w')
    for acc in accs:
        print(acc)
        fam = getFamily(acc)
        print(fam)
        fout.write('\t' + fam)
    fout.close()
    
def readFamFile(fileName):
    try:
        famFile = open(fileName, 'r')
    except FileNotFoundError:
        return {}
    fams = {}
    for line in famFile:
        cols = line.rstrip().split(',')
        fams[cols[0]] = cols[1]
    famFile.close()
    return fams
    
def getFamsForBlastAA():
    f = open('/Volumes/Data2/results/blast_aa_virus.csv', 'r')
    fout = open('/Volumes/Data2/results/blast_aa_virus_families.csv', 'w')
    fams = readFamFile('/Volumes/Data2/results/families.csv')
    famFile = open('/Volumes/Data2/results/families.csv', 'a')
    for line in f:
        cols = line.rstrip().split(',')
        acc = cols[9].rstrip('|')
        if '|' in acc:
            acc = acc.split('|')[1]
        if acc in fams:
            fam = fams[acc]
        else:
            fam = getFamily(acc)
            fams[acc] = fam
            famFile.write(acc + ',' + fam + '\n')
        fout.write(line.rstrip() + ',' + fam + '\n')
    f.close()
    fout.close()
    famFile.close()
    
#getFamsForBlastAA()

def getFamsForBlastV(csvFileName):
    f = open(csvFileName, 'r')
#    fout = open('/Volumes/Data2/results/blast_aa_virus_families.csv', 'w')
    fams = readFamFile('/Volumes/Data2/results/families.csv')
    famFile = open('/Volumes/Data2/results/families.csv', 'a')
    foundFams = []
    for line in f:
        cols = line.rstrip().split(',')
        acc = cols[9]
#        print(acc, end = ' -> ')
        if '|' in acc:
            acc = acc.split('|')[1]
#        print(acc)
        if acc in fams:
            fam = fams[acc]
        else:
#            print('Querying ' + acc)
            fam = getFamily(acc)
            fams[acc] = fam
            famFile.write(acc + ',' + fam + '\n')
#        fout.write(line.rstrip() + ',' + fam + '\n')
        foundFams.append(fam)
        if fam == 'Coronaviridae':
            print(acc, cols[0])
    f.close()
#    fout.close()
    famFile.close()
    
    foundFams = list(set(foundFams))
    foundFams.sort()
    for fam in foundFams:
        print(fam)
    
    return foundFams
    
#getFamsForBlastV('/Volumes/Data/aegypti/analyzed/Amacuzac-Mexico-10.LIN210A1625/spades_all/blast_scaffolds_all.csv')
#getFamsForBlastV('/Volumes/Data/aegypti/analyzed/Bangkok_Thailand_10.LIN210A1688/spades_all/blast_scaffolds_all.csv')
#getFamsForBlastV('/Volumes/Data/aegypti/analyzed/La_Lope-Gabon-10.LIN210A1646/spades_all/blast_scaffolds_all.csv')
#getFamsForBlastV('/Volumes/Data/aegypti/analyzed/Tahiti_FrenchPolynesia_10.LIN210A1708/spades_all/blast_scaffolds_all.csv')


def examineContigs():
    fam = {}
    dir = Path('.').resolve()
    files = [file for file in dir.iterdir() if ('_hits_features.xml' in file.name) or (('_hits.xml' in file.name) and not Path(str(file)[:-4] + '_features.xml').exists())]
    files.sort()
    for file in files:
        index = file.name.index('_all')
        specimen = file.name[:index]
        hits = {}
        tree = ET.parse(str(file))
        root = tree.getroot()
        for contig in root:   # one per virus
            name = contig.attrib['name']
            if name.count('_') == 6:
                name = '_'.join(name.split('_')[:-1])
            if name not in hits:
                hits[name] = []
            v = contig.find('virushit')
#            seqid, stitle = v.attrib['seqid'].strip('ref|'), v.attrib['stitle'].lstrip('|')  # should be same for all v
            seqid = v.attrib['seqid']            # should be same for all v
            if '|' in seqid:
                seqid = seqid.split('|')[1]
            stitle = v.attrib['stitle']
            if '|' in stitle:
                stitle = stitle.lstrip('|')
            if seqid not in fam:
                fam[seqid] = getFamily(seqid)
            length = 0
            bitscore = 0
            pident = 0
            evalue = float('inf')
            sstart1 = 100000
            sstart2 = 0
            for v in contig.findall('virushit'):
                evalue = min(evalue, float(v.find('evalue').text))
                bitscore = max(bitscore, float(v.find('bitscore').text))
                length += len(v.find('qseq').text)
                pident = max(pident, float(v.find('pident').text))
                sstart1 = min(sstart1, int(v.find('sstart').text), int(v.find('send').text))
                sstart2 = max(sstart2, int(v.find('sstart').text), int(v.find('send').text))
            hits[name].append((seqid, stitle, pident, length, evalue, bitscore, sstart1, sstart2))
        for name in hits:
            hits[name].sort()
        names = list(hits.keys())
        names.sort(key = lambda n: len(hits[n]))
        print('\n' + '*' * 40 + '\n' + specimen + '\n' + '*' * 40)
        for name in names:
            print('\n' + name)
            maxPident = 0
            maxLength = 0
            for t in hits[name]:
                maxPident = max(maxPident, t[2])
                maxLength = max(maxLength, t[3])
            print('   {0:<11}  {1:<40}  {2:<11}  {3:<25}  {4:<6}   {5:>4}  {6:>9}  {7:>5}'.format('ACCESSION', 'NAME', 'MIN-MAX POS', 'FAMILY', 'PIDENT', 'LEN', 'E-VALUE', 'BITS'))
            for t in hits[name]:
                t = t + (fam[t[0]],)
                print('   {0}  {1:<40}  {6:>5}-{7:<5}  {8:<25} '.format(*(t)), end = '')
                if t[2] == maxPident:
                    print('*{2:6.3f}*'.format(*(t)), end = '')
                else:
                    print(' {2:6.3f} '.format(*(t)), end = '')
                if t[3] == maxLength:
                    print(' *{3:*>4}*'.format(*(t)), end = '')
                else:
                    print('  {3:>4} '.format(*(t)), end = '')
                print(' {4:>9}  {5:5.1f}'.format(*(t)))

#examineContigs()


def complexity(dna, W = 0):
    """https://doi.org/10.1016/S0097-8485(99)00007-8"""
    
    dna = dna.upper()
    dna = dna.replace('-', '')
    
    if W == 0:
        W = int(math.log(len(dna), 4)) + 2
    
    C = 1
    for k in range(1, W + 1):
        kmers = set()
        for i in range(0, len(dna) - k + 1):
            kmers.add(dna[i:i+k])
        C *= (len(kmers) / (4 ** k))
    return C
    
#print(complexity('TCGCTGACCCATTATACAAAAGGTACGCAGTCACAGAACAAGTCTGCTCCCACTGTTTGTATGCATGCGGTTTCAGGATCTATTTCACTCCCCTCCCGGGGTTCTTT', 3))
    
def getComplexityForBlastAA():
    f = open('/Volumes/Data2/results/blast_aa_virus.csv', 'r')
    for line in f:
        cols = line.rstrip().split(',')
        seq1 = cols[3]
        print(str(complexity(seq1, 3)))
    f.close()
    
#getComplexityForBlastAA()

def evalueDistribution(fileName):
    tree = ET.parse(fileName)
    root = tree.getroot()
    for contig in root:   # one per virus
        v = contig.find('virushit')
#        seqid, stitle = v.attrib['seqid'].strip('ref|'), v.attrib['stitle'].lstrip('|')  # should be same for all v
        seqid = v.attrib['seqid']            # should be same for all v
        if '|' in seqid:
            seqid = seqid.split('|')[1]
        stitle = v.attrib['stitle']
        if '|' in stitle:
            stitle = stitle.lstrip('|')
        evalues = []
        lengths = []
        for v in contig.findall('vectorhitleft') + contig.findall('vectorhitright'):
            evalue = float(v.find('evalue').text)
            bitscore = float(v.find('bitscore').text)
            length = len(v.find('qseq').text)
            evalues.append(evalue)
            lengths.append(length)
        pyplot.hist(evalues)
        pyplot.show()
        
#evalueDistribution('/Volumes/Data2/results/specimens/Bangkok_Thailand_10.LIN210A1688/xml/Bangkok_Thailand_10.LIN210A1688_all_hits_features.xml')

#print(complexity('acacaccccca'))

def getDuplicateSpecies(fileName):
    allSpecies = {}
    for r in SeqIO.parse(fileName, 'fasta'):
        parts = r.description.split(' |')
        acc = parts[0]
        species = parts[1]
        if species not in allSpecies:
            allSpecies[species] = []
        allSpecies[species].append(acc)
        
    speciesNames = list(allSpecies.keys())
    speciesNames.sort(key = lambda s: len(allSpecies[s]), reverse = True)
    for species in speciesNames:
        print(species, len(allSpecies[species]))
        
#getDuplicateSpecies('/home/havill/data/blastdb/source/newdb/newdb_no_retro_or_unverified.fasta')

def getFamily(accessionID):
    """
    Query NCBI to get family of virus corresponding to an accession number.
    """
    
    Entrez.email = 'havill@denison.edu'
    try:
        handle = Entrez.efetch(db = 'nucleotide', id = accessionID, rettype = 'native', retmode = 'xml')
        result = handle.read()
    except:
        print('Exception:', sys.exc_info()[1])
        return ''
        
    root = ET.fromstring(result)
    lineage = root.find('.//OrgName_lineage')

    if lineage is not None:
        lineageParts = lineage.text.split('; ')
        if len(lineageParts) <= 6:    # 1-6
            return lineageParts[-1]
        elif len(lineageParts) <= 8:  # 7-8
            return lineageParts[6]
        else:
            return lineageParts[7]    # 9-
    else:
        return ''
    
FAMILY_CSV = '/home/havill/data2/results/families.csv'

def readFamFile():
    try:
        famFile = open(FAMILY_CSV, 'r')
    except FileNotFoundError:
        return {}
    fams = {}
    for line in famFile:
        cols = line.rstrip().split(',')
        fams[cols[0]] = cols[1]
    famFile.close()
    return fams
    
def addFamily(acc, name):
    famFile = open(FAMILY_CSV, 'a')
    famFile.write(acc + ',' + name + '\n')
    famFile.close()

def getSpecies(fileName):
    allSpecies = {}
    count = 0
    for r in SeqIO.parse(fileName, 'fasta'):
        count += 1
        parts = r.description.split(' |')
        acc = parts[0]
        species = parts[1]
        if species not in allSpecies:
            allSpecies[species] = []
        allSpecies[species].append(acc)
        
    speciesNames = list(allSpecies.keys())
    speciesNames.sort()
    fams = readFamFile()
    for species in speciesNames:
        seqid = allSpecies[species][0]
        if seqid not in fams:
            fam = getFamily(seqid)
            fams[seqid] = fam
            addFamily(seqid, fam)
        print(species + '\t' + fams[seqid] + '\t' + str(len(allSpecies[species])) + '\t' + '\t'.join(allSpecies[species]))
    print(count)

#getSpecies('/home/havill/data/blastdb/source/newdb/newdb_no_retro_or_unverified.fasta')


