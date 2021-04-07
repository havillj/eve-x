from Bio import Entrez
from lxml import etree as ET

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
        print(lineage.text)
        lineageParts = lineage.text.split('; ')
        if len(lineageParts) <= 6:    # 1-6
            return lineageParts[-1]
        elif len(lineageParts) <= 8:  # 7-8
            return lineageParts[6]
        else:
            return lineageParts[7]    # 9-
    else:
        return ''

#print(getFamily('NC_038287.1'))
#print(getFamily('NC_038512.1'))

def getFamsForResults():
    f = open('/home/havill/data/aegypti/results/HITS_all/results_all.txt', 'r')
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
    #print(fams)

accs = ['NC_040532.1', 
'NC_033701.1', 
'NC_040599.1', 
'NC_028241.1', 
'NC_038287.1', 
'NC_035132.1', 
'NC_038278.1', 
'NC_034542.1', 
'NC_034537.1', 
'NC_031268.1', 
'NC_033705.1', 
'NC_025342.1', 
'NC_039206.1', 
'NC_038285.1', 
'NC_038276.1', 
'NC_022580.1', 
'NC_025378.1', 
'NC_025384.1', 
'NC_038279.1', 
'NC_038284.1', 
'NC_034508.1', 
'NC_031236.1', 
'NC_031691.1', 
'NC_031305.1', 
'NC_038286.1', 
'NC_020803.1', 
'NC_034535.1', 
'NC_033103.1', 
'NC_033731.1', 
'NC_025394.1', 
'NC_028232.1', 
'NC_039202.1', 
'NC_031690.1', 
'NC_025364.1', 
'NC_034540.1', 
'NC_034538.1', 
'NC_034533.1', 
'NC_025340.1']

for acc in accs:
    fam = getFamily(acc)
    print(acc, fam)
