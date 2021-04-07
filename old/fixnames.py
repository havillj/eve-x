from Bio import Entrez
import sys
import os
from pathlib import Path
              
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

def fixNames(fn):
	f = open(fn, 'r')
	newfn = fn[:-4] + '_fixed.csv'
	new = open(newfn, 'w')
	
	header = f.readline()
	new.write(header)
	
	names = {}
	fixesNeeded = False
	
	for line in f:
		cols = line.split(',')
		if cols[3] == '':
			fixesNeeded = True
			acc = cols[2]
			if '.' in acc:
				acc = acc.split('.')[0]
			if acc in names:
				virusName = names[acc]
			else:
				virusName = getName(acc)
				if virusName == '':
					f.close()
					new.close()
					os.rename(newfn, newfn[:-4] + '_partial.csv')
					print('NCBI query failed.')
					return
				else:
					names[acc] = virusName
			cols[3] = virusName
			line = ','.join(cols)
		new.write(line)
		
	f.close()
	new.close()
	
	if not fixesNeeded:
		os.remove(newfn)
		print(fn + ': no missing names')
	else:
		print(fn + ': missing names saved in ' + newfn)
	
def main():
	if len(sys.argv) == 1:
		cwd = Path('.')
		csvFiles = [f for f in cwd.iterdir() if not f.is_dir() and f.name[-4:] == '.csv']
		for f in csvFiles:
			fixNames(f.name)
	else:
		fn = sys.argv[1]
		fixNames(fn)
	
main()
