import sys
import os.path
from pathlib import Path
import pysam
from Bio import Entrez

def checkEnds(filenameBAM, seq):
	bamfile = pysam.AlignmentFile(filenameBAM, 'rb', threads = 8)
	
	length = len(seq)
	for read in bamfile:
		if (read.query_sequence[:length] == seq) or (read.query_sequence[-length:] == seq):
			print(read.query_name, read.query_sequence)
			
	bamfile.close()
	
def main():
	checkEnds(sys.argv[1], 'CGGGTAACGTCGGTAG')  # CTCGTCGCCACTTTTA')
	
main()
