#!/usr/bin/env python3

import sys
import os.path
from pathlib import Path
import pysam

def getReads(filenameBAM, left, right):
	
	bamfile = pysam.AlignmentFile(filenameBAM, 'rb', threads = 8)
	index = filenameBAM.index('.bam')
	bamfile_out = pysam.AlignmentFile(filenameBAM[:index] + '_nonspanning.bam', 'wb', template = bamfile)
	for read in bamfile:
		print(read.reference_start, read.reference_end, end='')
		if read.reference_end is not None and (read.reference_start > left or read.reference_end < right):
			bamfile_out.write(read)
			print(' good')
		else:
			print()
	bamfile.close()
	bamfile_out.close()
	
def main():
#137741225 - 242
	if len(sys.argv) == 4:
		getReads(sys.argv[1], int(sys.argv[2]), int(sys.argv[3]))
		
main()
