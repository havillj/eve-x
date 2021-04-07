#!/usr/bin/env python3

import sys
import os.path
from pathlib import Path
import pysam

def getReads(filenameCSV, filenameBAM, left, right, getMates, virusName):

	csv = open(filenameCSV, 'r')
	queries = {}
	for line in csv:
		cols = line.split(',')
		query = cols[0]
# 		if query[-2] == '/':
# 			query = query[:-2]   # remove /1 or /2 to match BAM query_name
		start = int(cols[9])
		end = int(cols[10])
		if (start < end and start >= left and end <= right) or \
		   (end < start and end >= left and start <= right):
			queries[query] = ''
	csv.close()

	bamfile = pysam.AlignmentFile(filenameBAM, 'rb', threads = 8)
	index = filenameBAM.index('.sorted')
	bamfile_out = pysam.AlignmentFile(filenameBAM[:index] + '_' + virusName + '.bam', 'wb', template = bamfile, threads = 8)

	written = {}
	for read in bamfile:
		if read.is_read1:
			qn = read.query_name + '/1'
			qn_mate = read.query_name + '/2'
		elif read.is_read2:
			qn = read.query_name + '/2'
			qn_mate = read.query_name + '/1'
		else:
			qn = read.query_name
		read.query_name += '_VIRAL'
		if getMates:
			if (qn in queries or qn_mate in queries) and qn not in written:
				written[qn] = True
				bamfile_out.write(read)
		else:
			if qn in queries and qn not in written:
				bamfile_out.write(read)
				written[qn] = True

	bamfile.close()
	bamfile_out.close()

def main():
	if len(sys.argv) == 5:
		getReads(sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4]), True, 'cfav')

if __name__ == '__main__':
	main()
