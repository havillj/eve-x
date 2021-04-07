import sys
import os.path
from pathlib import Path
import pysam
from Bio import Entrez

def getMappedMates2(csvFileName, aedesFileName):
	"""
	Get mates of unmapped reads identified as viral.
	
	csvFileName:     name of csv file returned by BLAST
	aedesFileName:   name of original genomic BAM file
	
	"""
		
	# get all reads that mapped to viral genomes
	# store in a dictionary {query_name: (query_sequence, accession number), ...}
	
	viralQueries = {}
	
	viralCSV = open(csvFileName, 'r')
	for line in viralCSV:
		line = line.rstrip()
		cols = line.split(',')
		query = cols[0]
		readNum = 1
		if query[-2] == '/':
			readNum = int(query[-1])
			query = query[:-2]   # remove /1 or /2 to match BAM query_name
		if query not in viralQueries:
			viralQueries[query] = [None, None]
		viralQueries[query][readNum - 1] = line
	viralCSV.close()
	
	# mates for the viral reads from the aedes aegypti BAM file => new file "mapped_mates.bam"
	
	aedesSAM = pysam.AlignmentFile(aedesFileName, 'rb', threads = 8)
	newSAM = pysam.AlignmentFile('mapped_mates2.bam', 'wb', template = aedesSAM)
	mates = {}
	for read in aedesSAM:
		if read.query_name in viralQueries:
			if not read.is_unmapped and not read.is_supplementary and not read.is_secondary:
				newSAM.write(read)
				mates[read.query_name] = read
	newSAM.close()
	aedesSAM.close()
		
	# Sort and index mapped_mates.bam => mapped_mates_sorted.bam
	
	pysam.sort('-o', 'mapped_mates_sorted2.bam', 'mapped_mates2.bam')
	pysam.index('mapped_mates_sorted2.bam')
	
	# Write data to CSV file named results.csv
	
	written = []
	
	index = aedesFileName.index('.sorted')
	csv = open(aedesFileName[:index] + '_MATES_FROM_BLAST.csv', 'w')
	sortedBAM = pysam.AlignmentFile('mapped_mates_sorted2.bam', 'rb')
#	csv.write('Read ID,Aegypti Reference,Position\n')
	for read in sortedBAM:
		qn = read.query_name
		if read.is_read1:
			readNum = 1   # orig read is read 2
		else:
			readNum = 0   # orig read is read 1
		csv.write(viralQueries[qn][readNum] + ',' + mates[qn].reference_name + ',' + str(mates[qn].reference_start + 1) + ',' + str(mates[qn].reference_end + 1) + '\n')
		written.append(qn)
	sortedBAM.close()
	
	for qn in viralQueries:
		if qn not in written:
			if viralQueries[qn][0] is not None:
				csv.write(viralQueries[qn][0] + '\n')
			if viralQueries[qn][1] is not None:
				csv.write(viralQueries[qn][1] + '\n')
	csv.close()
	
def main():
	getMappedMates2(sys.argv[1], sys.argv[2])
	
main()
