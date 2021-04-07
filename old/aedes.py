#!/usr/bin/env python3

import sys
import os.path
from pathlib import Path
import pysam
from Bio import Entrez
              
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
		return ''

# def getUnmappedReads(filenameBAM):
# 	"""
# 	Get all reads that are unmapped but have mapped mates from a genomic BAM file.
# 	
# 	filenameBAM: file name of genomic BAM file
# 	
# 	"""
# 	
# 	bamfile = pysam.AlignmentFile(filenameBAM, 'rb', threads = 8)
# 	
# 	index = filenameBAM.index('.sorted')
# 	fastq = open(filenameBAM[:index] + '_unmapped_mates.fastq', 'w')
# 	
# 	count = 0
# 	count_unmapped = 0
# 	count_total = 0
# 	reads = {}
# 	for read in bamfile:
# 		count_total += 1
# 		if read.is_paired and not read.is_proper_pair:  # paired and not in proper pair
# 			if read.is_unmapped and not read.mate_is_unmapped:   # unmapped and mate is mapped
# 				seq = read.query_sequence
# 				if seq.count('N') / len(seq) <= 0.5:
# 					count += 1
# 					quality_string = ''
# 					for s in read.query_qualities:
# 						quality_string += chr(s + 0x21)
# 					fastq.write('@' + read.query_name + '\n' + seq + '\n+\n' + quality_string + '\n')
# 	bamfile.close()
# 	fastq.close()
# 	
# 	print('Total reads =', count_total)
# 	print('Paired but unmapped =', count, '({:.2f}%)'.format(100*count/count_total))
	
# def getAllUnmappedReads(filenameBAM):
# 	"""
# 	Get all unmapped reads from a genomic BAM file.
# 	
# 	filenameBAM: file name of genomic BAM file
# 	
# 	"""
# 	
# 	bamfile = pysam.AlignmentFile(filenameBAM, 'rb', threads = 8)
# 	
# 	try:
# 		index = filenameBAM.index('.sorted')
# 	except:
# 		index = filenameBAM.index('.bam')
# 		
# 	fastq = open(filenameBAM[:index] + '_unmapped.fastq', 'w')
# 	
# 	count = 0
# 	count_total = 0
# 	reads = {}
# 	for read in bamfile:
# 		count_total += 1
# 		if read.is_unmapped:
# 			seq = read.query_sequence
# 			if seq.count('N') / len(seq) <= 0.5:
# 				count += 1
# 				quality_string = ''
# 				for s in read.query_qualities:
# 					quality_string += chr(s + 0x21)
# 				fastq.write('@' + read.query_name + '\n' + seq + '\n+\n' + quality_string + '\n')
# 	bamfile.close()
# 	fastq.close()
# 	
# 	print('Total reads =', count_total)
# 	print('Unmapped reads =', count, '({:.2f}%)'.format(100*count/count_total))

def getAllUnmappedReads(filenameBAM):
	"""
	Get all unmapped reads from a genomic BAM file.
	
	filenameBAM: file name of genomic BAM file
	
	"""
	
	bamfile = pysam.AlignmentFile(filenameBAM, 'rb', threads = 8)
	
	try:
		index = filenameBAM.index('.sorted')
	except:
		index = filenameBAM.index('.bam')
		
	newBAM_fn = filenameBAM[:index] + '_unmapped.bam'
	newBAM = pysam.AlignmentFile(newBAM_fn, 'wb', template = bamfile)
	
	count = 0
	count_total = 0
	for read in bamfile:
		count_total += 1
		if read.is_unmapped:
			seq = read.query_sequence
			if seq.count('N') / len(seq) <= 0.5:
				count += 1
				newBAM.write(read)
	bamfile.close()
	newBAM.close()
	
	print('Total reads =', count_total)
	print('Unmapped reads =', count, '({:.2f}%)'.format(100*count/count_total))
	
#	fastq_text = pysam.fastq(newBAM_fn)
	fastq = open(filenameBAM[:index] + '_unmapped.fastq', 'w')
	fastq.write(pysam.fastq(newBAM_fn))
	fastq.close()
	
# def getAllUnmappedPairs2(filenameBAM):
# 	"""
# 	Get all unmapped pairs and write to 2 (forward and reverse) fastq files.
# 	
# 	filenameBAM: file name of genomic BAM file
# 	
# 	"""
# 	
# 	bamfile = pysam.AlignmentFile(filenameBAM, 'rb', threads = 8)
# 	
# 	index = filenameBAM.index('.sorted')
# 	fastq_f = open(filenameBAM[:index] + '_unmapped_pairs1.fastq', 'w')
# 	fastq_r = open(filenameBAM[:index] + '_unmapped_pairs2.fastq', 'w')
# 	fwd = {}
# 	rev = {}
# 	
# 	count = 0
# 	count_unmapped = 0
# 	count_total = 0
# 	reads = {}
# 	for read in bamfile:
# 		count_total += 1
# 		if read.is_unmapped:
# 			seq = read.query_sequence
# 			if seq.count('N') / len(seq) <= 0.5:
# 				count += 1
# 				quality_string = ''
# 				for s in read.query_qualities:
# 					quality_string += chr(s + 0x21)
# 				if read.is_read1:
# 					fwd[read.query_name] = (seq, quality_string)
# 				else:
# 					rev[read.query_name] = (seq, quality_string)
# 					
# 	for qn in fwd:
# 		if qn in rev:
# 			fastq_f.write('@' + qn + '\n' + fwd[qn][0] + '\n+\n' + fwd[qn][1] + '\n')
# 			fastq_r.write('@' + qn + '\n' + rev[qn][0] + '\n+\n' + rev[qn][1] + '\n')
# 					
# 	bamfile.close()
# 	fastq_f.close()
# 	fastq_r.close()
# 	
# 	print('Total reads =', count_total)
# 	print('Unmapped =', count, '({:.2f}%)'.format(100*count/count_total))
	
def getAllUnmappedPairs(filenameBAM):
	"""
	Get all pairs in which at least one read is unmapped, and write to 2 (forward and reverse) fastq files.
	
	filenameBAM: file name of genomic BAM file
	
	"""
	
	bamfile = pysam.AlignmentFile(filenameBAM, 'rb', threads = 8)
	
	index = filenameBAM.index('.sorted')
	fastq_f = open(filenameBAM[:index] + '_unmapped_pairs1.fastq', 'w')
	fastq_r = open(filenameBAM[:index] + '_unmapped_pairs2.fastq', 'w')
	fwd = {}
	rev = {}
	
	count = 0
	count_total = 0
	count_singles = 0
	reads = {}
	
	for read in bamfile:
		count_total += 1
		if read.is_paired and (read.is_unmapped or read.mate_is_unmapped):  ####
			seq = read.query_sequence
			if seq.count('N') / len(seq) <= 0.5:  # may omit a pair
				if read.is_unmapped:
					count += 1
				quality_string = ''
				for s in read.query_qualities:
					quality_string += chr(s + 0x21)
				if read.is_read1:
					fwd[read.query_name] = (seq, quality_string)
				else:
					rev[read.query_name] = (seq, quality_string)
		elif not read.is_paired and read.is_unmapped:
			count_singles += 1
					
	for qn in fwd:
		if qn in rev:
			fastq_f.write('@' + qn + '\n' + fwd[qn][0] + '\n+\n' + fwd[qn][1] + '\n')
			fastq_r.write('@' + qn + '\n' + rev[qn][0] + '\n+\n' + rev[qn][1] + '\n')
					
	bamfile.close()
	fastq_f.close()
	fastq_r.close()
	
	print('Total reads =', count_total)
	print('Unmapped paired-end reads =', count, '({:.2f}%)'.format(100*count/count_total))
	print('Unmapped single-end reads (omitted) =', count_singles)

def getMappedMates(viralMatchesDir, aedesFileName):
	"""
	Get mates of unmapped reads identified as viral by GD.
	
	viralMatchesDir: name of directory containing a subdirectory for each virus matched
	                 each subdirectory contains a file named consensus_alignment_sorted.bam
	aedesFileName:   name of original genomic BAM file
	
	"""
	
	if viralMatchesDir[-1] != '/':
		viralMatchesDir += '/'
		
	# get all reads that mapped to viral genomes
	# store in a dictionary {query_name: (query_sequence, accession number), ...}
	
	viralQueries = {}
	path = Path(viralMatchesDir)
	for accessionDir in path.iterdir():
		if accessionDir.is_dir():
			print(accessionDir.name)
	
			viralSAM = pysam.AlignmentFile(str(accessionDir) + '/consensus_alignment_sorted.bam', 'rb')
			for read in viralSAM:
				query = read.query_name
				if query[-2] == '/':
					query = query[:-2]   # remove /1 or /2 to match BAM query_name
				if query in viralQueries:
					print(query)
				viralQueries[query] = (read.query_sequence, accessionDir.name, read.reference_start + 1, read.reference_end + 1)
			viralSAM.close()
	
	# mates for the viral reads from the aedes aegypti BAM file => new file "mapped_mates.bam"
	
	aedesSAM = pysam.AlignmentFile(aedesFileName, 'rb', threads = 8)
	newSAM = pysam.AlignmentFile(viralMatchesDir + '/mapped_mates.bam', 'wb', template = aedesSAM)
	mates = {}
	for read in aedesSAM:
		if read.query_name in viralQueries:
			if not read.is_unmapped and not read.is_supplementary and not read.is_secondary:
				newSAM.write(read)
				mates[read.query_name] = read
	newSAM.close()
	aedesSAM.close()
		
	# Sort and index mapped_mates.bam => mapped_mates_sorted.bam
	
	pysam.sort('-o', viralMatchesDir + '/mapped_mates_sorted.bam', viralMatchesDir + '/mapped_mates.bam')
	pysam.index(viralMatchesDir + '/mapped_mates_sorted.bam')
	
	# Write data to CSV file named results.csv
	
	names = {}
	
	try:
		index = aedesFileName.index('.sorted')
	except:
		index = aedesFileName.index('.bam')
		
	csv = open(aedesFileName[:index] + '_RESULTS.csv', 'w')
	sortedBAM = pysam.AlignmentFile(viralMatchesDir + '/mapped_mates_sorted.bam', 'rb')
#	csv.write('Read ID,Read Seq,Virus Match,Virus Match Name,Mate Seq,Aegypti Reference,Position\n')
	csv.write('Read ID,Read Seq,Virus Match,Virus Match Name,Ref Start,Ref End,Mate Seq,Aegypti Reference,Ref Start,Ref End\n')
	for read in sortedBAM:
		qn = read.query_name
		acc = viralQueries[qn][1]
		
		if '.' in acc:
			acc = acc.split('.')[0]
		if acc in names:
			virusName = names[acc]
		else:
			virusName = getName(acc)
			if virusName != '':
				names[acc] = virusName
					
# 		if '.' in acc:
# 			virusName = getName(acc.split('.')[0])
# 		else:
# 			virusName = getName(acc)
#		csv.write(qn + ',' + viralQueries[qn][0] + ',' + acc + ','+ virusName + ',' + mates[qn].query_sequence + ',' + mates[qn].reference_name + ',' + str(mates[qn].reference_start + 1) + '\n')
		csv.write(qn + ',' + viralQueries[qn][0] + ',' + acc + ','+ virusName + ',' + str(viralQueries[qn][2]) + ',' + str(viralQueries[qn][3]) + ',' + mates[qn].query_sequence + ',' + mates[qn].reference_name + ',' + str(mates[qn].reference_start + 1) + ',' + str(mates[qn].reference_end + 1) + '\n')
	sortedBAM.close()
	csv.close()
	
def printUsage(command):
	print('Usage:')
	print('   python3 ' + command + ' -1 genomic_file.bam')
	print('or')
	print('   python3 ' + command + ' -1')
	print('or')
	print('   python3 ' + command + ' -1 -r')
	print('or')
	print('   python3 ' + command + ' -2 viral_matches_dir genomic_file.bam')

def processDir(dir):
	for file in dir.iterdir():
		if file.is_dir():
			continue
		fn = str(file)
		try:
			index = fn.index('.sorted.deduped.merged.bam')
			if fn[-4:] != '.bam':
				continue
		except:
			try:
				index = fn.index('.bam')
			except:
				continue
	
		if not Path(fn[:index] + '_unmapped.fastq').exists():
			print('Processing ' + fn + '...')
			getAllUnmappedReads(fn)
# 		else:
# 			if not Path(fn[:index] + '_unmapped_mates.fastq').exists():
# 				print('Processing ' + fn + '...')
# 				getUnmappedReads(fn)
				
def main():
	if len(sys.argv) < 2 or not(('-1' in sys.argv) or ('-2' in sys.argv)):
		printUsage(sys.argv[0])
		return 1
		
	if '-1' in sys.argv and len(sys.argv) >= 2:
		if len(sys.argv) == 3 and '-r' not in sys.argv:  # one named file
			getAllUnmappedReads(sys.argv[2])
			return 0
		else:
			cwd = Path('.')
			if len(sys.argv) == 2 and '-r' not in sys.argv:  # bam file(s) in cwd
				processDir(cwd)
			else:                                            # bam files in all subdirs of cwd
				dirs = [d for d in cwd.iterdir() if d.is_dir()]  
				for dir in dirs:
					processDir(dir.resolve())
			return 0
	elif '-2' in sys.argv and len(sys.argv) == 4:
		getMappedMates(sys.argv[2], sys.argv[3])
		return 0
	else:
		printUsage(sys.argv[0])
		return 1
	
main()
