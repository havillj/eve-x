import sys
import os.path
from pathlib import Path
import pysam

# def pairs(filename):
# 	samfile = pysam.AlignmentFile(filename, 'rb')
# 	samfile2 = pysam.AlignmentFile(filename, 'rb')
# 	
# 	print('Total reads =', samfile.mapped + samfile.unmapped)
# 	
# 	count = 0
# 	count_unmapped = 0
# 	count_total = 0
# 	for read in samfile:
# 		count_total += 1
# 		if not read.is_supplementary and read.is_paired and not read.is_proper_pair:
# 			count += 1
# 			print(count)
# #			print('\n' + read.query_name)
# #			print(read.seq)
# 			
# # 			print(read)
# # 			print()
# 			try:
# 				samfile2.mate(read)
# # 				print(samfile2.mate(read))
# 			except ValueError:
# 				count_unmapped += 1
# # 				print('Mate unmapped.')
# #			print('\n')
# 			
# #			return
# 			
# 	samfile.close()
# 	samfile2.close()
# 	
# 	print('Total reads =', count_total)
# 	print('Paired but not proper =', count)
# 	print('Paired, not proper, mate not mapped =', count_unmapped)


# stdin: samtools view 
# stdout: fastq file
# def pairs():
# 	count = 0
# 	count_unmapped = 0
# 	count_total = 0
# 	reads = {}
# 	for line in sys.stdin:
# 		cols = line.rstrip().split('\t')
# 		flag = int(cols[1])
# 		count_total += 1
# 		if flag & 3 == 1:  # paired and not in proper pair
# 			if flag & 12 == 4:   # unmapped and pair is mapped
# 				seq = cols[9]
# 				if seq.count('N') / len(seq) <= 0.5:
# 					count += 1
# 					qname = cols[0]
# #					rname = cols[3]
# #					cigar = cols[5]
# 					qual = cols[10]
# #					print('>' + qname + '_' + str(flag) + '_' + rname + '\n' + seq)
# 					print('@' + qname + '\n' + seq + '\n+\n' + qual)
# 	
# 	print('Total reads =', count_total)
# 	print('Paired but unmapped =', count)
# #	print('Paired, not proper, mate not mapped =', count_unmapped)

def pairs(filenameBAM):
	bamfile = pysam.AlignmentFile(filenameBAM, 'rb', threads = 8)
	fastq = open(filenameBAM + '_unmapped_mates.fastq', 'w')
	
	count = 0
	count_unmapped = 0
	count_total = 0
	reads = {}
	for read in bamfile:
		count_total += 1
		if read.is_paired and not read.is_proper_pair:  # paired and not in proper pair
			if read.is_unmapped and not read.mate_is_unmapped:   # unmapped and mate is mapped
				seq = read.query_sequence
				if seq.count('N') / len(seq) <= 0.5:
					count += 1
					quality_string = ''
					for s in read.query_qualities:
						quality_string += chr(s + 0x21)
					fastq.write('@' + read.query_name + '\n' + seq + '\n+\n' + quality_string + '\n')
	bamfile.close()
	fastq.close()
	
	print('Total reads =', count_total)
	print('Paired but unmapped =', count)


# def clean():
# 	header = input()
# 	while header:
# 		seq = input()
# 		plus = input()
# 		qual = input()
# #		blank = input()
# 		if seq.count('N') <= 5:
# 			print(header + '\n' + seq + '\n+\n' + qual)
# 		try:
# 			header = input()
# 		except EOFError:
# 			break
	
def getMappedMates(viralMatchesDir, aedesFileName):
	if viralMatchesDir[-1] != '/':
		viralMatchesDir += '/'
		
	# get all reads that mapped to viral genomes
	
	viralQueries = {}
	path = Path(viralMatchesDir)
	for accessionDir in path.iterdir():
		if accessionDir.is_dir():
			print(accessionDir.name)
	
			viralSAM = pysam.AlignmentFile(str(accessionDir) + '/consensus_alignment_sorted.bam', 'rb')
			for read in viralSAM:
				if read.query_name in viralQueries:
					print(read.query_name)
				viralQueries[read.query_name] = (read.query_sequence, accessionDir.name)
			viralSAM.close()
	
#	print(viralQueries)
	
	# get mates for the viral reads from the aedes aegypti BAM file
	
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
	
# 	for qn in viralQueries:
# 		print(qn + ',' + viralQueries[qn] + ',' + viralAccession + ',,' + mates[qn].query_sequence + ',' + mates[qn].reference_name + ',' + str(mates[qn].reference_start + 1))
	
	pysam.sort('-o', viralMatchesDir + '/mapped_mates_sorted.bam', viralMatchesDir + '/mapped_mates.bam')
	pysam.index(viralMatchesDir + '/mapped_mates_sorted.bam')
	
	csv = open(viralMatchesDir + '/results.csv', 'w')
	sortedBAM = pysam.AlignmentFile(viralMatchesDir + '/mapped_mates_sorted.bam', 'rb')
	csv.write('Read ID,Read Seq,Virus Match,Virus Match Name,Mate Seq,Aegypti Reference,Position\n')
	for read in sortedBAM:
		qn = read.query_name
		csv.write(qn + ',' + viralQueries[qn][0] + ',' + viralQueries[qn][1] + ',,' + mates[qn].query_sequence + ',' + mates[qn].reference_name + ',' + str(mates[qn].reference_start + 1) + '\n')
	sortedBAM.close()
	csv.close()
	

def main():
	if len(sys.argv) == 2:
		pairs(sys.argv[1])
	else:
		getMappedMates(sys.argv[1], sys.argv[2])
	
main()

# Mexico
# 317235754 total
# 101683138 are paired but not proper

# bin/samtools flagstat Amacuzac_Mexico-_12.2-116626.sorted.deduped.merged.bam
# 317235754 + 0 in total (QC-passed reads + QC-failed reads)
# 0 + 0 secondary
# 18385124 + 0 supplementary
# 34339111 + 0 duplicates
# 298206884 + 0 mapped (94.00% : N/A)
# 298850630 + 0 paired in sequencing
# 149425315 + 0 read1
# 149425315 + 0 read2
# 212185516 + 0 properly paired (71.00% : N/A)
# 274580490 + 0 with itself and mate mapped
# 5241270 + 0 singletons (1.75% : N/A)
# 39657766 + 0 with mate mapped to a different chr
# 15303750 + 0 with mate mapped to a different chr (mapQ>=5)