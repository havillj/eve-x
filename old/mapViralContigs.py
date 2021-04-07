import sys
import os
import pysam
from pathlib import Path
from multiprocessing import Process
from getreadsfromblast import getReads
from checkinsertion import getHits, consolidate

# VIRUS_NAME = 'cfav'
# VIRUS_DB = 'cfavdb'

# VIRUS_NAME = 'liao_ning'
# VIRUS_DB = 'liao_ning_db'

# VIRUS_NAME = 'xincheng'
# VIRUS_DB = 'xinchengdb'

# VIRUS_NAME = 'anphevirus'
# VIRUS_DB = 'anphevirusdb'

# VIRUS_NAME = 'horseradish'
# VIRUS_DB = 'horseradish_db'

# VIRUS_NAME = 'aus_anopheles'
# VIRUS_DB = 'aus_anophelesdb'

VIRUS_NAMES = ['cfav', 'liao_ning', 'xincheng', 'anphevirus', 'horseradish', 'aus_anopheles']
VIRUS_DBS = ['cfavdb', 'liao_ning_db', 'xinchengdb', 'anphevirusdb', 'horseradish_db', 'aus_anophelesdb']
VIRUS_NAME = ''
VIRUS_DB = ''

def doit(dir, filenameBAM):
	os.chdir(dir)
	try:
		index = filenameBAM.index('.sorted')
	except:
		index = filenameBAM.index('.bam')
	newBAM_fn = filenameBAM[:index] + '_unmapped.bam'

	print('\nAnalyzing ' + filenameBAM)

	# Get all unmapped reads.

	print('1/9: Getting unmapped reads...')
#	os.system('python3 /home/havill/bin/aedes.py -1 ' + filenameBAM)
#	os.system('samtools view -b -f 4 ' + filenameBAM + ' > ' + newBAM_fn)  # or python3 ~/bin/aedes.py -1 fn # => newBAM_fn

	if not Path(newBAM_fn).exists():
		bamfile = pysam.AlignmentFile(filenameBAM, 'rb', threads = 8)
		newBAM = pysam.AlignmentFile(newBAM_fn, 'wb', template = bamfile, threads = 8)

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
	else:
		print('Skipping - ' + newBAM_fn + ' exists.')

	# Convert unmapped reads to FASTA.

	print('2/9: Writing unmapped reads to a FASTA file...')
	newFASTA_fn = newBAM_fn[:-4] + '.fasta'
	if not Path(newFASTA_fn).exists():
		os.system('samtools fasta ' + newBAM_fn + ' > ' + newFASTA_fn)
	else:
		print('Skipping - ' + newFASTA_fn + ' exists.')

	# BLAST unmapped reads against viral genome.

	print('3/9: Searching for viral reads with BLAST...')
	os.system('blastn -query ' + newFASTA_fn + ' -db ' + VIRUS_DB + ' -num_threads 8 -task blastn -evalue 0.0001 -outfmt "10 qaccver saccver stitle pident length mismatch gapopen qstart qend sstart send evalue bitscore" -out blast_unmapped_' + VIRUS_NAME + '.csv')

	# Get matching reads in a BAM file.

	print('4/9: Getting viral sequences (plus paired-end mates)...')
	getReads('blast_unmapped_' + VIRUS_NAME + '.csv', filenameBAM, 0, 12000, True, VIRUS_NAME)
	viralBAM_fn = filenameBAM[:index] + '_' + VIRUS_NAME + '.bam'

	# Convert matching reads to FASTQ.

	print('5/9: Writing viral paired-end reads to FASTQ files...')
	os.system('mkdir spades_' + VIRUS_NAME)
	os.system('samtools fastq -1 spades_' + VIRUS_NAME + '/viral1.fastq -2 spades_' + VIRUS_NAME + '/viral2.fastq -s spades_' + VIRUS_NAME + '/viral_single.fastq ' + viralBAM_fn)

	# Assemble matching reads.

	print('6/9: Assembling reads with SPAdes...')
	if os.path.getsize('spades_' + VIRUS_NAME + '/viral_single.fastq') > 0:
		result = os.system('/home/havill/bin/SPAdes-3.13.1-Linux/bin/spades.py --pe1-1 spades_' + VIRUS_NAME + '/viral1.fastq --pe1-2 spades_' + VIRUS_NAME + '/viral2.fastq --s1 spades_' + VIRUS_NAME + '/viral_single.fastq -o spades_' + VIRUS_NAME)
	else:
		os.remove('spades_' + VIRUS_NAME + '/viral_single.fastq')
		result = os.system('/home/havill/bin/SPAdes-3.13.1-Linux/bin/spades.py --pe1-1 spades_' + VIRUS_NAME + '/viral1.fastq --pe1-2 spades_' + VIRUS_NAME + '/viral2.fastq -o spades_' + VIRUS_NAME)

	if Path('spades_' + VIRUS_NAME + '/scaffolds.fasta').exists():  # assembly succeeded
		# BLAST scaffolds against viral and aegypti genomes.
		print('7/9: BLASTing scaffolds against viral genome...')
		os.system('blastn -query spades_' + VIRUS_NAME + '/scaffolds.fasta -db ' + VIRUS_DB + ' -num_threads 8 -task blastn -evalue 0.0001 -outfmt "10 qseqid qstart qend qseq sstart send sseq evalue bitscore sseqid" -out spades_' + VIRUS_NAME + '/blast_scaffolds_' + VIRUS_NAME + '.csv')

		print('8/9: BLASTing scaffolds against Aedes aegypti genome...')
		os.system('blastn -query spades_' + VIRUS_NAME + '/scaffolds.fasta -db aegyptidb -num_threads 8 -outfmt "10 qseqid qstart qend qseq sseqid sstart send sseq evalue bitscore" -out spades_' + VIRUS_NAME + '/blast_scaffolds_aa.csv')

		# Combine BLAST results to find putative insertions.
		print('9/9: Combining BLAST results to locate putative insertions...')
		getHits(Path('.').resolve(), VIRUS_NAME)
	else:
		print('Assembly failed.  No results found.')

	print('Done.')

def process_cmdline():
	if (sys.argv[1] == '-a'):
		if len(sys.argv) >= 3:
			pattern = sys.argv[2]
		else:
			pattern = ''
		dir = Path('.').resolve()   # ~/data/aegypti/analyzed/
		# workers = []
		# for subdir in dir.iterdir():
		# 	if subdir.is_dir() and subdir.name[:len(pattern)] == pattern and \
		# 	  not Path(subdir / ('spades_' + VIRUS_NAME)).exists():
		# 		for file in subdir.iterdir():
		# 			if not file.is_dir() and (file.name[-26:] == '.sorted.deduped.merged.bam'):
		# 				p = Process(target = doit, args = (subdir.resolve(), file.name))
		# 				workers.append(p)
		# 				p.start()
		# 				break
		#
		# for p in workers:
		# 	p.join()

		print('Consolidating results...')
		os.chdir(dir)
		if not Path('HITS_' + VIRUS_NAME).exists():
			os.mkdir('HITS_' + VIRUS_NAME)
		#os.system('cp -v */spades_' + VIRUS_NAME + '/*viral_hits* ' + 'HITS_' + VIRUS_NAME)
		consolidate(VIRUS_NAME)
	elif len(sys.argv) == 2:
		filenameBAM = sys.argv[1]
		doit(Path('.').resolve(), filenameBAM)

def main():
	global VIRUS_NAME
	global VIRUS_DB

	for index in range(len(VIRUS_NAMES)):
		VIRUS_NAME = VIRUS_NAMES[index]
		VIRUS_DB = VIRUS_DBS[index]
		process_cmdline()

if __name__ == '__main__':
	main()
