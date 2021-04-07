import sys
import os
from pathlib import Path
from checkinsertion import getHits, consolidate

VIRUS_NAMES = ['cfav', 'liao_ning', 'xincheng', 'anphevirus', 'horseradish',
               'aus_anopheles']
VIRUS_DBS = ['cfavdb', 'liao_ning_db', 'xinchengdb', 'anphevirusdb',
             'horseradish_db', 'aus_anophelesdb']

AEGYPTI = ['ALL', 'Amacuzac_Mexico', 'Bangkok_Thailand', 'Cairns_Australia',
           'Cebu_City_Philippines', 'Cuanda_Angola', 'El_Dorado_Argentina',
		   'HoChiMin_Vietnam', 'La_Lope-Gabon', 'Maricopa_County_Arizona-USA',
		   'Paqueta-Brazil', 'Skukusa_South_Africa', 'Tahiti_FrenchPolynesia']

def doIt(aegyptiDir, virusName, virusDB):
	""" Assemble combined reads for one region and one virus."""

	print('Processing ' + aegyptiDir + ' for virus ' + virusName + '...')
	combinedDir = aegyptiDir + '_combined'
	if not Path(combinedDir).exists():
		os.system('mkdir ' + combinedDir)
	os.chdir(combinedDir)  # combined directory name

	if aegyptiDir == 'ALL':
		prefix = ''
	else:
		prefix = aegyptiDir[:6]

	if not Path('spades_' + virusName).exists():
		os.system('mkdir spades_' + virusName)
	else:
		print(combinedDir + '/spades_' + virusName + ' already exists.  Skipping...')
		os.chdir('..')
		return
	print('Copying files for assembly...')
	os.system('cat ../' + prefix + '*/spades_' + virusName + '/viral1.fastq > spades_' + virusName + '/viral1.fastq')
	os.system('cat ../' + prefix + '*/spades_' + virusName + '/viral2.fastq > spades_' + virusName + '/viral2.fastq')
	os.system('cat ../' + prefix + '*/spades_' + virusName + '/viral_single.fastq > spades_' + virusName + '/viral_single.fastq')

	print('6/9: Assembling reads with SPAdes...')
	if os.path.getsize('spades_' + virusName + '/viral_single.fastq') > 0:
		result = os.system('/home/havill/bin/SPAdes-3.13.1-Linux/bin/spades.py --pe1-1 spades_' + virusName + '/viral1.fastq --pe1-2 spades_' + virusName + '/viral2.fastq --s1 spades_' + virusName + '/viral_single.fastq -o spades_' + virusName)
	else:
		os.remove('spades_' + virusName + '/viral_single.fastq')
		result = os.system('/home/havill/bin/SPAdes-3.13.1-Linux/bin/spades.py --pe1-1 spades_' + virusName + '/viral1.fastq --pe1-2 spades_' + virusName + '/viral2.fastq -o spades_' + virusName)

	if Path('spades_' + virusName + '/scaffolds.fasta').exists():  # assembly succeeded
		# BLAST scaffolds against viral and aegypti genomes.
		print('7/9: BLASTing scaffolds against viral genome...')
		os.system('blastn -query spades_' + virusName + '/scaffolds.fasta -db ' + virusDB + ' -num_threads 8 -task blastn -evalue 0.0001 -outfmt "10 qseqid qstart qend qseq sstart send sseq evalue bitscore sseqid" -out spades_' + virusName + '/blast_scaffolds_' + virusName + '.csv')

		print('8/9: BLASTing scaffolds against Aedes aegypti genome...')
		os.system('blastn -query spades_' + virusName + '/scaffolds.fasta -db aegyptidb -num_threads 8 -outfmt "10 qseqid qstart qend qseq sseqid sstart send sseq evalue bitscore" -out spades_' + virusName + '/blast_scaffolds_aa.csv')

		# Combine BLAST results to find putative insertions.
		print('9/9: Combining BLAST results to locate putative insertions...')
		getHits(Path('.').resolve(), virusName)
	else:
		print('Assembly failed.  No results found.')

	print('Copying hits to ../HITS_' + virusName + '...')
	if not Path('../HITS_' + virusName).exists():
		os.mkdir('../HITS_' + virusName)
	os.system('cp -v spades_' + virusName + '/*viral_hits* ' + '../HITS_' + virusName)

	if not Path('../HITS_' + virusName + '_combined').exists():
		os.mkdir('../HITS_' + virusName + '_combined')
	os.system('cp -v spades_' + virusName + '/*viral_hits* ' + '../HITS_' + virusName + '_combined')

	print('Done.')
	os.chdir('..')

def main():
	# assumes cwd is data/egypti/analyzed
	#doIt(AEGYPTI[1], 'cfav', 'cfavdb')
	#consolidate('cfav_combined')
	#return

	for index in range(len(VIRUS_NAMES)):
		virusName = VIRUS_NAMES[index]
		virusDB = VIRUS_DBS[index]
		print('Beginning to process virus ' + virusName + '...')
		for aegypti in AEGYPTI:
			doIt(aegypti, virusName, virusDB)

		print('Consolidating results for ' + virusName + '...')
		consolidate(virusName + '_combined')

main()
