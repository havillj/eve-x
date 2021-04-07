import sys
import os
from pathlib import Path

def readCSV(csvFilename):
	csv = open(csvFilename, 'r')
	hits = {}
	for line in csv:
		line = line.rstrip()
		row = line.split(',')
		contig = row[0]
		if contig in hits:
			hits[contig].append(tuple(row[1:]))
		else:
			hits[contig] = [tuple(row[1:])]
	csv.close()

	return hits

def getHits(dir, virusName):
	csvVirusFilename = str(dir / ('spades_' + virusName + '/blast_scaffolds_' + virusName + '.csv'))
	csvAedesFilename = str(dir / ('spades_' + virusName + '/blast_scaffolds_aa.csv'))

	outFilenameBothEnds = str(dir / ('spades_' + virusName) / (dir.name + '_viral_hits2.txt'))
	outFilenameOneEnd = str(dir / ('spades_' + virusName) / (dir.name + '_viral_hits1.txt'))

	# os.system('blastn -query ' + str(dir) + \
    #   '/spadesCFAV/scaffolds.fasta -db aegyptidb -num_threads 8 -outfmt "10 qseqid qstart qend qseq sseqid sstart send sseq evalue bitscore" -out ' + \
    #   csvAedesFilename)
	#
	# os.system('blastn -query ' + str(dir) + \
	#   '/spadesCFAV/scaffolds.fasta -db cfavdb -num_threads 8 -task blastn -evalue 0.0001 -outfmt "10 qseqid qstart qend qseq sstart send sseq evalue bitscore" -out ' + \
	#   csvVirusFilename)

	virusHits = readCSV(csvVirusFilename)
	aedesHits = readCSV(csvAedesFilename)

	# Combine overlapping viral hits in each contig into long hits
	newVirusHits = {}
	for contig in virusHits:
		newVirusHits[contig] = []
		leftIndex = 0
		rightIndex = 0
		if len(virusHits[contig]) > 1:
			virusHits[contig].sort(key = lambda hit: int(hit[0]))
			virusHit = virusHits[contig][0]
			left = int(virusHit[0])
			right = int(virusHit[1])
			for index in range(1, len(virusHits[contig])):
				virusHit = virusHits[contig][index]
				nextLeft = int(virusHit[0])
				if nextLeft <= right + 1:
					right =int(virusHit[1])
					rightIndex = index
				else:
					newVirusHits[contig].append((leftIndex, rightIndex))
					leftIndex = index
					rightIndex = index
					left = int(virusHit[0])
					right = int(virusHit[1])
		newVirusHits[contig].append((leftIndex, rightIndex))

	# Combine virus and aedes hit results to locate putative insertions.

	outFileBothEnds = open(outFilenameBothEnds, 'w')
	outFileOneEnd = open(outFilenameOneEnd, 'w')
	summaryBothEnds = 'Summary of hits:\n   '
	summaryOneEnd = 'Summary of hits:\n   '

	for contig in virusHits:
		if contig in aedesHits:
			for leftIndex, rightIndex in newVirusHits[contig]:       # partition aedes hits into before and after viral hit
				virus_qstart = int(virusHits[contig][leftIndex][0])
				virus_qend = int(virusHits[contig][rightIndex][1])
				before = []
				after = []
				for index in range(len(aedesHits[contig])):
					aedesHit = aedesHits[contig][index]
					aedes_qstart = int(aedesHit[0])
					aedes_qend = int(aedesHit[1])
					aedes_chr = aedesHit[3]
					aedes_sstart = int(aedesHit[4])
					aedes_send = int(aedesHit[5])
					if aedes_qend < virus_qstart + 15:
						before.append((index, aedes_sstart, aedes_send, aedes_chr))
					elif aedes_qstart > virus_qend - 15:
						after.append((index, aedes_sstart, aedes_send, aedes_chr))

				beforeHits = []      # find pairs of before/after aedes hits that are within 10Kbp of each other
				afterHits = []
				for b in before:
					for a in after:
						if (b[3] == a[3]) and (abs(b[2] - a[1]) < 10000):
							beforeHits.append(aedesHits[contig][b[0]])
							afterHits.append(aedesHits[contig][a[0]])

				if (len(beforeHits) > 0):                               # write hits with before&after to one file
					outFileBothEnds.write('Contig: ' + contig + '\n')
					outFileBothEnds.write('\n   Viral hit:' + '\n')
					for index in range(leftIndex, rightIndex + 1):
						virusHit = virusHits[contig][index]
						outFileBothEnds.write('      ' + ', '.join(virusHit[:2]+virusHit[3:5]+virusHit[6:]) + '\n')
						summaryBothEnds += virusHit[3] + '-' + virusHit[4]
						if index != rightIndex:
							summaryBothEnds += ':'
						else:
							summaryBothEnds += '\n   '
					beforeHits = list(set(beforeHits))
					beforeHits.sort(key = lambda hit: int(hit[0]))
					afterHits = list(set(afterHits))
					afterHits.sort(key = lambda hit: int(hit[0]))
					outFileBothEnds.write('\n   Aedes aegypti hits:' + '\n')

					outFileBothEnds.write('      Left:  ')
					spaces = 0
					for beforeHit in beforeHits:
						outFileBothEnds.write(' ' * spaces + ', '.join(beforeHit[:2]+beforeHit[3:6]+beforeHit[7:]) + '\n')
						spaces = 13
					outFileBothEnds.write('\n      Right: ')
					spaces = 0
					for afterHit in afterHits:
						outFileBothEnds.write(' ' * spaces + ', '.join(afterHit[:2]+afterHit[3:6]+afterHit[7:]) + '\n')
						spaces = 13
					outFileBothEnds.write('\n')
				else:                                                  # write less good hits to another file
					outFileOneEnd.write('Contig: ' + contig + '\n')
					outFileOneEnd.write('\n   Viral hit:' + '\n')
					for index in range(leftIndex, rightIndex + 1):
						virusHit = virusHits[contig][index]
						outFileOneEnd.write('      ' + ', '.join(virusHit[:2]+virusHit[3:5]+virusHit[6:]) + '\n')
						summaryOneEnd += virusHit[3] + '-' + virusHit[4]
						if index != rightIndex:
							summaryOneEnd += ':'
						else:
							summaryOneEnd += '\n   '
					outFileOneEnd.write('\n   Aedes aegypti hits (one side only):' + '\n')
					if len(before) > 0:
						outFileOneEnd.write('      Left:  ')
						spaces = 0
						for b in before:
							beforeHit = aedesHits[contig][b[0]]
							outFileOneEnd.write(' ' * spaces + ', '.join(beforeHit[:2]+beforeHit[3:6]+beforeHit[7:]) + '\n')
							spaces = 13
					if len(after) > 0:
						outFileOneEnd.write('\n      Right: ')
						spaces = 0
						for a in after:
							afterHit = aedesHits[contig][a[0]]
							outFileOneEnd.write(' ' * spaces + ', '.join(afterHit[:2]+afterHit[3:6]+afterHit[7:]) + '\n')
							spaces = 13
					outFileOneEnd.write('\n')
	outFileBothEnds.write(summaryBothEnds + '\n')
	outFileOneEnd.write(summaryOneEnd + '\n')
	outFileBothEnds.close()
	outFileOneEnd.close()

def consolidate(virusName):
	dir = Path('HITS_' + virusName).resolve()
	viralHits1 = {}
	viralHits2 = {}
	allHits = []
	for file in dir.iterdir():
		if '_viral_hits' in file.name:
			if 'hits1' in file.name:
				number = 1
			else:
				number = 2
			index = file.name.index('_viral_hits')
			specimen = file.name[:index]
			f = open(str(file), 'r')
			line = f.readline()
			while line:
				if line == 'Summary of hits:\n':
					# print(specimen, end = ', ')
					if number == 1:
						viralHits1[specimen] = []
					else:
						viralHits2[specimen] = []
					line = f.readline()
					line = line.strip()
					while line:
						if '-' in line and ':' not in line:
							l, r = line.split('-')
							if int(l) > int(r):
								hit = (int(r), int(l))
							else:
								hit = (int(l), int(r))
						elif ':' in line:
							parts = line.split('-')
							if int(parts[0]) > int(parts[-1]):
								parts = parts[::-1]
								hit = [int(parts[0])]
								for index in range(1, len(parts) - 1):
									r, l = parts[index].split(':')
									hit.append(int(l))
									hit.append(int(r))
								hit.append(int(parts[-1]))
								hit = tuple(hit)
							else:
								hit = [int(parts[0])]
								for index in range(1, len(parts) - 1):
									r, l = parts[index].split(':')
									hit.append(int(r))
									hit.append(int(l))
								hit.append(int(parts[-1]))
								hit = tuple(hit)
						# print(hit, end = ', ')
						if number == 1:
							viralHits1[specimen].append(hit)
						else:
							viralHits2[specimen].append(hit)
						if hit not in allHits:
							allHits.append(hit)
						line = f.readline()
						line = line.strip()
					# print()
					break
				line = f.readline()
			f.close()

	out = open(str(dir / ('results_' + virusName + '.txt')), 'w')
	allHits.sort()                        # write (start, end) header
	hitCounts1 = {}
	hitCounts2 = {}
	out.write('\t')
	for hit in allHits:
		out.write(str(hit) + '\t')
		hitCounts1[hit] = 0
		hitCounts2[hit] = 0
	out.write('\n')
	out.write('\t')                       # write length header
	for hit in allHits:
		length = sum([abs(hit[i+1] - hit[i]) + 1 for i in range(0, len(hit), 2)])
		out.write(str(length) + '\t')
	out.write('\n')
	specimens = list(viralHits1.keys())
	specimens.sort()
	for specimen in specimens:
		out.write(specimen + '\t')
		for hit in allHits:
			if hit in viralHits2[specimen]:
				out.write('2\t')
				hitCounts2[hit] += 1
			elif hit in viralHits1[specimen]:
				out.write('1\t')
				hitCounts1[hit] += 1
			else:
				out.write('0\t')
		out.write('\n')
	out.write('Total 2\t')
	for hit in allHits:
		out.write(str(hitCounts2[hit]) + '\t')
	out.write('\n')
	out.write('Total 1\t')
	for hit in allHits:
		out.write(str(hitCounts1[hit]) + '\t')
	out.write('\n')
	out.write('Total all\t')
	for hit in allHits:
		out.write(str(hitCounts1[hit] + hitCounts2[hit]) + '\t')
	out.write('\n')
	out.close()

def main():
	if (len(sys.argv) == 3) and (sys.argv[1] == '-c'):
		virusName = sys.argv[2]
		if not Path('HITS_' + virusName).exists():
			os.mkdir('HITS_' + virusName)
		os.system('cp -v */spades_' + virusName + '/*viral_hits* ' + 'HITS_' + virusName)
		consolidate(virusName)
	elif len(sys.argv) == 3:
		dir = Path(sys.argv[1])
		dir = dir.resolve()
		virusName = sys.argv[2]
		getHits(dir, virusName)
	elif len(sys.argv) == 2:
		dir = Path('.').resolve()
		virusName = sys.argv[1]
		for subdir in dir.iterdir():
			if subdir.is_dir() and ('test' not in str(subdir)) and (subdir / 'spades_' + virusName / 'scaffolds.fasta').exists():
				print(subdir)
				getHits(subdir, virusName)

if __name__ == '__main__':
	main()
