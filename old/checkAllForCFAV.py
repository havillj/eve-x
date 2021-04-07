#!/usr/bin/env python3

from pathlib import Path
import pysam
import os

def main():
	print('specimen,323-3780,3914-4254,6313-6369,7311-7486,8549-8737,other')
	cwd = Path('.')
	for dir in cwd.iterdir():
		if not dir.is_dir():
			continue
		subdir1 = dir / 'gd2'
		subdir2 = dir / 'gd'
		if subdir1.exists():
			subdir = subdir1
		elif subdir2.exists():
			subdir = subdir2
		else:
			continue
		bamFile = subdir / 'consensus-BAM-files' / 'NC_001564.2_1' / 'consensus_alignment_sorted.bam'
		print(dir.name + ', ', end = '')
		if bamFile.exists():
#			print(str(bamFile))
			if not Path(str(bamFile) + '.bai').exists():
#				pysam.index(str(bamFile))
				os.system('samtools index ' + str(bamFile))
			count = pysam.view('-c', str(bamFile), 'NC_001564.2_1_consensus:0')
			counts = []
			counts.append(pysam.view('-c', str(bamFile), 'NC_001564.2_1_consensus:323-3780'))
			counts.append(pysam.view('-c', str(bamFile), 'NC_001564.2_1_consensus:3914-4254'))
			counts.append(pysam.view('-c', str(bamFile), 'NC_001564.2_1_consensus:6313-6369'))
			counts.append(pysam.view('-c', str(bamFile), 'NC_001564.2_1_consensus:7311-7486'))
			counts.append(pysam.view('-c', str(bamFile), 'NC_001564.2_1_consensus:8549-8737'))
			other = int(count.rstrip()) - sum([int(c.rstrip()) for c in counts])
			print(', '.join([c.rstrip() for c in counts]) + ', ' + str(other))
		else:
			print('0, 0, 0, 0, 0, 0')
	
main()
