from pathlib import Path
import os

# p = Path('.')
# for f in p.iterdir():
# 	try:
# 		index = f.name.index('.sorted.deduped.merged.bam')
# 	except:
# 		continue
# 	dir = f.name[:index]
# 	Path(dir).mkdir()
	
p = Path('.')
for f in p.iterdir():
	if f.is_dir():
		continue
	try:
		index = f.name.index('.sorted.deduped.merged.bam')
	except:
		index = f.name.index('_unmapped_mates.fastq')
	os.rename(f.name, f.name[:index] + '/' + f.name)
