import os
from pathlib import Path

def unzipFiles():
	cwd = Path('.')
	for file in cwd.iterdir():
		if file.is_dir():
			continue
		fn = str(file)
		if fn[-4:] == '.zip':
			os.system('unzip ' + fn)
	
	os.system('mkdir consensus-BAM-files')
	os.system('mv NC* consensus-BAM-files')
	os.system('rm *.zip')

def main():
	unzipFiles()
	
main()
