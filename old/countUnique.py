#!/usr/bin/env python3

import sys

def countUnique(filenameCSV):
	
	csv = open(filenameCSV, 'r')
	reads = {}
	for line in csv:
		cols = line.split(',')
		query = cols[0]
		reads[query] = ''
	csv.close()
	
	return len(reads)
	
def main():
	print(countUnique(sys.argv[1]))
	
main()
