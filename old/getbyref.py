import sys
import pysam

def getReads(filenameBAM):
	
	bamfile = pysam.AlignmentFile(filenameBAM, 'rb', threads = 8)
	
	index = filenameBAM.index('.sorted')
	newBAM = pysam.AlignmentFile(filenameBAM[:index] + '_scaff.bam', 'wb', template = bamfile)
	
	genome = ['NC_035107.1', 'NC_035108.1', 'NC_035109.1', 'NC_035159.1']
	for read in bamfile:
		id = read.reference_id
		refname = bamfile.get_reference_name(id)
		if not read.is_unmapped and refname[:2] == 'NW':  # not in genome:
			newBAM.write(read)
		
	bamfile.close()
	newBAM.close()
	
def main():
	getReads(sys.argv[1])
	
main()
