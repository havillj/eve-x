import sys

def fq2fa(fn):
	fq = open(fn, 'r')
	fa = open(fn[:-6] + '.fasta', 'w')
	
	desc = fq.readline()
	while desc != '':
		seq = fq.readline()
		plus = fq.readline()
		qual = fq.readline()
		
		fa.write('>' + desc[1:])
		fa.write(seq)
		
		desc = fq.readline()
		
	fq.close()
	fa.close()
	
def main():
	if len(sys.argv) > 1:
		fq2fa(sys.argv[1])
		
main()
