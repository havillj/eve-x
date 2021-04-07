import sys

def split(fn):
	fa = open(fn, 'r')
	index = fn.index('.fa')
	fa1 = open(fn[:index] + '1' + fn[index:], 'w')
	fa2 = open(fn[:index] + '2' + fn[index:], 'w')
	fa3 = None
	
	desc = fa.readline()
	while desc != '':
		seq = fa.readline()
		if desc[-3:-1] == '/1':
			fa1.write(desc[:-3] + '\n')
			fa1.write(seq)
		elif desc[-3:-1] == '/2':
			fa2.write(desc[:-3] + '\n')
			fa2.write(seq)
		else:
			if fa3 is None:
				fa3 = open(fn[:index] + '0' + fn[index:], 'w')
			fa3.write(desc)
			fa3.write(seq)
		
		desc = fa.readline()
		
	fa.close()
	fa1.close()
	fa2.close()
	if fa3 is not None:
		fa3.close()
	
def main():
	if len(sys.argv) > 1:
		split(sys.argv[1])
		
main()
