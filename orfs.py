import sys
sys.path.append('/home/havill/bin/modules')
import compbio
import orf

dna = compbio.readFASTA1(sys.argv[1])
orfs = orf.orf(dna, 20)
orf.printORFs(orfs)
for o in orfs:
	if o[0] < o[1]:
		print('(' + str(o[0]) + ', ' + str(o[1]) + ', 1),')
	else:
		print('(' + str(o[1]) + ', ' + str(o[0]) + ', -1),')
