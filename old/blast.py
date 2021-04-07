import os
import sys

os.system('echo ' + sys.argv[2]+ '| /usr/local/bin/blastn -task blastn -subject /home/havill/data/aegypti/genomes/' + sys.argv[1] + '.fasta')
