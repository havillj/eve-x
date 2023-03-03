from Bio import SeqIO

fileName = 'newdb.fasta'
retro = '../pararnavirae_plus_unverified.acc'
outFileName = 'newdb_no_retro_or_unverified.fasta'

retroFile = open(retro, 'r')
retroAccs = {}
for line in retroFile:
    acc = line.rstrip()
    if '.' in acc:
        acc = acc[:-2]
    retroAccs[acc] = 0
retroFile.close()
        
keepRecords = (r for r in SeqIO.parse(fileName, 'fasta') if r.description.split(' |')[0].split('.')[0] not in retroAccs)

SeqIO.write(keepRecords, outFileName, 'fasta')
