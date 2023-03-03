from Bio import SeqIO

fileName = 'rna_viruses.fasta'
retroCSV = 'retroviruses.csv'
outFileName = 'rna_viruses_no_retro.fasta'

retroFile = open(retroCSV, 'r')
retroNames = []
for line in retroFile:
    retroNames.append(line.rstrip())
retroFile.close()

records = list(SeqIO.parse(fileName, 'fasta'))

keepRecords = []
for r in records:
    species = r.description.split('|')[1]
    if species not in retroNames:
        keepRecords.append(r)

outputFile = open(outFileName, 'w')
SeqIO.write(keepRecords, outputFile, 'fasta')

print(len(records) - len(keepRecords))
