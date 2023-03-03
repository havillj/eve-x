# qseqid qstart qend qseq sstart send sseq evalue bitscore sseqid stitle pident

QSEQID = 0    # contig id
QSTART = 1    # start pos in contig
QEND = 2      # end pos in contig
QSEQ = 3      # sequence in contig
SSTART = 4    # start pos in virus genome
SEND = 5      # end pos in virus genome
SSEQ = 6      # sequence in virus genome
EVALUE = 7    # e-value of alignment
BITSCORE = 8  # bitscore of alignment
SSEQID = 9    # accession number of viral genome hit
STITLE = 10   # species name of viral genome hit
PIDENT = 11   # % identities in alignment

class VirusHit:

    def __init__(self, data):
        self.contig = data[QSEQID]
        self.qstart = int(data[QSTART])
        self.qend = int(data[QEND])
        self.qseq = data[QSEQ]
        self.sstart = int(data[SSTART])
        self.send = int(data[SEND])
        self.sseq = data[SSEQ]
        self.evalue = float(data[EVALUE])
        self.bitscore = float(data[BITSCORE])
        self.sseqid = data[SSEQID]
        self.stitle = data[STITLE]
        self.pident = float(data[PIDENT])
        
        assert self.qstart < self.qend
        self.qlength = self.qend - self.qstart + 1
        self.slength = abs(self.send - self.sstart) + 1
