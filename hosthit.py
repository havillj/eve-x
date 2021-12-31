# qseqid qstart qend qseq sseqid sstart send sseq evalue bitscore
    
H_QSEQID = 0    # contig id
H_QSTART = 1    # start pos in contig
H_QEND = 2      # end pos in contig
H_QSEQ = 3      # sequence in contig aligned to sequence in host genome
H_SSEQID = 4    # accession number of host genome hit
H_SSTART = 5    # start pos in host genome
H_SEND = 6      # end pos in host genome
H_SSEQ = 7      # sequence in host genome aligned to sequence in contig
H_EVALUE = 8    # e-value of alignment
H_BITSCORE = 9  # bitscore of alignment

class HostHit:
    
    def __init__(self, data):
        self.contig = data[H_QSEQID]
        self.qstart = int(data[H_QSTART])
        self.qend = int(data[H_QEND])
        self.qseq  = data[H_QSEQ]
        self.sseqid = data[H_SSEQID]
        self.sstart = int(data[H_SSTART])
        self.send = int(data[H_SEND])
        self.sseq = data[H_SSEQ]
        self.evalue = float(data[H_EVALUE])
        self.bitscore = float(data[H_BITSCORE])
