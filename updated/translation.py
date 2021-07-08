from Bio import SeqIO, Entrez, AlignIO
from Bio.Seq import Seq
from Bio.Align import MultipleSeqAlignment
from Bio.SeqFeature import SeqFeature, FeatureLocation
import os

SEQUENCES_DIR = "/Volumes/Data2/results/viruses/sequences"


#path = "/Volumes/Data2/specimens/Amacuzac-Mexico-13.LIN210A1627/results/sequences/Amacuzac-Mexico-13.LIN210A1627_hits_unaligned.fasta"

def getTranslation():
    hitsFiles = []
    paths = os.listdir(SEQUENCES_DIR)

    for fam in paths:
        path = f"/Volumes/Data2/results/viruses/sequences/{fam}/sequences_{fam}_per_contig_unaligned.fasta"
        outPath = f"/Volumes/Data2/results/viruses/sequences/{fam}/sequences_{fam}_per_contig_translated.fasta"
        with open(path, 'r') as reader, open(outPath, 'w') as outfile:        #count = 1
            for line in reader:
                if line[0] == ">":
                    startPos = []
                    endPos = []
                    title = line.split(".")[0] + "." + line.split(".")[1][0]
                    for x in line.split("_")[-1].split(","):    #get start and end positions in DNA
                        startPos.append(int(x.split("-")[0]))
                        endPos.append(int(x.split("-")[1]))
                else:
                    seq = Seq(line.strip().replace('-', 'N'))
                    rev = seq.reverse_complement()

                    for frame in range(3):
                        frameStart = []
                        end = []
                        revStart = []
                        revEnd = []
                        posStr = ""
                        revPosStr = ""
                        frameEnd = len(seq) - len(seq[frame:])%3

                        for x in range(len(startPos)):
                            frameStart.append(startPos[x] + frame)
                            end.append(endPos[x] - len(seq[frame:])%3)
                            posStr = posStr + str(frameStart[x]) + "-" + str(end[x]) + ","

                            revStart.append(endPos[x] - frame)
                            revEnd.append(startPos[x] + len(rev[frame:])%3)
                            revPosStr = revPosStr + str(revStart[x]) + "-" + str(revEnd[x]) + ","

                        tempSeq = seq[frame:frameEnd]
                        outfile.write(title + "_" + posStr[:-1] + "_(+)\n")
                        outfile.write(str(tempSeq.translate()) + "\n")

                        revSeq = rev[frame:frameEnd]
                        outfile.write(title + "_" + revPosStr[:-1] + "_(-)\n")
                        outfile.write(str(revSeq.translate()) + "\n")

getTranslation()
