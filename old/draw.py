from Bio.Graphics import GenomeDiagram
from Bio.SeqFeature import SeqFeature, FeatureLocation
from reportlab.lib import colors
from reportlab.lib.units import cm
import sys
from pathlib import Path

INTERVAL = 1000

def addFeature(gd_diagram, index, read):
	track = gd_diagram.get_tracks()[index]
	feature_set = track.get_sets()[0]

	feature = SeqFeature(FeatureLocation(read[0], read[1]), strand = +1)
	arrow_color = colors.blue
	
	feature_set.add_feature(feature, color = arrow_color, sigil = 'BOX', name = read[2][:12], label = 1, label_size = 6)
	
def drawChromosomes(dir):
	path = Path(dir).resolve().absolute()
	filenames = [f.name for f in path.iterdir() if f.name[-4:].lower() == '.csv']
	filenames.sort()
	
	chrLengths = (310827022, 474425716, 409777670)
	chrs =[GenomeDiagram.Diagram('Chromosome ' + str(index + 1), circular = 1) for index in range(3)]
	chrMap = {'NC_035107.1': 0, 'NC_035108.1': 1, 'NC_035109.1': 2}
	
	counts = [{}, {}, {}]
			
	aedesIndex = 0
	for filename in filenames:
		f = open(path.name + '/' + filename, 'r')
		header = f.readline()
		
		for chrNum in range(3):
			track = chrs[chrNum].new_track(aedesIndex + 1, name = filename[:-4], greytrack = 0, greytrack_labels = 0, scale_format = 'SInt', scale_largetick_interval = 100000000, scale_color = colors.green, scale_fontsize = 6)
			fs = track.new_set()
		
		for line in f:
			cols = line[:-1].split(',')
			virus = cols[3]
			if len(cols) == 7:
				length = len(cols[4])  
				chrAcc = cols[5]       
				pos = int(cols[6])     
			else:  # len(cols) == 10
				length = len(cols[6])  # 4
				chrAcc = cols[7]       # 5
				pos = int(cols[8])     # 6

			if chrAcc in chrMap:
				chrIndex = chrMap[chrAcc]
				addFeature(chrs[chrIndex], aedesIndex, (pos, pos + length - 1, virus))

				p = pos + length // 2
				p = int(round(p / INTERVAL) * INTERVAL)
				if p not in counts[chrIndex]:
					counts[chrIndex][p] = [aedesIndex]
				else:
					if aedesIndex not in counts[chrIndex][p]:
						counts[chrIndex][p].append(aedesIndex)
		f.close()
		aedesIndex += 1
	
	for chrIndex in range(len(chrs)):
		graph = chrs[chrIndex].new_track(aedesIndex + 1, axis_labels = 0, name = 'Counts', scale_format = 'SInt', scale_largetick_interval = 100000000, scale_color = colors.blue, scale_fontsize = 6)
		set = graph.new_set('graph')
		data = []
		for p in range(0, chrLengths[chrIndex], INTERVAL):
			if p in counts[chrIndex]:
				data.append((p, len(counts[chrIndex][p])))
			else:
				data.append((p, 0))
		set.new_graph(data, 'Counts', style = 'bar', color=colors.blue, center = 0)
	
		chrs[chrIndex].draw(format = 'circular', orientation = 'portrait', pagesize = (2*20.32*cm, 2*27.94*cm), start = 0, end = chrLengths[chrIndex])
		chrs[chrIndex].write(path.name + '/chr' + str(chrIndex + 1) + '.pdf', 'PDF')

def addFeatureAll(gd_diagram, index, read):
	track = gd_diagram.get_tracks()[index]
	feature_set = track.get_sets()[0]

	feature = SeqFeature(FeatureLocation(read[0], read[1]), strand = +1)
	arrow_color = colors.blue
	
	feature_set.add_feature(feature, color = arrow_color, sigil = 'BOX', name = read[2][:3], label = 1, label_size = 4)

def drawChromosomesAll(dir):
	path = Path(dir).resolve().absolute()
	filenames = [f.name for f in path.iterdir() if f.name[-4:].lower() == '.csv']
	filenames.sort()
	
	chrLengths = (310827022, 474425716, 409777670)
	chrs =[GenomeDiagram.Diagram('Chromosome ' + str(index + 1), circular = 1) for index in range(3)]
	chrMap = {'NC_035107.1': 0, 'NC_035108.1': 1, 'NC_035109.1': 2}
	
	counts = [{}, {}, {}]
			
	aedesIndex = 0
	for filename in filenames:
		f = open(path.name + '/' + filename, 'r')
		header = f.readline()
		
		locations = ['Mexico' , 'Thailand', 'Angola', 'Argentina', 'Gabon', 'USA', 'Brazil']
		locIndex = 0
		for location in locations:
			if location in filename:
				break
			locIndex += 1
		color = colors.Color(0, locIndex * 0.11, 0, 1)
		
		for chrNum in range(3):
			track = chrs[chrNum].new_track(aedesIndex + 1, name = locations[locIndex][:4], greytrack = True, greytrack_labels = 3, greytrack_fontsize = 4, scale_ticks = 0, scale_format = 'SInt', scale_largetick_interval = 100000000, scale_color = color, scale_fontsize = 4)
			fs = track.new_set()
		
		for line in f:
			cols = line[:-1].split(',')
			virus = cols[3]
			if len(cols) == 7:
				length = len(cols[4])  
				chrAcc = cols[5]       
				pos = int(cols[6])     
			else:  # len(cols) == 10
				length = len(cols[6])  # 4
				chrAcc = cols[7]       # 5
				pos = int(cols[8])     # 6

			if chrAcc in chrMap:
				chrIndex = chrMap[chrAcc]
				addFeatureAll(chrs[chrIndex], aedesIndex, (pos, pos + length - 1, virus))

				p = pos + length // 2
				p = int(round(p / INTERVAL) * INTERVAL)
				if p not in counts[chrIndex]:
					counts[chrIndex][p] = [aedesIndex]
				else:
					if aedesIndex not in counts[chrIndex][p]:
						counts[chrIndex][p].append(aedesIndex)
		f.close()
		aedesIndex += 1
	
	for chrIndex in range(len(chrs)):
		graph = chrs[chrIndex].new_track(aedesIndex + 1, axis_labels = 0, name = 'Counts', scale_format = 'SInt', scale_largetick_interval = 100000000, scale_color = colors.blue, scale_fontsize = 4)
		set = graph.new_set('graph')
		data = []
		for p in range(0, chrLengths[chrIndex], INTERVAL):
			if p in counts[chrIndex]:
				data.append((p, 2 * len(counts[chrIndex][p])))
			else:
				data.append((p, 0))
		set.new_graph(data, 'Counts', style = 'bar', color=colors.blue, center = 0)
	
		chrs[chrIndex].draw(format = 'circular', orientation = 'portrait', pagesize = (2*20.32*cm, 2*27.94*cm), start = 0, end = chrLengths[chrIndex])
		chrs[chrIndex].write(path.name + '/chr' + str(chrIndex + 1) + '.pdf', 'PDF')

def main():
	if sys.argv[1] == '-a':
		dir = sys.argv[2]
		drawChromosomesAll(dir)
	else:
		dir = sys.argv[1]
		drawChromosomes(dir)

main()
