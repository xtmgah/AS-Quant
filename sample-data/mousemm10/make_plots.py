import matplotlib
matplotlib.use("agg")
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np
import pandas as pd
import bisect
from bisect import bisect_left
import csv
import methods
import os

y_limit = 0

def bi_contains(lst, item):
    return bisect_left(lst, item)

def Generate_read_coverate_plot(ax, pathin, sample, chrom, geneID, start, end, startAll, endAll):
	bam_file_reader= open(pathin+sample+'/'+chrom+".txt", "rt")
	bam_read = csv.reader(bam_file_reader, delimiter="\t")
	bam_list = list(bam_read)
	position_row = [int(bam_list[i][1]) for i in range(len(bam_list))]

	ax = ax or plt.gca()
	x_formatter = matplotlib.ticker.ScalarFormatter(useOffset=False)
	x_formatter.set_scientific(False)
	ax.xaxis.set_major_formatter(x_formatter)
	ax.tick_params(axis='both', which='major', labelsize=7)

	pos1 = bi_contains(position_row, startAll)
	pos2 = bi_contains(position_row, endAll)
	if(int(bam_list[pos2][1]) != endAll):
		pos2 = pos2 - 1

	p = []
	c = []
	read = 0
	length = endAll - startAll + 1
	for t in range(length):
		p.append(t+startAll)
		c.append(0)
		
	for t in range(pos1, pos2+1):
		position = int(bam_list[t][1])
		read = int(bam_list[t][2])
		index = p.index(position)
		c[index] = read

	p = np.array(p)
	c = np.array(c)

	pos3 = bi_contains(p,start)
	pos4 = bi_contains(p,end)

	global y_limit
	m = max(c)
	if m > y_limit:
		y_limit = m

	caption = ax.fill_between(p,c, color="skyblue", alpha=0.9, label = sample)
	ax.fill_between(p[pos3:pos4+1],c[pos3:pos4+1], color="red")
	ax.legend(handles = [caption])
	ax.set_xlim(startAll, endAll)
	ax.autoscale(enable = True)

	return y_limit

def Generate_annotation_plot(ax, isoforms, exonCountList, exonStartList, exonEndList, start, end, startAll, endAll):
	ax = ax or plt.gca()
	ax.autoscale(enable = True)
	x_formatter = matplotlib.ticker.ScalarFormatter(useOffset=False)
	x_formatter.set_scientific(False)
	ax.xaxis.set_major_formatter(x_formatter)
	print("isoforms are: ",isoforms)
	ystart = 0
	height = 3
	for i in range(isoforms):
		if i>=15:
			print("15 isoforms of this Gene is plotted.")
			break;
		else:
			ax.hlines(y=(ystart+ystart+height)/2, xmin=startAll, xmax=endAll, linewidth=1, color='blue', alpha = 0.5, linestyle = '--')
			ecount = int(exonCountList[i])
			stList = exonStartList[i]
			enList = exonEndList[i]
			for p in range(ecount):
				ex_s = int(stList[p])
				width = int(enList[p]) - int(stList[p]) + 1
				#print(start,end, ex_s, int(enList[p]))
				
				rect = patches.Rectangle((ex_s,ystart), width, height, color = 'skyblue', alpha=0.9, fill = True)
				ax.add_patch(rect)
				
				if ((start >= ex_s) and (end <= int(enList[p]))):
					width = end - start + 1
					rect1 = patches.Rectangle((start,ystart), width, height, color = 'red', fill = True)
					ax.add_patch(rect1)

			ystart +=5

	ax.set_xlim(startAll, endAll)
	ax.autoscale(enable = True)
	ax.spines["top"].set_visible(False)
	ax.spines["right"].set_visible(False)
	ax.spines["left"].set_visible(False)
	ax.set_yticklabels([])
	ax.tick_params(left=False, axis='both', which='major', labelsize=7)

	return

def Take_user_inputs(samplenames, pathin, ChromDict):
	ChromDict = methods.MakeFullDictionary(ann_list, chromosomes)

	region = input("Enter the range: (chr:Gene:Start-End)	")
	#print("You entered " + str(region))
	chrom, geneID, rng = region.split(':')
	start, end = rng.split('-')
	print("Chromosome: ",chrom)
	print("Gene Name: ", geneID)
	print("Region Start: ", start)
	print("Region End: ", end)

	GeneDict = ChromDict[chrom]
	exList = GeneDict[geneID.upper()]

	mergedExList = methods.MergeIntervals(exList)
	startAll = int(exList[0].st)
	endAll = int(exList[-1].en)

	df = pd.DataFrame(ann_list)
	ann_tt = df.loc[df[2]==chrom]

	exonStartList = {}
	exonEndList = {}
	exonCountList = {}

	isoforms = 0
	for a_row in ann_tt.itertuples():
		geneName = a_row[2].strip()
		if geneName.upper() == geneID.upper():
			exonCount = int(a_row[9])
			exonCountList[isoforms] = exonCount
			exonStartList[isoforms] = a_row[10].split(',')
			exonEndList[isoforms] = a_row[11].split(',')
			isoforms+=1

	y_axis_height = 20
	title = ""+chrom+":"+start+"-"+end+"("+geneID+")"

	fig = plt.figure()

	ax1 = fig.add_subplot(3,1,1)
	ax1.set_title(title, color = "black")
	ax1.set_ylabel('Counts')
	y_limit = Generate_read_coverate_plot(ax1, pathin, samplenames[0], chrom, geneID, int(start), int(end), startAll, endAll)

	ax2 = fig.add_subplot(3,1,2)
	ax2.set_xlabel('Position')
	ax2.set_ylabel('Counts')
	y_limit = Generate_read_coverate_plot(ax2, pathin, samplenames[1], chrom, geneID, int(start), int(end), startAll, endAll)

	# set ylim for both axes after getting max ylim values
	ax1.set_ylim(0, y_limit)
	ax2.set_ylim(0, y_limit)
	
	ax3 = fig.add_subplot(3,1,3)
	ax3.set_ylabel('Isoforms')
	Generate_annotation_plot(ax3, isoforms, exonCountList, exonStartList, exonEndList, int(start), int(end), startAll, endAll)

	plotout = path+'Plots/'
	os.makedirs(plotout, exist_ok=True)
	plt.savefig(plotout+title+'.png')


######### Main starts here #################


path = '/home/fahmi/Desktop/AS-Quant/mousemm10/'			# the main directory
pathin = path+'input/'

chromosomes = ['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chrX','chrY']
samplenames = ['s1', 's2']

ann_file_reader= open(path+'annotation.csv', "rt")
ann_read = csv.reader(ann_file_reader, delimiter="\t")
ann_list = list(ann_read)

Take_user_inputs(samplenames, pathin, ann_list)

"""
Test input: 

chr10:TMEM254:81845986-81846014

chr8:EPHX2:27358480-27358527

chr20:RBM39:34328446-34328519

chr4:HNRPDL:83346715-83346820

"""

