import csv
from operator import attrgetter
import pandas as pd
import time
import bisect
from bisect import bisect_left
import count_pvalue

class EXON:
	def __init__(self):
		self.st = 0
		self.en = 0

class Stack:
	def __init__(self):
		self.items = []

	def size(self):
		return len(self.items)

	def isEmpty(self):
		return self.items == []

	def push(self, val):
		self.items.append(val)

	def top(self):
		if self.isEmpty():
			return None
		else:
			return self.items[self.size()-1]

	def pop(self):
		if self.isEmpty():
			return None
		else:
			return self.items.pop()


def InsertIntoOldChromDict(ChromDict, chr, exListNew, geneName):
	GeneDict = ChromDict[chr]
	if geneName not in GeneDict.keys():
		GeneDict[geneName] = exListNew
	else:
		exList = GeneDict[geneName]
		exList.extend(exListNew)
		GeneDict[geneName] = exList
	return GeneDict


def bi_contains(lst, item):
    return bisect_left(lst, item)

def MergeIntervals(inputlist):
	n = len(inputlist)
	inputlist.sort(key = attrgetter('st'), reverse = False)

	st = Stack()
	st.push(inputlist[0])

	for i in range(1,n):
		stacktop = st.top()
		if inputlist[i].st <= stacktop.en:
			st.pop()
			stacktop.en = max(stacktop.en,inputlist[i].en)
			st.push(stacktop)
		else:
			st.push(inputlist[i])

	mergedExList = []
	while(True):
		if st.size() == 0:
			break;
		stacktop = st.top()
		mergedExList.append(stacktop)
		st.pop()

	return mergedExList



def CountTotalReadCoverage(chrom, exList, bam_list, position_row):
	totalCount = 0
	for p in range(len(exList)):
		start = int(exList[p].st)
		end = int(exList[p].en)

		pos1 = bi_contains(position_row, start)
		pos2 = bi_contains(position_row, end)

		#print("len",len(bam_list))
		#print("hehehe:",start, end, pos1, pos2)

		if(pos1 < len(bam_list) and pos2 < len(bam_list)):
			if(int(bam_list[pos2][1]) != end):
				pos2 = pos2 - 1
			#print(pos1, pos2)
			#print(int(bam_list[pos1][1]), int(bam_list[pos2][1]))

			for t in range(pos1, pos2+1):
				read = int(bam_list[t][2])
				totalCount += read
	#print("Total Count: ", totalCount)
		
	return totalCount

def writeResult(chrom, geneID, start, end, newList, bam_list, position_row, RC, mergedExListLength, writer):
	targetRC = CountTotalReadCoverage(chrom, newList, bam_list, position_row)
	targetLength = end - start + 1
	#### Avoiding divide by zero error ##### 
	if targetLength == 0:
		averageTargetRC = 0
	else:
		averageTargetRC = targetRC/targetLength

	if mergedExListLength==targetLength:
		averageRCexcludingTarget = 0
	else:
		averageRCexcludingTarget = (RC-targetRC)/(mergedExListLength-targetLength)

	writer.writerow([chrom, geneID, start, end, targetRC, targetLength, RC, mergedExListLength, RC-targetRC, mergedExListLength-targetLength,  averageTargetRC, averageRCexcludingTarget])

def callA3SS(chrom, geneID, bam_list, position_row, t_row, RC, mergedExListLength, writer):
	newList = []
	e = EXON()
	strand = t_row[7]
	if strand == '+':
		e.st = int(t_row[5])
		e.en = int(t_row[3])-1
	else:
		e.st = int(t_row[4])+1
		e.en = int(t_row[6])
	newList.append(e)
	writeResult(chrom, geneID, e.st, e.en, newList, bam_list, position_row, RC, mergedExListLength, writer)

def callA5SS(chrom, geneID, bam_list, position_row, t_row, RC, mergedExListLength, writer):
	newList = []
	e = EXON()
	strand = t_row[7]
	if strand == '+':
		e.st = int(t_row[4])+1
		e.en = int(t_row[6])
	else:
		e.st = int(t_row[5])-1
		e.en = int(t_row[4])
	newList.append(e)
	writeResult(chrom, geneID, e.st, e.en, newList, bam_list, position_row, RC, mergedExListLength, writer)

def callMXE(chrom, geneID, bam_list, position_row, t_row, RC, mergedExListLength, writer):
	newList = []
	e1 = EXON()
	e1.st = int(t_row[3])
	e1.en = int(t_row[4])
	newList.append(e1)
	writeResult(chrom, geneID, e1.st, e1.en, newList, bam_list, position_row, RC, mergedExListLength, writer)

	newList = []
	e2 = EXON()
	e2.st = int(t_row[5])
	e2.en = int(t_row[6])
	newList.append(e2)
	writeResult(chrom, geneID, e2.st, e2.en, newList, bam_list, position_row, RC, mergedExListLength, writer)

def callSE_RI(chrom, geneID, bam_list, position_row, t_row, RC, mergedExListLength, writer):
	newList = []
	e = EXON()
	e.st = int(t_row[3])
	e.en = int(t_row[4])
	newList.append(e)
	writeResult(chrom, geneID, e.st, e.en, newList, bam_list, position_row, RC, mergedExListLength, writer)


def Generate(ChromDict, chromosomes, filename, sample, pathout, pathin):
	tt = time.time()
	target_file_reader= open(filename+'.csv', "rt")
	target_read = csv.reader(target_file_reader, delimiter="\t")

	target_list = list(target_read)

	checkDict = {}
	position_row = []
	with open(pathout+filename+"_"+sample+".csv",'w') as f:
		writer = csv.writer(f, delimiter='\t')
		writer.writerow(['Chrom', 'Gene Name', 'Exon Start', 'Exon End', 'Target Read Count', 'Target Length', 'Others: Read Count', 'Others: Length', 'Others: Read Count- target RC', 'Others: Length - exon length', 'Average Target Read Count(n)', 'Average Read Count All(N)'])

		for chrom in chromosomes:
			count = 0
			bam_file_reader= open(pathin+chrom+".txt", "rt")
			bam_read = csv.reader(bam_file_reader, delimiter="\t")
			bam_list = list(bam_read)
			position_row = [int(bam_list[i][1]) for i in range(len(bam_list))]
			
			df = pd.DataFrame(target_list)
			target_tt = df.loc[df[0]==chrom]
			for t_row in target_tt.itertuples():
				geneID = t_row[2].strip().upper()
				#print(geneID)
				if geneID not in checkDict.keys():
					GeneDict = ChromDict[chrom]
					exList = GeneDict[geneID]

					mergedExList = MergeIntervals(exList)
					"""
					for pp in range(len(mergedExList)):
						print(mergedExList[pp].st, mergedExList[pp].en)
					"""
					mergedExListLength = 0
					for p in range(len(mergedExList)):
						mergedExListLength += mergedExList[p].en - mergedExList[p].st + 1

					RC = CountTotalReadCoverage(chrom, mergedExList, bam_list, position_row)
					checkDict[geneID] = RC
				else:
					RC = checkDict[geneID]

				if(filename == 'SE' or filename == 'RI'):
					callSE_RI(chrom, geneID, bam_list, position_row, t_row, RC, mergedExListLength, writer)
				elif(filename == 'MXE'):				
					callMXE(chrom, geneID, bam_list, position_row, t_row, RC, mergedExListLength, writer)
				elif(filename == 'A5SS'):
					callA5SS(chrom, geneID, bam_list, position_row, t_row, RC, mergedExListLength, writer)
				elif(filename == 'A3SS'):
					callA3SS(chrom, geneID, bam_list, position_row, t_row, RC, mergedExListLength, writer)

	print("Elapsed time: ",(time.time()-tt)/60)
	f.close()


def MakeFullDictionary(ann_list, chromosomes):
	ChromDict = {}
	
	for a_row in ann_list:
		chrom = a_row[2].strip()
		geneName = a_row[1].strip().upper()

		exonCount = int(a_row[8])
		exonStartList = a_row[9].split(',')
		exonEndList = a_row[10].split(',')

		exList = []
		for i in range(exonCount):
			exonStart = int(exonStartList[i])
			exonEnd = int(exonEndList[i])
			newExon = EXON()
			newExon.st = exonStart
			newExon.en = exonEnd
			exList.append(newExon)
			
		if chrom in chromosomes:
			if chrom not in ChromDict.keys():
				GeneDict = {}
				GeneDict[geneName] = exList
			else:
				GeneDict = InsertIntoOldChromDict(ChromDict, chrom, exList, geneName)

			ChromDict[chrom] = GeneDict
			GeneDict_x = ChromDict[chrom]

	return ChromDict