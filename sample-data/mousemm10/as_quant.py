import csv
import os
from operator import attrgetter
import pandas as pd
import time
import bisect
import count_pvalue
import initial
import methods
import preprocess

startTime = time.time()


chromosomes = ['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chrX','chrY']
target_filenames = ['SE', 'RI', 'MXE', 'A3SS', 'A5SS']
samplenames = ['s1', 's2']

current_path = '/home/fahmi/Desktop/AS-Quant/mousemm10/'	# the main directory
input_path = current_path+'input/'			# input directory
output_path = current_path+'output/'					# output directory

os.makedirs(output_path, exist_ok=True)


for sample in samplenames:
	print("Running from as-quant")
	preprocess.SamtoText(input_path, sample, chromosomes)

ann_file_reader= open(current_path+'annotation.csv', "rt")
ann_read = csv.reader(ann_file_reader, delimiter="\t")
ann_list = list(ann_read)

ChromDict = methods.MakeFullDictionary(ann_list, chromosomes)

# call generate method 5 times for 5 different targets and with 2 samples

for filename in target_filenames:
	print("Filename: ",filename)
	for sample in samplenames:
		pathin = input_path+sample+'/'
		print("Sample name: ", sample)
		methods.Generate(ChromDict, chromosomes, filename, sample, output_path, pathin)


for filename in target_filenames:
	count_pvalue.Count_pvalue(filename, output_path, samplenames[0], samplenames[1])

totalTime = time.time() - startTime
print("Total program time is : ",totalTime)
