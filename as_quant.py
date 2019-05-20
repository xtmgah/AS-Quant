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
import sys

startTime = time.time()

chromosomes_h = ['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21', 'chr22','chrX','chrY']
chromosomes_m = ['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chrX','chrY']
target_filenames = ['SE', 'RI', 'MXE', 'A3SS', 'A5SS']

if(len(sys.argv)<5):
	print("Please provide all of the mandatory arguments. Example: $ python3 as_quant.py -s human s1 s2")
	sys.exit()


if len(sys.argv)==7:
	for ii in range(len(sys.argv)):
		if sys.argv[ii] == '-s' or sys.argv[ii] == '-S':
			species = sys.argv[ii+1]
		if sys.argv[ii] == '-o' or sys.argv[ii] == '-O':
			output_dir = sys.argv[ii+1]+'/'

elif len(sys.argv)==5:
	species = sys.argv[2]
	output_dir = 'Output/'
os.makedirs(output_dir, exist_ok=True)

input1_dir = sys.argv[-2]
input2_dir = sys.argv[-1]
print(output_dir, input1_dir, input2_dir)

input1 = input1_dir.split('/')
input2 = input2_dir.split('/')
s1_name = input1[-2:]
s2_name = input2[-2:]

s1_dir = ""
for i in range(len(input1)-2):
	s1_dir+=input1[i]+'/'

s2_dir = ""
for i in range(len(input1)-2):
	s2_dir+=input2[i]+'/'


samplenames = [s1_name[0], s2_name[0]]
bamfile_names = [s1_name[1], s2_name[1]]

#print("samplename and bamfiles: ",samplenames, bamfile_names)
#print("Prefix dir s1_dir s2_dir: ",s1_dir, s2_dir)

if species =='human':
	chromosomes = chromosomes_h
	inp = 'homo_sapiens/'
elif species == 'mouse':
	chromosomes = chromosomes_m
	inp = 'mus_musculus/'

for i in range(2):
	if i==0:
		input_path = s1_dir
	else:
		input_path = s2_dir
	preprocess.SamtoText(input_path, samplenames[i], bamfile_names[i], chromosomes)


ann_file_reader= open(inp+'annotation.csv', "rt")
ann_read = csv.reader(ann_file_reader, delimiter="\t")
ann_list = list(ann_read)

ChromDict = methods.MakeFullDictionary(ann_list, chromosomes)

# call generate method 5 times for 5 different target exons lists and with two input samples for each case
for filename in target_filenames:
	print("Taget Exon type: ",filename)
	for i in range(2):
		if i==0:
			input_path = s1_dir
		else:
			input_path = s2_dir
		pathin = input_path+samplenames[i]+'/'
		print("Sample name: ", samplenames[i])
		methods.Generate(ChromDict, chromosomes, filename, samplenames[i], inp, pathin, output_dir)

for filename in target_filenames:
	count_pvalue.Count_pvalue(filename, output_dir, samplenames[0], samplenames[1])

totalTime = time.time() - startTime
print("Total program time is : ",totalTime)

# python3 as_quant.py -s human /home/fahmi/Desktop/AS-Quant/mousemm10/input/s1/accepted_hits.bam /home/fahmi/Desktop/AS-Quant/mousemm10/input/s2/hits.bam