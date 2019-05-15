import csv
from operator import attrgetter
import pandas as pd
import time
from scipy.stats import chisquare

def Count_pvalue(filename, output_path, control, case):
	S1output_reader= open(output_path+filename+"_"+control+".csv", "rt", encoding='ascii')
	S1out = csv.reader(S1output_reader, delimiter="\t")
	S2output_reader= open(output_path+filename+"_"+case+".csv", "rt", encoding='ascii')
	S2out = csv.reader(S2output_reader, delimiter="\t")

	S1_list = list(S1out)
	S2_list = list(S2out)

	with open(output_path+filename+"_"+control+"_Vs_"+case+".csv",'w') as f:
		writer = csv.writer(f, dialect='excel',delimiter='\t')
		writer.writerow(['Chrom', 'Gene Name', 'Exon Start', 'Exon End', 'p-value', 'Ratio difference', 'Absolute Ratio difference', 'Chrom region', 'Target RC '+control, 'Target RC '+case, 'Target Exon Length', 'Others: Read Count '+control, 'Others: Length '+control, 'Others: Read Count '+case, 'Others: Length '+case, 'Average read coverage '+control, 'Average read coverage '+case, 'Chrom region Long'])
		for i in range(1, len(S1_list)):
			n1 = float(S1_list[i][10])
			n2 = float(S2_list[i][10])
			N1 = float(S1_list[i][11])
			N2 = float(S2_list[i][11])
			print(n1, n2, N1, N2)

			if N1!=0 and N2!=0:
				ratio_diff = (n1/N1) - (n2/N2)
			else:
				ratio_diff = 0
			
			abs_ratio_diff = abs(ratio_diff)
			if(abs_ratio_diff>=0.1):
				N1 = N1 + n1
				N2 = N2 + n2
				if (N1+N2)>0:
					P0 = (n1+n2)/(N1+N2)
				n10 = N1 * P0
				n20 = N2 * P0

				res = chisquare([n1, N1-n1, n2, N2-n2], f_exp=[n10, N1-n10, n20, N2-n20], ddof = 1)
				chrom_region_long = S1_list[i][0]+':'+S1_list[i][1]+":"+S1_list[i][2]+"-"+S1_list[i][3]
				chrom_region = S1_list[i][0]+':'+S1_list[i][2]+"-"+S1_list[i][3]
				writer.writerow([S1_list[i][0], S1_list[i][1], S1_list[i][2], S1_list[i][3], res[1], ratio_diff, abs_ratio_diff, chrom_region, S1_list[i][4], S2_list[i][4], S1_list[i][5], S1_list[i][6], S1_list[i][7], S2_list[i][6], S2_list[i][7], n1, n2, chrom_region_long])

	f.close()

