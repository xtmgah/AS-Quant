import os
import time


def SamtoText(path, sample, chromosomes):
	out_directory = path+sample+"/"
	os.makedirs(out_directory, exist_ok=True)

	cmd1 = "./samtools index "+path+sample+"/accepted_hits.bam"		## make samtools index accepted_hits.bam.bai
	os.system(cmd1)
	print("./Samtools index run completed.")

	for chrom in chromosomes:
		print("Start of ", chrom)
		tt = time.time()
		cmd2 = "./samtools view -b "+path+sample+"/accepted_hits.bam "+chrom+" -o "+out_directory+chrom+".bam"
		cmd3 = "./samtools pileup "+out_directory+chrom+".bam | cut -f 1,2,4 > "+out_directory+chrom+".txt"    ### Need to use pileup, not mpileup
		command = cmd2+";"+cmd3
		os.system(command)
		print("Samtools Time: ", time.time()-tt)
	print("Input text file generation completed")

