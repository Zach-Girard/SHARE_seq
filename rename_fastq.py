#!/usr/bin/env python

import sys
import os
import argparse
import gzip
import logging
import Colorer
import itertools
import subprocess
import sys
import gzip
import os
import argparse
import pandas as pd
import subprocess
import glob
from utils import *


current_file_base_name = __file__.split("/")[-1].split(".")[0]
def my_args():
	mainParser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

	mainParser.add_argument("-r1","--read1",  help="R1 fastq file", required=True)
	mainParser.add_argument("-r2","--read2",  help="R2 fastq file", required=True)
	mainParser.add_argument("--sample_ID",  help="sample_ID", required=True)
	mainParser.add_argument("--barcode_list",  help="list of barcodes 1", required=True)
	mainParser.add_argument("--error",  help="barcode 1 allowed number of mismatches", type=int,default=1)
	mainParser.add_argument("--spacer1_length",  help="spacer1_length", type=int, default=15)
	mainParser.add_argument("--barcode_length",  help="barcode_length", type=int, default=8)
	mainParser.add_argument("--spacer2_length",  help="spacer2_length", type=int, default=30)
	mainParser.add_argument("--spacer3_length",  help="spacer3_length", type=int, default=30)
	mainParser.add_argument("--revcomp",  help="reverse complement barcode sequence",action='store_true')


	##------- add parameters above ---------------------
	args = mainParser.parse_args()	
	return args


def output_to_fastq_gz(file_name,list_of_lines):
	f = gzip.open(file_name,"wt")
	[f.write(x) for x in list_of_lines]
	f.close()

def share_seq_rename_fastq(Read1,Read2,label, barcode_list, error, sp1_length,bc1_length,sp2_length,bc2_length,sp3_length,bc3_length):
	"""general utils for demultiplexing read in PE mode"""
	
	barcode_dict = get_barcode_dict(barcode_list) # used to count reads
	
	junk_list_R1 = []
	junk_list_R2 = []
	matched_list_R1 = []
	matched_list_R2 = []
	f1 = gzip.open(Read1, 'rt')
	f2 = gzip.open(Read2, 'rt')	
	logging.info("Allowed barcode error: %s"%(error))
	# edit distance is better than hamming distance
	bc3_kmer = k_mer_distance(bc3_length,error,barcode_list=barcode_list)
	bc2_kmer = bc3_kmer
	bc1_kmer = bc3_kmer
	logging.info("Finish enumerate possible kmers")
	bc1_start = sp1_length
	bc1_end = bc1_start+bc1_length
	bc2_start = bc1_end+sp2_length
	bc2_end = bc2_start+bc2_length
	bc3_start = bc2_end+sp3_length
	bc3_end = bc3_start+bc3_length
	
	bc3_count =0
	bc2_count =0
	bc1_count =0
	bc1_bc2_count =0
	bc1_bc3_count =0
	bc2_bc3_count =0
	line1 = f1.readline()
	line2 = f2.readline()
	count = 0
	total_valid_reads = 0
	while (line1):
		count +=1
		if count % 10000 == 0:
			logging.info("%s reads processed"%(count))
		name,barcode_r1 = line1.split()
		# print (line2)
		_ = line2.split()
		name_r1 = name + "1\n" 
		name_r2 = name + "2\n"
		barcode_r1 = barcode_r1.split(":")[-1]

		line1 = f1.readline()
		line2 = f2.readline()
		bc3 = barcode_r1[bc3_start:bc3_end]
		bc2 = barcode_r1[bc2_start:bc2_end]
		bc1 = barcode_r1[bc1_start:bc1_end]
		# print (barcode_r1,bc1,bc2,bc3)
		bc1_dist,barcode_1 = find_dist(bc1_kmer,bc1,barcode_list)
		bc2_dist,barcode_2 = find_dist(bc2_kmer,bc2,barcode_list)
		bc3_dist,barcode_3 = find_dist(bc3_kmer,bc3,barcode_list)
		
		if bc1_dist:
			bc1_count+=1
		if bc2_dist:
			bc2_count+=1
		if bc3_dist:
			bc3_count+=1
		if bc1_dist and bc2_dist:
			bc1_bc2_count += 1
		if bc1_dist and bc3_dist:
			bc1_bc3_count += 1
		if bc2_dist and bc3_dist:
			bc2_bc3_count += 1
		if bc1_dist and bc2_dist and bc3_dist :
			barcode_dict[barcode_1][barcode_2][barcode_3]+=1
			# first_line_r1 = '@' + ",".join(["NNNN",barcode_1,barcode_2,barcode_3]) + ',' + name_r1[1:]
			# first_line_r2 = '@' + ",".join(["NNNN",barcode_1,barcode_2,barcode_3])+ ',' + name_r2[1:]
			first_line_r1 = name + "_"+"".join([barcode_1,barcode_2,barcode_3]) + "\n"
			first_line_r2 = first_line_r1
			matched_list_R1.append(first_line_r1)
			matched_list_R2.append(first_line_r2)

			matched_list_R1.append(line1)
			matched_list_R2.append(line2)	
			
			third_line_r1 = f1.readline()
			third_line_r2 = f2.readline()
			matched_list_R1.append(third_line_r1)
			matched_list_R2.append(third_line_r2)	
	
			four_line_r1 = f1.readline()
			four_line_r2 = f2.readline()
			matched_list_R1.append(four_line_r1)
			matched_list_R2.append(four_line_r2)	
			total_valid_reads += 1

		else:
			junk_list_R1.append(name_r1)
			junk_list_R2.append(name_r2)
			junk_list_R1.append(line1)
			junk_list_R2.append(line2)
			
			third_line_r1 = f1.readline()
			third_line_r2 = f2.readline()
			junk_list_R1.append(third_line_r1)
			junk_list_R2.append(third_line_r2)	
	
			four_line_r1 = f1.readline()
			four_line_r2 = f2.readline()
			junk_list_R1.append(four_line_r1)
			junk_list_R2.append(four_line_r2)	
				
			
		line1 = f1.readline()
		line2 = f2.readline()

	junk_R1_output = label + ".junk.R1.fastq.gz"
	junk_R2_output = label + ".junk.R2.fastq.gz"
	output_to_fastq_gz(junk_R1_output,junk_list_R1)
	output_to_fastq_gz(junk_R2_output,junk_list_R2)
	
	matched_R1_output = label + ".matched.R1.fastq.gz"
	matched_R2_output = label + ".matched.R2.fastq.gz"
	output_to_fastq_gz(matched_R1_output,matched_list_R1)
	output_to_fastq_gz(matched_R2_output,matched_list_R2)
	df = dict3d_to_df(barcode_dict)
	df.to_csv(label + ".total_number_reads.tsv",sep="\t",index=False)
	print ("Sample: %s has %s BC1 %s BC2 %s BC3"%(label,bc1_count,bc2_count,bc3_count))
	print ("Sample: %s has %s BC1_2 %s BC1_3 %s BC2_3"%(label,bc1_bc2_count,bc1_bc3_count,bc2_bc3_count))
	print ("Sample: %s has %s total valid reads. Total read is %s"%(label,total_valid_reads,count))

def main():

	args = my_args()
	
	# read barcode
	logging.info("Reading barcode file")
	barcode_list = pd.read_csv(args.barcode_list,header=None)[0].tolist()
	barcode_list = [x.replace(" ","") for x in barcode_list]

	if args.revcomp:
		barcode_list = [revcomp(x) for x in barcode_list]


	share_seq_rename_fastq(args.read1,args.read2,args.sample_ID, barcode_list, args.error, args.spacer1_length,args.barcode_length,args.spacer2_length,args.barcode_length,args.spacer3_length,args.barcode_length)
	
if __name__ == "__main__":
	main()

























