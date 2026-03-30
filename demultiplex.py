#!/usr/bin/env python
import sys
import os
import argparse
import gzip
from Bio.Seq import Seq
from Bio import SeqIO
import logging
import Colorer
import itertools
import string
import pandas as pd
from utils import *

"""
delmultiplex undetermined R1 and R2 fastq by sample barcode.

Sample barcode can be provided as fasta or tsv


"""

current_file_base_name = __file__.split("/")[-1].split(".")[0]
def my_args():
	mainParser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	# mainParser.add_argument('-o',"--output",  help="output fastq file name", default="output.fastq.gz")	
	mainParser.add_argument('-r1',  help="input undetermined R1 fastq.gz", required=True)	
	mainParser.add_argument('-r2',  help="input undetermined R2 fastq.gz", required=True)	
	mainParser.add_argument('-b',"--barcode",  help="barcode file in fasta or table, csv or tsv format", required=True)	
	mainParser.add_argument('-n',"--num_mismatch",  help="number of mismatch allowed", default=1,type=int)	
	mainParser.add_argument("--revcomp",  help="revcomp barcode",action='store_true')

	##------- add parameters above ---------------------
	args = mainParser.parse_args()	
	return args




def main():

	args = my_args()
	
	# read barcode
	logging.info("Reading barcode file")
	try:
		barcode = read_barcode_table(args.barcode,rc = args.revcomp)
	except:
		barcode = read_fasta(args.barcode,rc = args.revcomp)
	# open files
	R1 = {}
	R2 = {}
	log_count = {}
	for b in barcode:
		R1[barcode[b]] = gzip.open("%s.R1.fastq.gz"%(b), 'wt')
		R2[barcode[b]] = gzip.open("%s.R2.fastq.gz"%(b), 'wt')
		log_count[barcode[b]] = 0
	log_count["unmatched"] = 0												
	# add mismatch
	if args.num_mismatch>0:
		for b in barcode:
			# print (b)
			for k in mismatch_sequence(barcode[b],args.num_mismatch):
				R1[k] = R1[barcode[b]]
				R2[k] = R2[barcode[b]]
	
	R1['unmatched'] = gzip.open("unmatched.R1.fastq.gz", 'wt')
	R2['unmatched'] = gzip.open("unmatched.R2.fastq.gz", 'wt')
	# read fastq
	# logging.info("Reading barcode file")
	n = 4
	count = 0
	# https://www.biostars.org/p/317524/
	with gzip.open(args.r1, 'rt') as fh:
		lines = []
		for line in fh:
			lines.append(line.strip())
			if len(lines) == n:
				# print (lines)
				
				# logging.info("%s reads processed"%(count))
				if count %100000 == 0:
					# print (count,"reads has been processed")
					logging.info("%s matched R1 reads processed"%(count))
				# @NB551526:91:HC27NAFX2:1:11101:17620:1055 1:N:0:CTGCCTAA
				# print (lines[0])
				myString = lines[0].split("+")[-1]
				# print (myString)
				try:
					count += 1
					R1[myString].write("\n".join(lines)+"\n")
					log_count[myString] += 1
				except:
					R1['unmatched'].write("\n".join(lines)+"\n")
					log_count['unmatched'] += 1
				lines=[]

	n = 4
	count = 0
	# https://www.biostars.org/p/317524/
	with gzip.open(args.r2, 'rt') as fh:
		lines = []
		for line in fh:
			lines.append(line.strip())
			if len(lines) == n:
				# print (lines)
				
				# logging.info("%s reads processed"%(count))
				if count %100000 == 0:
					# print (count,"reads has been processed")
					logging.info("%s matched R2 reads processed"%(count))
				# @NB551526:91:HC27NAFX2:1:11101:17620:1055 1:N:0:CTGCCTAA
				myString = lines[0].split("+")[-1]
				# print (myString)
				try:
					count += 1
					R2[myString].write("\n".join(lines)+"\n")
				except:
					R2['unmatched'].write("\n".join(lines)+"\n")
				lines=[]

	df = pd.DataFrame.from_dict(log_count,orient="index",columns=['Total_reads'])
	df = df.sort_values("Total_reads",ascending=False)
	print (df)
	df.to_csv("SHARE-seq.demultiplex.stats.tsv",sep="\t")
if __name__ == "__main__":
	main()

























