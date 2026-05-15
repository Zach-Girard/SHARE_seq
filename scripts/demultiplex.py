#!/usr/bin/env python
import argparse
import gzip
import logging
import itertools
import csv
import pandas as pd

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


def revcomp(seq):
	tab = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")
	return seq.translate(tab)[::-1]


def read_fasta(path, rc=False):
	"""Return ({sample_id: barcode_seq}, {}) from FASTA (no sample-type column)."""
	barcode = {}
	current_id = None
	current_seq = []
	with open(path) as f:
		for raw in f:
			line = raw.strip()
			if not line:
				continue
			if line.startswith(">"):
				if current_id is not None:
					seq = "".join(current_seq).upper()
					barcode[current_id] = revcomp(seq) if rc else seq
				current_id = line[1:].split()[0]
				current_seq = []
			else:
				current_seq.append(line)
	if current_id is not None:
		seq = "".join(current_seq).upper()
		barcode[current_id] = revcomp(seq) if rc else seq
	# FASTA has no table columns for sample type; caller may ignore second value.
	return barcode, {}


def read_barcode_table(path, rc=False):
	"""Return ({sample_id: barcode_seq}, {sample_id: sample_type}) from TSV/CSV table.
	Column 1 = sample name, column 2 = barcode (typical), column 3 = sample type when present."""
	rows = []
	with open(path, newline="") as f:
		sample = f.read(4096)
		f.seek(0)
		try:
			dialect = csv.Sniffer().sniff(sample)
			reader = csv.reader(f, dialect)
		except Exception:
			reader = csv.reader(f, delimiter="\t")
		for row in reader:
			row = [x.strip() for x in row if str(x).strip() != ""]
			if row:
				rows.append(row)
	if not rows:
		raise ValueError(f"No barcode rows found in {path}")

	def looks_like_barcode(s):
		u = s.upper()
		return len(u) >= 6 and all(c in "ACGTN" for c in u)

	barcode = {}
	# Column 1 = sample name (sid); column 3 = sample type (RNA/ATAC/etc.) when present.
	sample_type_by_sid = {}
	for row in rows:
		sample_type = ""
		if len(row) >= 3:
			sample_type = row[2].strip()
		if len(row) == 1:
			sid = row[0]
			seq = row[0]
		else:
			a, b = row[0], row[1]
			if looks_like_barcode(a) and not looks_like_barcode(b):
				seq, sid = a, b
			elif looks_like_barcode(b) and not looks_like_barcode(a):
				sid, seq = a, b
			else:
				# default: first column sample ID, second column barcode
				sid, seq = a, b
		seq = seq.upper()
		barcode[sid] = revcomp(seq) if rc else seq
		sample_type_by_sid[sid] = sample_type
	return barcode, sample_type_by_sid


def mismatch_sequence(seq, max_mismatch):
	"""Generate sequences within <= max_mismatch substitutions."""
	seq = seq.upper()
	alphabet = ("A", "C", "G", "T")
	out = {seq}
	if max_mismatch <= 0:
		return out
	positions = range(len(seq))
	for d in range(1, max_mismatch + 1):
		for idxs in itertools.combinations(positions, d):
			for repl in itertools.product(alphabet, repeat=d):
				s = list(seq)
				ok = False
				for i, base in zip(idxs, repl):
					if s[i] != base:
						s[i] = base
						ok = True
				if ok:
					out.add("".join(s))
	return out




def main():

	args = my_args()
	
	# read barcode (+ optional sample type from column 3 of table input)
	logging.info("Reading barcode file")
	sample_type_by_sid = {}
	try:
		barcode, sample_type_by_sid = read_barcode_table(args.barcode, rc=args.revcomp)
	except Exception:
		barcode, sample_type_by_sid = read_fasta(args.barcode, rc=args.revcomp)
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

	# Map observed index sequence -> sample name (column 1); include mismatch variants.
	seq_to_sid = {}
	for sid, seq in barcode.items():
		seq_to_sid[seq] = sid
		if args.num_mismatch > 0:
			for k in mismatch_sequence(seq, args.num_mismatch):
				seq_to_sid[k] = sid

	rows_out = []
	for seq_key, total in sorted(log_count.items(), key=lambda kv: kv[1], reverse=True):
		if seq_key == "unmatched":
			rows_out.append(
				{
					"Sample_Index": "unmatched",
					"Sample_Name": "unmatched",
					"Sample_Type": "",
					"Total_reads": total,
				}
			)
		else:
			sid = seq_to_sid.get(seq_key, "")
			rows_out.append(
				{
					"Sample_Index": seq_key,
					"Sample_Name": sid if sid else "",
					"Sample_Type": sample_type_by_sid.get(sid, "") if sid else "",
					"Total_reads": total,
				}
			)
	df = pd.DataFrame(rows_out)
	print(df)
	df.to_csv("SHARE-seq.demultiplex.stats.tsv", sep="\t", index=False)
if __name__ == "__main__":
	main()

























