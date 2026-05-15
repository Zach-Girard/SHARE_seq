#!/usr/bin/env python

import argparse
import gzip
import logging
import itertools
import csv
from collections import defaultdict


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


def revcomp(seq):
	tab = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")
	return seq.translate(tab)[::-1]


def hamming_distance(a, b):
	if len(a) != len(b):
		return 10**9
	return sum(ch1 != ch2 for ch1, ch2 in zip(a, b))


def _mismatch_neighbors(seq, max_mismatch):
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


def k_mer_distance(barcode_length, error, barcode_list):
	lookup = {}
	for bc in barcode_list:
		bc = bc.upper()
		if len(bc) != barcode_length:
			continue
		for variant in _mismatch_neighbors(bc, error):
			lookup.setdefault(variant, set()).add(bc)
	return lookup


def find_dist(kmer_lookup, query, barcode_list):
	query = query.upper()
	candidates = list(kmer_lookup.get(query, []))
	if not candidates:
		return (False, None)
	best = sorted(candidates, key=lambda bc: (hamming_distance(query, bc), bc))[0]
	return (True, best)


def read_fastq_record(handle):
	h = handle.readline()
	if not h:
		return None
	s = handle.readline()
	p = handle.readline()
	q = handle.readline()
	if not s or not p or not q:
		return None
	return (h, s, p, q)


def read_barcode_list(path):
	barcodes = []
	with open(path, newline="") as f:
		sample = f.read(4096)
		f.seek(0)
		try:
			dialect = csv.Sniffer().sniff(sample)
			reader = csv.reader(f, dialect)
		except Exception:
			reader = csv.reader(f, delimiter="\t")
		for row in reader:
			if not row:
				continue
			val = str(row[0]).strip().replace(" ", "")
			if val:
				barcodes.append(val)
	return barcodes

def share_seq_rename_fastq(Read1,Read2,label, barcode_list, error, sp1_length,bc1_length,sp2_length,bc2_length,sp3_length,bc3_length):
	"""general utils for demultiplexing read in PE mode"""

	barcode_counts = defaultdict(int)

	junk_R1_output = label + ".junk.R1.fastq.gz"
	junk_R2_output = label + ".junk.R2.fastq.gz"
	matched_R1_output = label + ".matched.R1.fastq.gz"
	matched_R2_output = label + ".matched.R2.fastq.gz"

	f1 = gzip.open(Read1, 'rt')
	f2 = gzip.open(Read2, 'rt')
	j1 = gzip.open(junk_R1_output, "wt")
	j2 = gzip.open(junk_R2_output, "wt")
	m1 = gzip.open(matched_R1_output, "wt")
	m2 = gzip.open(matched_R2_output, "wt")
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
	count = 0
	total_valid_reads = 0
	while True:
		rec1 = read_fastq_record(f1)
		rec2 = read_fastq_record(f2)
		if rec1 is None or rec2 is None:
			break
		count +=1
		if count % 10000 == 0:
			logging.info("%s reads processed"%(count))
		h1, s1, p1, q1 = rec1
		h2, s2, p2, q2 = rec2
		header_bits = h1.strip().split()
		name = header_bits[0]
		barcode_r1 = header_bits[-1].split(":")[-1] if len(header_bits) > 1 else ""
		name_r1 = name + "1\n"
		name_r2 = name + "2\n"

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
			barcode_counts[(barcode_1, barcode_2, barcode_3)] += 1
			first_line_r1 = name + "_"+"".join([barcode_1,barcode_2,barcode_3]) + "\n"
			first_line_r2 = first_line_r1
			m1.write(first_line_r1)
			m2.write(first_line_r2)
			m1.write(s1)
			m2.write(s2)
			m1.write(p1)
			m2.write(p2)
			m1.write(q1)
			m2.write(q2)
			total_valid_reads += 1

		else:
			j1.write(name_r1)
			j2.write(name_r2)
			j1.write(s1)
			j2.write(s2)
			j1.write(p1)
			j2.write(p2)
			j1.write(q1)
			j2.write(q2)

	f1.close()
	f2.close()
	j1.close()
	j2.close()
	m1.close()
	m2.close()

	with open(label + ".total_number_reads.tsv", "w", newline="") as out:
		w = csv.writer(out, delimiter="\t")
		w.writerow(["BC1", "BC2", "BC3", "Total_reads"])
		for (b1, b2, b3), c in sorted(barcode_counts.items(), key=lambda kv: kv[1], reverse=True):
			w.writerow([b1, b2, b3, c])

	print ("Sample: %s has %s BC1 %s BC2 %s BC3"%(label,bc1_count,bc2_count,bc3_count))
	print ("Sample: %s has %s BC1_2 %s BC1_3 %s BC2_3"%(label,bc1_bc2_count,bc1_bc3_count,bc2_bc3_count))
	print ("Sample: %s has %s total valid reads. Total read is %s"%(label,total_valid_reads,count))

def main():

	args = my_args()
	
	# read barcode
	logging.info("Reading barcode file")
	barcode_list = read_barcode_list(args.barcode_list)

	if args.revcomp:
		barcode_list = [revcomp(x) for x in barcode_list]


	share_seq_rename_fastq(args.read1,args.read2,args.sample_ID, barcode_list, args.error, args.spacer1_length,args.barcode_length,args.spacer2_length,args.barcode_length,args.spacer3_length,args.barcode_length)
	
if __name__ == "__main__":
	main()

























