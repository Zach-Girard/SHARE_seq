#!/usr/bin/env python3

import argparse
import gzip
import sys
from typing import List, Optional


def hamming_distance(a: str, b: str) -> int:
    if len(a) != len(b):
        return max(len(a), len(b))
    return sum(1 for x, y in zip(a, b) if x != y)


def load_whitelist(path: str):
    barcodes = []
    with open(path) as f:
        for line in f:
            bc = line.strip()
            if not bc:
                continue
            barcodes.append(bc)
    return barcodes


def open_fastq_read(path: str):
    """Open FASTQ for reading; gzip if path ends with .gz."""
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path, "rt")


def open_fastq_write(path: str):
    """Open FASTQ for writing; gzip if path ends with .gz."""
    return gzip.open(path, "wt") if path.endswith(".gz") else open(path, "w")


def match_eight_mer(mer: str, whitelist_8: List[str]) -> Optional[str]:
    """
    Match one 8 bp block to the whitelist: prefer exact, else any entry with Hamming distance 1.
    Returns the corrected whitelist sequence, or None if no match.
    """
    if len(mer) != 8:
        return None
    eightmers = [w for w in whitelist_8 if len(w) == 8]
    for w in eightmers:
        if w == mer:
            return w
    for w in eightmers:
        if hamming_distance(mer, w) == 1:
            return w
    return None


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Filter paired-end R1/R3 reads by matching the 24 bp barcode prefix on R3 "
            "(three 8 bp blocks) against an 8 bp whitelist file: each block must match "
            "some whitelist entry with at most 1 mismatch (same idea as STARsolo 1MM per block)."
        )
    )
    parser.add_argument("--r1", required=True, help="Input R1 FASTQ (.gz or plain).")
    parser.add_argument("--r3", required=True, help="Input R3 FASTQ (.gz or plain).")
    parser.add_argument(
        "--whitelist",
        required=True,
        help="Path to whitelist with one 8 bp barcode per line (e.g. barcodes_RC.txt).",
    )
    parser.add_argument(
        "--out-r1",
        required=True,
        help="Output filtered R1 FASTQ (use .gz for gzip).",
    )
    parser.add_argument(
        "--out-r3",
        required=True,
        help="Output filtered R3 FASTQ (use .gz for gzip).",
    )
    parser.add_argument(
        "--summary",
        required=True,
        help="Output summary text file with basic statistics.",
    )

    args = parser.parse_args()

    whitelist = load_whitelist(args.whitelist)
    if not whitelist:
        print(f"Whitelist file {args.whitelist} is empty or missing valid barcodes.", file=sys.stderr)
        sys.exit(1)

    total_reads = 0
    kept_reads = 0
    rejected_reads = 0

    with open_fastq_read(args.r1) as r1_in, open_fastq_read(args.r3) as r3_in, open_fastq_write(
        args.out_r1
    ) as r1_out, open_fastq_write(args.out_r3) as r3_out:
        while True:
            r1_lines = [r1_in.readline() for _ in range(4)]
            r3_lines = [r3_in.readline() for _ in range(4)]

            if not r1_lines[0] or not r3_lines[0]:
                break

            total_reads += 1

            seq_r3 = r3_lines[1].strip()
            if len(seq_r3) < 24:
                rejected_reads += 1
                continue

            bc24 = seq_r3[:24]
            b1, b2, b3 = bc24[0:8], bc24[8:16], bc24[16:24]

            w1 = match_eight_mer(b1, whitelist)
            w2 = match_eight_mer(b2, whitelist)
            w3 = match_eight_mer(b3, whitelist)

            if w1 is None or w2 is None or w3 is None:
                rejected_reads += 1
                continue

            corrected_bc24 = w1 + w2 + w3
            kept_reads += 1
            for line in r1_lines:
                r1_out.write(line)

            seq_full = r3_lines[1].strip()
            new_seq = corrected_bc24 + seq_full[24:]
            r3_out.write(r3_lines[0])
            r3_out.write(new_seq + "\n")
            r3_out.write(r3_lines[2])
            r3_out.write(r3_lines[3])

    with open(args.summary, "w") as sf:
        sf.write(f"Total reads processed:\t{total_reads}\n")
        sf.write(f"Reads kept (all 3x8bp matched within 1 mismatch):\t{kept_reads}\n")
        sf.write(f"Reads rejected:\t{rejected_reads}\n")
        if total_reads > 0:
            frac = kept_reads / total_reads
            sf.write(f"Fraction kept:\t{frac:.4f}\n")


if __name__ == "__main__":
    main()
