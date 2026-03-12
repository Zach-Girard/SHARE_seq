#!/usr/bin/env python3

import argparse
import gzip
import sys
from collections import Counter


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


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Filter paired-end R1/R3 reads by matching 24bp barcodes in R3 "
            "to a whitelist with up to 1 mismatch in each 8bp block."
        )
    )
    parser.add_argument("--r1", required=True, help="Input R1 FASTQ (gzipped).")
    parser.add_argument("--r3", required=True, help="Input R3 FASTQ (gzipped).")
    parser.add_argument(
        "--whitelist",
        required=True,
        help="Path to barcodes_RC.txt whitelist (one barcode per line).",
    )
    parser.add_argument(
        "--out-r1",
        required=True,
        help="Output filtered R1 FASTQ (gzipped).",
    )
    parser.add_argument(
        "--out-r3",
        required=True,
        help="Output filtered R3 FASTQ (gzipped).",
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

    # Open input and output FASTQs
    with gzip.open(args.r1, "rt") as r1_in, gzip.open(
        args.r3, "rt"
    ) as r3_in, gzip.open(args.out_r1, "wt") as r1_out, gzip.open(
        args.out_r3, "wt"
    ) as r3_out:
        while True:
            # Read one record (4 lines) from each FASTQ
            r1_lines = [r1_in.readline() for _ in range(4)]
            r3_lines = [r3_in.readline() for _ in range(4)]

            # Check for EOF
            if not r1_lines[0] or not r3_lines[0]:
                break

            total_reads += 1

            seq_r3 = r3_lines[1].strip()
            if len(seq_r3) < 24:
                rejected_reads += 1
                continue

            bc24 = seq_r3[:24]
            b1, b2, b3 = bc24[0:8], bc24[8:16], bc24[16:24]

            matched = False
            corrected_bc24 = None
            for wl in whitelist:
                if len(wl) < 24:
                    continue
                w1, w2, w3 = wl[0:8], wl[8:16], wl[16:24]
                if (
                    hamming_distance(b1, w1) <= 1
                    and hamming_distance(b2, w2) <= 1
                    and hamming_distance(b3, w3) <= 1
                ):
                    matched = True
                    corrected_bc24 = wl[:24]
                    break

            if matched:
                kept_reads += 1
                for line in r1_lines:
                    r1_out.write(line)
                # Overwrite the first 24bp of the R3 sequence with the matched whitelist barcode
                seq_full = r3_lines[1].strip()
                if corrected_bc24 is not None and len(seq_full) >= 24:
                    new_seq = corrected_bc24 + seq_full[24:]
                    r3_out.write(r3_lines[0])
                    r3_out.write(new_seq + "\n")
                    r3_out.write(r3_lines[2])
                    r3_out.write(r3_lines[3])
                else:
                    # Fallback: write original if something unexpected happens
                    for line in r3_lines:
                        r3_out.write(line)
            else:
                rejected_reads += 1

    with open(args.summary, "w") as sf:
        sf.write(f"Total reads processed:\t{total_reads}\n")
        sf.write(f"Reads kept (all 3x8bp matched within 1 mismatch):\t{kept_reads}\n")
        sf.write(f"Reads rejected:\t{rejected_reads}\n")
        if total_reads > 0:
            frac = kept_reads / total_reads
            sf.write(f"Fraction kept:\t{frac:.4f}\n")


if __name__ == "__main__":
    main()

