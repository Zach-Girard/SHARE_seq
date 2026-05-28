#!/usr/bin/env python3
"""
Build cutadapt barcode FASTA for sgRNA demultiplexing.

Reads sgrna_demux_barcodes.tsv (Sample_Name, Sample_Index, ...) and writes:

  >sgRNA_C1200
  ^GCAGAGTC

The FASTA header (sample name) becomes cutadapt -o {name}.fastq.gz output basename.
The sequence uses ^ for 5' anchored matching (cutadapt -g).
"""

from __future__ import annotations

import argparse
import csv
import re
import sys


def load_samples(path: str) -> list[tuple[str, str]]:
    samples: list[tuple[str, str]] = []
    with open(path, newline="", errors="replace") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        if reader.fieldnames:
            for row in reader:
                name = (row.get("Sample_Name") or row.get("sample_name") or "").strip()
                index = (row.get("Sample_Index") or row.get("sample_index") or "").strip()
                if name and index:
                    samples.append((name, index))
            if samples:
                return samples
        fh.seek(0)
        reader2 = csv.reader(fh, delimiter="\t")
        for row in reader2:
            if not row or row[0].startswith("#"):
                continue
            if row[0].lower() in ("sample_name", "sample"):
                continue
            if len(row) < 2:
                continue
            samples.append((row[0].strip(), row[1].strip()))
    return samples


def write_barcode_fa(path: str, samples: list[tuple[str, str]]) -> None:
    seen_names: set[str] = set()
    with open(path, "w") as fh:
        for name, index in samples:
            if name in seen_names:
                raise ValueError(f"Duplicate sample name in barcode table: {name}")
            seen_names.add(name)
            idx = index.upper()
            if not re.fullmatch(r"[ACGTN]+", idx):
                raise ValueError(f"Invalid Sample_Index for {name}: {index!r}")
            fh.write(f">{name}\n")
            fh.write(f"^{idx}\n")


def main() -> int:
    p = argparse.ArgumentParser(description="Build cutadapt barcode.fa for sgRNA demux.")
    p.add_argument("--barcode-table", required=True, help="sgrna_demux_barcodes.tsv")
    p.add_argument("-o", "--output", default="barcode.fa")
    args = p.parse_args()

    if not samples := load_samples(args.barcode_table):
        print(f"ERROR: no samples in {args.barcode_table}", file=sys.stderr)
        return 1

    write_barcode_fa(args.output, samples)
    print(f"Wrote {args.output} ({len(samples)} barcodes)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
