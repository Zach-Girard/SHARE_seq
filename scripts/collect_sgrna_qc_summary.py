#!/usr/bin/env python3
"""
Collect sgRNA demux and counting outputs into one QC summary TSV.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import os
import statistics
import sys
from typing import Dict, List, Optional


def parse_args():
    p = argparse.ArgumentParser(description="Collect sgRNA QC summary metrics.")
    p.add_argument("--manifest", required=True, help="sgRNA_run.tsv with fastq/sample/library columns")
    p.add_argument("--out-dir", required=True, help="sgRNA analysis output directory")
    p.add_argument("--out", default="sgrna_qc_summary.tsv", help="Output TSV")
    return p.parse_args()


def count_fastq_reads(path: str) -> Optional[int]:
    if not path or not os.path.isfile(path):
        return None
    try:
        opener = gzip.open if path.endswith(".gz") else open
        with opener(path, "rt", errors="replace") as fh:
            return sum(1 for _ in fh) // 4
    except Exception:
        return None


def load_matrix(path: str) -> tuple[int, int, int, float, float]:
    """Return cells, cells_with_gRNA, assigned_counts, mean_counts_per_cell, median_counts_per_cell."""
    if not os.path.isfile(path):
        return 0, 0, 0, 0.0, 0.0
    row_sums: List[int] = []
    with open(path, newline="", errors="replace") as fh:
        reader = csv.reader(fh)
        header = next(reader, None)
        if not header:
            return 0, 0, 0, 0.0, 0.0
        for row in reader:
            vals = []
            for val in row[1:]:
                try:
                    vals.append(int(float(val)))
                except Exception:
                    vals.append(0)
            row_sums.append(sum(vals))
    cells = len(row_sums)
    cells_with_grna = sum(1 for v in row_sums if v > 0)
    assigned = sum(row_sums)
    mean_counts = statistics.mean(row_sums) if row_sums else 0.0
    median_counts = statistics.median(row_sums) if row_sums else 0.0
    return cells, cells_with_grna, assigned, mean_counts, median_counts


def load_grna_counts(path: str) -> tuple[int, int, str, int]:
    """Return library_size, detected_gRNAs, top_gRNA, top_count."""
    if not os.path.isfile(path):
        return 0, 0, "", 0
    library_size = 0
    detected = 0
    top_grna = ""
    top_count = 0
    with open(path, newline="", errors="replace") as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            library_size += 1
            grna = row.get("sgRNA_sequence", "")
            try:
                count = int(float(row.get("count", "0") or 0))
            except Exception:
                count = 0
            if count > 0:
                detected += 1
            if count > top_count:
                top_grna = grna
                top_count = count
    return library_size, detected, top_grna, top_count


def fmt(value) -> str:
    if value is None:
        return ""
    if isinstance(value, float):
        return f"{value:.2f}"
    return str(value)


def main() -> int:
    args = parse_args()
    if not os.path.isfile(args.manifest):
        print(f"ERROR: manifest not found: {args.manifest}", file=sys.stderr)
        return 1

    out_dir = os.path.abspath(args.out_dir)
    rows: List[Dict[str, object]] = []
    with open(args.manifest, newline="", errors="replace") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            sample = (row.get("sample_name") or "").strip()
            if not sample:
                continue
            demux_fastq = (row.get("fastq") or "").strip()
            sample_dir = os.path.join(out_dir, sample)
            matched_fastq = os.path.join(sample_dir, f"{sample}.matched.R1.fastq.gz")
            matrix = os.path.join(sample_dir, f"final_{sample}.gRNA.count.csv")
            counts = os.path.join(sample_dir, "gRNA_counts_final.csv")

            cells, cells_with_grna, assigned, mean_counts, median_counts = load_matrix(matrix)
            library_size, detected, top_grna, top_count = load_grna_counts(counts)
            demux_reads = count_fastq_reads(demux_fastq)
            matched_reads = count_fastq_reads(matched_fastq)
            pct_cells_with_grna = (
                100.0 * cells_with_grna / cells
            ) if cells else None
            pct_guides_detected = (
                100.0 * detected / library_size
            ) if library_size else None

            rows.append(
                {
                    "sample_name": sample,
                    "experimental_group": row.get("experimental_group", ""),
                    "demux_reads": demux_reads,
                    "matched_reads": matched_reads,
                    "assigned_gRNA_reads": assigned,
                    "cells": cells,
                    "cells_with_gRNA": cells_with_grna,
                    "pct_cells_with_gRNA": pct_cells_with_grna,
                    "library_size": library_size,
                    "detected_gRNAs": detected,
                    "pct_gRNAs_detected": pct_guides_detected,
                    "mean_gRNA_counts_per_cell": mean_counts,
                    "median_gRNA_counts_per_cell": median_counts,
                    "top_gRNA": top_grna,
                    "top_gRNA_count": top_count,
                }
            )

    out_path = os.path.join(out_dir, args.out)
    fieldnames = [
        "sample_name",
        "experimental_group",
        "demux_reads",
        "matched_reads",
        "assigned_gRNA_reads",
        "cells",
        "cells_with_gRNA",
        "pct_cells_with_gRNA",
        "library_size",
        "detected_gRNAs",
        "pct_gRNAs_detected",
        "mean_gRNA_counts_per_cell",
        "median_gRNA_counts_per_cell",
        "top_gRNA",
        "top_gRNA_count",
    ]
    with open(out_path, "w", newline="") as out:
        writer = csv.DictWriter(out, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow({key: fmt(row.get(key)) for key in fieldnames})

    print(f"Wrote {out_path} ({len(rows)} samples)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
