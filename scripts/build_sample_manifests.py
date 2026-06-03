#!/usr/bin/env python3
"""
Parse the combined sample_barcode_file (input.tsv) and:

1. Write sgRNA.tsv (sgRNA sample metadata).
2. Write sgrna_demux_barcodes.tsv and sgrna_barcode.fa (cutadapt demux from shared Undetermined R1).
3. Write demux_barcodes.tsv for global Undetermined demux (RNA/ATAC only).

sgRNA row layout (5 columns; demux from shared sgRNA Undetermined FASTQs, not per-row paths):

  sample_name  Sample_Index  sgRNA  Experimental_Group  grna_library_csv

Example:

  sgRNA_C1200  gcagagtc  sgRNA  C1200  C1200_gRNA_library.csv
  sgRNA_C6991  gagcagca  sgRNA  C6991  C6991_gRNA_library.csv

Undetermined input is paired R1+R2 FASTQs via params (defaults: sgRNA_Undetermined_S0_R1_001.fastq.gz and R2 in RAW_FASTQ/).

RNA/ATAC rows (columns 1–4 minimum):

  sample_name  Sample_Index  RNA|ATAC  Experimental_Group
"""

from __future__ import annotations

import argparse
import csv
import os
import re
import sys
from typing import List, Optional


def split_cols(line: str) -> List[str]:
    return [c.strip() for c in re.split(r"\t|,", line)]


def is_header_row(cols: List[str]) -> bool:
    if not cols:
        return True
    c0 = cols[0].lower()
    c2 = cols[2].lower() if len(cols) > 2 else ""
    if c0 in ("sample", "sample_name", "sample_id"):
        return True
    if c2 in ("type", "sample_type"):
        return True
    return False


def normalize_type(value: str) -> str:
    u = (value or "").strip().upper()
    if u == "SGRNA":
        return "sgRNA"
    if u in ("RNA", "ATAC"):
        return u
    return ""


def looks_like_library_csv(value: str) -> bool:
    return (value or "").lower().endswith(".csv")


def resolve_path(project_dir: str, raw_fastq_dir: str, path_value: str) -> str:
    if not path_value:
        return ""
    if os.path.isabs(path_value):
        return os.path.abspath(path_value)
    for base in (project_dir, raw_fastq_dir, os.path.join(project_dir, raw_fastq_dir)):
        if not base:
            continue
        candidate = os.path.abspath(os.path.join(base, path_value))
        if os.path.exists(candidate):
            return candidate
    return os.path.abspath(os.path.join(project_dir, path_value))


def parse_row(cols: List[str]) -> Optional[dict]:
    if len(cols) < 3:
        return None
    sample = cols[0]
    if not sample:
        return None

    stype = normalize_type(cols[2])
    if not stype:
        return None

    if stype == "sgRNA":
        if len(cols) < 5:
            return None
        lib = cols[-1]
        if not looks_like_library_csv(lib):
            return None
        return {
            "sample": sample,
            "sample_index": cols[1],
            "sample_type": stype,
            "experimental_group": cols[3] if len(cols) > 3 else "",
            "grna_library_csv": lib,
        }

    return {
        "sample": sample,
        "sample_index": cols[1] if len(cols) > 1 else "",
        "sample_type": stype,
        "experimental_group": cols[3] if len(cols) > 3 else "",
        "grna_library_csv": "",
    }


def load_rows(path: str) -> List[dict]:
    rows: List[dict] = []
    with open(path, "r", errors="replace") as fh:
        for raw in fh:
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            cols = split_cols(line)
            if is_header_row(cols):
                continue
            parsed = parse_row(cols)
            if parsed:
                rows.append(parsed)
    return rows


def write_demux_barcodes(path: str, demux_rows: List[dict]) -> None:
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["Sample_Name", "Sample_Index", "Sample_Type", "Experimental_Group"])
        for row in demux_rows:
            w.writerow(
                [
                    row["sample"],
                    row["sample_index"],
                    row["sample_type"],
                    row.get("experimental_group", ""),
                ]
            )


def write_sgrna_barcode_fa(path: str, sgrna_rows: List[dict]) -> None:
    """cutadapt -g file:barcode.fa format: >Sample_Name then ^Sample_Index."""
    seen: set[str] = set()
    with open(path, "w") as fh:
        for row in sgrna_rows:
            name = row["sample"]
            if name in seen:
                raise ValueError(f"Duplicate sgRNA sample name: {name}")
            seen.add(name)
            idx = row["sample_index"].strip().upper()
            if not idx or not re.fullmatch(r"[ACGTN]+", idx):
                raise ValueError(f"Invalid Sample_Index for {name}: {row['sample_index']!r}")
            fh.write(f">{name}\n")
            fh.write(f"^{idx}\n")


def write_sgrna_manifest(
    path: str, sgrna_rows: List[dict], project_dir: str, raw_fastq_dir: str
) -> None:
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["sample_name", "sample_index", "experimental_group", "grna_library_csv"])
        for row in sgrna_rows:
            lib_abs = resolve_path(project_dir, raw_fastq_dir, row["grna_library_csv"])
            w.writerow(
                [
                    row["sample"],
                    row["sample_index"],
                    row.get("experimental_group", ""),
                    lib_abs,
                ]
            )


def main() -> int:
    p = argparse.ArgumentParser(description="Build sgRNA and demux manifests from input.tsv.")
    p.add_argument("--sample-barcode-file", required=True)
    p.add_argument("--project-dir", required=True)
    p.add_argument("--raw-fastq-dir", default="RAW_FASTQ")
    p.add_argument("--out-sgrna", default="sgRNA.tsv")
    p.add_argument("--out-sgrna-demux-barcodes", default="sgrna_demux_barcodes.tsv")
    p.add_argument("--out-sgrna-barcode-fa", default="sgrna_barcode.fa")
    p.add_argument("--out-demux-barcodes", default="demux_barcodes.tsv")
    args = p.parse_args()

    project_dir = os.path.abspath(args.project_dir)
    raw_fastq_dir = os.path.abspath(
        args.raw_fastq_dir
        if os.path.isabs(args.raw_fastq_dir)
        else os.path.join(project_dir, args.raw_fastq_dir)
    )
    barcode_path = args.sample_barcode_file
    if not os.path.isfile(barcode_path):
        if not os.path.isabs(barcode_path):
            for candidate in (
                os.path.join(project_dir, barcode_path),
                os.path.join(project_dir, raw_fastq_dir, os.path.basename(barcode_path)),
            ):
                if os.path.isfile(candidate):
                    barcode_path = candidate
                    break

    if not os.path.isfile(barcode_path):
        print(f"ERROR: sample barcode file not found: {barcode_path}", file=sys.stderr)
        return 1

    all_rows = load_rows(barcode_path)
    demux_rows = [r for r in all_rows if r["sample_type"] in ("RNA", "ATAC")]
    sgrna_rows = [r for r in all_rows if r["sample_type"] == "sgRNA"]

    for row in demux_rows:
        if not row.get("sample_index"):
            print(
                f"ERROR: RNA/ATAC sample {row['sample']} missing Sample_Index (column 2).",
                file=sys.stderr,
            )
            return 1

    for row in sgrna_rows:
        if not row.get("sample_index"):
            print(
                f"ERROR: sgRNA sample {row['sample']} missing Sample_Index (column 2).",
                file=sys.stderr,
            )
            return 1
        if not row.get("grna_library_csv"):
            print(
                f"ERROR: sgRNA sample {row['sample']} requires grna_library_csv (column 5).",
                file=sys.stderr,
            )
            return 1
        lib = resolve_path(project_dir, raw_fastq_dir, row["grna_library_csv"])
        if not os.path.isfile(lib):
            print(f"ERROR: gRNA library CSV not found for {row['sample']}: {lib}", file=sys.stderr)
            return 1

    write_demux_barcodes(args.out_demux_barcodes, demux_rows)
    write_sgrna_manifest(args.out_sgrna, sgrna_rows, project_dir, raw_fastq_dir)
    write_demux_barcodes(args.out_sgrna_demux_barcodes, sgrna_rows)
    if sgrna_rows:
        write_sgrna_barcode_fa(args.out_sgrna_barcode_fa, sgrna_rows)

    print(f"Wrote {args.out_demux_barcodes} ({len(demux_rows)} RNA/ATAC samples)")
    print(f"Wrote {args.out_sgrna_demux_barcodes} ({len(sgrna_rows)} sgRNA samples)")
    if sgrna_rows:
        print(f"Wrote {args.out_sgrna_barcode_fa} (cutadapt barcode.fa)")
    print(f"Wrote {args.out_sgrna} ({len(sgrna_rows)} sgRNA samples)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
