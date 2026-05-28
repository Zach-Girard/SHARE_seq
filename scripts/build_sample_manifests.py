#!/usr/bin/env python3
"""
Parse the combined sample_barcode_file (input.tsv) and:

1. Write sgRNA.tsv for sgRNA analysis (fastq, sample_name, grna_library_csv).
2. Write demux_barcodes.tsv containing only RNA/ATAC rows for undetermined demultiplexing.

sgRNA row layouts supported
---------------------------
Five columns (no demux index; library in column 2):

  sample_name  grna_library_csv  sgRNA  Experimental_Group  fastq

Six columns (optional demux index placeholder in column 2, e.g. NA):

  sample_name  Sample_Index  sgRNA  Experimental_Group  grna_library_csv  fastq

RNA/ATAC rows use the standard layout (columns 1–4 minimum):

  sample_name  Sample_Index  RNA|ATAC  Experimental_Group
"""

from __future__ import annotations

import argparse
import csv
import os
import re
import sys
from typing import List, Optional, Tuple


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
        if len(cols) >= 6:
            return {
                "sample": sample,
                "sample_index": cols[1],
                "sample_type": stype,
                "experimental_group": cols[3] if len(cols) > 3 else "",
                "grna_library_csv": cols[4],
                "fastq": cols[5],
            }
        if len(cols) >= 5:
            return {
                "sample": sample,
                "sample_index": "",
                "sample_type": stype,
                "experimental_group": cols[3],
                "grna_library_csv": cols[1],
                "fastq": cols[4],
            }
        return None

    # RNA / ATAC
    return {
        "sample": sample,
        "sample_index": cols[1] if len(cols) > 1 else "",
        "sample_type": stype,
        "experimental_group": cols[3] if len(cols) > 3 else "",
        "grna_library_csv": "",
        "fastq": "",
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


def write_sgrna_manifest(
    path: str, sgrna_rows: List[dict], project_dir: str, raw_fastq_dir: str
) -> None:
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["fastq", "sample_name", "grna_library_csv", "experimental_group"])
        for row in sgrna_rows:
            fastq_abs = resolve_path(project_dir, raw_fastq_dir, row["fastq"])
            lib_abs = resolve_path(project_dir, raw_fastq_dir, row["grna_library_csv"])
            w.writerow(
                [
                    fastq_abs,
                    row["sample"],
                    lib_abs,
                    row.get("experimental_group", ""),
                ]
            )


def main() -> int:
    p = argparse.ArgumentParser(description="Build sgRNA.tsv and demux_barcodes.tsv from input.tsv.")
    p.add_argument("--sample-barcode-file", required=True)
    p.add_argument("--project-dir", required=True)
    p.add_argument("--raw-fastq-dir", default="RAW_FASTQ")
    p.add_argument("--out-sgrna", default="sgRNA.tsv")
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
        if os.path.isabs(barcode_path):
            pass
        else:
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
        if not row.get("grna_library_csv") or not row.get("fastq"):
            print(
                f"ERROR: sgRNA sample {row['sample']} requires grna_library_csv and fastq.",
                file=sys.stderr,
            )
            return 1
        lib = resolve_path(project_dir, raw_fastq_dir, row["grna_library_csv"])
        fq = resolve_path(project_dir, raw_fastq_dir, row["fastq"])
        if not os.path.isfile(lib):
            print(f"ERROR: gRNA library CSV not found for {row['sample']}: {lib}", file=sys.stderr)
            return 1
        if not os.path.isfile(fq):
            print(f"ERROR: FASTQ not found for {row['sample']}: {fq}", file=sys.stderr)
            return 1

    write_demux_barcodes(args.out_demux_barcodes, demux_rows)
    write_sgrna_manifest(args.out_sgrna, sgrna_rows, project_dir, raw_fastq_dir)

    print(f"Wrote {args.out_demux_barcodes} ({len(demux_rows)} RNA/ATAC samples)")
    print(f"Wrote {args.out_sgrna} ({len(sgrna_rows)} sgRNA samples)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
