#!/usr/bin/env python3
"""
Build sgRNA_run.tsv for run_lsf / downstream steps after per-pool demultiplexing.

Maps each sample to demux/<sample_name>/<sample_name>.R1.fastq.gz.
"""

from __future__ import annotations

import argparse
import csv
import os
import sys


def main() -> int:
    p = argparse.ArgumentParser(description="Build post-demux sgRNA run manifest.")
    p.add_argument("--sgrna-manifest", required=True, help="sgRNA.tsv from BUILD_SAMPLE_MANIFESTS")
    p.add_argument("--demux-dir", required=True, help="Project demux/ directory")
    p.add_argument("--out", default="sgRNA_run.tsv")
    args = p.parse_args()

    if not os.path.isfile(args.sgrna_manifest):
        print(f"ERROR: not found: {args.sgrna_manifest}", file=sys.stderr)
        return 1

    demux_dir = os.path.abspath(args.demux_dir)
    rows_out = []
    with open(args.sgrna_manifest, newline="", errors="replace") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            sample = (row.get("sample_name") or "").strip()
            if not sample:
                continue
            demux_r1 = os.path.join(demux_dir, sample, f"{sample}.R1.fastq.gz")
            if not os.path.isfile(demux_r1):
                print(f"ERROR: demuxed R1 not found: {demux_r1}", file=sys.stderr)
                return 1
            rows_out.append(
                {
                    "fastq": demux_r1,
                    "sample_name": sample,
                    "grna_library_csv": row.get("grna_library_csv", ""),
                    "experimental_group": row.get("experimental_group", ""),
                }
            )

    if not rows_out:
        print("No sgRNA samples to write.")
        return 0

    with open(args.out, "w", newline="") as out:
        w = csv.DictWriter(
            out,
            fieldnames=["fastq", "sample_name", "grna_library_csv", "experimental_group"],
            delimiter="\t",
        )
        w.writeheader()
        w.writerows(rows_out)

    print(f"Wrote {args.out} ({len(rows_out)} samples)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
