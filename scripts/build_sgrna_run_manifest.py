#!/usr/bin/env python3
"""
Build sgRNA_run.tsv for rename + gRNA counting after cutadapt demultiplexing.

Maps each sample to sgRNA/demux/<sample_name>/<sample_name>.R1.fastq.gz.
"""

from __future__ import annotations

import argparse
import csv
import os
import sys
from typing import List, Optional


def _staged_demux_for_sample(sample: str, staged_paths: List[str]) -> Optional[str]:
    suffix = f"/{sample}/{sample}.R1.fastq.gz"
    for path in staged_paths:
        path = os.path.abspath(path)
        if path.endswith(suffix) or path.endswith(f"{sample}.R1.fastq.gz"):
            return path
    return None


def main() -> int:
    p = argparse.ArgumentParser(description="Build post-demux sgRNA run manifest.")
    p.add_argument("--sgrna-manifest", required=True, help="sgRNA.tsv from BUILD_SAMPLE_MANIFESTS")
    p.add_argument("--project-dir", required=True, help="Nextflow launch / publish directory")
    p.add_argument(
        "--demux-r1",
        action="append",
        default=[],
        help="Staged demux R1 paths from SGRNA_DEMULTIPLEX_CUTADAPT (for validation before publish)",
    )
    p.add_argument("--out", default="sgRNA_run.tsv")
    args = p.parse_args()

    if not os.path.isfile(args.sgrna_manifest):
        print(f"ERROR: not found: {args.sgrna_manifest}", file=sys.stderr)
        return 1

    project_dir = os.path.abspath(args.project_dir)
    demux_dir = os.path.join(project_dir, "sgRNA", "demux")
    staged = [os.path.abspath(p) for p in args.demux_r1 if p]
    rows_out = []
    with open(args.sgrna_manifest, newline="", errors="replace") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            sample = (row.get("sample_name") or "").strip()
            if not sample:
                continue
            demux_r1 = os.path.join(demux_dir, sample, f"{sample}.R1.fastq.gz")
            if not os.path.isfile(demux_r1) and not _staged_demux_for_sample(sample, staged):
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
