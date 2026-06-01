#!/usr/bin/env python3
"""
Stage sgRNA demux FASTQs into the analysis work directory and rewrite manifest paths.

Reads sgRNA_run.tsv (fastq column may point at projectDir/sgRNA/demux/...) and copies
each R1 into <out_dir>/demux/<sample>/<sample>.R1.fastq.gz so downstream steps use
local paths inside the Nextflow task work dir.
"""

from __future__ import annotations

import argparse
import csv
import os
import shutil
import sys


def main() -> int:
    p = argparse.ArgumentParser(description="Stage sgRNA demux FASTQs for SGRNA_ANALYSIS.")
    p.add_argument("--manifest", required=True, help="Input sgRNA_run.tsv")
    p.add_argument("--out-dir", default=".", help="Analysis work directory")
    p.add_argument("--out-manifest", default="sgRNA_run.local.tsv", help="Rewritten manifest")
    args = p.parse_args()

    if not os.path.isfile(args.manifest):
        print(f"ERROR: manifest not found: {args.manifest}", file=sys.stderr)
        return 1

    out_dir = os.path.abspath(args.out_dir)
    demux_root = os.path.join(out_dir, "demux")
    rows_out = []

    with open(args.manifest, newline="", errors="replace") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            sample = (row.get("sample_name") or "").strip()
            src = (row.get("fastq") or "").strip()
            if not sample or not src:
                continue
            if not os.path.isfile(src):
                print(f"ERROR: demuxed R1 not found: {src}", file=sys.stderr)
                return 1

            dest_dir = os.path.join(demux_root, sample)
            os.makedirs(dest_dir, exist_ok=True)
            dest = os.path.join(dest_dir, f"{sample}.R1.fastq.gz")
            if os.path.abspath(src) != os.path.abspath(dest):
                shutil.copy2(src, dest)

            rows_out.append(
                {
                    "fastq": dest,
                    "sample_name": sample,
                    "grna_library_csv": row.get("grna_library_csv", ""),
                    "experimental_group": row.get("experimental_group", ""),
                }
            )

    if not rows_out:
        print("ERROR: no sgRNA samples in manifest.", file=sys.stderr)
        return 1

    out_manifest = os.path.join(out_dir, args.out_manifest)
    with open(out_manifest, "w", newline="") as out:
        w = csv.DictWriter(
            out,
            fieldnames=["fastq", "sample_name", "grna_library_csv", "experimental_group"],
            delimiter="\t",
        )
        w.writeheader()
        w.writerows(rows_out)

    print(f"Staged {len(rows_out)} sample(s) under {demux_root}")
    print(f"Wrote {out_manifest}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
