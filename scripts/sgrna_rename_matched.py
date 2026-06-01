#!/usr/bin/env python3
"""
Run rename_fastq.py on cutadapt-demuxed R1 FASTQs to produce matched FASTQs.

Reads the `fastq` column from sgRNA_run.tsv (paths under sgRNA/demux/<sample>/).
Writes {sample_name}.matched.R1.fastq.gz under sgRNA/<sample_name>/.
"""

from __future__ import annotations

import argparse
import csv
import os
import shutil
import subprocess
import sys


def resolve_demux_fastq(fastq_manifest: str, sample_name: str) -> str:
    """Resolve demuxed R1 path from manifest fastq column."""
    if not fastq_manifest:
        raise FileNotFoundError(f"No fastq path in manifest for sample {sample_name}")
    path = os.path.abspath(fastq_manifest)
    if not os.path.isfile(path):
        raise FileNotFoundError(f"Demuxed R1 not found for {sample_name}: {path}")
    return path


def rename_sample(
    sample_id: str,
    input_fastq: str,
    barcode_list: str,
    rename_script: str,
    sample_out: str,
    error: int,
) -> None:
    if not os.path.isfile(input_fastq):
        raise FileNotFoundError(f"Input FASTQ not found for {sample_id}: {input_fastq}")
    if not os.path.isfile(barcode_list):
        raise FileNotFoundError(f"Barcode list not found: {barcode_list}")

    os.makedirs(sample_out, exist_ok=True)
    work_dir = os.path.dirname(input_fastq) or os.getcwd()

    cmd = [
        sys.executable,
        rename_script,
        "-r1",
        input_fastq,
        "-r2",
        input_fastq,
        "--sample_ID",
        sample_id,
        "--barcode_list",
        barcode_list,
        "--error",
        str(error),
    ]

    print(f"Renaming barcodes for {sample_id}: {' '.join(cmd)}", flush=True)
    subprocess.run(cmd, check=True, cwd=work_dir)

    for suffix in (
        ".matched.R1.fastq.gz",
        ".matched.R2.fastq.gz",
        ".junk.R1.fastq.gz",
        ".junk.R2.fastq.gz",
        ".total_number_reads.tsv",
    ):
        src = os.path.join(work_dir, f"{sample_id}{suffix}")
        if os.path.isfile(src):
            dest = os.path.join(sample_out, os.path.basename(src))
            if os.path.abspath(src) != os.path.abspath(dest):
                shutil.copy2(src, dest)


def main() -> int:
    p = argparse.ArgumentParser(description="SHARE-seq rename (matched FASTQs) for sgRNA samples.")
    p.add_argument("--manifest", required=True, help="sgRNA_run.tsv")
    p.add_argument("--out-dir", required=True, help="Base output directory (sgRNA/)")
    p.add_argument("--barcode-list", required=True, help="8bp barcode whitelist (e.g. barcode1.list)")
    p.add_argument(
        "--rename-script",
        default=None,
        help="Path to rename_fastq.py (default: alongside this script)",
    )
    p.add_argument("--error", type=int, default=1, help="Allowed barcode mismatches")
    args = p.parse_args()

    rename_script = args.rename_script or os.path.join(
        os.path.dirname(os.path.abspath(__file__)), "rename_fastq.py"
    )
    if not os.path.isfile(rename_script):
        print(f"ERROR: rename script not found: {rename_script}", file=sys.stderr)
        return 1

    out_dir = os.path.abspath(args.out_dir)

    n = 0
    with open(args.manifest, newline="", errors="replace") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            sample = (row.get("sample_name") or "").strip()
            if not sample:
                continue
            fastq = (row.get("fastq") or "").strip()
            try:
                input_r1 = resolve_demux_fastq(fastq, sample)
            except FileNotFoundError as exc:
                print(f"ERROR: {exc}", file=sys.stderr)
                return 1
            sample_out = os.path.join(out_dir, sample)
            rename_sample(
                sample,
                input_r1,
                os.path.abspath(args.barcode_list),
                rename_script,
                sample_out,
                args.error,
            )
            n += 1

    if n == 0:
        print("No sgRNA samples in manifest.")
    else:
        print(f"Rename complete for {n} sample(s).")
    return 0


if __name__ == "__main__":
    sys.exit(main())
