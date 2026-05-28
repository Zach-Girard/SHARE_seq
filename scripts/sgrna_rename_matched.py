#!/usr/bin/env python3
"""
Run rename_fastq.py (SHARE-seq step 2) on run_lsf outputs to produce matched FASTQs.

Expects per-sample inputs like {sample_name}_C3.fastq.gz (same file used for R1 and R2).
Writes {sample_name}.matched.R1.fastq.gz under sgRNA/<sample_name>/.
"""

from __future__ import annotations

import argparse
import csv
import os
import shutil
import subprocess
import sys


def infer_c3_fastq(sample_name: str, fastq_manifest: str, search_dirs: list[str]) -> str:
    """Locate {sample}_C3.fastq.gz produced by run_lsf.py."""
    name = f"{sample_name}_C3.fastq.gz"
    candidates: list[str] = []
    for d in search_dirs:
        if d:
            candidates.append(os.path.join(d, name))
    if fastq_manifest:
        candidates.append(os.path.join(os.path.dirname(fastq_manifest), name))

    seen: set[str] = set()
    for path in candidates:
        path = os.path.abspath(path)
        if path in seen:
            continue
        seen.add(path)
        if os.path.isfile(path):
            return path

    return os.path.abspath(candidates[0]) if candidates else name


def rename_sample(
    sample_id: str,
    c3_fastq: str,
    barcode_list: str,
    rename_script: str,
    sample_out: str,
    error: int,
) -> None:
    if not os.path.isfile(c3_fastq):
        raise FileNotFoundError(f"C3 FASTQ not found for {sample_id}: {c3_fastq}")
    if not os.path.isfile(barcode_list):
        raise FileNotFoundError(f"Barcode list not found: {barcode_list}")

    os.makedirs(sample_out, exist_ok=True)
    work_dir = os.path.dirname(c3_fastq) or os.getcwd()

    cmd = [
        sys.executable,
        rename_script,
        "-r1",
        c3_fastq,
        "-r2",
        c3_fastq,
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
    p.add_argument("--manifest", required=True, help="sgRNA.tsv")
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
    manifest = os.path.abspath(args.manifest)
    search_dirs = [os.getcwd(), out_dir, os.path.dirname(manifest)]

    n = 0
    with open(manifest, newline="", errors="replace") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            sample = (row.get("sample_name") or "").strip()
            if not sample:
                continue
            fastq = (row.get("fastq") or "").strip()
            c3 = infer_c3_fastq(sample, fastq, search_dirs)
            sample_out = os.path.join(out_dir, sample)
            try:
                rename_sample(
                    sample,
                    c3,
                    os.path.abspath(args.barcode_list),
                    rename_script,
                    sample_out,
                    args.error,
                )
            except FileNotFoundError as exc:
                print(f"ERROR: {exc}", file=sys.stderr)
                return 1
            n += 1

    if n == 0:
        print("No sgRNA samples in manifest.")
    else:
        print(f"Rename complete for {n} sample(s).")
    return 0


if __name__ == "__main__":
    sys.exit(main())
