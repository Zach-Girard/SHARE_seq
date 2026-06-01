#!/usr/bin/env python3
"""
sgRNA analysis entrypoint.

Stages sgRNA.tsv in the work directory, then runs a user-provided shell script
(scripts/sgrna_run.sh by default) that may invoke one or more commands.

Environment passed to the runner:
  SGRNA_MANIFEST  - absolute path to sgRNA.tsv in the work directory
  SGRNA_OUT_DIR   - absolute path to the output directory (sgRNA/)
  PROJECT_DIR     - Nextflow project directory
"""

from __future__ import annotations

import argparse
import csv
import os
import shutil
import subprocess
import sys


def _count_samples(manifest: str) -> int:
    with open(manifest, newline="", errors="replace") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        return sum(1 for row in reader if (row.get("sample_name") or "").strip())


def main() -> int:
    p = argparse.ArgumentParser(description="Run sgRNA analysis via a user shell script.")
    p.add_argument("--manifest", required=True, help="Path to sgRNA.tsv")
    p.add_argument("--out-dir", default="sgRNA", help="Output directory")
    p.add_argument(
        "--runner",
        required=True,
        help="Shell script with sgRNA commands (see scripts/sgrna_run.sh)",
    )
    p.add_argument("--project-dir", required=True, help="Nextflow project directory")
    p.add_argument(
        "--barcode-list",
        required=True,
        help="8bp barcode list for rename_fastq.py (--barcode_list)",
    )
    args = p.parse_args()

    if not os.path.isfile(args.manifest):
        print(f"ERROR: manifest not found: {args.manifest}", file=sys.stderr)
        return 1

    n = _count_samples(args.manifest)
    if n == 0:
        print("No sgRNA samples in manifest; nothing to run.")
        return 0

    os.makedirs(args.out_dir, exist_ok=True)

    manifest_abs = os.path.abspath(args.manifest)
    work_manifest = os.path.join(os.getcwd(), "sgRNA.tsv")
    if manifest_abs != os.path.abspath(work_manifest):
        shutil.copy2(manifest_abs, work_manifest)
    shutil.copy2(work_manifest, os.path.join(args.out_dir, "sgRNA.tsv"))

    runner = os.path.abspath(args.runner)
    if not os.path.isfile(runner):
        print(f"ERROR: sgRNA runner not found: {runner}", file=sys.stderr)
        return 1

    env = os.environ.copy()
    env["SGRNA_MANIFEST"] = os.path.abspath(work_manifest)
    env["SGRNA_OUT_DIR"] = os.path.abspath(args.out_dir)
    env["PROJECT_DIR"] = os.path.abspath(args.project_dir)
    env["SGRNA_BARCODE_LIST"] = os.path.abspath(args.barcode_list)

    print(f"Running sgRNA runner: {runner}", flush=True)
    subprocess.run(["bash", runner], check=True, env=env)
    print(f"sgRNA analysis complete for {n} sample(s); outputs under {args.out_dir}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
