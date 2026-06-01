#!/usr/bin/env python3
"""
sgRNA analysis entrypoint.

Runs scripts/sgrna_run.sh (or --runner) on sgRNA_run.tsv from BUILD_SGRNA_RUN_MANIFEST.

Environment passed to the runner:
  SGRNA_MANIFEST  - absolute path to sgRNA_run.tsv in the work directory
  SGRNA_OUT_DIR   - absolute path to the output directory (sgRNA/)
  PROJECT_DIR     - Nextflow project directory
"""

from __future__ import annotations

import argparse
import csv
import os
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

    out_dir = os.path.abspath(args.out_dir)
    project_dir = os.path.abspath(args.project_dir)
    stage_script = os.path.join(project_dir, "scripts", "stage_sgrna_demux_for_analysis.py")
    local_manifest = os.path.join(out_dir, "sgRNA_run.local.tsv")

    if os.path.isfile(stage_script):
        stage = subprocess.run(
            [
                sys.executable,
                stage_script,
                "--manifest",
                os.path.abspath(args.manifest),
                "--out-dir",
                out_dir,
                "--out-manifest",
                os.path.basename(local_manifest),
            ],
            check=False,
            capture_output=True,
            text=True,
        )
        if stage.returncode != 0:
            print(stage.stdout, file=sys.stderr)
            print(stage.stderr, file=sys.stderr)
            print(
                f"ERROR: failed to stage sgRNA demux FASTQs (exit {stage.returncode})",
                file=sys.stderr,
            )
            return stage.returncode
        work_manifest = local_manifest
    else:
        work_manifest = os.path.abspath(args.manifest)

    runner = (args.runner or "").strip()
    if not runner or runner.lower() == "null":
        runner = os.path.join(os.path.abspath(args.project_dir), "scripts", "sgrna_run.sh")
    runner = os.path.abspath(runner)
    if not os.path.isfile(runner):
        print(f"ERROR: sgRNA runner not found: {runner}", file=sys.stderr)
        return 1

    env = os.environ.copy()
    env["SGRNA_MANIFEST"] = os.path.abspath(work_manifest)
    env["SGRNA_OUT_DIR"] = out_dir
    env["PROJECT_DIR"] = project_dir
    env["SGRNA_BARCODE_LIST"] = os.path.abspath(args.barcode_list)

    print(f"Running sgRNA runner: {runner}", flush=True)
    print(f"Using manifest: {work_manifest}", flush=True)
    run = subprocess.run(["bash", runner], check=False, env=env, capture_output=True, text=True)
    if run.stdout:
        print(run.stdout, end="" if run.stdout.endswith("\n") else "\n", flush=True)
    if run.stderr:
        print(run.stderr, end="" if run.stderr.endswith("\n") else "\n", file=sys.stderr, flush=True)
    if run.returncode != 0:
        print(f"ERROR: sgRNA runner failed (exit {run.returncode})", file=sys.stderr)
        return run.returncode
    print(f"sgRNA analysis complete for {n} sample(s); outputs under {args.out_dir}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
