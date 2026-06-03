#!/usr/bin/env python3
"""
sgRNA analysis: stage demux FASTQs, SHARE-seq rename_fastq, match_gRNA, QC summary.

Input: sgRNA_run.tsv from BUILD_SGRNA_RUN_MANIFEST (fastq, fastq_r2, sample_name, ...).
For sgRNA, fastq_r2 is intentionally allowed to equal fastq (R1-only mode).
"""

from __future__ import annotations

import argparse
import csv
import glob
import gzip
import os
import shutil
import statistics
import subprocess
import sys
from typing import Dict, List, Optional

MANIFEST_FIELDS = [
    "fastq",
    "fastq_r2",
    "sample_name",
    "grna_library_csv",
    "experimental_group",
]

QC_SUMMARY_FIELDS = [
    "sample_name",
    "experimental_group",
    "assigned_gRNA_reads",
    "cells",
    "cells_with_gRNA",
    "library_size",
    "detected_gRNAs",
    "pct_gRNAs_detected",
    "mean_gRNA_counts_per_cell",
    "median_gRNA_counts_per_cell",
    "top_gRNA",
    "top_gRNA_count",
]


def _pipeline_env(share_seq_pipeline_dir: str) -> dict[str, str]:
    env = os.environ.copy()
    prev = env.get("PYTHONPATH", "")
    env["PYTHONPATH"] = share_seq_pipeline_dir + (os.pathsep + prev if prev else "")
    return env


def _require_utils(share_seq_pipeline_dir: str) -> None:
    utils_py = os.path.join(share_seq_pipeline_dir, "utils.py")
    if not os.path.isfile(utils_py):
        raise FileNotFoundError(
            f"utils.py not found: {utils_py}. Set --share-seq-pipeline-dir correctly."
        )


def _staged_demux_for_sample(
    sample: str, staged_paths: List[str], suffix: str
) -> Optional[str]:
    for path in staged_paths:
        path = os.path.abspath(path)
        if path.endswith(f"/{sample}/{sample}{suffix}") or path.endswith(f"{sample}{suffix}"):
            return path
    return None


def _resolve_src_fastq(
    sample: str,
    label: str,
    manifest_path: str,
    staged_paths: List[str],
    out_dir: str,
) -> str:
    if os.path.isfile(manifest_path):
        return os.path.abspath(manifest_path)
    suffix = ".R1.fastq.gz" if label == "R1" else ".R2.fastq.gz"
    staged = _staged_demux_for_sample(sample, staged_paths, suffix)
    if staged and os.path.isfile(staged):
        return staged
    fallback = os.path.join(out_dir, "demux", sample, f"{sample}{suffix}")
    if os.path.isfile(fallback):
        return os.path.abspath(fallback)
    raise FileNotFoundError(f"Demux {label} not found for {sample}: {manifest_path}")


def stage_demux_manifest(
    manifest: str,
    out_dir: str,
    local_manifest_name: str = "sgRNA_run.local.tsv",
    staged_r1: Optional[List[str]] = None,
    staged_r2: Optional[List[str]] = None,
) -> str:
    """Copy demux R1 into work dir and return path to rewritten manifest."""
    out_dir = os.path.abspath(out_dir)
    demux_root = os.path.join(out_dir, "demux")
    rows_out: List[dict[str, str]] = []
    extra_r1 = [os.path.abspath(p) for p in (staged_r1 or []) if p]
    extra_r2 = [os.path.abspath(p) for p in (staged_r2 or []) if p]

    with open(manifest, newline="", errors="replace") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            sample = (row.get("sample_name") or "").strip()
            src_r1 = (row.get("fastq") or "").strip()
            src_r2 = (row.get("fastq_r2") or "").strip() or src_r1
            if not sample or not src_r1:
                raise ValueError(
                    f"sample {sample or '?'} missing sample_name or fastq"
                )
            src_r1 = _resolve_src_fastq(sample, "R1", src_r1, extra_r1, out_dir)
            if src_r2 != src_r1:
                src_r2 = _resolve_src_fastq(sample, "R2", src_r2, extra_r2, out_dir)
            else:
                src_r2 = src_r1

            dest_dir = os.path.join(demux_root, sample)
            os.makedirs(dest_dir, exist_ok=True)
            dest_r1 = os.path.join(dest_dir, f"{sample}.R1.fastq.gz")
            if os.path.abspath(src_r1) != os.path.abspath(dest_r1):
                shutil.copy2(src_r1, dest_r1)

            rows_out.append(
                {
                    "fastq": dest_r1,
                    "fastq_r2": dest_r1,
                    "sample_name": sample,
                    "grna_library_csv": row.get("grna_library_csv", ""),
                    "experimental_group": row.get("experimental_group", ""),
                }
            )

    if not rows_out:
        raise ValueError("no sgRNA samples in manifest")

    out_manifest = os.path.join(out_dir, local_manifest_name)
    with open(out_manifest, "w", newline="") as out:
        w = csv.DictWriter(out, fieldnames=MANIFEST_FIELDS, delimiter="\t")
        w.writeheader()
        w.writerows(rows_out)

    print(f"Staged {len(rows_out)} sample(s) under {demux_root}", flush=True)
    print(f"Wrote {out_manifest}", flush=True)
    return out_manifest


def resolve_demux_fastq(path: str, sample_name: str, label: str) -> str:
    if not path:
        raise FileNotFoundError(f"No {label} in manifest for sample {sample_name}")
    resolved = os.path.abspath(path)
    if not os.path.isfile(resolved):
        raise FileNotFoundError(f"Demuxed {label} not found for {sample_name}: {resolved}")
    return resolved


def count_reads_gz(path: str) -> int:
    if not os.path.isfile(path):
        return 0
    try:
        with gzip.open(path, "rt") as fh:
            return sum(1 for _ in fh) // 4
    except OSError:
        return 0


def infer_matched_fastq(
    fastq_path: str, sample_name: str, search_dirs: List[str], out_dir: str
) -> str:
    """Locate {sample}.matched.R1.fastq.gz after rename_fastq."""
    matched = f"{sample_name}.matched.R1.fastq.gz"
    candidates: List[str] = [os.path.join(out_dir, sample_name, matched)]
    for d in search_dirs:
        candidates.append(os.path.join(d, matched))

    base = os.path.basename(fastq_path)
    if base.endswith(".R1.fastq.gz"):
        stem = base[: -len(".R1.fastq.gz")]
    elif base.endswith(".fastq.gz"):
        stem = base[: -len(".fastq.gz")]
    else:
        stem = os.path.splitext(base)[0]
    for name in (stem, sample_name):
        alt = f"{name}.matched.R1.fastq.gz"
        candidates.append(os.path.join(os.path.dirname(fastq_path), alt))
        for d in search_dirs:
            candidates.append(os.path.join(d, alt))

    seen: set[str] = set()
    for path in candidates:
        path = os.path.abspath(path)
        if path in seen:
            continue
        seen.add(path)
        if os.path.isfile(path):
            return path
    return os.path.abspath(candidates[0])


RENAME_OUTPUT_SUFFIXES = (
    ".matched.R1.fastq.gz",
    ".junk.R1.fastq.gz",
    ".total_number_reads.tsv",
)


def _rename_output_files(work_dir: str, sample_id: str) -> List[str]:
    """Paths produced by SHARE-seq rename in work_dir (demux/<sample>/)."""
    found: List[str] = []
    seen: set[str] = set()
    for suffix in RENAME_OUTPUT_SUFFIXES:
        path = os.path.join(work_dir, f"{sample_id}{suffix}")
        if os.path.isfile(path):
            abspath = os.path.abspath(path)
            if abspath not in seen:
                seen.add(abspath)
                found.append(abspath)
    for pattern in (
        "*.matched.R1.fastq.gz",
        "*.junk.R1.fastq.gz",
        "*.total_number_reads.tsv",
    ):
        for path in glob.glob(os.path.join(work_dir, pattern)):
            abspath = os.path.abspath(path)
            if abspath not in seen:
                seen.add(abspath)
                found.append(abspath)
    return found


def _find_matched_r1(work_dir: str, sample_out: str, sample_id: str) -> str:
    candidates: List[str] = [
        os.path.join(sample_out, f"{sample_id}.matched.R1.fastq.gz"),
        os.path.join(work_dir, f"{sample_id}.matched.R1.fastq.gz"),
    ]
    for directory in (sample_out, work_dir):
        candidates.extend(glob.glob(os.path.join(directory, "*.matched.R1.fastq.gz")))
    seen: set[str] = set()
    for path in candidates:
        path = os.path.abspath(path)
        if path in seen:
            continue
        seen.add(path)
        if os.path.isfile(path) and count_reads_gz(path) > 0:
            return path
    return os.path.abspath(candidates[0]) if candidates else ""


def rename_sample(
    sample_id: str,
    input_r1: str,
    input_r2: str,
    barcode_list: str,
    rename_script: str,
    share_seq_pipeline_dir: str,
    sample_out: str,
    error: int,
) -> None:
    _require_utils(share_seq_pipeline_dir)
    if not os.path.isfile(rename_script):
        raise FileNotFoundError(f"rename script not found: {rename_script}")

    os.makedirs(sample_out, exist_ok=True)
    work_dir = os.path.dirname(input_r1) or os.getcwd()
    env = _pipeline_env(share_seq_pipeline_dir)

    # R1-only sgRNA mode: pass demuxed R1 as both -r1 and -r2.
    cmd = [
        sys.executable,
        os.path.abspath(rename_script),
        "-r1",
        os.path.abspath(input_r1),
        "-r2",
        os.path.abspath(input_r1),
        "--sample_ID",
        sample_id,
        "--barcode_list",
        os.path.abspath(barcode_list),
        "--error",
        str(error),
    ]
    print(f"Step 1 rename_fastq [{sample_id}]: {' '.join(cmd)}", flush=True)
    subprocess.run(cmd, check=True, cwd=work_dir, env=env)

    # Ignore/delete R2 outputs from rename in sgRNA mode.
    for r2_name in (
        f"{sample_id}.matched.R2.fastq.gz",
        f"{sample_id}.junk.R2.fastq.gz",
    ):
        for d in (work_dir, sample_out):
            p = os.path.join(d, r2_name)
            if os.path.isfile(p):
                try:
                    os.remove(p)
                except OSError:
                    pass

    for src in _rename_output_files(work_dir, sample_id):
        dest = os.path.join(sample_out, os.path.basename(src))
        if os.path.abspath(src) != os.path.abspath(dest):
            shutil.copy2(src, dest)

    matched_r1 = _find_matched_r1(work_dir, sample_out, sample_id)
    if count_reads_gz(matched_r1) == 0:
        listing = ", ".join(sorted(os.listdir(work_dir))) if os.path.isdir(work_dir) else "?"
        raise RuntimeError(
            f"rename_fastq wrote no matched reads for {sample_id} "
            f"(work_dir={work_dir}: {listing}; check barcode headers in {input_r1})"
        )


def run_match_grna(
    match_script: str,
    pipeline_dir: str,
    fastq: str,
    library: str,
    matrix_out: str,
    start: int,
) -> None:
    _require_utils(pipeline_dir)
    if not os.path.isfile(match_script):
        raise FileNotFoundError(f"match_gRNA.py not found: {match_script}")

    env = _pipeline_env(pipeline_dir)
    cmd = [
        sys.executable,
        os.path.abspath(match_script),
        "-f",
        os.path.abspath(fastq),
        "-g",
        os.path.abspath(library),
        "-o",
        os.path.abspath(matrix_out),
        "-s",
        str(start),
    ]
    print(f"Step 2 match_gRNA: {' '.join(cmd)}", flush=True)
    subprocess.run(cmd, check=True, cwd=pipeline_dir, env=env)


def normalize_cell_matrix(matrix_path: str) -> None:
    rows: List[list] = []
    with open(matrix_path, newline="", errors="replace") as fh:
        reader = csv.reader(fh)
        header = next(reader, None)
        if not header:
            return
        first = header[0].strip()
        if first in ("cell_barcode", "cell"):
            if first == "cell_barcode":
                return
            header[0] = "cell_barcode"
        else:
            header = ["cell_barcode"] + header
        rows.append(header)
        for row in reader:
            if row:
                rows.append(row)
    with open(matrix_path, "w", newline="") as out:
        csv.writer(out).writerows(rows)


def write_grna_counts_final(matrix_path: str, out_path: str) -> None:
    totals: dict[str, int] = {}
    with open(matrix_path, newline="", errors="replace") as fh:
        reader = csv.reader(fh)
        header = next(reader, None)
        if not header or len(header) < 2:
            return
        for row in reader:
            if not row:
                continue
            for seq, val in zip(header[1:], row[1:]):
                try:
                    totals[seq] = totals.get(seq, 0) + int(float(val or 0))
                except (TypeError, ValueError):
                    continue
    with open(out_path, "w", newline="") as out:
        writer = csv.writer(out)
        writer.writerow(["sgRNA_sequence", "count"])
        for seq in sorted(totals):
            writer.writerow([seq, totals[seq]])


def load_matrix_stats(path: str) -> tuple[int, int, int, float, float]:
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
                except (TypeError, ValueError):
                    vals.append(0)
            row_sums.append(sum(vals))
    cells = len(row_sums)
    cells_with_grna = sum(1 for v in row_sums if v > 0)
    assigned = sum(row_sums)
    mean_counts = statistics.mean(row_sums) if row_sums else 0.0
    median_counts = statistics.median(row_sums) if row_sums else 0.0
    return cells, cells_with_grna, assigned, mean_counts, median_counts


def load_grna_counts_stats(path: str) -> tuple[int, int, str, int]:
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
            try:
                count = int(float(row.get("count", "0") or 0))
            except (TypeError, ValueError):
                count = 0
            if count > 0:
                detected += 1
            if count > top_count:
                top_grna = row.get("sgRNA_sequence", "")
                top_count = count
    return library_size, detected, top_grna, top_count


def _fmt_qc(value: object) -> str:
    if value is None:
        return ""
    if isinstance(value, float):
        return f"{value:.2f}"
    return str(value)


def process_sample(
    row: dict[str, str],
    out_dir: str,
    barcode_list: str,
    rename_script: str,
    match_script: str,
    pipeline_dir: str,
    rename_error: int,
    match_start: int,
) -> dict[str, object]:
    sample = (row.get("sample_name") or "").strip()
    grna_library = (row.get("grna_library_csv") or "").strip()
    input_r1 = resolve_demux_fastq((row.get("fastq") or "").strip(), sample, "fastq (R1)")
    raw_r2 = (row.get("fastq_r2") or "").strip()
    input_r2 = (
        resolve_demux_fastq(raw_r2, sample, "fastq_r2 (R2)")
        if raw_r2 and raw_r2 != (row.get("fastq") or "").strip()
        else input_r1
    )

    if not grna_library:
        raise ValueError(f"sample {sample} missing grna_library_csv")
    if not os.path.isfile(grna_library):
        raise FileNotFoundError(f"gRNA library not found: {grna_library}")

    sample_out = os.path.join(out_dir, sample)
    rename_sample(
        sample,
        input_r1,
        input_r2,
        barcode_list,
        rename_script,
        pipeline_dir,
        sample_out,
        rename_error,
    )

    search_dirs = [os.getcwd(), out_dir]
    matched = infer_matched_fastq((row.get("fastq") or "").strip(), sample, search_dirs, out_dir)
    if not os.path.isfile(matched):
        raise FileNotFoundError(
            f"matched FASTQ not found for {sample} after rename (expected under {sample_out})"
        )

    matrix_out = os.path.join(sample_out, f"final_{sample}.gRNA.count.csv")
    counts_out = os.path.join(sample_out, "gRNA_counts_final.csv")
    run_match_grna(match_script, pipeline_dir, matched, grna_library, matrix_out, match_start)
    normalize_cell_matrix(matrix_out)
    write_grna_counts_final(matrix_out, counts_out)

    cell_list_out = matrix_out + ".cell_with_gRNA.csv"
    if not os.path.isfile(cell_list_out):
        print(f"WARNING: expected {cell_list_out} from match_gRNA", file=sys.stderr)

    cells, cells_with_grna, assigned, mean_counts, median_counts = load_matrix_stats(matrix_out)
    library_size, detected, top_grna, top_count = load_grna_counts_stats(counts_out)
    pct_guides = (100.0 * detected / library_size) if library_size else None

    print(f"Completed {sample}", flush=True)
    return {
        "sample_name": sample,
        "experimental_group": row.get("experimental_group", ""),
        "assigned_gRNA_reads": assigned,
        "cells": cells,
        "cells_with_gRNA": cells_with_grna,
        "library_size": library_size,
        "detected_gRNAs": detected,
        "pct_gRNAs_detected": pct_guides,
        "mean_gRNA_counts_per_cell": mean_counts,
        "median_gRNA_counts_per_cell": median_counts,
        "top_gRNA": top_grna,
        "top_gRNA_count": top_count,
    }


def write_qc_summary(rows: List[dict[str, object]], out_path: str) -> None:
    with open(out_path, "w", newline="") as out:
        writer = csv.DictWriter(out, fieldnames=QC_SUMMARY_FIELDS, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow({key: _fmt_qc(row.get(key)) for key in QC_SUMMARY_FIELDS})
    print(f"Step 3 wrote {out_path} ({len(rows)} samples)", flush=True)


def load_manifest_rows(manifest: str) -> List[dict[str, str]]:
    with open(manifest, newline="", errors="replace") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def main() -> int:
    p = argparse.ArgumentParser(description="sgRNA demux staging, rename, match_gRNA, QC summary.")
    p.add_argument("--manifest", required=True, help="sgRNA_run.tsv")
    p.add_argument("--out-dir", default=".", help="Output directory (sgRNA/)")
    p.add_argument("--project-dir", required=True, help="Nextflow project directory")
    p.add_argument("--barcode-list", required=True, help="8bp barcode whitelist")
    p.add_argument(
        "--rename-fastq-script",
        required=True,
        help="SHARE_seq rename script (e.g. share_seq_step2_rename_fastq.py)",
    )
    p.add_argument("--match-grna-script", required=True, help="SHARE_seq match_gRNA.py")
    p.add_argument("--share-seq-pipeline-dir", required=True, help="SHARE_seq_pipeline (utils.py)")
    p.add_argument("--match-grna-start", type=int, default=43, help="match_gRNA -s/--start")
    p.add_argument("--rename-error", type=int, default=1, help="rename_fastq --error")
    p.add_argument("--qc-summary-out", default="sgrna_qc_summary.tsv", help="QC summary TSV name")
    p.add_argument(
        "--skip-stage",
        action="store_true",
        help="Use manifest paths as-is (do not copy demux into out-dir)",
    )
    p.add_argument(
        "--demux-r1",
        action="append",
        default=[],
        help="Staged demux R1 from SGRNA_DEMULTIPLEX (when manifest paths are project-dir only)",
    )
    p.add_argument(
        "--demux-r2",
        action="append",
        default=[],
        help="Staged demux R2 from SGRNA_DEMULTIPLEX",
    )
    args = p.parse_args()

    if not os.path.isfile(args.manifest):
        print(f"ERROR: manifest not found: {args.manifest}", file=sys.stderr)
        return 1

    out_dir = os.path.abspath(args.out_dir)
    os.makedirs(out_dir, exist_ok=True)
    pipeline_dir = os.path.abspath(args.share_seq_pipeline_dir)
    rename_script = os.path.abspath(args.rename_fastq_script)
    match_script = os.path.abspath(args.match_grna_script)
    barcode_list = os.path.abspath(args.barcode_list)

    try:
        if args.skip_stage:
            work_manifest = os.path.abspath(args.manifest)
        else:
            work_manifest = stage_demux_manifest(
                args.manifest,
                out_dir,
                staged_r1=args.demux_r1,
                staged_r2=args.demux_r2,
            )

        manifest_rows = load_manifest_rows(work_manifest)
        samples = [r for r in manifest_rows if (r.get("sample_name") or "").strip()]
        if not samples:
            print("No sgRNA samples in manifest; nothing to run.")
            return 0

        qc_rows: List[dict[str, object]] = []
        for row in samples:
            qc_rows.append(
                process_sample(
                    row,
                    out_dir,
                    barcode_list,
                    rename_script,
                    match_script,
                    pipeline_dir,
                    args.rename_error,
                    args.match_grna_start,
                )
            )

        summary_path = os.path.join(out_dir, args.qc_summary_out)
        write_qc_summary(qc_rows, summary_path)
        print(f"sgRNA analysis complete for {len(qc_rows)} sample(s); outputs under {out_dir}")
        return 0

    except (FileNotFoundError, ValueError, RuntimeError, subprocess.CalledProcessError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())
