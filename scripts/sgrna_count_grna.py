#!/usr/bin/env python3
"""
Assign reads to gRNA library sequences (fuzzy match) and build count matrices.

Single-sample mode:
  python3 sgrna_count_grna.py --grna-library lib.csv --fastq sample.matched.R1.fastq.gz --sample-id SAMPLE

Batch mode (reads sgRNA.tsv):
  python3 sgrna_count_grna.py --manifest sgRNA.tsv --out-dir sgRNA/
"""

from __future__ import annotations

import argparse
import collections
import csv
import gzip
import os
import sys

try:
    import regex
except ImportError:
    print("ERROR: pip package 'regex' is required (fuzzy matching).", file=sys.stderr)
    sys.exit(1)


def load_library(grna_library: str, max_mismatches: int) -> tuple[list[str], list[tuple[str, regex.Pattern]]]:
    grna_patterns: list[tuple[str, regex.Pattern]] = []
    unique_library_seqs: set[str] = set()

    with open(grna_library, newline="", errors="replace") as f:
        reader = csv.reader(f)
        next(reader, None)  # header
        for row in reader:
            if row and len(row) >= 2:
                seq = row[1].strip().upper()
                if seq and seq not in unique_library_seqs:
                    pattern = regex.compile(r"(?b)(%s){s<=%d}" % (seq, max_mismatches))
                    grna_patterns.append((seq, pattern))
                    unique_library_seqs.add(seq)

    library_ids = sorted(unique_library_seqs)
    return library_ids, grna_patterns


def infer_matched_fastq(
    fastq_path: str, sample_name: str, search_dirs: list[str], out_dir: str
) -> str:
    """Locate {sample}.matched.R1.fastq.gz from rename_fastq (under sgRNA/<sample>/)."""
    candidates: list[str] = []
    matched = f"{sample_name}.matched.R1.fastq.gz"
    candidates.append(os.path.join(out_dir, sample_name, matched))
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


def count_grna(
    grna_library: str,
    fastq_file: str,
    sample_id: str,
    out_dir: str,
    max_mismatches: int = 2,
) -> None:
    os.makedirs(out_dir, exist_ok=True)

    output_csv = os.path.join(out_dir, "gRNA_counts_final.csv")
    output_assigned_fastq = os.path.join(out_dir, "assigned_reads.fastq")
    output_cell_matrix = os.path.join(out_dir, f"final_{sample_id}.gRNA.count.csv")
    output_cell_list = os.path.join(out_dir, f"final_{sample_id}.gRNA.count.csv.cell_with_gRNA.csv")

    print(f"Compiling library {grna_library}...")
    library_ids, grna_patterns = load_library(grna_library, max_mismatches)
    total_gRNA_in_lib = len(library_ids)
    if total_gRNA_in_lib == 0:
        raise ValueError(f"No gRNA sequences loaded from {grna_library}")

    stats: collections.Counter = collections.Counter()
    final_counts: collections.Counter = collections.Counter()
    cell_gRNA_counts: collections.defaultdict = collections.defaultdict(
        lambda: collections.defaultdict(int)
    )

    print(f"Processing {fastq_file}...")
    with gzip.open(fastq_file, "rt") as f, open(output_assigned_fastq, "w") as f_out:
        while True:
            header = f.readline().strip()
            if not header:
                break

            seq_line = f.readline().strip()
            plus_line = f.readline().strip()
            qual_line = f.readline().strip()

            stats["total_reads"] += 1
            cell_barcode = header.split("_")[-1] if "_" in header else "Unknown"

            matches_for_this_read: dict[str, int] = {}
            for lib_seq, pattern in grna_patterns:
                m = pattern.search(seq_line)
                if m:
                    dist = sum(m.fuzzy_counts)
                    if lib_seq not in matches_for_this_read or dist < matches_for_this_read[lib_seq]:
                        matches_for_this_read[lib_seq] = dist

            if not matches_for_this_read:
                stats["unassigned_reads"] += 1
            else:
                best_score = min(matches_for_this_read.values())
                stats[f"mismatch_{best_score}"] += 1

                seqs_at_best_level = [
                    s for s, dist in matches_for_this_read.items() if dist == best_score
                ]

                if len(seqs_at_best_level) > 1:
                    stats["ambiguous_reads"] += 1
                else:
                    stats["assigned_reads"] += 1
                    target_seq = seqs_at_best_level[0]
                    final_counts[target_seq] += 1
                    cell_gRNA_counts[cell_barcode][target_seq] += 1
                    f_out.write(f"{header}\n{seq_line}\n{plus_line}\n{qual_line}\n")

            if stats["total_reads"] % 5000 == 0:
                print(f"Progress: {stats['total_reads']} reads...", end="\r")

    with open(output_cell_matrix, "w", newline="") as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(["cell_barcode"] + library_ids)
        for cb in sorted(cell_gRNA_counts.keys()):
            row = [cb] + [cell_gRNA_counts[cb].get(s, 0) for s in library_ids]
            writer.writerow(row)

    with open(output_cell_list, "w", newline="") as csvfile:
        writer = csv.writer(csvfile, quoting=csv.QUOTE_MINIMAL)
        for cb in sorted(cell_gRNA_counts.keys()):
            seqs_found = sorted(cell_gRNA_counts[cb].keys())
            writer.writerow([cb, ",".join(seqs_found)])

    with open(output_csv, "w", newline="") as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(["sgRNA_sequence", "count"])
        for s in library_ids:
            writer.writerow([s, final_counts.get(s, 0)])

    print("\n" + "=" * 45)
    print(f"sample\t\t\t\t{sample_id}")
    print(f"total_reads\t\t\t\t{stats['total_reads']}")
    print(f"assigned_reads\t\t\t\t{stats['assigned_reads']}")
    print(f"ambiguous_reads\t\t\t\t{stats['ambiguous_reads']}")
    print(f"unassigned_reads\t\t\t\t{stats['unassigned_reads']}")
    print(f"mismatch_0\t\t\t\t{stats['mismatch_0']}")
    print(f"mismatch_1\t\t\t\t{stats['mismatch_1']}")
    print(f"mismatch_2\t\t\t\t{stats['mismatch_2']}")
    print(f"total_unique_gRNA_in_lib\t\t{total_gRNA_in_lib}")
    print(f"gRNA_with_at_least_one_matched_barcode\t{len(final_counts)}")
    sum_mismatches = stats["mismatch_0"] + stats["mismatch_1"] + stats["mismatch_2"]
    sum_assigned_ambig = stats["assigned_reads"] + stats["ambiguous_reads"]
    print("-" * 45)
    print(f"Logic Check: {sum_assigned_ambig} (Assigned+Ambig) == {sum_mismatches} (Mismatches)")
    print("=" * 45)
    print(f"Done! Results under {out_dir}")


def run_from_manifest(manifest: str, out_dir: str, max_mismatches: int) -> int:
    search_dirs = [os.getcwd(), out_dir, os.path.dirname(os.path.abspath(manifest))]
    n = 0
    with open(manifest, newline="", errors="replace") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            sample = (row.get("sample_name") or "").strip()
            grna_library = (row.get("grna_library_csv") or "").strip()
            fastq = (row.get("fastq") or "").strip()
            if not sample:
                continue
            if not grna_library or not fastq:
                print(f"ERROR: sample {sample} missing grna_library_csv or fastq", file=sys.stderr)
                return 1
            if not os.path.isfile(grna_library):
                print(f"ERROR: gRNA library not found: {grna_library}", file=sys.stderr)
                return 1

            matched = infer_matched_fastq(fastq, sample, search_dirs, out_dir)
            if not os.path.isfile(matched):
                print(
                    f"ERROR: matched FASTQ not found for {sample} (tried {matched}). "
                    "Run run_lsf.py and sgrna_rename_matched.py first.",
                    file=sys.stderr,
                )
                return 1

            sample_out = os.path.join(out_dir, sample)
            count_grna(grna_library, matched, sample, sample_out, max_mismatches)
            n += 1

    if n == 0:
        print("No sgRNA samples in manifest.")
    return 0


def main() -> int:
    p = argparse.ArgumentParser(description="Count gRNAs from matched R1 FASTQ.")
    p.add_argument("--manifest", help="sgRNA.tsv for batch processing")
    p.add_argument("--out-dir", default=".", help="Base output directory (batch mode)")
    p.add_argument("--grna-library", help="gRNA library CSV (single-sample mode)")
    p.add_argument("--fastq", help="matched R1 FASTQ (single-sample mode)")
    p.add_argument("--sample-id", help="Sample name for output prefixes (single-sample mode)")
    p.add_argument("--max-mismatches", type=int, default=2)
    args = p.parse_args()

    if args.manifest:
        if not os.path.isfile(args.manifest):
            print(f"ERROR: manifest not found: {args.manifest}", file=sys.stderr)
            return 1
        return run_from_manifest(args.manifest, os.path.abspath(args.out_dir), args.max_mismatches)

    if not all([args.grna_library, args.fastq, args.sample_id]):
        print("ERROR: provide --manifest or (--grna-library, --fastq, --sample-id)", file=sys.stderr)
        return 1

    if not os.path.isfile(args.grna_library):
        print(f"ERROR: {args.grna_library} not found.", file=sys.stderr)
        return 1
    if not os.path.isfile(args.fastq):
        print(f"ERROR: {args.fastq} not found.", file=sys.stderr)
        return 1

    sample_out = os.path.join(args.out_dir, args.sample_id)
    count_grna(args.grna_library, args.fastq, args.sample_id, sample_out, args.max_mismatches)
    return 0


if __name__ == "__main__":
    sys.exit(main())
