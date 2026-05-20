#!/usr/bin/env python3
"""
Compare RNA (STARsolo filtered) and ATAC (ArchR) cell barcodes per experimental group.

Requires sample_barcode_file with columns:
  1 = sample name, 3 = RNA|ATAC, 4 = Experimental_Group (optional but required for overlap).
"""

from __future__ import annotations

import argparse
import csv
import os
import re
import sys
from collections import defaultdict
from typing import Dict, List, Optional, Set, Tuple

import matplotlib.pyplot as plt


def parse_args():
    p = argparse.ArgumentParser(description="RNA/ATAC cell barcode overlap per experimental group.")
    p.add_argument("--project-dir", required=True, help="Nextflow project / launch directory.")
    p.add_argument("--sample-barcode-file", required=True, help="Sample metadata TSV/CSV.")
    p.add_argument(
        "--star-alignment-mode",
        default="single",
        choices=["single", "paired"],
        help="STARsolo output root: STARsolo (single) or STARsolo_paired (paired).",
    )
    p.add_argument("--out-dir", default=".", help="Output directory (multiome_overlap/).")
    p.add_argument(
        "--atac-pre-counts-dir",
        default="",
        help="Directory of staged *.atac_cells.pre_dedup.counts.tsv files (from Nextflow).",
    )
    return p.parse_args()


def split_cols(line: str) -> list[str]:
    return [c.strip() for c in re.split(r"\t|,", line)]


def load_sample_metadata(path: str) -> Tuple[Dict[str, str], Dict[str, str]]:
    """Return sample -> type (RNA|ATAC), sample -> experimental group."""
    sample_type: Dict[str, str] = {}
    sample_group: Dict[str, str] = {}
    if not path or not os.path.isfile(path):
        return sample_type, sample_group

    with open(path, "r", errors="replace") as fh:
        for raw in fh:
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            cols = split_cols(line)
            if len(cols) < 3:
                continue
            sample = cols[0]
            stype = cols[2].upper()
            group = cols[3] if len(cols) >= 4 else ""
            if not sample:
                continue
            if sample.lower() in ("sample", "sample_name") or stype in ("TYPE", "SAMPLE_TYPE"):
                continue
            if stype in ("RNA", "ATAC"):
                sample_type[sample] = stype
            if group and group.lower() not in ("group", "experimental_group", "condition"):
                sample_group[sample] = group
    return sample_type, sample_group


def normalize_barcode(bc: str, length: int = 24) -> Optional[str]:
    if bc is None:
        return None
    s = str(bc).strip().upper()
    if not s:
        return None
    # ArchR / STAR may append sample suffix after '#'
    if "#" in s:
        s = s.split("#")[-1]
    s = re.sub(r"[^ACGTN]", "", s)
    if len(s) < length:
        return None
    if len(s) > length:
        s = s[:length]
    return s


def load_rna_barcodes(project_dir: str, sample: str, starsolo_root: str) -> Set[str]:
    rel = os.path.join(
        starsolo_root,
        sample,
        "Solo.out",
        "GeneFull",
        "filtered",
        "barcodes.tsv",
    )
    path = os.path.join(project_dir, rel)
    if not os.path.isfile(path):
        return set()
    out: Set[str] = set()
    with open(path, "r", errors="replace") as fh:
        for raw in fh:
            bc = normalize_barcode(raw.split("\t")[0].split()[0])
            if bc:
                out.add(bc)
    return out


def resolve_atac_pre_counts_path(
    project_dir: str,
    sample: str,
    atac_pre_counts_dir: Optional[str],
    staged_only: bool = False,
) -> Optional[str]:
    fname = f"{sample}.atac_cells.pre_dedup.counts.tsv"
    if atac_pre_counts_dir:
        staged = os.path.join(atac_pre_counts_dir, fname)
        if os.path.isfile(staged):
            return staged
        if staged_only:
            return None
    canonical = os.path.join(project_dir, "ATAC", sample, fname)
    if os.path.isfile(canonical):
        return canonical
    return None


def count_rows_in_counts_tsv(path: str) -> int:
    with open(path, "r", errors="replace") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        if not reader.fieldnames:
            return 0
        return sum(1 for _ in reader)


def load_atac_barcodes(
    project_dir: str,
    sample: str,
    atac_pre_counts_dir: Optional[str] = None,
    staged_only: bool = False,
) -> Tuple[Set[str], str]:
    path = resolve_atac_pre_counts_path(
        project_dir, sample, atac_pre_counts_dir, staged_only=staged_only
    )
    if not path:
        print(
            f"WARNING: missing {sample}.atac_cells.pre_dedup.counts.tsv "
            f"(re-run ESTIMATE_ATAC_CELLS; overlap uses pre-dedup ArchR cells only)",
            file=sys.stderr,
        )
        return set(), ""
    out: Set[str] = set()
    with open(path, "r", errors="replace") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        if not reader.fieldnames:
            return out, path
        bc_col = None
        for name in reader.fieldnames:
            if name and name.strip().lower() == "barcode":
                bc_col = name
                break
        if not bc_col:
            print(
                f"WARNING: no Barcode column in {path}; skipping sample {sample}",
                file=sys.stderr,
            )
            return out, path
        for row in reader:
            bc = normalize_barcode(row.get(bc_col, ""))
            if bc:
                out.add(bc)
    post_path = os.path.join(project_dir, "ATAC", sample, f"{sample}.atac_cells.counts.tsv")
    if os.path.isfile(post_path):
        post_rows = count_rows_in_counts_tsv(post_path)
        pre_rows = count_rows_in_counts_tsv(path)
        if len(out) == post_rows and pre_rows > post_rows:
            print(
                f"WARNING: {sample} unique pre-dedup barcodes ({len(out)}) equals post-dedup "
                f"row count ({post_rows}) but pre file has more rows ({pre_rows}); "
                "check barcode normalization.",
                file=sys.stderr,
            )
    print(
        f"ATAC {sample}: {len(out)} barcodes from pre-dedup {path}",
        file=sys.stderr,
    )
    return out, path


def write_shared_barcodes(path: str, barcodes: Set[str]) -> None:
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    with open(path, "w") as fh:
        for bc in sorted(barcodes):
            fh.write(bc + "\n")


def write_group_summary(path: str, row: dict) -> None:
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    fields = list(row.keys())
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t")
        w.writeheader()
        w.writerow(row)


def plot_overlap(rows: List[dict], out_png: str) -> None:
    if not rows:
        return
    groups = [r["Experimental_Group"] for r in rows]
    rna_n = [int(r["RNA_Cells"]) for r in rows]
    atac_n = [int(r["ATAC_Cells"]) for r in rows]
    shared_n = [int(r["Shared_Cells"]) for r in rows]

    x = range(len(groups))
    width = 0.25
    fig, ax = plt.subplots(figsize=(max(8, len(groups) * 1.4), 5))
    ax.bar([i - width for i in x], rna_n, width=width, label="RNA cells", color="#3b82f6")
    ax.bar(x, atac_n, width=width, label="ATAC cells (pre-dedup)", color="#f59e0b")
    ax.bar([i + width for i in x], shared_n, width=width, label="Shared barcodes", color="#10b981")
    ax.set_xticks(list(x))
    ax.set_xticklabels(groups, rotation=25, ha="right")
    ax.set_ylabel("Cell count")
    ax.set_title("RNA / ATAC cell barcode overlap by experimental group")
    ax.legend()
    fig.tight_layout()
    fig.savefig(out_png, dpi=200)
    plt.close(fig)


def main() -> int:
    args = parse_args()
    project_dir = os.path.abspath(args.project_dir)
    out_dir = os.path.abspath(args.out_dir)
    atac_pre_counts_dir = (args.atac_pre_counts_dir or "").strip()
    staged_only = False
    if atac_pre_counts_dir:
        atac_pre_counts_dir = os.path.abspath(atac_pre_counts_dir)
        staged_only = True
    os.makedirs(out_dir, exist_ok=True)

    starsolo_root = "STARsolo_paired" if args.star_alignment_mode == "paired" else "STARsolo"

    sample_type, sample_group = load_sample_metadata(args.sample_barcode_file)
    if not sample_group:
        print(
            "No Experimental_Group column (column 4) in sample_barcode_file; skipping overlap.",
            file=sys.stderr,
        )
        summary_path = os.path.join(out_dir, "overlap_by_group.tsv")
        with open(summary_path, "w", newline="") as fh:
            w = csv.writer(fh, delimiter="\t")
            w.writerow(
                [
                    "Experimental_Group",
                    "RNA_Samples",
                    "ATAC_Samples",
                    "ATAC_PreDedup_Files",
                    "RNA_Cells",
                    "ATAC_Cells",
                    "Shared_Cells",
                    "RNA_Only",
                    "ATAC_Only",
                    "Pct_RNA_in_ATAC",
                    "Pct_ATAC_in_RNA",
                    "Jaccard",
                ]
            )
            w.writerow(
                ["", "", "", "", "", "", "", "", "", "", "", "No experimental groups defined"]
            )
        return 0

    members_by_group: Dict[str, List[str]] = defaultdict(list)
    for sample, group in sample_group.items():
        members_by_group[group].append(sample)

    summary_rows: List[dict] = []

    for group in sorted(members_by_group):
        members = members_by_group[group]
        rna_samples = sorted(s for s in members if sample_type.get(s) == "RNA")
        atac_samples = sorted(s for s in members if sample_type.get(s) == "ATAC")

        rna_barcodes: Set[str] = set()
        for s in rna_samples:
            rna_barcodes |= load_rna_barcodes(project_dir, s, starsolo_root)

        atac_barcodes: Set[str] = set()
        atac_sources: List[str] = []
        for s in atac_samples:
            sample_barcodes, src = load_atac_barcodes(
                project_dir, s, atac_pre_counts_dir, staged_only=staged_only
            )
            atac_barcodes |= sample_barcodes
            if src:
                atac_sources.append(src)

        shared = rna_barcodes & atac_barcodes
        rna_only = rna_barcodes - atac_barcodes
        atac_only = atac_barcodes - rna_barcodes
        union = rna_barcodes | atac_barcodes

        pct_rna_in_atac = (100.0 * len(shared) / len(rna_barcodes)) if rna_barcodes else ""
        pct_atac_in_rna = (100.0 * len(shared) / len(atac_barcodes)) if atac_barcodes else ""
        jaccard = (len(shared) / len(union)) if union else ""

        row = {
            "Experimental_Group": group,
            "RNA_Samples": ",".join(rna_samples),
            "ATAC_Samples": ",".join(atac_samples),
            "ATAC_PreDedup_Files": ";".join(atac_sources),
            "RNA_Cells": len(rna_barcodes),
            "ATAC_Cells": len(atac_barcodes),
            "Shared_Cells": len(shared),
            "RNA_Only": len(rna_only),
            "ATAC_Only": len(atac_only),
            "Pct_RNA_in_ATAC": f"{pct_rna_in_atac:.2f}" if pct_rna_in_atac != "" else "",
            "Pct_ATAC_in_RNA": f"{pct_atac_in_rna:.2f}" if pct_atac_in_rna != "" else "",
            "Jaccard": f"{jaccard:.4f}" if jaccard != "" else "",
        }
        summary_rows.append(row)

        safe_group = re.sub(r"[^\w\-.]+", "_", group)
        group_dir = os.path.join(out_dir, safe_group)
        write_group_summary(os.path.join(group_dir, "overlap_summary.tsv"), row)
        if shared:
            write_shared_barcodes(os.path.join(group_dir, "shared_barcodes.txt"), shared)

    master_path = os.path.join(out_dir, "overlap_by_group.tsv")
    fieldnames = [
        "Experimental_Group",
        "RNA_Samples",
        "ATAC_Samples",
        "ATAC_PreDedup_Files",
        "RNA_Cells",
        "ATAC_Cells",
        "Shared_Cells",
        "RNA_Only",
        "ATAC_Only",
        "Pct_RNA_in_ATAC",
        "Pct_ATAC_in_RNA",
        "Jaccard",
    ]
    with open(master_path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
        w.writeheader()
        for row in summary_rows:
            w.writerow(row)

    plot_rows = [r for r in summary_rows if r.get("RNA_Cells") and r.get("ATAC_Cells")]
    if plot_rows:
        plot_overlap(plot_rows, os.path.join(out_dir, "overlap_by_group.png"))

    print(f"Wrote {master_path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
