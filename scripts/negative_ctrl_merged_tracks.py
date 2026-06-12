#!/usr/bin/env python3
"""
Post-pipeline: merged ATAC/RNA BAM + BigWig for cells with negative-control gRNAs.

Reads guide labels from the gRNA library CSV (default label: negative_ctrl), finds
matching cells in final_<sample>.gRNA.count.csv, and writes merged tracks.
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import re
import sys
from dataclasses import dataclass
from typing import Dict, List, Optional, Set, Tuple

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
if SCRIPT_DIR not in sys.path:
    sys.path.insert(0, SCRIPT_DIR)

from barcode_utils import normalize_barcode
from cell_overlap_by_group import (
    load_atac_barcodes,
    load_rna_barcodes,
    load_sample_metadata,
)
from grna_cell_tracks import (
    SampleSet,
    atac_bam_path,
    build_chrom_sizes,
    ensure_bam_index,
    extract_merged_bam,
    grna_count_matrix_path,
    make_bigwig,
    parse_nextflow_config,
    resolve_effective_genome_size,
    resolve_reference_fasta,
    rna_bam_path,
    write_barcode_list,
)

SEQUENCE_COL_CANDIDATES = (
    "sequence",
    "sgRNA_sequence",
    "sgrna_sequence",
    "grna_sequence",
    "guide_sequence",
    "grna",
    "sgrna",
    "protospacer",
)
NAME_COL_CANDIDATES = (
    "guide",
    "guide_id",
    "name",
    "id",
    "grna_id",
    "sgrna_id",
    "guide_name",
)


@dataclass
class NegativeGuide:
    name: str
    sequence: str
    label: str
    matrix_sequence: str = ""


@dataclass
class NegativeCtrlCell:
    cell_barcode: str
    total_count: int
    guide_hits: Dict[str, int]
    in_atac_call: bool = False
    in_rna_call: bool = False


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Merged ATAC/RNA BAM and BigWig for cells with negative-control gRNAs."
    )
    p.add_argument("--project-dir", required=True, help="Nextflow project / launch directory.")
    p.add_argument(
        "--sample-barcode-file",
        default="",
        help="Sample metadata TSV/CSV (required unless sample IDs are explicit).",
    )
    p.add_argument(
        "--experimental-group",
        default="",
        help="Experimental_Group to resolve sgRNA/ATAC/RNA samples.",
    )
    p.add_argument("--sgrna-sample", default="", help="Explicit sgRNA sample ID.")
    p.add_argument("--atac-sample", default="", help="Explicit ATAC sample ID.")
    p.add_argument("--rna-sample", default="", help="Explicit RNA sample ID.")
    p.add_argument(
        "--grna-library-csv",
        default="",
        help="gRNA library CSV (default: from manifests/sgRNA.tsv for the sgRNA sample).",
    )
    p.add_argument(
        "--control-label",
        default="negative_ctrl",
        help="Guide label in the library CSV marking negative controls (default: negative_ctrl).",
    )
    p.add_argument("--out-dir", required=True, help="Output directory.")
    p.add_argument(
        "--star-alignment-mode",
        default="single",
        choices=["single", "paired"],
        help="STARsolo output root: STARsolo or STARsolo_paired.",
    )
    p.add_argument(
        "--min-grna-count",
        type=int,
        default=1,
        help="Minimum total negative-control gRNA counts per cell (sum across controls).",
    )
    p.add_argument(
        "--modality",
        default="both",
        choices=["atac", "rna", "both"],
        help="Which modalities to extract (default: both).",
    )
    p.add_argument(
        "--require-modality-call",
        action="store_true",
        help="Keep only cells also called in ATAC (ArchR) and RNA (STARsolo filtered).",
    )
    p.add_argument("--skip-bigwig", action="store_true", help="Emit BAMs only.")
    p.add_argument(
        "--atac-bam-type",
        default="rmdup",
        choices=["rmdup", "mapped"],
        help="ATAC BAM source: post-dedup (default) or pre-dedup mapped.",
    )
    p.add_argument("--bigwig-bin-size", type=int, default=10, help="BigWig bin size (bp).")
    p.add_argument("--threads", type=int, default=4, help="Threads for samtools/bamCoverage.")
    p.add_argument("--reference-fasta", default="", help="Reference FASTA override.")
    p.add_argument("--effective-genome-size", type=int, default=0)
    p.add_argument(
        "--species-model",
        default="",
        choices=["", "human", "mouse", "hybrid"],
    )
    p.add_argument(
        "--human-genome-build",
        default="",
        choices=["", "hg19", "hg38"],
    )
    p.add_argument(
        "--nextflow-config",
        default="",
        help="Path to nextflow.config (default: <project-dir>/nextflow.config).",
    )
    return p.parse_args()


def normalize_sequence(seq: str) -> str:
    return re.sub(r"[^ACGTN]", "", (seq or "").strip().upper())


def sequence_matches_with_n(pattern: str, candidate: str) -> bool:
    """True when candidate matches pattern, treating N in pattern as any ACGT base."""
    pattern = normalize_sequence(pattern)
    candidate = normalize_sequence(candidate)
    if not pattern or not candidate or len(pattern) != len(candidate):
        return False
    for p, c in zip(pattern, candidate):
        if p == "N":
            if c not in "ACGT":
                return False
            continue
        if p != c:
            return False
    return True


def find_matrix_column_for_guide(
    pattern: str,
    headers: List[str],
) -> Optional[Tuple[int, str]]:
    """
    Map a library guide sequence to a count-matrix column index.

    Tries exact header match first, then N-aware matching (library N -> any base).
    """
    pattern_norm = normalize_sequence(pattern)
    if not pattern_norm:
        return None

    for idx, header in enumerate(headers):
        if idx == 0:
            continue
        candidate = normalize_sequence(header)
        if candidate == pattern_norm:
            return idx, candidate

    matches: List[Tuple[int, str]] = []
    for idx, header in enumerate(headers):
        if idx == 0 or not header:
            continue
        candidate = normalize_sequence(header)
        if sequence_matches_with_n(pattern_norm, candidate):
            matches.append((idx, candidate))

    if not matches:
        return None
    if len(matches) > 1:
        fixed_bases = sum(1 for base in pattern_norm if base != "N")
        matches.sort(
            key=lambda item: (
                -sum(
                    1
                    for p, c in zip(pattern_norm, item[1])
                    if p != "N" and p == c
                ),
                -fixed_bases,
                item[1],
            )
        )
        print(
            f"WARNING: guide {pattern_norm!r} matched multiple matrix columns; "
            f"using {matches[0][1]!r}",
            file=sys.stderr,
        )
    return matches[0]


def _pick_column(fieldnames: List[str], candidates: Tuple[str, ...]) -> Optional[str]:
    if not fieldnames:
        return None
    lower_map = {name.strip().lower(): name for name in fieldnames if name}
    for cand in candidates:
        if cand in lower_map:
            return lower_map[cand]
    return None


def _row_matches_control_label(row: Dict[str, str], control_label: str) -> Tuple[bool, str]:
    target = control_label.strip().lower()
    for key, value in row.items():
        if value is None or not key:
            continue
        val = str(value).strip().lower()
        if val == target:
            return True, str(value).strip()
    return False, ""


def load_negative_control_guides(
    library_path: str,
    control_label: str = "negative_ctrl",
) -> List[NegativeGuide]:
    if not os.path.isfile(library_path):
        raise FileNotFoundError(f"gRNA library not found: {library_path}")

    with open(library_path, newline="", errors="replace") as fh:
        reader = csv.DictReader(fh)
        if not reader.fieldnames:
            raise ValueError(f"Empty gRNA library: {library_path}")
        fieldnames = list(reader.fieldnames)
        rows = [row for row in reader if any((v or "").strip() for v in row.values())]

    if not rows:
        raise ValueError(f"No guide rows in gRNA library: {library_path}")

    seq_col = _pick_column(fieldnames, SEQUENCE_COL_CANDIDATES)
    name_col = _pick_column(fieldnames, NAME_COL_CANDIDATES)

    if not seq_col:
        raise ValueError(
            f"Could not find a sequence column in {library_path}. "
            f"Expected one of: {', '.join(SEQUENCE_COL_CANDIDATES)}"
        )

    guides: List[NegativeGuide] = []
    seen_sequences: Set[str] = set()
    for row in rows:
        matched, matched_label = _row_matches_control_label(row, control_label)
        if not matched:
            continue
        seq = normalize_sequence(row.get(seq_col, ""))
        if not seq:
            continue
        if seq in seen_sequences:
            continue
        seen_sequences.add(seq)
        name = (row.get(name_col, "") if name_col else "").strip() or matched_label or seq[:16]
        guides.append(NegativeGuide(name=name, sequence=seq, label=matched_label or control_label))

    if not guides:
        raise ValueError(
            f"No guides labeled {control_label!r} in {library_path}. "
            "Check --control-label or library CSV columns."
        )
    return guides


def resolve_grna_library_csv(
    project_dir: str,
    sgrna_sample: str,
    explicit_path: str,
) -> str:
    if explicit_path:
        for base in (
            explicit_path,
            os.path.join(project_dir, explicit_path),
            os.path.join(project_dir, "RAW_FASTQ", explicit_path),
        ):
            if os.path.isfile(base):
                return os.path.abspath(base)
        raise FileNotFoundError(f"gRNA library not found: {explicit_path}")

    manifest = os.path.join(project_dir, "manifests", "sgRNA.tsv")
    if os.path.isfile(manifest):
        with open(manifest, newline="", errors="replace") as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                if (row.get("sample_name") or "").strip() == sgrna_sample:
                    lib = (row.get("grna_library_csv") or "").strip()
                    if lib and os.path.isfile(lib):
                        return os.path.abspath(lib)

    barcode_file = os.path.join(project_dir, "RAW_FASTQ", "input.tsv")
    if os.path.isfile(barcode_file):
        sample_type, _ = load_sample_metadata(barcode_file)
        with open(barcode_file, "r", errors="replace") as fh:
            for raw in fh:
                cols = re.split(r"\t|,", raw.strip())
                if len(cols) >= 5 and cols[0].strip() == sgrna_sample:
                    lib_name = cols[-1].strip()
                    for base in (
                        project_dir,
                        os.path.join(project_dir, "RAW_FASTQ"),
                    ):
                        candidate = os.path.join(base, lib_name)
                        if os.path.isfile(candidate):
                            return os.path.abspath(candidate)

    raise FileNotFoundError(
        f"Could not resolve gRNA library for {sgrna_sample!r}. Pass --grna-library-csv."
    )


def resolve_samples_from_args(args: argparse.Namespace, project_dir: str) -> SampleSet:
    """Thin wrapper matching grna_cell_tracks.resolve_samples without grna-sequence arg."""
    ns = argparse.Namespace(
        sgrna_sample=args.sgrna_sample,
        atac_sample=args.atac_sample,
        rna_sample=args.rna_sample,
        sample_barcode_file=args.sample_barcode_file,
        experimental_group=args.experimental_group,
    )
    from grna_cell_tracks import resolve_samples

    return resolve_samples(ns, project_dir)


def map_guides_to_matrix_columns(
    headers: List[str],
    guides: List[NegativeGuide],
) -> Dict[str, int]:
    mapping: Dict[str, int] = {}
    missing: List[str] = []
    used_columns: Set[int] = set()
    for guide in guides:
        match = find_matrix_column_for_guide(guide.sequence, headers)
        if match is None:
            missing.append(guide.sequence)
            continue
        idx, matrix_seq = match
        if idx in used_columns:
            print(
                f"WARNING: matrix column {matrix_seq!r} already matched; "
                f"skipping duplicate guide {guide.name!r}",
                file=sys.stderr,
            )
            continue
        used_columns.add(idx)
        guide.matrix_sequence = matrix_seq
        mapping[guide.sequence] = idx
    if not mapping:
        raise ValueError(
            "None of the negative-control sequences were found as columns in the count matrix. "
            "Library sequences with N are matched to resolved matrix columns (N matches any ACGT). "
            f"Missing: {', '.join(missing[:5])}"
            + (f" (+{len(missing) - 5} more)" if len(missing) > 5 else "")
        )
    if missing:
        print(
            f"WARNING: {len(missing)} negative-control sequence(s) not in count matrix; skipping.",
            file=sys.stderr,
        )
    print(
        f"Matched {len(mapping)} negative-control guide(s) to count matrix columns",
        file=sys.stderr,
    )
    return mapping


def load_cells_with_negative_controls(
    matrix_path: str,
    col_indices: Dict[str, int],
    min_total_count: int,
    atac_calls: Set[str],
    rna_calls: Set[str],
    require_modality_call: bool,
) -> List[NegativeCtrlCell]:
    if not os.path.isfile(matrix_path):
        raise FileNotFoundError(f"sgRNA count matrix not found: {matrix_path}")

    hits: List[NegativeCtrlCell] = []
    with open(matrix_path, newline="", errors="replace") as fh:
        reader = csv.reader(fh)
        header = next(reader, None)
        if not header:
            raise ValueError(f"Empty count matrix: {matrix_path}")
        if header[0].strip().lower() not in ("cell_barcode", "barcode", "cell"):
            header = ["cell_barcode"] + header

        for row in reader:
            if not row:
                continue
            bc = normalize_barcode(row[0])
            if not bc:
                continue
            guide_hits: Dict[str, int] = {}
            total = 0
            for seq, idx in col_indices.items():
                if idx >= len(row):
                    continue
                try:
                    count = int(float(row[idx] or 0))
                except (TypeError, ValueError):
                    count = 0
                if count > 0:
                    guide_hits[seq] = count
                    total += count
            if total < min_total_count:
                continue
            in_atac = bc in atac_calls
            in_rna = bc in rna_calls
            if require_modality_call and (not in_atac or not in_rna):
                continue
            hits.append(
                NegativeCtrlCell(
                    cell_barcode=bc,
                    total_count=total,
                    guide_hits=guide_hits,
                    in_atac_call=in_atac,
                    in_rna_call=in_rna,
                )
            )

    hits.sort(key=lambda c: (-c.total_count, c.cell_barcode))
    return hits


def write_negative_guides(path: str, guides: List[NegativeGuide]) -> None:
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    with open(path, "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(["guide_name", "library_sequence", "matrix_sequence", "label"])
        for g in guides:
            if not g.matrix_sequence:
                continue
            writer.writerow([g.name, g.sequence, g.matrix_sequence, g.label])


def write_selected_cells(path: str, cells: List[NegativeCtrlCell]) -> None:
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    with open(path, "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(
            [
                "cell_barcode",
                "total_neg_ctrl_count",
                "neg_ctrl_sequences",
                "in_atac_call",
                "in_rna_call",
            ]
        )
        for cell in cells:
            seqs = ";".join(
                f"{seq}:{count}" for seq, count in sorted(cell.guide_hits.items())
            )
            writer.writerow(
                [
                    cell.cell_barcode,
                    cell.total_count,
                    seqs,
                    int(cell.in_atac_call),
                    int(cell.in_rna_call),
                ]
            )


def build_merged_outputs(
    project_dir: str,
    out_dir: str,
    samples: SampleSet,
    barcodes: List[str],
    modality: str,
    atac_bam_type: str,
    star_alignment_mode: str,
    skip_bigwig: bool,
    bigwig_bin_size: int,
    effective_genome_size: int,
    threads: int,
) -> dict:
    merged_dir = os.path.join(out_dir, "merged")
    os.makedirs(merged_dir, exist_ok=True)
    barcodes_path = os.path.join(merged_dir, "merged_barcodes.txt")
    write_barcode_list(barcodes_path, barcodes)

    merged: dict = {
        "n_cells": len(barcodes),
        "barcodes_file": barcodes_path,
        "atac_reads": 0,
        "rna_reads": 0,
        "atac_bam": "",
        "rna_bam": "",
        "atac_bw": "",
        "rna_bw": "",
        "warnings": [],
    }

    if not barcodes:
        merged["warnings"].append("no cells selected for merge")
        return merged

    if modality in ("atac", "both"):
        atac_src = atac_bam_path(project_dir, samples.atac_sample, atac_bam_type)
        ensure_bam_index(atac_src, threads)
        atac_out = os.path.join(merged_dir, "ATAC.negative_ctrl.merged.bam")
        reads, method = extract_merged_bam(atac_src, barcodes, atac_out, threads)
        if method == "stream":
            merged["warnings"].append("ATAC merged: used streaming CB/QNAME filter")
        merged["atac_reads"] = reads
        if reads > 0:
            merged["atac_bam"] = atac_out
            if not skip_bigwig:
                bw_out = os.path.join(merged_dir, "ATAC.negative_ctrl.merged.bw")
                if make_bigwig(
                    atac_out, bw_out, bigwig_bin_size, effective_genome_size, threads
                ):
                    merged["atac_bw"] = bw_out
        else:
            merged["warnings"].append("ATAC merged: zero reads")

    if modality in ("rna", "both"):
        rna_src = rna_bam_path(project_dir, samples.rna_sample, star_alignment_mode)
        ensure_bam_index(rna_src, threads)
        rna_out = os.path.join(merged_dir, "RNA.negative_ctrl.merged.bam")
        reads, method = extract_merged_bam(rna_src, barcodes, rna_out, threads)
        if method == "stream":
            merged["warnings"].append("RNA merged: used streaming CB/QNAME filter")
        merged["rna_reads"] = reads
        if reads > 0:
            merged["rna_bam"] = rna_out
            if not skip_bigwig:
                bw_out = os.path.join(merged_dir, "RNA.negative_ctrl.merged.bw")
                if make_bigwig(
                    rna_out, bw_out, bigwig_bin_size, effective_genome_size, threads
                ):
                    merged["rna_bw"] = bw_out
        else:
            merged["warnings"].append("RNA merged: zero reads")

    return merged


def main() -> int:
    args = parse_args()
    project_dir = os.path.abspath(args.project_dir)
    out_dir = os.path.abspath(args.out_dir)
    os.makedirs(out_dir, exist_ok=True)

    try:
        samples = resolve_samples_from_args(args, project_dir)
        library_path = resolve_grna_library_csv(
            project_dir, samples.sgrna_sample, args.grna_library_csv
        )
        guides = load_negative_control_guides(library_path, args.control_label)
    except (FileNotFoundError, ValueError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1

    matrix_path = grna_count_matrix_path(project_dir, samples.sgrna_sample)
    with open(matrix_path, newline="", errors="replace") as fh:
        header = next(csv.reader(fh), [])
    if header and header[0].strip().lower() not in ("cell_barcode", "barcode", "cell"):
        header = ["cell_barcode"] + header

    try:
        col_indices = map_guides_to_matrix_columns(header, guides)
    except ValueError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1

    starsolo_root = "STARsolo_paired" if args.star_alignment_mode == "paired" else "STARsolo"
    atac_calls, _ = load_atac_barcodes(project_dir, samples.atac_sample)
    rna_calls = load_rna_barcodes(project_dir, samples.rna_sample, starsolo_root)

    cells = load_cells_with_negative_controls(
        matrix_path,
        col_indices,
        args.min_grna_count,
        atac_calls,
        rna_calls,
        args.require_modality_call,
    )

    write_negative_guides(os.path.join(out_dir, "negative_ctrl_guides.tsv"), guides)
    write_selected_cells(os.path.join(out_dir, "selected_cells.tsv"), cells)

    if not cells:
        print("No cells matched negative-control gRNA criteria.", file=sys.stderr)
        manifest = {
            "control_label": args.control_label,
            "grna_library_csv": library_path,
            "n_negative_guides": len(guides),
            "n_cells": 0,
            "merged": {},
        }
        with open(os.path.join(out_dir, "run_manifest.json"), "w") as fh:
            json.dump(manifest, fh, indent=2)
        return 0

    config_path = args.nextflow_config or os.path.join(project_dir, "nextflow.config")
    nf_config = parse_nextflow_config(config_path)
    try:
        reference_fasta, species_key = resolve_reference_fasta(project_dir, args, nf_config)
        effective_genome_size = resolve_effective_genome_size(
            species_key, args.effective_genome_size
        )
    except (FileNotFoundError, ValueError) as exc:
        if args.skip_bigwig:
            reference_fasta = ""
            effective_genome_size = args.effective_genome_size or 2_864_785_220
        else:
            print(f"ERROR: {exc}", file=sys.stderr)
            return 1

    if not args.skip_bigwig and reference_fasta:
        build_chrom_sizes(reference_fasta, out_dir)

    barcodes = [c.cell_barcode for c in cells]
    print(
        f"Merging tracks for {len(cells)} cell(s) with {len(col_indices)} negative-control guide(s)",
        flush=True,
    )
    merged = build_merged_outputs(
        project_dir,
        out_dir,
        samples,
        barcodes,
        args.modality,
        args.atac_bam_type,
        args.star_alignment_mode,
        args.skip_bigwig,
        args.bigwig_bin_size,
        effective_genome_size,
        args.threads,
    )

    manifest = {
        "control_label": args.control_label,
        "grna_library_csv": library_path,
        "experimental_group": samples.experimental_group,
        "sgrna_sample": samples.sgrna_sample,
        "atac_sample": samples.atac_sample,
        "rna_sample": samples.rna_sample,
        "modality": args.modality,
        "require_modality_call": args.require_modality_call,
        "n_negative_guides": len(guides),
        "n_negative_guides_in_matrix": len(col_indices),
        "n_cells": len(cells),
        "reference_fasta": reference_fasta,
        "effective_genome_size": effective_genome_size,
        "merged": merged,
    }
    manifest_path = os.path.join(out_dir, "run_manifest.json")
    with open(manifest_path, "w") as fh:
        json.dump(manifest, fh, indent=2)

    print(f"Wrote {len(cells)} cell(s) to {out_dir}/selected_cells.tsv")
    print(f"Wrote merged tracks under {out_dir}/merged/")
    print(f"Wrote {manifest_path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
