#!/usr/bin/env python3
"""
Post-pipeline: merged ATAC/RNA BAM + BigWig for cells with negative-control gRNAs.

Reads guide labels from the gRNA library CSV listed in sample metadata
(column 5 for the sgRNA sample in the selected experimental group), finds
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
    resolve_sample_barcode_file,
    resolve_samples,
    rna_bam_path,
    write_barcode_list,
)

# SHARE-seq library default layout (no header):
#   col1 = guide name, col2 = sequence, col3 = label (e.g. negative_ctrl)
DEFAULT_NAME_COL = 0
DEFAULT_SEQUENCE_COL = 1
DEFAULT_LABEL_COL = 2

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
        required=True,
        help="Sample metadata TSV/CSV with gRNA library in column 5 for sgRNA rows.",
    )
    p.add_argument(
        "--experimental-group",
        required=True,
        help="Experimental_Group used to resolve sgRNA/ATAC/RNA sample IDs.",
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
    return p.parse_args()


def resolve_reference_from_nextflow(project_dir: str) -> Tuple[str, int, str]:
    """Reference FASTA and effective genome size from <project-dir>/nextflow.config."""
    config_path = os.path.join(project_dir, "nextflow.config")
    nf_config = parse_nextflow_config(config_path)
    nf_ref_args = argparse.Namespace(
        reference_fasta="",
        species_model="",
        human_genome_build="",
    )
    reference_fasta, species_key = resolve_reference_fasta(
        project_dir, nf_ref_args, nf_config
    )
    effective_genome_size = resolve_effective_genome_size(species_key, 0)
    return reference_fasta, effective_genome_size, config_path


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


def library_sequence_matches_matrix(library_seq: str, matrix_header: str) -> bool:
    """Match library protospacer to count-matrix column (exact, N-aware, or substring)."""
    lib = normalize_sequence(library_seq)
    mat = normalize_sequence(matrix_header)
    if not lib or not mat:
        return False
    if lib == mat or sequence_matches_with_n(lib, mat):
        return True
    # Matrix columns may contain the full construct; library col2 is often the protospacer.
    if len(lib) >= 8 and lib in mat:
        return True
    return False


def find_matrix_column_for_guide(
    pattern: str,
    headers: List[str],
) -> Optional[Tuple[int, str]]:
    """Map a library guide sequence (column 2) to a count-matrix column index."""
    pattern_norm = normalize_sequence(pattern)
    if not pattern_norm:
        return None

    for idx, header in enumerate(headers):
        if idx == 0:
            continue
        if library_sequence_matches_matrix(pattern_norm, header):
            return idx, normalize_sequence(header)

    return None


def _split_library_line(line: str) -> List[str]:
    return [c.strip() for c in re.split(r",|\t", line.strip())]


def _looks_like_header_row(cols: List[str]) -> bool:
    if len(cols) < 3:
        return False
    c0, c2 = cols[0].lower(), cols[2].lower()
    if c0 in ("guide", "guide_name", "name", "id", "grna_id"):
        return True
    if c2 in ("label", "type", "category", "control", "guide_type", "class"):
        return True
    return False


def _looks_like_data_row(cols: List[str], control_label: str) -> bool:
    """True for SHARE-seq rows like name,ACGT...,negative_ctrl."""
    if len(cols) < 3:
        return False
    if cols[2].strip().lower() == control_label.strip().lower():
        return True
    seq = normalize_sequence(cols[1])
    return len(seq) >= 8 and bool(re.fullmatch(r"[ACGTN]+", seq))


def parse_grna_library_rows(
    library_path: str,
) -> Tuple[List[Dict[str, str]], str]:
    """
    Parse gRNA library CSV.

    Default SHARE-seq layout (with or without header row):
      guide_name, sequence, label
    """
    with open(library_path, newline="", errors="replace") as fh:
        raw_lines = [ln.strip() for ln in fh if ln.strip() and not ln.startswith("#")]

    if not raw_lines:
        raise ValueError(f"Empty gRNA library: {library_path}")

    first_cols = _split_library_line(raw_lines[0])
    has_header = _looks_like_header_row(first_cols) and not _looks_like_data_row(
        first_cols, "negative_ctrl"
    )

    rows: List[Dict[str, str]] = []
    if has_header:
        with open(library_path, newline="", errors="replace") as fh:
            reader = csv.DictReader(fh)
            if not reader.fieldnames:
                raise ValueError(f"Empty gRNA library: {library_path}")
            fieldnames = list(reader.fieldnames)
            if len(fieldnames) < 3:
                raise ValueError(
                    f"gRNA library must have at least 3 columns (name, sequence, label). "
                    f"Found headers: {fieldnames}"
                )
            # SHARE-seq layout is always positional: col1=name, col2=sequence, col3=label.
            # Header names are ignored so a column named 'grna' in col1 cannot be mistaken
            # for the sequence column.
            name_col = fieldnames[DEFAULT_NAME_COL]
            seq_col = fieldnames[DEFAULT_SEQUENCE_COL]
            label_col = fieldnames[DEFAULT_LABEL_COL]
            for row in reader:
                if not any((v or "").strip() for v in row.values()):
                    continue
                rows.append(
                    {
                        "guide_name": (row.get(name_col) or "").strip(),
                        "sequence": (row.get(seq_col) or "").strip(),
                        "label": (row.get(label_col) or "").strip(),
                    }
                )
        mode = (
            f"3-column layout (col2=sequence, col3=label); "
            f"headers: [{name_col!r}, {seq_col!r}, {label_col!r}]"
        )
    else:
        start = 1 if _looks_like_header_row(first_cols) else 0
        for line in raw_lines[start:]:
            cols = _split_library_line(line)
            if len(cols) < 3:
                continue
            rows.append(
                {
                    "guide_name": cols[DEFAULT_NAME_COL],
                    "sequence": cols[DEFAULT_SEQUENCE_COL],
                    "label": cols[DEFAULT_LABEL_COL],
                }
            )
        mode = "positional columns: guide_name (col1), sequence (col2), label (col3)"

    if not rows:
        raise ValueError(f"No guide rows in gRNA library: {library_path}")
    print(f"Parsed gRNA library ({mode}): {len(rows)} row(s)", file=sys.stderr)
    return rows, mode


def load_negative_control_guides(
    library_path: str,
    control_label: str = "negative_ctrl",
) -> List[NegativeGuide]:
    if not os.path.isfile(library_path):
        raise FileNotFoundError(f"gRNA library not found: {library_path}")

    rows, _mode = parse_grna_library_rows(library_path)
    target = control_label.strip().lower()

    guides: List[NegativeGuide] = []
    seen_sequences: Set[str] = set()
    for row in rows:
        label_val = (row.get("label") or "").strip()
        if label_val.lower() != target:
            continue
        seq = normalize_sequence(row.get("sequence", ""))
        if not seq:
            continue
        if seq in seen_sequences:
            continue
        seen_sequences.add(seq)
        name = (row.get("guide_name") or "").strip() or seq[:16]
        guides.append(NegativeGuide(name=name, sequence=seq, label=label_val or control_label))

    if not guides:
        raise ValueError(
            f"No guides with label {control_label!r} in column 3 of {library_path}. "
            "Expected format: guide_name,sequence,label"
        )
    print(
        f"Loaded {len(guides)} negative-control guide(s) from library label column",
        file=sys.stderr,
    )
    example = guides[0]
    print(
        f"Example negative control: name={example.name!r}, "
        f"sequence={example.sequence!r}",
        file=sys.stderr,
    )
    return guides


def _resolve_library_path_value(project_dir: str, lib_value: str) -> Optional[str]:
    if not lib_value:
        return None
    if os.path.isabs(lib_value) and os.path.isfile(lib_value):
        return os.path.abspath(lib_value)
    for base in (
        lib_value,
        project_dir,
        os.path.join(project_dir, "RAW_FASTQ"),
    ):
        candidate = os.path.abspath(os.path.join(base, lib_value))
        if os.path.isfile(candidate):
            return candidate
    return None


def resolve_grna_library_path(project_dir: str, library_path: str) -> str:
    """Resolve --grna-library-csv to an existing file."""
    resolved = _resolve_library_path_value(project_dir, library_path)
    if resolved:
        return resolved
    raise FileNotFoundError(f"gRNA library not found: {library_path}")


def resolve_grna_library_from_sample_barcode(
    project_dir: str,
    sample_barcode_file: str,
    sgrna_sample: str,
) -> str:
    """Resolve gRNA library path from sample barcode file column 5."""
    barcode_path = resolve_sample_barcode_file(project_dir, sample_barcode_file)
    if not barcode_path or not os.path.isfile(barcode_path):
        raise FileNotFoundError(f"Sample barcode file not found: {sample_barcode_file}")

    with open(barcode_path, "r", errors="replace") as fh:
        for raw in fh:
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            cols = [c.strip() for c in re.split(r"\t|,", line)]
            if len(cols) < 5:
                continue
            if cols[0].lower() in ("sample", "sample_name"):
                continue
            if cols[0] != sgrna_sample:
                continue
            library_raw = cols[4]
            if not library_raw:
                raise ValueError(
                    f"Column 5 (gRNA library CSV) is empty for sgRNA sample {sgrna_sample!r} "
                    f"in {barcode_path}"
                )
            resolved = resolve_grna_library_path(project_dir, library_raw)
            print(
                f"Resolved gRNA library from sample barcode file: {library_raw} -> {resolved}",
                flush=True,
            )
            return resolved

    raise ValueError(
        f"Could not find sgRNA sample {sgrna_sample!r} with a column-5 gRNA library "
        f"in sample barcode file {barcode_path}"
    )


def resolve_samples_from_group(args: argparse.Namespace, project_dir: str) -> SampleSet:
    """Resolve sample IDs from sample barcode file and experimental group."""
    group_args = argparse.Namespace(
        sgrna_sample="",
        atac_sample="",
        rna_sample="",
        sample_barcode_file=args.sample_barcode_file,
        experimental_group=args.experimental_group,
    )
    return resolve_samples(group_args, project_dir)


def write_library_audit(
    out_dir: str,
    library_path: str,
    guides: List[NegativeGuide],
) -> str:
    """Record which library file was read."""
    audit_path = os.path.join(out_dir, "library_audit.tsv")
    os.makedirs(out_dir, exist_ok=True)
    with open(audit_path, "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(["library_path", library_path])
        writer.writerow([])
        writer.writerow(["guide_name", "library_sequence", "label"])
        for guide in guides:
            writer.writerow([guide.name, guide.sequence, guide.label])
    return audit_path


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
            "None of the negative-control sequences (library column 2) were found in the "
            "count matrix column headers. "
            "Check library_audit.tsv for the file/sequences actually read. "
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
        samples = resolve_samples_from_group(args, project_dir)
        library_path = resolve_grna_library_from_sample_barcode(
            project_dir,
            args.sample_barcode_file,
            samples.sgrna_sample,
        )
        print(f"Using gRNA library: {library_path}", flush=True)
        print(
            "Resolved samples from experimental group "
            f"{samples.experimental_group!r}: "
            f"sgRNA={samples.sgrna_sample}, ATAC={samples.atac_sample}, RNA={samples.rna_sample}",
            flush=True,
        )
        guides = load_negative_control_guides(library_path, args.control_label)
    except (FileNotFoundError, ValueError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1

    audit_path = write_library_audit(out_dir, library_path, guides)
    print(f"Wrote {audit_path}", flush=True)

    matrix_path = grna_count_matrix_path(project_dir, samples.sgrna_sample)
    print(f"Using count matrix: {matrix_path}", flush=True)
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
            "library_audit": audit_path,
            "experimental_group": samples.experimental_group,
            "sgrna_sample": samples.sgrna_sample,
            "atac_sample": samples.atac_sample,
            "rna_sample": samples.rna_sample,
            "n_negative_guides": len(guides),
            "n_cells": 0,
            "merged": {},
        }
        with open(os.path.join(out_dir, "run_manifest.json"), "w") as fh:
            json.dump(manifest, fh, indent=2)
        return 0

    try:
        reference_fasta, effective_genome_size, config_path = resolve_reference_from_nextflow(
            project_dir
        )
    except (FileNotFoundError, ValueError) as exc:
        if args.skip_bigwig:
            reference_fasta = ""
            effective_genome_size = 0
            config_path = os.path.join(project_dir, "nextflow.config")
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
        "library_audit": audit_path,
        "sgrna_sample": samples.sgrna_sample,
        "atac_sample": samples.atac_sample,
        "rna_sample": samples.rna_sample,
        "modality": args.modality,
        "require_modality_call": args.require_modality_call,
        "n_negative_guides": len(guides),
        "n_negative_guides_in_matrix": len(col_indices),
        "n_cells": len(cells),
        "nextflow_config": config_path,
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
