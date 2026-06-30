#!/usr/bin/env python3
"""
Post-pipeline: find cells containing a target gRNA and emit per-cell ATAC/RNA BAM + BigWig.

Requires a completed Nextflow run with sgRNA, ATAC, and RNA branches for the target group.
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import re
import subprocess
import sys
import tempfile
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import dataclass, field
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

EFFECTIVE_GENOME_SIZES = {
    "hg19": 2_864_785_220,
    "hg38": 2_913_022_398,
    "mm10": 2_652_783_500,
    "hybrid": 2_913_022_398 + 2_652_783_500,
}


@dataclass
class SampleSet:
    sgrna_sample: str
    atac_sample: str
    rna_sample: str
    experimental_group: str = ""


@dataclass
class CellHit:
    cell_barcode: str
    grna_count: int
    in_atac_call: bool = False
    in_rna_call: bool = False


@dataclass
class RunConfig:
    project_dir: str
    out_dir: str
    samples: SampleSet
    grna_sequence: str
    modality: str
    star_alignment_mode: str
    min_grna_count: int
    require_modality_call: bool
    max_cells: Optional[int]
    skip_bigwig: bool
    atac_bam_type: str
    bigwig_bin_size: int
    threads: int
    reference_fasta: str
    effective_genome_size: int
    enable_stream_fallback: bool = False
    chrom_sizes_path: str = field(default="")


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Extract per-cell ATAC/RNA BAM and BigWig for cells containing a target gRNA."
    )
    p.add_argument("--project-dir", required=True, help="Nextflow project / launch directory.")
    p.add_argument(
        "--sample-barcode-file",
        default="",
        help="Sample metadata TSV/CSV (required unless all sample IDs are explicit).",
    )
    p.add_argument(
        "--experimental-group",
        default="",
        help="Experimental_Group to resolve sgRNA/ATAC/RNA samples (column 4 of barcode file).",
    )
    p.add_argument("--sgrna-sample", default="", help="Explicit sgRNA sample ID.")
    p.add_argument("--atac-sample", default="", help="Explicit ATAC sample ID.")
    p.add_argument("--rna-sample", default="", help="Explicit RNA sample ID.")
    p.add_argument("--grna-sequence", required=True, help="Full sgRNA sequence (matrix column).")
    p.add_argument(
        "--out-dir",
        default="",
        help="Output directory. If omitted, use --output-name to place outputs under "
        "grna_tracks/<experimental_group>/<output_name>/.",
    )
    p.add_argument(
        "--output-name",
        default="",
        help="Custom run folder name under grna_tracks/<experimental_group>/ (used when --out-dir is omitted).",
    )
    p.add_argument(
        "--star-alignment-mode",
        default="single",
        choices=["single", "paired"],
        help="STARsolo output root: STARsolo or STARsolo_paired.",
    )
    p.add_argument("--min-grna-count", type=int, default=1, help="Minimum gRNA count per cell.")
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
    p.add_argument("--max-cells", type=int, default=0, help="Safety cap on cells processed (0 = no cap).")
    p.add_argument(
        "--skip-merged",
        action="store_true",
        help="Skip merged BAM/BigWig across all selected cells.",
    )
    p.add_argument("--skip-bigwig", action="store_true", help="Emit BAMs only.")
    p.add_argument(
        "--atac-bam-type",
        default="rmdup",
        choices=["rmdup", "mapped"],
        help="ATAC BAM source: post-dedup (default) or pre-dedup mapped.",
    )
    p.add_argument("--bigwig-bin-size", type=int, default=10, help="BigWig bin size (bp).")
    p.add_argument("--threads", type=int, default=4, help="Threads per samtools/bamCoverage job.")
    p.add_argument(
        "--reference-fasta",
        default="",
        help="Reference FASTA for chrom sizes (default: parse nextflow.config by species).",
    )
    p.add_argument(
        "--effective-genome-size",
        type=int,
        default=0,
        help="Effective genome size for RPGC normalization (default: inferred from species).",
    )
    p.add_argument(
        "--species-model",
        default="",
        choices=["", "human", "mouse", "hybrid"],
        help="Species for reference FASTA lookup when not using --reference-fasta.",
    )
    p.add_argument(
        "--human-genome-build",
        default="",
        choices=["", "hg19", "hg38"],
        help="Human build when species-model=human.",
    )
    p.add_argument(
        "--nextflow-config",
        default="",
        help="Path to nextflow.config (default: <project-dir>/nextflow.config).",
    )
    p.add_argument(
        "--jobs",
        type=int,
        default=4,
        help="Parallel cell workers (each uses --threads for subprocesses).",
    )
    p.add_argument(
        "--enable-stream-fallback",
        action="store_true",
        help="Enable full BAM streaming fallback when tag-based filtering finds no reads (slower).",
    )
    return p.parse_args()


def resolve_sample_barcode_file(project_dir: str, path_value: str) -> str:
    if not path_value:
        return ""
    if os.path.isabs(path_value) and os.path.isfile(path_value):
        return path_value
    for base in (
        project_dir,
        os.path.join(project_dir, "RAW_FASTQ"),
    ):
        candidate = os.path.join(base, path_value)
        if os.path.isfile(candidate):
            return os.path.abspath(candidate)
    return os.path.abspath(os.path.join(project_dir, path_value))


def parse_nextflow_config(config_path: str) -> Dict[str, str]:
    values: Dict[str, str] = {}
    if not os.path.isfile(config_path):
        return values
    key_re = re.compile(
        r"^\s*(species_model|human_genome_build|hg19_fasta|hg38_fasta|mm10_fasta|hybrid_fasta)\s*=\s*['\"]([^'\"]+)['\"]"
    )
    with open(config_path, "r", errors="replace") as fh:
        for line in fh:
            m = key_re.match(line)
            if m:
                values[m.group(1)] = m.group(2)
    return values


def resolve_reference_fasta(
    project_dir: str,
    args: argparse.Namespace,
    config: Dict[str, str],
) -> Tuple[str, str]:
    if args.reference_fasta:
        fasta = os.path.abspath(args.reference_fasta)
        if not os.path.isfile(fasta):
            raise FileNotFoundError(f"Reference FASTA not found: {fasta}")
        return fasta, _infer_species_key_from_fasta(fasta, config, args)

    species = (args.species_model or config.get("species_model") or "human").lower()
    if species == "human":
        build = (args.human_genome_build or config.get("human_genome_build") or "hg19").lower()
        key = build
        fasta_key = f"{build}_fasta"
    elif species == "mouse":
        key = "mm10"
        fasta_key = "mm10_fasta"
    elif species == "hybrid":
        key = "hybrid"
        fasta_key = "hybrid_fasta"
    else:
        raise ValueError(f"Unsupported species_model: {species}")

    fasta = config.get(fasta_key, "")
    if not fasta or not os.path.isfile(fasta):
        raise FileNotFoundError(
            f"Reference FASTA not found for {key} ({fasta_key}={fasta!r}). "
            "Pass --reference-fasta explicitly."
        )
    return fasta, key


def _infer_species_key_from_fasta(
    fasta: str,
    config: Dict[str, str],
    args: argparse.Namespace,
) -> str:
    for key in ("hg19_fasta", "hg38_fasta", "mm10_fasta", "hybrid_fasta"):
        if os.path.abspath(config.get(key, "")) == os.path.abspath(fasta):
            return key.replace("_fasta", "")
    if args.species_model == "human":
        return (args.human_genome_build or config.get("human_genome_build") or "hg19").lower()
    if args.species_model == "mouse":
        return "mm10"
    if args.species_model == "hybrid":
        return "hybrid"
    name = os.path.basename(fasta).lower()
    if "hg38" in name:
        return "hg38"
    if "hg19" in name:
        return "hg19"
    if "mm10" in name:
        return "mm10"
    return "hg19"


def resolve_effective_genome_size(species_key: str, override: int) -> int:
    if override > 0:
        return override
    size = EFFECTIVE_GENOME_SIZES.get(species_key)
    if size is None:
        raise ValueError(
            f"No default effective genome size for {species_key!r}; pass --effective-genome-size."
        )
    return size


def build_chrom_sizes(reference_fasta: str, out_dir: str) -> str:
    os.makedirs(out_dir, exist_ok=True)
    chrom_sizes = os.path.join(out_dir, "chrom.sizes")
    fai = reference_fasta + ".fai"
    if not os.path.isfile(fai):
        subprocess.run(
            ["samtools", "faidx", reference_fasta],
            check=True,
        )
    with open(fai, "r", errors="replace") as fh, open(chrom_sizes, "w") as out:
        for line in fh:
            cols = line.split("\t")
            if len(cols) >= 2:
                out.write(f"{cols[0]}\t{cols[1]}\n")
    return chrom_sizes


def resolve_samples(args: argparse.Namespace, project_dir: str) -> SampleSet:
    explicit = {
        "sgRNA": (args.sgrna_sample or "").strip(),
        "ATAC": (args.atac_sample or "").strip(),
        "RNA": (args.rna_sample or "").strip(),
    }
    if all(explicit.values()):
        return SampleSet(
            sgrna_sample=explicit["sgRNA"],
            atac_sample=explicit["ATAC"],
            rna_sample=explicit["RNA"],
            experimental_group=(args.experimental_group or "").strip(),
        )

    barcode_file = resolve_sample_barcode_file(project_dir, args.sample_barcode_file)
    if not barcode_file or not os.path.isfile(barcode_file):
        raise FileNotFoundError(
            "Provide --sample-barcode-file or explicit --sgrna-sample, --atac-sample, --rna-sample."
        )
    group = (args.experimental_group or "").strip()
    if not group:
        raise ValueError("Provide --experimental-group or explicit sample IDs for all modalities.")

    sample_type, sample_group = load_sample_metadata(barcode_file)
    members = [s for s, g in sample_group.items() if g == group]
    if not members:
        raise ValueError(f"No samples found for experimental group {group!r}.")

    by_type: Dict[str, List[str]] = {"sgRNA": [], "ATAC": [], "RNA": []}
    for sample in members:
        stype = sample_type.get(sample)
        if stype in by_type:
            by_type[stype].append(sample)

    resolved = {}
    for stype, samples in by_type.items():
        if len(samples) == 0:
            raise ValueError(f"Group {group!r} has no {stype} sample.")
        if len(samples) > 1:
            raise ValueError(
                f"Group {group!r} has multiple {stype} samples ({', '.join(samples)}); "
                f"pass --{stype.lower()}-sample explicitly."
            )
        resolved[stype] = samples[0]

    for stype, val in explicit.items():
        key = "sgRNA" if stype == "sgRNA" else stype
        if val:
            resolved[key] = val

    return SampleSet(
        sgrna_sample=resolved["sgRNA"],
        atac_sample=resolved["ATAC"],
        rna_sample=resolved["RNA"],
        experimental_group=group,
    )


def _sanitize_output_name(value: str) -> str:
    return re.sub(r"[^A-Za-z0-9._-]+", "_", (value or "").strip()).strip("._-")


def resolve_out_dir(args: argparse.Namespace, project_dir: str, samples: SampleSet) -> str:
    """
    Resolve output directory.

    Priority:
      1) --out-dir (explicit absolute/relative path)
      2) --output-name under grna_tracks/<experimental_group>/<output_name>
    """
    if args.out_dir:
        return os.path.abspath(args.out_dir)

    output_name = _sanitize_output_name(args.output_name)
    if not output_name:
        raise ValueError("Provide --out-dir or --output-name.")
    if not samples.experimental_group:
        raise ValueError(
            "--output-name requires --experimental-group to build group-specific output path."
        )
    return os.path.join(project_dir, "grna_tracks", samples.experimental_group, output_name)


def grna_count_matrix_path(project_dir: str, sgrna_sample: str) -> str:
    return os.path.join(
        project_dir,
        "sgRNA",
        sgrna_sample,
        f"final_{sgrna_sample}.gRNA.count.csv",
    )


def find_grna_column(headers: List[str], grna_sequence: str) -> Tuple[int, str]:
    target = grna_sequence.strip().upper()
    for idx, header in enumerate(headers):
        if header and header.strip().upper() == target:
            return idx, header
    suggestions = [
        h for h in headers[1:]
        if h and (target in h.upper() or h.upper() in target)
    ][:5]
    msg = f"gRNA sequence not found as a column in count matrix: {grna_sequence!r}"
    if suggestions:
        msg += f". Near-miss columns: {', '.join(suggestions)}"
    raise ValueError(msg)


def load_cells_for_grna(
    matrix_path: str,
    grna_sequence: str,
    min_count: int,
    atac_calls: Set[str],
    rna_calls: Set[str],
    require_modality_call: bool,
    max_cells: Optional[int],
) -> List[CellHit]:
    if not os.path.isfile(matrix_path):
        raise FileNotFoundError(f"sgRNA count matrix not found: {matrix_path}")

    by_barcode: Dict[str, CellHit] = {}
    duplicate_rows = 0
    with open(matrix_path, newline="", errors="replace") as fh:
        reader = csv.reader(fh)
        header = next(reader, None)
        if not header:
            raise ValueError(f"Empty count matrix: {matrix_path}")

        if header[0].strip().lower() not in ("cell_barcode", "barcode", "cell"):
            header = ["cell_barcode"] + header

        col_idx, _ = find_grna_column(header, grna_sequence)

        for row in reader:
            if not row or col_idx >= len(row):
                continue
            bc = normalize_barcode(row[0])
            if not bc:
                continue
            try:
                count = int(float(row[col_idx] or 0))
            except (TypeError, ValueError):
                count = 0
            if count < min_count:
                continue

            in_atac = bc in atac_calls
            in_rna = bc in rna_calls
            if require_modality_call and (not in_atac or not in_rna):
                continue

            existing = by_barcode.get(bc)
            if existing is not None:
                # Two matrix rows normalize to the same 24bp barcode; they would
                # collapse into a single cells/<barcode>/ directory. Keep the
                # higher count so per-cell outputs stay consistent with the manifest.
                duplicate_rows += 1
                if count > existing.grna_count:
                    existing.grna_count = count
                    existing.in_atac_call = in_atac
                    existing.in_rna_call = in_rna
                continue

            by_barcode[bc] = CellHit(
                cell_barcode=bc,
                grna_count=count,
                in_atac_call=in_atac,
                in_rna_call=in_rna,
            )

    if duplicate_rows:
        print(
            f"WARNING: collapsed {duplicate_rows} duplicate barcode row(s) into "
            f"{len(by_barcode)} unique cell(s) (kept max count per barcode).",
            file=sys.stderr,
        )

    hits = list(by_barcode.values())
    hits.sort(key=lambda c: (-c.grna_count, c.cell_barcode))
    if max_cells and len(hits) > max_cells:
        print(
            f"WARNING: {len(hits)} cells match; capping at --max-cells={max_cells}",
            file=sys.stderr,
        )
        hits = hits[:max_cells]
    return hits


def atac_bam_path(project_dir: str, sample: str, bam_type: str) -> str:
    suffix = "q30.rmdup.sorted.bam" if bam_type == "rmdup" else "q30.mapped.sorted.bam"
    return os.path.join(project_dir, "ATAC", sample, f"{sample}.{suffix}")


def rna_bam_path(project_dir: str, sample: str, star_alignment_mode: str) -> str:
    root = "STARsolo_paired" if star_alignment_mode == "paired" else "STARsolo"
    return os.path.join(project_dir, root, sample, "Aligned.sortedByCoord.out.bam")


def ensure_bam_index(bam_path: str, threads: int) -> None:
    bai = bam_path + ".bai"
    # Regenerate when missing OR stale: on re-runs the BAM is rewritten in place,
    # so a leftover index from a prior run points at invalid byte offsets and
    # makes downstream tools (e.g. bamCoverage) hit "Invalid BGZF header".
    if os.path.isfile(bai) and os.path.getmtime(bai) >= os.path.getmtime(bam_path):
        return
    if os.path.isfile(bai):
        try:
            os.remove(bai)
        except OSError:
            pass
    subprocess.run(
        ["samtools", "index", "-@", str(threads), bam_path],
        check=True,
    )


_QNAME_BC_RE = re.compile(r"[ACGTNacgtn]{24}")
# Pipeline ATAC BAMs and STARsolo BAMs store the 24bp cell barcode as CB:Z:<sequence>.
CB_TAG_PREFIX = "CB:Z:"


def cb_tag_field(barcode: str) -> str:
    return f"{CB_TAG_PREFIX}{barcode}"


def count_bam_reads(bam_path: str) -> int:
    if not os.path.isfile(bam_path) or os.path.getsize(bam_path) == 0:
        return 0
    result = subprocess.run(
        ["samtools", "view", "-c", bam_path],
        check=True,
        capture_output=True,
        text=True,
    )
    return int(result.stdout.strip() or 0)


def _cb_from_sam_optional_fields(fields: List[str]) -> Optional[str]:
    for field in fields:
        if field.startswith(CB_TAG_PREFIX):
            return normalize_barcode(field[len(CB_TAG_PREFIX) :])
    return None


def _cb_from_qname(qname: str) -> Optional[str]:
    candidate = qname.rsplit("_", 1)[-1] if "_" in qname else qname
    candidate = candidate.strip().split()[0].split("/")[0]
    matches = _QNAME_BC_RE.findall(candidate) or _QNAME_BC_RE.findall(qname)
    if not matches:
        return None
    return normalize_barcode(matches[-1])


def _sam_line_matches_barcode(sam_line: str, barcode: str) -> bool:
    return _sam_line_matches_barcodes(sam_line, {barcode})


def _sam_line_matches_barcodes(sam_line: str, barcodes: Set[str]) -> bool:
    if not barcodes:
        return False
    cols = sam_line.split("\t")
    if len(cols) < 11:
        return False
    cb = _cb_from_sam_optional_fields(cols[11:])
    if cb and cb in barcodes:
        return True
    qname_bc = _cb_from_qname(cols[0])
    if qname_bc and qname_bc in barcodes:
        return True
    qname = cols[0]
    return any(bc in qname for bc in barcodes)


def _write_sorted_bam_from_stream(
    bam_stream: bytes,
    out_bam: str,
    threads: int,
) -> int:
    if not bam_stream:
        open(out_bam, "wb").close()
        return 0
    subprocess.run(
        ["samtools", "sort", "-@", str(threads), "-o", out_bam, "-"],
        input=bam_stream,
        check=True,
        capture_output=True,
    )
    read_count = count_bam_reads(out_bam)
    if read_count > 0:
        subprocess.run(
            ["samtools", "index", "-@", str(threads), out_bam],
            check=True,
        )
    return read_count


def extract_cell_bam_by_tag(
    source_bam: str,
    barcode: str,
    out_bam: str,
    threads: int,
) -> Tuple[int, str]:
    """Filter reads with CB:Z:<barcode> via samtools -D / -d / -e."""
    os.makedirs(os.path.dirname(out_bam) or ".", exist_ok=True)
    attempts: List[Tuple[str, List[str]]] = []

    # samtools documents barcode lists via -D TAG:file (matches CB:Z values in file).
    with tempfile.NamedTemporaryFile(
        mode="w",
        suffix=".barcodes.txt",
        delete=False,
    ) as fh:
        fh.write(barcode + "\n")
        barcode_file = fh.name
    try:
        attempts.append(
            (
                "CB:-D",
                [
                    "samtools",
                    "view",
                    "-@",
                    str(threads),
                    "-b",
                    "-D",
                    f"CB:{barcode_file}",
                    source_bam,
                ],
            )
        )
        attempts.extend(
            [
                (
                    "CB:-d",
                    [
                        "samtools",
                        "view",
                        "-@",
                        str(threads),
                        "-b",
                        "-d",
                        f"CB:{barcode}",
                        source_bam,
                    ],
                ),
                (
                    "CB:-e",
                    [
                        "samtools",
                        "view",
                        "-@",
                        str(threads),
                        "-b",
                        "-e",
                        f'CB == "{barcode}"',
                        source_bam,
                    ],
                ),
            ]
        )
        for label, cmd in attempts:
            view = subprocess.run(cmd, capture_output=True)
            if view.returncode != 0:
                continue
            if view.stdout:
                return _write_sorted_bam_from_stream(view.stdout, out_bam, threads), label
            # Tag-filter command succeeded but no reads matched this barcode.
            # Treat as a valid empty result to avoid expensive full-file scanning.
            open(out_bam, "wb").close()
            return 0, label
    finally:
        try:
            os.remove(barcode_file)
        except OSError:
            pass
    return 0, ""


def extract_cell_bam_by_qname(
    source_bam: str,
    barcode: str,
    out_bam: str,
    threads: int,
) -> Tuple[int, str]:
    """
    Filter reads whose QNAME ends with _<barcode>.

    STARsolo RNA BAMs in this pipeline carry the 24bp cell barcode in the read
    name (e.g. ..._TACCGAGCTCGAAGTGGCCAATGT) rather than a CB:Z tag, so the CB
    tag filters find nothing. This uses a C-level samtools filter expression
    (single pass, no Python decoding).
    """
    os.makedirs(os.path.dirname(out_bam) or ".", exist_ok=True)
    bc = normalize_barcode(barcode) or barcode.strip().upper()
    cmd = [
        "samtools",
        "view",
        "-@",
        str(threads),
        "-b",
        "-e",
        f'qname =~ "_{bc}$"',
        source_bam,
    ]
    view = subprocess.run(cmd, capture_output=True)
    if view.returncode != 0:
        return 0, ""
    if view.stdout:
        return _write_sorted_bam_from_stream(view.stdout, out_bam, threads), "QNAME:-e"
    open(out_bam, "wb").close()
    return 0, "QNAME:-e"


def extract_cell_bam_by_stream(
    source_bam: str,
    barcode: str,
    out_bam: str,
    threads: int,
) -> int:
    """Scan BAM and keep reads whose CB tag or QNAME matches the cell barcode."""
    os.makedirs(os.path.dirname(out_bam) or ".", exist_ok=True)
    view = subprocess.Popen(
        ["samtools", "view", "-@", str(threads), "-h", source_bam],
        stdout=subprocess.PIPE,
    )
    assert view.stdout is not None
    header_lines: List[bytes] = []
    matched: List[bytes] = []
    for raw in view.stdout:
        if raw.startswith(b"@"):
            header_lines.append(raw)
            continue
        line = raw.decode("utf-8", errors="replace")
        if _sam_line_matches_barcode(line, barcode):
            matched.append(raw)
    view.wait()
    if view.returncode != 0:
        raise subprocess.CalledProcessError(view.returncode, "samtools view")

    if not matched:
        open(out_bam, "wb").close()
        return 0

    return _write_sorted_bam_from_stream(b"".join(header_lines + matched), out_bam, threads)


def extract_cell_bam(
    source_bam: str,
    barcode: str,
    out_bam: str,
    threads: int,
    enable_stream_fallback: bool = False,
    use_qname: bool = False,
) -> Tuple[int, str]:
    """
    Return (read_count, method_label).

    When use_qname is True (STARsolo RNA BAMs, which store the barcode in the
    read name rather than a CB:Z tag), skip the CB-tag filters entirely.
    """
    method = ""
    if not use_qname:
        # 1) CB:Z tag filters (ATAC BAMs carry CB tags).
        reads, method = extract_cell_bam_by_tag(source_bam, barcode, out_bam, threads)
        if reads > 0:
            return reads, method
    # 2) QNAME filter (STARsolo RNA BAMs store the barcode in the read name).
    reads, qmethod = extract_cell_bam_by_qname(source_bam, barcode, out_bam, threads)
    if reads > 0:
        return reads, qmethod
    # 3) Optional full Python scan (slow); off by default.
    if not enable_stream_fallback:
        open(out_bam, "wb").close()
        return 0, method or qmethod
    reads = extract_cell_bam_by_stream(source_bam, barcode, out_bam, threads)
    if reads > 0:
        return reads, "stream"
    return 0, ""


def write_barcode_list(path: str, barcodes: List[str]) -> None:
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    with open(path, "w") as fh:
        for bc in sorted(barcodes):
            fh.write(bc + "\n")


def extract_merged_bam_by_tag(
    source_bam: str,
    barcodes: List[str],
    out_bam: str,
    threads: int,
) -> Tuple[int, str]:
    """Filter source BAM to reads from any listed CB:Z barcode."""
    if not barcodes:
        return 0, ""
    os.makedirs(os.path.dirname(out_bam) or ".", exist_ok=True)
    with tempfile.NamedTemporaryFile(
        mode="w",
        suffix=".barcodes.txt",
        delete=False,
    ) as fh:
        for bc in barcodes:
            fh.write(bc + "\n")
        barcode_file = fh.name
    try:
        view = subprocess.run(
            [
                "samtools",
                "view",
                "-@",
                str(threads),
                "-b",
                "-D",
                f"CB:{barcode_file}",
                source_bam,
            ],
            capture_output=True,
        )
        if view.returncode == 0:
            if view.stdout:
                reads = _write_sorted_bam_from_stream(view.stdout, out_bam, threads)
                return reads, "CB:-D"
            open(out_bam, "wb").close()
            return 0, "CB:-D"
    finally:
        try:
            os.remove(barcode_file)
        except OSError:
            pass
    return 0, ""


def extract_merged_bam_by_stream(
    source_bam: str,
    barcodes: Set[str],
    out_bam: str,
    threads: int,
) -> int:
    """Scan BAM and keep reads matching any selected cell barcode."""
    if not barcodes:
        return 0
    os.makedirs(os.path.dirname(out_bam) or ".", exist_ok=True)
    view = subprocess.Popen(
        ["samtools", "view", "-@", str(threads), "-h", source_bam],
        stdout=subprocess.PIPE,
    )
    assert view.stdout is not None
    header_lines: List[bytes] = []
    matched: List[bytes] = []
    for raw in view.stdout:
        if raw.startswith(b"@"):
            header_lines.append(raw)
            continue
        line = raw.decode("utf-8", errors="replace")
        if _sam_line_matches_barcodes(line, barcodes):
            matched.append(raw)
    view.wait()
    if view.returncode != 0:
        raise subprocess.CalledProcessError(view.returncode, "samtools view")
    if not matched:
        open(out_bam, "wb").close()
        return 0
    return _write_sorted_bam_from_stream(b"".join(header_lines + matched), out_bam, threads)


def extract_merged_bam_by_qname(
    source_bam: str,
    barcodes: List[str],
    out_bam: str,
    threads: int,
    chunk_size: int = 400,
) -> Tuple[int, str]:
    """
    Merge reads whose QNAME ends with _<barcode> for any selected barcode.

    For STARsolo RNA BAMs (no CB:Z tag) the 24bp barcode lives in the read name.
    Uses C-level samtools filter expressions, chunked to keep the regex small,
    then concatenates the per-chunk hits. A read's QNAME contains exactly one
    barcode, so chunks are disjoint (no duplicate reads).
    """
    os.makedirs(os.path.dirname(out_bam) or ".", exist_ok=True)
    bcs = sorted({normalize_barcode(b) or b.strip().upper() for b in barcodes if b})
    if not bcs:
        open(out_bam, "wb").close()
        return 0, ""

    parts: List[str] = []
    cat_bam = f"{out_bam}.cat.bam"
    try:
        for i in range(0, len(bcs), chunk_size):
            chunk = bcs[i : i + chunk_size]
            expr = f'qname =~ "_({"|".join(chunk)})$"'
            view = subprocess.run(
                ["samtools", "view", "-@", str(threads), "-b", "-e", expr, source_bam],
                capture_output=True,
            )
            if view.returncode != 0:
                return 0, ""
            if view.stdout:
                part = f"{out_bam}.part{i // chunk_size}.bam"
                with open(part, "wb") as fh:
                    fh.write(view.stdout)
                parts.append(part)

        if not parts:
            open(out_bam, "wb").close()
            return 0, "QNAME:-e"

        if len(parts) == 1:
            source = parts[0]
        else:
            subprocess.run(["samtools", "cat", "-o", cat_bam, *parts], check=True)
            source = cat_bam

        subprocess.run(
            ["samtools", "sort", "-@", str(threads), "-o", out_bam, source],
            check=True,
        )
        reads = count_bam_reads(out_bam)
        if reads > 0:
            ensure_bam_index(out_bam, threads)
        return reads, "QNAME:-e"
    finally:
        for p in (*parts, cat_bam):
            if os.path.isfile(p):
                try:
                    os.remove(p)
                except OSError:
                    pass


def extract_merged_bam(
    source_bam: str,
    barcodes: List[str],
    out_bam: str,
    threads: int,
) -> Tuple[int, str]:
    """Return (read_count, method_label) for all cells combined."""
    barcode_set = set(barcodes)
    reads, method = extract_merged_bam_by_tag(source_bam, barcodes, out_bam, threads)
    if method:
        return reads, method
    reads = extract_merged_bam_by_stream(source_bam, barcode_set, out_bam, threads)
    if reads > 0:
        return reads, "stream"
    return 0, ""


def _build_merged_bam(
    cfg: RunConfig,
    label: str,
    cell_bams: List[str],
    source_bam: str,
    barcodes: List[str],
    out_bam: str,
) -> Tuple[int, List[str]]:
    """
    Build a merged BAM. Prefer merging the already-extracted per-cell BAMs
    (fast, avoids rescanning the full source BAM); fall back to a direct
    source extraction only when no per-cell BAMs are available.
    """
    warnings: List[str] = []
    reads = merge_cell_bams(cell_bams, out_bam, cfg.threads)
    if reads > 0:
        print(f"{label} merged: combined {reads} reads from per-cell BAMs", flush=True)
        return reads, warnings

    reads, method = extract_merged_bam(source_bam, barcodes, out_bam, cfg.threads)
    if method == "stream":
        warnings.append(f"{label} merged: used streaming CB/QNAME filter")
    if reads > 0:
        print(f"{label} merged: extracted {reads} reads from source BAM", flush=True)
    return reads, warnings


def build_merged_tracks(
    cfg: RunConfig,
    cells: List[CellHit],
    manifest_rows: List[dict],
) -> dict:
    """
    Build merged ATAC/RNA BAM and BigWig across all selected cells.

    Uses the same cell list as per-cell extraction (already filtered by
    --require-modality-call when that flag is set). The merged BAM is built by
    combining the per-cell BAMs produced earlier, which is much faster than
    rescanning the full source BAM and lets BigWig failures stay non-fatal.
    """
    merged_dir = os.path.join(cfg.out_dir, "merged")
    os.makedirs(merged_dir, exist_ok=True)
    barcodes = [c.cell_barcode for c in cells]
    barcodes_path = os.path.join(merged_dir, "merged_barcodes.txt")
    write_barcode_list(barcodes_path, barcodes)

    atac_cell_bams = [r.get("atac_bam", "") for r in manifest_rows]
    rna_cell_bams = [r.get("rna_bam", "") for r in manifest_rows]

    merged: dict = {
        "n_cells": len(barcodes),
        "barcodes_file": barcodes_path,
        "require_modality_call": cfg.require_modality_call,
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

    if cfg.modality in ("atac", "both"):
        atac_src = atac_bam_path(cfg.project_dir, cfg.samples.atac_sample, cfg.atac_bam_type)
        atac_out = os.path.join(merged_dir, "ATAC.merged.bam")
        reads, warns = _build_merged_bam(
            cfg, "ATAC", atac_cell_bams, atac_src, barcodes, atac_out
        )
        merged["warnings"].extend(warns)
        merged["atac_reads"] = reads
        if reads > 0:
            merged["atac_bam"] = atac_out
            if not cfg.skip_bigwig:
                bw_out = os.path.join(merged_dir, "ATAC.merged.bw")
                print("ATAC merged: generating BigWig", flush=True)
                ok, err = safe_make_bigwig(
                    atac_out,
                    bw_out,
                    cfg.bigwig_bin_size,
                    cfg.effective_genome_size,
                    cfg.threads,
                )
                if ok:
                    merged["atac_bw"] = bw_out
                elif err:
                    msg = f"ATAC merged BigWig skipped: {err}"
                    merged["warnings"].append(msg)
                    print(f"WARNING: {msg}", file=sys.stderr, flush=True)
        else:
            merged["warnings"].append("ATAC merged: zero reads")

    if cfg.modality in ("rna", "both"):
        rna_src = rna_bam_path(cfg.project_dir, cfg.samples.rna_sample, cfg.star_alignment_mode)
        rna_out = os.path.join(merged_dir, "RNA.merged.bam")
        reads, warns = _build_merged_bam(
            cfg, "RNA", rna_cell_bams, rna_src, barcodes, rna_out
        )
        merged["warnings"].extend(warns)
        merged["rna_reads"] = reads
        if reads > 0:
            merged["rna_bam"] = rna_out
            if not cfg.skip_bigwig:
                bw_out = os.path.join(merged_dir, "RNA.merged.bw")
                print("RNA merged: generating BigWig", flush=True)
                ok, err = safe_make_bigwig(
                    rna_out,
                    bw_out,
                    cfg.bigwig_bin_size,
                    cfg.effective_genome_size,
                    cfg.threads,
                )
                if ok:
                    merged["rna_bw"] = bw_out
                elif err:
                    msg = f"RNA merged BigWig skipped: {err}"
                    merged["warnings"].append(msg)
                    print(f"WARNING: {msg}", file=sys.stderr, flush=True)
        else:
            merged["warnings"].append("RNA merged: zero reads")

    return merged


def make_bigwig(
    bam_path: str,
    out_bw: str,
    bin_size: int,
    effective_genome_size: int,
    threads: int,
) -> bool:
    if count_bam_reads(bam_path) == 0:
        return False
    subprocess.run(
        [
            "bamCoverage",
            "-b",
            bam_path,
            "-o",
            out_bw,
            "--binSize",
            str(bin_size),
            "--normalizeUsing",
            "RPGC",
            "--effectiveGenomeSize",
            str(effective_genome_size),
            "--numberOfProcessors",
            str(threads),
        ],
        check=True,
    )
    return True


def safe_make_bigwig(
    bam_path: str,
    out_bw: str,
    bin_size: int,
    effective_genome_size: int,
    threads: int,
) -> Tuple[bool, str]:
    """Generate a BigWig without raising. Returns (success, error_message)."""
    try:
        ensure_bam_index(bam_path, threads)
        if make_bigwig(bam_path, out_bw, bin_size, effective_genome_size, threads):
            return True, ""
        return False, "zero reads"
    except subprocess.CalledProcessError as exc:
        return False, f"bamCoverage failed (exit {exc.returncode})"
    except OSError as exc:
        return False, f"bamCoverage error: {exc}"


def merge_cell_bams(cell_bams: List[str], out_bam: str, threads: int) -> int:
    """Merge already-extracted per-cell BAMs into one sorted, indexed BAM."""
    existing = [
        b for b in cell_bams if b and os.path.isfile(b) and os.path.getsize(b) > 0
    ]
    os.makedirs(os.path.dirname(out_bam) or ".", exist_ok=True)
    if not existing:
        open(out_bam, "wb").close()
        return 0
    subprocess.run(
        ["samtools", "merge", "-f", "-@", str(threads), out_bam, *existing],
        check=True,
    )
    reads = count_bam_reads(out_bam)
    if reads > 0:
        ensure_bam_index(out_bam, threads)
    return reads


def process_one_cell(
    cell: CellHit,
    cfg_dict: dict,
) -> dict:
    cfg = RunConfig(**cfg_dict)
    cell_dir = os.path.join(cfg.out_dir, "cells", cell.cell_barcode)
    os.makedirs(cell_dir, exist_ok=True)
    result = {
        "cell_barcode": cell.cell_barcode,
        "grna_count": cell.grna_count,
        "atac_reads": 0,
        "rna_reads": 0,
        "atac_bam": "",
        "rna_bam": "",
        "atac_bw": "",
        "rna_bw": "",
        "warnings": [],
    }

    if cfg.modality in ("atac", "both"):
        atac_src = atac_bam_path(cfg.project_dir, cfg.samples.atac_sample, cfg.atac_bam_type)
        atac_out = os.path.join(cell_dir, f"ATAC.{cell.cell_barcode}.bam")
        reads, method = extract_cell_bam(
            atac_src,
            cell.cell_barcode,
            atac_out,
            cfg.threads,
            cfg.enable_stream_fallback,
        )
        if method == "stream":
            result["warnings"].append("ATAC: used streaming CB/QNAME filter")
        result["atac_reads"] = reads
        result["atac_bam"] = atac_out if reads > 0 else ""
        if reads == 0:
            result["warnings"].append("ATAC: zero reads")
        elif not cfg.skip_bigwig:
            bw_out = os.path.join(cell_dir, f"ATAC.{cell.cell_barcode}.bw")
            ok, err = safe_make_bigwig(
                atac_out,
                bw_out,
                cfg.bigwig_bin_size,
                cfg.effective_genome_size,
                cfg.threads,
            )
            if ok:
                result["atac_bw"] = bw_out
            elif err:
                result["warnings"].append(f"ATAC BigWig skipped: {err}")

    if cfg.modality in ("rna", "both"):
        rna_src = rna_bam_path(cfg.project_dir, cfg.samples.rna_sample, cfg.star_alignment_mode)
        rna_out = os.path.join(cell_dir, f"RNA.{cell.cell_barcode}.bam")
        reads, method = extract_cell_bam(
            rna_src,
            cell.cell_barcode,
            rna_out,
            cfg.threads,
            cfg.enable_stream_fallback,
            use_qname=True,
        )
        if method == "stream":
            result["warnings"].append("RNA: used streaming CB/QNAME filter")
        result["rna_reads"] = reads
        result["rna_bam"] = rna_out if reads > 0 else ""
        if reads == 0:
            result["warnings"].append("RNA: zero reads")
        elif not cfg.skip_bigwig:
            bw_out = os.path.join(cell_dir, f"RNA.{cell.cell_barcode}.bw")
            ok, err = safe_make_bigwig(
                rna_out,
                bw_out,
                cfg.bigwig_bin_size,
                cfg.effective_genome_size,
                cfg.threads,
            )
            if ok:
                result["rna_bw"] = bw_out
            elif err:
                result["warnings"].append(f"RNA BigWig skipped: {err}")

    return result


def write_selected_cells(path: str, cells: List[CellHit]) -> None:
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    with open(path, "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(
            ["cell_barcode", "grna_count", "in_atac_call", "in_rna_call"]
        )
        for cell in cells:
            writer.writerow(
                [
                    cell.cell_barcode,
                    cell.grna_count,
                    int(cell.in_atac_call),
                    int(cell.in_rna_call),
                ]
            )


def validate_inputs(cfg: RunConfig) -> None:
    matrix = grna_count_matrix_path(cfg.project_dir, cfg.samples.sgrna_sample)
    if not os.path.isfile(matrix):
        raise FileNotFoundError(f"sgRNA count matrix not found: {matrix}")

    if cfg.modality in ("atac", "both"):
        atac_bam = atac_bam_path(cfg.project_dir, cfg.samples.atac_sample, cfg.atac_bam_type)
        if not os.path.isfile(atac_bam):
            raise FileNotFoundError(f"ATAC BAM not found: {atac_bam}")
        ensure_bam_index(atac_bam, cfg.threads)

    if cfg.modality in ("rna", "both"):
        rna_bam = rna_bam_path(cfg.project_dir, cfg.samples.rna_sample, cfg.star_alignment_mode)
        if not os.path.isfile(rna_bam):
            raise FileNotFoundError(f"RNA BAM not found: {rna_bam}")
        ensure_bam_index(rna_bam, cfg.threads)


def resolve_parallelism(args: argparse.Namespace) -> Tuple[int, int, int]:
    """
    Keep jobs x threads within allocated CPU slots.

    Returns (effective_threads, effective_jobs, allocated_cpus).
    """
    allocated = os.getenv("LSB_DJOB_NUMPROC", "").strip()
    try:
        allocated_cpus = int(allocated) if allocated else (os.cpu_count() or 1)
    except ValueError:
        allocated_cpus = os.cpu_count() or 1
    allocated_cpus = max(1, allocated_cpus)

    threads = max(1, int(args.threads))
    jobs = max(1, int(args.jobs))

    if threads > allocated_cpus:
        threads = allocated_cpus
    max_jobs_for_threads = max(1, allocated_cpus // threads)
    if jobs > max_jobs_for_threads:
        jobs = max_jobs_for_threads
    return threads, jobs, allocated_cpus


def main() -> int:
    args = parse_args()
    project_dir = os.path.abspath(args.project_dir)

    try:
        samples = resolve_samples(args, project_dir)
        out_dir = resolve_out_dir(args, project_dir, samples)
    except (FileNotFoundError, ValueError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1
    os.makedirs(out_dir, exist_ok=True)
    print(f"Using output directory: {out_dir}", flush=True)

    effective_threads, effective_jobs, allocated_cpus = resolve_parallelism(args)
    if effective_threads != args.threads or effective_jobs != args.jobs:
        print(
            "Adjusted parallelism to avoid CPU oversubscription: "
            f"--threads {args.threads} -> {effective_threads}, "
            f"--jobs {args.jobs} -> {effective_jobs} "
            f"(allocated CPUs: {allocated_cpus})",
            file=sys.stderr,
        )

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
            effective_genome_size = args.effective_genome_size or EFFECTIVE_GENOME_SIZES["hg19"]
            species_key = "hg19"
        else:
            print(f"ERROR: {exc}", file=sys.stderr)
            return 1

    chrom_sizes_path = ""
    if not args.skip_bigwig and reference_fasta:
        chrom_sizes_path = build_chrom_sizes(reference_fasta, out_dir)

    starsolo_root = "STARsolo_paired" if args.star_alignment_mode == "paired" else "STARsolo"
    atac_calls, _ = load_atac_barcodes(project_dir, samples.atac_sample)
    rna_calls = load_rna_barcodes(project_dir, samples.rna_sample, starsolo_root)

    matrix_path = grna_count_matrix_path(project_dir, samples.sgrna_sample)
    max_cells = args.max_cells if args.max_cells > 0 else None

    try:
        cells = load_cells_for_grna(
            matrix_path,
            args.grna_sequence,
            args.min_grna_count,
            atac_calls,
            rna_calls,
            args.require_modality_call,
            max_cells,
        )
    except (FileNotFoundError, ValueError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1

    if not cells:
        print("No cells matched the gRNA criteria.", file=sys.stderr)
        write_selected_cells(os.path.join(out_dir, "selected_cells.tsv"), [])
        return 0

    if len(cells) > 50:
        print(
            f"WARNING: processing {len(cells)} cells (each emits up to 4 files). "
            "This may take a while; use --max-cells to cap.",
            file=sys.stderr,
        )

    write_selected_cells(os.path.join(out_dir, "selected_cells.tsv"), cells)

    cfg = RunConfig(
        project_dir=project_dir,
        out_dir=out_dir,
        samples=samples,
        grna_sequence=args.grna_sequence,
        modality=args.modality,
        star_alignment_mode=args.star_alignment_mode,
        min_grna_count=args.min_grna_count,
        require_modality_call=args.require_modality_call,
        max_cells=max_cells,
        skip_bigwig=args.skip_bigwig,
        atac_bam_type=args.atac_bam_type,
        bigwig_bin_size=args.bigwig_bin_size,
        threads=effective_threads,
        reference_fasta=reference_fasta,
        effective_genome_size=effective_genome_size,
        enable_stream_fallback=args.enable_stream_fallback,
        chrom_sizes_path=chrom_sizes_path,
    )

    try:
        validate_inputs(cfg)
    except FileNotFoundError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1

    cfg_dict = {
        "project_dir": cfg.project_dir,
        "out_dir": cfg.out_dir,
        "samples": cfg.samples,
        "grna_sequence": cfg.grna_sequence,
        "modality": cfg.modality,
        "star_alignment_mode": cfg.star_alignment_mode,
        "min_grna_count": cfg.min_grna_count,
        "require_modality_call": cfg.require_modality_call,
        "max_cells": cfg.max_cells,
        "skip_bigwig": cfg.skip_bigwig,
        "atac_bam_type": cfg.atac_bam_type,
        "bigwig_bin_size": cfg.bigwig_bin_size,
        "threads": cfg.threads,
        "reference_fasta": cfg.reference_fasta,
        "effective_genome_size": cfg.effective_genome_size,
        "enable_stream_fallback": cfg.enable_stream_fallback,
        "chrom_sizes_path": cfg.chrom_sizes_path,
    }

    manifest_rows: List[dict] = []
    jobs = effective_jobs
    if jobs == 1:
        for cell in cells:
            manifest_rows.append(process_one_cell(cell, cfg_dict))
    else:
        with ProcessPoolExecutor(max_workers=jobs) as pool:
            futures = {
                pool.submit(process_one_cell, cell, cfg_dict): cell for cell in cells
            }
            for fut in as_completed(futures):
                manifest_rows.append(fut.result())

    manifest_rows.sort(key=lambda r: (-int(r["grna_count"]), r["cell_barcode"]))

    merged_tracks: dict = {}
    if not args.skip_merged:
        print(
            f"Building merged tracks for {len(cells)} cell(s)"
            + (" (require-modality-call)" if args.require_modality_call else ""),
            flush=True,
        )
        merged_tracks = build_merged_tracks(cfg, cells, manifest_rows)

    manifest = {
        "grna_sequence": args.grna_sequence,
        "experimental_group": samples.experimental_group,
        "sgrna_sample": samples.sgrna_sample,
        "atac_sample": samples.atac_sample,
        "rna_sample": samples.rna_sample,
        "modality": args.modality,
        "require_modality_call": args.require_modality_call,
        "threads": effective_threads,
        "jobs": effective_jobs,
        "enable_stream_fallback": args.enable_stream_fallback,
        "cells_processed": len(manifest_rows),
        "reference_fasta": reference_fasta,
        "effective_genome_size": effective_genome_size,
        "merged": merged_tracks,
        "results": manifest_rows,
    }
    manifest_path = os.path.join(out_dir, "run_manifest.json")
    with open(manifest_path, "w") as fh:
        json.dump(manifest, fh, indent=2)

    print(f"Wrote {len(cells)} cell(s) under {out_dir}/cells/")
    if merged_tracks:
        print(f"Wrote merged tracks under {out_dir}/merged/")
    print(f"Wrote {manifest_path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
