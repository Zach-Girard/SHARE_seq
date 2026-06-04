## SHARE-seq Nextflow Pipeline

Nextflow workflow for SHARE-seq: RNA (STARsolo), ATAC (BWA + ArchR), optional sgRNA (cutadapt demux + gRNA counting), and HTML QC with multiome barcode overlap.

### Pipeline overview

| Step | Process(es) | Main outputs |
|------|-------------|--------------|
| 0. Manifests | `BUILD_SAMPLE_MANIFESTS` | `manifests/` |
| 1. RNA/ATAC split + demux | `SPLIT_UNDETERMINED_FASTQ`, `DEMULTIPLEX`, `MERGE_DEMUX_CHUNKS`, `BUILD_DEMUX_STATS_FROM_MERGED` | `demux/<sample>/`, `SHARE-seq.demultiplex.stats.tsv` |
| 2. RNA/ATAC barcode QC | `RENAME_FASTQ`, `FASTQC_DEMUX` | `demux/<sample>/`, `fastqc_demux/<sample>/` |
| 3. sgRNA branch (optional) | `SGRNA_DEMULTIPLEX_CUTADAPT`, `BUILD_SGRNA_RUN_MANIFEST`, `SGRNA_ANALYSIS` | `sgRNA/demux/`, `sgRNA/<sample>/`, `sgRNA/sgrna_qc_summary.tsv` |
| 4. Sample routing | Column 3 of `sample_barcode_file` | — |
| 5. ATAC (ATAC only) | `BWA_INDEX`, `BWA_ALIGN_ATAC`, `ESTIMATE_ATAC_CELLS`, `MULTIQC_ATAC` | `ATAC/<sample>/`, `multiqc_atac/` |
| 6. Poly-T filter (RNA only) | `POLYT_FILTER` | `polyt_filtered/<sample>/` |
| 7. Trim (optional, RNA) | `TRIM_R1`, `TRIM_R2_PROTECTED`, `FASTQC_TRIMMED` | `trimmed/<sample>/`, `fastqc_trimmed/` |
| 8. Barcode prepend (RNA) | `PREPEND_HEADER_BARCODES` | `trimmed/` or `polyt_filtered/` |
| 9. STAR index (RNA) | `STAR_INDEX` | staged `STAR_index_selected/` |
| 10a. STARsolo single (RNA) | `STARSOLO_SINGLE` | `STARsolo/<sample>/` |
| 10b. STARsolo paired (RNA) | `BUILD_PAIRED_WHITELIST`, `STARSOLO_PAIRED` | `STARsolo_paired/<sample>/` |
| 11. QC + reports | `KNEE_PLOT`, `BARNYARD_PLOT`, `HYBRID_SPLIT_SPECIES`, `CELL_OVERLAP_BY_GROUP`, `BUILD_QC_HTML` | `multiome_overlap/`, `QC_Report*` |

Step 3 runs in parallel with steps 1–2 when the barcode file includes `sgRNA` rows. sgRNA-only runs skip steps 1–2 and 4–10.

**Key behaviour**

- **Read structure (RNA/ATAC after demux):** R1 = cDNA; R2 = UMI + PolyT + cDNA. `RENAME_FASTQ` (`scripts/rename_fastq.py`) error-corrects three 8bp round barcodes in R1 and writes the 24bp cell barcode into read headers.
- **Routing:** Column 1 = sample ID; column 3 = `RNA`, `ATAC`, or `sgRNA`. RNA → Poly-T → optional trim → barcode prepend → STARsolo. ATAC → BWA branch after demux QC. sgRNA → separate path ([below](#sgrna-branch)).
- **Experimental groups:** Optional column 4 (`Experimental_Group`) drives QC group cards and `CELL_OVERLAP_BY_GROUP` (shared 24bp barcodes across RNA STARsolo cells, ATAC ArchR pre-dedup `*.atac_cells.pre_dedup.counts.tsv`, and sgRNA `final_*.gRNA.count.csv`, including `RNA_ATAC_sgRNA_Shared`).
- **Poly-T / trim / STARsolo:** Poly-T always runs on RNA (R2 anchor). Trimming is off by default (`trim_reads = false`). STARsolo uses R1 + barcode-prepended R2; barcodes are pre-corrected at rename.
- **Hybrid extras:** Barnyard and species-split plots only when `species_model = hybrid`.
- **QC:** `QC_Report.html`, per-sample pages, assets, and bundle zip; ATAC/sgRNA/multiome sections appear when outputs exist.

See `main.nf` and `nextflow.config` for wiring and parameters.

### Environment

- Conda env: `environment.yml`; setup: `ENVIRONMENT_SETUP.md`.
- Add new dependencies to `environment.yml` and update the env before running.

### Version control

Commit `environment.yml` and `ENVIRONMENT_SETUP.md`; keep dependency changes in `environment.yml`.

### Reference and index configuration

Prebuilt paths in `nextflow.config` (FASTA, GTF, BWA prefix, STAR index):

| Set | Used when |
|-----|-----------|
| `hybrid_*` | `species_model = hybrid` |
| `mm10_*` | `species_model = mouse` |
| `hg19_*` / `hg38_*` | `species_model = human` and `human_genome_build` (default `hg19`) |

### Sample barcode file and inputs

Put inputs under `RAW_FASTQ/` (or set `params.raw_fastq`). Pass **filenames only** for FASTQs and the barcode file (resolved under `raw_fastq/`).

**`--sample_barcode_file` is required** for every run.

**RNA/ATAC rows** (≥3 columns):

| Col | Field |
|-----|--------|
| 1 | Sample ID |
| 2 | Sample index (demux) |
| 3 | `RNA` or `ATAC` |
| 4 | Optional `Experimental_Group` |

**sgRNA rows** (5 columns):

| Col | Field |
|-----|--------|
| 1 | Sample ID |
| 2 | Sample index (cutadapt) |
| 3 | `sgRNA` |
| 4 | `Experimental_Group` |
| 5 | gRNA library CSV filename (resolved under project dir or `RAW_FASTQ/`) |

**Undetermined FASTQs**

| Modality | Params | Default / auto-detect |
|----------|--------|------------------------|
| RNA/ATAC | `--undetermined_r1`, `--undetermined_r2` (both or neither) | `*Undetermined*R1*.fastq.gz` in `RAW_FASTQ/`, excluding `*gRNA*` |
| sgRNA | `--sgrna_undetermined_r1` | `sgRNA_Undetermined_S0_R1_001.fastq.gz` or `*gRNA*Undetermined*R1*.fastq.gz` |

sgRNA demux is **R1-only**; outputs go to `sgRNA/demux/<sample>/`, not `demux/`.

**Examples**

```tsv
Sample_Name	Sample_Index	Sample_Type	Experimental_Group
RNA_A	ACGTACGT	RNA	Group_1
ATAC_A	TGCATGCA	ATAC	Group_1
```

```tsv
sgRNA_C1200	gcagagtc	sgRNA	C1200	C1200_gRNA_library.csv
```

**Demux notes**

- `BUILD_SAMPLE_MANIFESTS` writes `demux_barcodes.tsv` (RNA/ATAC) and, if needed, `sgRNA.tsv`, `sgrna_demux_barcodes.tsv`, `sgrna_barcode.fa`.
- RNA/ATAC: chunked demux (`split_reads`, default 4M; `demux_max_forks`, default 8) via `scripts/demultiplex.py`.
- sgRNA: `scripts/sgrna_demux_cutadapt.sh` with `sgrna_cutadapt_jobs` (default 4).
- `RENAME_FASTQ` publishes `{sample}.matched.R1.fastq.gz` / `.R2.fastq.gz` under `demux/<sample>/`.

### ATAC branch

Runs on **matched** demux R1/R2 for `ATAC` samples.

1. **Index** — `BWA_INDEX` stages the configured BWA prefix as `BWA_index_selected/`.
2. **Align + filter** — `bwa mem -C`, then `samtools view -q 30 -F 4`. Intermediate `.q30.mapped.bam` is removed after sorting; retained pre-dedup artifact: `{sample}.q30.mapped.sorted.bam` (+ flagstat/idxstats/stats).
3. **CB tagging** — `CB:Z:<24bp>` from QNAME (with regex fallback); `*.cbtag_qc.tsv`.
4. **Dedup** — `sort -n` → `fixmate -m` → coordinate sort → `samtools markdup --barcode-tag CB -r` → `{sample}.q30.rmdup.sorted.bam` (+ metrics).
5. **Cells** — `ESTIMATE_ATAC_CELLS` (ArchR) on pre- and post-dedup BAMs → `*.atac_cells.summary.tsv`, `*.atac_cells.counts.tsv`, `*.atac_cells.pre_dedup.counts.tsv`, `*.archr_tagged.stats.tsv`.
6. **Reports** — `MULTIQC_ATAC` → `multiqc_atac/ATAC_MultiQC.html`; tables/plots in `BUILD_QC_HTML`.

ArchR thresholds: `atac_min_frags_for_cell`, `atac_min_tss_for_cell` in `nextflow.config`.

### sgRNA branch

Independent of RNA/ATAC demux (`demultiplex.py`). Requires external SHARE-seq scripts under `share_seq_pipeline_dir` (see `nextflow.config`).

| | RNA/ATAC | sgRNA |
|---|----------|-------|
| Demux | `demultiplex.py` (chunked R1+R2) | cutadapt on shared R1 (`sgrna_demux_cutadapt.sh`) |
| Demux dir | `demux/<sample>/` | `sgRNA/demux/<sample>/` |
| Rename | `scripts/rename_fastq.py` | `share_seq_step2_rename_fastq.py` via `sgrna_analyze.py` |
| Downstream | STARsolo or ATAC | `match_gRNA.py` → count matrix |

**Flow**

1. `SGRNA_DEMULTIPLEX_CUTADAPT` — cutadapt `--no-indels`, `--no-trim`, `-g file:barcode.fa`, `-j ${sgrna_cutadapt_jobs}`; per sample `{sample}.R1.fastq.gz`; unassigned `untrimmed_R1.fastq.gz`.
2. `BUILD_SGRNA_RUN_MANIFEST` — `sgRNA/sgRNA_run.tsv` (`fastq` and `fastq_r2` both point at demuxed R1).
3. `SGRNA_ANALYSIS` — stage R1; rename with **same file as `-r1` and `-r2`** (drops spurious `.matched.R2` / `.junk.R2`); `match_gRNA.py` (`match_grna_start`, default 43). Needs `utils.py` on `PYTHONPATH` (Conda: `python-levenshtein`, `seaborn`, etc.).

**Outputs:** `sgRNA/<sample>/final_<sample>.gRNA.count.csv`, `sgRNA/sgrna_qc_summary.tsv`, and related rename/match files.

**QC:** Combined sgRNA table, per-group KPIs, top-20 guides per sample; per-sample **sgRNA** tab; overlap column `RNA_ATAC_sgRNA_Shared` when matrices exist. sgRNA-only projects still get `QC_Report.html`.

Manual cutadapt-only job: `scripts/sgrna_split.lsf`.

### Barcode configuration

Default 8bp list: `barcodes_RC.txt` (`barcodes_rc = false` uses the file as written; `true` reverse-complements before use).

One whitelist drives `RENAME_FASTQ`, single-end STARsolo (`CB_UMI_Complex`), and the paired 24bp whitelist build.

```groovy
barcodes_8bp_file = 'my_barcodes.txt'
barcodes_rc       = true   // recommended if barcodes are in sequencing orientation
```

### Running the pipeline

| Profile | Executor | Typical use | Default resources |
|---------|----------|-------------|-------------------|
| (none) | LSF | Cluster; uses top-level `process {}` in config | 8 CPUs, 40 GB |
| `lsf` | LSF | Same, with per-process overrides in `nextflow.config` | Mostly 8 CPUs, 12 GB (STAR higher) |
| `local` | local | Workstation test | 4 CPUs, 16 GB |

`standard` profile only sets LSF and is equivalent to the implicit default.

```bash
# Cluster (default executor is LSF)
nextflow run main.nf

# Workstation
nextflow run main.nf -profile local

# Explicit LSF tuning
nextflow run main.nf -profile lsf
```

**LSF submit example** (`run_shareseq.lsf`):

```bash
nextflow run main.nf -profile lsf \
  --undetermined_r1 Undetermined_S0_R1_001.fastq.gz \
  --undetermined_r2 Undetermined_S0_R2_001.fastq.gz \
  --sample_barcode_file sample_barcodes.tsv
```

**Hybrid paired-end example:**

```bash
nextflow run main.nf -profile lsf \
  --species_model hybrid \
  --star_alignment_mode paired \
  --barcodes_8bp_file my_barcodes.txt \
  --barcodes_rc true
```

**Pre-run checklist**

- Barcode file in `RAW_FASTQ/` with valid `RNA` / `ATAC` / `sgRNA` rows.
- RNA/ATAC: undetermined R1+R2 (params or auto-detect); reference paths valid for `species_model` / `human_genome_build`.
- RNA: `star_alignment_mode`, optional `trim_reads`, barcode file settings.
- ATAC: BWA/STAR index selection matches sample species.
- sgRNA: shared Undetermined R1, gRNA library CSVs, `share_seq_pipeline_dir` and script paths.
- Conda env created from `environment.yml`.
