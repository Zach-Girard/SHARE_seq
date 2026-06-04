## SHARE-seq Nextflow Pipeline

Nextflow workflow for SHARE-seq: scRNA, scATAC, optional sgRNA, and HTML QC Reporting.

The following README.md is intended to provide details about the setup, configuration, and sample-type specific processes in this pipeline.  

First-time users should put intentional focus on the Environment Setup, Running the Pipeline, and Pre-Run Checklist Sections. 

See `main.nf` and `nextflow.config` for channel wiring and the full list of parameters.

### Pipeline Overview


| Step                       | Process(es)                                                                                      | Main outputs                                                     |
| -------------------------- | ------------------------------------------------------------------------------------------------ | ---------------------------------------------------------------- |
| 0. Manifests               | `BUILD_SAMPLE_MANIFESTS`                                                                         | `manifests/`                                                     |
| 1. RNA/ATAC split + demux  | `SPLIT_UNDETERMINED_FASTQ`, `DEMULTIPLEX`, `MERGE_DEMUX_CHUNKS`, `BUILD_DEMUX_STATS_FROM_MERGED` | `demux/<sample>/`, `SHARE-seq.demultiplex.stats.tsv`             |
| 2. RNA/ATAC barcode QC     | `RENAME_FASTQ`, `FASTQC_DEMUX`                                                                   | `demux/<sample>/`, `fastqc_demux/<sample>/`                      |
| 3. sgRNA branch (optional) | `SGRNA_DEMULTIPLEX_CUTADAPT`, `BUILD_SGRNA_RUN_MANIFEST`, `SGRNA_ANALYSIS`                       | `sgRNA/demux/`, `sgRNA/<sample>/`, `sgRNA/sgrna_qc_summary.tsv`  |
| 4. Sample routing          | Column 3 of `sample_barcode_file`                                                                | —                                                                |
| 5. ATAC (ATAC only)        | `BWA_INDEX`, `BWA_ALIGN_ATAC`, `ESTIMATE_ATAC_CELLS`, `MULTIQC_ATAC`                             | `ATAC/<sample>/`, `multiqc_atac/`                                |
| 6–10. RNA (RNA only)       | `POLYT_FILTER`, optional trim, `PREPEND_HEADER_BARCODES`, `STAR_INDEX`, `STARSOLO_*`             | `polyt_filtered/`, `trimmed/`, `STARsolo/` or `STARsolo_paired/` |
| 11. QC + reports           | `KNEE_PLOT`, `BARNYARD_PLOT`, `HYBRID_SPLIT_SPECIES`, `CELL_OVERLAP_BY_GROUP`, `BUILD_QC_HTML`   | `multiome_overlap/`, `QC_Report*`                                |


### Environment Setup

A conda environment has been assembled for this pipeline's dependencies. Use the files below to create this environment. 

- Conda environment: `environment.yml`; setup instructions: `ENVIRONMENT_SETUP.md`.

### Species reference and index file configuration

Prebuilt reference paths in `nextflow.config` (FASTA, GTF, BWA prefix, STAR index) have been set for this pipeline. Use `species_model` and `human_genome_build` parameters to change the default paths.


| Set                 | Used when                                                         |
| ------------------- | ----------------------------------------------------------------- |
| `hybrid_*`          | `species_model = hybrid`                                          |
| `mm10_*`            | `species_model = mouse`                                           |
| `hg19_*` / `hg38_*` | `species_model = human` and `human_genome_build` (default `hg19`) |


### Sample barcode file and inputs

Place input files under `RAW_FASTQ/` (or set `params.raw_fastq`). Pass **filenames only** for FASTQs and the sample barcode file; paths resolve under `RAW_FASTQ/` directory.

The `--sample_barcode_file` parameter is **required** for every run.

**RNA/ATAC rows** (≥3 columns):


| Col | Field                         |
| --- | ----------------------------- |
| 1   | Sample ID                     |
| 2   | Sample index (demux)          |
| 3   | `RNA` or `ATAC`               |
| 4   | Optional `Experimental_Group` |


**sgRNA rows** (5 columns):


| Col | Field                                                         |
| --- | ------------------------------------------------------------- |
| 1   | Sample ID                                                     |
| 2   | Sample index (cutadapt)                                       |
| 3   | `sgRNA`                                                       |
| 4   | `Experimental_Group`                                          |
| 5   | gRNA library CSV (resolved under project dir or `RAW_FASTQ/`) |


**Undetermined FASTQs**


| Modality | Params                                                     | Default / auto-detect                                                        |
| -------- | ---------------------------------------------------------- | ---------------------------------------------------------------------------- |
| RNA/ATAC | `--undetermined_r1`, `--undetermined_r2` (both or neither) | `*Undetermined*R1*.fastq.gz` in `RAW_FASTQ/`, excluding `*gRNA`*             |
| sgRNA    | `--sgrna_undetermined_r1`                                  | `sgRNA_Undetermined_S0_R1_001.fastq.gz` or `*gRNA*Undetermined*R1*.fastq.gz` |


**Example Sample Barcode File**

```tsv
Sample_Name Sample_Index  Sample_Type Experimental_Group  gRNA library CSV
RNA_A ACGTACGT  RNA Group_1
ATAC_A  TGCATGCA  ATAC  Group_1
sgRNA_A gcagagtc  sgRNA Group_1 Group_1_gRNA_library.csv
```

### Cell Barcode file configuration

Default 8bp cell barcode whitelist: `barcodes_RC.txt`. With `barcodes_rc = false`, barcodes are used as written; with `true`, each 8bp barcode is reverse-complemented before use.

The same whitelist drives `RENAME_FASTQ`, single-end STARsolo (`CB_UMI_Complex`), and the combinatorial 24bp whitelist for paired-end STARsolo.

```groovy
barcodes_8bp_file = 'my_barcodes.txt'
barcodes_rc       = true   // if barcodes are in sequencing orientation
```

### RNA Branch

RNA samples share steps 1–2 with ATAC (chunked demux, `scripts/rename_fastq.py`, FastQC on matched reads), then continue through Poly-T filtering and STARsolo.

**Read structure after demux:** R1 = cDNA; R2 = UMI + PolyT + cDNA. `RENAME_FASTQ` process error-corrects three 8bp cell barcodes in R1 and places the 24bp cell barcode in read headers.

**Processing**

1. **Poly-T filter** (`POLYT_FILTER`) — always run; R2 is the cutadapt anchor; R1 reads are synced. Outputs under `polyt_filtered/<sample>/` (`extracted` and `noPolyT` buckets).
2. **Trim** (optional) — `trim_reads = true` runs fastp on extracted R1 and on R2 with the first `umi_len` bases (UMI) protected (`TRIM_R1`, `TRIM_R2_PROTECTED`, `FASTQC_TRIMMED`). Default is `trim_reads = false`.
3. **Barcode prepend** (`PREPEND_HEADER_BARCODES`) — 24bp cell barcode from the R2 header is prepended to R2 sequence → `withBarcodes_`* in `trimmed/<sample>/` or `polyt_filtered/<sample>/`.
4. **STAR** — `STAR_INDEX` stages the configured index; alignment uses R1 (cDNA) + barcode-prepended R2.

**STARsolo modes** (`star_alignment_mode` in `nextflow.config`):


| Mode               | Process                                     | Output dir                  | Whitelist                                      |
| ------------------ | ------------------------------------------- | --------------------------- | ---------------------------------------------- |
| `single` (default) | `STARSOLO_SINGLE`                           | `STARsolo/<sample>/`        | `CB_UMI_Complex`                               |
| `paired`           | `BUILD_PAIRED_WHITELIST`, `STARSOLO_PAIRED` | `STARsolo_paired/<sample>/` | `CB_UMI_Simple` + 24bp combinatorial whitelist |


Cell Barcodes are already corrected after rename process, so STARsolo does not perform additional barcode mismatch correction.

**Hybrid** If `species_model = hybrid`: `BARNYARD_PLOT` and `HYBRID_SPLIT_SPECIES` processes run after STARsolo.

**Downstream QC:** knee plots, `Barcodes.stats`, `GeneFull Summary.csv`, and per-sample RNA pages in `QC_Report/<sample>/`.

### ATAC Branch

ATAC samples use matched demux R1/R2 from step 2, then the dedicated alignment branch (see [Pipeline Overview](#pipeline-overview) step 5).

1. **Index** — `BWA_INDEX` stages the configured BWA prefix as `BWA_index_selected/`.
2. **Align + filter** — `bwa mem -C`, then `samtools view -q 30 -F 4`. Pre-dedup artifact: `{sample}.q30.mapped.sorted.bam` (+ flagstat/idxstats/stats).
3. **CB tagging** — `CB:Z:<24bp>` from QNAME (regex fallback); `{sample}.cbtag_qc.tsv`.
4. **Dedup** — `sort -n` → `fixmate -m` → coordinate sort → `samtools markdup --barcode-tag CB -r` → `{sample}.q30.rmdup.sorted.bam`.
5. **Cells** — `ESTIMATE_ATAC_CELLS` (ArchR) on pre- and post-dedup BAMs → `*.atac_cells.summary.tsv`, `*.atac_cells.counts.tsv`, `*.atac_cells.pre_dedup.counts.tsv`, `*.archr_tagged.stats.tsv`.
6. **Reports** — `MULTIQC_ATAC` → `multiqc_atac/ATAC_MultiQC.html`; ATAC tables in `BUILD_QC_HTML`.

ArchR thresholds: `atac_min_frags_for_cell`, `atac_min_tss_for_cell` in `nextflow.config`. Multiome overlap uses pre-dedup ArchR cell barcodes (`*.atac_cells.pre_dedup.counts.tsv`).

### sgRNA Branch

Independent of RNA/ATAC demux (`demultiplex.py`). 


|            | RNA/ATAC                         | sgRNA                                                    |
| ---------- | -------------------------------- | -------------------------------------------------------- |
| Demux      | `demultiplex.py` (chunked R1+R2) | cutadapt on shared R1 (`sgrna_demux_cutadapt.sh`)        |
| Demux dir  | `demux/<sample>/`                | `sgRNA/demux/<sample>/`                                  |
| Rename     | `scripts/rename_fastq.py`        | `share_seq_step2_rename_fastq.py` via `sgrna_analyze.py` |
| Downstream | STARsolo or ATAC                 | `match_gRNA.py` → count matrix                           |


**Flow**

1. `SGRNA_DEMULTIPLEX_CUTADAPT` — Utilizes cutadapt `--no-indels`, `--no-trim`, `-g file:barcode.fa`, `-j ${sgrna_cutadapt_jobs}` (default 4); Processes per sample `{sample}.R1.fastq.gz`; unassigned reads == `untrimmed_R1.fastq.gz`.
2. `BUILD_SGRNA_RUN_MANIFEST` — `sgRNA/sgRNA_run.tsv` (`fastq` and `fastq_r2` both point at demuxed R1).
3. `SGRNA_ANALYSIS` — stage R1; rename with the same file as `-r1` and `-r2`; `match_gRNA.py` (`match_grna_start`, default 43).

**Outputs:** `sgRNA/<sample>/final_<sample>.gRNA.count.csv`, `sgRNA/sgrna_qc_summary.tsv`, and related rename/match artifacts.

**QC:** Combined sgRNA summary, per-group KPIs, top-20 guides per sample; `RNA_ATAC_sgRNA_Shared` in `multiome_overlap/` when matrices exist.

### Running the Pipeline

The pipeline uses the LSF executor and per-process resources in `nextflow.config` (queue `standard`; tune CPUs, memory, and queue for your site).

```bash
nextflow run main.nf \
  --undetermined_r1 Undetermined_S0_R1_001.fastq.gz \
  --undetermined_r2 Undetermined_S0_R2_001.fastq.gz \
  --sgrna_undetermined_r1 sgRNA_Undetermined_S0_R1_001.fastq.gz \
  --sample_barcode_file sample_barcodes.tsv
```

**LSF submit example** — create `run_shareseq.lsf` that activates the Conda environment, `cd` to the project root, runs the command above, then `bsub < run_shareseq.lsf`.

**Hybrid paired-end example:**

```bash
nextflow run main.nf \
  --undetermined_r1 Undetermined_S0_R1_001.fastq.gz \
  --undetermined_r2 Undetermined_S0_R2_001.fastq.gz \
  --sgrna_undetermined_r1 sgRNA_Undetermined_S0_R1_001.fastq.gz \
  --sample_barcode_file sample_barcodes.tsv \
  --species_model hybrid \
  --star_alignment_mode paired \
  --barcodes_8bp_file my_barcodes.txt \
  --barcodes_rc true
```

### Pre-run Checklist

- Conda env created from `environment.yml`.
- Sample Barcode file in `RAW_FASTQ/` with valid `RNA`, `ATAC`, and/or `sgRNA` rows.
- RNA/ATAC: undetermined R1+R2 files (params or auto-detect) in `RAW_FASTQ/`; reference paths set for `species_model` / `human_genome_build`.
- RNA: `star_alignment_mode`, optional `trim_reads`, cell barcode file in root and whitelist settings adjusted.
- ATAC: BWA/STAR index selection matches sample species.
- sgRNA: Undetermined R1 and gRNA library CSVs in `RAW_FASTQ/`.

