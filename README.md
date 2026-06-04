## SHARE-seq Nextflow Pipeline

Nextflow workflow for SHARE-seq: scRNA, scATAC, optional sgRNA, and HTML QC Report.

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


Column 3 routes each sample to the [RNA](#rna), [ATAC](#atac), or [sgRNA](#sgrna) branch. Optional column 4 (`Experimental_Group`) feeds QC group cards and multiome overlap (`CELL_OVERLAP_BY_GROUP`).

Final QC: `QC_Report.html`, per-sample pages under `QC_Report/`, assets, and `QC_Report_bundle.zip`.

See `main.nf` and `nextflow.config` for channel wiring and parameters.

### Environment

- Conda environment: `environment.yml`; setup instructions: `ENVIRONMENT_SETUP.md`.

### Reference and index configuration

Prebuilt reference paths in `nextflow.config` (FASTA, GTF, BWA prefix, STAR index):


| Set                 | Used when                                                         |
| ------------------- | ----------------------------------------------------------------- |
| `hybrid_*`          | `species_model = hybrid`                                          |
| `mm10_*`            | `species_model = mouse`                                           |
| `hg19_*` / `hg38_*` | `species_model = human` and `human_genome_build` (default `hg19`) |


### Sample barcode file and inputs

Place inputs under `RAW_FASTQ/` (or set `params.raw_fastq`). Pass **filenames only** for FASTQs and the barcode file; paths resolve under `raw_fastq/`.

`--sample_barcode_file` is **required** for every run.

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


**Examples**

```tsv
Sample_Name	Sample_Index	Sample_Type	Experimental_Group
RNA_A	ACGTACGT	RNA	Group_1
ATAC_A	TGCATGCA	ATAC	Group_1
```

```tsv
sgRNA_C1200	gcagagtc	sgRNA	C1200	C1200_gRNA_library.csv
```



### Barcode configuration

Default 8bp whitelist: `barcodes_RC.txt`. With `barcodes_rc = false`, barcodes are used as written; with `true`, each 8bp barcode is reverse-complemented before use.

The same whitelist drives `RENAME_FASTQ`, single-end STARsolo (`CB_UMI_Complex`), and the combinatorial 24bp whitelist for paired-end STARsolo.

```groovy
barcodes_8bp_file = 'my_barcodes.txt'
barcodes_rc       = true   // if barcodes are in sequencing orientation
```

### RNA

RNA samples share steps 1–2 with ATAC (chunked demux, `scripts/rename_fastq.py`, FastQC on matched reads), then continue through Poly-T filtering and STARsolo.

**Read structure after demux:** R1 = cDNA; R2 = UMI + PolyT + cDNA. `RENAME_FASTQ` error-corrects three 8bp round barcodes in R1 and places the 24bp cell barcode in read headers.

**Processing**

1. **Poly-T filter** (`POLYT_FILTER`) — always run; R2 is the cutadapt anchor; R1 is synced. Outputs under `polyt_filtered/<sample>/` (`extracted` and `noPolyT` buckets).
2. **Trim** (optional) — `trim_reads = true` runs fastp on extracted R1 and on R2 with the first `umi_len` bases (UMI) protected (`TRIM_R1`, `TRIM_R2_PROTECTED`, `FASTQC_TRIMMED`). Default is `trim_reads = false`.
3. **Barcode prepend** (`PREPEND_HEADER_BARCODES`) — 24bp cell barcode from the R2 header is prepended to R2 sequence → `withBarcodes_`* in `trimmed/<sample>/` or `polyt_filtered/<sample>/`.
4. **STAR** — `STAR_INDEX` stages the configured index; alignment uses R1 (cDNA) + barcode-prepended R2.

**STARsolo modes** (`star_alignment_mode` in `nextflow.config`):


| Mode               | Process                                     | Output dir                  | Whitelist                                      |
| ------------------ | ------------------------------------------- | --------------------------- | ---------------------------------------------- |
| `single` (default) | `STARSOLO_SINGLE`                           | `STARsolo/<sample>/`        | `CB_UMI_Complex`                               |
| `paired`           | `BUILD_PAIRED_WHITELIST`, `STARSOLO_PAIRED` | `STARsolo_paired/<sample>/` | `CB_UMI_Simple` + 24bp combinatorial whitelist |


Barcodes are already corrected at rename, so STARsolo does not perform additional barcode matching.

**Hybrid** (`species_model = hybrid`): `BARNYARD_PLOT` and `HYBRID_SPLIT_SPECIES` run after STARsolo.

**Downstream QC:** knee plots, `Barcodes.stats`, `GeneFull Summary.csv`, and per-sample RNA pages in `QC_Report/<sample>/`.

### ATAC

ATAC samples use matched demux R1/R2 from step 2, then the dedicated alignment branch (see [Pipeline Overview](#pipeline-overview) step 5).

1. **Index** — `BWA_INDEX` stages the configured BWA prefix as `BWA_index_selected/`.
2. **Align + filter** — `bwa mem -C`, then `samtools view -q 30 -F 4`. Pre-dedup artifact: `{sample}.q30.mapped.sorted.bam` (+ flagstat/idxstats/stats).
3. **CB tagging** — `CB:Z:<24bp>` from QNAME (regex fallback); `{sample}.cbtag_qc.tsv`.
4. **Dedup** — `sort -n` → `fixmate -m` → coordinate sort → `samtools markdup --barcode-tag CB -r` → `{sample}.q30.rmdup.sorted.bam`.
5. **Cells** — `ESTIMATE_ATAC_CELLS` (ArchR) on pre- and post-dedup BAMs → `*.atac_cells.summary.tsv`, `*.atac_cells.counts.tsv`, `*.atac_cells.pre_dedup.counts.tsv`, `*.archr_tagged.stats.tsv`.
6. **Reports** — `MULTIQC_ATAC` → `multiqc_atac/ATAC_MultiQC.html`; ATAC tables in `BUILD_QC_HTML`.

ArchR thresholds: `atac_min_frags_for_cell`, `atac_min_tss_for_cell` in `nextflow.config`. Multiome overlap uses pre-dedup ArchR cell barcodes (`*.atac_cells.pre_dedup.counts.tsv`).

### sgRNA

Independent of RNA/ATAC demux (`demultiplex.py`). Requires external SHARE-seq scripts under `share_seq_pipeline_dir` (see `nextflow.config`).


|            | RNA/ATAC                         | sgRNA                                                    |
| ---------- | -------------------------------- | -------------------------------------------------------- |
| Demux      | `demultiplex.py` (chunked R1+R2) | cutadapt on shared R1 (`sgrna_demux_cutadapt.sh`)        |
| Demux dir  | `demux/<sample>/`                | `sgRNA/demux/<sample>/`                                  |
| Rename     | `scripts/rename_fastq.py`        | `share_seq_step2_rename_fastq.py` via `sgrna_analyze.py` |
| Downstream | STARsolo or ATAC                 | `match_gRNA.py` → count matrix                           |


**Flow**

1. `SGRNA_DEMULTIPLEX_CUTADAPT` — cutadapt `--no-indels`, `--no-trim`, `-g file:barcode.fa`, `-j ${sgrna_cutadapt_jobs}` (default 4); per sample `{sample}.R1.fastq.gz`; unassigned `untrimmed_R1.fastq.gz`.
2. `BUILD_SGRNA_RUN_MANIFEST` — `sgRNA/sgRNA_run.tsv` (`fastq` and `fastq_r2` both point at demuxed R1).
3. `SGRNA_ANALYSIS` — stage R1; rename with the same file as `-r1` and `-r2`; `match_gRNA.py` (`match_grna_start`, default 43). Requires SHARE-seq `utils.py` on `PYTHONPATH`.

**Outputs:** `sgRNA/<sample>/final_<sample>.gRNA.count.csv`, `sgRNA/sgrna_qc_summary.tsv`, and related rename/match artifacts.

**QC:** Combined sgRNA summary, per-group KPIs, top-20 guides per sample; per-sample **sgRNA** tab only on sgRNA sample pages; `RNA_ATAC_sgRNA_Shared` in `multiome_overlap/` when matrices exist.

Manual cutadapt-only job: `scripts/sgrna_split.lsf`.

### Running the Pipeline

The pipeline uses the LSF executor and per-process resources in `nextflow.config` (queue `standard`; tune CPUs, memory, and queue for your site).

```bash
nextflow run main.nf \
  --undetermined_r1 Undetermined_S0_R1_001.fastq.gz \
  --undetermined_r2 Undetermined_S0_R2_001.fastq.gz \
  --sample_barcode_file sample_barcodes.tsv
```

**LSF submit example** — create `run_shareseq.lsf` that activates the Conda environment, `cd` to the project root, runs the command above, then `bsub < run_shareseq.lsf`.

**Hybrid paired-end example:**

```bash
nextflow run main.nf \
  --species_model hybrid \
  --star_alignment_mode paired \
  --barcodes_8bp_file my_barcodes.txt \
  --barcodes_rc true
```

**Pre-run checklist**

- Barcode file in `RAW_FASTQ/` with valid `RNA`, `ATAC`, and/or `sgRNA` rows.
- RNA/ATAC: undetermined R1+R2 (params or auto-detect); reference paths valid for `species_model` / `human_genome_build`.
- RNA: `star_alignment_mode`, optional `trim_reads`, barcode whitelist settings.
- ATAC: BWA/STAR index selection matches sample species.
- sgRNA: shared Undetermined R1, gRNA library CSVs, `share_seq_pipeline_dir` and script paths.
- Conda env created from `environment.yml`.

