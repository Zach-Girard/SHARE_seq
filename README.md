## SHARE-seq Nextflow Pipeline

### Pipeline overview

End-to-end SHARE-seq processing from raw FASTQs to STARsolo quantification and QC, with RNA/ATAC sample-type routing.

| Step | Process(es) | Output folder |
|------|-------------|---------------|
| 1. Split + parallel demultiplex | `SPLIT_UNDETERMINED_FASTQ`, `DEMULTIPLEX`, `MERGE_DEMUX_CHUNKS` | `RAW_FASTQ/`, `demux/<sample>/` |
| 2. Barcode validation + demux QC | `RENAME_FASTQ`, `FASTQC_DEMUX` | `demux/<sample>/`, `fastqc_demux/<sample>/` |
| 3. RNA/ATAC routing | workflow channel logic from `sample_barcode_file` column 3 | N/A |
| 4. ATAC alignment + filtering/QC (ATAC only) | `BWA_ALIGN_ATAC` | `ATAC/<sample>/` |
| 5. Poly-T filtering (RNA only) | `POLYT_FILTER` | `polyt_filtered/<sample>/` |
| 6. Trimming (optional, RNA only) | `TRIM_R1`, `TRIM_R2_PROTECTED`, `FASTQC_TRIMMED` | `trimmed/<sample>/`, `fastqc_trimmed/` |
| 7. Barcode prepend (RNA only) | `PREPEND_HEADER_BARCODES` | `trimmed/<sample>/` or `polyt_filtered/<sample>/` |
| 8. STAR genome index validation (RNA only) | `STAR_INDEX` | staged `STAR_index_selected/` symlink |
| 9a. Single-end alignment (RNA only) | `STARSOLO_SINGLE` | `STARsolo/<sample>/` |
| 9b. Paired-end alignment (RNA only) | `BUILD_PAIRED_WHITELIST`, `STARSOLO_PAIRED` | `STARsolo_paired/<sample>/` |
| 10. Downstream QC + reports | `KNEE_PLOT`, `BARNYARD_PLOT`, `HYBRID_SPLIT_SPECIES`, `CELL_OVERLAP_BY_GROUP`, `BUILD_QC_HTML` | STARsolo dirs, ATAC dirs, `multiome_overlap/`, `QC_Report*` |

**Key behaviour:**

- **Read structure after demux**: R1 = cDNA, R2 = UMI + PolyT + cDNA. The 24bp cell barcode (3×8bp round barcodes validated by `scripts/rename_fastq.py`) is in the read header.
- **Sample-type routing**: `sample_barcode_file` column 1 is sample ID and column 3 is sample type (`RNA` or `ATAC`). RNA samples continue through Poly-T/STARsolo; ATAC samples go through a dedicated BWA branch.
- **Optional experimental groups**: `sample_barcode_file` can include a 4th column (`Experimental_Group`). The QC Overview can use this to show one KPI card per group with combined group-level estimates (ATAC uses pre-dedup ArchR estimate when available), and `CELL_OVERLAP_BY_GROUP` compares shared 24bp cell barcodes between RNA (STARsolo filtered) and ATAC (ArchR pre-dedup cells from `*.atac_cells.pre_dedup.counts.tsv`) within each group.
- **ATAC branch**: `BWA_INDEX` validates and stages a prebuilt species/build-specific BWA index configured in `nextflow.config`. `BWA_ALIGN_ATAC` runs `bwa mem -C`, keeps mapped reads with `MAPQ >= 30`, adds `CB:Z` tags from QNAME, performs barcode-aware duplicate removal (`samtools markdup --barcode-tag CB -r`), and writes per-sample BAM + QC metrics (`flagstat`, `idxstats`, `stats`) under `ATAC/<sample>/`.
- **Poly-T filtering** is always run. R2 is the anchor (cutadapt matches UMI + PolyT pattern); R1 (cDNA) is synced. Reads are split into `extracted` and `noPolyT` buckets.
- **Trimming** (`trim_reads = true`) runs fastp on Poly-T–extracted R1 (standard) and R2 (first `umi_len` bases protected). Read-dropping filters are disabled to keep R1/R2 in sync. When trimming is off (default), Poly-T–extracted reads flow directly to barcode prepend.
- **Barcode prepend** extracts the 24bp cell barcode from R2's read header and prepends it to R2's sequence. The `withBarcodes_*` output is written to `trimmed/<sample>/` when trimming is enabled, or `polyt_filtered/<sample>/` when disabled.
- **STARsolo alignment** uses R1 (cDNA) + `withBarcodes_R2` (24bp CB + UMI + cDNA). Single-end mode uses `CB_UMI_Complex`; paired-end mode uses `CB_UMI_Simple` with a combinatorial 24bp whitelist. Barcodes are already error-corrected by `scripts/rename_fastq.py`, so no additional barcode matching is needed.
- **STAR genome index** uses a prebuilt index directory configured in `nextflow.config` (validated and staged for each run).
- **Barnyard plots and species splits** are generated only when `species_model = hybrid`.
- **QC report outputs** include `QC_Report.html`, `QC_Report/` (per-sample pages), `QC_Report_assets/`, and `QC_Report_bundle.zip`. Reports now include an ATAC section in the combined page and conditional per-sample ATAC blocks when ATAC outputs exist.

See `main.nf` for channel wiring and `nextflow.config` for all tunable parameters.

### Environment

- The master Conda environment is defined in `environment.yml`.
- Installation and usage instructions (without Homebrew) are in `ENVIRONMENT_SETUP.md`.

When adding new tools or libraries for this project, **always**:

1. Add them to `environment.yml`.
2. Recreate or update your local environment as described in `ENVIRONMENT_SETUP.md`.

### Version control

When you initialize this project as a Git repository and push to GitHub, ensure that:

- `environment.yml` and `ENVIRONMENT_SETUP.md` are committed.
- Any future changes to dependencies are made via `environment.yml` and committed as well.

### Reference and index configuration

This pipeline now uses fixed, prebuilt reference assets configured in `nextflow.config`:

- FASTA and GTF paths per reference set
- BWA index prefix path per reference set
- STAR index directory per reference set

Supported reference sets:

- `hybrid`
- `hg19`
- `hg38`
- `mm10`

Selection logic:

- `species_model=hybrid` -> uses `hybrid_*` paths
- `species_model=mouse` -> uses `mm10_*` paths
- `species_model=human` -> uses `hg19_*` or `hg38_*` based on `human_genome_build` (default: `hg19`)

### Raw FASTQ input (demultiplexing)

- Place undetermined R1/R2 FASTQ files and the sample barcode mapping file in the `RAW_FASTQ/` directory (configurable via `params.raw_fastq`).
- Specify filenames (not full paths) via `--undetermined_r1`, `--undetermined_r2`, and `--sample_barcode_file`; they are resolved relative to `raw_fastq/`.
- `sample_barcode_file` must provide at least 3 columns:
  - column 1 = sample ID
  - column 3 = sample type (`RNA` or `ATAC`)
- Optional column 4 = `Experimental_Group` (free text group label used in the QC report Overview).
- Example `input.tsv` (`sample_barcode_file`) with optional `Experimental_Group`:

```tsv
Sample_Name	Sample_Index	Sample_Type	Experimental_Group
RNA_A	ACGTACGT	RNA	Group_1
ATAC_A	TGCATGCA	ATAC	Group_1
RNA_B	GATTACAA	RNA	Group_2
```
- Undetermined FASTQs are first split into chunks (`split_reads`) and demultiplexed in parallel, then merged per sample.
- `RENAME_FASTQ` validates the three SHARE-seq round barcodes embedded in R1's sequence, rewrites headers with error-corrected barcodes, and writes per-sample outputs to `demux/<sample>/`.
- `FASTQC_DEMUX` writes per-sample reports to `fastqc_demux/<sample>/`.
- ATAC samples then run through `BWA_INDEX`/`BWA_ALIGN_ATAC` and write outputs to `ATAC/<sample>/`:
  - `<sample>.q30.rmdup.sorted.bam`
  - `<sample>.q30.rmdup.sorted.bam.bai`
  - `<sample>.q30.mapped.flagstat.txt`
  - `<sample>.q30.mapped.idxstats.txt`
  - `<sample>.q30.mapped.stats.txt`
  - `<sample>.q30.rmdup.flagstat.txt`
  - `<sample>.q30.rmdup.idxstats.txt`
  - `<sample>.q30.rmdup.stats.txt`
  - `<sample>.cbtag_qc.tsv` (CB-tagging QC: total alignments, tagged alignments, missing-CB alignments, tagged %)

### ATAC branch details (alignment, deduplication, QC, cell estimation)

This section documents the ATAC-specific logic in detail.

#### 1) Alignment and MAPQ filtering

- Input: demultiplexed `matched` R1/R2 FASTQs for samples labeled `ATAC` in `sample_barcode_file` column 3.
- `BWA_ALIGN_ATAC` runs:
  - `bwa mem -C` against the selected genome index (`human`, `mouse`, or `hybrid`)
  - `samtools view -q 30 -F 4` to retain mapped alignments with `MAPQ >= 30`
- Intermediate mapped BAM:
  - `ATAC/<sample>/<sample>.q30.mapped.bam`

#### 2) Barcode tagging before deduplication

- After alignment, the pipeline adds `CB:Z:<24bp_barcode>` tags to alignments before duplicate removal.
- The barcode is extracted from read QNAME (SHARE-seq header convention; fallback regex search for a 24bp `[ACGTN]{24}` token).
- This enables barcode-aware duplicate removal in the next step.
- Tagging QC is written to:
  - `ATAC/<sample>/<sample>.cbtag_qc.tsv`

#### 3) Barcode-aware deduplication strategy

- The pipeline uses:
  - `samtools sort -n` -> `samtools fixmate -m` -> coordinate sort -> `samtools markdup --barcode-tag CB -r`
- Duplicate identity is determined by standard paired-end alignment geometry (position/orientation/mate information) **within each barcode group** (`CB`), instead of collapsing identical coordinates across all cells.
- This is the intended strategy for single-cell ATAC to avoid over-collapsing fragments from different cells that share coordinates.

#### 4) Pre/post dedup ATAC QC metrics

- Pre-dedup metrics (on `.q30.mapped`):
  - `flagstat`, `idxstats`, `stats`
- Post-dedup metrics (on `.q30.rmdup`):
  - `flagstat`, `idxstats`, `stats`
- Combined and per-sample HTML report sections use these files to show:
  - reads pre/post dedup
  - mitochondrial fraction pre/post
  - retained reads %
  - insert-size summaries (via ATAC MultiQC where available)

#### 5) ArchR-based ATAC cell estimation

- `ESTIMATE_ATAC_CELLS` runs ArchR on both:
  - pre-dedup BAM: `ATAC/<sample>/<sample>.q30.mapped.sorted.bam`
  - post-dedup BAM: `ATAC/<sample>/<sample>.q30.rmdup.sorted.bam`
- ArchR flow:
  - `addArchRGenome(...)`
  - `createArrowFiles(minFrags = atac_min_frags_for_cell, minTSS = atac_min_tss_for_cell, bcTag = "CB")`
  - `ArchRProject(...)`
- Summary/count outputs:
  - `ATAC/<sample>/<sample>.atac_cells.summary.tsv`
  - `ATAC/<sample>/<sample>.atac_cells.counts.tsv`
  - `ATAC/<sample>/<sample>.archr_tagged.stats.tsv` (input alignment counts pre/post dedup)

#### 6) ATAC metrics surfaced in reports

- Combined report (`ATAC Key Summary by Sample`) includes:
  - estimated cells (ArchR)
  - median nFrags
  - median TSS enrichment
  - CB-tagged %
  - pre/post read and mt metrics
- Per-sample ATAC report includes:
  - pre/post dedup comparison table
  - ArchR summary metrics
  - CB-tagging QC rows
  - link to global ATAC MultiQC report
- **Dependencies**: Script dependencies come from `environment.yml` (no local `utils.py` required).

### Barcode configuration

- By default, the pipeline uses `barcodes_RC.txt` as the 8bp barcode list.
- A single configuration drives **all barcode-dependent steps**:
  - `RENAME_FASTQ` barcode validation (error-correction against whitelist).
  - Single-end STARsolo (`CB_UMI_Complex` whitelist).
  - Building the 24bp combinatorial whitelist for paired-end STARsolo.

To use **your own 8bp barcode file**:

1. Place your 8bp barcode list (one barcode per line) in the project directory.
2. In `nextflow.config`, set:

   ```groovy
   barcodes_8bp_file = 'my_barcodes.txt'  // your 8bp barcode list
   barcodes_rc       = true               // let the pipeline reverse complement them (recommended)
   ```

- If `barcodes_rc = true`, the pipeline will:
  - Reverse complement every 8bp barcode in `barcodes_8bp_file`.
  - Use this reverse-complemented list everywhere:
    - For `RENAME_FASTQ` barcode validation.
    - As the whitelist for single-end STARsolo.
    - As the input to build the 24bp paired-end STARsolo whitelist (all 3-barcode combinations).
- If `barcodes_rc = false`, the pipeline will use `barcodes_8bp_file` **as-is** in all of the above places.

### Running the pipeline

The pipeline supports two execution modes via profiles in `nextflow.config`:

| Profile   | Executor | Use case                          | Resources (default) |
|-----------|----------|------------------------------------|---------------------|
| `standard`| LSF      | Default; for HPC clusters         | 8 CPUs, 40 GB       |
| `local`   | local    | Workstation/laptop testing        | 4 CPUs, 16 GB       |
| `lsf`     | LSF      | Explicit LSF profile              | 8 CPUs, 40 GB       |

**Run locally** (e.g. on a Mac or workstation without LSF):

```bash
nextflow run main.nf -profile local
```

**Run on an LSF-based HPC cluster** (default):

```bash
nextflow run main.nf
# or explicitly:
nextflow run main.nf -profile lsf
```

- **LSF defaults** (tune to your site in `nextflow.config`):
  - **Queue**: `standard`
  - **Queue size**: `executor.queueSize = 200`
  - **ATAC overrides**: `BWA_INDEX` (index path validation/staging), `BWA_ALIGN_ATAC` (per-sample alignment/filter/QC)
  - **Other processes**: profile-specific overrides in `nextflow.config`

**One-time environment setup on the HPC:**

1. Load Conda (or Mamba) and create the environment:

   ```bash
   module load miniconda            # or your site's preferred module
   cd /path/to/NextFlow
   conda env create -f environment.yml
   ```

2. Ensure Nextflow is available (e.g. `module load nextflow` or via Conda in `environment.yml`).

**Submitting a run via LSF (bsub):**

Create a submission script, e.g. `run_shareseq.lsf`:

```bash
#!/bin/bash
#BSUB -J shareseq
#BSUB -q standard
#BSUB -n 8
#BSUB -M 51200
#BSUB -R "rusage[mem=51200]"
#BSUB -oo shareseq_%J.out

module load miniconda
conda activate nextflow-master

cd $LS_SUBCWD

nextflow run main.nf -profile lsf \
  --undetermined_r1 Undetermined_S0_R1_001.fastq.gz \
  --undetermined_r2 Undetermined_S0_R2_001.fastq.gz \
  --sample_barcode_file sample_barcodes.fasta
```

Submit with:

```bash
bsub < run_shareseq.lsf
```

**Hybrid / paired-end examples:**

```bash
nextflow run main.nf -profile lsf \
  --species_model hybrid \
  --star_alignment_mode paired \
  --barcodes_8bp_file my_barcodes.txt \
  --barcodes_rc true
```

**Checklist before running:**

- `RAW_FASTQ/` contains the undetermined R1/R2 and sample barcode file, and `--undetermined_r1`, `--undetermined_r2`, `--sample_barcode_file` are set.
- `species_model` is set to `human`, `mouse`, or `hybrid`.
- If `species_model=human`, set `human_genome_build` to `hg19` or `hg38` (default is `hg19`).
- The corresponding `*_fasta`, `*_gtf`, `*_bwa_prefix`, and `*_star_index` paths are valid in `nextflow.config`.
- `barcodes_8bp_file` and `barcodes_rc` are set correctly for your barcode scheme.
- `star_alignment_mode` is `single` or `paired` as needed.
- Set `trim_reads = true` to enable fastp trimming. R1 (cDNA) gets standard fastp; R2's first `umi_len` bases (UMI) are always protected.
- For ATAC samples, ensure `species_model`/`human_genome_build` resolve to the intended configured BWA prefix and STAR index.

