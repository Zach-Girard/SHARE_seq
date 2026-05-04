## SHARE-seq Nextflow Pipeline

### Pipeline overview

End-to-end SHARE-seq processing from raw FASTQs to STARsolo quantification and QC, with RNA/ATAC sample-type routing.

| Step | Process(es) | Output folder |
|------|-------------|---------------|
| 1. Split + parallel demultiplex | `SPLIT_UNDETERMINED_FASTQ`, `DEMULTIPLEX`, `MERGE_DEMUX_CHUNKS` | `RAW_FASTQ/`, `demux/<sample>/` |
| 2. Barcode validation + demux QC | `RENAME_FASTQ`, `FASTQC_DEMUX` | `demux/<sample>/`, `fastqc_demux/<sample>/` |
| 3. RNA/ATAC routing | workflow channel logic from `sample_barcode_file` column 3 | N/A |
| 4a. ATAC BWA index (ATAC only) | `BWA_INDEX` | `BWA_index_<species>/` |
| 4b. ATAC alignment + filtering/QC (ATAC only) | `BWA_ALIGN_ATAC` | `ATAC/<sample>/` |
| 5. Poly-T filtering (RNA only) | `POLYT_FILTER` | `polyt_filtered/<sample>/` |
| 6. Trimming (optional, RNA only) | `TRIM_R1`, `TRIM_R2_PROTECTED`, `FASTQC_TRIMMED` | `trimmed/<sample>/`, `fastqc_trimmed/` |
| 7. Barcode prepend (RNA only) | `PREPEND_HEADER_BARCODES` | `trimmed/<sample>/` or `polyt_filtered/<sample>/` |
| 8. STAR genome index (RNA only) | `STAR_INDEX` | `STAR_index_<species><len>bp/` |
| 9a. Single-end alignment (RNA only) | `STARSOLO_SINGLE` | `STARsolo/<sample>/` |
| 9b. Paired-end alignment (RNA only) | `BUILD_PAIRED_WHITELIST`, `STARSOLO_PAIRED` | `STARsolo_paired/<sample>/` |
| 10. Downstream QC + reports | `KNEE_PLOT`, `BARNYARD_PLOT`, `HYBRID_SPLIT_SPECIES`, `BUILD_QC_HTML` | STARsolo dirs, ATAC dirs, `QC_Report*` |

**Key behaviour:**

- **Read structure after demux**: R1 = cDNA, R2 = UMI + PolyT + cDNA. The 24bp cell barcode (3×8bp round barcodes validated by `rename_fastq.py`) is in the read header.
- **Sample-type routing**: `sample_barcode_file` column 1 is sample ID and column 3 is sample type (`RNA` or `ATAC`). RNA samples continue through Poly-T/STARsolo; ATAC samples go through a dedicated BWA branch.
- **ATAC branch**: `BWA_INDEX` builds a species-specific BWA index once and reuses it across runs using a genome fingerprint. `BWA_ALIGN_ATAC` runs `bwa mem -C`, keeps mapped reads with `MAPQ >= 30`, adds `CB:Z` tags from QNAME, performs barcode-aware duplicate removal (`samtools markdup --barcode-tag CB -r`), and writes per-sample BAM + QC metrics (`flagstat`, `idxstats`, `stats`) under `ATAC/<sample>/`.
- **Poly-T filtering** is always run. R2 is the anchor (cutadapt matches UMI + PolyT pattern); R1 (cDNA) is synced. Reads are split into `extracted` and `noPolyT` buckets.
- **Trimming** (`trim_reads = true`) runs fastp on Poly-T–extracted R1 (standard) and R2 (first `umi_len` bases protected). Read-dropping filters are disabled to keep R1/R2 in sync. When trimming is off (default), Poly-T–extracted reads flow directly to barcode prepend.
- **Barcode prepend** extracts the 24bp cell barcode from R2's read header and prepends it to R2's sequence. The `withBarcodes_*` output is written to `trimmed/<sample>/` when trimming is enabled, or `polyt_filtered/<sample>/` when disabled.
- **STARsolo alignment** uses R1 (cDNA) + `withBarcodes_R2` (24bp CB + UMI + cDNA). Single-end mode uses `CB_UMI_Complex`; paired-end mode uses `CB_UMI_Simple` with a combinatorial 24bp whitelist. Barcodes are already error-corrected by `rename_fastq.py`, so no additional barcode matching is needed.
- **STAR genome index** is reused across runs when species, read length, genome FASTA, and GTF all match (fingerprint-based).
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

### Genome assemblies

- Genome download and hybrid genome creation are handled by the script `Genomes/prepare_genomes.sh`.
- By default, it downloads the **GRCh38** (human) and **GRCm39** (mouse) primary assemblies from Ensembl and builds a hybrid FASTA.

To use it:

```bash
cd Genomes
chmod +x prepare_genomes.sh   # first time only
./prepare_genomes.sh
```

If you want to download **different genomes**, edit the URL variables at the top of `Genomes/prepare_genomes.sh` (e.g. point to different species or assembly builds) and then rerun the script.

The Nextflow pipeline has a `species_model` parameter (default: `human`) that can be set to:

- `human`
- `mouse`
- `hybrid` (mixed human–mouse model)

If you change `species_model`, make sure you have downloaded and prepared **matching genome and GTF files** under `Genomes/` and `GTF/` (for example, using the provided helper scripts or by customizing them).

### GTF annotations

- GTF annotation files for GRCh38 and GRCm39 are downloaded by the script `GTF/download_gtf.sh`.

To use it:

```bash
cd GTF
chmod +x download_gtf.sh   # first time only
./download_gtf.sh
```

This will place the human and mouse GTFs under `GTF/Homo_sapiens/` and `GTF/Mus_musculus/`.

### Raw FASTQ input (demultiplexing)

- Place undetermined R1/R2 FASTQ files and the sample barcode mapping file in the `RAW_FASTQ/` directory (configurable via `params.raw_fastq`).
- Specify filenames (not full paths) via `--undetermined_r1`, `--undetermined_r2`, and `--sample_barcode_file`; they are resolved relative to `raw_fastq/`.
- `sample_barcode_file` must provide at least 3 columns:
  - column 1 = sample ID
  - column 3 = sample type (`RNA` or `ATAC`)
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

- `ESTIMATE_ATAC_CELLS` uses the post-dedup BAM:
  - `ATAC/<sample>/<sample>.q30.rmdup.sorted.bam`
- It creates an ArchR-ready tagged BAM for Arrow creation:
  - `ATAC/<sample>/<sample>.archr_tagged.sorted.bam`
- ArchR flow:
  - `addArchRGenome(...)`
  - `createArrowFiles(minFrags = atac_min_frags_for_cell, minTSS = atac_min_tss_for_cell, bcTag = "CB")`
  - `ArchRProject(...)`
- Summary/count outputs:
  - `ATAC/<sample>/<sample>.atac_cells.summary.tsv`
  - `ATAC/<sample>/<sample>.atac_cells.counts.tsv`
  - `ATAC/<sample>/<sample>.archr_tagged.stats.tsv`

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
  - **ATAC overrides**: `BWA_INDEX` (high-memory one-time index build), `BWA_ALIGN_ATAC` (per-sample alignment/filter/QC)
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
- `Genomes/` and `GTF/` have been prepared using `Genomes/prepare_genomes.sh` and `GTF/download_gtf.sh`.
- `species_model` matches your genome/GTF setup (`human`, `mouse`, or `hybrid`).
- `barcodes_8bp_file` and `barcodes_rc` are set correctly for your barcode scheme.
- `star_alignment_mode` is `single` or `paired` as needed.
- Set `trim_reads = true` to enable fastp trimming. R1 (cDNA) gets standard fastp; R2's first `umi_len` bases (UMI) are always protected.
- For ATAC samples, ensure `species_model` matches the intended BWA genome (`human`, `mouse`, or `hybrid`) since BWA index reuse is species-specific.

