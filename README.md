## NextFlow project

This repository contains a Nextflow-based SHARE-seq workflow and supporting resources.

### Pipeline overview

This workflow implements an end-to-end SHARE-seq processing pipeline:

- **Input discovery**: finds raw FASTQs under `params.raw_fastq` (default `RAW_FASTQ/` or `demux/`).
- **Demultiplexing** (`DEMULTIPLEX_PLACEHOLDER`): currently a pass-through copy into `demux/`, designed to be replaced by a real SHARE-seq demux script.
- **Raw QC** (`FASTQC_RAW`): runs FastQC on each raw/demuxed FASTQ into `fastqc_raw/`.
- **Poly-T filtering** (`POLYT_FILTER`): splits R1/R2/R3 into `matched` vs `noPolyT` buckets under `polyt_filtered/`.
- **Trimming** (`TRIM_FASTQ`): trims Poly-T–matched R1/R3 only, producing `trimmed/` FASTQs with names like `sample.matched.trimmed.R1.fastq.gz` and `sample.matched.trimmed.R3.fastq.gz`, plus fastp reports.
- **Trimmed QC** (`FASTQC_TRIMMED`): runs FastQC on all trimmed FASTQs into `fastqc_trimmed/`.
- **Barcode prepending** (`ADD_R2_BARCODES_TO_R3`): extracts three barcodes from R2 (configured by `bc_coords`) and UMIs from R3 (`umi_len`), prepends them to R3, and writes `withBarcodes_*` into `trimmed/`.
- **STAR genome index** (`STAR_INDEX`): builds or reuses a STAR genome index per `species_model` and read length, using fingerprinting (species, read length, genome FASTA, GTF) to avoid redundant rebuilds.
- **STARsolo alignment** (`STARSOLO_SINGLE` / `STARSOLO_PAIRED`): runs STARsolo on trimmed R1 + withBarcodes_R3, producing per-sample `STARsolo/<sample>/` (and optionally `STARsolo_paired/<sample>/`) outputs.
- **Downstream QC** (`KNEE_PLOT`, `BARNYARD_PLOT`, `HYBRID_SPLIT_SPECIES`): generates knee plots and, for `species_model = hybrid`, barnyard collision plots and species-purity split summaries from STARsolo GeneFull outputs.

See `main.nf` for exact channel wiring and `nextflow.config` for tunable parameters.

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

### Raw FASTQ input

- Place all raw, undemultiplexed FASTQ files into the `RAW_FASTQ/` directory before running any demultiplexing or downstream workflows.
- This folder is not managed by any script; it is simply the **expected input location** for your sequencing data.

### Barcode configuration

- By default, the pipeline uses `barcodes_RC.txt` as the 8bp barcode list.
- A single configuration drives **all barcode-dependent steps**:
  - Single-end STARsolo.
  - Paired-end barcode matching.
  - Building the 24bp whitelist for paired-end STARsolo.

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
    - As the whitelist for single-end STARsolo.
    - As the whitelist for paired-end barcode matching.
    - As the input to build the 24bp paired-end STARsolo whitelist (all 3-barcode combinations).
- If `barcodes_rc = false`, the pipeline will use `barcodes_8bp_file` **as-is** in all of the above places.



### Running the pipeline

The pipeline supports two execution modes via profiles in `nextflow.config`:

| Profile   | Executor | Use case                          | Resources (default) |
|-----------|----------|------------------------------------|---------------------|
| `standard`| LSF      | Default; for HPC clusters         | 8 CPUs, 50 GB       |
| `local`   | local    | Workstation/laptop testing         | 4 CPUs, 16 GB       |
| `lsf`     | LSF      | Explicit LSF with queue/workDir    | 8 CPUs, 50 GB       |

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
  - **Queue**: `general`
  - **Resources per process**: 8 CPUs, 50 GB RAM
  - **Work directory**: `LSF_SCRATCH` if defined, otherwise `/scratch/$USER/nextflow_work`

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
#BSUB -q general
#BSUB -n 8
#BSUB -M 51200
#BSUB -R "rusage[mem=51200]"
#BSUB -oo shareseq_%J.out

module load miniconda
conda activate nextflow-master

cd $LS_SUBCWD

nextflow run main.nf -profile lsf \
  --species_model human \
  --star_alignment_mode single
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

- `RAW_FASTQ/` (or `params.raw_fastq`) contains your raw, undemultiplexed FASTQs.
- `Genomes/` and `GTF/` have been prepared using `Genomes/prepare_genomes.sh` and `GTF/download_gtf.sh`.
- `barcodes_8bp_file` and `barcodes_rc` are set correctly for your barcode scheme.

