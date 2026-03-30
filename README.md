## SHARE-seq Nextflow Pipeline

### Pipeline overview

End-to-end SHARE-seq processing from raw FASTQs to STARsolo quantification and QC.

| Step | Process(es) | Output folder |
|------|-------------|---------------|
| 1. Demultiplex + barcode validation | `DEMULTIPLEX`, `RENAME_FASTQ` | `demux/` |
| 2. Demux QC | `FASTQC_DEMUX` | `fastqc_demux/` |
| 3. Poly-T filtering | `POLYT_FILTER` | `polyt_filtered/` |
| 4. Trimming (optional) | `TRIM_R1`, `TRIM_R2_PROTECTED`, `FASTQC_TRIMMED` | `trimmed/`, `fastqc_trimmed/` |
| 5. Barcode prepend | `PREPEND_HEADER_BARCODES` | `trimmed/` or `polyt_filtered/` |
| 6. STAR genome index | `STAR_INDEX` | `STAR_index_<species><len>bp/` |
| 7a. Single-end alignment | `STARSOLO_SINGLE` | `STARsolo/<sample>/` |
| 7b. Paired-end alignment | `BUILD_PAIRED_WHITELIST`, `STARSOLO_PAIRED` | `STARsolo_paired/<sample>/` |
| 8. Downstream QC | `KNEE_PLOT`, `BARNYARD_PLOT`, `HYBRID_SPLIT_SPECIES` | inside alignment output folder |

**Key behaviour:**

- **Read structure after demux**: R1 = cDNA, R2 = UMI + PolyT + cDNA. The 24bp cell barcode (3×8bp round barcodes validated by `rename_fastq.py`) is in the read header.
- **Poly-T filtering** is always run. R2 is the anchor (cutadapt matches UMI + PolyT pattern); R1 (cDNA) is synced. Reads are split into `matched` and `noPolyT` buckets.
- **Trimming** (`trim_reads = true`) runs fastp on Poly-T–matched R1 (standard) and R2 (first `umi_len` bases protected). Read-dropping filters are disabled to keep R1/R2 in sync. When trimming is off (default), Poly-T–matched reads flow directly to barcode prepend.
- **Barcode prepend** extracts the 24bp cell barcode from R2's read header and prepends it to R2's sequence. The `withBarcodes_*` output is written to `trimmed/` when trimming is enabled, or `polyt_filtered/` when disabled.
- **STARsolo alignment** uses R1 (cDNA) + `withBarcodes_R2` (24bp CB + UMI + cDNA). Single-end mode uses `CB_UMI_Complex`; paired-end mode uses `CB_UMI_Simple` with a combinatorial 24bp whitelist. Barcodes are already error-corrected by `rename_fastq.py`, so no additional barcode matching is needed.
- **STAR genome index** is reused across runs when species, read length, genome FASTA, and GTF all match (fingerprint-based).
- **Barnyard plots and species splits** are generated only when `species_model = hybrid`.

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
- `DEMULTIPLEX` splits undetermined reads by sample index barcode (from the read header). Output goes to `demux/`.
- `RENAME_FASTQ` validates the three SHARE-seq round barcodes embedded in R1's sequence, rewrites headers with error-corrected barcodes, and outputs per-sample `<sample>.matched.R1.fastq.gz` / `.R2.fastq.gz` to `demux/`.
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
  - **Resources per process**: 8 CPUs, 40 GB RAM

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

