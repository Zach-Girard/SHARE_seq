## NextFlow project

This repository contains a Nextflow-based workflow and supporting resources.

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



### Running on an LSF-based HPC cluster

This project is designed to run on LSF clusters using the `lsf` profile in `nextflow.config`.

- **Recommended defaults** (tune to your site):
  - **Queue**: `general`.
  - **Resources per process**: `4` CPUs, `16 GB` RAM.
  - **Work directory**: `LSF_SCRATCH` if defined, otherwise `/scratch/$USER/nextflow_work`.

**One-time environment setup on the HPC:**

1. Load Conda (or Mamba) and create the environment:
```
bash
module load miniconda            # or your site's preferred module
cd /path/to/NextFlow
conda env create -f environment.yml
```

2. Make sure nextflow is available:
Either via a site module (e.g. module load nextflow), or
Installed in your Conda environment (already included in `environment.yml`).
Submitting a run via LSF (bsub):

Create a simple submission script, for example `run_shareseq.lsf`:
```
#!/bin/bash
#BSUB -J shareseq
#BSUB -q general
#BSUB -n 4
#BSUB -M 16000
#BSUB -R "rusage[mem=16000]"
#BSUB -oo shareseq_%J.out

module load miniconda
conda activate nextflow-master

cd $LS_SUBCWD

nextflow run main.nf -profile lsf \
  --species_model human \
  --star_alignment_mode single 
```

Submit with:
```
bsub < run_shareseq.lsf
```

Hybrid / paired-end examples:

```
nextflow run main.nf -profile lsf \
  --species_model hybrid \
  --star_alignment_mode paired \
  --barcodes_8bp_file my_barcodes.txt \
  --barcodes_rc true
```

Checklist before running on the cluster:

`RAW_FASTQ/` contains your raw, undemultiplexed FASTQs.
`Genomes/` and `GTF/` have been prepared using:
`Genomes/prepare_genomes.sh` (or your customized version).
`GTF/download_gtf.sh` (or your customized version).
`barcodes_8bp_file` and `barcodes_rc` are set correctly for your barcode scheme.

