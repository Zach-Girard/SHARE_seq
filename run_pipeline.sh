#!/bin/bash
# BSUB directives for the Nextflow driver job (submits and monitors child jobs).
# This job stays alive for the full pipeline run; child jobs get resources from nextflow.config.
#BSUB -J shareseq
#BSUB -n 1
#BSUB -M 2000
#BSUB -R "rusage[mem=2000MB]"
#BSUB -oo shareseq_%J.out

cd "$LS_SUBCWD"
nextflow run main.nf \
  --undetermined_r1 Undetermined_S0_001_R1.fastq.gz \
  --undetermined_r2 Undetermined_S0_001_R2.fastq.gz \
  --sample_barcode_file input.tsv \
  -resume
