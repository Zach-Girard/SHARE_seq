#!/bin/bash
# BSUB directives for the Nextflow driver job (submits and monitors child jobs).
# This job stays alive for the full pipeline run; child jobs get resources from nextflow.config.
#BSUB -J shareseq
#BSUB -q general
#BSUB -n 1
#BSUB -M 16000
#BSUB -R "rusage[mem=16000]"
#BSUB -oo shareseq_%J.out
#BSUB -eo shareseq_%J.err

cd "$LS_SUBCWD"
nextflow run main.nf -profile lsf --raw_fastq demux
