#!/bin/bash
# BSUB directives for the Nextflow driver job (submits and monitors child jobs).
# This job stays alive for the full pipeline run; child jobs get resources from nextflow.config.
#BSUB -J shareseq
#BSUB -n 1
#BSUB -M 2000
#BSUB -R "rusage[mem=2000MB]"
#BSUB -oo shareseq_%J.out
#BSUB -eo shareseq_%J.err

cd "$LS_SUBCWD"
nextflow run main.nf
