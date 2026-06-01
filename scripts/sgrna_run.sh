#!/usr/bin/env bash
# sgRNA analysis after cutadapt demux (input: sgRNA_run.tsv with demuxed R1 paths):
#   1. rename_fastq.py -> {sample}.matched.R1.fastq.gz
#   2. gRNA counting on matched R1
#
# Environment (set by scripts/sgrna_analyze.py):
#   SGRNA_MANIFEST      - sgRNA_run.tsv (copied as sgRNA.tsv in work dir)
#   SGRNA_OUT_DIR       - output directory for this task
#   PROJECT_DIR         - Nextflow project root
#   SGRNA_BARCODE_LIST  - effective 8bp barcode whitelist (from pipeline)

set -euo pipefail

if [[ -z "${SGRNA_BARCODE_LIST:-}" ]]; then
  echo "ERROR: SGRNA_BARCODE_LIST is not set" >&2
  exit 1
fi

echo "Step 1: rename_fastq.py on demuxed R1 -> {sample}.matched.R1.fastq.gz"
python3 "${PROJECT_DIR}/scripts/sgrna_rename_matched.py" \
  --manifest "${SGRNA_MANIFEST}" \
  --out-dir "${SGRNA_OUT_DIR}" \
  --barcode-list "${SGRNA_BARCODE_LIST}"

echo "Step 2: gRNA assignment and count matrices"
python3 "${PROJECT_DIR}/scripts/sgrna_count_grna.py" \
  --manifest "${SGRNA_MANIFEST}" \
  --out-dir "${SGRNA_OUT_DIR}" \
  --max-mismatches 2
