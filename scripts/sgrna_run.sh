#!/usr/bin/env bash
# sgRNA analysis:
#   1. run_lsf.py
#   2. rename_fastq.py (share_seq_step2) -> {sample}.matched.R1.fastq.gz
#   3. gRNA counting on matched R1
#
# Environment (set by scripts/sgrna_analyze.py):
#   SGRNA_MANIFEST      - sgRNA.tsv
#   SGRNA_OUT_DIR       - sgRNA/
#   PROJECT_DIR         - Nextflow project root
#   SGRNA_POOL          - LSF pool (default: CROP_gRNA)
#   SGRNA_BARCODE_LIST  - effective 8bp barcode whitelist (from pipeline)

set -euo pipefail

POOL="${SGRNA_POOL:-CROP_gRNA}"
if [[ -z "${SGRNA_BARCODE_LIST:-}" ]]; then
  echo "ERROR: SGRNA_BARCODE_LIST is not set" >&2
  exit 1
fi

echo "Step 1: run_lsf.py -f sgRNA.tsv -p ${POOL}"
run_lsf.py -f "${SGRNA_MANIFEST}" -p "${POOL}"

echo "Step 2: rename_fastq.py on {sample}_C3.fastq.gz -> {sample}.matched.R1.fastq.gz"
python3 "${PROJECT_DIR}/scripts/sgrna_rename_matched.py" \
  --manifest "${SGRNA_MANIFEST}" \
  --out-dir "${SGRNA_OUT_DIR}" \
  --barcode-list "${SGRNA_BARCODE_LIST}"

echo "Step 3: gRNA assignment and count matrices"
python3 "${PROJECT_DIR}/scripts/sgrna_count_grna.py" \
  --manifest "${SGRNA_MANIFEST}" \
  --out-dir "${SGRNA_OUT_DIR}" \
  --max-mismatches 2
