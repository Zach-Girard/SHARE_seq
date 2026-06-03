#!/usr/bin/env bash
# sgRNA demultiplex: cutadapt with anchored sample indices (same as split.lsf).
set -euo pipefail

usage() {
  echo "Usage: $0 --barcode-table TSV --input R1.fastq.gz --input-r2 R2.fastq.gz --out-dir demux [--error-rate 0.15]" >&2
  exit 1
}

BARCODE_TABLE=""
BARCODE_FA=""
INPUT_R1=""
INPUT_R2=""
OUT_DIR="sgRNA/demux"
ERROR_RATE="0.15"
JOBS=1

while [[ $# -gt 0 ]]; do
  case "$1" in
    --barcode-table) BARCODE_TABLE="$2"; shift 2 ;;
    --barcode-fa) BARCODE_FA="$2"; shift 2 ;;
    --input) INPUT_R1="$2"; shift 2 ;;
    --input-r2) INPUT_R2="$2"; shift 2 ;;
    --out-dir) OUT_DIR="$2"; shift 2 ;;
    --error-rate) ERROR_RATE="$2"; shift 2 ;;
    --jobs) JOBS="$2"; shift 2 ;;
    -h|--help) usage ;;
    *) echo "Unknown option: $1" >&2; usage ;;
  esac
done

[[ -n "${BARCODE_TABLE}" && -n "${INPUT_R1}" && -n "${INPUT_R2}" ]] || usage
[[ -f "${BARCODE_TABLE}" ]] || { echo "ERROR: barcode table not found: ${BARCODE_TABLE}" >&2; exit 1; }
[[ -f "${INPUT_R1}" ]] || { echo "ERROR: input R1 FASTQ not found: ${INPUT_R1}" >&2; exit 1; }
[[ -f "${INPUT_R2}" ]] || { echo "ERROR: input R2 FASTQ not found: ${INPUT_R2}" >&2; exit 1; }

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BUILD_FA="${SCRIPT_DIR}/build_sgrna_barcode_fa.py"
if [[ ! -f "${BUILD_FA}" ]]; then
  echo "ERROR: ${BUILD_FA} not found" >&2
  exit 1
fi

WORKDIR="$(pwd)"
LOCAL_BARCODE_FA="${WORKDIR}/barcode.fa"
mkdir -p "${OUT_DIR}"

if [[ -n "${BARCODE_FA}" && -f "${BARCODE_FA}" ]]; then
  cp "${BARCODE_FA}" "${LOCAL_BARCODE_FA}"
else
  python3 "${BUILD_FA}" --barcode-table "${BARCODE_TABLE}" -o "${LOCAL_BARCODE_FA}"
fi

demux_read() {
  local label="$1"
  local input_fastq="$2"
  local untrimmed="${OUT_DIR}/untrimmed_${label}.fastq.gz"
  echo "Running cutadapt (sgRNA demux ${label}) on ${input_fastq} with -j ${JOBS} ..."
  cutadapt \
    --no-indels \
    -j "${JOBS}" \
    -e "${ERROR_RATE}" \
    -g "file:${LOCAL_BARCODE_FA}" \
    --no-trim \
    --untrimmed-output "${untrimmed}" \
    -o "{name}.${label}.fastq.gz" \
    "${input_fastq}"
}

demux_read "R1" "${INPUT_R1}"
demux_read "R2" "${INPUT_R2}"

while IFS=$'\t' read -r sample_name sample_index _rest; do
  [[ "${sample_name}" == "Sample_Name" || "${sample_name}" == "sample_name" ]] && continue
  [[ -z "${sample_name}" ]] && continue
  sample_dir="${OUT_DIR}/${sample_name}"
  mkdir -p "${sample_dir}"
  for label in R1 R2; do
    src="${WORKDIR}/${sample_name}.${label}.fastq.gz"
    dest="${sample_dir}/${sample_name}.${label}.fastq.gz"
    if [[ -f "${src}" ]]; then
      cp -f "${src}" "${dest}"
    else
      echo "ERROR: no demux ${label} output for ${sample_name} (expected ${src})" >&2
      exit 1
    fi
  done
done < "${BARCODE_TABLE}"

echo "sgRNA cutadapt demux complete under ${OUT_DIR}/"
