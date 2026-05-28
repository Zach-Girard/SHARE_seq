#!/usr/bin/env bash
# sgRNA demultiplex: cutadapt with anchored sample indices (same as split.lsf).
set -euo pipefail

usage() {
  echo "Usage: $0 --barcode-table TSv --input R1.fastq.gz --out-dir demux [--error-rate 0.15]" >&2
  exit 1
}

BARCODE_TABLE=""
BARCODE_FA=""
INPUT_FASTQ=""
OUT_DIR="demux"
ERROR_RATE="0.15"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --barcode-table) BARCODE_TABLE="$2"; shift 2 ;;
    --barcode-fa) BARCODE_FA="$2"; shift 2 ;;
    --input) INPUT_FASTQ="$2"; shift 2 ;;
    --out-dir) OUT_DIR="$2"; shift 2 ;;
    --error-rate) ERROR_RATE="$2"; shift 2 ;;
    -h|--help) usage ;;
    *) echo "Unknown option: $1" >&2; usage ;;
  esac
done

[[ -n "${BARCODE_TABLE}" && -n "${INPUT_FASTQ}" ]] || usage
[[ -f "${BARCODE_TABLE}" ]] || { echo "ERROR: barcode table not found: ${BARCODE_TABLE}" >&2; exit 1; }
[[ -f "${INPUT_FASTQ}" ]] || { echo "ERROR: input FASTQ not found: ${INPUT_FASTQ}" >&2; exit 1; }

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BUILD_FA="${SCRIPT_DIR}/build_sgrna_barcode_fa.py"
if [[ ! -f "${BUILD_FA}" ]]; then
  echo "ERROR: ${BUILD_FA} not found" >&2
  exit 1
fi

WORKDIR="$(pwd)"
LOCAL_BARCODE_FA="${WORKDIR}/barcode.fa"
UNTRIMMED="${WORKDIR}/untrimmed.fastq.gz"

if [[ -n "${BARCODE_FA}" && -f "${BARCODE_FA}" ]]; then
  cp "${BARCODE_FA}" "${LOCAL_BARCODE_FA}"
else
  python3 "${BUILD_FA}" --barcode-table "${BARCODE_TABLE}" -o "${LOCAL_BARCODE_FA}"
fi

echo "Running cutadapt (sgRNA demux) on ${INPUT_FASTQ} ..."
cutadapt \
  --no-indels \
  -e "${ERROR_RATE}" \
  -g "file:${LOCAL_BARCODE_FA}" \
  --no-trim \
  --untrimmed-output "${UNTRIMMED}" \
  -o "{name}.fastq.gz" \
  "${INPUT_FASTQ}"

mkdir -p "${OUT_DIR}"
while IFS=$'\t' read -r sample_name sample_index _rest; do
  [[ "${sample_name}" == "Sample_Name" || "${sample_name}" == "sample_name" ]] && continue
  [[ -z "${sample_name}" ]] && continue
  sample_dir="${OUT_DIR}/${sample_name}"
  mkdir -p "${sample_dir}"
  src="${WORKDIR}/${sample_name}.fastq.gz"
  dest="${sample_dir}/${sample_name}.R1.fastq.gz"
  if [[ -f "${src}" ]]; then
    cp -f "${src}" "${dest}"
  else
    echo "WARNING: no demux output for ${sample_name} (expected ${src})" >&2
  fi
done < "${BARCODE_TABLE}"

echo "sgRNA cutadapt demux complete under ${OUT_DIR}/"
