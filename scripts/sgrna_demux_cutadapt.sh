#!/usr/bin/env bash
# sgRNA demultiplex: cutadapt with anchored sample indices (same as split.lsf).
set -euo pipefail

usage() {
  echo "Usage: $0 --barcode-table TSV --input R1.fastq.gz --out-dir demux" >&2
  exit 1
}

BARCODE_TABLE=""
BARCODE_FA=""
INPUT_R1=""
OUT_DIR="sgRNA/demux"
JOBS=1

while [[ $# -gt 0 ]]; do
  case "$1" in
    --barcode-table) BARCODE_TABLE="$2"; shift 2 ;;
    --barcode-fa) BARCODE_FA="$2"; shift 2 ;;
    --input) INPUT_R1="$2"; shift 2 ;;
    --out-dir) OUT_DIR="$2"; shift 2 ;;
    --jobs) JOBS="$2"; shift 2 ;;
    -h|--help) usage ;;
    *) echo "Unknown option: $1" >&2; usage ;;
  esac
done

[[ -n "${BARCODE_TABLE}" && -n "${INPUT_R1}" ]] || usage
[[ -f "${BARCODE_TABLE}" ]] || { echo "ERROR: barcode table not found: ${BARCODE_TABLE}" >&2; exit 1; }
[[ -f "${INPUT_R1}" ]] || { echo "ERROR: input R1 FASTQ not found: ${INPUT_R1}" >&2; exit 1; }

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

echo "Running cutadapt (sgRNA demux R1-only) with -j ${JOBS} ..."
cutadapt \
  --no-indels \
  -j "${JOBS}" \
  -g "file:${LOCAL_BARCODE_FA}" \
  --no-trim \
  --untrimmed-output "${OUT_DIR}/untrimmed_R1.fastq.gz" \
  -o "{name}.R1.fastq.gz" \
  "${INPUT_R1}"

while IFS=$'\t' read -r sample_name sample_index _rest; do
  [[ "${sample_name}" == "Sample_Name" || "${sample_name}" == "sample_name" ]] && continue
  [[ -z "${sample_name}" ]] && continue
  sample_dir="${OUT_DIR}/${sample_name}"
  mkdir -p "${sample_dir}"
  src="${WORKDIR}/${sample_name}.R1.fastq.gz"
  dest="${sample_dir}/${sample_name}.R1.fastq.gz"
  if [[ -f "${src}" ]]; then
    cp -f "${src}" "${dest}"
  else
    echo "ERROR: no demux R1 output for ${sample_name} (expected ${src})" >&2
    exit 1
  fi
done < "${BARCODE_TABLE}"

echo "sgRNA cutadapt demux complete under ${OUT_DIR}/"
