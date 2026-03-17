#!/usr/bin/env bash
# LSF BSUB directives (used when submitting: bsub < download_gtf.sh)
#BSUB -J download_gtf
#BSUB -n 8
#BSUB -R "rusage[mem=1GB]"
#BSUB -oo download_gtf_%J.out
#BSUB -eo download_gtf_%J.err

set -euo pipefail

cd "$(dirname "$0")"

mkdir -p Homo_sapiens Mus_musculus hybrid

HUMAN_GTF_URL="https://ftp.ensembl.org/pub/release-111/gtf/homo_sapiens/Homo_sapiens.GRCh38.111.gtf.gz"
MOUSE_GTF_URL="https://ftp.ensembl.org/pub/release-111/gtf/mus_musculus/Mus_musculus.GRCm39.111.gtf.gz"

HUMAN_OUT="Homo_sapiens/$(basename "$HUMAN_GTF_URL")"
MOUSE_OUT="Mus_musculus/$(basename "$MOUSE_GTF_URL")"
HYBRID_GTF_GZ="hybrid/hybrid_human_mouse.gtf.gz"

download_file() {
  local url="$1"
  local out="$2"
  local label="$3"

  if [ -f "$out" ]; then
    echo "${label}: ${out} already exists, skipping download."
    return 0
  fi

  echo "${label}: downloading to ${out}..."
  curl -kL -o "$out" "$url"
  echo "${label}: download complete."
}

echo "Starting parallel GTF downloads..."

download_file "$HUMAN_GTF_URL" "$HUMAN_OUT" "Human (GRCh38)" &
PID_HUMAN=$!

download_file "$MOUSE_GTF_URL" "$MOUSE_OUT" "Mouse (GRCm39)" &
PID_MOUSE=$!

wait "$PID_HUMAN"
wait "$PID_MOUSE"

echo "All GTF downloads finished."
echo "  Human GTF:  ${HUMAN_OUT}"
echo "  Mouse GTF:  ${MOUSE_OUT}"

echo "Building hybrid GTF with prefixed chromosome and IDs..."
rm -f "$HYBRID_GTF_GZ"

{
  gunzip -c "$HUMAN_OUT" \
    | awk 'BEGIN{OFS="\t"} \
           /^#/ {print; next} \
           { $1 = "hs_" $1; \
             gsub(/gene_id "([^"]+)"/, "gene_id \"hs_\\1\""); \
             gsub(/transcript_id "([^"]+)"/, "transcript_id \"hs_\\1\""); \
             print }'

  gunzip -c "$MOUSE_OUT" \
    | awk 'BEGIN{OFS="\t"} \
           /^#/ {print; next} \
           { $1 = "mm_" $1; \
             gsub(/gene_id "([^"]+)"/, "gene_id \"mm_\\1\""); \
             gsub(/transcript_id "([^"]+)"/, "transcript_id \"mm_\\1\""); \
             print }'
} | gzip -c > "$HYBRID_GTF_GZ"

echo "Hybrid GTF written to ${HYBRID_GTF_GZ}"

