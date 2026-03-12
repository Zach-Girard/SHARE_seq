#!/usr/bin/env bash
set -euo pipefail

cd "$(dirname "$0")"

mkdir -p GRCh38 GRCm39 hybrid

HUMAN_URL="https://ftp.ensembl.org/pub/current_fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"
MOUSE_URL="https://ftp.ensembl.org/pub/current_fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna.primary_assembly.fa.gz"

HUMAN_FASTA_GZ="GRCh38/$(basename "$HUMAN_URL")"
MOUSE_FASTA_GZ="GRCm39/$(basename "$MOUSE_URL")"
HYBRID_FASTA_GZ="hybrid/hybrid_human_mouse.fa.gz"

DECOMPRESS_HUMAN=("gunzip" "-c")
DECOMPRESS_MOUSE=("gunzip" "-c")
COMPRESS=("gzip" "-c")

if command -v pigz >/dev/null 2>&1; then
  echo "pigz detected – using parallel (de)compression."
  DECOMPRESS_HUMAN=("pigz" "-dc")
  DECOMPRESS_MOUSE=("pigz" "-dc")
  COMPRESS=("pigz" "-c")
fi

echo "Downloading GRCh38 primary assembly..."
if [ ! -f "$HUMAN_FASTA_GZ" ]; then
  curl -kL -o "$HUMAN_FASTA_GZ" "$HUMAN_URL"
else
  echo "  GRCh38 FASTA already exists, skipping download."
fi

echo "Downloading GRCm39 primary assembly..."
if [ ! -f "$MOUSE_FASTA_GZ" ]; then
  curl -kL -o "$MOUSE_FASTA_GZ" "$MOUSE_URL"
else
  echo "  GRCm39 FASTA already exists, skipping download."
fi

echo "Building hybrid genome FASTA with prefixed chromosome names (streaming)..."
rm -f "$HYBRID_FASTA_GZ"

{
  "${DECOMPRESS_HUMAN[@]}" "$HUMAN_FASTA_GZ" | sed 's/^>\([^ ]\+\)/>hs_\1/'
  "${DECOMPRESS_MOUSE[@]}" "$MOUSE_FASTA_GZ" | sed 's/^>\([^ ]\+\)/>mm_\1/'
} | "${COMPRESS[@]}" > "$HYBRID_FASTA_GZ"

echo "Done. Hybrid genome written to ${HYBRID_FASTA_GZ}"
