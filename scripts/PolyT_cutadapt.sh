#!/bin/bash

# Poly-T filtering for 2-read SHARE-seq FASTQs.
# Usage:
#   PolyT_cutadapt.sh R1.fastq.gz R2.fastq.gz [OUTDIR]
#
# R1 = cDNA, R2 = UMI + Poly-T + cDNA (cell barcode is in the read header).
# R2 is the anchor: cutadapt matches the UMI+Poly-T pattern in R2.
# R1 is synced based on the Poly-T–extracted / noPolyT IDs derived from R2.

set -euo pipefail

if [[ $# -lt 2 || $# -gt 3 ]]; then
  echo "Usage: $0 R1.fastq.gz R2.fastq.gz [OUTDIR]" >&2
  exit 1
fi

R1_IN="$1"
R2_IN="$2"
OUTDIR="${3:-polyt_filtered}"

mkdir -p "${OUTDIR}"

basename_r2="$(basename "$R2_IN")"
# Strip read designator and extension to get the sample label
LABEL="${basename_r2%.R2.fastq.gz}"
LABEL="${LABEL%_R2.fastq.gz}"

# --- 1. FILTER R2 (THE ANCHOR) ---
echo "Filtering R2 for Poly-T ($LABEL)..."
cutadapt \
  -g "^NNNNNNNNNNTTTTTT" \
  --no-trim \
  --overlap 16 \
  -e 0.2 \
  -o "${OUTDIR}/${LABEL}.extracted.R2.fastq.gz" \
  --untrimmed-output "${OUTDIR}/${LABEL}.noPolyT.R2.fastq.gz" \
  "$R2_IN"

# --- 2. EXTRACT ID LISTS ---
echo "Generating ID lists for $LABEL..."
python3 - "$LABEL" "$OUTDIR" <<'PY'
import gzip, sys

label, outdir = sys.argv[1], sys.argv[2]

def write_ids(in_path, out_path):
    with gzip.open(in_path, "rt") as fin, open(out_path, "w") as fout:
        while True:
            header = fin.readline()
            if not header:
                break
            seq = fin.readline()
            plus = fin.readline()
            qual = fin.readline()
            if not qual:
                break
            fout.write(header.split()[0] + "\n")

write_ids(f"{outdir}/{label}.extracted.R2.fastq.gz", f"{label}_extracted_ids.txt")
write_ids(f"{outdir}/{label}.noPolyT.R2.fastq.gz", f"{label}_noPolyT_ids.txt")
PY

# --- 3. SYNC R1 USING EXTRACTED/NOPOLYT ID SETS ---
python3 - "$LABEL" "$R1_IN" "$OUTDIR" <<'PY'
import gzip, sys

label, r1_path, outdir = sys.argv[1], sys.argv[2], sys.argv[3]

def load_ids(path):
    ids = set()
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line:
                ids.add(line)
    return ids

def filter_fastq(in_path, out_path, keep_ids):
    with gzip.open(in_path, "rt") as fin, gzip.open(out_path, "wt") as fout:
        while True:
            header = fin.readline()
            if not header:
                break
            seq = fin.readline()
            plus = fin.readline()
            qual = fin.readline()
            if not qual:
                break
            hid = header.split()[0]
            if hid in keep_ids:
                fout.write(header)
                fout.write(seq)
                fout.write(plus)
                fout.write(qual)

extracted_ids = load_ids(f"{label}_extracted_ids.txt")
no_polyt_ids = load_ids(f"{label}_noPolyT_ids.txt")

print(f"Creating {label}.extracted.R1.fastq.gz...")
filter_fastq(r1_path, f"{outdir}/{label}.extracted.R1.fastq.gz", extracted_ids)
print(f"Creating {label}.noPolyT.R1.fastq.gz...")
filter_fastq(r1_path, f"{outdir}/{label}.noPolyT.R1.fastq.gz", no_polyt_ids)
PY

# --- 4. CLEANUP ---
rm "${LABEL}_extracted_ids.txt" "${LABEL}_noPolyT_ids.txt"
echo "Poly-T filtering complete for $LABEL."
