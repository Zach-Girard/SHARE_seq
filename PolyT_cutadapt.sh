#!/bin/bash

# Simple Poly-T filtering script intended to be called by Nextflow.
# Usage:
#   PolyT_cutadapt.sh R1.fastq.gz R2.fastq.gz R3.fastq.gz
#
# R3 is treated as the anchor; R1 and R2 are synced based on the
# matched / no-PolyT IDs derived from R3.

set -euo pipefail

if [[ $# -lt 3 || $# -gt 4 ]]; then
  echo "Usage: $0 R1.fastq.gz R2.fastq.gz R3.fastq.gz [OUTDIR]" >&2
  exit 1
fi

R1_IN="$1"
R2_IN="$2"
R3_IN="$3"
OUTDIR="${4:-polyt_filtered}"

mkdir -p "${OUTDIR}"

basename_r3="$(basename "$R3_IN")"
LABEL="${basename_r3%_R3.fastq.gz}"

# --- 2. FILTER R3 (THE ANCHOR) ---
# This removes Poly-T from R3 and creates the two primary buckets.
echo "Filtering R3 for $LABEL..."
cutadapt \
  -g "^NNNNNNNNNNTTTTTT" \
  --no-trim \
  --overlap 16 \
  -e 0.2 \
  -o "${OUTDIR}/${LABEL}.matched.R3.fastq.gz" \
  --untrimmed-output "${OUTDIR}/${LABEL}.noPolyT.R3.fastq.gz" \
  "$R3_IN"

# --- 3. EXTRACT ID LISTS ---
echo "Generating ID lists for $LABEL..."
# Use Python instead of GNU sed/zcat to be portable across macOS/Linux.
python - "$LABEL" "$OUTDIR" <<'PY'
import gzip
import sys
label, outdir = sys.argv[1], sys.argv[2]

def write_ids(in_path, out_path):
    with gzip.open(in_path, "rt") as fin, open(out_path, "w") as fout:
        # FASTQ: 4 lines per record; header is line 1 of each block
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

write_ids(f"{outdir}/{label}.matched.R3.fastq.gz", f"{label}_matched_ids.txt")
write_ids(f"{outdir}/{label}.noPolyT.R3.fastq.gz", f"{label}_noPolyT_ids.txt")
PY

# --- 4. SYNC R1 AND R2 USING PYTHON (no zcat/awk dependencies) ---
python - "$LABEL" "$R1_IN" "$R2_IN" "$OUTDIR" <<'PY'
import gzip
import sys

label, r1_path, r2_path, outdir = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4]

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

matched_ids = load_ids(f"{label}_matched_ids.txt")
no_polyt_ids = load_ids(f"{label}_noPolyT_ids.txt")

print(f"Creating {label}.matched.R2.fastq.gz...")
filter_fastq(r2_path, f"{outdir}/{label}.matched.R2.fastq.gz", matched_ids)
print(f"Creating {label}.noPolyT.R2.fastq.gz...")
filter_fastq(r2_path, f"{outdir}/{label}.noPolyT.R2.fastq.gz", no_polyt_ids)

print(f"Creating {label}.matched.R1.fastq.gz...")
filter_fastq(r1_path, f"{outdir}/{label}.matched.R1.fastq.gz", matched_ids)
print(f"Creating {label}.noPolyT.R1.fastq.gz...")
filter_fastq(r1_path, f"{outdir}/{label}.noPolyT.R1.fastq.gz", no_polyt_ids)
PY

# --- 5. CLEANUP ---
rm "${LABEL}_matched_ids.txt" "${LABEL}_noPolyT_ids.txt"
echo "Workflow complete for $LABEL."
