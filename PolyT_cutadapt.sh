#!/bin/bash

#BSUB -J PolyT
#BSUB -o %J.out
#BSUB -e %J.err
#BSUB -n 16
#BSUB -R "rusage[mem=3GB]"
#BSUB -W 12:00


module load cutadapt/5.2

# --- 1. SET FILENAMES ---
LABEL="000H7YWLK_LP12066_TATCCTCT_L001"
R1_IN="${LABEL}_R1.fastq.gz"
R2_IN="${LABEL}_R2.fastq.gz"
R3_IN="${LABEL}_R3.fastq.gz"

# --- 2. FILTER R3 (THE ANCHOR) ---
# This removes Poly-T from R3 and creates the two primary buckets.
echo "Filtering R3 for $LABEL..."
cutadapt \
  -g "^NNNNNNNNNNTTTTTT" \
  --no-trim \
  --overlap 16 \
  -e 0.2 \
  -o "${LABEL}.matched.R3.fastq.gz" \
  --untrimmed-output "${LABEL}.noPolyT.R3.fastq.gz" \
  "$R3_IN"

# --- 3. EXTRACT ID LISTS ---
echo "Generating ID lists for $LABEL..."
# We use the newly created labeled files to get the IDs
zcat "${LABEL}.matched.R3.fastq.gz" | sed -n '1~4p' | awk '{print $1}' > "${LABEL}_matched_ids.txt"
zcat "${LABEL}.noPolyT.R3.fastq.gz" | sed -n '1~4p' | awk '{print $1}' > "${LABEL}_noPolyT_ids.txt"

# --- 4. SYNC R1 AND R2 USING AWK ---
# This function pulls the corresponding reads from R1/R2 based on the R3 results
sync_reads() {
    local ID_FILE=$1
    local INPUT_FASTQ=$2
    local OUTPUT_FASTQ=$3
    
    echo "Creating $OUTPUT_FASTQ..."
    zcat "$INPUT_FASTQ" | awk '
        NR==FNR {ids[$1]; next} 
        {
            header=$1
            if (header in ids) {
                print $0;        # Header
                getline; print;  # Sequence
                getline; print;  # Plus
                getline; print;  # Quality
            }
        }' "$ID_FILE" - | gzip > "$OUTPUT_FASTQ"
}

# Apply the sync to R2 (Creates labeled Matched and NoPolyT)
sync_reads "${LABEL}_matched_ids.txt" "$R2_IN" "${LABEL}.matched.R2.fastq.gz"
sync_reads "${LABEL}_noPolyT_ids.txt" "$R2_IN" "${LABEL}.noPolyT.R2.fastq.gz"

# Apply the sync to R1 (Creates labeled Matched and NoPolyT)
sync_reads "${LABEL}_matched_ids.txt" "$R1_IN" "${LABEL}.matched.R1.fastq.gz"
sync_reads "${LABEL}_noPolyT_ids.txt" "$R1_IN" "${LABEL}.noPolyT.R1.fastq.gz"

# --- 5. CLEANUP ---
rm "${LABEL}_matched_ids.txt" "${LABEL}_noPolyT_ids.txt"
echo "Workflow complete for $LABEL."
