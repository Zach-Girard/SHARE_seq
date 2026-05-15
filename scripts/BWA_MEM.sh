#!/bin/bash

#BSUB -J bwa-mem
#BSUB -o %J.out
#BSUB -e %J.err
#BSUB -n 16
#BSUB -R "rusage[mem=10GB]"
#BSUB -W 12:00



# Request ~32GB for indexing

bwa index -p hg38_bwa /research_jude/rgs01_jude/dept/HEM/common/sequencing/chenggrp/zgirard/Genomes/GRCh38.primary_assembly.genome.fa




# 1. Align R1 and R2
# -C flag: Appends the comment from the FASTQ header to the SAM record
# This ensures the barcode (after the space or underscore) is carried over.

bwa mem -t 16 -C hg38_bwa /research_jude/rgs01_jude/dept/HEM/common/sequencing/chenggrp/zgirard/MY_SHARE_CROP/SHARE_seq/220929_NovaSeq/NewSampleSheet/shareseq_yli11_2022-10-06/ATAC_C1200/ATAC_C1200.matched.R1.fastq.gz /research_jude/rgs01_jude/dept/HEM/common/sequencing/chenggrp/zgirard/MY_SHARE_CROP/SHARE_seq/220929_NovaSeq/NewSampleSheet/shareseq_yli11_2022-10-06/ATAC_C1200/ATAC_C1200.matched.R2.fastq.gz > ATAC_C1200.aligned.sam

# 2. Convert to BAM and Filter
# We keep only mapped, high-quality (MAPQ > 30) reads.

samtools view -bS -q 30 ATAC_C1200.aligned.sam | samtools sort -@ 4 -o ATAC_C1200.aligned_sorted.bam -

# 3. Index BAM

samtools index ATAC_C1200.aligned_sorted.bam
