nextflow.enable.dsl=2

/*
 * SHARE-seq processing pipeline.
 *
 * Steps:
 *   1. Demultiplex by sample index barcode, then validate SHARE-seq round barcodes
 *   2. FastQC on demultiplexed reads
 *   3. Poly-T filtering (R2 anchor, R1 synced)
 *   4. Optional fastp trimming (R1 standard; R2 with UMI protection)
 *   5. Barcode prepending (24bp cell barcode from R2 header → R2 sequence)
 *   6. STAR genome index (built or reused via fingerprint)
 *   7. STARsolo alignment (single-end CB_UMI_Complex or paired-end CB_UMI_Simple)
 *   8. QC: knee plots, barnyard plots (hybrid), species-purity splits (hybrid)
 *
 * Requires:
 *   - RAW_FASTQ/ directory with undetermined R1/R2 fastq.gz and sample barcode file
 *   - Genomes/ and GTF/ prepared via helper scripts
 *   - Python dependencies from environment.yml (no local utils.py required)
 */

params.genomes_dir         = params.genomes_dir         ?: 'Genomes'
params.gtf_dir             = params.gtf_dir             ?: 'GTF'
params.umi_len             = params.umi_len             ?: 10
params.total_bc_len        = params.total_bc_len        ?: 24
params.species_model       = params.species_model       ?: 'human'  // 'human', 'mouse', or 'hybrid'
// 'single' = STARsolo CB_UMI_Complex on R1 + withBarcodes_R2 (outputs under STARsolo/)
// 'paired' = STARsolo CB_UMI_Simple with 24bp combinatorial whitelist (outputs under STARsolo_paired/)
// Both modes use trimmed or untrimmed reads depending on params.trim_reads.
params.star_alignment_mode = params.star_alignment_mode ?: 'single' // 'single' or 'paired'
params.barcodes_8bp_file   = params.barcodes_8bp_file   ?: 'barcodes_RC.txt'
params.barcodes_rc         = (params.barcodes_rc in [true, 'true'])
// Optional adapter/quality trimming with fastp (default off).
// R1 is trimmed normally; R2 protects the first `umi_len` bases (UMI) from any trimming.
params.trim_reads          = (params.trim_reads in [true, 'true'])

// Demultiplexing inputs — expected in params.raw_fastq directory
params.raw_fastq           = params.raw_fastq           ?: 'RAW_FASTQ'
params.undetermined_r1     = params.undetermined_r1     ?: null   // R1 fastq.gz filename (inside raw_fastq/)
params.undetermined_r2     = params.undetermined_r2     ?: null   // R2 fastq.gz filename (inside raw_fastq/)
params.sample_barcode_file = params.sample_barcode_file ?: null   // sample index barcode mapping (inside raw_fastq/)
params.demux_mismatches    = params.demux_mismatches    ?: 1      // allowed mismatches for sample index matching
params.split_reads         = params.split_reads         ?: 1000000
params.split_fastq_bin     = params.split_fastq_bin     ?: '/home/yli11/HemTools/bin/splitFastq'

// Derive a single "effective" 8bp whitelist file used everywhere:
//  - RENAME_FASTQ barcode validation (error-correction against whitelist)
//  - Single-end STARsolo (CB_UMI_Complex whitelist)
//  - Building 24bp paired-end STARsolo whitelist
def effectiveCbWhitelistPath = params.barcodes_8bp_file
if (params.barcodes_rc) {
    def inFile  = new File(params.barcodes_8bp_file)
    def outFile = new File((inFile.parent ?: '.') + "/${inFile.name}.rc8bp")

    outFile.withWriter { w ->
        inFile.eachLine { line ->
            def bc = line.trim()
            if (bc) {
                def rc = bc.reverse().tr('ACGTacgtnN', 'TGCAtgcanN')
                w << rc << '\n'
            }
        }
    }
    effectiveCbWhitelistPath = outFile.path
}

// Nextflow tasks run in per-task `work/` directories; make this absolute so STAR
// can always read it via the filesystem regardless of cwd.
effectiveCbWhitelistPath = file(effectiveCbWhitelistPath).toAbsolutePath().toString()

// Log key configuration
log.info "Using Genomes directory      : ${params.genomes_dir}"
log.info "Using GTF directory          : ${params.gtf_dir}"
log.info "Species model                : ${params.species_model}"
log.info "STARsolo alignment mode      : ${params.star_alignment_mode}"
log.info "UMI length (R2)              : ${params.umi_len}"
log.info "Total barcode length         : ${params.total_bc_len}"
log.info "User 8bp barcode file        : ${params.barcodes_8bp_file}"
log.info "Reverse-complement barcodes? : ${params.barcodes_rc}"
log.info "Effective 8bp whitelist file : ${effectiveCbWhitelistPath}"
log.info "Trim reads (fastp)           : ${params.trim_reads}"

/*
 * Channel definitions
 */

def rawDir = params.raw_fastq
log.info "RAW_FASTQ directory           : ${rawDir}"
log.info "Undetermined FASTQs           : explicit --undetermined_r1/--undetermined_r2 or auto-detect *Undetermined*R1*.fastq.gz in ${rawDir}"
def sampleBarcodePath = params.sample_barcode_file ? "${rawDir}/${params.sample_barcode_file}" : 'not set'
log.info "Sample barcode file           : ${sampleBarcodePath}"
log.info "Split reads per chunk         : ${params.split_reads}"
log.info "splitFastq executable         : ${params.split_fastq_bin}"

def loadSampleTypes = { File barcodePath ->
    def sampleTypes = [:]
    barcodePath.eachLine { raw ->
        def line = raw?.trim()
        if (!line || line.startsWith('#')) {
            return
        }
        def cols = line.split(/\t|,/, -1).collect { it.trim() }
        if (cols.size() < 3) {
            return
        }
        def sample = cols[0]
        def sampleType = cols[2].toUpperCase()
        if (!sample || sample.equalsIgnoreCase('sample') || sampleType in ['TYPE', 'SAMPLE_TYPE']) {
            return
        }
        if (sampleType in ['RNA', 'ATAC']) {
            sampleTypes[sample] = sampleType
        }
    }
    return sampleTypes
}

/*
 * Processes
 */

// Step 1: Split undetermined R1/R2 by sample index barcode in the read header.
// Produces per-sample <name>.R1.fastq.gz / <name>.R2.fastq.gz.
process SPLIT_UNDETERMINED_FASTQ {
    tag { pair_id }

    // Default stage-in is symlink; splitFastq may write next to the input path and hit RAW_FASTQ instead of work/.
    stageInMode 'copy'

    publishDir "${projectDir}/${rawDir}", mode: 'copy', overwrite: true

    input:
    tuple val(pair_id), path(r1_undetermined), path(r2_undetermined), path(barcode_file)

    // Capture only split FASTQ chunk files (avoid recursive directory artifacts).
    output:
    tuple val(pair_id), path("split_r1/*.fastq.gz"), path("split_r2/*.fastq.gz")

    script:
    def r1stem = r1_undetermined.name.replaceFirst(/(\.fastq|\.fq)(\.gz)?$/, '')
    def r2stem = r2_undetermined.name.replaceFirst(/(\.fastq|\.fq)(\.gz)?$/, '')
    """
    set -euo pipefail
    module load gcc
    mkdir -p split_r1 split_r2
    ${params.split_fastq_bin} -n ${params.split_reads} -i ${r1_undetermined} -o split_r1/Split_${r1stem}
    ${params.split_fastq_bin} -n ${params.split_reads} -i ${r2_undetermined} -o split_r2/Split_${r2stem}
    # splitFastq may emit plain .fastq (not .gz); normalize to .fastq.gz for demultiplex.py.
    for f in split_r1/*.fastq split_r2/*.fastq; do
      [ -e "\$f" ] || continue
      gzip -f "\$f"
    done
    """
}

process DEMULTIPLEX {
    tag { demux_chunk_id }

    publishDir "${projectDir}/demux", mode: 'copy', overwrite: true

    input:
    tuple val(demux_chunk_id), path(r1_undetermined), path(r2_undetermined), path(barcode_file)

    output:
    path "*.R1.fastq.gz", emit: demux_r1
    path "*.R2.fastq.gz", emit: demux_r2
    path "SHARE-seq.demultiplex.stats.tsv", emit: stats

    """
    python "${projectDir}/demultiplex.py" \\
      -r1 ${r1_undetermined} \\
      -r2 ${r2_undetermined} \\
      -b ${barcode_file} \\
      -n ${params.demux_mismatches} \\
      --revcomp

    rm -f unmatched.R1.fastq.gz unmatched.R2.fastq.gz
    """
}

process MERGE_DEMUX_CHUNKS {
    tag { sample_id }

    publishDir "${projectDir}/demux/${sample_id}", mode: 'copy', overwrite: true

    input:
    tuple val(sample_id), path(r1_parts, stageAs: 'r1_parts??/*'), path(r2_parts, stageAs: 'r2_parts??/*')

    output:
    tuple val(sample_id), path("${sample_id}.R1.fastq.gz"), path("${sample_id}.R2.fastq.gz")

    """
    cat r1_parts*/* > ${sample_id}.R1.fastq.gz
    cat r2_parts*/* > ${sample_id}.R2.fastq.gz
    """
}

// Step 2: Validate three SHARE-seq round barcodes in R1 sequence, rewrite headers
// with matched barcodes, and split into matched/junk output pairs.
// Python dependencies are provided via environment.yml.
process RENAME_FASTQ {
    tag { sample_id }

    publishDir "${projectDir}/demux/${sample_id}", mode: 'copy', overwrite: true

    input:
    tuple val(sample_id), path(r1_demux), path(r2_demux), path(barcode_list)

    output:
    path "${sample_id}.matched.R1.fastq.gz", emit: matched_r1
    path "${sample_id}.matched.R2.fastq.gz", emit: matched_r2
    path "${sample_id}.junk.R1.fastq.gz"
    path "${sample_id}.junk.R2.fastq.gz"
    path "${sample_id}.total_number_reads.tsv", emit: rename_stats

    """
    python "${projectDir}/rename_fastq.py" \\
      -r1 ${r1_demux} \\
      -r2 ${r2_demux} \\
      --sample_ID ${sample_id} \\
      --barcode_list ${barcode_list} \\
      --error 1
    """
}

process FASTQC_DEMUX {
    tag { sample_id }

    publishDir "${projectDir}/fastqc_demux/${sample_id}", mode: 'copy', overwrite: true

    input:
    tuple val(sample_id), path(fastq)

    output:
    path "fastqc_demux/*_fastqc.*", emit: demux_reports

    """
    mkdir -p fastqc_demux
    fastqc -o fastqc_demux -f fastq ${fastq}
    """
}

process POLYT_FILTER {
    tag { sample_id }

    publishDir "${projectDir}/polyt_filtered/${sample_id}", mode: 'copy', overwrite: true

    input:
    tuple val(sample_id), path(r1_fastq), path(r2_fastq)

    output:
    tuple val(sample_id), path("polyt_filtered/*.fastq.gz"), emit: polyt_outputs

    """
    mkdir -p polyt_filtered
    bash "${projectDir}/PolyT_cutadapt.sh" \\
      ${r1_fastq} \\
      ${r2_fastq} \\
      polyt_filtered
    """
}

process TRIM_R1 {
    tag { sample_id }

    publishDir "${projectDir}/trimmed/${sample_id}", mode: 'copy', overwrite: true

    input:
    tuple val(sample_id), path(fastq)

    output:
    tuple val(sample_id), path("trimmed/*.fastq.gz"), emit: trimmed_r1
    path "trimmed/*.fastp.json"
    path "trimmed/*.fastp.html"

    script:
    def stem = fastq.name.replaceFirst(/\.(fastq|fq)(\.gz)?$/, '')
    def outBase = stem.replaceFirst(/\.(R[123])$/, '.trimmed.$1')
    """
    mkdir -p trimmed
    fastp \\
      -i ${fastq} \\
      -o trimmed/${outBase}.fastq.gz \\
      -w ${task.cpus} \\
      --disable_quality_filtering \\
      --disable_length_filtering \\
      -j trimmed/${outBase}.fastp.json \\
      -h trimmed/${outBase}.fastp.html
    """
}

process TRIM_R2_PROTECTED {
    tag { sample_id }

    publishDir "${projectDir}/trimmed/${sample_id}", mode: 'copy', overwrite: true

    input:
    tuple val(sample_id), path(fastq), val(protect_len)

    output:
    tuple val(sample_id), path("trimmed/*.fastq.gz"), emit: trimmed_r2
    path "trimmed/*.fastp.json"
    path "trimmed/*.fastp.html"

    script:
    def stem = fastq.name.replaceFirst(/\.(fastq|fq)(\.gz)?$/, '')
    def outBase = stem.replaceFirst(/\.(R[123])$/, '.trimmed.$1')
    """
    mkdir -p trimmed

    python3 - "${fastq}" "${protect_len}" <<'PY'
import gzip, sys, os

inpath = sys.argv[1]
protect = int(sys.argv[2])

opener = gzip.open if inpath.endswith('.gz') else open

with opener(inpath, 'rt') as fq, \\
     open('_protected.fastq', 'w') as pf, \\
     open('_suffix.fastq', 'w') as sf:
    while True:
        header = fq.readline()
        if not header:
            break
        seq  = fq.readline().rstrip('\\n')
        plus = fq.readline()
        qual = fq.readline().rstrip('\\n')

        pf.write(header)
        pf.write(seq[:protect] + '\\n')
        pf.write(plus)
        pf.write(qual[:protect] + '\\n')

        sf.write(header)
        sf.write(seq[protect:] + '\\n')
        sf.write(plus)
        sf.write(qual[protect:] + '\\n')
PY

    fastp \\
      -i _suffix.fastq \\
      -o _suffix_trimmed.fastq.gz \\
      -w ${task.cpus} \\
      --disable_quality_filtering \\
      --disable_length_filtering \\
      -j trimmed/${outBase}.fastp.json \\
      -h trimmed/${outBase}.fastp.html

    python3 - "_protected.fastq" "_suffix_trimmed.fastq.gz" "trimmed/${outBase}.fastq.gz" <<'PY'
import gzip, sys

prot_path = sys.argv[1]
trim_path = sys.argv[2]
out_path  = sys.argv[3]

with open(prot_path, 'r') as pf, \\
     gzip.open(trim_path, 'rt') as tf, \\
     gzip.open(out_path, 'wt') as out:
    while True:
        p_header = pf.readline()
        if not p_header:
            break
        p_seq  = pf.readline().rstrip('\\n')
        _      = pf.readline()
        p_qual = pf.readline().rstrip('\\n')

        t_header = tf.readline()
        if not t_header:
            out.write(p_header)
            out.write(p_seq + '\\n')
            out.write('+\\n')
            out.write(p_qual + '\\n')
            continue
        t_seq  = tf.readline().rstrip('\\n')
        _      = tf.readline()
        t_qual = tf.readline().rstrip('\\n')

        out.write(p_header)
        out.write(p_seq + t_seq + '\\n')
        out.write('+\\n')
        out.write(p_qual + t_qual + '\\n')
PY

    rm -f _protected.fastq _suffix.fastq
    """
}

process FASTQC_TRIMMED {
    tag { file(fastq).name }

    publishDir "${projectDir}", mode: 'copy', overwrite: true

    input:
    path fastq

    output:
    path "fastqc_trimmed/*_fastqc.*", emit: trimmed_reports

    """
    mkdir -p fastqc_trimmed
    fastqc -o fastqc_trimmed -f fastq ${fastq}
    """
}

// Extract the 24bp cell barcode from the R2 read header (appended by rename_fastq.py
// as @readname_<BC1><BC2><BC3>) and prepend it to the R2 sequence/quality lines.
// Result: each R2 read becomes 24bp_CB + UMI + cDNA.
process PREPEND_HEADER_BARCODES {
    tag { sample_id }

    publishDir "${projectDir}", mode: 'copy', overwrite: true

    input:
    tuple val(sample_id), path(r2_fastq), val(out_dir), val(bc_len)

    output:
    tuple val(sample_id), path("${out_dir}/withBarcodes_*"), emit: r2_with_barcodes

    script:
    def outName = "withBarcodes_${r2_fastq.name}"
    """
    mkdir -p ${out_dir}

    python3 - "${r2_fastq}" "${out_dir}/${outName}" "${bc_len}" <<'PY'
import gzip, sys

in_path  = sys.argv[1]
out_path = sys.argv[2]
bc_len   = int(sys.argv[3])

opener_r = gzip.open if in_path.endswith('.gz') else open
opener_w = gzip.open if out_path.endswith('.gz') else open

with opener_r(in_path, 'rt') as fin, opener_w(out_path, 'wt') as fout:
    while True:
        header = fin.readline()
        if not header:
            break
        seq  = fin.readline().rstrip('\\n')
        plus = fin.readline()
        qual = fin.readline().rstrip('\\n')

        barcode = header.rstrip('\\n').split('_')[-1][:bc_len]
        bc_qual = 'I' * len(barcode)

        fout.write(header)
        fout.write(barcode + seq + '\\n')
        fout.write(plus)
        fout.write(bc_qual + qual + '\\n')
PY
    """
}

process DETERMINE_READ_LENGTH {
    tag { file(r1_fastq).name }

    input:
    path r1_fastq

    output:
    stdout emit: read_length

    """
    python - "${r1_fastq}" <<'PY'
import gzip, sys

path = sys.argv[1]
opener = gzip.open if path.endswith(".gz") else open
with opener(path, "rt") as f:
    f.readline()
    seq = f.readline().rstrip("\\n\\r")
if not seq:
    raise SystemExit(f"Failed to determine read length from {path}")
print(len(seq), end="")
PY
    """
}

process STAR_INDEX {
    tag { params.species_model }

    input:
    tuple val(read_length), path(barcoded_r2_dependency)

    output:
    path "STAR_index_*", emit: star_index

    script:
    def proj = projectDir
    def genomes = params.genomes_dir
    def gtf = params.gtf_dir
    def species = params.species_model

    """
    set -euo pipefail

    readLen=${read_length}
    overhang=\$((readLen - 1))

    case "${species}" in
      mouse)  genomeFasta="${proj}/${genomes}/GRCm39/Mus_musculus.GRCm39.dna.primary_assembly.fa.gz"
              gtfFile="${proj}/${gtf}/Mus_musculus/Mus_musculus.GRCm39.111.gtf.gz"
              indexDir="STAR_index_mouse\${readLen}bp" ;;
      hybrid) genomeFasta="${proj}/${genomes}/hybrid/hybrid_human_mouse.fa.gz"
              gtfFile="${proj}/${gtf}/hybrid/hybrid_human_mouse.gtf.gz"
              indexDir="STAR_index_hybrid\${readLen}bp" ;;
      *)      genomeFasta="${proj}/${genomes}/GRCh38/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"
              gtfFile="${proj}/${gtf}/Homo_sapiens/Homo_sapiens.GRCh38.111.gtf.gz"
              indexDir="STAR_index_human\${readLen}bp" ;;
    esac

    # Build a fingerprint so we can safely reuse an existing index only when
    # species/read length/genome fasta/GTF are unchanged.
    python - "\${genomeFasta}" "\${gtfFile}" "${species}" "\${readLen}" > current_index_fingerprint.json <<'PY'
import json
import os
import sys

genome, gtf, species, read_len = sys.argv[1:5]

def file_meta(path):
    st = os.stat(path)
    return {
        "path": os.path.realpath(path),
        "size": st.st_size,
        "mtime": int(st.st_mtime),
    }

fp = {
    "species_model": species,
    "read_length": int(read_len),
    "genome_fasta": file_meta(genome),
    "gtf": file_meta(gtf),
}
print(json.dumps(fp, sort_keys=True))
PY

    REUSE_DIR="${proj}/\${indexDir}"
    REUSE_FP="\${REUSE_DIR}/.index_fingerprint.json"

    if [[ -d "\${REUSE_DIR}" && -f "\${REUSE_FP}" ]] && cmp -s current_index_fingerprint.json "\${REUSE_FP}"; then
      echo "Reusing existing STAR index: \${REUSE_DIR}"
      ln -s "\${REUSE_DIR}" "\${indexDir}"
      exit 0
    fi

    mkdir -p "\${indexDir}"

    # STAR genomeGenerate expects an uncompressed FASTA (and works best with plain-text GTF).
    if [[ "\${genomeFasta}" == *.gz ]]; then
      gzip -dc "\${genomeFasta}" > genome.fa
    else
      ln -sf "\${genomeFasta}" genome.fa
    fi

    if [[ "\${gtfFile}" == *.gz ]]; then
      gzip -dc "\${gtfFile}" > annotations.gtf
    else
      ln -sf "\${gtfFile}" annotations.gtf
    fi

    STAR --runMode genomeGenerate \\
         --runThreadN ${task.cpus} \\
         --genomeDir \${indexDir} \\
         --genomeFastaFiles genome.fa \\
         --sjdbGTFfile annotations.gtf \\
         --sjdbOverhang \${overhang}

    # Save fingerprint in generated index and mirror index to project root for reuse.
    cp current_index_fingerprint.json "\${indexDir}/.index_fingerprint.json"
    TMP_REUSE_DIR="\${REUSE_DIR}.tmp"
    rm -rf "\${TMP_REUSE_DIR}" "\${REUSE_DIR}"
    cp -R "\${indexDir}" "\${TMP_REUSE_DIR}"
    mv "\${TMP_REUSE_DIR}" "\${REUSE_DIR}"
    """
}

process STARSOLO_SINGLE {
    tag { sample_id }

    publishDir "${projectDir}", mode: 'copy', overwrite: true

    input:
    tuple val(sample_id), path(r1_fastq), path(r2_with_barcodes), path(star_index_dir)

    output:
    path "STARsolo/${sample_id}/", emit: starsolo_out

    """
    mkdir -p STARsolo/${sample_id}
    STAR \\
      --genomeDir ${star_index_dir} \\
      --runThreadN ${task.cpus} \\
      --readFilesCommand zcat \\
      --readFilesIn ${r1_fastq} ${r2_with_barcodes} \\
      --soloType CB_UMI_Complex \\
      --soloCBwhitelist ${effectiveCbWhitelistPath} ${effectiveCbWhitelistPath} ${effectiveCbWhitelistPath} \\
      --soloCBmatchWLtype 1MM \\
      --soloCBposition 0_0_0_7 0_8_0_15 0_16_0_23 \\
      --soloUMIposition 0_24_0_${24 + params.umi_len - 1} \\
      --soloBarcodeReadLength 0 \\
      --soloFeatures Gene GeneFull \\
      --soloStrand Unstranded \\
      --outSAMtype BAM SortedByCoordinate \\
      --outFileNamePrefix STARsolo/${sample_id}/
    """
}

process KNEE_PLOT {
    tag { sample_id }

    publishDir "${projectDir}", mode: 'copy', overwrite: true

    input:
    tuple val(sample_id), path(starsolo_dir), val(qc_base)

    output:
    path "${qc_base}/${sample_id}/${sample_id}_knee_plot.png", emit: knee_plot

    """
    GENEFULL_DIR="${starsolo_dir}/Solo.out/GeneFull"
    SUMMARY_CSV="\${GENEFULL_DIR}/Summary.csv"

    if [ ! -f "\${SUMMARY_CSV}" ]; then
      echo "Summary.csv not found for sample ${sample_id} at \${SUMMARY_CSV}" >&2
      exit 1
    fi

    # Extract "Estimated Number of Cells" (2nd CSV column) without awk regex literals (/.../)
    EST=\$(python - "\${SUMMARY_CSV}" <<'PY'
import csv, sys
path = sys.argv[1]
with open(path, newline="") as f:
    for row in csv.reader(f):
        if not row:
            continue
        if row[0].strip() == "Estimated Number of Cells":
            print((row[1] if len(row) > 1 else "").strip())
            raise SystemExit(0)
raise SystemExit(1)
PY
    )

    if [ -z "\${EST}" ]; then
      echo "Could not extract Estimated Number of Cells from \${SUMMARY_CSV}" >&2
      exit 1
    fi

    python "${projectDir}/Visualization_scripts/Knee_plot.py" \\
      --genefull-dir "\${GENEFULL_DIR}" \\
      --estimated-cells "\${EST}" \\
      --output "${qc_base}/${sample_id}/${sample_id}_knee_plot.png"
    """
}

process BARNYARD_PLOT {
    tag { sample_id }

    publishDir "${projectDir}", mode: 'copy', overwrite: true

    input:
    tuple val(sample_id), path(starsolo_dir), val(qc_base)

    output:
    path "${qc_base}/${sample_id}/80%_collision_plot.png", emit: barnyard80
    path "${qc_base}/${sample_id}/90%_collision_plot.png", emit: barnyard90

    """
    FILTERED_DIR="${starsolo_dir}/Solo.out/GeneFull/filtered"

    if [ ! -d "\${FILTERED_DIR}" ]; then
      echo "Filtered GeneFull directory not found for sample ${sample_id} at \${FILTERED_DIR}" >&2
      exit 1
    fi

    # 80% purity barnyard plot
    python "${projectDir}/Visualization_scripts/BarnyardPlot.py" \\
      --filtered-dir "\${FILTERED_DIR}" \\
      --human-threshold 0.8 \\
      --mouse-threshold 0.2 \\
      --collision-low 0.2 \\
      --collision-high 0.8 \\
      --title "Human-Mouse Collision (80% Purity)" \\
      --output "${qc_base}/${sample_id}/80%_collision_plot.png"

    # 90% purity barnyard plot
    python "${projectDir}/Visualization_scripts/BarnyardPlot.py" \\
      --filtered-dir "\${FILTERED_DIR}" \\
      --human-threshold 0.9 \\
      --mouse-threshold 0.1 \\
      --collision-low 0.1 \\
      --collision-high 0.9 \\
      --title "Human-Mouse Collision (90% Purity)" \\
      --output "${qc_base}/${sample_id}/90%_collision_plot.png"
    """
}

process HYBRID_SPLIT_SPECIES {
    tag { sample_id }

    publishDir "${projectDir}", mode: 'copy', overwrite: true

    input:
    tuple val(sample_id), path(starsolo_dir), val(qc_base)

    output:
    path "${qc_base}/${sample_id}/species_split_purity_0.9", emit: split_09
    path "${qc_base}/${sample_id}/species_split_purity_0.85", emit: split_085
    path "${qc_base}/${sample_id}/species_split_purity_0.8", emit: split_08

    """
    FILTERED_DIR="${starsolo_dir}/Solo.out/GeneFull/filtered"

    if [ ! -d "\${FILTERED_DIR}" ]; then
      echo "Filtered GeneFull directory not found for sample ${sample_id} at \${FILTERED_DIR}" >&2
      exit 1
    fi

    # Purity 0.9
    python "${projectDir}/Split_Species_By_Purity.py" \\
      --input "\${FILTERED_DIR}" \\
      --output "${qc_base}/${sample_id}/species_split_purity_0.9" \\
      --purity 0.9

    # Purity 0.85
    python "${projectDir}/Split_Species_By_Purity.py" \\
      --input "\${FILTERED_DIR}" \\
      --output "${qc_base}/${sample_id}/species_split_purity_0.85" \\
      --purity 0.85

    # Purity 0.8
    python "${projectDir}/Split_Species_By_Purity.py" \\
      --input "\${FILTERED_DIR}" \\
      --output "${qc_base}/${sample_id}/species_split_purity_0.8" \\
      --purity 0.8
    """
}

process BUILD_PAIRED_WHITELIST {
    tag "paired_whitelist"

    input:
    val(dummy)

    output:
    path "whitelist_paired.txt", emit: paired_whitelist_file

    """
    python "${projectDir}/Build_Paired_Whitelist.py" \\
      --barcodes "${effectiveCbWhitelistPath}" \\
      --output whitelist_paired.txt
    """
}

process BUILD_QC_HTML {
    tag "qc_report"

    publishDir "${projectDir}", mode: 'copy', overwrite: true

    input:
    val(report_trigger)

    output:
    path "QC_Report.html", emit: qc_html
    path "QC_Report_assets/*", emit: qc_assets
    path "QC_Report/**", emit: qc_report_dir
    path "QC_Report_bundle.zip", emit: qc_bundle

    """
    python3 - "${projectDir}" <<'PY'
import glob
import html
import os
import csv
import sys
import shutil
import zipfile

proj = sys.argv[1]
out_path = "QC_Report.html"
assets_dir = "QC_Report_assets"
per_sample_dir = "QC_Report"
os.makedirs(assets_dir, exist_ok=True)
os.makedirs(per_sample_dir, exist_ok=True)

def rel_list(pattern):
    return sorted([
        os.path.relpath(p, proj)
        for p in glob.glob(os.path.join(proj, pattern))
        if os.path.isfile(p)
    ])

def safe_asset_name(rel_path):
    name = rel_path.replace("/", "__")
    return name.replace(" ", "_")

def stage_asset(rel_path):
    src = os.path.join(proj, rel_path)
    if not os.path.isfile(src):
        return None
    dest = os.path.join(assets_dir, safe_asset_name(rel_path))
    shutil.copy2(src, dest)
    return dest

def stage_asset_in(rel_path, dest_dir):
    src = os.path.join(proj, rel_path)
    if not os.path.isfile(src):
        return None
    os.makedirs(dest_dir, exist_ok=True)
    dest = os.path.join(dest_dir, safe_asset_name(rel_path))
    shutil.copy2(src, dest)
    return dest

def read_table_preview(rel_path, max_rows=8):
    abs_path = os.path.join(proj, rel_path)
    if not os.path.exists(abs_path):
        return "<p><em>Missing file</em></p>"
    rows = []
    delim = "," if rel_path.endswith(".csv") else "\\t"
    try:
        with open(abs_path, newline="") as fh:
            reader = csv.reader(fh, delimiter=delim)
            for i, row in enumerate(reader):
                rows.append([html.escape(x) for x in row])
                if i + 1 >= max_rows:
                    break
    except Exception as e:
        return f"<p><em>Could not parse preview: {html.escape(str(e))}</em></p>"
    if not rows:
        return "<p><em>File is empty</em></p>"
    cells = []
    for ridx, row in enumerate(rows):
        tag = "th" if ridx == 0 else "td"
        cells.append("<tr>" + "".join(f"<{tag}>{c}</{tag}>" for c in row) + "</tr>")
    return "<table>" + "".join(cells) + "</table>"

def read_text_preview(rel_path, max_lines=80):
    abs_path = os.path.join(proj, rel_path)
    if not os.path.exists(abs_path):
        return "<p><em>Missing file</em></p>"
    try:
        lines = []
        with open(abs_path, "r", errors="replace") as fh:
            for i, line in enumerate(fh):
                lines.append(html.escape(line.rstrip("\\n")))
                if i + 1 >= max_lines:
                    break
        if not lines:
            return "<p><em>File is empty</em></p>"
        return "<pre>" + "\\n".join(lines) + "</pre>"
    except Exception as e:
        return f"<p><em>Could not read file: {html.escape(str(e))}</em></p>"

def links_block(title, paths):
    if not paths:
        return f"<h3>{html.escape(title)}</h3><p><em>No files found.</em></p>"
    items = []
    for p in paths:
        asset = stage_asset(p)
        if asset is None:
            continue
        items.append(f'<li><a href="{html.escape(asset)}">{html.escape(p)}</a></li>')
    if not items:
        return f"<h3>{html.escape(title)}</h3><p><em>No readable files found.</em></p>"
    return f"<h3>{html.escape(title)}</h3><ul>{''.join(items)}</ul>"

def text_files_block(title, paths, max_lines=80):
    if not paths:
        return f"<h3>{html.escape(title)}</h3><p><em>No files found.</em></p>"
    chunks = [f"<h3>{html.escape(title)}</h3>"]
    for p in paths:
        chunks.append(f"<h4>{html.escape(p)}</h4>")
        chunks.append(read_text_preview(p, max_lines=max_lines))
    return "".join(chunks)

def table_files_block(title, paths, max_rows=12):
    if not paths:
        return f"<h3>{html.escape(title)}</h3><p><em>No files found.</em></p>"
    chunks = [f"<h3>{html.escape(title)}</h3>"]
    for p in paths:
        chunks.append(f"<h4>{html.escape(p)}</h4>")
        chunks.append(read_table_preview(p, max_rows=max_rows))
    return "".join(chunks)

def image_files_block(title, paths):
    if not paths:
        return f"<h3>{html.escape(title)}</h3><p><em>No files found.</em></p>"
    chunks = [f"<h3>{html.escape(title)}</h3>"]
    for p in paths:
        asset = stage_asset(p)
        if asset is None:
            continue
        chunks.append(f"<h4>{html.escape(p)}</h4>")
        chunks.append(f'<img src="{html.escape(asset)}" alt="{html.escape(p)}" style="max-width: 1200px; width: 100%; border: 1px solid #ddd; margin-bottom: 14px;" />')
    return "".join(chunks)

demux_total = rel_list("demux/*.total_number_reads.tsv")
demux_stats = rel_list("demux/SHARE-seq.demultiplex.stats.tsv")
fastqc_html = rel_list("fastqc_demux/*_fastqc.html")

starsolo_logs = rel_list("STARsolo/*/Log.final.out") + rel_list("STARsolo_paired/*/Log.final.out")
knee_plots = rel_list("STARsolo/*/*_knee_plot.png") + rel_list("STARsolo_paired/*/*_knee_plot.png")
barcodes_stats = rel_list("STARsolo/*/Solo.out/Barcodes.stats") + rel_list("STARsolo_paired/*/Solo.out/Barcodes.stats")
summary_csv = rel_list("STARsolo/*/Solo.out/GeneFull/Summary.csv") + rel_list("STARsolo_paired/*/Solo.out/GeneFull/Summary.csv")
barnyard = rel_list("STARsolo/*/*collision_plot.png") + rel_list("STARsolo_paired/*/*collision_plot.png")

parts = []
parts.append('''<!doctype html>
<html>
<head>
  <meta charset="utf-8">
  <title>SHARE-seq QC Report</title>
  <style>
    body { font-family: Arial, sans-serif; margin: 24px; line-height: 1.45; }
    h1, h2, h3 { margin-top: 1.1em; }
    code { background: #f3f3f3; padding: 2px 4px; border-radius: 3px; }
    table { border-collapse: collapse; margin: 8px 0 18px 0; }
    th, td { border: 1px solid #ccc; padding: 4px 8px; font-size: 12px; }
  </style>
</head>
<body>
<h1>SHARE-seq QC and Visualization Report</h1>
<p>Generated from pipeline outputs in <code>__PROJECT_DIR__</code>.</p>
'''.replace("__PROJECT_DIR__", html.escape(proj)))

parts.append("<h2>Demultiplexing</h2>")
if demux_stats:
    parts.append("<h3>Demultiplex Stats Preview</h3>")
    parts.append(read_table_preview(demux_stats[0]))
else:
    parts.append("<p><em>No demultiplex stats file found.</em></p>")

parts.append("<h2>FastQC (Demultiplexed)</h2>")
parts.append(links_block("fastqc_demux/*_fastqc.html", fastqc_html))

parts.append("<h2>STARsolo QC</h2>")
parts.append(text_files_block("STARsolo*/<sample>/Log.final.out", starsolo_logs, max_lines=120))
parts.append(image_files_block("STARsolo*/<sample>/*_knee_plot.png", knee_plots))
parts.append(text_files_block("STARsolo*/<sample>/Solo.out/Barcodes.stats", barcodes_stats, max_lines=120))
parts.append(table_files_block("STARsolo*/<sample>/Solo.out/GeneFull/Summary.csv", summary_csv, max_rows=20))

parts.append("<h2>Hybrid Barnyard Plots (if applicable)</h2>")
parts.append(image_files_block("STARsolo*/<sample>/*collision_plot.png", barnyard))

parts.append("</body></html>")

with open(out_path, "w") as out:
    out.write("\\n".join(parts))

def sample_from_starsolo_path(rel_path):
    bits = rel_path.split("/")
    if len(bits) >= 3 and bits[0] in ("STARsolo", "STARsolo_paired"):
        return bits[1]
    return None

all_sample_candidates = set()
for p in starsolo_logs + knee_plots + barcodes_stats + summary_csv + barnyard:
    s = sample_from_starsolo_path(p)
    if s:
        all_sample_candidates.add(s)

for sample in sorted(all_sample_candidates):
    sample_root = os.path.join(per_sample_dir, sample)
    sample_assets = os.path.join(sample_root, "assets")
    os.makedirs(sample_assets, exist_ok=True)

    sample_logs = [p for p in starsolo_logs if sample_from_starsolo_path(p) == sample]
    sample_knee = [p for p in knee_plots if sample_from_starsolo_path(p) == sample]
    sample_barcodes = [p for p in barcodes_stats if sample_from_starsolo_path(p) == sample]
    sample_summary = [p for p in summary_csv if sample_from_starsolo_path(p) == sample]
    sample_barnyard = [p for p in barnyard if sample_from_starsolo_path(p) == sample]

    def sample_image_block(title, paths):
        if not paths:
            return f"<h3>{html.escape(title)}</h3><p><em>No files found.</em></p>"
        chunks = [f"<h3>{html.escape(title)}</h3>"]
        for p in paths:
            asset = stage_asset_in(p, sample_assets)
            if asset is None:
                continue
            rel_asset = os.path.relpath(asset, sample_root)
            chunks.append(f"<h4>{html.escape(p)}</h4>")
            chunks.append(f'<img src="{html.escape(rel_asset)}" alt="{html.escape(p)}" style="max-width: 1200px; width: 100%; border: 1px solid #ddd; margin-bottom: 14px;" />')
        return "".join(chunks)

    sample_parts = []
    sample_parts.append(f'''<!doctype html>
<html>
<head>
  <meta charset="utf-8">
  <title>SHARE-seq QC Report - {html.escape(sample)}</title>
  <style>
    body {{ font-family: Arial, sans-serif; margin: 24px; line-height: 1.45; }}
    h1, h2, h3 {{ margin-top: 1.1em; }}
    code {{ background: #f3f3f3; padding: 2px 4px; border-radius: 3px; }}
    pre {{ background: #f8f8f8; border: 1px solid #ddd; padding: 10px; overflow-x: auto; }}
    table {{ border-collapse: collapse; margin: 8px 0 18px 0; }}
    th, td {{ border: 1px solid #ccc; padding: 4px 8px; font-size: 12px; }}
  </style>
</head>
<body>
<h1>SHARE-seq QC Report: {html.escape(sample)}</h1>
<p>Project path: <code>{html.escape(proj)}</code></p>
''')
    sample_parts.append(text_files_block("Log.final.out", sample_logs, max_lines=120))
    sample_parts.append(sample_image_block("Knee plots", sample_knee))
    sample_parts.append(text_files_block("Barcodes.stats", sample_barcodes, max_lines=120))
    sample_parts.append(table_files_block("GeneFull Summary.csv", sample_summary, max_rows=20))
    sample_parts.append(sample_image_block("Barnyard plots (if hybrid)", sample_barnyard))
    sample_parts.append("</body></html>")

    with open(os.path.join(sample_root, "index.html"), "w") as sf:
        sf.write("\\n".join(sample_parts))

with zipfile.ZipFile("QC_Report_bundle.zip", "w", compression=zipfile.ZIP_DEFLATED) as zf:
    zf.write(out_path, arcname=out_path)
    for root, _, files in os.walk(assets_dir):
        for f in files:
            p = os.path.join(root, f)
            zf.write(p, arcname=p)
    for root, _, files in os.walk(per_sample_dir):
        for f in files:
            p = os.path.join(root, f)
            zf.write(p, arcname=p)
PY
    """
}

process STARSOLO_PAIRED {
    tag { sample_id }

    publishDir "${projectDir}", mode: 'copy', overwrite: true

    input:
    tuple val(sample_id), path(r1_paired), path(r2_paired), path(star_index_dir), path(paired_whitelist)

    output:
    path "STARsolo_paired/${sample_id}/", emit: starsolo_paired_out

    """
    mkdir -p STARsolo_paired/${sample_id}
    STAR \\
      --genomeDir ${star_index_dir} \\
      --runThreadN ${task.cpus} \\
      --readFilesCommand zcat \\
      --readFilesIn ${r1_paired} ${r2_paired} \\
      --soloType CB_UMI_Simple \\
      --soloCBwhitelist ${paired_whitelist} \\
      --soloCBmatchWLtype Exact \\
      --soloCBstart 1 \\
      --soloCBlen 24 \\
      --soloUMIstart 25 \\
      --soloUMIlen ${params.umi_len} \\
      --soloBarcodeReadLength 0 \\
      --soloBarcodeMate 2 \\
      --clip5pNbases 0 ${24 + params.umi_len + 15} \\
      --outFilterMatchNmin 20 \\
      --outFilterMatchNminOverLread 0 \\
      --outFilterScoreMinOverLread 0 \\
      --alignEndsType Local \\
      --soloFeatures Gene GeneFull \\
      --soloStrand Unstranded \\
      --outSAMtype BAM SortedByCoordinate \\
      --outFileNamePrefix STARsolo_paired/${sample_id}/
    """
}

/*
 * Main workflow
 */

workflow {
    main:
    // Validate and build undetermined input pairs here for parser compatibility across Nextflow versions.
    if (!params.sample_barcode_file) {
        error "Demultiplexing requires --sample_barcode_file"
    }
    if ((params.undetermined_r1 && !params.undetermined_r2) || (!params.undetermined_r1 && params.undetermined_r2)) {
        error "Provide both --undetermined_r1 and --undetermined_r2, or neither to auto-detect Undetermined FASTQ pairs."
    }

    def barcodeFile = file("${rawDir}/${params.sample_barcode_file}")
    if (!barcodeFile.exists()) {
        error "Sample barcode file not found: ${barcodeFile}"
    }
    def sampleTypeMap = loadSampleTypes(barcodeFile as File)
    if (!sampleTypeMap) {
        error "No valid sample types found in ${barcodeFile}. Expected column 1 = sample name and column 3 = RNA or ATAC."
    }
    def nRNA = sampleTypeMap.findAll { k, v -> v == 'RNA' }.size()
    def nATAC = sampleTypeMap.findAll { k, v -> v == 'ATAC' }.size()
    log.info "Sample types loaded            : RNA=${nRNA}, ATAC=${nATAC}"

    if (params.undetermined_r1 && params.undetermined_r2) {
        Channel
            .of(tuple(
                params.undetermined_r1.replaceFirst(/(\.fastq|\.fq)(\.gz)?$/, ''),
                file("${rawDir}/${params.undetermined_r1}"),
                file("${rawDir}/${params.undetermined_r2}"),
                barcodeFile
            ))
            .set { ch_undetermined_pairs }
    } else {
        Channel
            .fromPath("${rawDir}/*Undetermined*R1*.fastq.gz")
            .map { r1 ->
                def r2name = r1.name.replaceFirst(/R1/, 'R2')
                def r2 = file("${r1.parent}/${r2name}")
                if (!r2.exists()) {
                    error "Missing R2 pair for ${r1}: expected ${r2}"
                }
                tuple(r1.name.replaceFirst(/(\.fastq|\.fq)(\.gz)?$/, ''), r1, r2, barcodeFile)
            }
            .ifEmpty {
                error "No Undetermined R1 FASTQs found in ${rawDir}. Provide --undetermined_r1/--undetermined_r2 or place files matching *Undetermined*R1*.fastq.gz in RAW_FASTQ."
            }
            .set { ch_undetermined_pairs }
    }

    // Step 1: Split undetermined FASTQs in RAW_FASTQ so demultiplex can fan out in parallel.
    SPLIT_UNDETERMINED_FASTQ(ch_undetermined_pairs)

    SPLIT_UNDETERMINED_FASTQ.out
        .flatMap { pair_id, r1_split_files, r2_split_files ->
            def r1List = (r1_split_files instanceof List) ? r1_split_files : [r1_split_files]
            def r2List = (r2_split_files instanceof List) ? r2_split_files : [r2_split_files]
            def r1Sorted = r1List.toList().sort { it.name }
            def r2Sorted = r2List.toList().sort { it.name }
            if (r1Sorted.size() != r2Sorted.size()) {
                error "Unequal split chunk counts for ${pair_id}: R1=${r1Sorted.size()} R2=${r2Sorted.size()}"
            }
            def chunks = []
            for (int i = 0; i < r1Sorted.size(); i++) {
                chunks << tuple("${pair_id}__chunk${i + 1}", r1Sorted[i], r2Sorted[i], barcodeFile)
            }
            return chunks
        }
        .set { ch_demux_input }

    // Step 2: Demultiplex split chunks by sample index barcode in parallel
    DEMULTIPLEX(ch_demux_input)

    // Merge chunk-level demultiplex output back to one R1/R2 per sample for downstream steps
    DEMULTIPLEX.out.demux_r1
        .flatten()
        .map { f -> tuple(f.name.replaceFirst(/\.R1\.fastq\.gz$/, ''), f) }
        .join(
            DEMULTIPLEX.out.demux_r2.flatten()
                .map { f -> tuple(f.name.replaceFirst(/\.R2\.fastq\.gz$/, ''), f) }
        )
        .groupTuple()
        .map { sample_id, r1_parts, r2_parts ->
            tuple(sample_id, r1_parts.sort { it.name }, r2_parts.sort { it.name })
        }
        .set { ch_demux_chunk_grouped }

    MERGE_DEMUX_CHUNKS(ch_demux_chunk_grouped)

    // Step 3: Pair per-sample R1/R2 and validate SHARE-seq round barcodes
    def cb_whitelist = file(effectiveCbWhitelistPath)
    MERGE_DEMUX_CHUNKS.out
        .map { sample_id, r1, r2 -> tuple(sample_id, r1, r2, cb_whitelist) }
        .set { ch_demux_pairs }

    RENAME_FASTQ(ch_demux_pairs)

    // Attach sample type to matched outputs using sample_barcode_file (col1 sample, col3 RNA/ATAC).
    def ch_matched_r1_typed = RENAME_FASTQ.out.matched_r1
        .map { f ->
            def sample = f.name.replaceFirst(/\.matched\.R1\.fastq\.gz$/, '')
            def sampleType = sampleTypeMap[sample]
            if (!sampleType) {
                error "Sample ${sample} not found in sample_barcode_file with RNA/ATAC type."
            }
            tuple(sample, sampleType, f)
        }

    def ch_matched_r2_typed = RENAME_FASTQ.out.matched_r2
        .map { f ->
            def sample = f.name.replaceFirst(/\.matched\.R2\.fastq\.gz$/, '')
            def sampleType = sampleTypeMap[sample]
            if (!sampleType) {
                error "Sample ${sample} not found in sample_barcode_file with RNA/ATAC type."
            }
            tuple(sample, sampleType, f)
        }

    def ch_matched_pairs_typed = ch_matched_r1_typed
        .join(ch_matched_r2_typed)
        .map { sample_id, sample_type, r1, _, r2 -> tuple(sample_id, sample_type, r1, r2) }

    // ATAC and RNA both run through FASTQC_DEMUX (up to first five processes).
    def ch_input_fastq = ch_matched_pairs_typed
        .flatMap { sample_id, sample_type, r1, r2 -> [tuple(sample_id, r1), tuple(sample_id, r2)] }

    /*
     * Normalize sample identifiers across different FASTQ naming conventions.
     * This is critical because downstream `.join()` operations require exact key matches.
     */
    def normalizeSampleId = { String filenameOrBase ->
        def s = filenameOrBase
        // drop common FASTQ extensions
        s = s.replaceFirst(/(\.fastq|\.fq)(\.gz)?$/, '')
        // drop STARsolo/barcode-prefix artifacts
        s = s.replaceFirst(/^withBarcodes_/, '')
        // drop pipeline-generated tags
        s = s.replaceAll(/(\.matched|\.extracted|\.trimmed|_trimmed)$/, '')
        s = s.replaceAll(/(\.matched|\.extracted|\.trimmed|_trimmed)/, '')
        // drop read designators (common Illumina + simpler forms)
        s = s.replaceFirst(/_R[123](_\d+)?$/, '')
        s = s.replaceFirst(/\.R[123]$/, '')
        // drop residual lane/read counters sometimes left behind
        s = s.replaceFirst(/_\d+$/, '')
        return s
    }

    FASTQC_DEMUX(ch_input_fastq)

    // Only RNA samples continue through the full workflow.
    def ch_rna_pairs = ch_matched_pairs_typed
        .filter { sample_id, sample_type, r1, r2 -> sample_type == 'RNA' }
        .map { sample_id, sample_type, r1, r2 -> tuple(sample_id, r1, r2) }

    // Split demuxed FASTQs: R1 = cDNA, R2 = UMI + cDNA
    def ch_r1_demux = ch_rna_pairs
        .map { sample_id, r1, r2 -> tuple(sample_id, r1) }

    def ch_r2_demux = ch_rna_pairs
        .map { sample_id, r1, r2 -> tuple(sample_id, r2) }

    // Poly-T filtering: R2 is the anchor (UMI + PolyT + cDNA), R1 (cDNA) is synced
    def ch_r1_r2_for_polyt = ch_r1_demux
        .join(ch_r2_demux)
        .map { sample_id, r1, r2 -> tuple(sample_id, r1, r2) }
        .ifEmpty {
            log.warn "No R1/R2 FASTQ pairs found for Poly-T filtering."
            Channel.empty()
        }

    def poly_ch = POLYT_FILTER(ch_r1_r2_for_polyt)
    def ch_polyt_fastq = poly_ch.polyt_outputs
        .flatMap { sample_id, files ->
            def fList = (files instanceof List) ? files : [files]
            fList.collect { f -> tuple(sample_id, f) }
        }
    // R1 = cDNA, R2 = UMI + PolyT + cDNA
    def ch_polyt_r1 = ch_polyt_fastq
        .filter { sample_id, f -> f.name.contains('extracted.R1') }
    def ch_polyt_r2 = ch_polyt_fastq
        .filter { sample_id, f -> f.name.contains('extracted.R2') }

    // Optional fastp trimming (runs before barcode prepend).
    // R1 (cDNA): standard fastp. R2 (UMI+cDNA): first umi_len bp protected.
    // Read-dropping filters are disabled in both to keep R1/R2 in sync.
    def ch_r1_for_downstream
    def ch_r2_for_barcode_prepend
    def barcode_out_dir

    if (params.trim_reads) {
        log.info "Trimming enabled: R1 via fastp; R2 with first ${params.umi_len}bp protected."
        barcode_out_dir = 'trimmed'

        TRIM_R1(ch_polyt_r1)
        ch_r1_for_downstream = TRIM_R1.out.trimmed_r1

        ch_polyt_r2
            .map { sample_id, r2 -> tuple(sample_id, r2, params.umi_len) }
            .set { ch_r2_for_trim }
        TRIM_R2_PROTECTED(ch_r2_for_trim)
        ch_r2_for_barcode_prepend = TRIM_R2_PROTECTED.out.trimmed_r2

        FASTQC_TRIMMED(
            TRIM_R1.out.trimmed_r1
                .mix(TRIM_R2_PROTECTED.out.trimmed_r2)
                .map { sample_id, fastq -> fastq }
        )
    } else {
        log.info "Trimming disabled: using Poly-T–extracted reads directly."
        barcode_out_dir = 'polyt_filtered'
        ch_r1_for_downstream = ch_polyt_r1
        ch_r2_for_barcode_prepend = ch_polyt_r2
    }

    // Read length from R1 (cDNA read) drives STAR index sjdbOverhang
    def ch_r1_for_index = ch_r1_for_downstream.map { sample_id, r1 -> r1 }.take(1)
    DETERMINE_READ_LENGTH(ch_r1_for_index)
    def ch_read_length = DETERMINE_READ_LENGTH.out.read_length.map { it.trim() }

    // Barcode prepend: extract 24bp barcode from R2 header and prepend to R2 sequence.
    // Only R2 is needed; the barcode is already in the header from rename_fastq.py.
    ch_r2_for_barcode_prepend
        .map { sample_id, r2 -> tuple(sample_id, r2, barcode_out_dir, params.total_bc_len) }
        .set { ch_r2_for_barcode }

    PREPEND_HEADER_BARCODES(ch_r2_for_barcode)

    ch_read_length
        .combine(PREPEND_HEADER_BARCODES.out.r2_with_barcodes.map { sample_id, r2wb -> r2wb }.take(1))
        .set { ch_star_index_input }
    STAR_INDEX(ch_star_index_input)

    // Single-end STARsolo (CB_UMI_Complex): R1 (cDNA) + withBarcodes_R2 (24bp CB + UMI + cDNA)
    if( params.star_alignment_mode == 'single' ) {
        log.info "Running single-end STARsolo alignment."

        ch_r1_for_downstream
            .map { sample_id, r1 -> tuple(sample_id, r1) }
            .set { ch_r1_for_align }

        PREPEND_HEADER_BARCODES.out.r2_with_barcodes
            .map { sample_id, r2wb -> tuple(sample_id, r2wb) }
            .set { ch_r2_barcode }

        ch_r1_for_align
            .join(ch_r2_barcode)
            .combine(STAR_INDEX.out.star_index)
            .map { sample_id, r1, r2wb, idx -> tuple(sample_id, r1, r2wb, idx) }
            .set { ch_starsolo_single }

        STARSOLO_SINGLE(ch_starsolo_single)

        STARSOLO_SINGLE.out.starsolo_out
            .map { dir ->
                def sample = dir.baseName  // STARsolo/<sample>/
                tuple(sample, dir, 'STARsolo')
            }
            .set { ch_starsolo_for_hybrid_qc }

        KNEE_PLOT(ch_starsolo_for_hybrid_qc)

        def ch_report_trigger_single = KNEE_PLOT.out.knee_plot
        if( params.species_model == 'hybrid' ) {
            BARNYARD_PLOT(ch_starsolo_for_hybrid_qc)
            HYBRID_SPLIT_SPECIES(ch_starsolo_for_hybrid_qc)
            ch_report_trigger_single = ch_report_trigger_single
                .mix(BARNYARD_PLOT.out.barnyard80)
                .mix(BARNYARD_PLOT.out.barnyard90)
        }
        ch_report_trigger_single.collect().map { 1 }.set { ch_build_qc_report_single }
        BUILD_QC_HTML(ch_build_qc_report_single)
    // Paired-end STARsolo (CB_UMI_Simple): R1 (cDNA) + withBarcodes_R2 (24bp CB + UMI + cDNA).
    // Barcodes are already error-corrected by rename_fastq.py, so no barcode matching needed.
    } else if( params.star_alignment_mode == 'paired' ) {
        log.info "Running paired-end STARsolo CB_UMI_Simple alignment."

        ch_r1_for_downstream
            .map { sample_id, r1 -> tuple(sample_id, r1) }
            .set { ch_r1_for_paired }

        PREPEND_HEADER_BARCODES.out.r2_with_barcodes
            .map { sample_id, r2wb -> tuple(sample_id, r2wb) }
            .set { ch_r2wb_for_paired }

        ch_r1_for_paired
            .join(ch_r2wb_for_paired)
            .combine(STAR_INDEX.out.star_index)
            .map { sample_id, r1, r2wb, idx -> tuple(sample_id, r1, r2wb, idx) }
            .set { ch_starsolo_paired }

        // 24bp whitelist: all 8bp^3 barcode combinations for Exact CB matching
        Channel
            .of(1)
            .set { ch_build_whitelist }
        BUILD_PAIRED_WHITELIST(ch_build_whitelist)

        ch_starsolo_paired
            .combine(BUILD_PAIRED_WHITELIST.out.paired_whitelist_file)
            .set { ch_starsolo_paired_with_wl }

        STARSOLO_PAIRED(ch_starsolo_paired_with_wl)

        STARSOLO_PAIRED.out.starsolo_paired_out
            .map { dir ->
                def sample = dir.baseName
                tuple(sample, dir, 'STARsolo_paired')
            }
            .set { ch_starsolo_paired_for_qc }

        KNEE_PLOT(ch_starsolo_paired_for_qc)

        def ch_report_trigger_paired = KNEE_PLOT.out.knee_plot
        if( params.species_model == 'hybrid' ) {
            BARNYARD_PLOT(ch_starsolo_paired_for_qc)
            HYBRID_SPLIT_SPECIES(ch_starsolo_paired_for_qc)
            ch_report_trigger_paired = ch_report_trigger_paired
                .mix(BARNYARD_PLOT.out.barnyard80)
                .mix(BARNYARD_PLOT.out.barnyard90)
        }
        ch_report_trigger_paired.collect().map { 1 }.set { ch_build_qc_report_paired }
        BUILD_QC_HTML(ch_build_qc_report_paired)
    } else {
        log.info "Unknown star_alignment_mode = ${params.star_alignment_mode}; skipping STARsolo alignment."
    }
}

