nextflow.enable.dsl=2

/*
 * SHARE-seq processing pipeline.
 *
 * Steps:
 *   1. Demultiplex (placeholder pass-through)
 *   2. FastQC on raw reads
 *   3. Poly-T filtering (R1/R2/R3)
 *   4. Optional fastp trimming (R1 full; R3 with UMI protection)
 *   5. Barcode prepending (R2 barcodes → R3)
 *   6. STAR genome index (built or reused via fingerprint)
 *   7. STARsolo alignment (single-end CB_UMI_Complex or paired-end CB_UMI_Simple)
 *   8. QC: knee plots, barnyard plots (hybrid), species-purity splits (hybrid)
 *
 * Requires:
 *   - params.raw_fastq directory with R1/R2/R3 *.fastq.gz files
 *   - Genomes/ and GTF/ prepared via helper scripts
 */

params.raw_fastq           = params.raw_fastq           ?: 'RAW_FASTQ'
params.genomes_dir         = params.genomes_dir         ?: 'Genomes'
params.gtf_dir             = params.gtf_dir             ?: 'GTF'
params.umi_len             = params.umi_len             ?: 10
params.bc_coords           = params.bc_coords           ?: '15-23,53-61,91-99'
params.species_model       = params.species_model       ?: 'human'  // 'human', 'mouse', or 'hybrid'
// 'single' = STARsolo CB_UMI_Complex on R1 + withBarcodes_R3 (outputs under STARsolo/)
// 'paired' = Paired_Barcode_Matcher then STARsolo CB_UMI_Simple (outputs under STARsolo_paired/)
// Both modes use trimmed or untrimmed reads depending on params.trim_reads.
params.star_alignment_mode = params.star_alignment_mode ?: 'single' // 'single' or 'paired'
params.barcodes_8bp_file   = params.barcodes_8bp_file   ?: 'barcodes_RC.txt'
params.barcodes_rc         = (params.barcodes_rc in [true, 'true'])
// Optional adapter/quality trimming with fastp (default off).
// R1 is trimmed normally; R3 protects the first `umi_len` bases (UMI) from any trimming.
params.trim_reads          = (params.trim_reads in [true, 'true'])

// Derive a single "effective" 8bp whitelist file used everywhere:
//  - Single-end STARsolo
//  - Paired-end barcode matching
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
log.info "Using RAW_FASTQ directory    : ${params.raw_fastq}"
log.info "Using Genomes directory      : ${params.genomes_dir}"
log.info "Using GTF directory          : ${params.gtf_dir}"
log.info "Species model                : ${params.species_model}"
log.info "STARsolo alignment mode      : ${params.star_alignment_mode}"
log.info "UMI length (R3)              : ${params.umi_len}"
log.info "Barcode coords (R2)          : ${params.bc_coords}"
log.info "User 8bp barcode file        : ${params.barcodes_8bp_file}"
log.info "Reverse-complement barcodes? : ${params.barcodes_rc}"
log.info "Effective 8bp whitelist file : ${effectiveCbWhitelistPath}"
log.info "Trim reads (fastp)           : ${params.trim_reads}"

/*
 * Channel definitions
 */

// Resolve raw_fastq directory and scan for *.fastq.gz files
def rawDirFs = new File(params.raw_fastq)
log.info "Resolved RAW_FASTQ path          : ${rawDirFs.absolutePath}"
if (!rawDirFs.isDirectory()) {
    log.warn "RAW_FASTQ is not a directory: ${rawDirFs.absolutePath}"
}

Channel
    .fromList(
        ((rawDirFs.isDirectory()
            ? ((rawDirFs.listFiles() ?: []).findAll { it.isFile() && it.name.endsWith('.fastq.gz') })
            : []) as List)
            .collect { it.toString() }
            .sort()
    )
    .ifEmpty { error "No FASTQ files found in: ${params.raw_fastq}" }
    .map { p -> file(p) }
    .set { ch_raw_fastq }

/*
 * Processes
 */

// Placeholder: copies raw FASTQs into demux/. Replace with real demux logic.
process DEMULTIPLEX_PLACEHOLDER {
    tag { file(raw_fastq).name }

    publishDir "${projectDir}", mode: 'copy', overwrite: true

    input:
    path raw_fastq

    output:
    path "demux/${raw_fastq.name}", emit: demuxed_fastq

    """
    mkdir -p demux
    cp ${raw_fastq} demux/
    """
}

process FASTQC_RAW {
    tag { file(fastq).name }

    publishDir "${projectDir}", mode: 'copy', overwrite: true

    input:
    path fastq

    output:
    path "fastqc_raw/${fastq.simpleName}_fastqc.*", emit: raw_reports

    """
    mkdir -p fastqc_raw
    fastqc -o fastqc_raw -f fastq ${fastq}
    """
}

process POLYT_FILTER {
    tag { file(r3_fastq).name }

    publishDir "${projectDir}", mode: 'copy', overwrite: true

    input:
    tuple path(r1_fastq), path(r2_fastq), path(r3_fastq)

    output:
    // Emit individual Poly-T–filtered FASTQs (matched + noPolyT) for downstream use
    path "polyt_filtered/*.fastq.gz", emit: polyt_outputs

    """
    mkdir -p polyt_filtered
    bash "${projectDir}/PolyT_cutadapt.sh" \\
      ${r1_fastq} \\
      ${r2_fastq} \\
      ${r3_fastq} \\
      polyt_filtered
    """
}

process TRIM_R1 {
    tag { file(fastq).name }

    publishDir "${projectDir}", mode: 'copy', overwrite: true

    input:
    path fastq

    output:
    path "trimmed/*.fastq.gz", emit: trimmed_r1
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

process TRIM_R3_PROTECTED {
    tag { file(fastq).name }

    publishDir "${projectDir}", mode: 'copy', overwrite: true

    input:
    tuple path(fastq), val(protect_len)

    output:
    path "trimmed/*.fastq.gz", emit: trimmed_r3
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

process ADD_R2_BARCODES_TO_R3 {
    tag { file(r3_fastq).name }

    publishDir "${projectDir}", mode: 'copy', overwrite: true

    input:
    tuple path(r2_fastq), path(r3_fastq), val(out_dir)

    output:
    path "${out_dir}/withBarcodes_*", emit: r3_with_barcodes

    """
    mkdir -p ${out_dir}
    python "${projectDir}/Read3_Barcode_Addition.py" \\
      -r3 ${r3_fastq} \\
      -r2 ${r2_fastq} \\
      -bc ${params.bc_coords} \\
      -umi 0-${params.umi_len}
    mv withBarcodes_* ${out_dir}/
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
    tuple val(read_length), path(barcoded_r3_dependency)

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
    tuple val(sample_id), path(r1_fastq), path(r3_with_barcodes), path(star_index_dir)

    output:
    path "STARsolo/${sample_id}/", emit: starsolo_out

    """
    mkdir -p STARsolo/${sample_id}
    STAR \\
      --genomeDir ${star_index_dir} \\
      --runThreadN ${task.cpus} \\
      --readFilesCommand zcat \\
      --readFilesIn ${r1_fastq} ${r3_with_barcodes} \\
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

process PAIRED_BARCODE_MATCH {
    tag { sample_id }

    input:
    tuple val(sample_id), path(r1_fastq), path(r3_with_barcodes)

    output:
    tuple val(sample_id), path("paired/${sample_id}_R1.paired.fastq.gz"), path("paired/${sample_id}_R3.paired.fastq.gz"), path("paired/${sample_id}_barcode_match_summary.txt"), emit: paired_reads
    path "paired/${sample_id}_R1.paired.fastq.gz", emit: r1_paired
    path "paired/${sample_id}_R3.paired.fastq.gz", emit: r3_paired
    path "paired/${sample_id}_barcode_match_summary.txt", emit: summary

    """
    mkdir -p paired
    python "${projectDir}/Paired_Barcode_Matcher.py" \\
      --r1 "${r1_fastq}" \\
      --r3 "${r3_with_barcodes}" \\
      --whitelist "${effectiveCbWhitelistPath}" \\
      --out-r1 "paired/${sample_id}_R1.paired.fastq.gz" \\
      --out-r3 "paired/${sample_id}_R3.paired.fastq.gz" \\
      --summary "paired/${sample_id}_barcode_match_summary.txt"
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

process STARSOLO_PAIRED {
    tag { sample_id }

    publishDir "${projectDir}", mode: 'copy', overwrite: true

    input:
    tuple val(sample_id), path(r1_paired), path(r3_paired), path(star_index_dir), path(paired_whitelist)

    output:
    path "STARsolo_paired/${sample_id}/", emit: starsolo_paired_out

    """
    mkdir -p STARsolo_paired/${sample_id}
    STAR \\
      --genomeDir ${star_index_dir} \\
      --runThreadN ${task.cpus} \\
      --readFilesCommand zcat \\
      --readFilesIn ${r1_paired} ${r3_paired} \\
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
    def ch_input_fastq = ch_raw_fastq
    DEMULTIPLEX_PLACEHOLDER(ch_input_fastq)
    ch_input_fastq = DEMULTIPLEX_PLACEHOLDER.out.demuxed_fastq

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
        s = s.replaceAll(/(\.matched|\.trimmed|_trimmed)$/, '')
        s = s.replaceAll(/(\.matched|\.trimmed|_trimmed)/, '')
        // drop read designators (common Illumina + simpler forms)
        s = s.replaceFirst(/_R[123](_\d+)?$/, '')
        s = s.replaceFirst(/\.R[123]$/, '')
        // drop residual lane/read counters sometimes left behind
        s = s.replaceFirst(/_\d+$/, '')
        return s
    }

    FASTQC_RAW(ch_input_fastq)

    // Build R1/R2/R3 triplets from demuxed FASTQs for Poly-T filtering
    def ch_r1_demux = ch_input_fastq
        .filter { f -> f.name.contains('_R1') }
        .map { f -> tuple(f.baseName.replaceFirst('_R1',''), f) }

    def ch_r2_demux = ch_input_fastq
        .filter { f -> f.name.contains('_R2') }
        .map { f -> tuple(f.baseName.replaceFirst('_R2',''), f) }

    def ch_r3_demux = ch_input_fastq
        .filter { f -> f.name.contains('_R3') }
        .map { f -> tuple(f.baseName.replaceFirst('_R3',''), f) }

    def ch_r1_r2_r3_for_polyt = ch_r1_demux
        .join(ch_r2_demux)
        .join(ch_r3_demux)
        .map { sample_id, r1, r2, r3 -> tuple(r1, r2, r3) }
        .ifEmpty {
            log.warn "No R1/R2/R3 FASTQ triplets found in demux/ for Poly-T filtering."
            Channel.empty()
        }

    def poly_ch = POLYT_FILTER(ch_r1_r2_r3_for_polyt)
    def ch_polyt_fastq = poly_ch.polyt_outputs.flatten()
    def ch_polyt_r1 = ch_polyt_fastq.filter { it.name.contains('matched.R1') }
    def ch_polyt_r2 = ch_polyt_fastq.filter { it.name.contains('matched.R2') }
    def ch_polyt_r3 = ch_polyt_fastq.filter { it.name.contains('matched.R3') }

    // Optional fastp trimming (runs before barcode prepend).
    // When enabled: R1 gets standard fastp; R3 protects the first umi_len bp (UMI).
    // Read-dropping filters are disabled in both to keep R1/R3 in sync.
    def ch_r1_for_downstream
    def ch_r3_for_barcode_prepend
    def barcode_out_dir

    if (params.trim_reads) {
        log.info "Trimming enabled: R1 via fastp; R3 with first ${params.umi_len}bp protected."
        barcode_out_dir = 'trimmed'

        TRIM_R1(ch_polyt_r1)
        ch_r1_for_downstream = TRIM_R1.out.trimmed_r1

        ch_polyt_r3
            .map { r3 -> tuple(r3, params.umi_len) }
            .set { ch_r3_for_trim }
        TRIM_R3_PROTECTED(ch_r3_for_trim)
        ch_r3_for_barcode_prepend = TRIM_R3_PROTECTED.out.trimmed_r3

        FASTQC_TRIMMED(TRIM_R1.out.trimmed_r1.mix(TRIM_R3_PROTECTED.out.trimmed_r3))
    } else {
        log.info "Trimming disabled: using Poly-T–matched reads directly."
        barcode_out_dir = 'polyt_filtered'
        ch_r1_for_downstream = ch_polyt_r1
        ch_r3_for_barcode_prepend = ch_polyt_r3
    }

    // R2 is never trimmed — used only for barcode extraction
    def ch_r2_source = ch_polyt_r2

    // Read length from R1 drives STAR index sjdbOverhang
    def ch_r1_for_index = ch_r1_for_downstream.take(1)
    DETERMINE_READ_LENGTH(ch_r1_for_index)
    def ch_read_length = DETERMINE_READ_LENGTH.out.read_length.map { it.trim() }

    // Barcode prepend: R2 barcodes → R3. Output dir matches the R3 source folder.
    def ch_r3_source = ch_r3_for_barcode_prepend

    ch_r2_source
        .map { r2 -> tuple(normalizeSampleId(r2.name), r2) }
        .join(
            ch_r3_source.map { r3 -> tuple(normalizeSampleId(r3.name), r3) }
        )
        .map { sample_id, r2, r3 -> tuple(r2, r3, barcode_out_dir) }
        .set { ch_r2_r3_pairs }

    ADD_R2_BARCODES_TO_R3(ch_r2_r3_pairs)

    ch_read_length
        .combine(ADD_R2_BARCODES_TO_R3.out.r3_with_barcodes.take(1))
        .set { ch_star_index_input }
    STAR_INDEX(ch_star_index_input)

    // Single-end STARsolo (CB_UMI_Complex): R1 + withBarcodes_R3
    if( params.star_alignment_mode == 'single' ) {
        log.info "Running single-end STARsolo alignment."

        def ch_r1_source_for_align = ch_r1_for_downstream

        ch_r1_source_for_align
            .map { r1 ->
                def sample = normalizeSampleId(r1.name)
                tuple(sample, r1)
            }
            .set { ch_r1_for_align }

        ADD_R2_BARCODES_TO_R3.out.r3_with_barcodes
            .map { r3 -> 
                def sample = normalizeSampleId(r3.name)
                tuple(sample, r3)
            }
            .set { ch_r3_barcode }

        ch_r1_for_align
            .join(ch_r3_barcode)
            .map { sample_id, r1, r3 -> tuple(sample_id, r1, r3) }
            .combine(STAR_INDEX.out.star_index)
            .map { sample_id, r1, r3, idx -> tuple(sample_id, r1, r3, idx) }
            .set { ch_starsolo_single }

        STARSOLO_SINGLE(ch_starsolo_single)

        STARSOLO_SINGLE.out.starsolo_out
            .map { dir ->
                def sample = dir.baseName  // STARsolo/<sample>/
                tuple(sample, dir, 'STARsolo')
            }
            .set { ch_starsolo_for_hybrid_qc }

        KNEE_PLOT(ch_starsolo_for_hybrid_qc)

        if( params.species_model == 'hybrid' ) {
            BARNYARD_PLOT(ch_starsolo_for_hybrid_qc)
            HYBRID_SPLIT_SPECIES(ch_starsolo_for_hybrid_qc)
        }
    // Paired-end STARsolo (CB_UMI_Simple): barcode-match R1+withBarcodes_R3,
    // then align matched pairs with a 24bp combinatorial whitelist.
    } else if( params.star_alignment_mode == 'paired' ) {
        log.info "Running paired-end path: barcode matching (3x8bp, 1MM) then STARsolo CB_UMI_Simple."

        def ch_r1_source_for_align_paired = ch_r1_for_downstream

        ch_r1_source_for_align_paired
            .map { r1 -> tuple(normalizeSampleId(r1.name), r1) }
            .set { ch_r1_for_paired }

        ADD_R2_BARCODES_TO_R3.out.r3_with_barcodes
            .map { r3 -> tuple(normalizeSampleId(r3.name), r3) }
            .set { ch_r3_barcode_for_paired }

        ch_r1_for_paired
            .join(ch_r3_barcode_for_paired)
            .map { sample_id, r1, r3 -> tuple(sample_id, r1, r3) }
            .set { ch_pairs_for_matching }

        PAIRED_BARCODE_MATCH(ch_pairs_for_matching)

        PAIRED_BARCODE_MATCH.out.paired_reads
            .map { row ->
                def (sample_id, r1p, r3p, _summary) = row
                tuple(sample_id, r1p, r3p)
            }
            .combine(STAR_INDEX.out.star_index)
            .map { sample_id, r1p, r3p, idx -> tuple(sample_id, r1p, r3p, idx) }
            .ifEmpty {
                log.warn "No paired reads available for STARSOLO_PAIRED (barcode matching channel is empty)."
                Channel.empty()
            }
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

        if( params.species_model == 'hybrid' ) {
            BARNYARD_PLOT(ch_starsolo_paired_for_qc)
            HYBRID_SPLIT_SPECIES(ch_starsolo_paired_for_qc)
        }
    } else {
        log.info "Unknown star_alignment_mode = ${params.star_alignment_mode}; skipping STARsolo alignment."
    }
}

