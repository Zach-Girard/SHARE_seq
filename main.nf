nextflow.enable.dsl=2

/*
 * SHARE-seq processing pipeline.
 *
 * Steps:
 *   1. Split undetermined FASTQs into read-count chunks, then demultiplex chunk pairs in parallel
 *   2. Merge chunk-level demultiplex output back to per-sample FASTQ pairs
 *   3. Validate SHARE-seq round barcodes (rename_fastq.py) and run FastQC on matched reads
 *   4. Route by sample type from sample_barcode_file (column 3):
 *      - ATAC: stop after FASTQC_DEMUX
 *      - RNA : continue through full RNA workflow below
 *   5. Poly-T filtering (R2 anchor, R1 synced) for RNA samples
 *   6. Optional fastp trimming (R1 standard; R2 with UMI protection)
 *   7. Barcode prepending (24bp cell barcode from R2 header → R2 sequence)
 *   8. STAR genome index (built or reused via fingerprint)
 *   9. STARsolo alignment (single-end CB_UMI_Complex or paired-end CB_UMI_Simple)
 *  10. QC: knee plots, barnyard plots (hybrid), species-purity splits (hybrid), multiome cell overlap by group, HTML report outputs
 *
 * Requires:
 *   - RAW_FASTQ/ directory with undetermined R1/R2 fastq.gz and sample barcode file
 *   - sample_barcode_file column 1 = sample name, column 3 = sample type (RNA/ATAC), optional column 4 = Experimental_Group
 *   - Reference FASTA/GTF and BWA/STAR indexes configured in nextflow.config
 *   - Python dependencies from environment.yml (no local utils.py required)
 */

params.umi_len             = params.umi_len             ?: 10
params.total_bc_len        = params.total_bc_len        ?: 24
params.species_model       = params.species_model       ?: 'human'  // 'human', 'mouse', or 'hybrid'
params.human_genome_build  = params.human_genome_build  ?: 'hg19'   // when species_model=human: 'hg19' or 'hg38'
// 'single' = STARsolo CB_UMI_Complex on R1 + withBarcodes_R2 (outputs under STARsolo/)
// 'paired' = STARsolo CB_UMI_Simple with 24bp combinatorial whitelist (outputs under STARsolo_paired/)
// Both modes use trimmed or untrimmed reads depending on params.trim_reads.
params.star_alignment_mode = params.star_alignment_mode ?: 'single' // 'single' or 'paired'
params.barcodes_8bp_file   = params.barcodes_8bp_file   ?: 'barcodes_RC.txt'
params.barcodes_rc         = (params.barcodes_rc in [true, 'true'])
// Optional adapter/quality trimming with fastp (default off).
// R1 is trimmed normally; R2 protects the first `umi_len` bases (UMI) from any trimming.
params.trim_reads          = (params.trim_reads in [true, 'true'])
params.atac_min_frags_for_cell = params.atac_min_frags_for_cell ?: 1000
params.atac_min_tss_for_cell   = params.atac_min_tss_for_cell ?: 4

// Demultiplexing inputs — expected in params.raw_fastq directory
params.raw_fastq           = params.raw_fastq           ?: 'RAW_FASTQ'
params.undetermined_r1     = params.undetermined_r1     ?: null   // R1 fastq.gz filename (inside raw_fastq/)
params.undetermined_r2     = params.undetermined_r2     ?: null   // R2 fastq.gz filename (inside raw_fastq/)
params.sample_barcode_file = params.sample_barcode_file ?: null   // sample index barcode mapping (inside raw_fastq/)
params.demux_mismatches    = params.demux_mismatches    ?: 1      // allowed mismatches for sample index matching
params.split_reads         = params.split_reads         ?: 1000000
params.split_fastq_bin     = params.split_fastq_bin     ?: '/home/yli11/HemTools/bin/splitFastq'

// Fixed reference + index paths (no genome/GTF/index building in pipeline)
params.hybrid_fasta      = params.hybrid_fasta      ?: null
params.hybrid_gtf        = params.hybrid_gtf        ?: null
params.hybrid_bwa_prefix = params.hybrid_bwa_prefix ?: null
params.hybrid_star_index = params.hybrid_star_index ?: null

params.hg19_fasta      = params.hg19_fasta      ?: null
params.hg19_gtf        = params.hg19_gtf        ?: null
params.hg19_bwa_prefix = params.hg19_bwa_prefix ?: null
params.hg19_star_index = params.hg19_star_index ?: null

params.hg38_fasta      = params.hg38_fasta      ?: null
params.hg38_gtf        = params.hg38_gtf        ?: null
params.hg38_bwa_prefix = params.hg38_bwa_prefix ?: null
params.hg38_star_index = params.hg38_star_index ?: null

params.mm10_fasta      = params.mm10_fasta      ?: null
params.mm10_gtf        = params.mm10_gtf        ?: null
params.mm10_bwa_prefix = params.mm10_bwa_prefix ?: null
params.mm10_star_index = params.mm10_star_index ?: null

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

def resolveReferenceConfig = {
    def species = (params.species_model ?: 'human').toString().toLowerCase()
    def humanBuild = (params.human_genome_build ?: 'hg19').toString().toLowerCase()
    def key = (species == 'human') ? humanBuild : species
    if (!(key in ['hybrid', 'hg19', 'hg38', 'mm10'])) {
        error "Unsupported reference selection. species_model=${params.species_model}, human_genome_build=${params.human_genome_build}. Expected hybrid | hg19 | hg38 | mm10."
    }
    def cfg = [
        key        : key,
        fasta      : params["${key}_fasta"],
        gtf        : params["${key}_gtf"],
        bwa_prefix : params["${key}_bwa_prefix"],
        star_index : params["${key}_star_index"],
    ]
    cfg.each { k, v ->
        if (k != 'key' && (!v || v.toString().trim() == '')) {
            error "Missing required reference parameter: ${key}_${k}"
        }
    }
    return cfg
}

def selectedReferenceConfig = resolveReferenceConfig()
def selectedBwaPrefixName = new File(selectedReferenceConfig.bwa_prefix.toString()).name
def bwaRequiredExts = ['.amb', '.ann', '.bwt', '.pac', '.sa']
bwaRequiredExts.each { ext ->
    def p = file("${selectedReferenceConfig.bwa_prefix}${ext}")
    if (!p.exists()) {
        error "Missing BWA index component: ${p}"
    }
}
if (!file(selectedReferenceConfig.star_index.toString()).exists()) {
    error "STAR index directory not found: ${selectedReferenceConfig.star_index}"
}
if (!file(selectedReferenceConfig.fasta.toString()).exists()) {
    error "Reference FASTA not found: ${selectedReferenceConfig.fasta}"
}
if (!file(selectedReferenceConfig.gtf.toString()).exists()) {
    error "Reference GTF not found: ${selectedReferenceConfig.gtf}"
}

// Log key configuration
log.info "Species model                : ${params.species_model}"
log.info "Human genome build           : ${params.human_genome_build}"
log.info "Selected reference set       : ${selectedReferenceConfig.key}"
log.info "Reference FASTA              : ${selectedReferenceConfig.fasta}"
log.info "Reference GTF                : ${selectedReferenceConfig.gtf}"
log.info "BWA index prefix             : ${selectedReferenceConfig.bwa_prefix}"
log.info "STAR index directory         : ${selectedReferenceConfig.star_index}"
log.info "STARsolo alignment mode      : ${params.star_alignment_mode}"
log.info "UMI length (R2)              : ${params.umi_len}"
log.info "Total barcode length         : ${params.total_bc_len}"
log.info "User 8bp barcode file        : ${params.barcodes_8bp_file}"
log.info "Reverse-complement barcodes? : ${params.barcodes_rc}"
log.info "Effective 8bp whitelist file : ${effectiveCbWhitelistPath}"
log.info "Trim reads (fastp)           : ${params.trim_reads}"
log.info "ATAC min frags for cell      : ${params.atac_min_frags_for_cell}"
log.info "ATAC min TSS ratio for cell  : ${params.atac_min_tss_for_cell}"

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

def loadSampleTypes = { barcodePathObj ->
    File barcodePath = (barcodePathObj instanceof File) ? barcodePathObj : new File(barcodePathObj.toString())
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

// Step 1: Split undetermined input FASTQs into smaller chunks for parallel demultiplexing.
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

    // Chunk-level outputs stay in work/; merged stats are published by MERGE_DEMUX_STATS.

    input:
    tuple val(demux_chunk_id), path(r1_undetermined), path(r2_undetermined), path(barcode_file)

    // Chunk-level outputs stay in work/; final sample-level demux outputs are published by MERGE_DEMUX_CHUNKS/RENAME_FASTQ.
    output:
    path "*.R1.fastq.gz", emit: demux_r1
    path "*.R2.fastq.gz", emit: demux_r2

    """
    python3 "${projectDir}/scripts/demultiplex.py" \\
      -r1 ${r1_undetermined} \\
      -r2 ${r2_undetermined} \\
      -b ${barcode_file} \\
      -n ${params.demux_mismatches} \\
      --revcomp
    rm -f SHARE-seq.demultiplex.stats.tsv

    rm -f unmatched.R1.fastq.gz unmatched.R2.fastq.gz
    """
}

process BUILD_DEMUX_STATS_FROM_MERGED {
    tag "demux_stats_from_merged"

    publishDir "${projectDir}/demux", mode: 'copy', overwrite: true

    input:
    path(r1_fastqs)
    path(barcode_file)

    output:
    path "SHARE-seq.demultiplex.stats.tsv", emit: merged_stats

    """
    python3 - <<'PY'
import csv
import gzip
import os
import subprocess

tab = chr(9)
barcode_file = "${barcode_file}"

def read_barcode_table(path):
    rows = []
    with open(path, newline="") as f:
        sample = f.read(4096)
        f.seek(0)
        try:
            dialect = csv.Sniffer().sniff(sample)
            reader = csv.reader(f, dialect)
        except Exception:
            reader = csv.reader(f, delimiter=tab)
        for row in reader:
            row = [x.strip() for x in row if str(x).strip() != ""]
            if row:
                rows.append(row)
    barcode = {}
    sample_type = {}
    for row in rows:
        if len(row) < 2:
            continue
        sid = row[0]
        seq = row[1].upper()
        stype = row[2].strip() if len(row) >= 3 else ""
        barcode[sid] = seq
        sample_type[sid] = stype
    return barcode, sample_type

def count_reads_gz(path):
    # Fast line-count using zcat/wc when available, fallback to Python stream count.
    try:
        out = subprocess.check_output(f"zcat '{path}' | wc -l", shell=True, text=True).strip()
        lines = int(out.split()[0])
        return lines // 4
    except Exception:
        n = 0
        with gzip.open(path, "rt") as fh:
            for _ in fh:
                n += 1
        return n // 4

barcode_by_sample, type_by_sample = read_barcode_table(barcode_file)

rows = []
for fp in sorted([p for p in os.listdir(".") if p.endswith(".R1.fastq.gz")]):
    sample = fp.replace(".R1.fastq.gz", "")
    total = count_reads_gz(fp)
    rows.append((
        barcode_by_sample.get(sample, ""),
        sample,
        type_by_sample.get(sample, ""),
        total
    ))

rows.sort(key=lambda x: x[3], reverse=True)

with open("SHARE-seq.demultiplex.stats.tsv", "w", newline="") as out:
    w = csv.writer(out, delimiter=tab)
    w.writerow(["Sample_Index", "Sample_Name", "Sample_Type", "Total_reads"])
    for r in rows:
        w.writerow(r)
PY
    """
}

process MERGE_DEMUX_CHUNKS {
    tag { sample_id }

    // Final per-sample demultiplexed FASTQs.
    publishDir "${projectDir}/demux/${sample_id}", mode: 'copy', overwrite: true

    input:
    tuple val(sample_id), path(r1_parts, stageAs: 'r1_parts??/*'), path(r2_parts, stageAs: 'r2_parts??/*')

    output:
    tuple val(sample_id), path("${sample_id}.R1.fastq.gz"), path("${sample_id}.R2.fastq.gz")

    """
    cat r1_parts*/* > "${sample_id}.R1.fastq.gz"
    cat r2_parts*/* > "${sample_id}.R2.fastq.gz"
    """
}

// Validate three SHARE-seq round barcodes in R1 sequence, rewrite headers
// with matched barcodes, and split into matched/junk output pairs.
// Python dependencies are provided via environment.yml.
process RENAME_FASTQ {
    tag { sample_id }

    // Keep rename outputs grouped under demux/<sample_id>/.
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
    python3 "${projectDir}/scripts/rename_fastq.py" \\
      -r1 ${r1_demux} \\
      -r2 ${r2_demux} \\
      --sample_ID ${sample_id} \\
      --barcode_list ${barcode_list} \\
      --error 1
    """
}

process FASTQC_DEMUX {
    tag { sample_id }

    // Per-sample FastQC outputs.
    publishDir "${projectDir}/fastqc_demux/${sample_id}", mode: 'copy', overwrite: true

    input:
    tuple val(sample_id), path(fastq)

    output:
    path "*_fastqc.*", emit: demux_reports

    """
    fastqc -o . -f fastq ${fastq}
    """
}

process POLYT_FILTER {
    tag { sample_id }

    // RNA-only branch: per-sample Poly-T outputs under polyt_filtered/<sample_id>/.
    publishDir "${projectDir}/polyt_filtered/${sample_id}", mode: 'copy', overwrite: true

    input:
    tuple val(sample_id), path(r1_fastq), path(r2_fastq)

    output:
    tuple val(sample_id), path("*extracted*.fastq.gz"), emit: polyt_outputs

    """
    bash "${projectDir}/scripts/PolyT_cutadapt.sh" \\
      ${r1_fastq} \\
      ${r2_fastq} \\
      .
    """
}

process TRIM_R1 {
    tag { sample_id }

    // Per-sample trimmed outputs under trimmed/<sample_id>/.
    publishDir "${projectDir}/trimmed/${sample_id}", mode: 'copy', overwrite: true

    input:
    tuple val(sample_id), path(fastq)

    output:
    tuple val(sample_id), path("*.fastq.gz"), emit: trimmed_r1
    path "*.fastp.json"
    path "*.fastp.html"

    script:
    def stem = fastq.name.replaceFirst(/\.(fastq|fq)(\.gz)?$/, '')
    def outBase = stem.replaceFirst(/\.(R[123])$/, '.trimmed.$1')
    """
    fastp \\
      -i ${fastq} \\
      -o ${outBase}.fastq.gz \\
      -w ${task.cpus} \\
      --disable_quality_filtering \\
      --disable_length_filtering \\
      -j ${outBase}.fastp.json \\
      -h ${outBase}.fastp.html
    """
}

process TRIM_R2_PROTECTED {
    tag { sample_id }

    // Per-sample trimmed outputs under trimmed/<sample_id>/.
    publishDir "${projectDir}/trimmed/${sample_id}", mode: 'copy', overwrite: true

    input:
    tuple val(sample_id), path(fastq), val(protect_len)

    output:
    tuple val(sample_id), path("*.fastq.gz"), emit: trimmed_r2
    path "*.fastp.json"
    path "*.fastp.html"

    script:
    def stem = fastq.name.replaceFirst(/\.(fastq|fq)(\.gz)?$/, '')
    def outBase = stem.replaceFirst(/\.(R[123])$/, '.trimmed.$1')
    """
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
      -j ${outBase}.fastp.json \\
      -h ${outBase}.fastp.html

    python3 - "_protected.fastq" "_suffix_trimmed.fastq.gz" "${outBase}.fastq.gz" <<'PY'
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

    // Route withBarcodes outputs to either:
    //   - polyt_filtered/<sample_id>/ (trim disabled)
    //   - trimmed/<sample_id>/        (trim enabled)
    publishDir "${projectDir}", mode: 'copy', overwrite: true, saveAs: { filename -> "${out_dir}/${sample_id}/${filename}" }

    input:
    tuple val(sample_id), path(r2_fastq), val(out_dir), val(bc_len)

    output:
    tuple val(sample_id), path("withBarcodes_*"), emit: r2_with_barcodes

    script:
    def outName = "withBarcodes_${r2_fastq.name}"
    """
    python3 - "${r2_fastq}" "${outName}" "${bc_len}" <<'PY'
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

process STAR_INDEX {
    tag { params.species_model }

    input:
    path barcoded_r2_dependency

    output:
    path "STAR_index_selected", emit: star_index

    script:
    def starIndexSource = selectedReferenceConfig.star_index.toString()

    """
    set -euo pipefail
    [ -d "${starIndexSource}" ] || { echo "ERROR: STAR index directory not found: ${starIndexSource}" >&2; exit 1; }
    [ -f "${starIndexSource}/Genome" ] || { echo "ERROR: STAR index seems invalid (missing Genome file): ${starIndexSource}" >&2; exit 1; }
    ln -s "${starIndexSource}" STAR_index_selected
    """
}

process BWA_INDEX {
    tag { params.species_model }

    input:
    val(dummy)

    output:
    path "BWA_index_selected", emit: bwa_index

    script:
    def bwaPrefix = selectedReferenceConfig.bwa_prefix.toString()
    def bwaIndexDir = new File(bwaPrefix).parent
    """
    set -euo pipefail
    [ -d "${bwaIndexDir}" ] || { echo "ERROR: BWA index directory not found: ${bwaIndexDir}" >&2; exit 1; }
    ln -s "${bwaIndexDir}" BWA_index_selected
    """
}

process BWA_ALIGN_ATAC {
    tag { sample_id }

    publishDir "${projectDir}", mode: 'copy', overwrite: true

    input:
    tuple val(sample_id), path(r1_fastq), path(r2_fastq), val(prefix_name), path(bwa_index_dir)

    output:
    path "ATAC/${sample_id}/", emit: atac_align_out

    """
    set -euo pipefail
    command -v bwa >/dev/null 2>&1 || { echo "ERROR: bwa not found in PATH. Install bwa in environment.yml / activate env."; exit 127; }
    command -v samtools >/dev/null 2>&1 || { echo "ERROR: samtools not found in PATH. Install samtools in environment.yml / activate env."; exit 127; }

    mkdir -p ATAC/${sample_id}
    BWA_PREFIX="${bwa_index_dir}/${prefix_name}"

    bwa mem -t ${task.cpus} -C "\${BWA_PREFIX}" "${r1_fastq}" "${r2_fastq}" \
      | samtools view -@ ${task.cpus} -b -q 30 -F 4 - \
      > ATAC/${sample_id}/${sample_id}.q30.mapped.bam

    # Add CB:Z:<24bp_barcode> from QNAME so duplicate removal is barcode-aware.
    samtools view -@ ${task.cpus} -h ATAC/${sample_id}/${sample_id}.q30.mapped.bam \
      | python3 -c 'import re,sys
pat=re.compile(r"[ACGTNacgtn]{24}")
for raw in sys.stdin:
    if raw.startswith("@"):
        sys.stdout.write(raw); continue
    cols=raw.rstrip("\\n").split("\\t")
    if len(cols)<11:
        continue
    qname=cols[0]
    candidate=qname.rsplit("_",1)[-1] if "_" in qname else qname
    candidate=candidate.strip().split()[0].split("/")[0]
    matches=pat.findall(candidate) or pat.findall(qname)
    if not matches:
        continue
    cb=matches[-1].upper()
    cols=[x for x in cols if not x.startswith("CB:Z:")]
    cols.append("CB:Z:"+cb)
    sys.stdout.write("\\t".join(cols)+"\\n")
' \
      | samtools view -@ ${task.cpus} -b -o ATAC/${sample_id}/${sample_id}.q30.mapped.cbtag.bam -

    total_alignments=\$(samtools view -@ ${task.cpus} -c ATAC/${sample_id}/${sample_id}.q30.mapped.bam)
    cbtag_alignments=\$(samtools view -@ ${task.cpus} -c ATAC/${sample_id}/${sample_id}.q30.mapped.cbtag.bam)
    missing_cb_alignments=\$(( total_alignments - cbtag_alignments ))
    cbtag_pct=\$(python3 - "\${total_alignments}" "\${cbtag_alignments}" <<'PY'
import sys
total = float(sys.argv[1])
kept = float(sys.argv[2])
print(f"{(100.0*kept/total):.4f}" if total > 0 else "0.0000")
PY
    )
    {
      echo -e "Metric\tValue"
      echo -e "TotalAlignments\t\${total_alignments}"
      echo -e "CBTaggedAlignments\t\${cbtag_alignments}"
      echo -e "MissingCBAlignments\t\${missing_cb_alignments}"
      echo -e "CBTaggedPct\t\${cbtag_pct}"
    } > ATAC/${sample_id}/${sample_id}.cbtag_qc.tsv

    samtools sort -@ ${task.cpus} -n \
      -o ATAC/${sample_id}/${sample_id}.q30.namesort.bam \
      ATAC/${sample_id}/${sample_id}.q30.mapped.cbtag.bam

    samtools fixmate -@ ${task.cpus} -m \
      ATAC/${sample_id}/${sample_id}.q30.namesort.bam \
      ATAC/${sample_id}/${sample_id}.q30.fixmate.bam

    samtools sort -@ ${task.cpus} \
      -o ATAC/${sample_id}/${sample_id}.q30.mapped.sorted.bam \
      ATAC/${sample_id}/${sample_id}.q30.fixmate.bam

    # Pre-dedup snapshot (after MAPQ filter, before duplicate removal), from coordinate-sorted BAM.
    samtools index -@ ${task.cpus} ATAC/${sample_id}/${sample_id}.q30.mapped.sorted.bam
    samtools flagstat ATAC/${sample_id}/${sample_id}.q30.mapped.sorted.bam > ATAC/${sample_id}/${sample_id}.q30.mapped.flagstat.txt
    samtools idxstats ATAC/${sample_id}/${sample_id}.q30.mapped.sorted.bam > ATAC/${sample_id}/${sample_id}.q30.mapped.idxstats.txt
    samtools stats ATAC/${sample_id}/${sample_id}.q30.mapped.sorted.bam > ATAC/${sample_id}/${sample_id}.q30.mapped.stats.txt

    samtools markdup -@ ${task.cpus} -r --barcode-tag CB \
      ATAC/${sample_id}/${sample_id}.q30.mapped.sorted.bam \
      ATAC/${sample_id}/${sample_id}.q30.rmdup.sorted.bam

    samtools index -@ ${task.cpus} ATAC/${sample_id}/${sample_id}.q30.rmdup.sorted.bam
    samtools flagstat ATAC/${sample_id}/${sample_id}.q30.rmdup.sorted.bam > ATAC/${sample_id}/${sample_id}.q30.rmdup.flagstat.txt
    samtools idxstats ATAC/${sample_id}/${sample_id}.q30.rmdup.sorted.bam > ATAC/${sample_id}/${sample_id}.q30.rmdup.idxstats.txt
    samtools stats ATAC/${sample_id}/${sample_id}.q30.rmdup.sorted.bam > ATAC/${sample_id}/${sample_id}.q30.rmdup.stats.txt

    rm -f \
      ATAC/${sample_id}/${sample_id}.q30.mapped.bam \
      ATAC/${sample_id}/${sample_id}.q30.mapped.cbtag.bam \
      ATAC/${sample_id}/${sample_id}.q30.namesort.bam \
      ATAC/${sample_id}/${sample_id}.q30.fixmate.bam \
      ATAC/${sample_id}/${sample_id}.q30.mapped.sorted.bai
    """
}

process MULTIQC_ATAC {
    tag "atac_multiqc"

    publishDir "${projectDir}", mode: 'copy', overwrite: true

    input:
    path(atac_dirs)

    output:
    path "multiqc_atac/ATAC_MultiQC.html", emit: multiqc_report
    path "multiqc_atac/ATAC_MultiQC_data", emit: multiqc_data

    """
    set -euo pipefail
    command -v multiqc >/dev/null 2>&1 || { echo "ERROR: multiqc not found in PATH. Install multiqc in environment.yml / activate env."; exit 127; }

    mkdir -p multiqc_atac
    multiqc --force --outdir multiqc_atac --filename ATAC_MultiQC.html .
    """
}

process ESTIMATE_ATAC_CELLS {
    tag { sample_id }

    publishDir "${projectDir}", mode: 'copy', overwrite: true

    input:
    tuple val(sample_id), path(atac_dir), path(cb_whitelist)

    output:
    path "ATAC/${sample_id}/${sample_id}.atac_cells.summary.tsv", emit: atac_cell_summary
    path "ATAC/${sample_id}/${sample_id}.atac_cells.pre_dedup.counts.tsv", emit: atac_cell_counts_pre
    path "ATAC/${sample_id}/${sample_id}.atac_cells.counts.tsv", emit: atac_cell_counts
    path "ATAC/${sample_id}/${sample_id}.archr_tagged.stats.tsv", emit: archr_tagged_stats

    """
    set -euo pipefail
    command -v samtools >/dev/null 2>&1 || { echo "ERROR: samtools not found in PATH."; exit 127; }
    command -v Rscript >/dev/null 2>&1 || { echo "ERROR: Rscript not found in PATH. Install r-base/r-archr in environment.yml / activate env."; exit 127; }

    mkdir -p ATAC/${sample_id}
    BAM_POST="${atac_dir}/${sample_id}.q30.rmdup.sorted.bam"
    BAM_PRE="${atac_dir}/${sample_id}.q30.mapped.sorted.bam"
    [ -f "\$BAM_POST" ] || { echo "Missing post-dedup BAM: \$BAM_POST" >&2; exit 1; }
    [ -f "\$BAM_PRE" ] || { echo "Missing pre-dedup BAM: \$BAM_PRE" >&2; exit 1; }

    case "${params.species_model}" in
      mouse)  archrGenome="mm10" ;;
      human)  case "${params.human_genome_build}" in
                hg38) archrGenome="hg38" ;;
                *)    archrGenome="hg19" ;;
              esac ;;
      *)      archrGenome="hg19" ;;
    esac

    {
      echo -e "Metric\tValue"
      echo -e "InputAlignmentsPreDedup\t\$(samtools view -@ ${task.cpus} -c \"\$BAM_PRE\")"
      echo -e "InputAlignmentsPostDedup\t\$(samtools view -@ ${task.cpus} -c \"\$BAM_POST\")"
    } > ATAC/${sample_id}/${sample_id}.archr_tagged.stats.tsv

    Rscript - "${sample_id}" "\$BAM_PRE" "\$BAM_POST" \\
        "${params.atac_min_frags_for_cell}" "${params.atac_min_tss_for_cell}" "\${archrGenome}" "${task.cpus}" \\
        "ATAC/${sample_id}/${sample_id}.atac_cells.summary.tsv" \\
        "ATAC/${sample_id}/${sample_id}.atac_cells.pre_dedup.counts.tsv" \\
        "ATAC/${sample_id}/${sample_id}.atac_cells.counts.tsv" <<'RSCRIPT'
suppressPackageStartupMessages(library(ArchR))

args <- commandArgs(trailingOnly = TRUE)
sample_id <- args[[1]]
bam_pre <- args[[2]]
bam_post <- args[[3]]
min_frags <- as.numeric(args[[4]])
min_tss <- as.numeric(args[[5]])
archr_genome <- args[[6]]
threads <- as.integer(args[[7]])
out_summary <- args[[8]]
out_counts_pre <- args[[9]]
out_counts_post <- args[[10]]

write_cell_counts <- function(cell_col, out_path) {
  pick_col <- function(candidates) {
    for (candidate in candidates) {
      if (candidate %in% colnames(cell_col)) {
        return(candidate)
      }
    }
    NA_character_
  }
  frag_col <- pick_col(c("nFrags", "NFrags", "nfrags"))
  tss_col <- pick_col(c("TSSEnrichment", "TSS.enrichment", "TSSRatio"))
  cell_names <- rownames(cell_col)
  barcodes <- sub("^.*#", "", cell_names)
  fragments <- if (!is.na(frag_col)) cell_col[[frag_col]] else rep(NA_real_, nrow(cell_col))
  tss <- if (!is.na(tss_col)) cell_col[[tss_col]] else rep(NA_real_, nrow(cell_col))
  counts <- data.frame(
    Barcode = barcodes,
    Fragments = fragments,
    TSSEnrichment = tss,
    stringsAsFactors = FALSE
  )
  utils::write.table(counts, file = out_path, sep = "\t", quote = FALSE, row.names = FALSE)
  nrow(cell_col)
}

if (archr_genome == "hg38" && !requireNamespace("BSgenome.Hsapiens.UCSC.hg38", quietly = TRUE)) {
  stop("Missing package BSgenome.Hsapiens.UCSC.hg38. Install via conda dependency bioconductor-bsgenome.hsapiens.ucsc.hg38.")
}
if (archr_genome == "hg19" && !requireNamespace("BSgenome.Hsapiens.UCSC.hg19", quietly = TRUE)) {
  stop("Missing package BSgenome.Hsapiens.UCSC.hg19. Install via conda dependency bioconductor-bsgenome.hsapiens.ucsc.hg19.")
}
if (archr_genome == "mm10" && !requireNamespace("BSgenome.Mmusculus.UCSC.mm10", quietly = TRUE)) {
  stop("Missing package BSgenome.Mmusculus.UCSC.mm10. Install via conda dependency bioconductor-bsgenome.mmusculus.ucsc.mm10.")
}

addArchRLocking(locking = TRUE)
addArchRThreads(threads = threads)
addArchRGenome(archr_genome)

run_archr <- function(sample_name, bam_path, out_dir) {
  inputFiles <- bam_path
  names(inputFiles) <- sample_name
  arrow_files <- createArrowFiles(
    inputFiles = inputFiles,
    sampleNames = names(inputFiles),
    minFrags = min_frags,
    minTSS = min_tss,
    addTileMat = TRUE,
    addGeneScoreMat = TRUE,
    bcTag = "CB",
    force = TRUE
  )
  proj <- ArchRProject(
    ArrowFiles = arrow_files,
    outputDirectory = out_dir,
    copyArrows = TRUE
  )
  list(arrow_files = arrow_files, cell_col = as.data.frame(getCellColData(proj)))
}

pre_res <- run_archr(sample_id, bam_pre, paste0(sample_id, "_ArchR_preDedup"))
post_res <- run_archr(sample_id, bam_post, paste0(sample_id, "_ArchR_postDedup"))

summarize_archr_cell_col <- function(cell_col) {
  pick_col <- function(candidates) {
    for (candidate in candidates) {
      if (candidate %in% colnames(cell_col)) {
        return(candidate)
      }
    }
    NA_character_
  }
  n <- nrow(cell_col)
  if (n == 0) {
    return(list(
      median_fragments = 0,
      median_tss = 0,
      median_reads_in_tss = 0,
      median_reads_in_prom = 0,
      median_prom_ratio = 0,
      median_reads_in_black = 0,
      median_black_ratio = 0
    ))
  }
  frag_col <- pick_col(c("nFrags", "NFrags", "nfrags"))
  tss_col <- pick_col(c("TSSEnrichment", "TSS.enrichment", "TSSRatio"))
  reads_in_tss_col <- pick_col(c("ReadsInTSS", "readsInTSS"))
  reads_in_prom_col <- pick_col(c("ReadsInPromoter", "readsInPromoter", "ReadsInPromoters"))
  prom_ratio_col <- pick_col(c("PromoterRatio", "promoterRatio", "FRIP.Promoter", "FRIP_Promoter"))
  reads_in_black_col <- pick_col(c("ReadsInBlacklist", "readsInBlacklist"))
  black_ratio_col <- pick_col(c("BlacklistRatio", "blacklistRatio", "FRIP.Blacklist", "FRIP_Blacklist"))
  fragments <- if (!is.na(frag_col)) cell_col[[frag_col]] else rep(NA_real_, n)
  tss <- if (!is.na(tss_col)) cell_col[[tss_col]] else rep(NA_real_, n)
  reads_in_tss <- if (!is.na(reads_in_tss_col)) cell_col[[reads_in_tss_col]] else rep(NA_real_, n)
  reads_in_prom <- if (!is.na(reads_in_prom_col)) cell_col[[reads_in_prom_col]] else rep(NA_real_, n)
  prom_ratio <- if (!is.na(prom_ratio_col)) cell_col[[prom_ratio_col]] else rep(NA_real_, n)
  reads_in_black <- if (!is.na(reads_in_black_col)) cell_col[[reads_in_black_col]] else rep(NA_real_, n)
  black_ratio <- if (!is.na(black_ratio_col)) cell_col[[black_ratio_col]] else rep(NA_real_, n)
  list(
    median_fragments = if (length(fragments) && any(!is.na(fragments))) stats::median(fragments, na.rm = TRUE) else 0,
    median_tss = if (length(tss) && any(!is.na(tss))) stats::median(tss, na.rm = TRUE) else 0,
    median_reads_in_tss = if (length(reads_in_tss) && any(!is.na(reads_in_tss))) stats::median(reads_in_tss, na.rm = TRUE) else 0,
    median_reads_in_prom = if (length(reads_in_prom) && any(!is.na(reads_in_prom))) stats::median(reads_in_prom, na.rm = TRUE) else 0,
    median_prom_ratio = if (length(prom_ratio) && any(!is.na(prom_ratio))) stats::median(prom_ratio, na.rm = TRUE) else 0,
    median_reads_in_black = if (length(reads_in_black) && any(!is.na(reads_in_black))) stats::median(reads_in_black, na.rm = TRUE) else 0,
    median_black_ratio = if (length(black_ratio) && any(!is.na(black_ratio))) stats::median(black_ratio, na.rm = TRUE) else 0
  )
}

pre_m <- summarize_archr_cell_col(pre_res[["cell_col"]])
post_m <- summarize_archr_cell_col(post_res[["cell_col"]])

write_cell_counts(pre_res[["cell_col"]], out_counts_pre)
estimated_cells_post <- write_cell_counts(post_res[["cell_col"]], out_counts_post)
estimated_cells_pre <- nrow(pre_res[["cell_col"]])
estimated_cells <- estimated_cells_post

summary <- data.frame(
  Metric = c(
    "Sample",
    "Method",
    "ArchRGenome",
    "MinFragmentsForCell",
    "MinTSSRatioForCell",
    "EstimatedCellsPreDedup",
    "EstimatedCells",
    "MedianFragmentsPerEstimatedCellPreDedup",
    "MedianFragmentsPerEstimatedCell",
    "MedianTSSEnrichmentPerEstimatedCellPreDedup",
    "MedianTSSEnrichmentPerEstimatedCell",
    "MedianReadsInTSSPerEstimatedCellPreDedup",
    "MedianReadsInTSSPerEstimatedCell",
    "MedianReadsInPromoterPerEstimatedCellPreDedup",
    "MedianReadsInPromoterPerEstimatedCell",
    "MedianPromoterRatioPerEstimatedCellPreDedup",
    "MedianPromoterRatioPerEstimatedCell",
    "MedianReadsInBlacklistPerEstimatedCellPreDedup",
    "MedianReadsInBlacklistPerEstimatedCell",
    "MedianBlacklistRatioPerEstimatedCellPreDedup",
    "MedianBlacklistRatioPerEstimatedCell",
    "ArchRArrowFilePreDedup",
    "ArchRArrowFilePostDedup"
  ),
  Value = c(
    sample_id,
    "ArchR createArrowFiles",
    archr_genome,
    min_frags,
    min_tss,
    estimated_cells_pre,
    estimated_cells,
    pre_m[["median_fragments"]],
    post_m[["median_fragments"]],
    pre_m[["median_tss"]],
    post_m[["median_tss"]],
    pre_m[["median_reads_in_tss"]],
    post_m[["median_reads_in_tss"]],
    pre_m[["median_reads_in_prom"]],
    post_m[["median_reads_in_prom"]],
    pre_m[["median_prom_ratio"]],
    post_m[["median_prom_ratio"]],
    pre_m[["median_reads_in_black"]],
    post_m[["median_reads_in_black"]],
    pre_m[["median_black_ratio"]],
    post_m[["median_black_ratio"]],
    paste(pre_res[["arrow_files"]], collapse = ","),
    paste(post_res[["arrow_files"]], collapse = ",")
  ),
  stringsAsFactors = FALSE
)
utils::write.table(summary, file = out_summary, sep = "\t", quote = FALSE, row.names = FALSE)
RSCRIPT
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
    EST=\$(python3 - "\${SUMMARY_CSV}" <<'PY'
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

    python3 "${projectDir}/scripts/Knee_plot.py" \\
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
    python3 "${projectDir}/scripts/BarnyardPlot.py" \\
      --filtered-dir "\${FILTERED_DIR}" \\
      --human-threshold 0.8 \\
      --mouse-threshold 0.2 \\
      --collision-low 0.2 \\
      --collision-high 0.8 \\
      --title "Human-Mouse Collision (80% Purity)" \\
      --output "${qc_base}/${sample_id}/80%_collision_plot.png"

    # 90% purity barnyard plot
    python3 "${projectDir}/scripts/BarnyardPlot.py" \\
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
    python3 "${projectDir}/scripts/Split_Species_By_Purity.py" \\
      --input "\${FILTERED_DIR}" \\
      --output "${qc_base}/${sample_id}/species_split_purity_0.9" \\
      --purity 0.9

    # Purity 0.85
    python3 "${projectDir}/scripts/Split_Species_By_Purity.py" \\
      --input "\${FILTERED_DIR}" \\
      --output "${qc_base}/${sample_id}/species_split_purity_0.85" \\
      --purity 0.85

    # Purity 0.8
    python3 "${projectDir}/scripts/Split_Species_By_Purity.py" \\
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
    python3 "${projectDir}/scripts/Build_Paired_Whitelist.py" \\
      --barcodes "${effectiveCbWhitelistPath}" \\
      --output whitelist_paired.txt
    """
}

process CELL_OVERLAP_BY_GROUP {
    tag "multiome_overlap"

    publishDir "${projectDir}/multiome_overlap", mode: 'copy', overwrite: true

    input:
    val(trigger)
    path(sample_barcode_file)

    output:
    path "overlap_by_group.tsv", emit: overlap_summary
    path "overlap_by_group.png", emit: overlap_plot, optional: true
    path "*/overlap_summary.tsv", emit: overlap_per_group, optional: true
    path "*/shared_barcodes.txt", emit: overlap_shared_barcodes, optional: true

    """
    set -euo pipefail
    python3 "${projectDir}/scripts/cell_overlap_by_group.py" \\
      --project-dir "${projectDir}" \\
      --sample-barcode-file "${sample_barcode_file}" \\
      --star-alignment-mode "${params.star_alignment_mode}" \\
      --out-dir "."
    """
}

process BUILD_QC_HTML {
    tag "qc_report"

    publishDir "${projectDir}", mode: 'copy', overwrite: true

    input:
    tuple val(report_trigger), path(report_input_files), path(report_builder_script), path(sample_barcode_file)

    output:
    path "QC_Report.html", emit: qc_html
    path "QC_Report_assets/*", emit: qc_assets
    path "QC_Report", emit: qc_report_dir
    path "QC_Report_bundle.zip", emit: qc_bundle

    """
    python3 "${report_builder_script}" \\
        "${projectDir}" \\
        "${params.species_model}" \\
        "${params.star_alignment_mode}" \\
        "${sample_barcode_file}"
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
    def sampleTypeMap = loadSampleTypes(barcodeFile)
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
            r1Sorted.eachWithIndex { r1f, i ->
                chunks << tuple("${pair_id}__chunk${i + 1}", r1f, r2Sorted[i], barcodeFile)
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
    BUILD_DEMUX_STATS_FROM_MERGED(
        MERGE_DEMUX_CHUNKS.out.map { sample_id, r1, r2 -> r1 }.collect(),
        barcodeFile
    )

    // Step 3: Pair per-sample R1/R2 and validate SHARE-seq round barcodes
    def cb_whitelist = file(effectiveCbWhitelistPath)
    MERGE_DEMUX_CHUNKS.out
        .map { sample_id, r1, r2 -> tuple(sample_id, r1, r2, cb_whitelist) }
        .set { ch_demux_pairs }

    RENAME_FASTQ(ch_demux_pairs)

    // Attach sample type to matched outputs using sample_barcode_file (col1 sample, col3 RNA/ATAC; optional col4 Experimental_Group for reporting).
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

    // ATAC branch: align matched R1/R2 with BWA MEM, then MAPQ filter + duplicate removal.
    def ch_atac_pairs_for_trigger = ch_matched_pairs_typed
        .filter { sample_id, sample_type, r1, r2 -> sample_type == 'ATAC' }
        .map { sample_id, sample_type, r1, r2 -> tuple(sample_id, r1, r2) }

    def ch_atac_pairs_for_align = ch_matched_pairs_typed
        .filter { sample_id, sample_type, r1, r2 -> sample_type == 'ATAC' }
        .map { sample_id, sample_type, r1, r2 -> tuple(sample_id, r1, r2) }

    ch_atac_pairs_for_trigger
        .take(1)
        .map { 1 }
        .set { ch_bwa_index_trigger }

    BWA_INDEX(ch_bwa_index_trigger)

    def bwaPrefixName = selectedBwaPrefixName

    ch_atac_pairs_for_align
        .combine(BWA_INDEX.out.bwa_index)
        .map { sample_id, r1, r2, bwa_index_dir ->
            tuple(sample_id, r1, r2, bwaPrefixName, bwa_index_dir)
        }
        .set { ch_atac_bwa_input }

    BWA_ALIGN_ATAC(ch_atac_bwa_input)

    BWA_ALIGN_ATAC.out.atac_align_out
        .map { dir -> tuple(dir.baseName, dir, cb_whitelist) }
        .set { ch_atac_cell_estimate_input }
    ESTIMATE_ATAC_CELLS(ch_atac_cell_estimate_input)

    BWA_ALIGN_ATAC.out.atac_align_out
        .collect()
        .filter { dirs -> dirs && dirs.size() > 0 }
        .set { ch_atac_multiqc_input }
    MULTIQC_ATAC(ch_atac_multiqc_input)

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

    def ch_trim_completion = Channel.empty()

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
        ch_trim_completion = FASTQC_TRIMMED.out.trimmed_reports
    } else {
        log.info "Trimming disabled: using Poly-T–extracted reads directly."
        barcode_out_dir = 'polyt_filtered'
        ch_r1_for_downstream = ch_polyt_r1
        ch_r2_for_barcode_prepend = ch_polyt_r2
    }

    // Barcode prepend: extract 24bp barcode from R2 header and prepend to R2 sequence.
    // Only R2 is needed; the barcode is already in the header from rename_fastq.py.
    ch_r2_for_barcode_prepend
        .map { sample_id, r2 -> tuple(sample_id, r2, barcode_out_dir, params.total_bc_len) }
        .set { ch_r2_for_barcode }

    PREPEND_HEADER_BARCODES(ch_r2_for_barcode)

    PREPEND_HEADER_BARCODES.out.r2_with_barcodes
        .map { sample_id, r2wb -> r2wb }
        .take(1)
        .set { ch_star_index_input }
    STAR_INDEX(ch_star_index_input)

    // Single-end STARsolo (CB_UMI_Complex): R1 (cDNA) + withBarcodes_R2 (24bp CB + UMI + cDNA)
    def ch_report_inputs = BUILD_DEMUX_STATS_FROM_MERGED.out.merged_stats
        .mix(RENAME_FASTQ.out.rename_stats)
        .mix(FASTQC_DEMUX.out.demux_reports)
        .mix(ch_trim_completion)
        .mix(BWA_ALIGN_ATAC.out.atac_align_out)
        .mix(ESTIMATE_ATAC_CELLS.out.atac_cell_summary)
        .mix(ESTIMATE_ATAC_CELLS.out.atac_cell_counts)
        .mix(ESTIMATE_ATAC_CELLS.out.atac_cell_counts_pre)
        .mix(ESTIMATE_ATAC_CELLS.out.archr_tagged_stats)
        .mix(MULTIQC_ATAC.out.multiqc_report)
        .mix(MULTIQC_ATAC.out.multiqc_data)

    def ch_report_barrier = BUILD_DEMUX_STATS_FROM_MERGED.out.merged_stats
        .mix(RENAME_FASTQ.out.rename_stats)
        .mix(FASTQC_DEMUX.out.demux_reports)
        .mix(POLYT_FILTER.out.polyt_outputs)
        .mix(ch_trim_completion)
        .mix(PREPEND_HEADER_BARCODES.out.r2_with_barcodes)
        .mix(STAR_INDEX.out.star_index)
        .mix(BWA_ALIGN_ATAC.out.atac_align_out)
        .mix(ESTIMATE_ATAC_CELLS.out.atac_cell_summary)
        .mix(ESTIMATE_ATAC_CELLS.out.atac_cell_counts)
        .mix(ESTIMATE_ATAC_CELLS.out.atac_cell_counts_pre)
        .mix(ESTIMATE_ATAC_CELLS.out.archr_tagged_stats)
        .mix(MULTIQC_ATAC.out.multiqc_report)
        .mix(MULTIQC_ATAC.out.multiqc_data)

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

        ch_report_barrier = ch_report_barrier
            .mix(STARSOLO_SINGLE.out.starsolo_out)
            .mix(KNEE_PLOT.out.knee_plot)
        ch_report_inputs = ch_report_inputs
            .mix(STARSOLO_SINGLE.out.starsolo_out)
            .mix(KNEE_PLOT.out.knee_plot)
        if( params.species_model == 'hybrid' ) {
            BARNYARD_PLOT(ch_starsolo_for_hybrid_qc)
            HYBRID_SPLIT_SPECIES(ch_starsolo_for_hybrid_qc)
            ch_report_barrier = ch_report_barrier
                .mix(BARNYARD_PLOT.out.barnyard80)
                .mix(BARNYARD_PLOT.out.barnyard90)
                .mix(HYBRID_SPLIT_SPECIES.out.split_09)
                .mix(HYBRID_SPLIT_SPECIES.out.split_085)
                .mix(HYBRID_SPLIT_SPECIES.out.split_08)
            ch_report_inputs = ch_report_inputs
                .mix(BARNYARD_PLOT.out.barnyard80)
                .mix(BARNYARD_PLOT.out.barnyard90)
                .mix(HYBRID_SPLIT_SPECIES.out.split_09)
                .mix(HYBRID_SPLIT_SPECIES.out.split_085)
                .mix(HYBRID_SPLIT_SPECIES.out.split_08)
        }
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

        ch_report_barrier = ch_report_barrier
            .mix(BUILD_PAIRED_WHITELIST.out.paired_whitelist_file)
            .mix(STARSOLO_PAIRED.out.starsolo_paired_out)
            .mix(KNEE_PLOT.out.knee_plot)
        ch_report_inputs = ch_report_inputs
            .mix(STARSOLO_PAIRED.out.starsolo_paired_out)
            .mix(KNEE_PLOT.out.knee_plot)
        if( params.species_model == 'hybrid' ) {
            BARNYARD_PLOT(ch_starsolo_paired_for_qc)
            HYBRID_SPLIT_SPECIES(ch_starsolo_paired_for_qc)
            ch_report_barrier = ch_report_barrier
                .mix(BARNYARD_PLOT.out.barnyard80)
                .mix(BARNYARD_PLOT.out.barnyard90)
                .mix(HYBRID_SPLIT_SPECIES.out.split_09)
                .mix(HYBRID_SPLIT_SPECIES.out.split_085)
                .mix(HYBRID_SPLIT_SPECIES.out.split_08)
            ch_report_inputs = ch_report_inputs
                .mix(BARNYARD_PLOT.out.barnyard80)
                .mix(BARNYARD_PLOT.out.barnyard90)
                .mix(HYBRID_SPLIT_SPECIES.out.split_09)
                .mix(HYBRID_SPLIT_SPECIES.out.split_085)
                .mix(HYBRID_SPLIT_SPECIES.out.split_08)
        }
    } else {
        log.info "Unknown star_alignment_mode = ${params.star_alignment_mode}; skipping STARsolo alignment."
    }

    // RNA/ATAC cell barcode overlap per Experimental_Group (column 4 of sample_barcode_file).
    def ch_overlap_trigger = ch_report_barrier.collect().map { 1 }
    CELL_OVERLAP_BY_GROUP(ch_overlap_trigger, barcodeFile)
    ch_report_barrier = ch_report_barrier.mix(CELL_OVERLAP_BY_GROUP.out.overlap_summary)
    // Only top-level overlap files go to BUILD_QC_HTML (per-group TSVs share basenames and collide).
    // Per-group outputs are published under projectDir/multiome_overlap/ by CELL_OVERLAP_BY_GROUP.
    ch_report_inputs = ch_report_inputs
        .mix(CELL_OVERLAP_BY_GROUP.out.overlap_summary)
        .mix(CELL_OVERLAP_BY_GROUP.out.overlap_plot)

    def ch_report_done = ch_report_barrier
        .collect()
        .map { trigger_items ->
            tuple(1, true)
        }

    def ch_report_files = ch_report_inputs
        .flatMap { x ->
            if (x == null) {
                return []
            }
            if (x instanceof java.nio.file.Path || x instanceof File) {
                return [x]
            }
            if (x instanceof List || x.getClass().isArray()) {
                return x.findAll { it instanceof java.nio.file.Path || it instanceof File }
            }
            return []
        }
        .collect()
        .map { files ->
            def uniq = files.unique { it.toString() }
            tuple(1, uniq)
        }

    ch_report_done
        .join(ch_report_files)
        .map { k, done, report_items -> tuple(1, report_items, file("${projectDir}/scripts/build_qc_report.py"), barcodeFile) }
        .set { ch_build_qc_report }
    BUILD_QC_HTML(ch_build_qc_report)
}

