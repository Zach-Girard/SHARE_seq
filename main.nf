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
 *  10. QC: knee plots, barnyard plots (hybrid), species-purity splits (hybrid), HTML report outputs
 *
 * Requires:
 *   - RAW_FASTQ/ directory with undetermined R1/R2 fastq.gz and sample barcode file
 *   - sample_barcode_file column 1 = sample name, column 3 = sample type (RNA/ATAC)
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
    python "${projectDir}/demultiplex.py" \\
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
    bash "${projectDir}/PolyT_cutadapt.sh" \\
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

process BWA_INDEX {
    tag { params.species_model }

    input:
    val(dummy)

    output:
    path "BWA_index_*", emit: bwa_index

    script:
    def proj = projectDir
    def genomes = params.genomes_dir
    def species = params.species_model
    def prefixName = (species == 'mouse') ? 'mm39_bwa' : (species == 'hybrid' ? 'hybrid_bwa' : 'hg38_bwa')
    """
    set -euo pipefail
    command -v bwa >/dev/null 2>&1 || { echo "ERROR: bwa not found in PATH. Install bwa in environment.yml / activate env."; exit 127; }

    case "${species}" in
      mouse)  genomeFasta="${proj}/${genomes}/GRCm39/Mus_musculus.GRCm39.dna.primary_assembly.fa.gz"
              indexDir="BWA_index_mouse"
              prefixName="mm39_bwa" ;;
      hybrid) genomeFasta="${proj}/${genomes}/hybrid/hybrid_human_mouse.fa.gz"
              indexDir="BWA_index_hybrid"
              prefixName="hybrid_bwa" ;;
      *)      genomeFasta="${proj}/${genomes}/GRCh38/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"
              indexDir="BWA_index_human"
              prefixName="hg38_bwa" ;;
    esac

    python - "\${genomeFasta}" "${species}" > current_bwa_index_fingerprint.json <<'PY'
import json
import os
import sys

genome, species = sys.argv[1:3]
st = os.stat(genome)
fp = {
    "species_model": species,
    "genome_fasta": {
        "path": os.path.realpath(genome),
        "size": st.st_size,
        "mtime": int(st.st_mtime),
    }
}
print(json.dumps(fp, sort_keys=True))
PY

    REUSE_DIR="${proj}/\${indexDir}"
    REUSE_FP="\${REUSE_DIR}/.bwa_index_fingerprint.json"

    if [[ -d "\${REUSE_DIR}" && -f "\${REUSE_FP}" ]] && cmp -s current_bwa_index_fingerprint.json "\${REUSE_FP}"; then
      echo "Reusing existing BWA index: \${REUSE_DIR}"
      cp -R "\${REUSE_DIR}" "\${indexDir}"
      exit 0
    fi

    mkdir -p "\${indexDir}"
    if [[ "\${genomeFasta}" == *.gz ]]; then
      gzip -dc "\${genomeFasta}" > genome.fa
    else
      ln -sf "\${genomeFasta}" genome.fa
    fi

    bwa index -p "\${indexDir}/\${prefixName}" genome.fa
    cp current_bwa_index_fingerprint.json "\${indexDir}/.bwa_index_fingerprint.json"

    TMP_REUSE_DIR="\${REUSE_DIR}.tmp"
    rm -rf "\${TMP_REUSE_DIR}" "\${REUSE_DIR}"
    cp -R "\${indexDir}" "\${TMP_REUSE_DIR}"
    mv "\${TMP_REUSE_DIR}" "\${REUSE_DIR}"
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

    samtools sort -@ ${task.cpus} -n \
      -o ATAC/${sample_id}/${sample_id}.q30.namesort.bam \
      ATAC/${sample_id}/${sample_id}.q30.mapped.bam

    samtools fixmate -@ ${task.cpus} -m \
      ATAC/${sample_id}/${sample_id}.q30.namesort.bam \
      ATAC/${sample_id}/${sample_id}.q30.fixmate.bam

    samtools sort -@ ${task.cpus} \
      -o ATAC/${sample_id}/${sample_id}.q30.possort.bam \
      ATAC/${sample_id}/${sample_id}.q30.fixmate.bam

    # Pre-dedup snapshot (after MAPQ filter, before duplicate removal), from coordinate-sorted BAM.
    samtools index -@ ${task.cpus} ATAC/${sample_id}/${sample_id}.q30.possort.bam
    samtools flagstat ATAC/${sample_id}/${sample_id}.q30.possort.bam > ATAC/${sample_id}/${sample_id}.q30.mapped.flagstat.txt
    samtools idxstats ATAC/${sample_id}/${sample_id}.q30.possort.bam > ATAC/${sample_id}/${sample_id}.q30.mapped.idxstats.txt
    samtools stats ATAC/${sample_id}/${sample_id}.q30.possort.bam > ATAC/${sample_id}/${sample_id}.q30.mapped.stats.txt

    samtools markdup -@ ${task.cpus} -r \
      ATAC/${sample_id}/${sample_id}.q30.possort.bam \
      ATAC/${sample_id}/${sample_id}.q30.rmdup.sorted.bam

    samtools index -@ ${task.cpus} ATAC/${sample_id}/${sample_id}.q30.rmdup.sorted.bam
    samtools flagstat ATAC/${sample_id}/${sample_id}.q30.rmdup.sorted.bam > ATAC/${sample_id}/${sample_id}.q30.rmdup.flagstat.txt
    samtools idxstats ATAC/${sample_id}/${sample_id}.q30.rmdup.sorted.bam > ATAC/${sample_id}/${sample_id}.q30.rmdup.idxstats.txt
    samtools stats ATAC/${sample_id}/${sample_id}.q30.rmdup.sorted.bam > ATAC/${sample_id}/${sample_id}.q30.rmdup.stats.txt

    rm -f \
      ATAC/${sample_id}/${sample_id}.q30.mapped.bam \
      ATAC/${sample_id}/${sample_id}.q30.namesort.bam \
      ATAC/${sample_id}/${sample_id}.q30.fixmate.bam \
      ATAC/${sample_id}/${sample_id}.q30.possort.bam \
      ATAC/${sample_id}/${sample_id}.q30.possort.bam.bai
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
    path "ATAC/${sample_id}/${sample_id}.atac_cells.counts.tsv", emit: atac_cell_counts

    """
    set -euo pipefail
    command -v samtools >/dev/null 2>&1 || { echo "ERROR: samtools not found in PATH."; exit 127; }

    BAM="${atac_dir}/${sample_id}.q30.rmdup.sorted.bam"
    [ -f "$BAM" ] || { echo "Missing BAM: $BAM" >&2; exit 1; }

    case "${params.species_model}" in
      mouse)  gtfFile="${projectDir}/${params.gtf_dir}/Mus_musculus/Mus_musculus.GRCm39.111.gtf.gz" ;;
      hybrid) gtfFile="${projectDir}/${params.gtf_dir}/hybrid/hybrid_human_mouse.gtf.gz" ;;
      *)      gtfFile="${projectDir}/${params.gtf_dir}/Homo_sapiens/Homo_sapiens.GRCh38.111.gtf.gz" ;;
    esac

    samtools view -@ ${task.cpus} -f 64 -F 2304 "$BAM" \
      | python3 - "${sample_id}" "${cb_whitelist}" "${gtfFile}" "${params.atac_min_frags_for_cell}" "${params.atac_min_tss_for_cell}" \
          "ATAC/${sample_id}/${sample_id}.atac_cells.summary.tsv" \
          "ATAC/${sample_id}/${sample_id}.atac_cells.counts.tsv" <<'PY'
import sys
import statistics
import gzip
import bisect

sample_id, wl_path, gtf_path, min_frags_s, min_tss_s, out_summary, out_counts = sys.argv[1:8]
min_frags = int(min_frags_s)
min_tss = float(min_tss_s)

valid = set()
with open(wl_path, "r", errors="replace") as fh:
    for line in fh:
        s = line.strip()
        if not s:
            continue
        valid.add(s.split()[0])

def load_tss_positions(path):
    opener = gzip.open if path.endswith(".gz") else open
    out = {}
    with opener(path, "rt", errors="replace") as fh:
        for raw in fh:
            if not raw or raw.startswith("#"):
                continue
            cols = raw.rstrip("\n").split("\t")
            if len(cols) < 9 or cols[2] != "gene":
                continue
            chrom = cols[0]
            try:
                start = int(cols[3])
                end = int(cols[4])
            except Exception:
                continue
            strand = cols[6]
            tss = start if strand == "+" else end
            out.setdefault(chrom, []).append(tss)
    for chrom in out:
        out[chrom].sort()
    return out

def classify_tss_hit(tss_list, pos):
    if not tss_list:
        return (False, False)
    i = bisect.bisect_left(tss_list, pos)
    left = i - 1
    right = i
    is_tss = False
    is_flank = False
    while left >= 0 and (pos - tss_list[left]) <= 2000:
        d = abs(pos - tss_list[left])
        if d <= 100:
            is_tss = True
        elif 1900 <= d <= 2000:
            is_flank = True
        left -= 1
    n = len(tss_list)
    while right < n and (tss_list[right] - pos) <= 2000:
        d = abs(pos - tss_list[right])
        if d <= 100:
            is_tss = True
        elif 1900 <= d <= 2000:
            is_flank = True
        right += 1
    return (is_tss, is_flank)

tss_by_chrom = load_tss_positions(gtf_path)
counts = {}
tss_counts = {}
flank_counts = {}
for raw in sys.stdin:
    cols = raw.rstrip("\\n").split("\\t")
    if len(cols) < 12:
        continue
    chrom = cols[2]
    if chrom == "*" or chrom not in tss_by_chrom:
        continue
    try:
        pos = int(cols[3])
    except Exception:
        continue
    cb = None
    for tag in cols[11:]:
        if tag.startswith("CB:Z:"):
            cb = tag[5:]
            break
    if not cb:
        continue
    if valid and cb not in valid:
        continue
    counts[cb] = counts.get(cb, 0) + 1
    is_tss, is_flank = classify_tss_hit(tss_by_chrom[chrom], pos)
    if is_tss:
        tss_counts[cb] = tss_counts.get(cb, 0) + 1
    if is_flank:
        flank_counts[cb] = flank_counts.get(cb, 0) + 1

rows = []
for bc, n in counts.items():
    tss = tss_counts.get(bc, 0)
    flank = flank_counts.get(bc, 0)
    tss_ratio = (tss + 1.0) / (flank + 1.0)
    pass_frags = n >= min_frags
    pass_tss = tss_ratio >= min_tss
    rows.append((bc, n, tss, flank, tss_ratio, pass_frags, pass_tss))

passing_both = [r for r in rows if r[5] and r[6]]
passing_frags_only = [r for r in rows if r[5]]
with open(out_counts, "w") as out:
    out.write("Barcode\\tFragments\\tTSSReads\\tFlankReads\\tTSSRatio\\tPassMinFrags\\tPassMinTSS\\tPassCell\\n")
    for bc, n, tss, flank, tss_ratio, pass_frags, pass_tss in sorted(rows, key=lambda x: x[1], reverse=True):
        out.write(f"{bc}\\t{n}\\t{tss}\\t{flank}\\t{tss_ratio:.6f}\\t{1 if pass_frags else 0}\\t{1 if pass_tss else 0}\\t{1 if (pass_frags and pass_tss) else 0}\\n")

with open(out_summary, "w") as out:
    out.write("Metric\\tValue\\n")
    out.write(f"Sample\\t{sample_id}\\n")
    out.write(f"MinFragmentsForCell\\t{min_frags}\\n")
    out.write(f"MinTSSRatioForCell\\t{min_tss}\\n")
    out.write(f"BarcodesWithFragments\\t{len(counts)}\\n")
    out.write(f"EstimatedCellsMinFragsOnly\\t{len(passing_frags_only)}\\n")
    out.write(f"EstimatedCells\\t{len(passing_both)}\\n")
    out.write(f"MedianFragmentsPerBarcode\\t{statistics.median(counts.values()) if counts else 0}\\n")
    out.write(f"MedianFragmentsPerEstimatedCell\\t{statistics.median([r[1] for r in passing_both]) if passing_both else 0}\\n")
PY
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
    tuple val(report_trigger), path(knee_plot_files)

    output:
    path "QC_Report.html", emit: qc_html
    path "QC_Report_assets/*", emit: qc_assets
    path "QC_Report", emit: qc_report_dir
    path "QC_Report_bundle.zip", emit: qc_bundle

    """
    python3 - "${projectDir}" ${knee_plot_files} <<'PY'
import glob
import html
import os
import csv
import sys
import shutil
import zipfile
import re
import datetime
import statistics

proj = sys.argv[1]
knee_plot_inputs = [p for p in sys.argv[2:] if os.path.isfile(p)]
species_model = "${params.species_model}"
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

def rel_list_recursive(pattern):
    return sorted([
        os.path.relpath(p, proj)
        for p in glob.glob(os.path.join(proj, pattern), recursive=True)
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
    delim = "," if rel_path.lower().endswith(".csv") else chr(9)
    try:
        with open(abs_path, newline="") as fh:
            reader = csv.reader(fh, delimiter=delim)
            for i, row in enumerate(reader):
                rows.append([html.escape(x) for x in row])
                if max_rows is not None and i + 1 >= max_rows:
                    break
    except Exception as e:
        return f"<p><em>Could not parse table: {html.escape(str(e))}</em></p>"
    if not rows:
        return "<p><em>File is empty</em></p>"
    cells = []
    for ridx, row in enumerate(rows):
        tag = "th" if ridx == 0 else "td"
        cells.append("<tr>" + "".join(f"<{tag}>{c}</{tag}>" for c in row) + "</tr>")
    return "<table>" + "".join(cells) + "</table>"

def load_demux_sample_names(rel_paths):
    names = set()
    delim_tsv = chr(9)
    for rel in rel_paths:
        p = os.path.join(proj, rel)
        if not os.path.isfile(p):
            continue
        delim = "," if rel.lower().endswith(".csv") else delim_tsv
        try:
            with open(p, newline="") as fh:
                reader = csv.reader(fh, delimiter=delim)
                rows = list(reader)
        except Exception:
            continue
        if not rows:
            continue
        header = rows[0]
        try:
            name_idx = header.index("Sample_Name")
        except ValueError:
            continue
        for r in rows[1:]:
            if len(r) > name_idx and r[name_idx].strip():
                names.add(r[name_idx].strip())
    return names

def demux_stats_html_for_sample(rel_paths, sample):
    delim_tsv = chr(9)
    for rel in rel_paths:
        p = os.path.join(proj, rel)
        if not os.path.isfile(p):
            continue
        delim = "," if rel.lower().endswith(".csv") else delim_tsv
        try:
            with open(p, newline="") as fh:
                reader = csv.reader(fh, delimiter=delim)
                all_rows = list(reader)
        except Exception:
            continue
        if not all_rows:
            continue
        header = all_rows[0]
        try:
            name_idx = header.index("Sample_Name")
        except ValueError:
            continue
        body = [r for r in all_rows[1:] if len(r) > name_idx and r[name_idx] == sample]
        if not body:
            continue
        esc_rows = [[html.escape(x) for x in row] for row in [header] + body]
        cells = []
        for ridx, row in enumerate(esc_rows):
            tag = "th" if ridx == 0 else "td"
            cells.append("<tr>" + "".join(f"<{tag}>{c}</{tag}>" for c in row) + "</tr>")
        return "<table>" + "".join(cells) + "</table>"
    return f"<p><em>No demultiplex stats row for <code>{html.escape(sample)}</code>.</em></p>"

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

def image_gallery_block(title, paths):
    if not paths:
        return f"<h3>{html.escape(title)}</h3><p><em>No files found.</em></p>"
    chunks = [f"<h3>{html.escape(title)}</h3>", '<div class="img-grid">']
    for p in paths:
        asset = stage_asset(p)
        if asset is None:
            continue
        chunks.append(
            '<figure class="img-card">'
            f'<img src="{html.escape(asset)}" alt="{html.escape(p)}" />'
            f'<figcaption>{html.escape(p)}</figcaption>'
            '</figure>'
        )
    chunks.append("</div>")
    return "".join(chunks)

def sample_from_report_path(rel_path):
    bits = rel_path.split("/")
    if len(bits) >= 3 and bits[0] in ("STARsolo", "STARsolo_paired"):
        return bits[1]
    if len(bits) >= 3 and bits[0] == "ATAC":
        return bits[1]
    return None

def path_matches_sample(rel_path, sample):
    if not rel_path or not sample:
        return False
    bits = rel_path.split("/")
    if sample in bits:
        return True
    base = os.path.basename(rel_path)
    patterns = [
        f"{sample}.",
        f"{sample}_",
        f"{sample}-",
        f"{sample}R1",
        f"{sample}R2",
    ]
    return any(tok in base for tok in patterns)

def parse_starsolo_log_metrics(rel_path):
    wanted = [
        "Number of input reads",
        "Uniquely mapped reads %",
        "% of reads mapped to multiple loci",
        "% of reads unmapped: other",
        "% of reads unmapped: too short",
    ]
    out = {k: "" for k in wanted}
    abs_path = os.path.join(proj, rel_path)
    if not os.path.isfile(abs_path):
        return out
    try:
        with open(abs_path, "r", errors="replace") as fh:
            for line in fh:
                if "|" not in line:
                    continue
                left, right = line.split("|", 1)
                key = left.strip()
                val = right.strip()
                if key in out:
                    out[key] = val
    except Exception:
        return out
    return out

def starsolo_summary_table(log_paths):
    if not log_paths:
        return "<p><em>No STARsolo Log.final.out files found.</em></p>"
    headers = [
        "Sample",
        "Number of input reads",
        "Uniquely mapped reads %",
        "% of reads mapped to multiple loci",
        "% of reads unmapped: other",
        "% of reads unmapped: too short",
    ]
    rows = []
    for p in sorted(log_paths):
        bits = p.split("/")
        if len(bits) >= 3 and bits[0] in ("STARsolo", "STARsolo_paired"):
            sample = bits[1]
        else:
            sample = os.path.basename(os.path.dirname(p))
        m = parse_starsolo_log_metrics(p)
        rows.append([
            sample,
            m["Number of input reads"],
            m["Uniquely mapped reads %"],
            m["% of reads mapped to multiple loci"],
            m["% of reads unmapped: other"],
            m["% of reads unmapped: too short"],
        ])
    cells = []
    cells.append("<tr>" + "".join(f"<th>{html.escape(h)}</th>" for h in headers) + "</tr>")
    for row in rows:
        cells.append("<tr>" + "".join(f"<td>{html.escape(str(c))}</td>" for c in row) + "</tr>")
    return "<table>" + "".join(cells) + "</table>"

def _parse_pct(value):
    if value is None:
        return 0.0
    s = str(value).strip().replace("%", "").replace(",", "")
    try:
        return float(s)
    except Exception:
        return 0.0

def _parse_number(value):
    if value is None:
        return None
    s = str(value).strip().replace(",", "").replace("%", "")
    if s == "":
        return None
    try:
        return float(s)
    except Exception:
        return None

def _fmt_int(v):
    try:
        return f"{int(round(v)):,}"
    except Exception:
        return "N/A"

def _fmt_float(v, nd=1, suffix=""):
    try:
        return f"{float(v):.{nd}f}{suffix}"
    except Exception:
        return "N/A"

def starsolo_metrics_by_sample(log_paths):
    out = {}
    for p in sorted(log_paths):
        sample = sample_from_report_path(p) or os.path.basename(os.path.dirname(p))
        out[sample] = parse_starsolo_log_metrics(p)
    return out

def summary_metrics_by_sample(paths):
    out = {}
    for p in sorted(paths):
        sample = sample_from_report_path(p) or os.path.basename(os.path.dirname(p))
        out[sample] = parse_summary_csv_metrics(p)
    return out

def overview_cards_html(sample_names, starsolo_by_sample, summary_by_sample):
    unique_vals = []
    med_reads_cell_vals = []
    est_cells_vals = []
    for s in sample_names:
        sm = starsolo_by_sample.get(s, {})
        unique = _parse_pct(sm.get("Uniquely mapped reads %", ""))
        if unique > 0:
            unique_vals.append(unique)
        sy = summary_by_sample.get(s, {})
        mrc = _parse_number(sy.get("Median Reads per Cell", ""))
        if mrc is not None:
            med_reads_cell_vals.append(mrc)
        est = _parse_number(sy.get("Estimated Number of Cells", ""))
        if est is not None:
            est_cells_vals.append(est)
    n_samples = len(sample_names)
    med_unique = statistics.median(unique_vals) if unique_vals else None
    med_reads_cell = statistics.median(med_reads_cell_vals) if med_reads_cell_vals else None
    total_est_cells = sum(est_cells_vals) if est_cells_vals else None
    generated = datetime.datetime.now().strftime("%Y-%m-%d %H:%M")

    cards = []
    cards.append('<div class="kpi-grid">')
    cards.append(f'<div class="kpi-card"><div class="kpi-label">Samples</div><div class="kpi-value">{_fmt_int(n_samples)}</div></div>')
    cards.append(f'<div class="kpi-card"><div class="kpi-label">Median Unique Mapping</div><div class="kpi-value">{_fmt_float(med_unique, 1, "%")}</div></div>')
    cards.append(f'<div class="kpi-card"><div class="kpi-label">Median Reads per Cell</div><div class="kpi-value">{_fmt_int(med_reads_cell)}</div></div>')
    cards.append(f'<div class="kpi-card"><div class="kpi-label">Total Estimated Cells</div><div class="kpi-value">{_fmt_int(total_est_cells)}</div></div>')
    cards.append('</div>')
    cards.append(f'<p class="meta-line"><strong>Generated:</strong> {html.escape(generated)} | <strong>Project:</strong> <code>{html.escape(proj)}</code></p>')
    return "".join(cards)

def sample_directory_table(sample_names, starsolo_by_sample, summary_by_sample):
    if not sample_names:
        return "<p><em>No samples detected.</em></p>"
    rows = []
    rows.append("<tr><th>Sample</th><th>Uniquely mapped reads %</th><th>Median Reads per Cell</th><th>Estimated Number of Cells</th><th>Sample report</th></tr>")
    for s in sample_names:
        sm = starsolo_by_sample.get(s, {})
        sy = summary_by_sample.get(s, {})
        unique = sm.get("Uniquely mapped reads %", "")
        mrc = sy.get("Median Reads per Cell", "")
        enc = sy.get("Estimated Number of Cells", "")
        link = f'QC_Report/{s}/index.html'
        rows.append(
            "<tr>"
            f"<td>{html.escape(s)}</td>"
            f"<td>{html.escape(unique)}</td>"
            f"<td>{html.escape(mrc)}</td>"
            f"<td>{html.escape(enc)}</td>"
            f'<td><a href="{html.escape(link)}">Open</a></td>'
            "</tr>"
        )
    return "<table>" + "".join(rows) + "</table>"

def sample_sidebar_links(sample_names):
    if not sample_names:
        return "<p><em>No sample reports</em></p>"
    chunks = ['<div class="sample-side-title">Sample Reports</div>', '<div class="sample-side-links">']
    for s in sample_names:
        link = f'./QC_Report/{s}/index.html'
        chunks.append(f'<a href="{html.escape(link)}">{html.escape(s)}</a>')
    chunks.append('</div>')
    return "".join(chunks)

def alignment_summary_chart(log_paths):
    if not log_paths:
        return "<p><em>No STARsolo Log.final.out files found.</em></p>"
    rows = []
    for p in sorted(log_paths):
        bits = p.split("/")
        if len(bits) >= 3 and bits[0] in ("STARsolo", "STARsolo_paired"):
            sample = bits[1]
        else:
            sample = os.path.basename(os.path.dirname(p))
        m = parse_starsolo_log_metrics(p)
        unique = _parse_pct(m.get("Uniquely mapped reads %", ""))
        multi = _parse_pct(m.get("% of reads mapped to multiple loci", ""))
        unmapped = max(0.0, 100.0 - unique - multi)
        total = unique + multi + unmapped
        if total <= 0:
            unique_w = multi_w = unmapped_w = 0.0
        else:
            unique_w = (unique / total) * 100.0
            multi_w = (multi / total) * 100.0
            unmapped_w = (unmapped / total) * 100.0
        rows.append((sample, unique, multi, unmapped, unique_w, multi_w, unmapped_w))

    chunks = []
    chunks.append('<div class="align-legend">')
    chunks.append('<span><i class="swatch unique"></i>Uniquely Mapped</span>')
    chunks.append('<span><i class="swatch multi"></i>Multi-mapped</span>')
    chunks.append('<span><i class="swatch unmapped"></i>Unmapped</span>')
    chunks.append('</div>')
    chunks.append('<div class="align-chart">')
    for sample, unique, multi, unmapped, unique_w, multi_w, unmapped_w in rows:
        chunks.append('<div class="align-row">')
        chunks.append(f'<div class="align-sample">{html.escape(sample)}</div>')
        chunks.append('<div class="align-bar">')
        chunks.append(f'<span class="seg unique" style="width:{unique_w:.4f}%"></span>')
        chunks.append(f'<span class="seg multi" style="width:{multi_w:.4f}%"></span>')
        chunks.append(f'<span class="seg unmapped" style="width:{unmapped_w:.4f}%"></span>')
        chunks.append('</div>')
        chunks.append(
            f'<div class="align-values">U: {unique:.2f}% | M: {multi:.2f}% | Un: {unmapped:.2f}%</div>'
        )
        chunks.append('</div>')
    chunks.append('</div>')
    return "".join(chunks)

def parse_barcodes_stats_metrics(rel_path):
    metrics = {}
    abs_path = os.path.join(proj, rel_path)
    if not os.path.isfile(abs_path):
        return metrics
    try:
        with open(abs_path, "r", errors="replace") as fh:
            for raw in fh:
                line = raw.strip()
                if not line:
                    continue
                key = None
                val = None
                if "\t" in line:
                    parts = [x.strip() for x in line.split("\t") if x.strip()]
                    if len(parts) >= 2:
                        key = parts[0]
                        val = parts[-1]
                elif "|" in line:
                    left, right = line.split("|", 1)
                    key = left.strip()
                    val = right.strip()
                elif ":" in line:
                    left, right = line.split(":", 1)
                    key = left.strip()
                    val = right.strip()
                else:
                    parts = re.split(r"\s{2,}", line)
                    if len(parts) >= 2:
                        key = parts[0].strip()
                        val = parts[-1].strip()
                if key and val:
                    metrics[key] = val
    except Exception:
        return {}
    return metrics

def barcodes_stats_summary_table(paths):
    if not paths:
        return "<p><em>No Barcodes.stats files found.</em></p>"
    sample_metrics = []
    metric_order = []
    metric_seen = set()
    for p in sorted(paths):
        bits = p.split("/")
        if len(bits) >= 3 and bits[0] in ("STARsolo", "STARsolo_paired"):
            sample = bits[1]
        else:
            sample = os.path.basename(os.path.dirname(p))
        metrics = parse_barcodes_stats_metrics(p)
        sample_metrics.append((sample, metrics))
        for k in metrics.keys():
            if k not in metric_seen:
                metric_seen.add(k)
                metric_order.append(k)
    sample_names = [s for s, _ in sample_metrics]
    headers = ["Metric"] + sample_names
    cells = []
    cells.append("<tr>" + "".join(f"<th>{html.escape(h)}</th>" for h in headers) + "</tr>")
    for metric in metric_order:
        row = [metric] + [metrics.get(metric, "") for _, metrics in sample_metrics]
        cells.append("<tr>" + "".join(f"<td>{html.escape(str(c))}</td>" for c in row) + "</tr>")
    return "<table>" + "".join(cells) + "</table>"

def parse_summary_csv_metrics(rel_path):
    metrics = {}
    abs_path = os.path.join(proj, rel_path)
    if not os.path.isfile(abs_path):
        return metrics
    try:
        with open(abs_path, newline="") as fh:
            reader = csv.reader(fh)
            for row in reader:
                if not row:
                    continue
                key = row[0].strip()
                val = (row[1].strip() if len(row) > 1 else "")
                if key:
                    metrics[key] = val
    except Exception:
        return {}
    return metrics

def summary_csv_key_metrics_table(paths):
    wanted_metrics = [
        "Number of Reads",
        "Reads With Valid Barcodes",
        "Reads Mapped to Genome: Unique",
        "Estimated Number of Cells",
        "Median Reads per Cell",
        "Median UMI per Cell",
        "Total GeneFull Detected",
    ]
    if not paths:
        return "<p><em>No Summary.csv files found.</em></p>"
    sample_metrics = []
    for p in sorted(paths):
        bits = p.split("/")
        if len(bits) >= 3 and bits[0] in ("STARsolo", "STARsolo_paired"):
            sample = bits[1]
        else:
            sample = os.path.basename(os.path.dirname(p))
        metrics = parse_summary_csv_metrics(p)
        sample_metrics.append((sample, metrics))
    headers = ["Metric"] + [s for s, _ in sample_metrics]
    cells = []
    cells.append("<tr>" + "".join(f"<th>{html.escape(h)}</th>" for h in headers) + "</tr>")
    for metric in wanted_metrics:
        row_html = [f"<td>{html.escape(metric)}</td>"]
        for _, m in sample_metrics:
            val = m.get(metric, "")
            if metric == "Reads Mapped to Genome: Unique":
                pct = _parse_pct(val)
                if pct >= 85.0:
                    cls = "qc-good"
                elif pct >= 70.0:
                    cls = "qc-warn"
                else:
                    cls = "qc-bad"
                row_html.append(f'<td class="{cls}">{html.escape(str(val))}</td>')
            else:
                row_html.append(f"<td>{html.escape(str(val))}</td>")
        cells.append("<tr>" + "".join(row_html) + "</tr>")
    return "<table>" + "".join(cells) + "</table>"

def parse_flagstat_metrics(rel_path):
    out = {
        "in_total": "",
        "mapped_pct": "",
        "properly_paired_pct": "",
        "duplicates_pct": "",
        "singletons_pct": "",
        "mate_diff_chr_count": "",
    }
    abs_path = os.path.join(proj, rel_path)
    if not os.path.isfile(abs_path):
        return out
    try:
        with open(abs_path, "r", errors="replace") as fh:
            for line in fh:
                s = line.strip()
                if " in total" in s and "+" in s:
                    out["in_total"] = s.split("+", 1)[0].strip()
                elif " mapped (" in s:
                    m = re.search(r'\\(([^)]*%)', s)
                    out["mapped_pct"] = m.group(1) if m else ""
                elif " properly paired (" in s:
                    m = re.search(r'\\(([^)]*%)', s)
                    out["properly_paired_pct"] = m.group(1) if m else ""
                elif " duplicates" in s and "(" in s:
                    m = re.search(r'\\(([^)]*%)', s)
                    out["duplicates_pct"] = m.group(1) if m else ""
                elif " singletons (" in s:
                    m = re.search(r'\\(([^)]*%)', s)
                    out["singletons_pct"] = m.group(1) if m else ""
                elif " with mate mapped to a different chr" in s and "(" not in s and "+" in s:
                    out["mate_diff_chr_count"] = s.split("+", 1)[0].strip()
    except Exception:
        return out
    return out

def parse_idxstats_metrics(rel_path):
    out = {
        "total_mapped": 0,
        "mt_mapped": 0,
        "autosomal_mapped": 0,
        "x_mapped": 0,
        "y_mapped": 0,
    }
    abs_path = os.path.join(proj, rel_path)
    if not os.path.isfile(abs_path):
        return out
    try:
        with open(abs_path, "r", errors="replace") as fh:
            for line in fh:
                cols = line.rstrip("\\n").split("\\t")
                if len(cols) < 3:
                    continue
                chrom = cols[0].strip()
                try:
                    mapped = int(cols[2])
                except Exception:
                    continue
                out["total_mapped"] += mapped
                c = chrom.upper()
                if c in ("MT", "CHRM"):
                    out["mt_mapped"] += mapped
                elif c == "X" or c == "CHRX":
                    out["x_mapped"] += mapped
                elif c == "Y" or c == "CHRY":
                    out["y_mapped"] += mapped
                elif re.fullmatch(r'(CHR)?([1-9]|1[0-9]|2[0-2])', c):
                    out["autosomal_mapped"] += mapped
    except Exception:
        return out
    return out

def parse_atac_cell_summary(rel_path):
    out = {}
    abs_path = os.path.join(proj, rel_path)
    if not os.path.isfile(abs_path):
        return out
    try:
        with open(abs_path, "r", errors="replace") as fh:
            for i, line in enumerate(fh):
                if i == 0:
                    continue
                cols = line.rstrip("\n").split("\t", 1)
                if len(cols) == 2:
                    out[cols[0].strip()] = cols[1].strip()
    except Exception:
        return out
    return out

def sample_name_matches(sample, candidate):
    if not sample or not candidate:
        return False
    s = str(sample).strip()
    c = str(candidate).strip()
    if c == s:
        return True
    s_l = s.lower()
    c_l = c.lower()
    if c_l == s_l:
        return True

    # MultiQC often stores full paths / filenames; compare against basename too.
    c_base = os.path.basename(c_l)

    def _norm(x):
        y = x
        y = y.replace(".fastq.gz", "").replace(".fq.gz", "")
        y = y.replace(".bam", "").replace(".txt", "")
        y = y.replace(".q30.rmdup.sorted", "").replace(".q30.rmdup", "")
        y = y.replace(".q30.possort", "").replace(".q30.mapped", "")
        y = y.replace(".flagstat", "").replace(".idxstats", "").replace(".stats", "")
        return y

    s_n = _norm(s_l)
    c_n = _norm(c_l)
    c_base_n = _norm(c_base)

    # Allow exact and containment matches after normalization.
    if s_n == c_n or s_n == c_base_n:
        return True
    if s_n in c_n or s_n in c_base_n:
        return True
    if c_n in s_n or c_base_n in s_n:
        return True
    return False

def multiqc_insert_size_for_sample(rel_paths, sample):
    pre_insert = ""
    post_insert = ""
    for rel in rel_paths:
        abs_path = os.path.join(proj, rel)
        if not os.path.isfile(abs_path):
            continue
        try:
            with open(abs_path, newline="") as fh:
                reader = csv.reader(fh, delimiter="\t")
                all_rows = list(reader)
        except Exception:
            continue
        if len(all_rows) < 2:
            continue
        header = all_rows[0]
        try:
            ins_idx = header.index("samtools_stats-insert_size_average")
        except ValueError:
            continue
        for r in all_rows[1:]:
            if not r:
                continue
            name = r[0] if len(r) > 0 else ""
            if sample_name_matches(sample, name):
                ins = (r[ins_idx] if ins_idx < len(r) else "").strip()
                if not ins:
                    continue
                name_l = str(name).lower()
                if ".q30.rmdup" in name_l:
                    post_insert = ins
                elif ".q30" in name_l:
                    pre_insert = ins
    return [pre_insert, post_insert]

def _fmt_pct(v):
    try:
        return f"{float(v):.2f}%"
    except Exception:
        return ""

def _to_int(s):
    try:
        return int(str(s).replace(",", "").strip())
    except Exception:
        return None

def atac_alignment_summary_table(flagstat_paths, idxstats_paths):
    if not flagstat_paths and not idxstats_paths:
        return "<p><em>No ATAC alignment metrics found.</em></p>"
    flag_by_sample = { (sample_from_report_path(p) or os.path.basename(os.path.dirname(p))): p for p in flagstat_paths }
    idx_by_sample = { (sample_from_report_path(p) or os.path.basename(os.path.dirname(p))): p for p in idxstats_paths }
    samples = sorted(set(list(flag_by_sample.keys()) + list(idx_by_sample.keys())))
    rows = []
    for sample in samples:
        f = parse_flagstat_metrics(flag_by_sample[sample]) if sample in flag_by_sample else {}
        i = parse_idxstats_metrics(idx_by_sample[sample]) if sample in idx_by_sample else {}
        total_mapped = i.get("total_mapped", 0)
        mt_frac = (100.0 * i.get("mt_mapped", 0) / total_mapped) if total_mapped > 0 else 0.0
        rows.append([
            sample,
            f.get("in_total", ""),
            f.get("properly_paired_pct", ""),
            f.get("singletons_pct", ""),
            f.get("mate_diff_chr_count", ""),
            _fmt_pct(mt_frac),
            f"{i.get('autosomal_mapped', 0):,}",
            f"{i.get('x_mapped', 0):,}",
            f"{i.get('y_mapped', 0):,}",
        ])
    cells = []
    cells.append(
        "<tr>"
        "<th>Sample</th>"
        "<th>Total reads (flagstat)</th>"
        "<th>Properly paired %</th>"
        "<th>Singleton %</th>"
        "<th>Mate on different chr</th>"
        "<th>Mitochondrial fraction</th>"
        "<th>Autosomal mapped</th>"
        "<th>X mapped</th>"
        "<th>Y mapped</th>"
        "</tr>"
    )
    for row in rows:
        cells.append("<tr>" + "".join(f"<td>{html.escape(str(c))}</td>" for c in row) + "</tr>")
    return "<table>" + "".join(cells) + "</table>"

def atac_key_summary_table(flagstat_prededup_paths, idxstats_prededup_paths, flagstat_rmdup_paths, idxstats_rmdup_paths, atac_cell_summary_paths):
    samples = sorted(set(
        [sample_from_report_path(p) or os.path.basename(os.path.dirname(p)) for p in (flagstat_prededup_paths + idxstats_prededup_paths + flagstat_rmdup_paths + idxstats_rmdup_paths)]
    ))
    if not samples:
        return "<p><em>No ATAC alignment metrics found.</em></p>"

    pre_flag = { (sample_from_report_path(p) or os.path.basename(os.path.dirname(p))): p for p in flagstat_prededup_paths }
    pre_idx = { (sample_from_report_path(p) or os.path.basename(os.path.dirname(p))): p for p in idxstats_prededup_paths }
    post_flag = { (sample_from_report_path(p) or os.path.basename(os.path.dirname(p))): p for p in flagstat_rmdup_paths }
    post_idx = { (sample_from_report_path(p) or os.path.basename(os.path.dirname(p))): p for p in idxstats_rmdup_paths }
    cell_sum = { (sample_from_report_path(p) or os.path.basename(os.path.dirname(p))): p for p in atac_cell_summary_paths }

    rows = []
    for s in samples:
        pf = parse_flagstat_metrics(pre_flag[s]) if s in pre_flag else {}
        pi = parse_idxstats_metrics(pre_idx[s]) if s in pre_idx else {}
        rf = parse_flagstat_metrics(post_flag[s]) if s in post_flag else {}
        ri = parse_idxstats_metrics(post_idx[s]) if s in post_idx else {}
        cs = parse_atac_cell_summary(cell_sum[s]) if s in cell_sum else {}

        pre_total = _to_int(pf.get("in_total", ""))
        post_total = _to_int(rf.get("in_total", ""))
        retained_pct = (100.0 * post_total / pre_total) if (pre_total and post_total is not None and pre_total > 0) else None
        removed_pct = (100.0 - retained_pct) if retained_pct is not None else None

        pre_mt = (100.0 * pi.get("mt_mapped", 0) / pi.get("total_mapped", 1)) if pi.get("total_mapped", 0) > 0 else 0.0
        post_mt = (100.0 * ri.get("mt_mapped", 0) / ri.get("total_mapped", 1)) if ri.get("total_mapped", 0) > 0 else 0.0

        rows.append([
            s,
            cs.get("EstimatedCells", ""),
            pf.get("in_total", ""),
            pf.get("properly_paired_pct", ""),
            _fmt_pct(pre_mt),
            rf.get("in_total", ""),
            rf.get("properly_paired_pct", ""),
            _fmt_pct(post_mt),
            _fmt_pct(retained_pct) if retained_pct is not None else "",
            _fmt_pct(removed_pct) if removed_pct is not None else "",
        ])

    cells = []
    cells.append(
        "<tr>"
        "<th>Sample</th>"
        "<th>Estimated cells</th>"
        "<th>Reads pre-dedup</th>"
        "<th>Properly paired % pre</th>"
        "<th>MT % pre</th>"
        "<th>Reads post-dedup</th>"
        "<th>Properly paired % post</th>"
        "<th>MT % post</th>"
        "<th>Reads retained</th>"
        "<th>Reads removed</th>"
        "</tr>"
    )
    for row in rows:
        cells.append("<tr>" + "".join(f"<td>{html.escape(str(c))}</td>" for c in row) + "</tr>")
    return "<table>" + "".join(cells) + "</table>"

def atac_sample_metric_cards(sample, flagstat_paths, idxstats_paths):
    flag = parse_flagstat_metrics(flagstat_paths[0]) if flagstat_paths else {}
    idx = parse_idxstats_metrics(idxstats_paths[0]) if idxstats_paths else {}
    total_mapped = idx.get("total_mapped", 0)
    mt_frac = (100.0 * idx.get("mt_mapped", 0) / total_mapped) if total_mapped > 0 else 0.0
    cards = []
    cards.append('<div class="kpi-grid">')
    cards.append(f'<div class="kpi-card"><div class="kpi-label">Sample</div><div class="kpi-value">{html.escape(sample)}</div></div>')
    cards.append(f'<div class="kpi-card"><div class="kpi-label">Total reads</div><div class="kpi-value">{html.escape(flag.get("in_total", "N/A"))}</div></div>')
    cards.append(f'<div class="kpi-card"><div class="kpi-label">Properly paired</div><div class="kpi-value">{html.escape(flag.get("properly_paired_pct", "N/A"))}</div></div>')
    cards.append(f'<div class="kpi-card"><div class="kpi-label">Singletons</div><div class="kpi-value">{html.escape(flag.get("singletons_pct", "N/A"))}</div></div>')
    cards.append(f'<div class="kpi-card"><div class="kpi-label">Mate diff chr</div><div class="kpi-value">{html.escape(flag.get("mate_diff_chr_count", "N/A"))}</div></div>')
    cards.append(f'<div class="kpi-card"><div class="kpi-label">Mitochondrial fraction</div><div class="kpi-value">{_fmt_pct(mt_frac)}</div></div>')
    cards.append('</div>')
    return "".join(cards)

demux_total = rel_list("demux/*.total_number_reads.tsv")
demux_stats = sorted(set(
    rel_list("demux/SHARE-seq.demultiplex.stats.tsv")
    + rel_list_recursive("demux/**/SHARE-seq.demultiplex.stats.tsv")
))
fastqc_html = sorted(set(
    rel_list("fastqc_demux/*/*_fastqc.html") +
    rel_list("fastqc_demux/*/fastqc_demux/*_fastqc.html") +
    rel_list("fastqc_demux/*_fastqc.html")
))
fastqc_trimmed_html = sorted(set(
    rel_list("fastqc_trimmed/*_fastqc.html") +
    rel_list_recursive("fastqc_trimmed/**/*_fastqc.html")
))

starsolo_logs = rel_list("STARsolo/*/Log.final.out") + rel_list("STARsolo_paired/*/Log.final.out")
knee_plots = sorted(set(
    rel_list("STARsolo/*/*_knee_plot.png") +
    rel_list("STARsolo_paired/*/*_knee_plot.png") +
    rel_list_recursive("STARsolo/**/*_knee_plot.png") +
    rel_list_recursive("STARsolo_paired/**/*_knee_plot.png")
))
barcodes_stats = rel_list("STARsolo/*/Solo.out/Barcodes.stats") + rel_list("STARsolo_paired/*/Solo.out/Barcodes.stats")
summary_csv = rel_list("STARsolo/*/Solo.out/GeneFull/Summary.csv") + rel_list("STARsolo_paired/*/Solo.out/GeneFull/Summary.csv")
barnyard = rel_list("STARsolo/*/*collision_plot.png") + rel_list("STARsolo_paired/*/*collision_plot.png")
atac_flagstat_rmdup = sorted(set(
    rel_list("ATAC/*/*.q30.rmdup.flagstat.txt") +
    rel_list("ATAC/*/*.flagstat.txt")
))
atac_idxstats_rmdup = sorted(set(
    rel_list("ATAC/*/*.q30.rmdup.idxstats.txt") +
    rel_list("ATAC/*/*.idxstats.txt")
))
atac_stats_rmdup = sorted(set(
    rel_list("ATAC/*/*.q30.rmdup.stats.txt") +
    rel_list("ATAC/*/*.stats.txt")
))
atac_flagstat_prededup = rel_list("ATAC/*/*.q30.mapped.flagstat.txt")
atac_idxstats_prededup = rel_list("ATAC/*/*.q30.mapped.idxstats.txt")
atac_stats_prededup = rel_list("ATAC/*/*.q30.mapped.stats.txt")
atac_cell_summary = rel_list("ATAC/*/*.atac_cells.summary.tsv")
atac_multiqc_report = sorted(set(
    rel_list("multiqc_atac/ATAC_MultiQC.html") +
    rel_list("ATAC_MultiQC.html")
))
atac_multiqc_tables = sorted(set(
    rel_list("multiqc_atac/ATAC_MultiQC_data/multiqc_general_stats.txt") +
    rel_list("multiqc_atac/ATAC_MultiQC_data/multiqc_samtools_flagstat.txt") +
    rel_list("multiqc_atac/ATAC_MultiQC_data/multiqc_samtools_stats.txt") +
    rel_list("ATAC_MultiQC_data/multiqc_general_stats.txt") +
    rel_list("ATAC_MultiQC_data/multiqc_samtools_flagstat.txt") +
    rel_list("ATAC_MultiQC_data/multiqc_samtools_stats.txt")
))
sample_names = sorted(set(
    [sample_from_report_path(p) for p in (starsolo_logs + barcodes_stats + summary_csv + knee_plots + barnyard + atac_flagstat_rmdup + atac_idxstats_rmdup + atac_stats_rmdup + atac_flagstat_prededup + atac_idxstats_prededup + atac_stats_prededup) if sample_from_report_path(p)]
    + list(load_demux_sample_names(demux_stats))
))
starsolo_by_sample = starsolo_metrics_by_sample(starsolo_logs)
summary_by_sample = summary_metrics_by_sample(summary_csv)

parts = []
parts.append('''<!doctype html>
<html>
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>SHARE-seq QC Report</title>
  <style>
    body {
      font-family: "Inter", "Segoe UI", Roboto, "Helvetica Neue", Arial, sans-serif;
      font-size: 14px;
      color: #000000;
      margin: 24px auto;
      max-width: 1400px;
      line-height: 1.55;
      padding: 0 12px;
      background: #ffffff;
      scroll-padding-top: 86px;
    }
    h1, h2, h3, h4 {
      color: #000000;
      margin-top: 1.15em;
      margin-bottom: 0.45em;
      font-weight: 700;
      letter-spacing: 0.1px;
    }
    h1 { font-size: 30px; margin-top: 0.3em; }
    h2 { font-size: 21px; border-bottom: 1px solid #dfe1df; padding-bottom: 4px; }
    h3 { font-size: 16px; font-weight: 700; }
    h4 { font-size: 13px; font-weight: 600; color: #000000; }
    p { margin: 0.45em 0 0.8em 0; }
    em { color: #000000; }
    code { background: #ffffff; color: #000000; border: 1px solid #dfe1df; padding: 2px 5px; border-radius: 4px; font-size: 0.92em; }
    table {
      border-collapse: collapse;
      margin: 10px 0 20px 0;
      width: 100%;
      background: #fff;
      box-shadow: 0 1px 2px rgba(141, 0, 52, 0.08);
      display: table;
    }
    th, td {
      border: 1px solid #dfe1df;
      padding: 6px 10px;
      font-size: 12px;
      vertical-align: top;
      white-space: nowrap;
    }
    th {
      background: #ffffff;
      font-weight: 700;
      text-align: left;
      color: #000000;
    }
    td:first-child { font-weight: 600; color: #000000; }
    tr:nth-child(even) td { background: #ffffff; }
    .tabs {
      display: flex;
      flex-wrap: wrap;
      gap: 8px;
      gap: 8px;
      padding-left: 10px;
      padding-right: 8px;
    }
    .tabs a {
      text-decoration: none;
      border: 1px solid rgba(255, 255, 255, 0.6);
      background: rgba(255, 255, 255, 0.14);
      color: #ffffff;
      padding: 7px 12px;
      border-radius: 6px;
      font-size: 12px;
      font-weight: 600;
    }
    .tabs a:hover { background: rgba(255, 255, 255, 0.26); border-color: #ffffff; }
    .top-banner {
      background: #d11947;
      color: #ffffff;
      border-radius: 10px;
      padding: 12px 14px 12px 14px;
      margin: 0 0 14px 0;
      box-shadow: 0 2px 6px rgba(141, 0, 52, 0.22);
    }
    .top-banner h1, .top-banner h2, .top-banner h3, .top-banner h4,
    .top-banner p, .top-banner em, .top-banner code { color: #ffffff; }
    .top-banner code {
      background: rgba(255, 255, 255, 0.14);
      border: 1px solid rgba(255, 255, 255, 0.45);
    }
    section {
      padding: 10px 14px 6px 14px;
      margin: 8px 0 14px 0;
      background: #ffffff;
      color: #000000;
      border: 1px solid #dfe1df;
      border-radius: 10px;
      box-shadow: 0 1px 2px rgba(141, 0, 52, 0.06);
      scroll-margin-top: 90px;
    }
    section h2, section h3, section h4, section p, section em, section a { color: #000000; }
    section h2 { border-bottom: 1px solid #dfe1df; }
    section code {
      background: #ffffff;
      border: 1px solid #dfe1df;
      color: #000000;
    }
    section table, section th, section td {
      background: #ffffff;
      color: #000000;
    }
    .img-grid {
      display: grid;
      grid-template-columns: repeat(auto-fit, minmax(360px, 1fr));
      gap: 12px;
      margin: 8px 0 18px 0;
    }
    .img-card {
      margin: 0;
      border: 1px solid #dfe1df;
      border-radius: 6px;
      padding: 8px;
      background: #fff;
    }
    .img-card img {
      width: 100%;
      height: auto;
      max-height: 260px;
      object-fit: contain;
      display: block;
    }
    .img-card figcaption {
      font-size: 11px;
      margin-top: 6px;
      word-break: break-word;
    }
    .align-legend { display: flex; gap: 14px; flex-wrap: wrap; margin: 6px 0 10px 0; font-size: 12px; }
    .align-legend .swatch {
      display: inline-block;
      width: 10px;
      height: 10px;
      border-radius: 2px;
      margin-right: 5px;
      vertical-align: middle;
    }
    .align-chart { margin: 6px 0 16px 0; }
    .align-row {
      display: grid;
      grid-template-columns: 180px minmax(280px, 1fr) 260px;
      gap: 10px;
      align-items: center;
      margin-bottom: 8px;
    }
    .align-sample { font-size: 12px; word-break: break-word; }
    .align-bar {
      display: flex;
      height: 14px;
      border: 1px solid #d0d0d0;
      border-radius: 3px;
      overflow: hidden;
      background: #f5f5f5;
    }
    .align-bar .seg { height: 100%; display: block; }
    .align-values { font-size: 11px; color: #000000; font-weight: 600; }
    .unique { background: #2e7d32; }
    .multi { background: #1976d2; }
    .unmapped { background: #d32f2f; }
    .qc-good { background: #e6f4ea; color: #1b5e20; font-weight: 700; }
    .qc-warn { background: #fff8e1; color: #8a6d00; font-weight: 700; }
    .qc-bad { background: #fdecea; color: #b71c1c; font-weight: 700; }
    .kpi-grid {
      display: grid;
      grid-template-columns: repeat(auto-fit, minmax(210px, 1fr));
      gap: 10px;
      margin: 8px 0 10px 0;
    }
    .kpi-card {
      border: 1px solid #dfe1df;
      border-radius: 8px;
      background: #ffffff;
      padding: 10px 12px;
    }
    .kpi-label {
      font-size: 12px;
      font-weight: 600;
      color: #000000;
      margin-bottom: 2px;
    }
    .kpi-value {
      font-size: 24px;
      font-weight: 700;
      color: #000000;
      line-height: 1.2;
    }
    .meta-line { margin-top: 4px; color: #000000; }
    .main-layout {
      display: grid;
      grid-template-columns: 230px minmax(0, 1fr);
      gap: 14px;
      align-items: start;
    }
    .sample-sidebar {
      position: sticky;
      top: 74px;
      align-self: start;
      background: #d11947;
      border: 1px solid #8d0034;
      border-radius: 10px;
      padding: 10px;
      color: #ffffff;
      box-shadow: 0 2px 6px rgba(141, 0, 52, 0.22);
      max-height: calc(100vh - 74px);
      overflow-y: auto;
    }
    .tabs-wrap {
      position: sticky;
      top: 0;
      z-index: 1000;
      background: #d11947;
      border-radius: 0 0 8px 8px;
      padding: 8px 0 6px 0;
      margin-bottom: 14px;
      box-shadow: 0 2px 4px rgba(141, 0, 52, 0.18);
    }
    .sample-side-title {
      font-size: 12px;
      font-weight: 700;
      color: #ffffff;
      margin-bottom: 8px;
      text-transform: uppercase;
      letter-spacing: 0.3px;
    }
    .sample-side-links {
      display: flex;
      flex-direction: column;
      gap: 6px;
    }
    .sample-side-links a {
      text-decoration: none;
      font-size: 12px;
      color: #ffffff;
      border: 1px solid rgba(255, 255, 255, 0.6);
      border-radius: 6px;
      padding: 6px 8px;
      background: rgba(255, 255, 255, 0.14);
      font-weight: 700;
    }
    .sample-side-links a:hover { background: rgba(255, 255, 255, 0.26); border-color: #ffffff; }
    .main-content { min-width: 0; }
    .starsolo-block { margin: 0 0 18px 0; }
    @media (max-width: 1100px) {
      .main-layout { grid-template-columns: 1fr; }
      .sample-sidebar { position: static; max-height: none; }
    }
  </style>
</head>
<body>
<div class="top-banner">
<h1>SHARE-seq QC and Visualization Report</h1>
<p>Generated from pipeline outputs in <code>__PROJECT_DIR__</code>.</p>
</div>
<div class="tabs-wrap">
<div class="tabs">
  <a href="#sec-overview">Overview</a>
  <a href="#sec-demux">Demultiplexing</a>
  <a href="#sec-fastqc">FastQC (Demultiplexed)</a>
  <a href="#sec-atac">ATAC QC</a>
  <a href="#sec-starsolo">RNA QC</a>
__BARNYARD_TAB__
</div>
</div>
<div class="main-layout">
<aside class="sample-sidebar">
__SAMPLE_SIDEBAR__
</aside>
<div class="main-content">
'''.replace("__PROJECT_DIR__", html.escape(proj)))
parts[-1] = parts[-1].replace("__SAMPLE_SIDEBAR__", sample_sidebar_links(sample_names))
parts[-1] = parts[-1].replace(
    "__BARNYARD_TAB__",
    '  <a href="#sec-barnyard">Hybrid Barnyard</a>' if species_model == 'hybrid' else ""
)

parts.append('<section id="sec-overview">')
parts.append("<h2>Overview</h2>")
parts.append(overview_cards_html(sample_names, starsolo_by_sample, summary_by_sample))
parts.append("<h3>Sample Directory</h3>")
parts.append(sample_directory_table(sample_names, starsolo_by_sample, summary_by_sample))
parts.append("</section>")

parts.append('<section id="sec-demux">')
parts.append("<h2>Demultiplexing</h2>")
if demux_stats:
    parts.append("<h3>Demultiplex Stats</h3>")
    for dsp in demux_stats:
        if len(demux_stats) > 1:
            parts.append(f"<h4>{html.escape(dsp)}</h4>")
        parts.append(read_table_preview(dsp, max_rows=None))
else:
    parts.append("<p><em>No demultiplex stats file found.</em></p>")
parts.append("</section>")

parts.append('<section id="sec-fastqc">')
parts.append("<h2>FastQC (Demultiplexed)</h2>")
parts.append(links_block("fastqc_demux/<sample>/*_fastqc.html", fastqc_html))
parts.append("</section>")

parts.append('<section id="sec-atac">')
parts.append("<h2>ATAC QC</h2>")
parts.append(
    "<p><strong>What these files mean:</strong> "
    "<code>.q30.mapped</code> metrics are captured after MAPQ>=30 filtering and before duplicate removal. "
    "<code>.q30.rmdup</code> metrics are captured after duplicate removal. "
    "Use the difference to estimate duplicate burden and retained unique signal.</p>"
)
parts.append('<div class="starsolo-block"><h3>ATAC Key Summary by Sample</h3>')
parts.append(atac_key_summary_table(atac_flagstat_prededup, atac_idxstats_prededup, atac_flagstat_rmdup, atac_idxstats_rmdup, atac_cell_summary))
parts.append("</div>")
parts.append("</section>")

parts.append('<section id="sec-starsolo">')
parts.append("<h2>RNA QC</h2>")
parts.append('<div class="starsolo-block"><h3>Summary.csv Key Metrics by Sample</h3>')
parts.append(summary_csv_key_metrics_table(summary_csv))
parts.append("</div>")
parts.append('<div class="starsolo-block"><h3>Alignment Summary Bar Chart</h3>')
parts.append(alignment_summary_chart(starsolo_logs))
parts.append("</div>")
parts.append('<div class="starsolo-block">')
parts.append(image_gallery_block("STARsolo*/<sample>/*_knee_plot.png", knee_plots))
parts.append("</div>")
parts.append("</section>")

if species_model == 'hybrid':
    parts.append('<section id="sec-barnyard">')
    parts.append("<h2>Hybrid Barnyard Plots</h2>")
    parts.append(image_files_block("STARsolo*/<sample>/*collision_plot.png", barnyard))
    parts.append("</section>")

parts.append("</div></div>")
parts.append("</body></html>")

with open(out_path, "w") as out:
    out.write("\\n".join(parts))

def sample_from_starsolo_path(rel_path):
    bits = rel_path.split("/")
    if len(bits) >= 3 and bits[0] in ("STARsolo", "STARsolo_paired"):
        return bits[1]
    return None

def sample_from_atac_path(rel_path):
    bits = rel_path.split("/")
    if len(bits) >= 3 and bits[0] == "ATAC":
        return bits[1]
    return None

all_sample_candidates = set()
for p in starsolo_logs + knee_plots + barcodes_stats + summary_csv + barnyard + atac_flagstat_rmdup + atac_idxstats_rmdup + atac_stats_rmdup + atac_flagstat_prededup + atac_idxstats_prededup + atac_stats_prededup:
    s = sample_from_starsolo_path(p)
    if s:
        all_sample_candidates.add(s)
    a = sample_from_atac_path(p)
    if a:
        all_sample_candidates.add(a)
all_sample_candidates.update(load_demux_sample_names(demux_stats))

for sample in sorted(all_sample_candidates):
    sample_root = os.path.join(per_sample_dir, sample)
    sample_assets = os.path.join(sample_root, "assets")
    os.makedirs(sample_assets, exist_ok=True)

    sample_logs = [p for p in starsolo_logs if sample_from_starsolo_path(p) == sample]
    sample_knee = [p for p in knee_plots if sample_from_starsolo_path(p) == sample]
    sample_barcodes = [p for p in barcodes_stats if sample_from_starsolo_path(p) == sample]
    sample_summary = [p for p in summary_csv if sample_from_starsolo_path(p) == sample]
    sample_barnyard = [p for p in barnyard if sample_from_starsolo_path(p) == sample]
    sample_fastqc_demux = [p for p in fastqc_html if path_matches_sample(p, sample)]
    sample_fastqc_trimmed = [p for p in fastqc_trimmed_html if path_matches_sample(p, sample)]
    sample_atac_flagstat_rmdup = [p for p in atac_flagstat_rmdup if sample_from_atac_path(p) == sample]
    sample_atac_idxstats_rmdup = [p for p in atac_idxstats_rmdup if sample_from_atac_path(p) == sample]
    sample_atac_stats_rmdup = [p for p in atac_stats_rmdup if sample_from_atac_path(p) == sample]
    sample_atac_flagstat_prededup = [p for p in atac_flagstat_prededup if sample_from_atac_path(p) == sample]
    sample_atac_idxstats_prededup = [p for p in atac_idxstats_prededup if sample_from_atac_path(p) == sample]
    sample_atac_stats_prededup = [p for p in atac_stats_prededup if sample_from_atac_path(p) == sample]
    sample_atac_multiqc_tables = atac_multiqc_tables
    has_atac = bool(sample_atac_flagstat_rmdup or sample_atac_idxstats_rmdup or sample_atac_stats_rmdup or sample_atac_flagstat_prededup or sample_atac_idxstats_prededup or sample_atac_stats_prededup)
    has_rna = bool(sample_logs or sample_knee or sample_barcodes or sample_summary or sample_barnyard)
    sample_tabs = [("sec-demux", "Demultiplexing")]
    if has_atac:
        sample_tabs.append(("sec-atac", "ATAC QC"))
    if has_rna:
        sample_tabs.append(("sec-rna", "RNA QC"))

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

    def sample_links_block(title, paths):
        if not paths:
            return f"<h3>{html.escape(title)}</h3><p><em>No files found.</em></p>"
        items = []
        for p in paths:
            asset = stage_asset_in(p, sample_assets)
            if asset is None:
                continue
            rel_asset = os.path.relpath(asset, sample_root)
            items.append(f'<li><a href="{html.escape(rel_asset)}">{html.escape(p)}</a></li>')
        if not items:
            return f"<h3>{html.escape(title)}</h3><p><em>No readable files found.</em></p>"
        return f"<h3>{html.escape(title)}</h3><ul>{''.join(items)}</ul>"

    def prefer_path(paths, must_contain):
        if not paths:
            return None
        preferred = [p for p in paths if must_contain in p]
        pool = preferred if preferred else paths
        return sorted(pool)[0] if pool else None

    def filter_paths(paths, must_contain):
        if not paths:
            return []
        return sorted([p for p in paths if must_contain in p])

    def atac_pre_post_table(sample):
        pre_flag_path = prefer_path(sample_atac_flagstat_prededup, ".q30.mapped.flagstat.txt")
        pre_idx_path = prefer_path(sample_atac_idxstats_prededup, ".q30.mapped.idxstats.txt")
        post_flag_path = prefer_path(sample_atac_flagstat_rmdup, ".q30.rmdup.flagstat.txt")
        post_idx_path = prefer_path(sample_atac_idxstats_rmdup, ".q30.rmdup.idxstats.txt")

        pre_flag = parse_flagstat_metrics(pre_flag_path) if pre_flag_path else {}
        pre_idx = parse_idxstats_metrics(pre_idx_path) if pre_idx_path else {}
        post_flag = parse_flagstat_metrics(post_flag_path) if post_flag_path else {}
        post_idx = parse_idxstats_metrics(post_idx_path) if post_idx_path else {}
        pre_insert_avg, post_insert_avg = multiqc_insert_size_for_sample(sample_atac_multiqc_tables, sample)

        pre_total = _to_int(pre_flag.get("in_total", ""))
        post_total = _to_int(post_flag.get("in_total", ""))
        retained_pct = (100.0 * post_total / pre_total) if (pre_total and post_total is not None and pre_total > 0) else None
        removed_pct = (100.0 - retained_pct) if retained_pct is not None else None
        pre_mt = (100.0 * pre_idx.get("mt_mapped", 0) / pre_idx.get("total_mapped", 1)) if pre_idx.get("total_mapped", 0) > 0 else 0.0
        post_mt = (100.0 * post_idx.get("mt_mapped", 0) / post_idx.get("total_mapped", 1)) if post_idx.get("total_mapped", 0) > 0 else 0.0

        rows = [
            ("Total reads", pre_flag.get("in_total", ""), post_flag.get("in_total", "")),
            ("Properly paired %", pre_flag.get("properly_paired_pct", ""), post_flag.get("properly_paired_pct", "")),
            ("Singleton %", pre_flag.get("singletons_pct", ""), post_flag.get("singletons_pct", "")),
            ("Mate on different chr", pre_flag.get("mate_diff_chr_count", ""), post_flag.get("mate_diff_chr_count", "")),
            ("Insert size avg", pre_insert_avg, post_insert_avg),
            ("Mitochondrial fraction", _fmt_pct(pre_mt), _fmt_pct(post_mt)),
            ("Autosomal mapped", f"{pre_idx.get('autosomal_mapped', 0):,}", f"{post_idx.get('autosomal_mapped', 0):,}"),
            ("X mapped", f"{pre_idx.get('x_mapped', 0):,}", f"{post_idx.get('x_mapped', 0):,}"),
            ("Y mapped", f"{pre_idx.get('y_mapped', 0):,}", f"{post_idx.get('y_mapped', 0):,}"),
            ("Reads retained", "", _fmt_pct(retained_pct) if retained_pct is not None else ""),
            ("Reads removed", "", _fmt_pct(removed_pct) if removed_pct is not None else ""),
        ]
        cells = ['<tr><th>Metric</th><th>Pre-dedup (.q30.mapped)</th><th>Post-dedup (.q30.rmdup)</th></tr>']
        for m, pre_v, post_v in rows:
            cells.append(f"<tr><td>{html.escape(str(m))}</td><td>{html.escape(str(pre_v))}</td><td>{html.escape(str(post_v))}</td></tr>")
        return "<table>" + "".join(cells) + "</table>"

    sample_parts = []
    sample_parts.append(f'''<!doctype html>
<html>
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>SHARE-seq QC Report - {html.escape(sample)}</title>
  <style>
    body {{
      font-family: "Inter", "Segoe UI", Roboto, "Helvetica Neue", Arial, sans-serif;
      font-size: 14px;
      color: #000000;
      margin: 24px auto;
      max-width: 1200px;
      line-height: 1.55;
      padding: 0 12px;
      background: #ffffff;
      scroll-padding-top: 78px;
    }}
    h1, h2, h3 {{
      color: #000000;
      margin-top: 1.15em;
      margin-bottom: 0.45em;
      font-weight: 700;
      letter-spacing: 0.1px;
    }}
    h1 {{ font-size: 28px; margin-top: 0.3em; }}
    h2 {{ font-size: 20px; border-bottom: 1px solid #dfe1df; padding-bottom: 4px; }}
    h3 {{ font-size: 16px; }}
    code {{ background: #ffffff; color: #000000; border: 1px solid #dfe1df; padding: 2px 5px; border-radius: 4px; font-size: 0.92em; }}
    pre {{
      background: #ffffff;
      border: 1px solid #dfe1df;
      padding: 10px;
      overflow-x: auto;
      line-height: 1.45;
      font-size: 12px;
    }}
    table {{
      border-collapse: collapse;
      margin: 10px 0 20px 0;
      width: 100%;
      background: #fff;
      box-shadow: 0 1px 2px rgba(141, 0, 52, 0.08);
    }}
    th, td {{
      border: 1px solid #dfe1df;
      padding: 6px 10px;
      font-size: 12px;
      vertical-align: top;
    }}
    th {{ background: #ffffff; font-weight: 700; text-align: left; color: #000000; }}
    tr:nth-child(even) td {{ background: #ffffff; }}
    .top-banner {{
      background: #d11947;
      color: #ffffff;
      border-radius: 10px;
      padding: 12px 14px 12px 14px;
      margin: 0 0 14px 0;
      box-shadow: 0 2px 6px rgba(141, 0, 52, 0.22);
    }}
    .top-banner h1, .top-banner p, .top-banner code {{ color: #ffffff; }}
    .top-banner code {{
      background: rgba(255, 255, 255, 0.14);
      border: 1px solid rgba(255, 255, 255, 0.45);
    }}
    section {{
      padding: 10px 14px 8px 14px;
      margin: 8px 0 14px 0;
      background: #ffffff;
      border: 1px solid #dfe1df;
      border-radius: 10px;
      scroll-margin-top: 90px;
    }}
    .tabs-wrap {{
      position: sticky;
      top: 0;
      z-index: 1000;
      background: #d11947;
      border-radius: 0 0 8px 8px;
      padding: 8px 0 6px 0;
      margin-bottom: 14px;
      box-shadow: 0 2px 4px rgba(141, 0, 52, 0.18);
    }}
    .tabs {{
      display: flex;
      gap: 8px;
      align-items: center;
      flex-wrap: wrap;
      margin: 0 6px;
    }}
    .tabs a {{
      text-decoration: none;
      font-size: 13px;
      font-weight: 700;
      color: #ffffff;
      border: 1px solid rgba(255, 255, 255, 0.6);
      border-radius: 8px;
      padding: 6px 10px;
      background: rgba(255, 255, 255, 0.14);
    }}
    .tabs a:hover {{
      background: rgba(255, 255, 255, 0.26);
      border-color: #ffffff;
    }}
  </style>
</head>
<body>
<div class="top-banner">
<h1>SHARE-seq QC Report: {html.escape(sample)}</h1>
<p>Project path: <code>{html.escape(proj)}</code></p>
</div>
''')
    sample_parts.append('<div class="tabs-wrap"><div class="tabs">')
    for sec_id, sec_label in sample_tabs:
        sample_parts.append(f'<a href="#{html.escape(sec_id)}">{html.escape(sec_label)}</a>')
    sample_parts.append("</div></div>")
    sample_parts.append('<section id="sec-demux">')
    sample_parts.append("<h2>Demultiplexing</h2>")
    sample_parts.append(demux_stats_html_for_sample(demux_stats, sample))
    sample_parts.append(sample_links_block("FastQC HTML (demultiplexed)", sample_fastqc_demux))
    sample_parts.append(sample_links_block("FastQC HTML (trimmed)", sample_fastqc_trimmed))
    sample_parts.append("</section>")
    if has_atac:
        sample_parts.append('<section id="sec-atac">')
        sample_parts.append("<h2>ATAC QC</h2>")
        sample_parts.append(atac_pre_post_table(sample))
        sample_parts.append(sample_links_block("ATAC MultiQC report", atac_multiqc_report))
        sample_parts.append("</section>")
    if has_rna:
        sample_parts.append('<section id="sec-rna">')
        sample_parts.append("<h2>RNA QC</h2>")
        sample_parts.append(text_files_block("Log.final.out", sample_logs, max_lines=120))
        sample_parts.append(sample_image_block("Knee plots", sample_knee))
        sample_parts.append(text_files_block("Barcodes.stats", sample_barcodes, max_lines=120))
        sample_parts.append(table_files_block("GeneFull Summary.csv", sample_summary, max_rows=20))
        sample_parts.append(sample_image_block("Barnyard plots (if hybrid)", sample_barnyard))
        sample_parts.append("</section>")
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

    def bwaPrefixName = (params.species_model == 'mouse') ? 'mm39_bwa' : (params.species_model == 'hybrid' ? 'hybrid_bwa' : 'hg38_bwa')

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
    def ch_report_barrier = BUILD_DEMUX_STATS_FROM_MERGED.out.merged_stats
        .mix(RENAME_FASTQ.out.rename_stats)
        .mix(FASTQC_DEMUX.out.demux_reports)
        .mix(POLYT_FILTER.out.polyt_outputs)
        .mix(ch_trim_completion)
        .mix(DETERMINE_READ_LENGTH.out.read_length)
        .mix(PREPEND_HEADER_BARCODES.out.r2_with_barcodes)
        .mix(STAR_INDEX.out.star_index)
        .mix(BWA_ALIGN_ATAC.out.atac_align_out)
        .mix(ESTIMATE_ATAC_CELLS.out.atac_cell_summary)
        .mix(ESTIMATE_ATAC_CELLS.out.atac_cell_counts)
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
        if( params.species_model == 'hybrid' ) {
            BARNYARD_PLOT(ch_starsolo_for_hybrid_qc)
            HYBRID_SPLIT_SPECIES(ch_starsolo_for_hybrid_qc)
            ch_report_barrier = ch_report_barrier
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
        if( params.species_model == 'hybrid' ) {
            BARNYARD_PLOT(ch_starsolo_paired_for_qc)
            HYBRID_SPLIT_SPECIES(ch_starsolo_paired_for_qc)
            ch_report_barrier = ch_report_barrier
                .mix(BARNYARD_PLOT.out.barnyard80)
                .mix(BARNYARD_PLOT.out.barnyard90)
                .mix(HYBRID_SPLIT_SPECIES.out.split_09)
                .mix(HYBRID_SPLIT_SPECIES.out.split_085)
                .mix(HYBRID_SPLIT_SPECIES.out.split_08)
        }
    } else {
        log.info "Unknown star_alignment_mode = ${params.star_alignment_mode}; skipping STARsolo alignment."
    }

    ch_report_barrier
        .collect()
        .map { trigger_items ->
            def knee_items = trigger_items.findAll { x ->
                x != null && x.toString().endsWith("_knee_plot.png")
            }
            tuple(1, knee_items)
        }
        .set { ch_build_qc_report }
    BUILD_QC_HTML(ch_build_qc_report)
}

