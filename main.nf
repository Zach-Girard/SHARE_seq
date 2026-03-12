nextflow.enable.dsl=2

/*
 * Minimal starter pipeline.
 * Assumes:
 *  - RAW_FASTQ/ contains raw FASTQ files
 *  - Genomes/ and GTF/ have been prepared using the helper scripts
 */

params.raw_fastq           = params.raw_fastq           ?: 'RAW_FASTQ'
params.genomes_dir         = params.genomes_dir         ?: 'Genomes'
params.gtf_dir             = params.gtf_dir             ?: 'GTF'
params.skip_demux          = params.skip_demux          ?: false
params.skip_trim           = params.skip_trim           ?: false
params.skip_polyt          = params.skip_polyt          ?: false
params.umi_len             = params.umi_len             ?: 10
params.bc_coords           = params.bc_coords           ?: '15-23,53-61,91-99'
params.species_model       = params.species_model       ?: 'human'  // 'human', 'mouse', or 'hybrid'
params.star_alignment_mode = params.star_alignment_mode ?: 'single' // 'single' or 'paired'
params.barcodes_8bp_file   = params.barcodes_8bp_file   ?: 'barcodes_RC.txt'
params.barcodes_rc         = (params.barcodes_rc in [true, 'true'])

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

// Log key configuration
log.info "Using RAW_FASTQ directory    : ${params.raw_fastq}"
log.info "Using Genomes directory      : ${params.genomes_dir}"
log.info "Using GTF directory          : ${params.gtf_dir}"
log.info "Species model                : ${params.species_model}"
log.info "STARsolo alignment mode      : ${params.star_alignment_mode}"
log.info "Skip demultiplexing          : ${params.skip_demux}"
log.info "Skip adapter trimming        : ${params.skip_trim}"
log.info "Skip Poly-T filtering        : ${params.skip_polyt}"
log.info "UMI length (R3)              : ${params.umi_len}"
log.info "Barcode coords (R2)          : ${params.bc_coords}"
log.info "User 8bp barcode file        : ${params.barcodes_8bp_file}"
log.info "Reverse-complement barcodes? : ${params.barcodes_rc}"
log.info "Effective 8bp whitelist file : ${effectiveCbWhitelistPath}"

/*
 * Channel definitions
 */

Channel
    .fromPath("${params.raw_fastq}/*.fastq.gz", checkIfExists: true)
    .ifEmpty { error "No FASTQ files found in: ${params.raw_fastq}" }
    .set { ch_raw_fastq }

/*
 * Processes
 */

// Placeholder demultiplexing step: currently just passes FASTQs through.
// Later, replace this with a real SHARE-seq demux script.
process DEMULTIPLEX_PLACEHOLDER {
    tag { file(raw_fastq).name }

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

    input:
    path fastq

    output:
    path "fastqc_raw/${fastq.simpleName}_fastqc.*", emit: raw_reports

    """
    mkdir -p fastqc_raw
    fastqc -o fastqc_raw -f fastq ${fastq}
    """
}

process TRIM_FASTQ {
    tag { file(fastq).name }

    input:
    path fastq

    output:
    path "trimmed/${fastq.simpleName}.trimmed.fastq.gz", emit: trimmed_fastq
    path "trimmed/${fastq.simpleName}.fastp.json"
    path "trimmed/${fastq.simpleName}.fastp.html"

    """
    mkdir -p trimmed
    fastp \\
      -i ${fastq} \\
      -o trimmed/${fastq.simpleName}.trimmed.fastq.gz \\
      -w ${task.cpus} \\
      -j trimmed/${fastq.simpleName}.fastp.json \\
      -h trimmed/${fastq.simpleName}.fastp.html
    """
}

process FASTQC_TRIMMED {
    tag { file(fastq).name }

    input:
    path fastq

    output:
    path "fastqc_trimmed/${fastq.simpleName}_fastqc.*", emit: trimmed_reports

    """
    mkdir -p fastqc_trimmed
    fastqc -o fastqc_trimmed -f fastq ${fastq}
    """
}

process POLYT_FILTER {
    tag { file(r3_fastq).name }

    input:
    path r3_fastq

    output:
    path "polyt_filtered/*", emit: polyt_outputs

    """
    mkdir -p polyt_filtered
    bash PolyT_cutadapt.sh ${r3_fastq}
    """
}

process ADD_R2_BARCODES_TO_R3 {
    tag { file(r3_fastq).name }

    input:
    tuple path(r2_fastq), path(r3_fastq)

    output:
    path "withBarcodes_*", emit: r3_with_barcodes

    """
    python Read3_Barcode_Addition.py \\
      -r3 ${r3_fastq} \\
      -r2 ${r2_fastq} \\
      -bc ${params.bc_coords} \\
      -umi 0-${params.umi_len}
    """
}

process DETERMINE_READ_LENGTH {
    tag { file(r1_fastq).name }

    input:
    path r1_fastq

    output:
    path "read_length.txt", emit: read_length

    """
    # Determine read length from the second line of the FASTQ file
    len=$(zcat ${r1_fastq} 2>/dev/null | sed -n '2p' | awk '{print length($0)}' || \\
         sed -n '2p' ${r1_fastq} | awk '{print length($0)}')
    echo ${len} > read_length.txt
    """
}

process STAR_INDEX {
    tag { params.species_model }

    input:
    path read_length_file

    output:
    path "STAR_index_*", emit: star_index

    script:
    def readLen   = new File(read_length_file.toString()).text.trim().toInteger()
    def overhang  = readLen - 1

    def (genomeFasta, gtfFile, indexDir) = {
        switch (params.species_model) {
            case 'mouse':
                [
                  "${params.genomes_dir}/GRCm39/Mus_musculus.GRCm39.dna.primary_assembly.fa.gz",
                  "${params.gtf_dir}/Mus_musculus/Mus_musculus.GRCm39.111.gtf.gz",
                  "STAR_index_mouse${readLen}bp"
                ]
            case 'hybrid':
                [
                  "${params.genomes_dir}/hybrid/hybrid_human_mouse.fa.gz",
                  "${params.gtf_dir}/hybrid/hybrid_human_mouse.gtf.gz",
                  "STAR_index_hybrid${readLen}bp"
                ]
            case 'human':
            default:
                [
                  "${params.genomes_dir}/GRCh38/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz",
                  "${params.gtf_dir}/Homo_sapiens/Homo_sapiens.GRCh38.111.gtf.gz",
                  "STAR_index_human${readLen}bp"
                ]
        }
    }()

    """
    mkdir -p ${indexDir}
    STAR --runMode genomeGenerate \\
         --runThreadN ${task.cpus} \\
         --genomeDir ${indexDir} \\
         --genomeFastaFiles ${genomeFasta} \\
         --sjdbGTFfile ${gtfFile} \\
         --sjdbOverhang ${overhang}
    """
}

process STARSOLO_SINGLE {
    tag { sample_id }

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
      --soloUMIposition 0_24_0_33 \\
      --soloBarcodeReadLength 0 \\
      --soloFeatures Gene GeneFull \\
      --soloStrand Unstranded \\
      --outSAMtype BAM SortedByCoordinate \\
      --outFileNamePrefix STARsolo/${sample_id}/
    """
}

process KNEE_PLOT {
    tag { sample_id }

    input:
    tuple val(sample_id), path(starsolo_dir)

    output:
    path "STARsolo/${sample_id}/${sample_id}_knee_plot.png", emit: knee_plot

    """
    GENEFULL_DIR="${starsolo_dir}/Solo.out/GeneFull"
    SUMMARY_CSV="${GENEFULL_DIR}/Summary.csv"

    if [ ! -f "${SUMMARY_CSV}" ]; then
      echo "Summary.csv not found for sample ${sample_id} at ${SUMMARY_CSV}" >&2
      exit 1
    fi

    EST=$(awk -F',' 'BEGIN{OFS=","} $1 ~ /^Estimated Number of Cells/ {gsub(/^[ \t]+|[ \t]+$/,"",$2); print $2; exit}' "${SUMMARY_CSV}")

    if [ -z "${EST}" ]; then
      echo "Could not extract Estimated Number of Cells from ${SUMMARY_CSV}" >&2
      exit 1
    fi

    python "${projectDir}/Visualization_scripts/Knee_plot.py" \\
      --genefull-dir "${GENEFULL_DIR}" \\
      --estimated-cells "${EST}" \\
      --output "STARsolo/${sample_id}/${sample_id}_knee_plot.png"
    """
}

process BARNYARD_PLOT {
    tag { sample_id }

    input:
    tuple val(sample_id), path(starsolo_dir)

    output:
    path "STARsolo/${sample_id}/80%_collision_plot.png", emit: barnyard80
    path "STARsolo/${sample_id}/90%_collision_plot.png", emit: barnyard90

    """
    FILTERED_DIR="${starsolo_dir}/Solo.out/GeneFull/filtered"

    if [ ! -d "${FILTERED_DIR}" ]; then
      echo "Filtered GeneFull directory not found for sample ${sample_id} at ${FILTERED_DIR}" >&2
      exit 1
    fi

    # 80% purity barnyard plot
    python "${projectDir}/Visualization_scripts/BarnyardPlot.py" \\
      --filtered-dir "${FILTERED_DIR}" \\
      --human-threshold 0.8 \\
      --mouse-threshold 0.2 \\
      --collision-low 0.2 \\
      --collision-high 0.8 \\
      --title "Human-Mouse Collision (80% Purity)" \\
      --output "STARsolo/${sample_id}/80%_collision_plot.png"

    # 90% purity barnyard plot
    python "${projectDir}/Visualization_scripts/BarnyardPlot.py" \\
      --filtered-dir "${FILTERED_DIR}" \\
      --human-threshold 0.9 \\
      --mouse-threshold 0.1 \\
      --collision-low 0.1 \\
      --collision-high 0.9 \\
      --title "Human-Mouse Collision (90% Purity)" \\
      --output "STARsolo/${sample_id}/90%_collision_plot.png"
    """
}

process HYBRID_SPLIT_SPECIES {
    tag { sample_id }

    input:
    tuple val(sample_id), path(starsolo_dir)

    output:
    path "STARsolo/${sample_id}/species_split_purity_0.9", emit: split_09
    path "STARsolo/${sample_id}/species_split_purity_0.85", emit: split_085
    path "STARsolo/${sample_id}/species_split_purity_0.8", emit: split_08

    """
    FILTERED_DIR="${starsolo_dir}/Solo.out/GeneFull/filtered"

    if [ ! -d "${FILTERED_DIR}" ]; then
      echo "Filtered GeneFull directory not found for sample ${sample_id} at ${FILTERED_DIR}" >&2
      exit 1
    fi

    # Purity 0.9
    python "${projectDir}/Split_Species_By_Purity.py" \\
      --input "${FILTERED_DIR}" \\
      --output "STARsolo/${sample_id}/species_split_purity_0.9" \\
      --purity 0.9

    # Purity 0.85
    python "${projectDir}/Split_Species_By_Purity.py" \\
      --input "${FILTERED_DIR}" \\
      --output "STARsolo/${sample_id}/species_split_purity_0.85" \\
      --purity 0.85

    # Purity 0.8
    python "${projectDir}/Split_Species_By_Purity.py" \\
      --input "${FILTERED_DIR}" \\
      --output "STARsolo/${sample_id}/species_split_purity_0.8" \\
      --purity 0.8
    """
}

process PAIRED_BARCODE_MATCH {
    tag { sample_id }

    input:
    tuple val(sample_id), path(r1_fastq), path(r3_with_barcodes)

    output:
    path "paired/${sample_id}_R1.paired.fastq.gz", emit: r1_paired
    path "paired/${sample_id}_R3.paired.fastq.gz", emit: r3_paired
    path "paired/${sample_id}_barcode_match_summary.txt", emit: summary

    """
    mkdir -p paired
    python "${projectDir}/Paired_Barcode_Matcher.py" \\
      --r1 ${r1_fastq} \\
      --r3 ${r3_with_barcodes} \\
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
      --barcodes ${effectiveCbWhitelistPath} \\
      --output whitelist_paired.txt
    """
}

process STARSOLO_PAIRED {
    tag { sample_id }

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
      --soloUMIlen 10 \\
      --soloBarcodeReadLength 0 \\
      --soloBarcodeMate 2 \\
      --clip5pNbases 0 49 \\
      --outFilterMatchNmin 20 \\
      --outFilterMatchNminOverLread 0 \\
      --outFilterScoreMinOverLread 0 \\
      --alignEndsType Local \\
      --soloFeatures Gene GeneFull \\
      --soloStrand Unstranded \\
      --outSAMtype BAM SortedByCoordinate \\
      --outFileNamePrefix STARsolo_paired/${sample_id}/${sample_id}_paired_
    """
}

/*
 * Main workflow
 */

workflow {
    main:
    def ch_input_fastq = ch_raw_fastq

    if( !params.skip_demux ) {
        DEMULTIPLEX_PLACEHOLDER(ch_input_fastq)
        ch_input_fastq = DEMULTIPLEX_PLACEHOLDER.out.demuxed_fastq
    } else {
        log.info "Skipping demultiplexing step (using input FASTQs as already demultiplexed)."
    }

    // FastQC on raw (or demuxed) FASTQs
    FASTQC_RAW(ch_input_fastq)

    // Optional Poly-T filtering on R3 (UMI+PolyT+cDNA) FASTQs
    if( !params.skip_polyt ) {
        // Assumes that after demultiplexing you have three FASTQs per sample:
        //  R1 = cDNA, R2 = barcodes+linkers, R3 = UMI + PolyT + cDNA
        // and that R3 files are named with '_R3' in the filename.
        Channel
            .fromPath("demux/*_R3*.fastq.gz", checkIfExists: false)
            .tap { ch ->
                if( ch.empty ) log.warn "No R3 FASTQ files found in demux/ for Poly-T filtering."
            }
            .set { ch_r3_for_polyt }

        POLYT_FILTER(ch_r3_for_polyt)
    } else {
        log.info "Skipping Poly-T filtering step on R3."
    }

    // Decide which R2/R3 FASTQs to use for barcode prepending:
    //  - If trimming ran: always use trimmed R2/R3 (these will be from Poly-T–filtered or demux, depending on skip_polyt).
    //  - Else if Poly-T filtering ran: use Poly-T–filtered R2/R3.
    //  - Else: fall back to demultiplexed R2/R3.
    def ch_r2_source
    def ch_r3_source

    if( !params.skip_trim ) {
        log.info "Using trimmed R2/R3 for barcode prepending."
        ch_r2_source = Channel.fromPath("trimmed/*_R2*.trimmed.fastq.gz", checkIfExists: false)
        ch_r3_source = Channel.fromPath("trimmed/*_R3*.trimmed.fastq.gz", checkIfExists: false)
    } else if( !params.skip_polyt ) {
        log.info "Using Poly-T–filtered R2/R3 for barcode prepending."
        ch_r2_source = Channel.fromPath("polyt_filtered/*_R2*.fastq.gz", checkIfExists: false)
        ch_r3_source = Channel.fromPath("polyt_filtered/*_R3*.fastq.gz", checkIfExists: false)
    } else {
        log.info "Using demultiplexed R2/R3 for barcode prepending."
        ch_r2_source = Channel.fromPath("demux/*_R2*.fastq.gz", checkIfExists: false)
        ch_r3_source = Channel.fromPath("demux/*_R3*.fastq.gz", checkIfExists: false)
    }

    ch_r2_source
        .map { r2 -> tuple(r2.baseName.replaceFirst('_R2',''), r2) }
        .join(
            ch_r3_source.map { r3 -> tuple(r3.baseName.replaceFirst('_R3',''), r3) }
        )
        .map { sample_id, r2, r3 -> tuple(r2, r3) }
        .set { ch_r2_r3_pairs }

    ADD_R2_BARCODES_TO_R3(ch_r2_r3_pairs)

    // Decide which FASTQs to trim:
    //  - If Poly-T filtering ran: trim the Poly-T–filtered FASTQs (R1/R2/R3).
    //  - Else: trim the demultiplexed FASTQs.
    def ch_for_trim
    if( !params.skip_polyt ) {
        log.info "Trimming Poly-T–filtered FASTQs (R1/R2/R3)."
        ch_for_trim = Channel.fromPath("polyt_filtered/*.fastq.gz", checkIfExists: false)
    } else {
        log.info "Trimming demultiplexed FASTQs (R1/R2/R3)."
        ch_for_trim = Channel.fromPath("demux/*.fastq.gz", checkIfExists: false)
    }

    // Optional trimming, then FastQC on trimmed FASTQs
    if( !params.skip_trim ) {
        TRIM_FASTQ(ch_for_trim)
        FASTQC_TRIMMED(TRIM_FASTQ.out.trimmed_fastq)
    } else {
        log.info "Skipping adapter/quality trimming."
    }

    // Determine read length from R1 to drive STAR index sjdbOverhang
    def ch_r1_for_index
    if( !params.skip_trim ) {
        log.info "Using trimmed R1 to determine read length for STAR index."
        ch_r1_for_index = Channel.fromPath("trimmed/*_R1*.trimmed.fastq.gz", checkIfExists: false)
    } else if( !params.skip_polyt ) {
        log.info "Using Poly-T–filtered R1 to determine read length for STAR index."
        ch_r1_for_index = Channel.fromPath("polyt_filtered/*_R1*.fastq.gz", checkIfExists: false)
    } else {
        log.info "Using demultiplexed R1 to determine read length for STAR index."
        ch_r1_for_index = Channel.fromPath("demux/*_R1*.fastq.gz", checkIfExists: false)
    }

    ch_r1_for_index
        .take(1)
        .set { ch_r1_single }

    DETERMINE_READ_LENGTH(ch_r1_single)
    STAR_INDEX(DETERMINE_READ_LENGTH.out.read_length)

    // Single-end STARsolo alignment using (trimmed, Poly-T-matched) R1 and withBarcodes_R3 when available
    if( params.star_alignment_mode == 'single' ) {
        log.info "Running single-end STARsolo alignment."

        // Choose R1 source for alignment:
        //  - If trimming ran: use trimmed R1.
        //  - Else if Poly-T filtering ran: use Poly-T–filtered R1.
        //  - Else: use demultiplexed R1.
        def ch_r1_source_for_align
        if( !params.skip_trim ) {
            log.info "Using trimmed R1 for alignment."
            ch_r1_source_for_align = Channel.fromPath("trimmed/*_R1*.trimmed.fastq.gz", checkIfExists: false)
        } else if( !params.skip_polyt ) {
            log.info "Using Poly-T–filtered R1 for alignment."
            ch_r1_source_for_align = Channel.fromPath("polyt_filtered/*_R1*.fastq.gz", checkIfExists: false)
        } else {
            log.info "Using demultiplexed R1 for alignment."
            ch_r1_source_for_align = Channel.fromPath("demux/*_R1*.fastq.gz", checkIfExists: false)
        }

        ch_r1_source_for_align
            .map { r1 -> tuple(r1.baseName.replaceFirst('_R1',''), r1) }
            .set { ch_r1_for_align }

        ADD_R2_BARCODES_TO_R3.out.r3_with_barcodes
            .map { r3 -> 
                def bn = r3.baseName.replaceFirst('withBarcodes_','')
                def sample = bn.replaceFirst('_R3','')
                tuple(sample, r3)
            }
            .set { ch_r3_barcode }

        ch_r1_for_align
            .join(ch_r3_barcode)
            .map { sample_id, r1, r3 -> tuple(sample_id, r1, r3) }
            .combine(STAR_INDEX.out.star_index) { triple, idx ->
                def (sample_id, r1, r3) = triple
                tuple(sample_id, r1, r3, idx)
            }
            .set { ch_starsolo_single }

        STARSOLO_SINGLE(ch_starsolo_single)

        // Generate per-sample GeneFull knee plots from STARsolo outputs
        STARSOLO_SINGLE.out.starsolo_out
            .map { dir ->
                def sample = dir.baseName  // STARsolo/<sample>/
                tuple(sample, dir)
            }
            .set { ch_starsolo_for_hybrid_qc }

        KNEE_PLOT(ch_starsolo_for_hybrid_qc)

        // For hybrid mixed-species model, also generate barnyard plots and split species by purity
        if( params.species_model == 'hybrid' ) {
            BARNYARD_PLOT(ch_starsolo_for_hybrid_qc)
            HYBRID_SPLIT_SPECIES(ch_starsolo_for_hybrid_qc)
        }
    } else if( params.star_alignment_mode == 'paired' ) {
        log.info "Running paired-end barcode matching; STARsolo paired-end alignment is currently a placeholder."

        // Choose R1 source for paired-end alignment (same logic as single-end)
        def ch_r1_source_for_align_paired
        if( !params.skip_trim ) {
            log.info "Using trimmed R1 for paired alignment."
            ch_r1_source_for_align_paired = Channel.fromPath("trimmed/*_R1*.trimmed.fastq.gz", checkIfExists: false)
        } else if( !params.skip_polyt ) {
            log.info "Using Poly-T–filtered R1 for paired alignment."
            ch_r1_source_for_align_paired = Channel.fromPath("polyt_filtered/*_R1*.fastq.gz", checkIfExists: false)
        } else {
            log.info "Using demultiplexed R1 for paired alignment."
            ch_r1_source_for_align_paired = Channel.fromPath("demux/*_R1*.fastq.gz", checkIfExists: false)
        }

        ch_r1_source_for_align_paired
            .map { r1 -> tuple(r1.baseName.replaceFirst('_R1',''), r1) }
            .set { ch_r1_for_paired }

        // Use withBarcodes_R3 output as the R3 source for matching
        ADD_R2_BARCODES_TO_R3.out.r3_with_barcodes
            .map { r3 ->
                def bn = r3.baseName.replaceFirst('withBarcodes_','')
                def sample = bn.replaceFirst('_R3','')
                tuple(sample, r3)
            }
            .set { ch_r3_barcode_for_paired }

        ch_r1_for_paired
            .join(ch_r3_barcode_for_paired)
            .map { sample_id, r1, r3 -> tuple(sample_id, r1, r3) }
            .set { ch_pairs_for_matching }

        // Run barcode matching to filter read pairs
        PAIRED_BARCODE_MATCH(ch_pairs_for_matching)

        // Combine matched pairs with STAR index for paired STARsolo
        PAIRED_BARCODE_MATCH.out.r1_paired
            .map { r1p ->
                def sample = r1p.baseName.replaceFirst('_R1\\.paired','')
                tuple(sample, r1p)
            }
            .set { ch_r1_paired_by_sample }

        PAIRED_BARCODE_MATCH.out.r3_paired
            .map { r3p ->
                def sample = r3p.baseName.replaceFirst('_R3\\.paired','')
                tuple(sample, r3p)
            }
            .set { ch_r3_paired_by_sample }

        def ch_starsolo_paired = ch_r1_paired_by_sample
            .join(ch_r3_paired_by_sample)
            .map { sample_id, r1p, r3p -> tuple(sample_id, r1p, r3p) }
            .combine(STAR_INDEX.out.star_index) { triple, idx ->
                def (sample_id, r1p, r3p) = triple
                tuple(sample_id, r1p, r3p, idx)
            }
            .set { ch_starsolo_paired }

        // Build paired-end 24bp whitelist for STARsolo using effective 8bp barcodes
        Channel
            .of(1)
            .set { ch_build_whitelist }
        BUILD_PAIRED_WHITELIST(ch_build_whitelist)

        // Broadcast whitelist to all samples and run paired-end STARsolo
        BUILD_PAIRED_WHITELIST.out.paired_whitelist_file
            .map { wl -> tuple(wl) }
            .set { ch_whitelist_broadcast }

        ch_starsolo_paired
            .combine(ch_whitelist_broadcast) { triple, wl ->
                def (sample_id, r1p, r3p, idx) = triple
                tuple(sample_id, r1p, r3p, idx, wl)
            }
            .set { ch_starsolo_paired_with_wl }

        STARSOLO_PAIRED(ch_starsolo_paired_with_wl)
    } else {
        log.info "Unknown star_alignment_mode = ${params.star_alignment_mode}; skipping STARsolo alignment."
    }
}

