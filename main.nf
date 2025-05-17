#!/usr/bin/env nextflow

// Import modules
include { KRAKEN2 } from './modules/local/kraken2/main'
include { BRACKEN } from './modules/local/bracken/main'
include { KRAKEN2KRONA; KRONA_PLOT } from './modules/local/krona/main'
include { GENERATE_SUMMARY; ALPHA_DIVERSITY; BETA_DIVERSITY } from './modules/local/diversity/main'
include { ARGS_OAP } from './modules/local/args_oap/main' // Import ARGs_OAP module
include { FASTP } from './modules/local/fastp/main' // Import FASTP module
include { FASTQC as FASTQC_RAW; FASTQC as FASTQC_TRIMMED; FASTQC as FASTQC_PHIX_REMOVED } from './modules/local/fastqc/main'
include { BOWTIE2_BUILD_PHIX; BOWTIE2_REMOVAL_PHIX } from './modules/local/bowtie2/main'
include { MEGAHIT } from './modules/local/megahit/main' // Import MEGAHIT module

// Define parameters
params.reads_dir = "$baseDir/reads"  // Thư mục chứa tất cả raw reads
params.outdir = null                 // Thư mục đầu ra PHẢI được chỉ định
params.krakendb = "$baseDir/krakendb"
params.kraken_threads = 8
params.bracken_threads = 8
params.bracken_threshold = 10
params.bracken_length = 150
params.bracken_level = ["P","C","O","F","G","S"] // Danh sách các cấp độ phân loại
params.read_pattern = "*_{1,2}*.{fq,fastq}.gz"  // Pattern mở rộng để tìm nhiều kiểu tên file hơn
params.run_args_oap = true  // Flag to enable/disable ARGs_OAP analysis

// Thêm tham số Docker
params.enable_docker = true  // Mặc định bật Docker
params.container_kraken2 = 'quay.io/biocontainers/kraken2:2.1.2--pl5321h9f5acd7_2'
params.container_bracken = 'quay.io/biocontainers/bracken:2.7--py39hc16433a_0'
params.container_krona = 'quay.io/biocontainers/krona:2.8--pl5262hdfd78af_2'
params.container_args_oap = 'quay.io/biocontainers/args_oap:3.2.4--pyhdfd78af_0' // Docker container for ARGs_OAP
params.container_fastp = 'quay.io/biocontainers/fastp:0.24.1--heae3180_0' // Docker container for FASTP
params.container_fastqc = 'quay.io/biocontainers/fastqc:0.12.1--hdfd78af_0' // Docker container for FASTQC
params.container_bowtie2 = 'quay.io/biocontainers/bowtie2:2.5.4--he96a11b_5' // Docker container for Bowtie2
params.container_megahit = 'quay.io/biocontainers/megahit:1.2.9--h2e03b76_1' // Docker container for MEGAHIT

// Thêm tham số cho Bowtie2 và PhiX removal
params.phix_fasta = "$baseDir/assets/GCA_002596845.1_ASM259684v1_genomic.fna.gz" // Đường dẫn đến file FASTA của PhiX
params.bowtie2_build_threads = params.kraken_threads ?: 8 // Số luồng cho bowtie2-build
params.bowtie2_align_threads = params.kraken_threads ?: 8 // Số luồng cho bowtie2 alignment

// Thêm tham số cho MEGAHIT
params.megahit_threads = 8 // Số luồng cho MEGAHIT

// Workflow definition
workflow {
    // Kiểm tra tham số trước khi thực hiện workflow
    // Kiểm tra tham số outdir có được cung cấp hay không
    if (params.outdir == null) {
        log.error "Thư mục đầu ra (outdir) không được chỉ định!"
        log.error "Vui lòng chỉ định đường dẫn đầu ra bằng tham số --outdir"
        exit 1
    }

    // Chuyển đổi bracken_level nếu nó là chuỗi
    def bracken_levels_list = params.bracken_level instanceof String ?
                        [params.bracken_level] : params.bracken_level

    log.info """\
             K R A K E N 2   A N D   B R A C K E N   P I P E L I N E    
             ===================================
             Reads Directory    : ${params.reads_dir}
             Read pattern       : ${params.read_pattern}
             Kraken2 Database   : ${params.krakendb}
             Kraken2 threads    : ${params.kraken_threads}
             Bracken threads    : ${params.bracken_threads}
             Bracken threshold  : ${params.bracken_threshold}
             Bracken read length: ${params.bracken_length}
             Bracken tax levels : ${bracken_levels_list.join(', ')}
             ARGs_OAP enabled   : ${params.run_args_oap}
             Output directory   : ${params.outdir}
             Using Docker       : ${params.enable_docker}
             Kraken2 container  : ${params.container_kraken2}
             Bracken container  : ${params.container_bracken}
             Krona container    : ${params.container_krona}
             ARGs_OAP container : ${params.container_args_oap}
             FASTP container    : ${params.container_fastp}
             FASTQC container   : ${params.container_fastqc}
             PhiX FASTA         : ${params.phix_fasta}
             Bowtie2 build threads: ${params.bowtie2_build_threads}
             Bowtie2 align threads: ${params.bowtie2_align_threads}
             Bowtie2 container  : ${params.container_bowtie2}
             MEGAHIT threads    : ${params.megahit_threads}
             MEGAHIT container  : ${params.container_megahit}
             """
             .stripIndent()

    // Kiểm tra database Kraken2 có tồn tại không
    def krakendb_path = file(params.krakendb)
    if( !krakendb_path.exists() ) {
        log.error "Kraken2 database not found at: ${params.krakendb}"
        log.error "Please provide the correct path using --krakendb parameter"
        exit 1
    }

    // Kiểm tra thư mục chứa reads có tồn tại không
    def reads_dir_path = file(params.reads_dir)
    if( !reads_dir_path.exists() ) {
        log.error "Reads directory not found at: ${params.reads_dir}"
        log.error "Please provide the correct path using --reads_dir parameter"
        exit 1
    }

    // Kiểm tra file PhiX FASTA có tồn tại không
    def phix_fasta_path = file(params.phix_fasta)
    if( !phix_fasta_path.exists() ) {
        log.error "PhiX FASTA file not found at: ${params.phix_fasta}"
        log.error "Please provide the correct path using --phix_fasta parameter, or ensure the default file exists at assets/GCA_002596845.1_ASM259684v1_genomic.fna.gz"
        exit 1
    }

    // Tạo channel từ các file reads trong thư mục
    read_pairs_ch = Channel
        .fromFilePairs("${params.reads_dir}/${params.read_pattern}", checkIfExists: true)
        .map { sample_id, reads -> tuple(["id": sample_id], reads) }  // Sửa meta map với trường id
        .ifEmpty { error "No read files found in ${params.reads_dir} with pattern ${params.read_pattern}" }

    // Run FastQC for initial quality control
    FASTQC_RAW(read_pairs_ch)

    // Run FASTP for quality control and trimming
    FASTP(read_pairs_ch)
    
    // Run FastQC on trimmed reads to assess quality after trimming
    FASTQC_TRIMMED(FASTP.out.trimmed_reads)
    
    // Build PhiX index using Bowtie2
    phix_fasta_ch = file(params.phix_fasta)
    BOWTIE2_BUILD_PHIX(phix_fasta_ch, params.bowtie2_build_threads)

    // Remove PhiX reads using Bowtie2
    BOWTIE2_REMOVAL_PHIX(FASTP.out.trimmed_reads, BOWTIE2_BUILD_PHIX.out.index_dir, params.bowtie2_align_threads)

    // Run FastQC on PhiX-filtered reads
    FASTQC_PHIX_REMOVED(BOWTIE2_REMOVAL_PHIX.out.reads)
    
    // Run MEGAHIT for metagenomic assembly
    MEGAHIT(BOWTIE2_REMOVAL_PHIX.out.reads, params.megahit_threads)

    // Run Kraken2 on the trimmed and PhiX-filtered reads
    BOWTIE2_REMOVAL_PHIX.out.reads
        .map { meta, reads -> tuple(meta.id, reads) }  // Chuyển đổi meta map sang sample_id
        .set { kraken_input_ch }
    KRAKEN2(kraken_input_ch)
    
    // Tạo channel với tất cả các cấp độ phân loại
    def bracken_levels_list_ch = params.bracken_level instanceof String ?
                        [params.bracken_level] : params.bracken_level
    bracken_ch = KRAKEN2.out.report.combine(Channel.fromList(bracken_levels_list_ch))
        .map { sample_id, kraken_report, level -> 
            tuple(sample_id, kraken_report, level)
        }
    
    // Chạy Bracken cho mỗi cấp độ phân loại
    BRACKEN(bracken_ch)
    
    // Chạy KRAKEN2KRONA để chuyển đổi báo cáo Kraken2 sang định dạng Krona
    KRAKEN2.out.report
        .map { sample_id, kraken_report -> tuple(sample_id, kraken_report) }
        .set { kraken_reports_ch }
    KRAKEN2KRONA(kraken_reports_ch)
    
    // Tạo HTML Krona charts từ text files
    KRONA_PLOT(KRAKEN2KRONA.out.krona_text)
    
    // Nhóm các báo cáo Bracken theo cấp độ phân loại (cho summary reports)
    BRACKEN.out.reports
        .map { sample_id, level, bracken_report -> tuple(level, bracken_report) }
        .groupTuple()
        .set { grouped_bracken_reports_ch }
    
    // Tạo báo cáo tổng hợp cho từng cấp độ phân loại từ báo cáo Bracken
    GENERATE_SUMMARY(grouped_bracken_reports_ch)
    
    // Nhóm các báo cáo Bracken cho phân tích đa dạng sinh học
    BRACKEN.out.reports
        .map { sample_id, level, bracken_report -> tuple(level, bracken_report) }
        .groupTuple()
        .set { grouped_bracken_reports_ch_for_diversity }
    
    // Tính toán đa dạng sinh học Alpha từ báo cáo Bracken
    ALPHA_DIVERSITY(grouped_bracken_reports_ch_for_diversity)
    
    // Tính toán đa dạng sinh học Beta từ báo cáo Bracken
    BETA_DIVERSITY(grouped_bracken_reports_ch_for_diversity)
    
    // Run ARGs_OAP if enabled
    if (params.run_args_oap) {
        FASTP.out.trimmed_reads
            .map { meta, reads -> tuple(meta.id, reads) }  // Chuyển đổi meta map sang sample_id
            .set { args_oap_input_ch }
        ARGS_OAP(args_oap_input_ch)
    }
}