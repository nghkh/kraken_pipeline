#!/usr/bin/env nextflow

// Import modules
include { KRAKEN2 } from './modules/local/kraken2/main'
include { BRACKEN } from './modules/local/bracken/main'
include { KRAKEN2KRONA; KRONA_PLOT } from './modules/local/krona/main'
include { GENERATE_SUMMARY; ALPHA_DIVERSITY; BETA_DIVERSITY } from './modules/local/diversity/main'

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

// Thêm tham số Docker
params.enable_docker = true  // Mặc định bật Docker
params.container_kraken2 = 'quay.io/biocontainers/kraken2:2.1.2--pl5321h9f5acd7_2'
params.container_bracken = 'quay.io/biocontainers/bracken:2.7--py39hc16433a_0'
params.container_krona = 'quay.io/biocontainers/krona:2.8--pl5262hdfd78af_2'

// Kiểm tra tham số outdir có được cung cấp hay không
if (params.outdir == null) {
    log.error "Thư mục đầu ra (outdir) không được chỉ định!"
    log.error "Vui lòng chỉ định đường dẫn đầu ra bằng tham số --outdir"
    exit 1
}

// Chuyển đổi bracken_level nếu nó là chuỗi
def bracken_levels = params.bracken_level instanceof String ? 
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
         Bracken tax levels : ${bracken_levels.join(', ')}
         Output directory   : ${params.outdir}
         Using Docker       : ${params.enable_docker}
         Kraken2 container  : ${params.container_kraken2}
         Bracken container  : ${params.container_bracken}
         Krona container    : ${params.container_krona}
         """
         .stripIndent()

// Kiểm tra database Kraken2 có tồn tại không
krakendb_path = file(params.krakendb)
if( !krakendb_path.exists() ) {
    log.error "Kraken2 database not found at: ${params.krakendb}"
    log.error "Please provide the correct path using --krakendb parameter"
    exit 1
}

// Kiểm tra thư mục chứa reads có tồn tại không
reads_dir_path = file(params.reads_dir)
if( !reads_dir_path.exists() ) {
    log.error "Reads directory not found at: ${params.reads_dir}"
    log.error "Please provide the correct path using --reads_dir parameter"
    exit 1
}

// Tạo channel từ các file reads trong thư mục
Channel
    .fromFilePairs("${params.reads_dir}/${params.read_pattern}", checkIfExists: true)
    .ifEmpty { error "No read files found in ${params.reads_dir} with pattern ${params.read_pattern}" }
    .set { read_pairs_ch }

// Workflow definition
workflow {
    // Run Kraken2
    KRAKEN2(read_pairs_ch)
    
    // Tạo channel với tất cả các cấp độ phân loại
    bracken_ch = KRAKEN2.out.report.combine(Channel.fromList(bracken_levels))
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
}