// Process for fastp quality control and preprocessing
process FASTP {
    tag "$meta.id"
    publishDir "${params.outdir}/${meta.id}/fastp", mode: 'copy'
    
    container params.container_fastp
    
    input:
    tuple val(meta), path(reads)
    
    output:
    tuple val(meta), path("${meta.id}_{1,2}.fastp.fastq.gz"), emit: trimmed_reads
    path "${meta.id}.fastp.html", emit: html_report
    path "${meta.id}.fastp.json", emit: json_report
    path "${meta.id}.fastp.log", emit: log_file
    path "versions.yml", emit: versions
    
    script:
    // Trích xuất tên mẫu từ meta map
    def prefix = meta.id
    def save_merged = false // Mặc định không lưu merged file
    def merge_fastq = save_merged ? "-m --merged_out ${prefix}.merged.fastq.gz" : ''
    def args = "--qualified_quality_phred 20 --length_required 50" // Thêm các tham số mặc định
    def adapter_list = ""
    def fail_fastq = ""
    
    """
    [ ! -f  ${prefix}_1.fastq.gz ] && ln -sf ${reads[0]} ${prefix}_1.fastq.gz
    [ ! -f  ${prefix}_2.fastq.gz ] && ln -sf ${reads[1]} ${prefix}_2.fastq.gz
    
    fastp \\
        --in1 ${prefix}_1.fastq.gz \\
        --in2 ${prefix}_2.fastq.gz \\
        --out1 ${prefix}_1.fastp.fastq.gz \\
        --out2 ${prefix}_2.fastp.fastq.gz \\
        --json ${prefix}.fastp.json \\
        --html ${prefix}.fastp.html \\
        $adapter_list \\
        $fail_fastq \\
        $merge_fastq \\
        --thread $task.cpus \\
        --detect_adapter_for_pe \\
        $args \\
        2> ${prefix}.fastp.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastp: \$(fastp --version 2>&1 | sed -e "s/fastp //g")
    END_VERSIONS
    """
}