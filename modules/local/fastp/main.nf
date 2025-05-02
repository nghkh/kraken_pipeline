// Process for fastp quality control and preprocessing
process FASTP {
    tag "$sample_id"
    publishDir "${params.outdir}/$sample_id/fastp", mode: 'copy'
    
    container params.container_fastp
    
    input:
    tuple val(sample_id), path(reads)
    
    output:
    tuple val(sample_id), path("${sample_id}_R{1,2}.trimmed.fastq.gz"), emit: trimmed_reads
    path "${sample_id}_fastp.html", emit: html_report
    path "${sample_id}_fastp.json", emit: json_report
    
    script:
    """
    fastp -i ${reads[0]} -I ${reads[1]} \
        -o ${sample_id}_R1.trimmed.fastq.gz -O ${sample_id}_R2.trimmed.fastq.gz \
        -h ${sample_id}_fastp.html -j ${sample_id}_fastp.json \
        --detect_adapter_for_pe \
        --qualified_quality_phred 20 \
        --length_required 50 \
        --thread ${task.cpus}
    """
}