// Process for Bracken abundance estimation
process BRACKEN {
    tag "$sample_id - $level"
    publishDir "${params.outdir}/$sample_id/bracken", mode: 'copy'
    
    // Thêm container cho Bracken
    container params.container_bracken
    
    input:
    tuple val(sample_id), path(kraken_report), val(level)
    
    output:
    tuple val(sample_id), val(level), path("${sample_id}.bracken.${level}.report"), emit: reports
    tuple val(sample_id), val(level), path("${sample_id}.bracken.${level}.output"), emit: outputs
    
    script:
    """
    bracken -d ${params.krakendb} \
        -i ${kraken_report} \
        -o ${sample_id}.bracken.${level}.output \
        -w ${sample_id}.bracken.${level}.report \
        -r ${params.bracken_length} \
        -l ${level} \
        -t ${params.bracken_threshold}
    """
}