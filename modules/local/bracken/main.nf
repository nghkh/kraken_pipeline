// Process for Bracken abundance estimation
process BRACKEN {
    tag "$sample_id - $level"
    publishDir "${params.outdir}/$sample_id/bracken", mode: 'copy'
    
    container params.container_bracken
    
    input:
    tuple val(sample_id), path(kraken_report), val(level)
    
    output:
    tuple val(sample_id), val(level), path("${sample_id}.bracken.${level}.report"), emit: reports
    
    script:
    """
    # Run Bracken with the output file name
    bracken -d ${params.krakendb} -i ${kraken_report} -o ${sample_id}.bracken.${level}.output -r ${params.bracken_length} -l ${level} -t ${params.bracken_threshold}
    
    # Rename the output file to the expected format for Nextflow
    mv ${sample_id}.bracken.${level}.output ${sample_id}.bracken.${level}.report
    """
}