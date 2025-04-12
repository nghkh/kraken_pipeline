// Process for Kraken2 taxonomic classification
process KRAKEN2 {
    tag "$sample_id"
    publishDir "${params.outdir}/$sample_id/kraken2", mode: 'copy'
    
    // Thêm container cho Kraken2
    container params.container_kraken2
    
    input:
    tuple val(sample_id), path(reads)
    
    output:
    tuple val(sample_id), path("${sample_id}.kraken2.report"), emit: report
    path "${sample_id}.kraken2.output"
    
    script:
    """
    kraken2 --db ${params.krakendb} \
        --threads ${params.kraken_threads} \
        --paired \
        --output ${sample_id}.kraken2.output \
        --report ${sample_id}.kraken2.report \
        ${reads[0]} ${reads[1]}
    """
}