// Process for converting Kraken reports to Krona format
process KRAKEN2KRONA {
    tag "$sample_id"
    publishDir "${params.outdir}/$sample_id/krona", mode: 'copy'
    
    container params.container_bracken
    
    input:
    tuple val(sample_id), path(kraken_report)
    
    output:
    tuple val(sample_id), path("${sample_id}.krona.txt"), emit: krona_text
    
    script:
    """
    python3 ${baseDir}/bin/kreport2krona.py -r ${kraken_report} -o ${sample_id}.krona.txt --intermediate-ranks
    """
}