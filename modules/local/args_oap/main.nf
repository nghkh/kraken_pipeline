// ARGs_OAP module for antibiotic resistance gene detection
process ARGS_OAP {
    tag "${sample_id}"
    label 'args_oap'
    
    // Use Docker container instead of conda
    container params.container_args_oap
    
    publishDir "${params.outdir}/args_oap", mode: 'copy'
    
    input:
    tuple val(sample_id), path(reads)
    
    output:
    tuple val(sample_id), path("${sample_id}_args_oap_results"), emit: results
    path "${sample_id}_args_oap_results/summary.txt", emit: summary
    path "${sample_id}_args_oap_summary.log", emit: log
    
    script:
    def args = task.ext.args ?: ''
    def prefix = sample_id
    
    """
    # Create directories for input and results
    mkdir -p ${prefix}_input
    mkdir -p ${prefix}_args_oap_results
    
    # Copy input files to input directory
    cp ${reads.join(' ')} ${prefix}_input/
    
    # Run ARGs-OAP stage one
    args_oap stage_one -i ${prefix}_input -o ${prefix}_args_oap_results -f fq.gz -t ${task.cpus}
    
    # Run ARGs-OAP stage two
    args_oap stage_two -i ${prefix}_args_oap_results -t ${task.cpus}
    
    # Create a summary log
    echo "ARGs-OAP analysis completed for sample: ${prefix}" > ${prefix}_args_oap_summary.log
    echo "Date: \$(date)" >> ${prefix}_args_oap_summary.log
    
    # If ARGs-OAP doesn't produce a summary.txt, we create a simple one
    if [ ! -f "${prefix}_args_oap_results/summary.txt" ]; then
        echo "ARGs-OAP analysis summary for ${prefix}" > ${prefix}_args_oap_results/summary.txt
        echo "Results directory: ${prefix}_args_oap_results" >> ${prefix}_args_oap_results/summary.txt
        echo "Stage one and stage two completed successfully" >> ${prefix}_args_oap_results/summary.txt
    fi
    """
}