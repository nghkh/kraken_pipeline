// Process for converting Kraken reports to Krona format
process KRAKEN2KRONA {
    tag "$sample_id"
    label 'process_low'
    publishDir "${params.outdir}/$sample_id/krona", mode: 'copy'
    
    container params.container_bracken
    
    input:
    tuple val(sample_id), path(kraken_report)
    
    output:
    tuple val(sample_id), path("${sample_id}.krona.txt"), emit: krona_text
    path "versions.yml", emit: versions
    
    script:
    def prefix = task.ext.prefix ?: "${sample_id}"
    def VERSION = '2.1.2' // Version of the kreport2krona conversion script
    """
    python3 ${baseDir}/bin/kreport2krona.py -r ${kraken_report} -o ${prefix}.krona.txt --intermediate-ranks
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        kreport2krona: $VERSION
    END_VERSIONS
    """
}

// Process for generating interactive Krona HTML charts
process KRONA_PLOT {
    tag "$sample_id"
    label 'process_low'
    publishDir "${params.outdir}/$sample_id/krona", mode: 'copy'
    
    container params.container_krona
    
    input:
    tuple val(sample_id), path(krona_text)
    
    output:
    tuple val(sample_id), path("${sample_id}.krona.html"), emit: html
    path "versions.yml", emit: versions
    
    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${sample_id}"
    def VERSION = '2.8' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    # Try to use ktImportText
    if command -v ktImportText &>/dev/null; then
        ktImportText $args -o ${prefix}.krona.html ${krona_text}
    else
        echo "ktImportText not found, trying ktImportTaxonomy"
        if command -v ktImportTaxonomy &>/dev/null; then
            ktImportTaxonomy $args -o ${prefix}.krona.html ${krona_text}
        else
            echo "No Krona import tools found. Checking all available tools in the container:"
            find /usr -type f -executable -name "kt*" 2>/dev/null
            
            # Last resort - try specific paths where KronaTools might be installed
            if [ -f "/usr/local/bin/ktImportText" ]; then
                /usr/local/bin/ktImportText $args -o ${prefix}.krona.html ${krona_text}
            elif [ -f "/opt/conda/bin/ktImportText" ]; then
                /opt/conda/bin/ktImportText $args -o ${prefix}.krona.html ${krona_text}
            else
                echo "Cannot find KronaTools. Creating fallback HTML."
                cat > ${prefix}.krona.html << EOF
<!DOCTYPE html>
<html>
<head><title>Krona fallback for ${prefix}</title></head>
<body>
<h1>Krona visualization could not be generated</h1>
<p>The Krona tools could not be found in the container.</p>
</body>
</html>
EOF
            fi
        fi
    fi
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        krona: $VERSION
    END_VERSIONS
    """
}