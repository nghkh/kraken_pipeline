// Process for generating summary report
process GENERATE_SUMMARY {
    tag "$level"
    publishDir "${params.outdir}/summary", mode: 'copy'
    
    // Sử dụng container của Bracken cho quá trình tạo báo cáo tổng hợp
    container params.container_bracken
    
    input:
    tuple val(level), path(bracken_reports)
    
    output:
    path "bracken_summary_${level}.txt"
    
    script:
    """
    echo "Sample\tTaxonomic_Level\tTotal_Reads\tClassified_Reads\tUnclassified_Reads" > bracken_summary_${level}.txt
    
    for report in ${bracken_reports}; do
        sample=\$(basename \$report | cut -d. -f1)
        
        # Kiểm tra xem có thông tin total_reads không
        if grep -q "total_reads" \$report; then
            total=\$(grep -A1 "total_reads" \$report | tail -n1 | awk '{print \$2}')
            # Kiểm tra xem total có phải số hợp lệ không
            if [[ -z "\$total" || "\$total" == *[^0-9]* ]]; then
                echo "\$sample\t${level}\tN/A\tN/A\tN/A" >> bracken_summary_${level}.txt
                continue
            fi
            
            # Kiểm tra xem có thông tin total_classified không
            if grep -q "total_classified" \$report; then
                classified=\$(grep -A1 "total_classified" \$report | tail -n1 | awk '{print \$2}')
                # Kiểm tra xem classified có phải số hợp lệ không
                if [[ -z "\$classified" || "\$classified" == *[^0-9]* ]]; then
                    echo "\$sample\t${level}\t\$total\tN/A\tN/A" >> bracken_summary_${level}.txt
                    continue
                fi
                
                # Tính toán số reads không được phân loại
                unclassified=\$((\$total - \$classified))
                echo "\$sample\t${level}\t\$total\t\$classified\t\$unclassified" >> bracken_summary_${level}.txt
            else
                echo "\$sample\t${level}\t\$total\tN/A\tN/A" >> bracken_summary_${level}.txt
            fi
        else
            echo "\$sample\t${level}\tN/A\tN/A\tN/A" >> bracken_summary_${level}.txt
        fi
    done
    """
}

// Process for calculating Alpha Diversity
process ALPHA_DIVERSITY {
    tag "$level"
    publishDir "${params.outdir}/diversity/alpha", mode: 'copy'
    
    container params.container_bracken
    
    input:
    tuple val(level), path(bracken_reports)
    
    output:
    path "alpha_diversity_${level}.tsv"
    
    script:
    """
    # Create header for output file
    echo -e "Sample\tTaxonomic_Level\tShannon\tBerger_Parker\tSimpson\tInverse_Simpson\tFisher" > alpha_diversity_${level}.tsv
    
    # Process each Bracken report
    for report in ${bracken_reports}; do
        sample=\$(basename \$report | cut -d. -f1)
        
        # Calculate Shannon diversity
        shannon=\$(python ${baseDir}/bin/alpha_diversity.py -f \$report -a Sh)
        shannon_value=\$(echo \$shannon | grep -oP "Shannon's diversity: \\K[0-9.]+" || echo "NA")
        
        # Calculate Berger-Parker diversity
        bp=\$(python ${baseDir}/bin/alpha_diversity.py -f \$report -a BP)
        bp_value=\$(echo \$bp | grep -oP "Berger-parker's diversity: \\K[0-9.]+" || echo "NA")
        
        # Calculate Simpson diversity
        simpson=\$(python ${baseDir}/bin/alpha_diversity.py -f \$report -a Si)
        simpson_value=\$(echo \$simpson | grep -oP "Simpson's index of diversity: \\K[0-9.]+" || echo "NA")
        
        # Calculate Inverse Simpson diversity
        inv_simpson=\$(python ${baseDir}/bin/alpha_diversity.py -f \$report -a ISi)
        inv_simpson_value=\$(echo \$inv_simpson | grep -oP "Simpson's Reciprocal Index: \\K[0-9.]+" || echo "NA")
        
        # Calculate Fisher's alpha (may take longer)
        fisher=\$(python ${baseDir}/bin/alpha_diversity.py -f \$report -a F || echo "Fisher's alpha...loading NA")
        fisher_value=\$(echo \$fisher | grep -oP "Fisher's index: \\K[0-9.]+" || echo "NA")
        
        # Add results to output file
        echo -e "\${sample}\t${level}\t\${shannon_value}\t\${bp_value}\t\${simpson_value}\t\${inv_simpson_value}\t\${fisher_value}" >> alpha_diversity_${level}.tsv
    done
    """
}

// Process for calculating Beta Diversity
process BETA_DIVERSITY {
    tag "$level"
    publishDir "${params.outdir}/diversity/beta", mode: 'copy'
    
    container params.container_bracken
    
    input:
    tuple val(level), path(bracken_reports)
    
    output:
    path "beta_diversity_${level}.tsv"
    
    script:
    """
    # Run beta diversity analysis
    python ${baseDir}/bin/beta_diversity.py --input ${bracken_reports} --type bracken --level ${level} > beta_diversity_${level}.tsv
    """
}