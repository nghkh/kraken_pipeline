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
    errorStrategy 'retry'
    maxRetries 3
    
    container params.container_bracken
    
    input:
    tuple val(level), path(bracken_reports)
    
    output:
    path "alpha_diversity_${level}.txt"
    
    script:
    """
    # Create header for output file
    echo -e "Sample\tTaxonomic_Level\tShannon\tBerger_Parker\tSimpson\tInverse_Simpson\tFisher" > alpha_diversity_${level}.txt
    
    # Process each Bracken report
    for report in ${bracken_reports}; do
        sample=\$(basename \$report | cut -d. -f1)
        echo -n "\${sample}\t${level}\t" >> alpha_diversity_${level}.txt
        
        # Calculate Shannon diversity
        python3 ${baseDir}/bin/alpha_diversity.py -f \$report -a Sh > shannon_output.txt
        shannon=\$(cat shannon_output.txt | grep "Shannon's diversity:" | awk '{print \$NF}')
        if [[ -z "\$shannon" ]]; then shannon="NA"; fi
        echo -n "\${shannon}\t" >> alpha_diversity_${level}.txt
        
        # Calculate Berger-Parker diversity
        python3 ${baseDir}/bin/alpha_diversity.py -f \$report -a BP > bp_output.txt
        bp=\$(cat bp_output.txt | grep "Berger-parker's diversity:" | awk '{print \$NF}')
        if [[ -z "\$bp" ]]; then bp="NA"; fi
        echo -n "\${bp}\t" >> alpha_diversity_${level}.txt
        
        # Calculate Simpson diversity
        python3 ${baseDir}/bin/alpha_diversity.py -f \$report -a Si > simpson_output.txt
        simpson=\$(cat simpson_output.txt | grep "Simpson's index of diversity:" | awk '{print \$NF}')
        if [[ -z "\$simpson" ]]; then simpson="NA"; fi
        echo -n "\${simpson}\t" >> alpha_diversity_${level}.txt
        
        # Calculate Inverse Simpson diversity
        python3 ${baseDir}/bin/alpha_diversity.py -f \$report -a ISi > invsimpson_output.txt
        inv_simpson=\$(cat invsimpson_output.txt | grep "Simpson's Reciprocal Index:" | awk '{print \$NF}')
        if [[ -z "\$inv_simpson" ]]; then inv_simpson="NA"; fi
        echo -n "\${inv_simpson}\t" >> alpha_diversity_${level}.txt
        
        # Calculate Fisher's alpha
        python3 ${baseDir}/bin/alpha_diversity.py -f \$report -a F > fisher_output.txt
        fisher=\$(cat fisher_output.txt | grep "Fisher's index:" | awk '{print \$NF}')
        if [[ -z "\$fisher" ]]; then fisher="NA"; fi
        echo "\${fisher}" >> alpha_diversity_${level}.txt
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
    path "beta_diversity_${level}.txt"
    
    script:
    bracken_files = bracken_reports.join(' ')
    """
    # Ensure numpy is installed
    pip install numpy scipy
    
    # Run beta diversity analysis using Bracken reports with Python script directly
    # Pass each bracken file individually to the script with -i flag
    # Specify the taxonomy level parameter to match the input level
    python3 ${baseDir}/bin/beta_diversity.py -i ${bracken_files} --type bracken --level ${level} > beta_diversity_${level}.txt
    """
}