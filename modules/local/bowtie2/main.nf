process BOWTIE2_BUILD_PHIX {
    tag "$fasta.name"
    label 'process_medium'
    label 'bowtie2'

    conda (params.enable_conda ? "bioconda::bowtie2=2.5.3" : null)
    container params.container_bowtie2 ?: 'quay.io/biocontainers/bowtie2:2.5.4--he96a11b_5'

    input:
    path fasta
    val threads

    output:
    path "bt2_phix_index_dir" , emit: index_dir // Thư mục chứa Bowtie2 index
    path "versions.yml"       , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def fasta_filename = fasta.name // To avoid problems with special characters in fasta path for tag
    """
    echo "Input FASTA file: ${fasta}"
    ls -l ${fasta}
    
    mkdir -p bt2_phix_index_dir
    echo "Created index directory:"
    ls -l bt2_phix_index_dir
    
    bowtie2-build \\
        --threads ${threads} \\
        ${fasta} \\
        "bt2_phix_index_dir/bt2_index_base"
        
    echo "After index building:"
    ls -l bt2_phix_index_dir/

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bowtie2: \$(bowtie2 --version 2>&1 | head -n 1 | sed 's/^.*bowtie2-align-s version //;s/ .*//' | sed 's/^bowtie2-align-l version //;s/ .*//' | sed 's/^bowtie2 version //;s/ .*//')
    END_VERSIONS
    """
}

process BOWTIE2_REMOVAL_PHIX {
    tag "$meta.id"
    label 'process_high'
    label 'bowtie2'

    conda (params.enable_conda ? "bioconda::bowtie2=2.5.3" : null)
    container params.container_bowtie2 ?: 'quay.io/biocontainers/bowtie2:2.5.4--he96a11b_5'

    input:
    tuple val(meta), path(reads) // reads là một list: [read1, read2]
    path index_dir               // Thư mục chứa Bowtie2 index từ BOWTIE2_BUILD_PHIX
    val threads

    output:
    tuple val(meta), path("*.unmapped*.fastq.gz") , emit: reads
    path  "*.mapped*.read_ids.txt", optional:true , emit: read_ids
    tuple val(meta), path("*.bowtie2.log")        , emit: log
    path "versions.yml"                           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "Index directory contents:"
    ls -l ${index_dir}/
    echo "Input reads:"
    ls -l ${reads[0]} ${reads[1]}
    
    bowtie2 \\
        --threads ${threads} \\
        -x ${index_dir}/bt2_index_base \\
        -1 ${reads[0]} \\
        -2 ${reads[1]} \\
        --un-conc-gz ${prefix}.unmapped.%.fastq.gz \\
        --al-conc-gz ${prefix}.mapped.%.fastq.gz \\
        --very-sensitive-local \\
        -S /dev/null 2> ${prefix}.bowtie2.log
        
    echo "After bowtie2 alignment:"
    ls -l ${prefix}*

    # Tạo file read_ids từ mapped reads nếu có
    if [ -f "${prefix}.mapped.1.fastq.gz" ] && [ -f "${prefix}.mapped.2.fastq.gz" ]; then
        zcat ${prefix}.mapped.1.fastq.gz | awk 'NR%4==1 {print substr(\$1,2)}' > ${prefix}.mapped.read_ids.txt
    fi

    # Đảm bảo file output unmapped tồn tại ngay cả khi không có reads nào không khớp
    if [ ! -f "${prefix}.unmapped.1.fastq.gz" ]; then
        touch "${prefix}.unmapped.1.fastq.gz"
    fi
    if [ ! -f "${prefix}.unmapped.2.fastq.gz" ]; then
        touch "${prefix}.unmapped.2.fastq.gz"
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bowtie2: \$(bowtie2 --version 2>&1 | head -n 1 | sed 's/^.*bowtie2-align-s version //;s/ .*//' | sed 's/^bowtie2-align-l version //;s/ .*//' | sed 's/^bowtie2 version //;s/ .*//')
    END_VERSIONS
    """
}