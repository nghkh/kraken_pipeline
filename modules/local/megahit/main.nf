#!/usr/bin/env nextflow

process MEGAHIT {
    tag "MEGAHIT on ${meta.id}"
    
    label 'megahit'
    
    publishDir "${params.outdir}/megahit/${meta.id}", mode: 'copy', saveAs: { filename ->
        if (filename.endsWith("debug.log")) return "logs/$filename"
        else return filename
    }
    
    input:
    tuple val(meta), path(reads)
    val threads
    
    output:
    tuple val(meta), path("${meta.id}_assembly"), emit: assembly_dir
    tuple val(meta), path("${meta.id}_assembly/final.contigs.fa"), optional: true, emit: contigs
    tuple val(meta), path("${meta.id}_assembly/*.log"), emit: logs
    path "megahit_debug.log", emit: debug
    
    script:
    def read_args = reads.size() == 2 ? "-1 ${reads[0]} -2 ${reads[1]}" : "--12 ${reads[0]}"
    """
    # Tạo file log debug
    echo "=== MEGAHIT DEBUG LOG ===" > megahit_debug.log
    echo "Date: \$(date)" >> megahit_debug.log
    echo "Sample ID: ${meta.id}" >> megahit_debug.log
    echo "Container: \$HOSTNAME" >> megahit_debug.log
    
    # Kiểm tra phiên bản MEGAHIT
    echo "=== MEGAHIT VERSION ===" >> megahit_debug.log
    megahit --version >> megahit_debug.log 2>&1 || echo "Command not found" >> megahit_debug.log
    
    # Kiểm tra files input
    echo "=== INPUT FILES ===" >> megahit_debug.log
    echo "Number of read files: ${reads.size()}" >> megahit_debug.log
    ls -la ${reads.join(' ')} >> megahit_debug.log 2>&1
    
    # Kiểm tra kích thước files
    for file in ${reads.join(' ')}; do
        echo "File size of \$file: \$(du -h \$file | cut -f1)" >> megahit_debug.log
    done
    
    # Xóa thư mục đầu ra nếu đã tồn tại
    if [ -d "${meta.id}_assembly" ]; then
        rm -rf "${meta.id}_assembly"
        echo "Removed existing output directory" >> megahit_debug.log
    fi
    
    # Chạy MEGAHIT và ghi log
    echo "=== RUNNING MEGAHIT ===" >> megahit_debug.log
    echo "Command: megahit $read_args -t $threads -o ${meta.id}_assembly --out-prefix ${meta.id}" >> megahit_debug.log
    
    megahit $read_args \\
        -t $threads \\
        -o ${meta.id}_assembly \\
        --out-prefix ${meta.id} \\
        --continue \\
        --min-contig-len 200 \\
        2> >(tee -a megahit_debug.log >&2) 
    
    # Kiểm tra mã thoát
    exit_code=\$?
    echo "MEGAHIT exit code: \$exit_code" >> megahit_debug.log
    
    # Kiểm tra thư mục và files đầu ra
    echo "=== OUTPUT FILES ===" >> megahit_debug.log
    if [ -d "${meta.id}_assembly" ]; then
        ls -la ${meta.id}_assembly/ >> megahit_debug.log 2>&1
    else
        echo "Output directory not created" >> megahit_debug.log
    fi
    
    # Kiểm tra kết quả sau khi chạy xong
    if [ ! -f "${meta.id}_assembly/final.contigs.fa" ]; then
        echo "WARNING: final.contigs.fa không được tìm thấy! Tạo file rỗng để đảm bảo luồng workflow tiếp tục" >> megahit_debug.log
        mkdir -p ${meta.id}_assembly
        touch ${meta.id}_assembly/final.contigs.fa
    else
        echo "Success: file final.contigs.fa được tạo thành công" >> megahit_debug.log
        echo "Size: \$(du -h ${meta.id}_assembly/final.contigs.fa | cut -f1)" >> megahit_debug.log
    fi
    """
}