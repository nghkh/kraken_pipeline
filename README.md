# Nextflow Kraken2 và Bracken Workflow

Workflow phân tích phân loại vi khuẩn sử dụng Kraken2 và Bracken với Nextflow.

## Tổng quan

Pipeline này thực hiện:
1. Phân loại các reads với Kraken2
2. Phân tích ước lượng abundance với Bracken
3. Tạo báo cáo tổng hợp

## Yêu cầu

- Nextflow (≥ 21.10.0)
- Kraken2
- Bracken
- Database Kraken2 đã được tạo sẵn

## Cách sử dụng

### 1. Chuẩn bị thư mục chứa raw reads

Chuẩn bị một thư mục chứa tất cả các file raw reads (paired-end) của các mẫu. File reads phải có định dạng tên như sau:
- `sample_name_R1.fastq.gz` (forward reads)
- `sample_name_R2.fastq.gz` (reverse reads)

Trong đó:
- `sample_name`: Tên mẫu, sẽ được sử dụng làm tên thư mục con trong kết quả đầu ra.

### 2. Chạy workflow

Chạy với thư mục reads và database đã có sẵn:

```bash
nextflow run main.nf --reads_dir /đường/dẫn/đến/thư_mục_reads --krakendb /đường/dẫn/đến/database
```

Hoặc sử dụng file cấu hình có sẵn:

```bash
nextflow run main.nf -c params/kraken_params.config
```

Bạn có thể điều chỉnh pattern để tìm đúng kiểu file read của bạn:

```bash
nextflow run main.nf --reads_dir /đường/dẫn/reads --read_pattern "*_{1,2}.fastq.gz"
```

Bạn có thể điều chỉnh các tham số khác:

```bash
nextflow run main.nf --reads_dir /đường/dẫn/reads --krakendb /đường/dẫn/db --outdir kết_quả --kraken_threads 16 --bracken_level G
```

### 3. Tiếp tục chạy từ lần chạy trước

Nếu quá trình chạy bị gián đoạn hoặc có lỗi, bạn có thể sử dụng tham số `--resume` để tiếp tục từ vị trí bị dừng:

```bash
nextflow run main.nf --reads_dir /đường/dẫn/đến/thư_mục_reads --krakendb /đường/dẫn/đến/database --outdir /đường/dẫn/đến/kết_quả --resume
```

Khi sử dụng tham số `--resume`, Nextflow sẽ:
- Bỏ qua các bước đã hoàn thành thành công
- Chỉ chạy lại những bước bị lỗi hoặc chưa được thực hiện
- Sử dụng lại kết quả tạm thời đã lưu trong thư mục `work/`

### 4. Các tham số chính

| Tham số | Mô tả | Giá trị mặc định |
|---------|-------|-----------------|
| `--reads_dir` | Thư mục chứa tất cả file reads | `$baseDir/reads` |
| `--read_pattern` | Pattern để tìm file reads | `*_R{1,2}*.fastq.gz` |
| `--krakendb` | Đường dẫn đến database Kraken2 | `/path/to/your/existing/krakendb` |
| `--outdir` | Thư mục đầu ra | `$baseDir/results` |
| `--kraken_threads` | Số luồng cho Kraken2 | `8` |
| `--bracken_threads` | Số luồng cho Bracken | `8` |
| `--bracken_threshold` | Ngưỡng phân loại | `10` |
| `--bracken_length` | Độ dài read | `150` |
| `--bracken_level` | Cấp độ phân loại (S=species, G=genus) | `S` |
| `--resume` | Tiếp tục từ lần chạy trước đó | không bắt buộc |

## Kết quả

Kết quả được tổ chức theo cấu trúc:

```
results/
├── sample1/
│   ├── kraken2/
│   │   ├── sample1.kraken2.output
│   │   └── sample1.kraken2.report
│   └── bracken/
│       ├── sample1.bracken.output
│       └── sample1.bracken.report
├── sample2/
│   ├── kraken2/
│   │   ├── sample2.kraken2.output
│   │   └── sample2.kraken2.report
│   └── bracken/
│       ├── sample2.bracken.output
│       └── sample2.bracken.report
└── summary/
    └── bracken_summary.txt
```