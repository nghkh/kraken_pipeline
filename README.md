# Kraken Pipeline

A Nextflow pipeline for metagenomic analysis using Kraken2 and Bracken for taxonomic classification and abundance estimation.

## Overview

This pipeline processes paired-end metagenomic sequencing data through the following steps:

1. **Kraken2**: Taxonomic classification of reads using a reference database
2. **Bracken**: Abundance estimation at multiple taxonomic levels
3. **Krona**: Visualization of taxonomic distributions
4. **Diversity Analysis**: Calculation of alpha and beta diversity metrics

## Requirements

- [Nextflow](https://www.nextflow.io/) (v20.10.0 or later)
- [Docker](https://www.docker.com/) (optional but recommended)
- [Kraken2](https://ccb.jhu.edu/software/kraken2/) database
- [Python](https://www.python.org/) (v3.6 or later)
- [pip](https://pip.pypa.io/en/stable/) (Python package installer)
- [NumPy](https://numpy.org/) (required for diversity analysis)

## Installation

1. Clone this repository:
   ```bash
   git clone https://github.com/nghkh/kraken_pipeline.git
   cd kraken_pipeline
   ```

2. Install Nextflow (if not already installed):
   ```bash
   curl -s https://get.nextflow.io | bash
   ```

3. Set up a Kraken2 database:
   ```bash
   ./bin/setup_krakendb.sh [database_directory]
   ```

## Usage

```bash
nextflow run main.nf --reads_dir [path/to/reads] --outdir [path/to/output] --krakendb [path/to/krakendb]
```

### Input Data

The pipeline expects paired-end sequencing reads in gzipped FASTQ format. By default, it looks for files matching the pattern `*_{1,2}*.{fq,fastq}.gz` in the specified `reads_dir`.

### Example

```bash
nextflow run main.nf --reads_dir /path/to/reads --outdir results --krakendb /path/to/krakendb
```

## Parameters

| Parameter | Description | Default |
|-----------|-------------|---------|
| `--reads_dir` | Directory containing raw reads | `$baseDir/reads` |
| `--outdir` | Output directory (required) | `null` |
| `--krakendb` | Path to Kraken2 database | `$baseDir/krakendb` |
| `--kraken_threads` | Number of threads for Kraken2 | `8` |
| `--bracken_threads` | Number of threads for Bracken | `8` |
| `--bracken_threshold` | Minimum number of reads required for a taxonomic classification | `10` |
| `--bracken_length` | Read length used for Bracken | `150` |
| `--bracken_level` | Taxonomic levels for Bracken analysis | `["P","C","O","F","G","S"]` |
| `--read_pattern` | File pattern to match read files | `*_{1,2}*.{fq,fastq}.gz` |
| `--enable_docker` | Enable Docker support | `true` |
| `--container_kraken2` | Kraken2 Docker container | `quay.io/biocontainers/kraken2:2.1.2--pl5321h9f5acd7_2` |
| `--container_bracken` | Bracken Docker container | `quay.io/biocontainers/bracken:2.7--py39hc16433a_0` |

## Output

The pipeline generates the following output structure:

```
[outdir]/
├── [sample_id]/
│   ├── kraken2/
│   │   ├── [sample_id].kraken2.output  # Raw Kraken2 output
│   │   └── [sample_id].kraken2.report  # Kraken2 report
│   ├── bracken/
│   │   ├── [sample_id].bracken.[level].output  # Bracken abundance estimates
│   │   └── [sample_id].bracken.[level].report  # Bracken report
│   └── krona/
│       └── [sample_id].krona.html  # Interactive Krona visualization
├── summary/
│   └── [level]/
│       └── summary_[level].tsv  # Combined abundance table
└── diversity/
    ├── alpha/
    │   └── alpha_diversity_[level].tsv  # Alpha diversity metrics
    └── beta/
        ├── beta_diversity_[level].tsv   # Beta diversity distance matrix
        └── beta_diversity_[level].pdf   # PCoA plot visualization
```

## Taxonomic Levels

The pipeline supports the following taxonomic levels for Bracken analysis:
- `P`: Phylum
- `C`: Class
- `O`: Order
- `F`: Family
- `G`: Genus
- `S`: Species

## Diversity Metrics

### Alpha Diversity
- Shannon index
- Simpson index
- Chao1 richness estimator

### Beta Diversity
- Bray-Curtis dissimilarity
- Jaccard distance
- UniFrac distance (if phylogenetic information is available)

## Citations

If you use this pipeline, please cite the following tools:

- **Kraken2**: Wood, D.E., Lu, J. & Langmead, B. Improved metagenomic analysis with Kraken 2. Genome Biol 20, 257 (2019).
- **Bracken**: Lu J, Breitwieser FP, Thielen P, Salzberg SL. Bracken: estimating species abundance in metagenomics data. PeerJ Computer Science 3:e104, 2017.
- **Krona**: Ondov BD, Bergman NH, Phillippy AM. Interactive metagenomic visualization in a Web browser. BMC Bioinformatics. 2011 Sep 30;12:385.
- **Nextflow**: Di Tommaso, P., Chatzou, M., Floden, E. W., Barja, P. P., Palumbo, E., & Notredame, C. (2017). Nextflow enables reproducible computational workflows. Nature Biotechnology, 35(4), 316-319.

## License

This pipeline is licensed under [MIT License](LICENSE).

## Contact

For questions or issues, please submit an issue on GitHub or contact:

Email: khoa.bach01@gmail.com
GitHub: https://github.com/nghkh/kraken_pipeline