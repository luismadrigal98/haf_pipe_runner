# haf_pipe_runner

Python-based wrapper to run the haf-pipe described in 

```
"Accurate Allele Frequencies from Ultra-low Coverage Pool-Seq Samples in Evolve-and-Resequence Experiments."
Tilk S, Bergland A, Goodman A, Schmidt P, Petrov D, Greenblum S.
G3: Genes|Genomes|Genetics, 9(12), 4159 LP - 4168, 2019. doi:10.1534/g3.119.400755
```

which relies on HARP.

```
"Maximum Likelihood Estimation of Frequencies of Known Haplotypes from Pooled Sequence Data."
Kessner D, Turner T, Novembre J.
Molecular Biology and Evolution 30 (5): 1145–58, 2013. doi:10.1093/molbev/mst016
```

## Files in this Repository

- `haf_pipe_runner.py` - HAF-pipe wrapper for single directories or individual BAM files
- `README.md` - This documentation file

## Overview

This wrapper currently supports processing BAM files in a single directory or individual BAM files. It includes parallel processing capabilities and automatic chromosome detection from BAM headers.

## Quick Start

### Basic Usage

```bash
# Process all BAM files in a single directory
python3 haf_pipe_runner.py \
    --input_vcf variants.vcf \
    --bam_file_or_dir /path/to/bam/directory \
    --haf_wrapper /path/to/haf_wrapper.sh \
    --reference_fasta /path/to/reference.fa \
    --parallel

# Process a single BAM file
python3 haf_pipe_runner.py \
    --input_vcf variants.vcf \
    --bam_file_or_dir /path/to/sample.bam \
    --haf_wrapper /path/to/haf_wrapper.sh \
    --reference_fasta /path/to/reference.fa
```

### For Multiple Chromosome Directories

If you have multiple directories (each containing BAM files for different chromosomes), you have two options:

#### Option 1: Simple Bash Loop (Recommended)

```bash
# For directories like Chr_01/, Chr_02/, etc.
for dir in /path/to/chromosome/dirs/Chr_*/; do
    echo "Processing directory: $dir"
    python3 haf_pipe_runner.py \
        --input_vcf variants.vcf \
        --bam_file_or_dir "$dir" \
        --haf_wrapper /path/to/haf_wrapper.sh \
        --reference_fasta /path/to/reference.fa \
        --parallel
done
```

#### Option 2: Create a Custom Batch Script

Create your own batch script that calls `haf_pipe_runner.py` for each directory. This gives you more control over error handling, logging, and parallel processing across directories.

## Features

- **HAF-pipe Integration**: Seamless wrapper for the haf-pipe tool
- **Chromosome Detection**: Automatic chromosome name extraction from BAM headers
- **Parallel Processing**: Multi-threaded execution for faster processing of multiple BAM files
- **Flexible Input**: Support for single BAM files or directories containing multiple BAM files
- **Error Handling**: Robust error handling and logging
- **Configurable Parameters**: Customizable window size, number of sites, and other HAF-pipe parameters

## Requirements

- Python 3.6+
- `samtools` (for BAM file processing)
- `haf-pipe` tool and its dependencies
- Your HAF-pipe wrapper script

## Installation

1. Clone this repository:
   ```bash
   git clone https://github.com/luismadrigal98/haf_pipe_runner.git
   cd haf_pipe_runner
   ```

2. Ensure dependencies are installed:
   ```bash
   # Install samtools (varies by system)
   sudo apt-get install samtools  # Ubuntu/Debian
   # or
   conda install samtools         # Conda
   # or
   brew install samtools          # macOS
   ```

3. Make sure your HAF-pipe installation is working and accessible

## Usage Examples

### Single BAM File

```bash
python3 haf_pipe_runner.py \
    --input_vcf my_variants.vcf \
    --bam_file_or_dir sample.bam \
    --haf_wrapper /opt/haf-pipe/wrapper.sh \
    --reference_fasta reference.fa
```

### Directory with Multiple BAM Files

```bash
python3 haf_pipe_runner.py \
    --input_vcf my_variants.vcf \
    --bam_file_or_dir /data/bam_files/ \
    --haf_wrapper /opt/haf-pipe/wrapper.sh \
    --reference_fasta reference.fa \
    --parallel
```

### Processing Multiple Chromosome Directories

If you have separate directories for each chromosome, use a simple bash loop:

## Command Line Options

### Common Options (All Scripts)

- `--input_vcf` / `--iv`: Path to input VCF file (required)
- `--haf_wrapper`: Path to HAF-pipe wrapper script (required)
- `--reference_fasta` / `-r`: Path to reference FASTA file (required)
- `--SNP_table` / `-s`: Output SNP table name (optional)
- `--window_size` / `-w`: Window size in kb (default: 20)
- `--nsites` / `-n`: Number of sites for imputation (default: 20)
- `--parallel`: Enable parallel processing
- `--max_workers`: Maximum number of parallel workers
- `--continue_on_error`: Continue processing if one file/directory fails

### Single Directory Script (`haf_pipe_runner.py`)

- `--bam_file_or_dir` / `--bfd`: BAM file or directory containing BAM files
- `--chrom-wise`: Enable chromosome-wise mode (default: True)

### Batch Script (`batch_haf_runner.py`)

- `base_directory`: Base directory to search for subdirectories (positional)
- `--directories`: Explicit list of directories to process
- `--dir_pattern`: Pattern to filter directory names
- `--haf_runner_script`: Path to haf_pipe_runner.py (default: ./haf_pipe_runner.py)

### Enhanced Script (`haf_pipe_runner_enhanced.py`)

- `--bam_file_or_dir`: Single BAM file or directory (original behavior)
- `--bam_directories`: List of directories containing BAM files
- `--base_directory`: Base directory for auto-discovery
- `--dir_pattern`: Pattern for directory filtering
- `--chrom-wise` / `--no-chrom-wise`: Control chromosome-wise mode

## Performance Tips

1. **Use Parallel Processing**: Always enable `--parallel` for multiple files/directories
2. **Optimize Worker Count**: Use `--max_workers` to control resource usage
3. **Directory Patterns**: Use `--dir_pattern` to avoid processing unwanted directories
4. **Error Recovery**: Use `--continue_on_error` for robust batch processing
5. **Resource Monitoring**: Monitor CPU and memory usage during large batch jobs

## Troubleshooting

### Common Issues

1. **"samtools not found"**: Ensure samtools is installed and in your PATH
2. **"HAF-pipe wrapper not found"**: Check the path to your wrapper script
3. **"No BAM files found"**: Verify directory structure and BAM file extensions
4. **Permission errors**: Ensure read/write permissions for input/output directories

### Debug Mode

Add verbose logging by modifying the logging level:

```python
logging.basicConfig(level=logging.DEBUG)
```

### Testing

Test with a small dataset first:

```bash
# Test with single directory
python3 haf_pipe_runner.py 
    --input_vcf test_variants.vcf 
    --bam_file_or_dir test_bam_dir/ 
    --haf_wrapper wrapper.sh 
    --reference_fasta test_ref.fa

# Test batch processing with dry run (if implemented)
python3 batch_haf_runner.py test_chr_dirs/ 
    --input_vcf test_variants.vcf 
    --haf_wrapper wrapper.sh 
    --reference_fasta test_ref.fa
```

## Contributing

1. Fork the repository
2. Create a feature branch (`git checkout -b feature/amazing-feature`)
3. Commit your changes (`git commit -m 'Add amazing feature'`)
4. Push to the branch (`git push origin feature/amazing-feature`)
5. Open a Pull Request

## License

This project is licensed under the MIT License - see the LICENSE file for details.

## Citation

If you use this wrapper in your research, please cite the original HAF-pipe paper:

```
Tilk S, Bergland A, Goodman A, Schmidt P, Petrov D, Greenblum S.
"Accurate Allele Frequencies from Ultra-low Coverage Pool-Seq Samples in Evolve-and-Resequence Experiments."
G3: Genes|Genomes|Genetics, 9(12), 4159 LP - 4168, 2019. doi:10.1534/g3.119.400755
```

## Support

For questions or issues, please:

1. Check the troubleshooting section above
2. Search existing GitHub issues
3. Create a new issue with detailed information about your problem

---

**Note**: This wrapper assumes that chromosome names in BAM files match those in the VCF and reference FASTA files. Ensure consistent naming across all input files.
