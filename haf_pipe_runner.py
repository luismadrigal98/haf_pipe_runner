"""
This Python program is an wrapper for executing a command-line tool called "haf-pipe".

@author: Luis Javier Madrigal-Roca

@date: 2025-09-04

"""

import subprocess
import os
import sys
import argparse
import logging

# Set up logger
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

def extract_chromosome_name(bam_file):
    """Extract chromosome name from BAM file using samtools."""
    try:
        result = subprocess.run(['samtools', 'view', '-H', bam_file], capture_output=True, text=True, check=True)
        for line in result.stdout.splitlines():
            if line.startswith('@SQ'):
                fields = line.split('\t')
                for field in fields:
                    if field.startswith('SN:'):
                        return field[3:]  # Return chromosome name after 'SN:'
    except subprocess.CalledProcessError as e:
        logger.error(f"Error reading BAM file {bam_file}: {e}")
    return None

def main():

    # Parse command-line arguments
    parser = argparse.ArgumentParser(description="Run haf-pipe with specified parameters.")

    parser.add_argument('--input_vcf', '--iv', type=str, required=True, help='Path to the input VCF file.')

    # This script can work with individual bam files, or a directory containing multiple bam files (in such case, the haf_pipeline will be run for each of the individual bam files)
    parser.add_argument('--bam_file_or_dir', '--bfd', type=str, required=True, help='Path to a BAM file or directory containing BAM files.')
    
    parser.add_argument('--SNP_table', '-s', type=str, required=False, help='Name of the SNP table to be generated. Default is input VCF name with .SNP_table.txt suffix.')
    parser.add_argument('--haf_wrapper', type=str, required=True, help='Path to the bash haf-pipe wrapper.')

    # Arguments specific to haf-pipe
    parser.add_argument('--reference_fasta', '-r', type=str, required=True, help='Path to the reference FASTA file.')
    parser.add_argument('--window_size', '-w', type=int, default=20, help='Window size in kb for haf-pipe. Default is 20.')
    parser.add_argument('--nsites', '-n', type=int, default=20, help="Number of sites to consider in the imputation step.")
    parser.add_argument('--chrom-wise', type=bool, default=True, help='Flag to indicate if haf-pipe should be run in chromosome-wise mode.')

    # Overall control
    parser.add_argument('--parallel', action='store_true', help='Enable parallel processing for faster execution.')
    parser.add_argument('--max_workers', type=int, default=None, help='Maximum number of parallel workers (default: auto-detect).')

    # Print help if no arguments are provided
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)
    
    args = parser.parse_args()

    input_vcf = args.input_vcf
    bam_file_or_dir = args.bam_file_or_dir
    SNP_table = args.SNP_table if args.SNP_table else os.path.splitext(os.path.basename(input_vcf))[0] + '.SNP_table.txt'
    haf_wrapper = args.haf_wrapper
    reference_fasta = args.reference_fasta
    window_size = args.window_size
    nsites = args.nsites
    chrom_wise = args.chrom_wise
    parallel = args.parallel
    max_workers = args.max_workers if args.max_workers is not None else os.cpu_count()

    # Check if input files and directories exist
    if not os.path.exists(input_vcf):
        logger.error(f"Input VCF file '{input_vcf}' does not exist.")
        sys.exit(1)

    if not os.path.exists(bam_file_or_dir):
        logger.error(f"BAM file or directory '{bam_file_or_dir}' does not exist.")
        sys.exit(1)
    
    if not os.path.exists(haf_wrapper):
        logger.error(f"haf-pipe wrapper '{haf_wrapper}' does not exist.")
        sys.exit(1)

    if not os.path.exists(reference_fasta):
        logger.error(f"Reference FASTA file '{reference_fasta}' does not exist.")
        sys.exit(1)

    # Construct the haf-pipe command for a BAM file
    if os.path.isdir(bam_file_or_dir):
        bam_files = [os.path.join(bam_file_or_dir, f) for f in os.listdir(bam_file_or_dir) if f.endswith('.bam')]
    else:
        bam_files = [bam_file_or_dir]

    if not parallel:
        for bam_file in bam_files:
            # This script assumes that the chromosome name in the BAM file matches that in the VCF and reference FASTA
            # So, we can use the bam file to extract the chromosome name if needed

            cmd = [
                haf_wrapper,
                '--input_vcf', input_vcf,
                '--bam_file', bam_file,
                '--SNP_table', SNP_table,
                '--reference_fasta', reference_fasta,
                '--window_size', str(window_size),
                '--nsites', str(nsites)
                ]

            if chrom_wise:
                try:
                    chr_name = extract_chromosome_name(bam_file)
                    cmd.append(['--chrom', chr_name])
                except Exception as e:
                    logger.error(f"Failed to extract chromosome name from {bam_file}: {e}")

            logger.info(f"Running haf-pipe command: {' '.join(cmd)}")
            subprocess.run(cmd, check=True)
    else:
        from concurrent.futures import ThreadPoolExecutor, as_completed

        def run_haf_pipe(bam_file):
            cmd = [
                haf_wrapper,
                '--input_vcf', input_vcf,
                '--bam_file', bam_file,
                '--SNP_table', SNP_table,
                '--reference_fasta', reference_fasta,
                '--window_size', str(window_size),
                '--nsites', str(nsites),
            ]

            if chrom_wise:
                try:
                    chr_name = extract_chromosome_name(bam_file)
                    cmd.append(['--chrom', chr_name])
                except Exception as e:
                    logger.error(f"Failed to extract chromosome name from {bam_file}: {e}")

            logger.info(f"Running haf-pipe command: {' '.join(cmd)}")
            subprocess.run(cmd, check=True)

        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            future_to_bam = {executor.submit(run_haf_pipe, bam): bam for bam in bam_files}
            for future in as_completed(future_to_bam):
                bam = future_to_bam[future]
                try:
                    future.result()
                except Exception as exc:
                    logger.error(f"haf-pipe generated an exception for BAM file {bam}: {exc}")

if __name__ == "__main__":
    main()