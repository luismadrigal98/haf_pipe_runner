"""
Enhanced HAF-pipe Runner

This enhanced version can handle both single direc    parser.add_argument('--encoding', '-e', type=str, default='sanger', choices=['sanger', 'illumina'],
                       help='Base quality encoding in BAM files (default: sanger)')
    parser.add_argument('--chromosome', type=str, required=False,
                       help='Manually specify target chromosome (e.g., Chr_01). If not provided, will auto-detect from BAM file paths.')
    parser.add_argument('--no-chrom-wise', dest='chrom_wise', action='store_false', default=True,
                       help='Disable chromosome-wise mode (default: enabled).')es and multiple directories
containing chromosome-specific BAM files.

@author: Luis Javier Madrigal-Roca
@date: 2025-09-04
"""

import subprocess
import os
import sys
import argparse
import logging
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor, as_completed
from multiprocessing import cpu_count

# Set up logger
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

# Import utility functions from processing_utilities.py
from src.processing_utilities import *

def main():
    # Parse command-line arguments
    parser = argparse.ArgumentParser(
        description="Run haf-pipe with specified parameters. Can handle single directories or multiple directories.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Single directory (original behavior)
  python3 haf_pipe_runner_enhanced.py --input_vcf variants.vcf --bam_file_or_dir /path/to/bam/dir --haf_wrapper wrapper.sh --reference_fasta ref.fa
  
  # Multiple directories (new feature)
  python3 haf_pipe_runner_enhanced.py --input_vcf variants.vcf --bam_directories /path/Chr_01 /path/Chr_02 /path/Chr_03 --haf_wrapper wrapper.sh --reference_fasta ref.fa
  
  # Auto-discover directories with pattern
  python3 haf_pipe_runner_enhanced.py --input_vcf variants.vcf --base_directory /path/to/chromosomes --dir_pattern "Chr_" --haf_wrapper wrapper.sh --reference_fasta ref.fa --parallel
        """
    )

    parser.add_argument('--input_vcf', '--iv', type=str, required=True, 
                       help='Path to the input VCF file.')

    # BAM file/directory specification (mutually exclusive)
    bam_group = parser.add_mutually_exclusive_group(required=True)
    bam_group.add_argument('--bam_file_or_dir', '--bfd', type=str, 
                          help='Path to a BAM file or single directory containing BAM files.')
    bam_group.add_argument('--bam_directories', nargs='+', 
                          help='List of directories containing BAM files (one per chromosome).')
    bam_group.add_argument('--base_directory', type=str,
                          help='Base directory to search for subdirectories containing BAM files.')
    
    # Directory pattern for auto-discovery
    parser.add_argument('--dir_pattern', type=str,
                       help='Pattern to filter directory names when using --base_directory (e.g., "Chr_", "chromosome_")')
    
    parser.add_argument('--SNP_table', '-s', type=str, required=False, 
                       help='Name of the SNP table to be generated. Default is input VCF name with .SNP_table.txt suffix.')
    parser.add_argument('--haf_wrapper', type=str, required=True, 
                       help='Path to the bash haf-pipe wrapper.')
    parser.add_argument('--haf_maindir', type=str, required=False,
                       help='Directory in which HAF-pipe is located. Default: directory of the haf_wrapper script.')
    parser.add_argument('--logfile', type=str, required=False,
                       help='Name of file to write HAF-pipe log to. Default: auto-generated timestamp-based name.')
    
    # Output directory
    parser.add_argument('--output_dir', '-o', type=str, default='haf_pipe_output', 
                       help='Directory to store haf-pipe output files. Default is "haf_pipe_output".')

    # Arguments specific to haf-pipe
    parser.add_argument('--reference_fasta', '-r', type=str, required=True, 
                       help='Path to the reference FASTA file.')
    parser.add_argument('--window_size', '-w', type=int, default=20, 
                       help='Window size in kb for haf-pipe. Default is 20.')
    parser.add_argument('--nsites', '-n', type=int, default=20, 
                       help="Number of sites to consider in the imputation step.")
    parser.add_argument('--encoding', '-e', type=str, default='illumina', choices=['sanger', 'illumina'],
                       help='Base quality encoding in BAM files (default: illumina)')
    parser.add_argument('--no-chrom-wise', dest='chrom_wise', action='store_false', default=True,
                       help='Disable chromosome-wise mode (default: enabled).')

    # Overall control
    parser.add_argument('--parallel', action='store_true', 
                       help='Enable parallel processing for faster execution.')
    parser.add_argument('--max_workers', type=int, default=None, 
                       help='Maximum number of parallel workers (default: auto-detect).')
    parser.add_argument('--continue_on_error', action='store_true',
                       help='Continue processing other files if one fails.')

    # Print help if no arguments are provided
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)
    
    args = parser.parse_args()

    # Set default SNP_table name
    if not args.SNP_table:
        args.SNP_table = os.path.splitext(os.path.basename(args.input_vcf))[0] + '.SNP_table.txt'

    # Set default HAF maindir if not provided
    if not args.haf_maindir:
        args.haf_maindir = os.path.dirname(os.path.abspath(args.haf_wrapper))
        logger.info(f"Using HAF-pipe main directory: {args.haf_maindir}")

    # Check if input files exist
    if not os.path.exists(args.input_vcf):
        logger.error(f"Input VCF file '{args.input_vcf}' does not exist.")
        sys.exit(1)

    if not os.path.exists(args.haf_wrapper):
        logger.error(f"haf-pipe wrapper '{args.haf_wrapper}' does not exist.")
        sys.exit(1)

    if not os.path.exists(args.reference_fasta):
        logger.error(f"Reference FASTA file '{args.reference_fasta}' does not exist.")
        sys.exit(1)

    if not os.path.exists(args.haf_maindir):
        logger.error(f"HAF-pipe main directory '{args.haf_maindir}' does not exist.")
        sys.exit(1)

    # Get BAM files based on input method
    if args.bam_file_or_dir:
        # Original single directory/file behavior
        if not os.path.exists(args.bam_file_or_dir):
            logger.error(f"BAM file or directory '{args.bam_file_or_dir}' does not exist.")
            sys.exit(1)
        
        if os.path.isdir(args.bam_file_or_dir):
            bam_files = [f for f in Path(args.bam_file_or_dir).glob('*.bam')]
        else:
            bam_files = [Path(args.bam_file_or_dir)]
            
    elif args.bam_directories:
        # Multiple explicit directories
        directories = []
        for d in args.bam_directories:
            if not os.path.exists(d):
                logger.error(f"Directory '{d}' does not exist.")
                if not args.continue_on_error:
                    sys.exit(1)
                continue
            directories.append(Path(d))
        bam_files = get_bam_files_from_directories(directories)
        
    elif args.base_directory:
        # Auto-discover directories
        if not os.path.exists(args.base_directory):
            logger.error(f"Base directory '{args.base_directory}' does not exist.")
            sys.exit(1)
        
        directories = find_bam_directories(args.base_directory, args.dir_pattern)
        if not directories:
            logger.error("No directories with BAM files found.")
            sys.exit(1)
        bam_files = get_bam_files_from_directories(directories)

    if not bam_files:
        logger.error("No BAM files found to process.")
        sys.exit(1)

    logger.info(f"Found {len(bam_files)} BAM files to process")

    # Group BAM files by chromosome for chromosome-wise processing
    if args.chrom_wise:
        if args.chromosome:
            # Manual chromosome specified - all BAM files belong to this chromosome
            bam_groups = {args.chromosome: bam_files}
            logger.info(f"Using manually specified chromosome {args.chromosome} for all BAM files")
        else:
            # Auto-detect chromosomes from BAM file paths
            bam_groups = group_bam_files_by_chromosome(bam_files)
            if not bam_groups:
                logger.error("No valid chromosome groups found. Disabling chromosome-wise mode.")
                args.chrom_wise = False
            else:
                logger.info(f"Grouped BAM files into {len(bam_groups)} chromosome groups: {list(bam_groups.keys())}")
    
    if not args.chrom_wise:
        # For non-chromosome-wise processing, treat all BAM files as one group
        bam_groups = {"all": bam_files}

    # Set up parallel processing parameters
    max_workers = args.max_workers or min(len(bam_files), cpu_count())

    # Process each chromosome group
    all_results = []
    
    for chromosome, chromosome_bam_files in bam_groups.items():
        logger.info(f"\n=== Processing chromosome group: {chromosome} ({len(chromosome_bam_files)} BAM files) ===")
        
        # Create chromosome-specific SNP table
        if args.chrom_wise and chromosome != "all":
            # Set target chromosome for this group
            args.target_chromosome = chromosome
            snp_table_name = f"{os.path.splitext(args.SNP_table)[0]}_{chromosome}.txt"
        else:
            snp_table_name = args.SNP_table
            args.target_chromosome = None
        
        # Store original SNP table name and set chromosome-specific one
        original_snp_table = args.SNP_table
        args.SNP_table = snp_table_name
        
        # Create SNP table for this chromosome
        snp_table, success, error = run_haf_pipe_SNP_table_and_imputation(args)
        if not success:
            logger.error(f"Failed to create SNP table for chromosome {chromosome}: {error}")
            if not args.continue_on_error:
                sys.exit(1)
            continue
        else:
            logger.info(f"SNP table created successfully for chromosome {chromosome}: {snp_table}")
            
        # Validate SNP table before processing BAM files
        valid, validation_error = validate_snp_table_for_bam_processing(snp_table)
        if not valid:
            logger.error(f"SNP table validation failed for chromosome {chromosome}: {validation_error}")
            if not args.continue_on_error:
                sys.exit(1)
            continue

        # Process BAM files for this chromosome
        chromosome_results = process_bam_files_for_chromosome(chromosome_bam_files, args, max_workers)
        all_results.extend(chromosome_results)
        
        # Restore original SNP table name for next iteration
        args.SNP_table = original_snp_table

    # Print summary
    successful = [r for r in all_results if r[1]]
    failed = [r for r in all_results if not r[1]]
    
    logger.info(f"\nProcessing complete!")
    logger.info(f"Successfully processed: {len(successful)} BAM files")
    if failed:
        logger.error(f"Failed to process: {len(failed)} BAM files")
        for bam_file, _, error in failed:
            logger.error(f"  - {bam_file}: {error}")
    
    if failed and not args.continue_on_error:
        sys.exit(1)

if __name__ == "__main__":
    main()
