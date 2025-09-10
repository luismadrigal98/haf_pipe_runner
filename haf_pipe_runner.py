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
from src.postprocessing_utilities import run_full_postprocessing, generate_summary_report

def main():
    # Parse command-line arguments
    parser = argparse.ArgumentParser(
        description="Run haf-pipe with specified parameters. Can handle single directories or multiple directories.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Single directory (original behavior)
  python3 haf_pipe_runner.py --input_vcf variants.vcf --bam_file_or_dir /path/to/bam/dir --haf_wrapper wrapper.sh --reference_fasta ref.fa
  
  # Multiple directories (new feature)
  python3 haf_pipe_runner.py --input_vcf variants.vcf --bam_directories /path/Chr_01 /path/Chr_02 /path/Chr_03 --haf_wrapper wrapper.sh --reference_fasta ref.fa
  
  # Auto-discover directories with pattern and SLURM processing
  python3 haf_pipe_runner.py --input_vcf variants.vcf --base_directory /path/to/chromosomes --dir_pattern "Chr_" --haf_wrapper wrapper.sh --reference_fasta ref.fa --slurm
  
  # With post-processing to consolidate results
  python3 haf_pipe_runner.py --input_vcf variants.vcf --base_directory /path/to/chromosomes --haf_wrapper wrapper.sh --reference_fasta ref.fa --slurm --run_postprocessing
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
    parser.add_argument('--chromosome', type=str, required=False,
                       help='Manually specify target chromosome (e.g., Chr_01). If not provided, will auto-detect from BAM file paths.')
    parser.add_argument('--window_size', '-w', type=int, default=20, 
                       help='Window size in kb for haf-pipe. Default is 20.')
    parser.add_argument('--nsites', '-n', type=int, default=20, 
                       help="Number of sites to consider in the imputation step.")
    parser.add_argument('--encoding', '-e', type=str, default='illumina', choices=['sanger', 'illumina'],
                       help='Base quality encoding in BAM files (default: illumina)')
    parser.add_argument('--imputation_method', '--impmethod', '-i', type=str, default='simpute', 
                       choices=['simpute', 'npute', 'none'],
                       help='Imputation method for HAF-pipe (default: simpute). Options: simpute, npute, none')
    parser.add_argument('--mincalls', '-m', type=int, default=2,
                       help='Keep only sites with at least this many ref|alt calls in SNP table (default: 2)')
    parser.add_argument('--keephets', '-k', action='store_true',
                       help='Keep heterozygous calls as ambiguous bases rather than treating them as missing')
    parser.add_argument('--no-chrom-wise', dest='chrom_wise', action='store_false', default=True,
                       help='Disable chromosome-wise mode (default: enabled).')

    # Overall control
    parser.add_argument('--parallel', action='store_true', 
                       help='Enable parallel processing for faster execution.')
    parser.add_argument('--slurm', action='store_true',
                       help='Submit jobs to SLURM cluster instead of local parallel processing.')
    parser.add_argument('--slurm_partition', type=str, default='sixhour,eeb,kelly,kucg',
                       help='SLURM partition(s) to use (default: sixhour,eeb,kelly,kucg)')
    parser.add_argument('--slurm_email', type=str, default='madrigalrocalj@ku.edu',
                       help='Email for SLURM notifications')
    parser.add_argument('--slurm_time', type=str, default='06:00:00',
                       help='SLURM time limit (default: 06:00:00)')
    parser.add_argument('--slurm_mem', type=str, default='15g',
                       help='SLURM memory per CPU (default: 15g)')
    parser.add_argument('--max_workers', type=int, default=None, 
                       help='Maximum number of parallel workers (default: auto-detect).')
    parser.add_argument('--continue_on_error', action='store_true',
                       help='Continue processing other files if one fails.')
    parser.add_argument('--keep_temp_files', action='store_true',
                       help='Keep all intermediate files for debugging (default: cleanup intermediate, keep final results)')
    parser.add_argument('--keep_all_files', action='store_true',
                       help='Keep all files including intermediate processing files')
    parser.add_argument('--temp_dir_per_bam', type=bool, default=True,
                       help='Create separate temporary directory for each BAM file (default: enabled)')
    parser.add_argument('--run_postprocessing', action='store_true',
                       help='Run post-processing to consolidate results into unified outputs')
    parser.add_argument('--consolidated_dir', type=str, default='consolidated_results',
                       help='Directory name for consolidated post-processing results (default: consolidated_results)')

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
            logger.info(f"Found {len(bam_files)} BAM files in directory: {args.bam_file_or_dir}")
        else:
            bam_files = [Path(args.bam_file_or_dir)]
            logger.info(f"Processing single BAM file: {args.bam_file_or_dir}")
            
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

    # Set up parallel processing parameters
    max_workers = args.max_workers or min(len(bam_files), cpu_count())

    # Process BAM files
    results = []
    
    if args.slurm:
        logger.info(f"Submitting {len(bam_files)} BAM files to SLURM cluster")
        results = submit_slurm_jobs(bam_files, args)
        
    elif args.parallel and len(bam_files) > 1:
        logger.info(f"Processing {len(bam_files)} BAM files in parallel with {max_workers} workers")
        
        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            future_to_bam = {executor.submit(run_haf_pipe_complete, bam, args): bam 
                           for bam in bam_files}
            
            for future in as_completed(future_to_bam):
                bam = future_to_bam[future]
                try:
                    result = future.result()
                    results.append(result)
                except Exception as exc:
                    error_msg = f"BAM file {bam} generated an exception: {exc}"
                    logger.error(error_msg)
                    results.append((bam, False, error_msg))
                    
                    if not args.continue_on_error:
                        logger.error("Stopping due to error (use --continue_on_error to continue)")
                        # Cancel remaining futures
                        for f in future_to_bam:
                            f.cancel()
                        break
    else:
        logger.info(f"Processing {len(bam_files)} BAM files sequentially")
        
        for i, bam_file in enumerate(bam_files, 1):
            logger.info(f"Processing BAM file {i}/{len(bam_files)}: {bam_file.name}")
            result = run_haf_pipe_complete(bam_file, args)
            results.append(result)
            
            if not result[1] and not args.continue_on_error:
                logger.error("Stopping due to error (use --continue_on_error to continue)")
                break

    # Print summary
    successful = [r for r in results if r[1]]
    failed = [r for r in results if not r[1]]
    
    logger.info(f"\nProcessing complete!")
    logger.info(f"Successfully processed: {len(successful)} BAM files")
    if failed:
        logger.error(f"Failed to process: {len(failed)} BAM files")
        for bam_file, _, error in failed:
            logger.error(f"  - {bam_file}: {error}")
    
    # Run post-processing if requested
    if args.run_postprocessing and successful:
        logger.info("\nStarting post-processing to consolidate results...")
        try:
            run_full_postprocessing(args.output_dir, args.consolidated_dir)
            generate_summary_report(args.output_dir, args.consolidated_dir)
            logger.info("Post-processing completed successfully!")
        except Exception as e:
            logger.error(f"Post-processing failed: {e}")
    elif args.run_postprocessing and not successful:
        logger.warning("Skipping post-processing due to no successful results")
    
    if failed and not args.continue_on_error:
        sys.exit(1)

if __name__ == "__main__":
    main()
