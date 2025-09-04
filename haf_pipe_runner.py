"""
Enhanced HAF-pipe Runner

This enhanced version can handle both single directories and multiple directories
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


def find_bam_directories(base_path, pattern=None):
    """Find directories containing BAM files."""
    bam_directories = []
    base_path = Path(base_path)
    
    if not base_path.exists():
        logger.error(f"Base path {base_path} does not exist")
        return []
    
    for item in base_path.iterdir():
        if item.is_dir():
            # Check if pattern matches (if provided)
            if pattern and not item.name.lower().startswith(pattern.lower().replace('*', '')):
                continue
                
            # Check if directory contains BAM files
            bam_files = list(item.glob('*.bam'))
            if bam_files:
                bam_directories.append(item)
                logger.info(f"Found {len(bam_files)} BAM files in {item}")
    
    return sorted(bam_directories)


def get_bam_files_from_directories(directories):
    """Get all BAM files from multiple directories."""
    all_bam_files = []
    for directory in directories:
        bam_files = [f for f in Path(directory).glob('*.bam')]
        all_bam_files.extend(bam_files)
        logger.info(f"Directory {directory}: {len(bam_files)} BAM files")
    return all_bam_files


def run_haf_pipe_single(bam_file, args):
    """Run haf-pipe for a single BAM file."""
    try:
        cmd = [
            args.haf_wrapper,
            '--input_vcf', args.input_vcf,
            '--bam_file', str(bam_file),
            '--SNP_table', args.SNP_table,
            '--reference_fasta', args.reference_fasta,
            '--window_size', str(args.window_size),
            '--nsites', str(args.nsites)
        ]

        if args.chrom_wise:
            try:
                chr_name = extract_chromosome_name(str(bam_file))
                if chr_name:
                    cmd.extend(['--chrom', chr_name])
            except Exception as e:
                logger.error(f"Failed to extract chromosome name from {bam_file}: {e}")

        logger.info(f"Running haf-pipe for {bam_file}")
        logger.debug(f"Command: {' '.join(cmd)}")
        
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        logger.info(f"Successfully processed {bam_file}")
        return (bam_file, True, None)
        
    except subprocess.CalledProcessError as e:
        error_msg = f"Error processing {bam_file}: {e}\nStdout: {e.stdout}\nStderr: {e.stderr}"
        logger.error(error_msg)
        return (bam_file, False, error_msg)
    except Exception as e:
        error_msg = f"Unexpected error processing {bam_file}: {e}"
        logger.error(error_msg)
        return (bam_file, False, error_msg)


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

    # Arguments specific to haf-pipe
    parser.add_argument('--reference_fasta', '-r', type=str, required=True, 
                       help='Path to the reference FASTA file.')
    parser.add_argument('--window_size', '-w', type=int, default=20, 
                       help='Window size in kb for haf-pipe. Default is 20.')
    parser.add_argument('--nsites', '-n', type=int, default=20, 
                       help="Number of sites to consider in the imputation step.")
    parser.add_argument('--chrom-wise', dest='chrom_wise', action='store_true', default=True,
                       help='Flag to indicate if haf-pipe should be run in chromosome-wise mode (default: True).')
    parser.add_argument('--no-chrom-wise', dest='chrom_wise', action='store_false',
                       help='Disable chromosome-wise mode.')

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

    # Set up parallel processing parameters
    max_workers = args.max_workers or min(len(bam_files), cpu_count())

    # Process BAM files
    results = []
    
    if args.parallel and len(bam_files) > 1:
        logger.info(f"Processing {len(bam_files)} BAM files in parallel with {max_workers} workers")
        
        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            future_to_bam = {executor.submit(run_haf_pipe_single, bam, args): bam 
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
        
        for bam_file in bam_files:
            result = run_haf_pipe_single(bam_file, args)
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
    
    if failed and not args.continue_on_error:
        sys.exit(1)


if __name__ == "__main__":
    main()
