"""
Utility functions for processing BAM files and running haf-pipe.

Includes functions to extract chromosome names, find BAM directories,
and run haf-pipe on individual BAM files.

"""

import subprocess
import os
import sys
import argparse
import logging
import re
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor, as_completed
from multiprocessing import cpu_count

# Set up logger
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

def extract_chromosome_from_path(bam_file_path):
    """
    Extract chromosome name from BAM file path or directory.
    Works with directory names like Chr_01, Chr_02, etc.
    Or filenames containing chromosome information.
    """
    path = Path(bam_file_path)
    
    # Consolidated chromosome patterns - covers all common formats
    chr_patterns = [
        r'(?:Chr|chr|chromosome)_?(\d+)',  # Chr_01, chr01, chromosome_01, etc.
        r'(?:CHR|CHROM)_?(\d+)',          # CHR01, CHROM_01, etc.
    ]
    
    # Check both directory name and filename
    search_strings = [path.parent.name, path.name]
    
    logger.debug(f"Extracting chromosome from path: {bam_file_path}")
    logger.debug(f"Directory name: {path.parent.name}")
    logger.debug(f"Filename: {path.name}")
    
    for search_str in search_strings:
        logger.debug(f"Searching in: '{search_str}'")
        for pattern in chr_patterns:
            match = re.search(pattern, search_str, re.IGNORECASE)
            if match:
                chr_num = match.group(1).zfill(2)  # Zero-pad to 2 digits
                result = f"Chr_{chr_num}"
                logger.debug(f"Found chromosome match: {result}")
                return result
    
    logger.warning(f"Could not extract chromosome name from path: {bam_file_path}")
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

def group_bam_files_by_chromosome(bam_files):
    """Group BAM files by chromosome for chromosome-wise processing."""
    bam_groups = {}
    
    for bam_file in bam_files:
        chromosome = extract_chromosome_from_path(str(bam_file))
        if chromosome:
            if chromosome not in bam_groups:
                bam_groups[chromosome] = []
            bam_groups[chromosome].append(bam_file)
            logger.debug(f"Added {bam_file.name} to chromosome group {chromosome}")
        else:
            logger.warning(f"Could not determine chromosome for BAM file: {bam_file}")
    
    return bam_groups

def process_bam_files_for_chromosome(bam_files, args, max_workers):
    """Process a list of BAM files for a specific chromosome."""
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
        
        for i, bam_file in enumerate(bam_files, 1):
            logger.info(f"Processing BAM file {i}/{len(bam_files)}: {bam_file.name}")
            result = run_haf_pipe_single(bam_file, args)
            results.append(result)
            
            if not result[1] and not args.continue_on_error:
                logger.error("Stopping due to error (use --continue_on_error to continue)")
                break
    
    return results

def validate_bam_file(bam_file):
    """Validate that a BAM file exists and has a proper index."""
    bam_path = Path(bam_file)
    if not bam_path.exists():
        return False, f"BAM file does not exist: {bam_file}"
    
    # Check for BAM index (.bai or .csi)
    bai_path = bam_path.with_suffix('.bam.bai')
    csi_path = bam_path.with_suffix('.bam.csi')
    
    if not bai_path.exists() and not csi_path.exists():
        logger.warning(f"No index found for BAM file: {bam_file}")
        # Don't fail, just warn - samtools can create index if needed
    
    return True, None

def run_haf_pipe_SNP_table_and_imputation(args):
    """
    Run haf-pipe for SNP table generation and imputation (tasks 1,2).
    
    This function is specifically for tasks involving SNP table generation and imputation. This will
    speed up processing in parallel of several BAMs, avoiding the overhead of multiple SNP table generations.

    @param args: Parsed command-line arguments
    @return: Tuple (SNP_table_file, success, error_message)

    """
    try:
        # Ensure output directory exists
        output_dir = Path(args.output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        
        # Construct full path for SNP table within output directory
        snp_table_path = output_dir / args.SNP_table
        
        cmd = [
            'sh', args.haf_wrapper,
            '--tasks', "1,2",  # Only tasks 1 and 2 for SNP table and imputation
            '--maindir', args.haf_maindir,  # Specify HAF-pipe main directory
            '--vcf', args.input_vcf,
            '--snptable', str(snp_table_path),
            '--impmethod', 'simpute',  # Add imputation method for task 2
            '--nsites', str(args.nsites),
            '--outdir', args.output_dir
        ]
        
        # Add logfile if specified
        if hasattr(args, 'logfile') and args.logfile:
            cmd.extend(['--logfile', args.logfile])
        
        # Add chromosome if specified for chromosome-wise processing
        if args.chrom_wise and hasattr(args, 'target_chromosome') and args.target_chromosome:
            cmd.extend(['--chrom', args.target_chromosome])
            logger.info(f"Using target chromosome for SNP table generation: {args.target_chromosome}")
        elif args.chrom_wise:
            logger.error("Chromosome-wise mode enabled but no target chromosome specified. This should not happen.")
            return (None, False, "No target chromosome specified for chromosome-wise processing")
        else:
            logger.info("Running SNP table generation without chromosome specification (non-chromosome-wise mode)")

        logger.info(f"Generating the SNP table for {args.input_vcf}")
        logger.debug(f"Command: {' '.join(cmd)}")
        
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        
        # Check if SNP table was created
        if not snp_table_path.exists():
            error_msg = f"SNP table was not created: {snp_table_path}"
            logger.error(error_msg)
            return (None, False, error_msg)
        
        else:
            logger.info(f"Successfully processed {args.input_vcf}")
            # Update args to use the full path for subsequent BAM processing
            args.SNP_table = str(snp_table_path)
            return (str(snp_table_path), True, None)
        
    except subprocess.CalledProcessError as e:
        error_msg = f"Error processing {args.input_vcf}: {e}\nStdout: {e.stdout}\nStderr: {e.stderr}"
        logger.error(error_msg)
        return (args.input_vcf, False, error_msg)
    
    except Exception as e:
        error_msg = f"Unexpected error processing {args.input_vcf}: {e}"
        logger.error(error_msg)
        return (args.input_vcf, False, error_msg)

def validate_snp_table_for_bam_processing(snp_table_path):
    """Validate that the SNP table exists and is readable for BAM processing."""
    snp_path = Path(snp_table_path)
    
    if not snp_path.exists():
        return False, f"SNP table not found: {snp_table_path}"
    
    if not snp_path.is_file():
        return False, f"SNP table path is not a file: {snp_table_path}"
    
    try:
        # Try to read the first few lines to validate format
        with open(snp_path, 'r') as f:
            lines = f.readlines()
            if len(lines) < 2:
                return False, f"SNP table appears to be empty or malformed: {snp_table_path}"
            
        logger.info(f"SNP table validated: {snp_table_path} ({len(lines)} lines)")
        return True, None
        
    except Exception as e:
        return False, f"Error reading SNP table {snp_table_path}: {e}"

def run_haf_pipe_single(bam_file, args):
    """Run haf-pipe for a single BAM file (tasks 3,4)."""
    try:
        cmd = [
            'sh', args.haf_wrapper,
            '--tasks', "3,4",  # Only tasks 3 and 4 for haplotype and allele frequencies
            '--maindir', args.haf_maindir,  # Specify HAF-pipe main directory
            '--snptable', args.SNP_table,  # Use the pre-generated SNP table
            '--bamfile', str(bam_file),
            '--refseq', args.reference_fasta,
            '--winsize', str(args.window_size),
            '--encoding', args.encoding,  # Use user-specified encoding
            '--outdir', args.output_dir
        ]
        
        # Add logfile if specified
        if hasattr(args, 'logfile') and args.logfile:
            cmd.extend(['--logfile', args.logfile])

        if args.chrom_wise:
            chr_name = extract_chromosome_from_path(str(bam_file))
            if chr_name:
                cmd.extend(['--chrom', chr_name])
                logger.info(f"Using chromosome: {chr_name} for BAM file: {bam_file}")
            else:
                logger.warning(f"Could not determine chromosome for {bam_file}, proceeding without --chrom flag")

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