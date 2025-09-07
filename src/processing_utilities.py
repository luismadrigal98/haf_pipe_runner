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
    
    for search_str in search_strings:
        for pattern in chr_patterns:
            match = re.search(pattern, search_str, re.IGNORECASE)
            if match:
                chr_num = match.group(1).zfill(2)  # Zero-pad to 2 digits
                return f"Chr_{chr_num}"
    
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