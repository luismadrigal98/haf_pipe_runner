"""
Utilities for processing BAM files and running haf-pipe.

@author: Luis Javier Madrigal-Roca
@date:2025-09-07

"""

import subprocess
import os
import sys
import argparse
import logging
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor, as_completed
from multiprocessing import cpu_count

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