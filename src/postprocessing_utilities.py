"""
Post-processing utilities for HAF-pipe runner results

This module provides functions to collect, consolidate, and combine
HAF-pipe outputs from multiple BAM files into unified results.

@author: Luis Javier Madrigal-Roca
@date: 2025-09-08
"""

import os
import pandas as pd
import numpy as np
from pathlib import Path
import logging
import re
from typing import List, Dict, Tuple, Optional
import glob

logger = logging.getLogger(__name__)

def extract_sample_info(bam_name: str) -> Dict[str, str]:
    """
    Extract sample information (year, replicate, chromosome) from BAM filename.
    
    Expected format: Mguttatus_bcftools_variants_Chr##_YYYY_X.sorted
    where ## is chromosome, YYYY is year, X is replicate (A or B)
    
    Args:
        bam_name: BAM filename without extension
        
    Returns:
        Dictionary with chromosome, year, replicate info
    """
    # Pattern to match the expected filename format
    pattern = r'.*Chr(\d+)_(\d{4})_([AB])\.sorted'
    match = re.match(pattern, bam_name)
    
    if match:
        chromosome = f"Chr_{match.group(1).zfill(2)}"  # Ensure 2-digit format
        year = match.group(2)
        replicate = match.group(3)
        
        return {
            'chromosome': chromosome,
            'year': year,
            'replicate': replicate,
            'sample_id': f"{year}_{replicate}"
        }
    else:
        logger.warning(f"Could not parse sample info from filename: {bam_name}")
        return {
            'chromosome': 'unknown',
            'year': 'unknown',
            'replicate': 'unknown',
            'sample_id': bam_name
        }

def collect_allele_counts_files(output_dir: str) -> List[Tuple[str, str]]:
    """
    Collect all .alleleCts files from temp directories.
    
    Args:
        output_dir: Base output directory containing BAM output subdirectories
        
    Returns:
        List of tuples (file_path, bam_name)
    """
    allele_files = []
    
    # Search for all .alleleCts files in temp directories
    pattern = os.path.join(output_dir, "*_output", "temp_*", "*.alleleCts")
    files = glob.glob(pattern)
    
    for file_path in files:
        # Extract BAM name from directory structure
        # Path format: output_dir/BAM_NAME_output/temp_BAM_NAME/file.alleleCts
        dir_parts = Path(file_path).parts
        bam_output_dir = None
        for part in dir_parts:
            if part.endswith('_output'):
                bam_output_dir = part
                break
        
        if bam_output_dir:
            bam_name = bam_output_dir.replace('_output', '')
            allele_files.append((file_path, bam_name))
        else:
            logger.warning(f"Could not extract BAM name from path: {file_path}")
    
    logger.info(f"Found {len(allele_files)} allele count files")
    return allele_files

def collect_numeric_files(output_dir: str) -> List[Tuple[str, str]]:
    """
    Collect all .numeric files from temp directories.
    
    Args:
        output_dir: Base output directory containing BAM output subdirectories
        
    Returns:
        List of tuples (file_path, bam_name)
    """
    numeric_files = []
    
    # Search for all .numeric files in temp directories
    pattern = os.path.join(output_dir, "*_output", "temp_*", "*.numeric")
    files = glob.glob(pattern)
    
    for file_path in files:
        # Extract BAM name from directory structure
        dir_parts = Path(file_path).parts
        bam_output_dir = None
        for part in dir_parts:
            if part.endswith('_output'):
                bam_output_dir = part
                break
        
        if bam_output_dir:
            bam_name = bam_output_dir.replace('_output', '')
            numeric_files.append((file_path, bam_name))
        else:
            logger.warning(f"Could not extract BAM name from path: {file_path}")
    
    logger.info(f"Found {len(numeric_files)} numeric files")
    return numeric_files

def load_allele_counts(file_path: str, sample_info: Dict[str, str]) -> pd.DataFrame:
    """
    Load and process an allele counts file.
    
    Args:
        file_path: Path to .alleleCts file
        sample_info: Sample information dictionary
        
    Returns:
        DataFrame with allele count data
    """
    try:
        df = pd.read_csv(file_path)
        
        # Add sample information
        for key, value in sample_info.items():
            df[key] = value
            
        # Calculate allele frequency
        df['total_counts'] = df['refCt'] + df['altCt']
        df['alt_freq'] = np.where(df['total_counts'] > 0, 
                                 df['altCt'] / df['total_counts'], 
                                 0.0)
        df['ref_freq'] = np.where(df['total_counts'] > 0, 
                                 df['refCt'] / df['total_counts'], 
                                 0.0)
        
        return df
        
    except Exception as e:
        logger.error(f"Error loading allele counts from {file_path}: {e}")
        return pd.DataFrame()

def combine_allele_counts(output_dir: str, consolidated_dir: str = "consolidated_results") -> str:
    """
    Combine all allele count files into a single consolidated file.
    
    Args:
        output_dir: Base output directory
        consolidated_dir: Directory name for consolidated results
        
    Returns:
        Path to consolidated allele counts file
    """
    # Create consolidated results directory
    consol_path = os.path.join(output_dir, consolidated_dir)
    os.makedirs(consol_path, exist_ok=True)
    
    # Collect all allele count files
    allele_files = collect_allele_counts_files(output_dir)
    
    if not allele_files:
        logger.warning("No allele count files found")
        return None
    
    all_data = []
    
    for file_path, bam_name in allele_files:
        sample_info = extract_sample_info(bam_name)
        df = load_allele_counts(file_path, sample_info)
        
        if not df.empty:
            all_data.append(df)
            logger.info(f"Loaded {len(df)} variants from {bam_name}")
    
    if not all_data:
        logger.error("No valid allele count data found")
        return None
    
    # Combine all data
    combined_df = pd.concat(all_data, ignore_index=True)
    
    # Sort by chromosome, position, year, replicate
    combined_df = combined_df.sort_values(['chromosome', 'pos', 'year', 'replicate'])
    
    # Save consolidated file
    output_file = os.path.join(consol_path, "consolidated_allele_counts.csv")
    combined_df.to_csv(output_file, index=False)
    
    logger.info(f"Consolidated allele counts saved to: {output_file}")
    logger.info(f"Total variants across all samples: {len(combined_df)}")
    
    return output_file

def create_allele_frequency_matrix(output_dir: str, consolidated_dir: str = "consolidated_results") -> str:
    """
    Create a matrix of allele frequencies across years and replicates.
    
    Args:
        output_dir: Base output directory
        consolidated_dir: Directory name for consolidated results
        
    Returns:
        Path to allele frequency matrix file
    """
    # First create consolidated allele counts
    consol_file = combine_allele_counts(output_dir, consolidated_dir)
    
    if not consol_file:
        return None
    
    # Load consolidated data
    df = pd.read_csv(consol_file)
    
    # Create pivot table for alternative allele frequencies
    alt_freq_matrix = df.pivot_table(
        index=['chromosome', 'pos', 'ref', 'alt'],
        columns='sample_id',
        values='alt_freq',
        fill_value=0.0
    )
    
    # Create pivot table for reference allele frequencies
    ref_freq_matrix = df.pivot_table(
        index=['chromosome', 'pos', 'ref', 'alt'],
        columns='sample_id',
        values='ref_freq',
        fill_value=1.0
    )
    
    # Create pivot table for total counts
    count_matrix = df.pivot_table(
        index=['chromosome', 'pos', 'ref', 'alt'],
        columns='sample_id',
        values='total_counts',
        fill_value=0
    )
    
    # Save matrices
    consol_path = os.path.join(output_dir, consolidated_dir)
    
    alt_freq_file = os.path.join(consol_path, "allele_frequency_matrix_alt.csv")
    ref_freq_file = os.path.join(consol_path, "allele_frequency_matrix_ref.csv")
    count_file = os.path.join(consol_path, "allele_count_matrix.csv")
    
    alt_freq_matrix.to_csv(alt_freq_file)
    ref_freq_matrix.to_csv(ref_freq_file)
    count_matrix.to_csv(count_file)
    
    logger.info(f"Allele frequency matrices saved:")
    logger.info(f"  Alternative alleles: {alt_freq_file}")
    logger.info(f"  Reference alleles: {ref_freq_file}")
    logger.info(f"  Read counts: {count_file}")
    
    return alt_freq_file

def create_vcf_like_output(output_dir: str, consolidated_dir: str = "consolidated_results") -> str:
    """
    Create a VCF-like output with allele frequencies as sample columns.
    
    Args:
        output_dir: Base output directory
        consolidated_dir: Directory name for consolidated results
        
    Returns:
        Path to VCF-like output file
    """
    consol_file = combine_allele_counts(output_dir, consolidated_dir)
    
    if not consol_file:
        return None
    
    df = pd.read_csv(consol_file)
    
    # Create VCF-like format
    vcf_data = []
    
    # Group by position to create VCF records
    for (chrom, pos, ref, alt), group in df.groupby(['chromosome', 'pos', 'ref', 'alt']):
        record = {
            'CHROM': chrom,
            'POS': pos,
            'ID': '.',
            'REF': ref,
            'ALT': alt,
            'QUAL': '.',
            'FILTER': 'PASS',
            'INFO': f"AC={group['altCt'].sum()};AN={group['total_counts'].sum()}"
        }
        
        # Add sample columns with allele frequencies
        for _, row in group.iterrows():
            sample_col = f"{row['sample_id']}_AF"
            record[sample_col] = f"{row['alt_freq']:.4f}"
            
            count_col = f"{row['sample_id']}_DP"
            record[count_col] = row['total_counts']
        
        vcf_data.append(record)
    
    # Convert to DataFrame and save
    vcf_df = pd.DataFrame(vcf_data)
    
    # Sort by chromosome and position
    vcf_df = vcf_df.sort_values(['CHROM', 'POS'])
    
    consol_path = os.path.join(output_dir, consolidated_dir)
    output_file = os.path.join(consol_path, "consolidated_variants_vcf_like.tsv")
    
    vcf_df.to_csv(output_file, sep='\t', index=False)
    
    logger.info(f"VCF-like output saved to: {output_file}")
    logger.info(f"Total variant positions: {len(vcf_df)}")
    
    return output_file

def copy_important_files(output_dir: str, consolidated_dir: str = "consolidated_results"):
    """
    Copy important result files to consolidated directory before cleanup.
    
    Args:
        output_dir: Base output directory
        consolidated_dir: Directory name for consolidated results
    """
    consol_path = os.path.join(output_dir, consolidated_dir)
    os.makedirs(consol_path, exist_ok=True)
    
    # Create subdirectory for individual files
    individual_dir = os.path.join(consol_path, "individual_results")
    os.makedirs(individual_dir, exist_ok=True)
    
    # Copy all .alleleCts files
    allele_files = collect_allele_counts_files(output_dir)
    for file_path, bam_name in allele_files:
        dest_file = os.path.join(individual_dir, f"{bam_name}.alleleCts")
        os.system(f"cp '{file_path}' '{dest_file}'")
        logger.info(f"Copied {bam_name}.alleleCts")
    
    # Copy all .numeric files
    numeric_files = collect_numeric_files(output_dir)
    for file_path, bam_name in numeric_files:
        dest_file = os.path.join(individual_dir, f"{bam_name}.numeric")
        os.system(f"cp '{file_path}' '{dest_file}'")
        logger.info(f"Copied {bam_name}.numeric")
    
    logger.info(f"Individual result files copied to: {individual_dir}")

def run_full_postprocessing(output_dir: str, consolidated_dir: str = "consolidated_results"):
    """
    Run complete post-processing pipeline.
    
    Args:
        output_dir: Base output directory containing HAF-pipe results
        consolidated_dir: Directory name for consolidated results
    """
    logger.info("Starting post-processing pipeline...")
    
    # Copy important files first
    copy_important_files(output_dir, consolidated_dir)
    
    # Create consolidated outputs
    consol_file = combine_allele_counts(output_dir, consolidated_dir)
    
    if consol_file:
        # Create frequency matrices
        matrix_file = create_allele_frequency_matrix(output_dir, consolidated_dir)
        
        # Create VCF-like output
        vcf_file = create_vcf_like_output(output_dir, consolidated_dir)
        
        logger.info("Post-processing pipeline completed successfully!")
        logger.info(f"Results available in: {os.path.join(output_dir, consolidated_dir)}")
    else:
        logger.error("Post-processing pipeline failed")

def generate_summary_report(output_dir: str, consolidated_dir: str = "consolidated_results") -> str:
    """
    Generate a summary report of the processing results.
    
    Args:
        output_dir: Base output directory
        consolidated_dir: Directory name for consolidated results
        
    Returns:
        Path to summary report file
    """
    consol_path = os.path.join(output_dir, consolidated_dir)
    
    # Check if consolidated file exists
    consol_file = os.path.join(consol_path, "consolidated_allele_counts.csv")
    
    if not os.path.exists(consol_file):
        logger.warning("No consolidated file found for summary")
        return None
    
    df = pd.read_csv(consol_file)
    
    # Generate summary statistics
    summary = []
    summary.append("HAF-pipe Processing Summary Report")
    summary.append("=" * 40)
    summary.append("")
    
    # Sample information
    samples = df.groupby(['year', 'replicate']).size()
    summary.append("Samples processed:")
    for (year, rep), count in samples.items():
        summary.append(f"  {year}_{rep}: {count} variants")
    summary.append("")
    
    # Chromosome information
    chroms = df.groupby('chromosome').size()
    summary.append("Chromosomes processed:")
    for chrom, count in chroms.items():
        summary.append(f"  {chrom}: {count} variants")
    summary.append("")
    
    # Variant statistics
    total_variants = len(df.groupby(['chromosome', 'pos', 'ref', 'alt']))
    summary.append(f"Total unique variant positions: {total_variants}")
    summary.append(f"Total variant-sample combinations: {len(df)}")
    summary.append("")
    
    # Allele frequency statistics
    alt_freq_stats = df['alt_freq'].describe()
    summary.append("Alternative allele frequency statistics:")
    for stat, value in alt_freq_stats.items():
        summary.append(f"  {stat}: {value:.4f}")
    summary.append("")
    
    # Coverage statistics
    coverage_stats = df['total_counts'].describe()
    summary.append("Coverage depth statistics:")
    for stat, value in coverage_stats.items():
        summary.append(f"  {stat}: {value:.1f}")
    
    # Save summary report
    report_file = os.path.join(consol_path, "processing_summary.txt")
    with open(report_file, 'w') as f:
        f.write('\n'.join(summary))
    
    logger.info(f"Summary report saved to: {report_file}")
    return report_file
