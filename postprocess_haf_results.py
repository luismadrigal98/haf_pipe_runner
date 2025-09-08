#!/usr/bin/env python3
"""
Standalone Post-processing Script for HAF-pipe Results

This script can be run independently to consolidate HAF-pipe results
from multiple BAM files into unified outputs.

Usage:
    python3 postprocess_haf_results.py --output_dir /path/to/haf_pipe_output
    python3 postprocess_haf_results.py --output_dir /path/to/haf_pipe_output --consolidated_dir my_results

@author: Luis Javier Madrigal-Roca
@date: 2025-09-08
"""

import argparse
import sys
import os
import logging
from pathlib import Path

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

# Import post-processing utilities
from src.postprocessing_utilities import (
    run_full_postprocessing, 
    generate_summary_report,
    collect_allele_counts_files,
    collect_numeric_files
)

def main():
    parser = argparse.ArgumentParser(
        description="Post-process HAF-pipe results to create consolidated outputs",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Basic post-processing
  python3 postprocess_haf_results.py --output_dir haf_pipe_output
  
  # Custom consolidated directory name
  python3 postprocess_haf_results.py --output_dir haf_pipe_output --consolidated_dir final_results
  
  # Just check what files would be processed
  python3 postprocess_haf_results.py --output_dir haf_pipe_output --dry_run
        """
    )
    
    parser.add_argument('--output_dir', '-o', type=str, required=True,
                       help='Directory containing HAF-pipe output subdirectories')
    parser.add_argument('--consolidated_dir', type=str, default='consolidated_results',
                       help='Directory name for consolidated results (default: consolidated_results)')
    parser.add_argument('--dry_run', action='store_true',
                       help='Show what files would be processed without actually processing them')
    parser.add_argument('--verbose', '-v', action='store_true',
                       help='Enable verbose logging')
    
    args = parser.parse_args()
    
    # Set logging level
    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)
    
    # Check if output directory exists
    if not os.path.exists(args.output_dir):
        logger.error(f"Output directory '{args.output_dir}' does not exist")
        sys.exit(1)
    
    # Check what files are available
    logger.info(f"Scanning for HAF-pipe results in: {args.output_dir}")
    
    allele_files = collect_allele_counts_files(args.output_dir)
    numeric_files = collect_numeric_files(args.output_dir)
    
    if not allele_files:
        logger.error("No .alleleCts files found in the output directory")
        logger.info("Expected structure: output_dir/BAM_NAME_output/temp_BAM_NAME/*.alleleCts")
        sys.exit(1)
    
    logger.info(f"Found {len(allele_files)} allele count files")
    logger.info(f"Found {len(numeric_files)} numeric files")
    
    if args.dry_run:
        logger.info("\nDry run - files that would be processed:")
        logger.info("\nAllele count files:")
        for file_path, bam_name in allele_files:
            logger.info(f"  {bam_name}: {file_path}")
        
        logger.info("\nNumeric files:")
        for file_path, bam_name in numeric_files:
            logger.info(f"  {bam_name}: {file_path}")
        
        logger.info(f"\nConsolidated results would be saved to: {os.path.join(args.output_dir, args.consolidated_dir)}")
        return
    
    # Run post-processing
    logger.info("Starting post-processing pipeline...")
    
    try:
        run_full_postprocessing(args.output_dir, args.consolidated_dir)
        generate_summary_report(args.output_dir, args.consolidated_dir)
        
        # Show final results
        consol_path = os.path.join(args.output_dir, args.consolidated_dir)
        logger.info(f"\nPost-processing completed successfully!")
        logger.info(f"Results saved to: {consol_path}")
        
        # List generated files
        if os.path.exists(consol_path):
            logger.info("\nGenerated files:")
            for file in sorted(os.listdir(consol_path)):
                file_path = os.path.join(consol_path, file)
                if os.path.isfile(file_path):
                    logger.info(f"  {file}")
                else:
                    logger.info(f"  {file}/ (directory)")
        
    except Exception as e:
        logger.error(f"Post-processing failed: {e}")
        if args.verbose:
            import traceback
            traceback.print_exc()
        sys.exit(1)

if __name__ == "__main__":
    main()
