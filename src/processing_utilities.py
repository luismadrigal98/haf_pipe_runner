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

def run_haf_pipe_complete(bam_file, args):
    """
    Run complete haf-pipe process for a single BAM file (all tasks 1,2,3,4).
    This creates per-BAM isolation with separate temporary directories.
    """
    try:
        # Validate BAM file and index
        is_valid, message = validate_bam_file(bam_file)
        if not is_valid:
            return (bam_file, False, f"BAM validation failed: {message}")
        
        # Create BAM-specific output directory if temp_dir_per_bam is enabled
        if args.temp_dir_per_bam:
            bam_name = Path(bam_file).stem
            # Create single BAM output directory - no temp subdirectory for now
            bam_main_output_dir = Path(args.output_dir) / f"{bam_name}_output"
            bam_main_output_dir.mkdir(parents=True, exist_ok=True)
            
            working_output_dir = str(bam_main_output_dir)
            logger.info(f"Using BAM-specific directory: {working_output_dir}")
        else:
            working_output_dir = args.output_dir
            bam_main_output_dir = Path(args.output_dir)
        
        # Create BAM-specific SNP table name
        bam_name = Path(bam_file).stem
        bam_snp_table = f"{bam_name}_{args.SNP_table}"
        bam_snp_table_path = Path(working_output_dir) / bam_snp_table
        
        # Convert all paths to absolute paths for HAF-pipe
        abs_working_dir = os.path.abspath(working_output_dir)
        abs_bam_file = os.path.abspath(bam_file)
        abs_vcf_file = os.path.abspath(args.input_vcf)
        abs_ref_file = os.path.abspath(args.reference_fasta)
        abs_snp_table = os.path.abspath(bam_snp_table_path)
        abs_logfile = os.path.abspath(os.path.join(working_output_dir, f"HAFpipe-{bam_name}.log"))
        
        cmd = [
            'sh', args.haf_wrapper,
            '--tasks', "1,2,3,4",  # All tasks for complete processing
            '--maindir', args.haf_maindir,
            '--vcf', abs_vcf_file,
            '--bamfile', abs_bam_file,
            '--snptable', abs_snp_table,
            '--refseq', abs_ref_file,
            '--winsize', str(args.window_size),
            '--nsites', str(args.nsites),
            '--encoding', args.encoding,
            '--impmethod', args.imputation_method,  # Correct parameter name
            '--mincalls', str(args.mincalls),
            '--outdir', abs_working_dir,
            '--logfile', abs_logfile
        ]
        
        # Add parallel harp processing if enabled
        if getattr(args, 'use_parallel_harp', False) and getattr(args, 'harp_parallel_jobs', 1) > 1:
            # Check if we have the parallel script
            parallel_script = Path(__file__).parent.parent / "infer_haplotype_freqs_parallel.sh"
            if parallel_script.exists():
                # Set environment variable to use parallel processing
                parallel_script_abs = os.path.abspath(parallel_script)
                logger.info(f"Using parallel harp processing with {args.harp_parallel_jobs} jobs")
                logger.info(f"Parallel script: {parallel_script_abs}")
                # We'll modify the HAF-pipe wrapper to use our parallel script if this env var is set
                env_vars_msg = f"HARP_PARALLEL_SCRIPT={parallel_script_abs}, HARP_PARALLEL_JOBS={args.harp_parallel_jobs}"
                logger.debug(f"Environment variables for parallel processing: {env_vars_msg}")
            else:
                logger.warning(f"Parallel harp script not found at {parallel_script}, using sequential processing")
        
        # Add keephets flag if specified
        if args.keephets:
            cmd.extend(['--keephets'])
        
        # Add chromosome if chromosome-wise processing is enabled
        if args.chrom_wise:
            if args.chromosome:
                # Manual chromosome specified
                cmd.extend(['--chrom', args.chromosome])
                logger.info(f"Using manually specified chromosome: {args.chromosome} for BAM file: {bam_file}")
            else:
                # Auto-detect chromosome from BAM file path
                chr_name = extract_chromosome_from_path(str(bam_file))
                if chr_name:
                    cmd.extend(['--chrom', chr_name])
                    logger.info(f"Using auto-detected chromosome: {chr_name} for BAM file: {bam_file}")
                else:
                    logger.warning(f"Could not determine chromosome for {bam_file}, proceeding without --chrom flag")

        logger.info(f"Running complete haf-pipe for {bam_file}")
        logger.debug(f"Command: {' '.join(cmd)}")
        logger.debug(f"Working directory: {abs_working_dir}")
        
        # Store current directory and change to working directory
        original_cwd = os.getcwd()
        
        try:
            # Change to the working directory so HAF-pipe can find files correctly
            os.chdir(abs_working_dir)
            
            # Set locale environment variables to fix harp locale issues
            env = os.environ.copy()
            env['LC_ALL'] = 'C'
            env['LANG'] = 'C'
            
            # Also set PATH to ensure harp and other tools are found
            if 'PATH' in env:
                env['PATH'] = env['PATH']
            
            logger.info(f"Running HAF-pipe with locale fix (LC_ALL=C, LANG=C)")
            logger.debug(f"Environment variables: LC_ALL={env.get('LC_ALL')}, LANG={env.get('LANG')}")
            
            result = subprocess.run(cmd, capture_output=True, text=True, check=True, env=env)
            
            # Log what files were actually produced
            logger.info(f"HAF-pipe completed for {bam_file}")
            output_files = list(Path(abs_working_dir).glob('*'))
            logger.debug(f"Output files produced: {[f.name for f in output_files if f.is_file()]}")
            
            # Check for expected output files based on HAF-pipe naming convention
            if args.chromosome:
                expected_freqs = Path(abs_working_dir) / f"{Path(abs_bam_file).name}.{args.chromosome}.freqs"
            else:
                chr_name = extract_chromosome_from_path(str(bam_file))
                expected_freqs = Path(abs_working_dir) / f"{Path(abs_bam_file).name}.{chr_name}.freqs"
            
            if expected_freqs.exists():
                logger.info(f"Found expected freqs file: {expected_freqs}")
            else:
                logger.warning(f"Expected freqs file not found: {expected_freqs}")
                # Look for any .freqs files
                freqs_files = list(Path(abs_working_dir).glob('*.freqs'))
                if freqs_files:
                    logger.info(f"Found freqs files: {[f.name for f in freqs_files]}")
                else:
                    logger.warning("No .freqs files found at all!")
                    
        finally:
            # Always restore the original working directory
            os.chdir(original_cwd)
        
        # Files are already in the right place, no need to copy
        
        # Clean up intermediate files if requested
        if args.temp_dir_per_bam and not args.keep_temp_files and not getattr(args, 'keep_all_files', False):
            removed, kept = cleanup_intermediate_files(working_output_dir, bam_name)
            logger.info(f"Cleaned up {len(removed)} intermediate files, preserved {len(kept)} final results")
        
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

def submit_slurm_jobs(bam_files, args):
    """Submit SLURM jobs for each BAM file."""
    import tempfile
    import time
    
    logger.info(f"Preparing to submit {len(bam_files)} SLURM jobs")
    
    # Create SLURM scripts directory
    slurm_dir = Path(args.output_dir) / "slurm_scripts"
    slurm_dir.mkdir(parents=True, exist_ok=True)
    
    job_ids = []
    results = []
    
    for i, bam_file in enumerate(bam_files, 1):
        bam_name = Path(bam_file).stem
        
        # Create SLURM script for this BAM file
        slurm_script = create_slurm_script(bam_file, args, slurm_dir)
        
        # Submit job
        try:
            cmd = ['sbatch', str(slurm_script)]
            result = subprocess.run(cmd, capture_output=True, text=True, check=True)
            
            # Extract job ID from sbatch output
            job_id = result.stdout.strip().split()[-1]
            job_ids.append(job_id)
            
            logger.info(f"Submitted job {job_id} for BAM {i}/{len(bam_files)}: {bam_file.name}")
            
            # Brief pause to avoid overwhelming the scheduler
            time.sleep(0.1)
            
        except subprocess.CalledProcessError as e:
            error_msg = f"Failed to submit SLURM job for {bam_file}: {e.stderr}"
            logger.error(error_msg)
            results.append((bam_file, False, error_msg))
            
            if not args.continue_on_error:
                logger.error("Stopping job submission due to error")
                break
    
    if job_ids:
        logger.info(f"Successfully submitted {len(job_ids)} SLURM jobs: {', '.join(job_ids)}")
        logger.info("Use 'squeue -u $USER' to monitor job status")
        logger.info("Job outputs will be in the slurm_scripts directory")
        
        # Return placeholder results - actual results will be in SLURM output files
        for bam_file in bam_files[:len(job_ids)]:
            results.append((bam_file, True, "Submitted to SLURM"))
    
    return results

def create_slurm_script(bam_file, args, slurm_dir):
    """Create a SLURM script for processing a single BAM file."""
    bam_name = Path(bam_file).stem
    script_path = slurm_dir / f"haf_pipe_{bam_name}.sh"
    
    # Determine chromosome for this BAM
    chromosome_arg = ""
    if args.chrom_wise:
        if args.chromosome:
            chromosome_arg = f"--chromosome {args.chromosome}"
        else:
            chr_name = extract_chromosome_from_path(str(bam_file))
            if chr_name:
                chromosome_arg = f"--chromosome {chr_name}"
    
    script_content = f"""#!/bin/bash
#SBATCH --job-name=haf_pipe_{bam_name}
#SBATCH --output={slurm_dir}/haf_pipe_{bam_name}.out
#SBATCH --error={slurm_dir}/haf_pipe_{bam_name}.err
#SBATCH --partition={args.slurm_partition}
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task={getattr(args, 'harp_parallel_jobs', 1)}
#SBATCH --mem-per-cpu={args.slurm_mem}
#SBATCH --time={args.slurm_time}
#SBATCH --mail-user={args.slurm_email}
#SBATCH --mail-type=FAIL

# Load environment
module load conda
eval "$(conda shell.bash hook)"
conda activate haf-pipe

# Fix locale for harp
export LC_ALL=C
export LANG=C

# Change to working directory
cd {os.getcwd()}

# Run haf-pipe for single BAM file
python {sys.argv[0]} \\
    --input_vcf {args.input_vcf} \\
    --bam_file_or_dir {bam_file} \\
    --haf_wrapper {args.haf_wrapper} \\
    --haf_maindir {args.haf_maindir} \\
    --reference_fasta {args.reference_fasta} \\
    --output_dir {args.output_dir} \\
    --window_size {args.window_size} \\
    --nsites {args.nsites} \\
    --encoding {args.encoding} \\
    --imputation_method {args.imputation_method} \\
    --mincalls {args.mincalls} \\
    {"--keephets" if args.keephets else ""} \\
    --logfile HAFpipe-{bam_name}.log \\
    {chromosome_arg} \\
    {"--keep_temp_files" if args.keep_temp_files else ""} \\
    {"--no-chrom-wise" if not args.chrom_wise else ""}

echo "Job completed for {bam_file}"
"""
    
    with open(script_path, 'w') as f:
        f.write(script_content)
    
    # Make script executable
    script_path.chmod(0o755)
    
    return script_path

def copy_results_to_main_output(temp_dir, bam_output_dir, bam_name):
    """
    Copy important results from temp directory to BAM output directory.
    This preserves final outputs when temp directories are cleaned up.
    
    Args:
        temp_dir: Path to temp directory (e.g., BAM_NAME_output/temp_BAM_NAME/)
        bam_output_dir: Path to BAM output directory (e.g., BAM_NAME_output/)  
        bam_name: Name of BAM file (without extension)
    """
    import shutil
    
    temp_path = Path(temp_dir)
    output_path = Path(bam_output_dir)
    
    # Ensure output directory exists
    output_path.mkdir(parents=True, exist_ok=True)
    
    # Copy important output files (HAF-pipe specific naming)
    # Priority files that should always be preserved:
    priority_extensions = ['.freqs', '.afSite']
    
    # Other useful files to preserve:
    other_extensions = ['.txt', '.simpute', '.npute', '.numeric', '.alleleCts']
    
    all_extensions = priority_extensions + other_extensions
    
    copied_files = []
    missing_priority = []
    
    for file_path in temp_path.iterdir():
        if file_path.is_file():
            # Check if this is an important file
            if any(file_path.suffix == ext for ext in all_extensions):
                # For priority files, use clean naming
                if file_path.suffix in priority_extensions:
                    dest_name = f"{bam_name}{file_path.suffix}"
                    dest_path = output_path / dest_name
                else:
                    # For other files, keep more descriptive names
                    dest_path = output_path / file_path.name
                
                try:
                    shutil.copy2(file_path, dest_path)
                    copied_files.append(dest_path.name)
                    logger.info(f"Preserved: {dest_path.name}")
                except Exception as e:
                    logger.warning(f"Failed to copy {file_path.name}: {e}")
    
    # Check for missing priority files
    for ext in priority_extensions:
        if not any(f.endswith(ext) for f in copied_files):
            missing_priority.append(ext)
    
    if missing_priority:
        logger.warning(f"Missing important output files for {bam_name}: {missing_priority}")
    
    logger.info(f"Copied {len(copied_files)} important files from temp directory to {output_path}")
    return copied_files

def cleanup_temp_directory(temp_dir):
    """Remove temporary directory and its contents."""
    import shutil
    
    try:
        shutil.rmtree(temp_dir)
        logger.debug(f"Cleaned up temporary directory: {temp_dir}")
    except Exception as e:
        logger.warning(f"Failed to cleanup temp directory {temp_dir}: {e}")

def cleanup_intermediate_files(output_dir, bam_name):
    """
    Remove intermediate files while keeping final results.
    
    Removes:
    - SNP tables (.txt)
    - Imputation files (.simpute, .npute) 
    - Numeric files (.numeric)
    - Allele count files (.alleleCts)
    - Log files (.log)
    
    Keeps:
    - Final frequency files (.freqs)
    - Final allele site files (.afSite)
    """
    output_path = Path(output_dir)
    
    # Files to remove (intermediate processing files)
    cleanup_extensions = ['.txt', '.simpute', '.npute', '.numeric', '.alleleCts', '.log']
    
    # Files to keep (final results)
    keep_extensions = ['.freqs', '.afSite']
    
    removed_files = []
    kept_files = []
    
    for file_path in output_path.iterdir():
        if file_path.is_file():
            if any(file_path.suffix == ext for ext in cleanup_extensions):
                try:
                    file_path.unlink()
                    removed_files.append(file_path.name)
                    logger.debug(f"Removed intermediate file: {file_path.name}")
                except Exception as e:
                    logger.warning(f"Failed to remove {file_path.name}: {e}")
            elif any(file_path.suffix == ext for ext in keep_extensions):
                kept_files.append(file_path.name)
                logger.debug(f"Kept final result: {file_path.name}")
    
    logger.info(f"Cleanup complete for {bam_name}: removed {len(removed_files)} intermediate files, kept {len(kept_files)} final results")
    
    if not kept_files:
        logger.warning(f"No final result files (.freqs, .afSite) found for {bam_name}")
    
    return removed_files, kept_files

def check_bam_completion_status(bam_file, output_dir):
    """
    Check if a BAM file has been successfully processed.
    
    Args:
        bam_file: Path to BAM file
        output_dir: Base output directory
        
    Returns:
        tuple: (is_complete, status_message, expected_files)
    """
    bam_name = Path(bam_file).stem
    bam_output_dir = Path(output_dir) / f"{bam_name}_output"
    
    if not bam_output_dir.exists():
        return False, "Output directory not found", []
    
    # Check for expected output files
    chr_name = extract_chromosome_from_path(str(bam_file))
    if not chr_name:
        return False, "Cannot determine chromosome", []
    
    expected_files = {
        'freqs': bam_output_dir / f"{Path(bam_file).name}.{chr_name}.freqs",
        'afSite': bam_output_dir / f"{Path(bam_file).name}.{chr_name}.afSite",
        'log': bam_output_dir / f"HAFpipe-{bam_name}.log"
    }
    
    existing_files = []
    missing_files = []
    
    for file_type, file_path in expected_files.items():
        if file_path.exists() and file_path.stat().st_size > 0:
            existing_files.append(file_type)
        else:
            missing_files.append(file_type)
    
    if not missing_files:
        return True, "All expected files present", existing_files
    elif existing_files:
        return False, f"Partially complete (missing: {', '.join(missing_files)})", existing_files
    else:
        return False, "No output files found", []

def filter_incomplete_bams(bam_files, output_dir):
    """
    Filter BAM files to only those that need (re)processing.
    
    Args:
        bam_files: List of BAM file paths
        output_dir: Base output directory
        
    Returns:
        tuple: (incomplete_bams, complete_bams, status_summary)
    """
    incomplete_bams = []
    complete_bams = []
    status_summary = []
    
    for bam_file in bam_files:
        is_complete, status, files = check_bam_completion_status(bam_file, output_dir)
        
        status_summary.append({
            'bam': Path(bam_file).name,
            'complete': is_complete,
            'status': status,
            'files': files
        })
        
        if is_complete:
            complete_bams.append(bam_file)
            logger.info(f"✓ {Path(bam_file).name}: Complete")
        else:
            incomplete_bams.append(bam_file)
            logger.info(f"⚠ {Path(bam_file).name}: {status}")
    
    return incomplete_bams, complete_bams, status_summary

def validate_bam_file(bam_file):
    """Validate that a BAM file exists and has a proper index."""
    bam_path = Path(bam_file)
    if not bam_path.exists():
        return False, f"BAM file does not exist: {bam_file}"
    
    # Check for BAM index (.bai) - HAF-pipe expects .bam.bai format
    bai_path = Path(str(bam_path) + '.bai')
    alt_bai_path = bam_path.with_suffix('.bai')  # Some tools create .bai without .bam
    
    if not bai_path.exists():
        if alt_bai_path.exists():
            # Create symlink from .bai to .bam.bai
            try:
                bai_path.symlink_to(alt_bai_path)
                logger.info(f"Created BAM index symlink: {bai_path}")
            except Exception as e:
                logger.warning(f"Could not create BAM index symlink: {e}")
                return False, f"BAM index file missing and could not create symlink: {bam_file}.bai"
        else:
            return False, f"BAM index file missing: {bam_file}.bai (required by HAF-pipe)"
    
    return True, "BAM file and index validated"
    """Validate that a BAM file exists and has a proper index."""
    bam_path = Path(bam_file)
    if not bam_path.exists():
        return False, f"BAM file does not exist: {bam_file}"
    
    # Check for BAM index (.bai) - HAF-pipe expects .bam.bai format
    bai_path = Path(str(bam_path) + '.bai')
    alt_bai_path = bam_path.with_suffix('.bai')  # Some tools create .bai without .bam
    
    if not bai_path.exists():
        if alt_bai_path.exists():
            # Create symlink from .bai to .bam.bai
            try:
                bai_path.symlink_to(alt_bai_path)
                logger.info(f"Created BAM index symlink: {bai_path}")
            except Exception as e:
                logger.warning(f"Could not create BAM index symlink: {e}")
                return False, f"BAM index file missing and could not create symlink: {bam_file}.bai"
        else:
            return False, f"BAM index file missing: {bam_file}.bai (required by HAF-pipe)"
    
    return True, "BAM file and index validated"

def check_haf_pipe_dependencies():
    """Check if HAF-pipe dependencies are available."""
    required_tools = ['harp', 'tabix', 'bgzip', 'Rscript']
    missing_tools = []
    
    for tool in required_tools:
        try:
            subprocess.run(['which', tool], check=True, capture_output=True)
        except subprocess.CalledProcessError:
            missing_tools.append(tool)
    
    if missing_tools:
        return False, f"Missing required tools: {', '.join(missing_tools)}"
    
    return True, "All dependencies found"
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