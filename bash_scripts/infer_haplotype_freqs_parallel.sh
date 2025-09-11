#!/bin/bash
# Parallel version of infer_haplotype_freqs.sh using GNU parallel or xargs
# This processes multiple harp windows simultaneously

usage() {
    echo "usage: infer_haplotype_freqs_parallel.sh"
    echo "        [ -j --jobs ]       number of parallel jobs (default: 4)"
    echo "        [ -o --outdir ]     output directory"
    echo "        [ -d --maindir ]    HAF-pipe main directory"
    echo "        [ -s --snptable ]   SNP table file"
    echo "        [ -m --method ]     imputation method"
    echo "        [ -b --bamfile ]    BAM file"
    echo "        [ -r --refseq ]     reference fasta"
    echo "        [ -e --encoding ]   base quality encoding (illumina/sanger)"
    echo "        [ -w --winsize ]    window size in kb (default: 20)"
}

# Set defaults
maindir=$(dirname "$0")/..
wins=20
method=''
encoding="illumina"
jobs=4  # Default to 4 parallel jobs

# Parse arguments
while [ "$1" != "" ]; do
    case $1 in
        -j | --jobs )           shift; jobs=$1 ;;
        -d | --maindir )        shift; maindir=$1 ;;
        -o | --outdir )         shift; outdir=$1 ;;
        -s | --snptable )       shift; snptable=$1 ;;
        -m | --method )         shift; method="."$1 ;;
        -b | --bamfile )        shift; bamfile=$1 ;;
        -r | --refseq )         shift; refseq=$1 ;;
        -e | --encoding )       shift; encoding=$1 ;;
        -w | --winsize )        shift; wins=$1 ;;
        -h | --help )           usage; exit 1 ;;
        * ) echo "unknown flag $1"; usage; exit 1 ;;
    esac
    shift
done

# Apply locale fix
export LC_ALL=C
export LANG=C

snptable=${snptable}${method}
if [ -z "$outdir" ]; then outdir=$(dirname $bamfile); fi
outfile=$outdir/$(basename $bamfile)

# Create index if needed
if [ ! -f ${snptable}.idx ]; then  
    ${maindir}/scripts/index_snp_table $snptable 50000
fi

echo "Processing with $jobs parallel jobs"

# Window parameters
likewindow=$(( $wins * 10000 ))
likestep=$(( $wins * 5000 ))
freqwindow=$(( $wins * 1000 ))
freqstep=$(( $wins * 100 ))

# Get chromosome info
chrom=$(head -1 $snptable | cut -f1 -d',')
chrStart=1
chrEnd=$(tail -n1 $snptable | cut -d',' -f1)

echo "Processing chromosome $chrom from $chrStart to $chrEnd with ${wins}kb windows"

# Function to process a single window
process_window() {
    local start=$1
    local chrom=$2
    local chrEnd=$3
    local bamfile=$4
    local refseq=$5
    local snptable=$6
    local outfile=$7
    local encoding=$8
    local likewindow=$9
    local freqwindow=${10}
    local freqstep=${11}
    
    # Apply locale fix in subshell
    export LC_ALL=C
    export LANG=C
    
    # Calculate window end
    local stop=$((start + likewindow))
    if [ $stop -gt $chrEnd ]; then stop=$chrEnd; fi
    
    local window_stem="${outfile}.${chrom}_${start}_${stop}"
    
    # Build encoding flag
    local encoding_flag=""
    if [ "$encoding" = "illumina" ]; then
        encoding_flag="-I"
    fi
    
    echo "Processing window ${chrom}:${start}-${stop}"
    
    # Run harp like
    if harp like \
        -b "$bamfile" \
        --refseq "$refseq" \
        --snps "$snptable" \
        -r "${chrom}:${start}-${stop}" \
        --stem "$window_stem" \
        $encoding_flag >/dev/null 2>&1; then
        
        # Run harp freq
        if harp freq \
            -b "$bamfile" \
            --refseq "$refseq" \
            --snps "$snptable" \
            -r "${chrom}:${start}-${stop}" \
            --stem "$window_stem" \
            --window_step $freqstep \
            --window_width $freqwindow \
            $encoding_flag >/dev/null 2>&1; then
            
            echo "✓ Completed window ${chrom}:${start}-${stop}"
            
            # Clean up intermediate files
            rm -rf "${window_stem}.output" 2>/dev/null
            rm -f "${window_stem}.hlk" 2>/dev/null
            
            return 0
        else
            echo "✗ harp freq failed for window ${chrom}:${start}-${stop}"
            return 1
        fi
    else
        echo "✗ harp like failed for window ${chrom}:${start}-${stop}"
        return 1
    fi
}

# Export function and variables for parallel execution
export -f process_window
export chrom chrEnd bamfile refseq snptable outfile encoding likewindow freqwindow freqstep

# Generate list of window start positions
window_starts=$(seq $chrStart $likestep $chrEnd)
total_windows=$(echo "$window_starts" | wc -l)

echo "Processing $total_windows windows with $jobs parallel jobs..."

# Process windows in parallel using xargs
echo "$window_starts" | xargs -n 1 -P $jobs -I {} bash -c '
    process_window {} "$chrom" "$chrEnd" "$bamfile" "$refseq" "$snptable" "$outfile" "$encoding" "$likewindow" "$freqwindow" "$freqstep"
'

# Check results
successful_freqs=$(ls ${outfile}.${chrom}_*.freqs 2>/dev/null | wc -l)
echo "Successfully processed $successful_freqs out of $total_windows windows"

if [ $successful_freqs -gt 0 ]; then
    echo "Combining frequency files..."
    
    # Combine and sort all frequency files
    cat ${outfile}.${chrom}_*.freqs | tr ' ' '\t' | sort -k2g | tr '\t' ' ' > ${outfile}.${chrom}.freqs
    
    # Clean up individual files
    rm ${outfile}.${chrom}_*.freqs
    
    final_lines=$(wc -l < ${outfile}.${chrom}.freqs)
    echo "Final frequency file created with $final_lines lines: ${outfile}.${chrom}.freqs"
else
    echo "No frequency files were created successfully"
    exit 1
fi
