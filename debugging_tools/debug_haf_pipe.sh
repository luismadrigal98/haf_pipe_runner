#!/bin/bash

# Debug version of HAF-pipe run with full output visibility
# This script runs the same commands as our Python wrapper but with debug output

# Apply locale fix
export LC_ALL=C
export LANG=C

# Debug parameters - modify these to match your actual run
BAMFILE="/path/to/your/bamfile.bam"
VCF="/path/to/your/variants.vcf"
REFERENCE="/path/to/your/reference.fasta"
HAF_MAINDIR="/run/user/1000/gvfs/sftp:host=hpc.crc.ku.edu/kuhpc/home/l338m483/bin/haf_pipe/HAFpipe-line"
OUTPUT_DIR="/mnt/1692B2EF92B2D28B/Ongoing_projects/haf_pipe_runner/debug_output"
CHROMOSOME="Chr_01"  # Adjust as needed

# Create output directory
mkdir -p "$OUTPUT_DIR"
cd "$OUTPUT_DIR"

# Generate SNP table name
BAMNAME=$(basename "$BAMFILE" .bam)
SNP_TABLE="${OUTPUT_DIR}/${BAMNAME}_variants.SNP_table.txt"

echo "=== DEBUG HAF-PIPE RUN ==="
echo "BAM file: $BAMFILE"
echo "VCF file: $VCF"
echo "Reference: $REFERENCE"
echo "HAF maindir: $HAF_MAINDIR"
echo "Output dir: $OUTPUT_DIR"
echo "Chromosome: $CHROMOSOME"
echo "SNP table: $SNP_TABLE"
echo "=========================="

# Task 1: Create SNP table
echo -e "\n=== TASK 1: Creating SNP table ==="
${HAF_MAINDIR}/scripts/make_SNPtable_from_vcf.sh \
    -v "$VCF" \
    -c "$CHROMOSOME" \
    -s "$SNP_TABLE" \
    --mincalls 2

echo "Task 1 exit code: $?"
echo "Files created:"
ls -la "${SNP_TABLE}"* 2>/dev/null || echo "No SNP table files found"

# Task 2: Impute SNP table
echo -e "\n=== TASK 2: Imputing SNP table ==="
${HAF_MAINDIR}/scripts/impute_SNPtable.sh "$SNP_TABLE"

echo "Task 2 exit code: $?"
echo "Files created:"
ls -la "${SNP_TABLE}.simpute"* 2>/dev/null || echo "No imputed files found"

# Task 3: Infer haplotype frequencies (MODIFIED FOR DEBUG)
echo -e "\n=== TASK 3: Inferring haplotype frequencies (DEBUG VERSION) ==="

# Check if required files exist
IMPUTED_SNP_TABLE="${SNP_TABLE}.simpute"
echo "Checking required files:"
echo "BAM file exists: $(test -f "$BAMFILE" && echo "YES" || echo "NO")"
echo "BAM index exists: $(test -f "${BAMFILE}.bai" && echo "YES" || echo "NO")"
echo "Reference exists: $(test -f "$REFERENCE" && echo "YES" || echo "NO")"
echo "SNP table exists: $(test -f "$IMPUTED_SNP_TABLE" && echo "YES" || echo "NO")"
echo "SNP table index exists: $(test -f "${IMPUTED_SNP_TABLE}.idx" && echo "YES" || echo "NO")"

# Create index if missing
if [ ! -f "${IMPUTED_SNP_TABLE}.idx" ]; then
    echo "Creating SNP table index..."
    ${HAF_MAINDIR}/scripts/index_snp_table "$IMPUTED_SNP_TABLE" 50000
    echo "Index creation exit code: $?"
fi

# Manual harp commands with full output (NO /dev/null redirection)
echo -e "\n=== MANUAL HARP COMMANDS (FULL OUTPUT) ==="

# Get chromosome info from SNP table
chrom=$(head -1 "$IMPUTED_SNP_TABLE" | cut -f1 -d',')
chrStart=1
chrEnd=$(tail -n1 "$IMPUTED_SNP_TABLE" | cut -d',' -f1)
outfile="$OUTPUT_DIR/$(basename "$BAMFILE")"

echo "Chromosome: $chrom"
echo "Start: $chrStart"
echo "End: $chrEnd"
echo "Output stem: $outfile"

# Window parameters (same as script)
wins=20  # Use smaller window for testing
likewindow=$(( $wins * 10000 ))
likestep=$(( $wins * 5000 ))
freqwindow=$(( $wins * 1000 ))
freqstep=$(( $wins * 100 ))

echo "Like window: $likewindow"
echo "Like step: $likestep"
echo "Freq window: $freqwindow"
echo "Freq step: $freqstep"

# Test with just the first window
start=$chrStart
stop=$(( $start + $likewindow ))
if [ $stop -gt $chrEnd ]; then stop=$chrEnd; fi

echo -e "\n=== Testing first window: ${chrom}:${start}-${stop} ==="

# Run harp like with full output
echo "Running harp like..."
harp like \
    -b "$BAMFILE" \
    --refseq "$REFERENCE" \
    --snps "$IMPUTED_SNP_TABLE" \
    -r "${chrom}:${start}-${stop}" \
    --stem "${outfile}.${chrom}_${start}_${stop}" \
    -I

echo "Harp like exit code: $?"
echo "Files created by harp like:"
ls -la "${outfile}.${chrom}_${start}_${stop}"* 2>/dev/null || echo "No harp like files found"

# Run harp freq with full output
echo -e "\nRunning harp freq..."
harp freq \
    -b "$BAMFILE" \
    --refseq "$REFERENCE" \
    --snps "$IMPUTED_SNP_TABLE" \
    -r "${chrom}:${start}-${stop}" \
    --stem "${outfile}.${chrom}_${start}_${stop}" \
    --window_step $freqstep \
    --window_width $freqwindow \
    -I

echo "Harp freq exit code: $?"
echo "Files created by harp freq:"
ls -la "${outfile}.${chrom}_${start}_${stop}"* 2>/dev/null || echo "No harp freq files found"

echo -e "\n=== DEBUG RUN COMPLETE ==="
echo "Check the output above for any error messages from harp commands"