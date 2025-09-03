#!/bin/bash

# SRA to Assembly Pipeline Script
# This script downloads FASTQ files from SRA and runs quality control, trimming, assembly, and assessment

set -e  # Exit on any error

# Function to print colored output
print_status() {
    echo -e "\033[1;32m[$(date '+%Y-%m-%d %H:%M:%S')] $1\033[0m"
}

print_error() {
    echo -e "\033[1;31m[$(date '+%Y-%m-%d %H:%M:%S')] ERROR: $1\033[0m"
}

# Function to check if command exists
check_command() {
    if ! command -v "$1" &> /dev/null; then
        print_error "$1 is not installed or not in PATH"
        exit 1
    fi
}

# Function to display usage
usage() {
    echo "Usage: $0 -a <SRA_ACCESSION> [-o <OUTPUT_DIR>] [-t <THREADS>]"
    echo ""
    echo "Options:"
    echo "  -a    SRA accession number (required)"
    echo "  -o    Output directory (default: current directory)"
    echo "  -t    Number of threads (default: 4)"
    echo "  -h    Show this help message"
    echo ""
    echo "Example:"
    echo "  $0 -a SRR1234567 -o ./analysis -t 8"
    exit 1
}

# Default values
OUTPUT_DIR="."
THREADS=4
SRA_ACCESSION=""

# Parse command line arguments
while getopts "a:o:t:h" opt; do
    case $opt in
        a)
            SRA_ACCESSION="$OPTARG"
            ;;
        o)
            OUTPUT_DIR="$OPTARG"
            ;;
        t)
            THREADS="$OPTARG"
            ;;
        h)
            usage
            ;;
        \?)
            echo "Invalid option: -$OPTARG" >&2
            usage
            ;;
    esac
done

# Check if SRA accession is provided
if [ -z "$SRA_ACCESSION" ]; then
    print_error "SRA accession is required"
    usage
fi

# Check if required tools are installed
print_status "Checking required tools..."
check_command "fastq-dump"
check_command "fastqc"
check_command "fastp"
check_command "megahit"
# check_command "quast"

# Create output directory
mkdir -p "$OUTPUT_DIR"
cd "$OUTPUT_DIR"

# Create subdirectories
mkdir -p raw_data qc_reports trimmed_data assembly assembly_qc

print_status "Starting pipeline for SRA accession: $SRA_ACCESSION"
print_status "Output directory: $(pwd)"
print_status "Using $THREADS threads"

# Step 1: Download FASTQ files from SRA
print_status "Step 1: Downloading FASTQ files from SRA..."
fastq-dump --split-files --gzip --outdir raw_data "$SRA_ACCESSION"

# Check if files were downloaded
if [ ! -f "raw_data/${SRA_ACCESSION}_1.fastq.gz" ]; then
    print_error "Failed to download FASTQ files or single-end data detected"
    # Try single-end download
    print_status "Attempting single-end download..."
    fastq-dump --gzip --outdir raw_data "$SRA_ACCESSION"
    
    if [ ! -f "raw_data/${SRA_ACCESSION}.fastq.gz" ]; then
        print_error "Failed to download FASTQ files"
        exit 1
    fi
    
    SINGLE_END=true
    R1_FILE="raw_data/${SRA_ACCESSION}.fastq.gz"
else
    SINGLE_END=false
    R1_FILE="raw_data/${SRA_ACCESSION}_1.fastq.gz"
    R2_FILE="raw_data/${SRA_ACCESSION}_2.fastq.gz"
fi

print_status "Download completed successfully"

# Step 2: Run FastQC on raw data
print_status "Step 2: Running FastQC on raw data..."
if [ "$SINGLE_END" = true ]; then
    fastqc -o qc_reports -t "$THREADS" "$R1_FILE"
else
    fastqc -o qc_reports -t "$THREADS" "$R1_FILE" "$R2_FILE"
fi
print_status "FastQC analysis completed"

# Step 3: Run fastp for quality trimming and filtering
print_status "Step 3: Running fastp for quality trimming..."
if [ "$SINGLE_END" = true ]; then
    fastp \
        -i "$R1_FILE" \
        -o "trimmed_data/${SRA_ACCESSION}_trimmed.fastq.gz" \
        -h "qc_reports/${SRA_ACCESSION}_fastp.html" \
        -j "qc_reports/${SRA_ACCESSION}_fastp.json" \
        --thread "$THREADS" \
        --qualified_quality_phred 20 \
        --unqualified_percent_limit 40 \
        --length_required 50
    
    TRIMMED_R1="trimmed_data/${SRA_ACCESSION}_trimmed.fastq.gz"
else
    fastp \
        -i "$R1_FILE" \
        -I "$R2_FILE" \
        -o "trimmed_data/${SRA_ACCESSION}_1_trimmed.fastq.gz" \
        -O "trimmed_data/${SRA_ACCESSION}_2_trimmed.fastq.gz" \
        -h "qc_reports/${SRA_ACCESSION}_fastp.html" \
        -j "qc_reports/${SRA_ACCESSION}_fastp.json" \
        --thread "$THREADS" \
        --qualified_quality_phred 20 \
        --unqualified_percent_limit 40 \
        --length_required 50
    
    TRIMMED_R1="trimmed_data/${SRA_ACCESSION}_1_trimmed.fastq.gz"
    TRIMMED_R2="trimmed_data/${SRA_ACCESSION}_2_trimmed.fastq.gz"
fi
print_status "Quality trimming completed"

# Step 4: Run FastQC on trimmed data
print_status "Step 4: Running FastQC on trimmed data..."
if [ "$SINGLE_END" = true ]; then
    fastqc -o qc_reports -t "$THREADS" "$TRIMMED_R1"
else
    fastqc -o qc_reports -t "$THREADS" "$TRIMMED_R1" "$TRIMMED_R2"
fi
print_status "Post-trimming FastQC analysis completed"

# Step 5: Run MEGAHIT for assembly
print_status "Step 5: Running MEGAHIT for genome assembly..."
if [ "$SINGLE_END" = true ]; then
    megahit \
        -r "$TRIMMED_R1" \
        -o "assembly/${SRA_ACCESSION}_megahit" \
        -t "$THREADS" \
        --k-list 21,29,39,59 \
        --min-contig-len 200
else
    megahit \
        -1 "$TRIMMED_R1" \
        -2 "$TRIMMED_R2" \
        -o "assembly/${SRA_ACCESSION}_megahit" \
        -t "$THREADS" \
        --k-list 21,29,39,59 \
        --min-contig-len 200
fi

# Copy final contigs to a standard location
cp "assembly/${SRA_ACCESSION}_megahit/final.contigs.fa" "assembly/${SRA_ACCESSION}_contigs.fa"
print_status "Genome assembly completed"

# Step 6: Run QUAST for assembly quality assessment

print_status "Step 6: Running QUAST for assembly quality assessment..."
cd "./quast-5.3.0"
quast.py \
    "../assembly/${SRA_ACCESSION}_contigs.fa" \
    -o "../assembly_qc/${SRA_ACCESSION}_quast" \
    -t "$2" \
    --min-contig 200

print_status "Assembly quality assessment completed"
cd ../output

# Generate summary report
print_status "Generating summary report..."
SUMMARY_FILE="pipeline_summary_${SRA_ACCESSION}.txt"

cat > "$SUMMARY_FILE" << EOF
SRA to Assembly Pipeline Summary
===============================
SRA Accession: $SRA_ACCESSION
Date: $(date)
Threads used: $THREADS

Files generated:
- Raw FASTQ files: raw_data/
- Quality control reports: qc_reports/
- Trimmed FASTQ files: trimmed_data/
- Assembly: assembly/${SRA_ACCESSION}_contigs.fa
- Assembly QC: assembly_qc/${SRA_ACCESSION}_quast/

Pipeline steps completed:
1. ✓ Downloaded FASTQ files from SRA
2. ✓ Quality control with FastQC (raw data)
3. ✓ Quality trimming with fastp
4. ✓ Quality control with FastQC (trimmed data)
5. ✓ Genome assembly with MEGAHIT
6. ✓ Assembly quality assessment with QUAST

Key files to review:
- FastQC reports: qc_reports/*.html
- fastp report: qc_reports/${SRA_ACCESSION}_fastp.html
- Assembly statistics: assembly_qc/${SRA_ACCESSION}_quast/report.html
- Final contigs: assembly/${SRA_ACCESSION}_contigs.fa
EOF

print_status "Pipeline completed successfully!"
print_status "Summary report saved to: $SUMMARY_FILE"
print_status "Check the HTML reports in qc_reports/ and assembly_qc/ directories"

echo ""
echo "Quick assembly stats:"
if command -v seqkit &> /dev/null; then
    seqkit stats "assembly/${SRA_ACCESSION}_contigs.fa"
else
    echo "Install seqkit for detailed sequence statistics"
    echo "Contig file: assembly/${SRA_ACCESSION}_contigs.fa"
    echo "Number of contigs: $(grep -c "^>" "assembly/${SRA_ACCESSION}_contigs.fa")"
fi
