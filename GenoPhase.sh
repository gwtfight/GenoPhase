#!/bin/bash

# GenoPhase.sh - Wrapper script for GenoPhase pipeline
#Author: Northeast Forestry University 
#Date: January , 2025

set -e  # Exit on error

# Check if required software is installed
command -v python3 >/dev/null 2>&1 || { echo "Python 3 is required but not installed. Aborting."; exit 1; }
command -v minimap2 >/dev/null 2>&1 || { echo "minimap2 is required but not installed. Aborting."; exit 1; }
command -v samtools >/dev/null 2>&1 || { echo "samtools is required but not installed. Aborting."; exit 1; }
command -v hifiasm >/dev/null 2>&1 || { echo "hifiasm is required but not installed. Aborting."; exit 1; }

# Default values
THREADS=16
MODE="full"
MIN_BASE_QUALITY=20
MIN_DEPTH=10
MAX_DEPTH=100
MIN_MAF=0.25
MAX_EDGE_DISTANCE=16000
MIN_SHARED_READS=3

# Parse command line arguments
function show_help {
    echo "Usage: $0 [options]"
    echo ""
    echo "GenoPhase: A tool for haplotype-resolved genome assembly using HiFi reads"
    echo ""
    echo "Options:"
    echo "  --hifi-reads FILE       Input HiFi reads in FASTA/FASTQ format (required)"
    echo "  --reference FILE        Reference genome in FASTA format (optional)"
    echo "  --output-dir DIR        Output directory (required)"
    echo "  --mode MODE             Pipeline mode: full|reference|phase|assemble (default: full)"
    echo "  --threads N             Number of threads to use (default: 16)"
    echo "  --min-base-quality N    Minimum base quality for SNP detection (default: 20)"
    echo "  --min-depth N           Minimum read depth for SNP detection (default: 10)"
    echo "  --max-depth N           Maximum read depth for SNP detection (default: 100)"
    echo "  --min-maf FLOAT         Minimum minor allele frequency for SNP detection (default: 0.25)"
    echo "  --max-edge-distance N   Maximum distance between SNPs to create an edge (default: 16000)"
    echo "  --min-shared-reads N    Minimum number of shared reads to create an edge (default: 3)"
    echo "  --force-remap           Force remapping of reads even if BAM file exists"
    echo "  --help                  Show this help message and exit"
    exit 1
}

# Parse arguments
while [[ $# -gt 0 ]]; do
    case "$1" in
        --hifi-reads)
            HIFI_READS="$2"
            shift 2
            ;;
        --reference)
            REFERENCE="$2"
            shift 2
            ;;
        --output-dir)
            OUTPUT_DIR="$2"
            shift 2
            ;;
        --mode)
            MODE="$2"
            shift 2
            ;;
        --threads)
            THREADS="$2"
            shift 2
            ;;
        --min-base-quality)
            MIN_BASE_QUALITY="$2"
            shift 2
            ;;
        --min-depth)
            MIN_DEPTH="$2"
            shift 2
            ;;
        --max-depth)
            MAX_DEPTH="$2"
            shift 2
            ;;
        --min-maf)
            MIN_MAF="$2"
            shift 2
            ;;
        --max-edge-distance)
            MAX_EDGE_DISTANCE="$2"
            shift 2
            ;;
        --min-shared-reads)
            MIN_SHARED_READS="$2"
            shift 2
            ;;
        --force-remap)
            FORCE_REMAP="--force-remap"
            shift
            ;;
        --help)
            show_help
            ;;
        *)
            echo "Unknown option: $1"
            show_help
            ;;
    esac
done

# Check required arguments
if [ -z "$HIFI_READS" ]; then
    echo "Error: --hifi-reads is required"
    show_help
fi

if [ -z "$OUTPUT_DIR" ]; then
    echo "Error: --output-dir is required"
    show_help
fi

# Prepare reference option
if [ -n "$REFERENCE" ]; then
    REFERENCE_OPT="--reference $REFERENCE"
else
    REFERENCE_OPT=""
fi

# Prepare force remap option
if [ -n "$FORCE_REMAP" ]; then
    FORCE_REMAP_OPT="--force-remap"
else
    FORCE_REMAP_OPT=""
fi

# Build command
CMD="python3 GenoPhase.py \
    --hifi-reads $HIFI_READS \
    $REFERENCE_OPT \
    --output-dir $OUTPUT_DIR \
    --mode $MODE \
    --threads $THREADS \
    --min-base-quality $MIN_BASE_QUALITY \
    --min-depth $MIN_DEPTH \
    --max-depth $MAX_DEPTH \
    --min-minor-allele-freq $MIN_MAF \
    --max-edge-distance $MAX_EDGE_DISTANCE \
    --min-shared-reads $MIN_SHARED_READS \
    $FORCE_REMAP_OPT"

# Execute command
echo "Running GenoPhase with the following command:"
echo "$CMD"
eval "$CMD"
