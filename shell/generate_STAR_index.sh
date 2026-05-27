#!/usr/bin/env bash

set -euo pipefail

# ==============================================================================
# Default settings
# ==============================================================================

USE_MODULES="false"
STAR_MODULE="STAR/2.7.11b-GCC-13.2.0"

RUN_THREAD_N=""
RUN_MODE="genomeGenerate"
GENOME_DIR=""
GENOME_FASTA_FILES=""
SJDB_GTF_FILE=""
LIMIT_GENOME_GENERATE_RAM=""
GENOME_SA_INDEX_N_BASES=""
SJDB_OVERHANG=""
GENOME_CHR_BIN_N_BITS=""

# ==============================================================================
# Argument parser
# ==============================================================================

while [[ $# -gt 0 ]]; do
    case "$1" in
        --runThreadN)
            RUN_THREAD_N="$2"
            shift 2
            ;;
        --runMode)
            RUN_MODE="$2"
            shift 2
            ;;
        --genomeDir)
            GENOME_DIR="$2"
            shift 2
            ;;
        --genomeFastaFiles)
            GENOME_FASTA_FILES="$2"
            shift 2
            ;;
        --sjdbGTFfile)
            SJDB_GTF_FILE="$2"
            shift 2
            ;;
        --limitGenomeGenerateRAM)
            LIMIT_GENOME_GENERATE_RAM="$2"
            shift 2
            ;;
        --genomeSAindexNbases)
            GENOME_SA_INDEX_N_BASES="$2"
            shift 2
            ;;
        --sjdbOverhang)
            SJDB_OVERHANG="$2"
            shift 2
            ;;
        --genomeChrBinNbits)
            GENOME_CHR_BIN_N_BITS="$2"
            shift 2
            ;;
        --use-modules)
            USE_MODULES="$2"
            shift 2
            ;;
        --star-module)
            STAR_MODULE="$2"
            shift 2
            ;;
        -h|--help)
            cat <<EOF
Usage:
  generate_STAR_index.sh \\
    --runThreadN <threads> \\
    --runMode genomeGenerate \\
    --genomeDir <star_index_output_dir> \\
    --genomeFastaFiles <reference.fa> \\
    --sjdbGTFfile <reference.gtf> \\
    --limitGenomeGenerateRAM <bytes> \\
    --genomeSAindexNbases <int> \\
    --sjdbOverhang <int> \\
    --genomeChrBinNbits <int> \\
    --use-modules true|false \\
    --star-module <module_name>
EOF
            exit 0
            ;;
        *)
            echo "ERROR: Unknown argument: $1" >&2
            exit 1
            ;;
    esac
done

# ==============================================================================
# Checks
# ==============================================================================

if [[ -z "$RUN_THREAD_N" ]]; then
    echo "ERROR: Missing --runThreadN" >&2
    exit 1
fi

if [[ -z "$RUN_MODE" ]]; then
    echo "ERROR: Missing --runMode" >&2
    exit 1
fi

if [[ -z "$GENOME_DIR" ]]; then
    echo "ERROR: Missing --genomeDir" >&2
    exit 1
fi

if [[ -z "$GENOME_FASTA_FILES" ]]; then
    echo "ERROR: Missing --genomeFastaFiles" >&2
    exit 1
fi

if [[ -z "$SJDB_GTF_FILE" ]]; then
    echo "ERROR: Missing --sjdbGTFfile" >&2
    exit 1
fi

if [[ -z "$LIMIT_GENOME_GENERATE_RAM" ]]; then
    echo "ERROR: Missing --limitGenomeGenerateRAM" >&2
    exit 1
fi

if [[ -z "$GENOME_SA_INDEX_N_BASES" ]]; then
    echo "ERROR: Missing --genomeSAindexNbases" >&2
    exit 1
fi

if [[ -z "$SJDB_OVERHANG" ]]; then
    echo "ERROR: Missing --sjdbOverhang" >&2
    exit 1
fi

if [[ -z "$GENOME_CHR_BIN_N_BITS" ]]; then
    echo "ERROR: Missing --genomeChrBinNbits" >&2
    exit 1
fi

if [[ ! -f "$GENOME_FASTA_FILES" ]]; then
    echo "ERROR: FASTA file not found:" >&2
    echo "$GENOME_FASTA_FILES" >&2
    exit 1
fi

if [[ ! -f "$SJDB_GTF_FILE" ]]; then
    echo "ERROR: GTF file not found:" >&2
    echo "$SJDB_GTF_FILE" >&2
    exit 1
fi

# ==============================================================================
# Load STAR
# ==============================================================================

if [[ "$USE_MODULES" == "true" ]]; then
    echo "Loading STAR module: $STAR_MODULE"

    module purge
    module load "$STAR_MODULE"
else
    echo "Not using environment modules. Expecting STAR to be available in PATH."
fi

if ! command -v STAR >/dev/null 2>&1; then
    echo "ERROR: STAR executable not found in PATH." >&2
    exit 1
fi

# ==============================================================================
# Generate STAR index
# ==============================================================================

now="$(date +"%T")"
echo "$now   Starting STAR genome generation"
echo "Reference FASTA:        $GENOME_FASTA_FILES"
echo "Reference GTF:          $SJDB_GTF_FILE"
echo "STAR index dir:         $GENOME_DIR"
echo "Threads:                $RUN_THREAD_N"
echo "Run mode:               $RUN_MODE"
echo "SA index bases:         $GENOME_SA_INDEX_N_BASES"
echo "sjdbOverhang:           $SJDB_OVERHANG"
echo "Max RAM:                $LIMIT_GENOME_GENERATE_RAM"
echo "Use modules:            $USE_MODULES"
echo "STAR module:            $STAR_MODULE"

mkdir -p "$GENOME_DIR"

if [[ -z "$GENOME_DIR" || "$GENOME_DIR" == "/" ]]; then
    echo "ERROR: Refusing to clean unsafe GENOME_DIR: '$GENOME_DIR'" >&2
    exit 1
fi

if [[ -n "$(find "$GENOME_DIR" -mindepth 1 -maxdepth 1 -print -quit)" ]]; then
    echo "WARNING: STAR index folder is not empty. Previous results will be overwritten:" >&2
    echo "$GENOME_DIR" >&2

    find "$GENOME_DIR" -mindepth 1 -maxdepth 1 -exec rm -rf {} +
fi

cd "$GENOME_DIR"

STAR \
    --runThreadN "$RUN_THREAD_N" \
    --runMode "$RUN_MODE" \
    --genomeDir "$GENOME_DIR" \
    --genomeFastaFiles "$GENOME_FASTA_FILES" \
    --sjdbGTFfile "$SJDB_GTF_FILE" \
    --limitGenomeGenerateRAM "$LIMIT_GENOME_GENERATE_RAM" \
    --genomeSAindexNbases "$GENOME_SA_INDEX_N_BASES" \
    --sjdbOverhang "$SJDB_OVERHANG" \
    --genomeChrBinNbits "$GENOME_CHR_BIN_N_BITS"

now="$(date +"%T")"
echo "$now   Finished STAR genome generation"