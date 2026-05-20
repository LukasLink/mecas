#!/bin/bash
#SBATCH -J gen_STAR_index    # Job Name                # chmod +x /home/link/Amplicon_barcode_analysis/Lukas_Pipeline/binned_PCR_amplicon_UMI_analysis/auxiliary_scripts/generate_star_index.sh
#SBATCH -A lsteinme             # profile of the group                # sbatch /home/link/Amplicon_barcode_analysis/Lukas_Pipeline/binned_PCR_amplicon_UMI_analysis/auxiliary_scripts/generate_star_index.sh
#SBATCH --mem 128g               # Total memory required for the job #  --dependency=afterok:46151850
#SBATCH -N 1                    # Number of nodes
#SBATCH -n 12                   # Number of CPUs
#SBATCH -t 04:00:00             # Runtime until the job is forcefully canceled
#SBATCH --qos normal 
#SBATCH -o /g/steinmetz/link/logs/log_%x_%A.out
#SBATCH -e /g/steinmetz/link/logs//log_%x_%A.err
#SBATCH --mail-type=BEGIN,END,FAIL        	# notifications for job start, done & fail
#SBATCH --mail-user=lukas.link@embl.de      # send-to address     # notifications for job done & fail

set -euo pipefail

# ==============================================================================
# STAR settings
# ==============================================================================

SJDB_OVERHANG=121              # This is the length read - 1
Genome_SA_index_N_Bases=10     # Calculated by Combined_pipeline_support.rmd
NUM_THREADS=10                 # The number of threads used by STAR; should NOT be higher than cpus-per-task above
MAX_MEM=100000000000           # The number of bytes available to STAR; should NOT be higher than mem set above

# ==============================================================================
# Paths
# ==============================================================================

REF_DIR="/g/steinmetz/link/Amplicon_barcode_analysis/library/ref_CM_specific_library"
STAR_INDEX="/g/steinmetz/link/Amplicon_barcode_analysis/library/star_index_CM_specific_library_Rlenght_122"

REFERENCE_FASTA="${REF_DIR}/reference.fa"
REFERENCE_GTF="${REF_DIR}/reference.gtf"

# ==============================================================================
# Checks
# ==============================================================================

mkdir -p /g/steinmetz/link/logs

if [[ ! -f "$REFERENCE_FASTA" ]]; then
    echo "ERROR: FASTA file not found:"
    echo "$REFERENCE_FASTA"
    exit 1
fi

if [[ ! -f "$REFERENCE_GTF" ]]; then
    echo "ERROR: GTF file not found:"
    echo "$REFERENCE_GTF"
    exit 1
fi

mkdir -p "$STAR_INDEX"

# Avoid accidentally mixing old and new STAR index files
if [[ -n "$(find "$STAR_INDEX" -mindepth 1 -print -quit)" ]]; then
    echo "ERROR: STAR index folder is not empty:"
    echo "$STAR_INDEX"
    echo ""
    echo "Please remove or rename the existing folder before regenerating the index."
    exit 1
fi

# ==============================================================================
# Generate STAR index
# ==============================================================================

now="$(date +"%T")"
echo "$now   Starting STAR genome generation"
echo "Reference FASTA: $REFERENCE_FASTA"
echo "Reference GTF:   $REFERENCE_GTF"
echo "STAR index dir:  $STAR_INDEX"
echo "Threads:         $NUM_THREADS"
echo "SA index bases:  $Genome_SA_index_N_Bases"
echo "sjdbOverhang:    $SJDB_OVERHANG"
echo "Max RAM:         $MAX_MEM"

module purge
ml STAR/2.7.11b-GCC-13.2.0

STAR \
    --runThreadN "$NUM_THREADS" \
    --runMode genomeGenerate \
    --genomeDir "$STAR_INDEX" \
    --genomeFastaFiles "$REFERENCE_FASTA" \
    --sjdbGTFfile "$REFERENCE_GTF" \
    --limitGenomeGenerateRAM "$MAX_MEM" \
    --genomeSAindexNbases "$Genome_SA_index_N_Bases" \
    --sjdbOverhang "$SJDB_OVERHANG"

now="$(date +"%T")"
echo "$now   Finished STAR genome generation"