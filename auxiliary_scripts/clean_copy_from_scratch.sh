#!/usr/bin/env bash
#SBATCH -J clean_copy
#SBATCH -A lsteinme                                     
#SBATCH --mem=40g
#SBATCH -N 1
#SBATCH --cpus-per-task=1
#SBATCH -t 04:00:00
#SBATCH --qos normal
#SBATCH -o /g/steinmetz/link/logs/log_%x_%j.out
#SBATCH -e /g/steinmetz/link/logs/log_%x_%j.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=lukas.link@embl.de

set -euo pipefail

SRC="/scratch/link/Amplicon_barcode_analysis"
DST="/g/steinmetz/link/Amplicon_barcode_analysis"

timestamp() { date +"%T"; }
log() { echo "$(timestamp)   $*"; }

if [[ ! -d "$SRC" ]]; then
  echo "ERROR: source directory does not exist: $SRC" >&2
  exit 1
fi

case "$DST/" in
  "$SRC"/*)
    echo "ERROR: destination must not be inside source directory." >&2
    echo "SRC: $SRC" >&2
    echo "DST: $DST" >&2
    exit 1
    ;;
esac

mkdir -p "$DST"

log "Copying everything from:"
log "  SRC: $SRC"
log "  DST: $DST"

rsync -a --human-readable --stats --info=progress2 -- "$SRC"/ "$DST"/

log "Copy completed successfully."
# sbatch ~/Amplicon_barcode_analysis/Lukas_Pipeline/binned_PCR_amplicon_UMI_analysis/auxiliary_scripts/clean_copy_from_scratch.sh