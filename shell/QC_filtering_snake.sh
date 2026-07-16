#!/usr/bin/env bash
set -euo pipefail

timestamp() { date +"%T"; }
log() { echo "$(timestamp)   $*"; }
die() { echo "$(timestamp)   ERROR: $*" >&2; exit 1; }

INPUT_FASTQ=""
OUTPUT_FASTQ=""
MIN_QUAL=""
QUAL_OFFSET=""
MIN_LENGTH=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --input-fastq)
      INPUT_FASTQ="${2:-}"
      shift 2
      ;;
    --output-fastq)
      OUTPUT_FASTQ="${2:-}"
      shift 2
      ;;
    --min-qual)
      MIN_QUAL="${2:-}"
      shift 2
      ;;
    --qual-offset)
      QUAL_OFFSET="${2:-}"
      shift 2
      ;;
    --min-length)
      MIN_LENGTH="${2:-}"
      shift 2
      ;;
    *)
      die "Unknown argument: $1"
      ;;
  esac
done

[[ -n "$INPUT_FASTQ" ]] || die "--input-fastq is required"
[[ -n "$OUTPUT_FASTQ" ]] || die "--output-fastq is required"
[[ -n "$MIN_QUAL" ]] || die "--min-qual is required"
[[ -n "$QUAL_OFFSET" ]] || die "--qual-offset is required"
[[ -n "$MIN_LENGTH" ]] || die "--min-length is required"

[[ -f "$INPUT_FASTQ" ]] || die "Input FASTQ does not exist: $INPUT_FASTQ"

[[ "$MIN_QUAL" =~ ^[0-9]+$ ]] || die "--min-qual must be an integer"
[[ "$QUAL_OFFSET" =~ ^[0-9]+$ ]] || die "--qual-offset must be an integer"
[[ "$MIN_LENGTH" =~ ^[0-9]+$ ]] || die "--min-length must be an integer"

mkdir -p "$(dirname "$OUTPUT_FASTQ")"

command -v seqtk >/dev/null 2>&1 || die "seqtk not found on PATH"

log "Input:  $INPUT_FASTQ"
log "Output: $OUTPUT_FASTQ"
log "min_qual: $MIN_QUAL"
log "qual_offset: $QUAL_OFFSET"
log "min_length: $MIN_LENGTH"

seqtk seq \
  -q"$MIN_QUAL" \
  -Q"$QUAL_OFFSET" \
  -L"$MIN_LENGTH" \
  -n N \
  "$INPUT_FASTQ" | gzip > "$OUTPUT_FASTQ"

log "QC filtering complete."