#!/usr/bin/env bash
set -euo pipefail

timestamp() { date +"%T"; }
log() { echo "$(timestamp)   $*"; }
die() { echo "$(timestamp)   ERROR: $*" >&2; exit 1; }

MODE=""
INPUT_FASTQ=""
OUTPUT_FASTQ=""
INPUT_FASTQ_R1=""
INPUT_FASTQ_R2=""
OUTPUT_FASTQ_R1=""
OUTPUT_FASTQ_R2=""
MIN_QUAL=""
MIN_LENGTH=""
MIN_LENGTH_R1=""
MIN_LENGTH_R2=""
JSON_REPORT=""
READ_COUNT_OUTPUT=""
THREADS="1"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --mode)
      MODE="${2:-}"
      shift 2
      ;;
    --input-fastq)
      INPUT_FASTQ="${2:-}"
      shift 2
      ;;
    --output-fastq)
      OUTPUT_FASTQ="${2:-}"
      shift 2
      ;;
    --input-fastq-r1)
      INPUT_FASTQ_R1="${2:-}"
      shift 2
      ;;
    --input-fastq-r2)
      INPUT_FASTQ_R2="${2:-}"
      shift 2
      ;;
    --output-fastq-r1)
      OUTPUT_FASTQ_R1="${2:-}"
      shift 2
      ;;
    --output-fastq-r2)
      OUTPUT_FASTQ_R2="${2:-}"
      shift 2
      ;;
    --min-qual)
      MIN_QUAL="${2:-}"
      shift 2
      ;;
    --min-length)
      MIN_LENGTH="${2:-}"
      shift 2
      ;;
    --min-length-r1)
      MIN_LENGTH_R1="${2:-}"
      shift 2
      ;;
    --min-length-r2)
      MIN_LENGTH_R2="${2:-}"
      shift 2
      ;;
    --json-report)
      JSON_REPORT="${2:-}"
      shift 2
      ;;
    --read-count-output)
      READ_COUNT_OUTPUT="${2:-}"
      shift 2
      ;;
    --threads)
      THREADS="${2:-}"
      shift 2
      ;;
    *)
      die "Unknown argument: $1"
      ;;
  esac
done

[[ "$MODE" == "SE" || "$MODE" == "PE" ]] || die "--mode must be SE or PE"
[[ -n "$MIN_QUAL" ]] || die "--min-qual is required"
[[ -n "$JSON_REPORT" ]] || die "--json-report is required"
[[ -n "$READ_COUNT_OUTPUT" ]] || die "--read-count-output is required"

[[ "$MIN_QUAL" =~ ^[0-9]+$ ]] || die "--min-qual must be an integer"
[[ "$THREADS" =~ ^[1-9][0-9]*$ ]] || die "--threads must be a positive integer"

command -v cutadapt >/dev/null 2>&1 || die "cutadapt not found on PATH"
command -v python >/dev/null 2>&1 || die "python not found on PATH"

mkdir -p "$(dirname "$JSON_REPORT")"
mkdir -p "$(dirname "$READ_COUNT_OUTPUT")"

if [[ "$MODE" == "SE" ]]; then
  [[ -n "$INPUT_FASTQ" ]] || die "--input-fastq is required in SE mode"
  [[ -n "$OUTPUT_FASTQ" ]] || die "--output-fastq is required in SE mode"
  [[ -n "$MIN_LENGTH" ]] || die "--min-length is required in SE mode"

  [[ -f "$INPUT_FASTQ" ]] || die "Input FASTQ does not exist: $INPUT_FASTQ"
  [[ "$MIN_LENGTH" =~ ^[1-9][0-9]*$ ]] || die "--min-length must be a positive integer"

  mkdir -p "$(dirname "$OUTPUT_FASTQ")"

  log "Mode: SE"
  log "Input: $INPUT_FASTQ"
  log "Output: $OUTPUT_FASTQ"
  log "min_qual: $MIN_QUAL"
  log "min_length: $MIN_LENGTH"
  log "threads: $THREADS"

  cutadapt \
    --cores "$THREADS" \
    --quality-cutoff "$MIN_QUAL" \
    --minimum-length "$MIN_LENGTH" \
    --output "$OUTPUT_FASTQ" \
    --json "$JSON_REPORT" \
    "$INPUT_FASTQ"

else
  [[ -n "$INPUT_FASTQ_R1" ]] || die "--input-fastq-r1 is required in PE mode"
  [[ -n "$INPUT_FASTQ_R2" ]] || die "--input-fastq-r2 is required in PE mode"
  [[ -n "$OUTPUT_FASTQ_R1" ]] || die "--output-fastq-r1 is required in PE mode"
  [[ -n "$OUTPUT_FASTQ_R2" ]] || die "--output-fastq-r2 is required in PE mode"
  [[ -n "$MIN_LENGTH_R1" ]] || die "--min-length-r1 is required in PE mode"
  [[ -n "$MIN_LENGTH_R2" ]] || die "--min-length-r2 is required in PE mode"

  [[ -f "$INPUT_FASTQ_R1" ]] || die "R1 input FASTQ does not exist: $INPUT_FASTQ_R1"
  [[ -f "$INPUT_FASTQ_R2" ]] || die "R2 input FASTQ does not exist: $INPUT_FASTQ_R2"
  [[ "$MIN_LENGTH_R1" =~ ^[1-9][0-9]*$ ]] || die "--min-length-r1 must be a positive integer"
  [[ "$MIN_LENGTH_R2" =~ ^[1-9][0-9]*$ ]] || die "--min-length-r2 must be a positive integer"

  mkdir -p "$(dirname "$OUTPUT_FASTQ_R1")"
  mkdir -p "$(dirname "$OUTPUT_FASTQ_R2")"

  log "Mode: PE"
  log "R1 input: $INPUT_FASTQ_R1"
  log "R2 input: $INPUT_FASTQ_R2"
  log "R1 output: $OUTPUT_FASTQ_R1"
  log "R2 output: $OUTPUT_FASTQ_R2"
  log "min_qual: $MIN_QUAL"
  log "R1 min_length: $MIN_LENGTH_R1"
  log "R2 min_length: $MIN_LENGTH_R2"
  log "threads: $THREADS"

  cutadapt \
    --cores "$THREADS" \
    --quality-cutoff "$MIN_QUAL" \
    --minimum-length "${MIN_LENGTH_R1}:${MIN_LENGTH_R2}" \
    --pair-filter any \
    --output "$OUTPUT_FASTQ_R1" \
    --paired-output "$OUTPUT_FASTQ_R2" \
    --json "$JSON_REPORT" \
    "$INPUT_FASTQ_R1" \
    "$INPUT_FASTQ_R2"
fi

python - "$JSON_REPORT" "$READ_COUNT_OUTPUT" <<'PY'
import json
import sys
from pathlib import Path

json_path = Path(sys.argv[1])
count_path = Path(sys.argv[2])

with json_path.open() as handle:
    report = json.load(handle)

try:
    output_count = int(report["read_counts"]["output"])
except (KeyError, TypeError, ValueError) as error:
    raise SystemExit(f"Could not extract read_counts.output from {json_path}: {error}")

count_path.write_text(f"{output_count}\n")
PY

[[ -s "$JSON_REPORT" ]] || die "Cutadapt JSON report was not created: $JSON_REPORT"
[[ -s "$READ_COUNT_OUTPUT" ]] || die "Read-count output was not created: $READ_COUNT_OUTPUT"

log "Retained reads/read pairs: $(cat "$READ_COUNT_OUTPUT")"
log "QC filtering complete."