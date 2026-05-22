#!/usr/bin/env bash
set -euo pipefail

################################################################################
# QC filter FASTQ files using seqtk
#
# Works in:
#   - local mode
#   - SLURM non-array mode
#   - SLURM array mode
#
# Expected manifest columns:
#   pipeline_name
#   symlink_file
#
# The script appends nothing to the manifest itself.
# R should later add manifest$qc_filtered_paths based on the output folder.
################################################################################

timestamp() {
  date +"%T"
}

log() {
  echo "$(timestamp)   $*"
}

die() {
  echo "$(timestamp)   ERROR: $*" >&2
  exit 1
}

usage() {
  cat <<'EOF'
Usage:
  qc_filter_reads.sh \
    --manifest <standardized_fastq_manifest.tsv> \
    --output-dir <QC_filtered_folder> \
    --min-qual <int> \
    --qual-offset <int> \
    --min-length <int> \
    [--use-modules true|false] \
    [--seqtk-module <module_name>]

Required:
  --manifest       TSV manifest from prepare_fastq_inputs()
  --output-dir     Folder for QC-filtered FASTQ files
  --min-qual       Minimum base quality for seqtk -q
  --qual-offset    Quality offset for seqtk -Q
  --min-length     Minimum read length for seqtk -L

Optional:
  --use-modules    Whether to load seqtk via environment modules. Default: false
  --seqtk-module   Module name for seqtk. Default: seqtk/1.3-GCC-11.2.0

Example:
  qc_filter_reads.sh \
    --manifest standardized_fastq_manifest.tsv \
    --output-dir QC_filtered \
    --min-qual 20 \
    --qual-offset 20 \
    --min-length 71 \
    --use-modules true \
    --seqtk-module seqtk/1.3-GCC-11.2.0
EOF
}

################################################################################
# Defaults
################################################################################

MANIFEST=""
OUTPUT_DIR=""
MIN_QUAL=""
QUAL_OFFSET=""
MIN_LENGTH=""
USE_MODULES="false"
SEQTK_MODULE="seqtk/1.3-GCC-11.2.0"

################################################################################
# Parse arguments
################################################################################

while [[ $# -gt 0 ]]; do
  case "$1" in
    --manifest)
      MANIFEST="${2:-}"
      shift 2
      ;;
    --output-dir)
      OUTPUT_DIR="${2:-}"
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
    --use-modules)
      USE_MODULES="${2:-}"
      shift 2
      ;;
    --seqtk-module)
      SEQTK_MODULE="${2:-}"
      shift 2
      ;;
    --help|-h)
      usage
      exit 0
      ;;
    *)
      die "Unknown argument: $1"
      ;;
  esac
done

################################################################################
# Validate arguments
################################################################################

[[ -n "$MANIFEST" ]] || die "--manifest is required"
[[ -n "$OUTPUT_DIR" ]] || die "--output-dir is required"
[[ -n "$MIN_QUAL" ]] || die "--min-qual is required"
[[ -n "$QUAL_OFFSET" ]] || die "--qual-offset is required"
[[ -n "$MIN_LENGTH" ]] || die "--min-length is required"

[[ -f "$MANIFEST" ]] || die "Manifest does not exist: $MANIFEST"

[[ "$MIN_QUAL" =~ ^[0-9]+$ ]] || die "--min-qual must be an integer"
[[ "$QUAL_OFFSET" =~ ^[0-9]+$ ]] || die "--qual-offset must be an integer"
[[ "$MIN_LENGTH" =~ ^[0-9]+$ ]] || die "--min-length must be an integer"

if [[ "$USE_MODULES" != "true" && "$USE_MODULES" != "false" ]]; then
  die "--use-modules must be true or false"
fi

mkdir -p "$OUTPUT_DIR"

################################################################################
# Load seqtk
################################################################################

if [[ "$USE_MODULES" == "true" ]]; then
  log "Loading seqtk module: $SEQTK_MODULE"

  # Some local machines do not have the module command.
  if ! command -v module >/dev/null 2>&1; then
    die "USE_MODULES=true but 'module' command is not available"
  fi

  module purge
  module load "$SEQTK_MODULE"
fi

command -v seqtk >/dev/null 2>&1 || die "seqtk not found on PATH"

################################################################################
# Determine array/local task settings
################################################################################

if [[ -n "${SLURM_ARRAY_TASK_ID:-}" ]]; then
  ARRAY_TASK_ID="$SLURM_ARRAY_TASK_ID"
  ARRAY_TASK_MIN="$SLURM_ARRAY_TASK_MIN"
  ARRAY_TASK_MAX="$SLURM_ARRAY_TASK_MAX"
  ARRAY_STEP="${SLURM_ARRAY_TASK_STEP:-1}"
else
  # Local mode or non-array SLURM mode:
  # one pseudo-task processes all samples.
  ARRAY_TASK_ID=0
  ARRAY_TASK_MIN=0
  ARRAY_TASK_MAX=0
  ARRAY_STEP=1
fi

if [[ "$ARRAY_STEP" -ne 1 ]]; then
  die "This script expects a contiguous Slurm array with step size 1, e.g. --array=0-2."
fi

zero_based_task_index=$(( ARRAY_TASK_ID - ARRAY_TASK_MIN ))
total_array_tasks=$(( ARRAY_TASK_MAX - ARRAY_TASK_MIN + 1 ))

if (( zero_based_task_index < 0 || zero_based_task_index >= total_array_tasks )); then
  die "Invalid array task index calculation."
fi

log "Array/local task ID: $ARRAY_TASK_ID"
log "Array task min: $ARRAY_TASK_MIN"
log "Array task max: $ARRAY_TASK_MAX"
log "Array step: $ARRAY_STEP"
log "Zero-based task index: $zero_based_task_index"
log "Total array/local tasks: $total_array_tasks"

################################################################################
# Read manifest
################################################################################

declare -a sample_names=()
declare -a input_fastqs=()

header="$(head -n 1 "$MANIFEST")"
IFS=$'\t' read -r -a header_cols <<< "$header"

pipeline_col=-1
symlink_file_col=-1

for idx in "${!header_cols[@]}"; do
  case "${header_cols[$idx]}" in
    pipeline_name)
      pipeline_col="$idx"
      ;;
    symlink_file)
      symlink_file_col="$idx"
      ;;
  esac
done

if (( pipeline_col < 0 )); then
  die "Manifest is missing required column: pipeline_name"
fi

if (( symlink_file_col < 0 )); then
  die "Manifest is missing required column: symlink_file"
fi

while IFS=$'\t' read -r -a fields; do
  [[ "${fields[$pipeline_col]:-}" == "pipeline_name" ]] && continue
  [[ -z "${fields[$pipeline_col]:-}" ]] && continue

  sample_names+=("${fields[$pipeline_col]}")
  input_fastqs+=("${fields[$symlink_file_col]}")
done < <(tail -n +2 "$MANIFEST")

total_samples="${#sample_names[@]}"

if [[ "$total_samples" -eq 0 ]]; then
  die "Manifest contains no samples: $MANIFEST"
fi

log "Total samples: $total_samples"

################################################################################
# Main execution: stride through FASTQ files
################################################################################

processed_any=false

for (( sample_index=zero_based_task_index; sample_index<total_samples; sample_index+=total_array_tasks )); do
  sample_name="${sample_names[$sample_index]}"
  input_fastq="${input_fastqs[$sample_index]}"

  if [[ ! -f "$input_fastq" ]]; then
    die "Input FASTQ does not exist for $sample_name: $input_fastq"
  fi

  filename="$(basename "$input_fastq")"
  output_fastq="$OUTPUT_DIR/$filename"

  log "Task $ARRAY_TASK_ID processing manifest index $sample_index: $sample_name"
  log "  input:  $input_fastq"
  log "  output: $output_fastq"

  seqtk seq \
    -q"$MIN_QUAL" \
    -Q"$QUAL_OFFSET" \
    -L"$MIN_LENGTH" \
    -n N \
    "$input_fastq" | gzip > "$output_fastq"

  processed_any=true
done

if [[ "$processed_any" == "false" ]]; then
  log "No samples assigned to task $ARRAY_TASK_ID. Finishing successfully."
else
  log "Task $ARRAY_TASK_ID completed all assigned samples."
fi