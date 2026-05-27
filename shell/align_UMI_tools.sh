#!/usr/bin/env bash
set -euo pipefail

################################################################################
# align_UMI_tools.sh
#
# Modes:
#   DATA_TYPE=reads
#     QC-filtered FASTQ -> umi_tools extract/cut -> STAR mapping -> BAM index -> idxstats
#
#   DATA_TYPE=umis
#     QC-filtered FASTQ -> umi_tools extract -> STAR mapping -> BAM index
#     -> umi_tools dedup -> dedup BAM index -> idxstats
#
# Expected manifest columns:
#   pipeline_name
#   qc_filtered_paths
#
# Array behavior:
#   - local / non-array SLURM: one pseudo-task processes all samples
#   - SLURM array: each task processes every N-th sample
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
  align_UMI_tools.sh \
    --manifest <manifest.tsv> \
    --output-dir <output_folder> \
    --star-index-folder <STAR_index_folder> \
    --data-type <reads|umis> \
    --umi-regex <regex> \
    --threads <int> \
    --use-modules true|false \
    --star-module <module_name> \
    --samtools-module <module_name> \
    --umi-tools-module <module_name>

Required:
  --manifest             Manifest with pipeline_name and qc_filtered_paths
  --output-dir           Main output folder
  --star-index-folder    Existing STAR genome index folder
  --data-type            reads or umis
  --threads              Threads for STAR
  --umi-regex            Cutting the Read down to just the sgRNA (and finding UMIs)

Module options:
  --use-modules          true or false
  --star-module          STAR module name
  --samtools-module      SAMtools module name
  --umi-tools-module     UMI-tools module name

Example:
  align_UMI_tools.sh \
    --manifest standardized_fastq_manifest.tsv \
    --output-dir /path/to/output \
    --star-index-folder /path/to/star_index \
    --data-type umis \
    --umi-regex '^(?P<umi_1>.{10})' \
    --threads 10 \
    --use-modules true \
    --star-module STAR/2.7.11b-GCC-13.2.0 \
    --samtools-module SAMtools/1.21-GCC-13.3.0 \
    --umi-tools-module UMI-tools/1.1.2-foss-2021b-Python-3.9.6
EOF
}

################################################################################
# Defaults
################################################################################

MANIFEST=""
OUTPUT_DIR=""
STAR_INDEX_FOLDER=""
DATA_TYPE=""
UMI_REGEX=""
THREADS="1"
LIMIT_BAM_SORT_RAM=""

USE_MODULES="false"
STAR_MODULE="STAR/2.7.11b-GCC-13.2.0"
SAMTOOLS_MODULE="SAMtools/1.21-GCC-13.3.0"
UMI_TOOLS_MODULE="UMI-tools/1.1.2-foss-2021b-Python-3.9.6"

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
    --star-index-folder)
      STAR_INDEX_FOLDER="${2:-}"
      shift 2
      ;;
    --data-type)
      DATA_TYPE="${2:-}"
      shift 2
      ;;
    --umi-regex)
      UMI_REGEX="${2:-}"
      shift 2
      ;;
    --threads)
      THREADS="${2:-}"
      shift 2
      ;;
    --limit-bam-sort-ram)
      LIMIT_BAM_SORT_RAM="${2:-}"
      shift 2
      ;;
    --use-modules)
      USE_MODULES="${2:-}"
      shift 2
      ;;
    --star-module)
      STAR_MODULE="${2:-}"
      shift 2
      ;;
    --samtools-module)
      SAMTOOLS_MODULE="${2:-}"
      shift 2
      ;;
    --umi-tools-module)
      UMI_TOOLS_MODULE="${2:-}"
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
[[ -n "$STAR_INDEX_FOLDER" ]] || die "--star-index-folder is required"
[[ -n "$DATA_TYPE" ]] || die "--data-type is required"
[[ -n "$THREADS" ]] || die "--threads is required"
[[ -n "$LIMIT_BAM_SORT_RAM" ]] || die "--limit-bam-sort-ram is required"
[[ "$LIMIT_BAM_SORT_RAM" =~ ^[0-9]+$ ]] || die "--limit-bam-sort-ram must be an integer"

[[ -f "$MANIFEST" ]] || die "Manifest does not exist: $MANIFEST"
[[ -d "$STAR_INDEX_FOLDER" ]] || die "STAR index folder does not exist: $STAR_INDEX_FOLDER"

[[ "$DATA_TYPE" == "reads" || "$DATA_TYPE" == "umis" ]] || die "--data-type must be reads or umis"
[[ "$THREADS" =~ ^[0-9]+$ ]] || die "--threads must be an integer"
(( THREADS >= 1 )) || die "--threads must be >= 1"

if [[ -z "$UMI_REGEX" ]]; then
  die "--umi-regex is required because umi_tools extract is always run before STAR"
fi

if [[ "$USE_MODULES" != "true" && "$USE_MODULES" != "false" ]]; then
  die "--use-modules must be true or false"
fi

################################################################################
# Output folders
################################################################################

UMI_EXTRACTED="$OUTPUT_DIR/UMI_extracted"
MAPPED="$OUTPUT_DIR/mapped"
DEDUP="$OUTPUT_DIR/dedup"
LOGS="$OUTPUT_DIR/logs"

mkdir -p "$UMI_EXTRACTED" "$MAPPED" "$DEDUP" "$LOGS"

################################################################################
# Load/check tools
################################################################################

load_module_if_needed() {
  local module_name="$1"

  if [[ "$USE_MODULES" == "true" ]]; then
    if ! command -v module >/dev/null 2>&1; then
      die "USE_MODULES=true but 'module' command is not available"
    fi

    module load "$module_name"
  fi
}

purge_modules_if_needed() {
  if [[ "$USE_MODULES" == "true" ]]; then
    module purge
  fi
}

################################################################################
# Determine array/local task settings
################################################################################

if [[ -n "${SLURM_ARRAY_TASK_ID:-}" ]]; then
  ARRAY_TASK_ID="$SLURM_ARRAY_TASK_ID"
  ARRAY_TASK_MIN="$SLURM_ARRAY_TASK_MIN"
  ARRAY_TASK_MAX="$SLURM_ARRAY_TASK_MAX"
  ARRAY_STEP="${SLURM_ARRAY_TASK_STEP:-1}"
else
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
log "Zero-based task index: $zero_based_task_index"
log "Total array/local tasks: $total_array_tasks"

################################################################################
# Read manifest using column-name lookup
################################################################################

declare -a sample_names=()
declare -a input_fastqs=()

header="$(head -n 1 "$MANIFEST")"
IFS=$'\t' read -r -a header_cols <<< "$header"

pipeline_col=-1
qc_fastq_col=-1

for idx in "${!header_cols[@]}"; do
  case "${header_cols[$idx]}" in
    pipeline_name)
      pipeline_col="$idx"
      ;;
    qc_filtered_paths)
      qc_fastq_col="$idx"
      ;;
  esac
done

if (( pipeline_col < 0 )); then
  die "Manifest is missing required column: pipeline_name"
fi

if (( qc_fastq_col < 0 )); then
  die "Manifest is missing required column: qc_filtered_paths"
fi

while IFS=$'\t' read -r -a fields; do
  [[ -z "${fields[$pipeline_col]:-}" ]] && continue

  sample_names+=("${fields[$pipeline_col]}")
  input_fastqs+=("${fields[$qc_fastq_col]}")
done < <(tail -n +2 "$MANIFEST")

total_samples="${#sample_names[@]}"

if [[ "$total_samples" -eq 0 ]]; then
  die "Manifest contains no samples: $MANIFEST"
fi

log "Total samples: $total_samples"

################################################################################
# Per-sample functions
################################################################################

run_umi_extraction_one_sample() {
  local sample_name="$1"
  local input_fastq="$2"

  local output_fastq="$UMI_EXTRACTED/${sample_name}.fastq.gz"

  log "Starting UMI extraction for $sample_name"
  log "  input:  $input_fastq"
  log "  output: $output_fastq"

  umi_tools extract \
    --stdin "$input_fastq" \
    --stdout "$output_fastq" \
    --extract-method=regex \
    --bc-pattern="$UMI_REGEX"

  log "Finished UMI extraction for $sample_name"
}

run_star_mapping_one_sample() {
  local sample_name="$1"
  local input_fastq="$2"

  log "Starting STAR mapping for $sample_name"
  log "  input: $input_fastq"

  STAR \
    --runThreadN "$THREADS" \
    --genomeDir "$STAR_INDEX_FOLDER" \
    --readFilesCommand zcat \
    --readFilesIn "$input_fastq" \
    --outFileNamePrefix "$MAPPED/${sample_name}_" \
    --limitBAMsortRAM "$LIMIT_BAM_SORT_RAM" \
    --outSAMtype BAM SortedByCoordinate

  log "Finished STAR mapping for $sample_name"
}

run_bam_index_one_sample() {
  local bam_file="$1"

  log "Indexing BAM: $bam_file"
  samtools index "$bam_file"
}

run_umi_dedup_one_sample() {
  local sample_name="$1"
  local mapped_bam="$2"

  local dedup_bam="$DEDUP/${sample_name}_dedup.bam"
  local dedup_log="$DEDUP/${sample_name}.dedup.log"
  local output_stats="$DEDUP/${sample_name}"

  log "Starting UMI deduplication for $sample_name"
  log "  input:  $mapped_bam"
  log "  output: $dedup_bam"

  umi_tools dedup \
    --stdin="$mapped_bam" \
    --log="$dedup_log" \
    --output-stats="$output_stats" \
    --method adjacency \
    -S "$dedup_bam"

  log "Finished UMI deduplication for $sample_name"
}

run_idxstats_one_sample() {
  local bam_file="$1"
  local output_dir="$2"

  local filename
  filename="$(basename "$bam_file" .bam)"

  log "Running idxstats for $bam_file"

  samtools idxstats "$bam_file" > "$output_dir/${filename}_idxstats.txt"
}

################################################################################
# Main execution
################################################################################

processed_any=false

for (( sample_index=zero_based_task_index; sample_index<total_samples; sample_index+=total_array_tasks )); do
  sample_name="${sample_names[$sample_index]}"
  qc_fastq="${input_fastqs[$sample_index]}"

  if [[ ! -f "$qc_fastq" ]]; then
    die "QC-filtered FASTQ does not exist for $sample_name: $qc_fastq"
  fi

  log "Task $ARRAY_TASK_ID processing manifest index $sample_index: $sample_name"

  # Always run umi_tools extract first.
  # For DATA_TYPE=reads, this is used to trim/cut reads down to the sgRNA sequence.
  # For DATA_TYPE=umis, this also extracts UMIs into the read name.
  purge_modules_if_needed
  load_module_if_needed "$UMI_TOOLS_MODULE"
  command -v umi_tools >/dev/null 2>&1 || die "umi_tools not found on PATH"
  run_umi_extraction_one_sample "$sample_name" "$qc_fastq"

  extracted_fastq="$UMI_EXTRACTED/${sample_name}.fastq.gz"

  if [[ ! -f "$extracted_fastq" ]]; then
    die "UMI-extracted FASTQ does not exist for $sample_name: $extracted_fastq"
  fi

  # Map extracted/trimmed reads to STAR reference.
  purge_modules_if_needed
  load_module_if_needed "$STAR_MODULE"
  command -v STAR >/dev/null 2>&1 || die "STAR not found on PATH"
  run_star_mapping_one_sample "$sample_name" "$extracted_fastq"

  mapped_bam="$MAPPED/${sample_name}_Aligned.sortedByCoord.out.bam"

  if [[ ! -f "$mapped_bam" ]]; then
    die "Mapped BAM does not exist for $sample_name: $mapped_bam"
  fi

  purge_modules_if_needed
  load_module_if_needed "$SAMTOOLS_MODULE"
  command -v samtools >/dev/null 2>&1 || die "samtools not found on PATH"
  run_bam_index_one_sample "$mapped_bam"

  if [[ "$DATA_TYPE" == "reads" ]]; then
    run_idxstats_one_sample "$mapped_bam" "$MAPPED"

  elif [[ "$DATA_TYPE" == "umis" ]]; then
    purge_modules_if_needed
    load_module_if_needed "$UMI_TOOLS_MODULE"
    command -v umi_tools >/dev/null 2>&1 || die "umi_tools not found on PATH"
    run_umi_dedup_one_sample "$sample_name" "$mapped_bam"

    dedup_bam="$DEDUP/${sample_name}_dedup.bam"

    if [[ ! -f "$dedup_bam" ]]; then
      die "Deduplicated BAM does not exist for $sample_name: $dedup_bam"
    fi

    purge_modules_if_needed
    load_module_if_needed "$SAMTOOLS_MODULE"
    command -v samtools >/dev/null 2>&1 || die "samtools not found on PATH"
    run_bam_index_one_sample "$dedup_bam"
    run_idxstats_one_sample "$dedup_bam" "$DEDUP"
  fi

  processed_any=true
done

if [[ "$processed_any" == "false" ]]; then
  log "No samples assigned to task $ARRAY_TASK_ID. Finishing successfully."
else
  log "Task $ARRAY_TASK_ID completed all assigned samples."
fi

################################################################################
# Clean TSV files from UMI-tools output
################################################################################

if [[ "$DATA_TYPE" == "umis" ]]; then
  log "Removing UMI-tools TSV stats files from dedup folder"

  find "$DEDUP" -maxdepth 1 -type f -name "*.tsv" -delete

  log "Finished removing TSV files"
fi

log "align_UMI_tools processing complete."