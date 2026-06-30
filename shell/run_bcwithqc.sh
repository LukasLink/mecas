#!/usr/bin/env bash
set -euo pipefail

timestamp() { date +"%T"; }
log() { echo "$(timestamp)   $*"; }
die() { echo "$(timestamp)   ERROR: $*" >&2; exit 1; }

################################################################################
# Defaults
################################################################################

MANIFEST=""
OUTPUT_DIR=""
BCWITHQC_CONFIG=""
STAR_INDEX=""
BCWITHQC_BIN="bcwithqc"

THREADS="1"
KEEP_INTERMEDIARY="true"
EXISTING_RESULTS_MODE="override"
VERBOSITY="-vv"

USE_MODULES="false"
STAR_MODULE="STAR/2.7.11b-GCC-13.2.0"

################################################################################
# Parse arguments
################################################################################

while [[ $# -gt 0 ]]; do
  case "$1" in
    --bcwithqc-bin)
      BCWITHQC_BIN="${2:-}"
      shift 2
      ;;
    --manifest)
      MANIFEST="${2:-}"
      shift 2
      ;;
    --output-dir)
      OUTPUT_DIR="${2:-}"
      shift 2
      ;;
    --bcwithqc-config)
      BCWITHQC_CONFIG="${2:-}"
      shift 2
      ;;
    --star-index-folder)
      STAR_INDEX="${2:-}"
      shift 2
      ;;
    --bcwithqc-bin)
      BCWITHQC_BIN="${2:-}"
      shift 2
      ;;
    --threads)
      THREADS="${2:-}"
      shift 2
      ;;
    --keep-intermediary)
      KEEP_INTERMEDIARY="${2:-}"
      shift 2
      ;;
    --existing-results-mode)
      EXISTING_RESULTS_MODE="${2:-}"
      shift 2
      ;;
    --verbosity)
      VERBOSITY="${2:-}"
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
    --help|-h)
      cat <<EOF
Usage:
  run_bcwithqc.sh \\
    --manifest <bcwithqc_manifest.tsv> \\
    --output-dir <output_folder> \\
    --bcwithqc-config <config.json> \\
    --bcwithqc-bin <bcwithqc executable> \\
    --threads <int> \\
    --use-modules true|false \\
    --star-module <module_name>
EOF
      exit 0
      ;;
    *)
      die "Unknown argument: $1"
      ;;
  esac
done

################################################################################
# Checks
################################################################################

[[ -n "$MANIFEST" ]] || die "--manifest is required"
[[ -n "$OUTPUT_DIR" ]] || die "--output-dir is required"
[[ -n "$BCWITHQC_CONFIG" ]] || die "--bcwithqc-config is required"
[[ -n "$STAR_INDEX" ]] || die "--star-index-folder is required"
[[ -n "$THREADS" ]] || die "--threads is required"

[[ -f "$MANIFEST" ]] || die "Manifest does not exist: $MANIFEST"
[[ -f "$BCWITHQC_CONFIG" ]] || die "bcwithqc config does not exist: $BCWITHQC_CONFIG"
[[ -d "$STAR_INDEX" ]] || die "STAR index folder does not exist: $STAR_INDEX"
[[ "$THREADS" =~ ^[0-9]+$ ]] || die "--threads must be an integer"
(( THREADS >= 1 )) || die "--threads must be >= 1"

if [[ "$USE_MODULES" != "true" && "$USE_MODULES" != "false" ]]; then
  die "--use-modules must be true or false"
fi

case "$EXISTING_RESULTS_MODE" in
  override|abort|skip) ;;
  *) die "--existing-results-mode must be one of: override, abort, skip" ;;
esac

if ! command -v "$BCWITHQC_BIN" >/dev/null 2>&1 && [[ ! -x "$BCWITHQC_BIN" ]]; then
  die "bcwithqc executable not found or not executable: $BCWITHQC_BIN"
fi

################################################################################
# Output folders
################################################################################

RUNS_DIR="$OUTPUT_DIR/bcwithqc_output"
LOGS_DIR="$OUTPUT_DIR/logs"

mkdir -p "$RUNS_DIR" "$LOGS_DIR"

################################################################################
# Helpers
################################################################################

dir_is_nonempty() {
  local d="$1"
  [[ -d "$d" ]] || return 1
  [[ -n "$(find "$d" -mindepth 1 -maxdepth 1 -print -quit)" ]]
}

handle_existing_run_output_dir() {
  local run_output_dir="$1"
  local run_name="$2"

  if dir_is_nonempty "$run_output_dir"; then
    case "$EXISTING_RESULTS_MODE" in
      override)
        log "[$run_name] Existing output directory found; deleting because mode is override: $run_output_dir"
        rm -rf -- "$run_output_dir"
        ;;
      abort)
        die "[$run_name] Output directory is non-empty: $run_output_dir"
        ;;
      skip)
        log "[$run_name] Existing output directory found; skipping because mode is skip: $run_output_dir"
        return 1
        ;;
    esac
  fi

  return 0
}

load_star_if_needed() {
  if [[ "$USE_MODULES" == "true" ]]; then
    if ! command -v module >/dev/null 2>&1; then
      die "USE_MODULES=true but module command is not available"
    fi

    module purge
    module load "$STAR_MODULE"
  fi

  command -v STAR >/dev/null 2>&1 || die "STAR not found on PATH"
}

clean_python_env_for_bcwithqc() {
  unset PYTHONPATH
  unset PYTHONHOME
  export PYTHONNOUSERSITE=1
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

log "Output folder: $OUTPUT_DIR"
log "bcwithqc config: $BCWITHQC_CONFIG"
log "STAR index: $STAR_INDEX"
log "bcwithqc executable: $BCWITHQC_BIN"
log "Threads: $THREADS"
log "Manifest: $MANIFEST"
log "Array/local task ID: $ARRAY_TASK_ID"
log "Zero-based task index: $zero_based_task_index"
log "Total array/local tasks: $total_array_tasks"

################################################################################
# Read manifest
################################################################################

declare -a sample_names=()
declare -a input_dirs=()

header="$(head -n 1 "$MANIFEST")"
IFS=$'\t' read -r -a header_cols <<< "$header"

sample_col=-1
input_dir_col=-1

for idx in "${!header_cols[@]}"; do
  case "${header_cols[$idx]}" in
    pipeline_name|sample_name)
      sample_col="$idx"
      ;;
    bcwithqc_input_dir|input_dir)
      input_dir_col="$idx"
      ;;
  esac
done

if (( sample_col < 0 )); then
  die "Manifest is missing required column: pipeline_name or sample_name"
fi

if (( input_dir_col < 0 )); then
  die "Manifest is missing required column: bcwithqc_input_dir or input_dir"
fi

while IFS=$'\t' read -r -a fields; do
  [[ -z "${fields[$sample_col]:-}" ]] && continue

  sample_names+=("${fields[$sample_col]}")
  input_dirs+=("${fields[$input_dir_col]}")
done < <(tail -n +2 "$MANIFEST")

total_samples="${#sample_names[@]}"

if [[ "$total_samples" -eq 0 ]]; then
  die "Manifest contains no samples: $MANIFEST"
fi

log "Total prepared samples: $total_samples"
################################################################################
# One bcwithqc run
################################################################################

run_bcwithqc_one_input_dir() {
  local input_dir="$1"
  local run_name="$2"

  local run_output_dir="$RUNS_DIR/$run_name"
  local star_dir="$run_output_dir/STAR_files"

  if ! handle_existing_run_output_dir "$run_output_dir" "$run_name"; then
    return 0
  fi

  mkdir -p "$run_output_dir" "$star_dir"

  log "Starting run: $run_name"
  log "Run input dir: $input_dir"
  log "Run output dir: $run_output_dir"

  ##############################################################################
  # 1. bcwithqc preprocess
  ##############################################################################

  log "[$run_name] Starting bcwithqc preprocess"
  
  log "BCWITHQC_BIN: $BCWITHQC_BIN"
  log "PYTHONPATH before cleanup: ${PYTHONPATH:-<unset>}"
  log "PYTHONHOME before cleanup: ${PYTHONHOME:-<unset>}"
  log "PYTHONNOUSERSITE before cleanup: ${PYTHONNOUSERSITE:-<unset>}"
  
  clean_python_env_for_bcwithqc
  
  log "PYTHONPATH after cleanup: ${PYTHONPATH:-<unset>}"
  log "PYTHONHOME after cleanup: ${PYTHONHOME:-<unset>}"
  log "PYTHONNOUSERSITE after cleanup: ${PYTHONNOUSERSITE:-<unset>}"
  
  preprocess_cmd=(
    "$BCWITHQC_BIN" "preprocess"
    "$input_dir"
    "--config=$BCWITHQC_CONFIG"
    "--output-dir=$run_output_dir"
    "--threads=$THREADS"
    "--output-format-bam"
  )

  if [[ -n "$VERBOSITY" ]]; then
    preprocess_cmd+=("$VERBOSITY")
  fi

  "${preprocess_cmd[@]}" 2>&1 | tee "$LOGS_DIR/${run_name}_01_preprocess.log"

  log "[$run_name] Finished bcwithqc preprocess"

  ##############################################################################
  # Find preprocessed FASTQ
  ##############################################################################

  mapfile -t preprocessed_bams < <(
    find "$run_output_dir" -maxdepth 1 \
      -name "*Aligned.out.bam" \
      -type f | sort
  )
  
  if [[ "${#preprocessed_bams[@]}" -eq 0 ]]; then
    die "[$run_name] No preprocessed BAM found in: $run_output_dir"
  fi
  
  if [[ "${#preprocessed_bams[@]}" -gt 1 ]]; then
    printf '%s\n' "${preprocessed_bams[@]}"
    die "[$run_name] Expected exactly one preprocessed BAM per run."
  fi
  
  local preprocessed_bam="${preprocessed_bams[0]}"
  log "[$run_name] Preprocessed BAM: $preprocessed_bam"

  ##############################################################################
  # 2. STAR alignment
  ##############################################################################
  # This was depracted when STAR was removed as obligatory for bcwithqc. 

  # log "[$run_name] Starting STAR alignment"
  # 
  # load_star_if_needed
  # 
  # local star_prefix="$star_dir/${run_name}_"
  # 
  # star_cmd=(
  #   STAR
  #   --runThreadN "$THREADS"
  #   --genomeDir "$STAR_INDEX"
  #   --readFilesIn "$preprocessed_fastq"
  #   --outFileNamePrefix "$star_prefix"
  #   --outFilterMultimapNmax 1
  #   --outSAMtype BAM Unsorted
  #   --outSAMattributes NH HI AS nM GX GN
  # )
  # 
  # if [[ "$preprocessed_fastq" == *.gz ]]; then
  #   star_cmd+=(--readFilesCommand zcat)
  # fi
  # 
  # "${star_cmd[@]}" 2>&1 | tee "$LOGS_DIR/${run_name}_02_STAR.log"
  # 
  # local aligned_bam="${star_prefix}Aligned.out.bam"
  # 
  # [[ -f "$aligned_bam" ]] || die "[$run_name] Expected STAR BAM not found: $aligned_bam"
  # 
  # log "[$run_name] Finished STAR alignment"
  # log "[$run_name] STAR BAM: $aligned_bam"

  ##############################################################################
  # 3. bcwithqc count
  ##############################################################################

  log "[$run_name] Starting bcwithqc count"
  
  log "BCWITHQC_BIN: $BCWITHQC_BIN"
  log "PYTHONPATH before cleanup: ${PYTHONPATH:-<unset>}"
  log "PYTHONHOME before cleanup: ${PYTHONHOME:-<unset>}"
  log "PYTHONNOUSERSITE before cleanup: ${PYTHONNOUSERSITE:-<unset>}"
  
  clean_python_env_for_bcwithqc
  
  log "PYTHONPATH after cleanup: ${PYTHONPATH:-<unset>}"
  log "PYTHONHOME after cleanup: ${PYTHONHOME:-<unset>}"
  log "PYTHONNOUSERSITE after cleanup: ${PYTHONNOUSERSITE:-<unset>}"
  
  count_cmd=(
    "$BCWITHQC_BIN" "count"
    "$run_output_dir"
    "--STAR-output-dir=$run_output_dir"
    "--config=$BCWITHQC_CONFIG"
    "--output-dir=$run_output_dir"
    "--threads=$THREADS"
  )

  if [[ "$KEEP_INTERMEDIARY" == "true" ]]; then
    count_cmd+=("--keep-intermediary")
  fi

  if [[ -n "$VERBOSITY" ]]; then
    count_cmd+=("$VERBOSITY")
  fi

  "${count_cmd[@]}" 2>&1 | tee "$LOGS_DIR/${run_name}_03_count.log"

  log "[$run_name] Finished bcwithqc count"
  log "[$run_name] Completed run"
}

################################################################################
# Main execution
################################################################################

processed_any=false

for (( sample_index=zero_based_task_index; sample_index<total_samples; sample_index+=total_array_tasks )); do
  sample_name="${sample_names[$sample_index]}"
  input_dir="${input_dirs[$sample_index]}"

  [[ -d "$input_dir" ]] || die "Prepared input directory missing for $sample_name: $input_dir"

  log "Task $ARRAY_TASK_ID processing manifest index $sample_index: $sample_name"
  run_bcwithqc_one_input_dir "$input_dir" "$sample_name"
  processed_any=true
done

if [[ "$processed_any" == "false" ]]; then
  log "No samples assigned to task $ARRAY_TASK_ID. Finishing successfully."
else
  log "Task $ARRAY_TASK_ID completed all assigned samples."
fi

log "bcwithqc processing complete."