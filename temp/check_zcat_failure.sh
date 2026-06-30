#!/usr/bin/env bash
set -euo pipefail

JOB_ID="${1:-55281059}"
LOG_DIR="${2:-}"

echo "============================================================"
echo "Checking Slurm array job: $JOB_ID"
echo "Date: $(date)"
echo "Host: $(hostname)"
echo "============================================================"
echo

if ! command -v sacct >/dev/null 2>&1; then
  echo "ERROR: sacct not found on PATH."
  exit 1
fi

OUTDIR="slurm_array_check_${JOB_ID}"
mkdir -p "$OUTDIR"

RAW_SACCT="$OUTDIR/sacct_raw.tsv"
TASK_SUMMARY="$OUTDIR/task_summary.tsv"
NODE_SUMMARY="$OUTDIR/node_summary.tsv"
FAILED_TASKS="$OUTDIR/failed_tasks.tsv"
LOG_HITS="$OUTDIR/log_failure_signatures.tsv"

echo "Collecting sacct information..."

sacct \
  -j "$JOB_ID" \
  --parsable2 \
  --noheader \
  --format=JobIDRaw,JobID,State,ExitCode,NodeList,Elapsed,Start,End,MaxRSS,AllocCPUS \
  > "$RAW_SACCT"

echo "Raw sacct output written to:"
echo "  $RAW_SACCT"
echo

awk -F'|' '
BEGIN {
  OFS="\t"
  print "JobIDRaw","JobID","ArrayTaskID","State","ExitCode","NodeList","Elapsed","Start","End","MaxRSS","AllocCPUS"
}
$2 ~ /_[0-9]+$/ {
  task_id=$2
  sub(/^.*_/, "", task_id)
  print $1,$2,task_id,$3,$4,$5,$6,$7,$8,$9,$10
}
' "$RAW_SACCT" > "$TASK_SUMMARY"

echo "Per-task summary written to:"
echo "  $TASK_SUMMARY"
echo

echo "============================================================"
echo "Task state counts"
echo "============================================================"
awk -F'\t' 'NR>1 {count[$4]++} END {for (s in count) print s, count[s]}' "$TASK_SUMMARY" | sort
echo

echo "============================================================"
echo "Exit-code counts"
echo "============================================================"
awk -F'\t' 'NR>1 {count[$5]++} END {for (e in count) print e, count[e]}' "$TASK_SUMMARY" | sort
echo

echo "============================================================"
echo "Node summary"
echo "============================================================"

awk -F'\t' '
BEGIN {
  OFS="\t"
}
NR == 1 { next }
{
  node=$6
  state=$4
  exitcode=$5

  total[node]++

  if (state !~ /^COMPLETED/) {
    failed[node]++
  }

  state_count[node,state]++
  exit_count[node,exitcode]++
}
END {
  print "NodeList","TotalTasks","NonCompletedTasks","States","ExitCodes"

  for (node in total) {
    states=""
    exits=""

    for (key in state_count) {
      split(key, a, SUBSEP)
      if (a[1] == node) {
        states = states a[2] ":" state_count[key] ";"
      }
    }

    for (key in exit_count) {
      split(key, a, SUBSEP)
      if (a[1] == node) {
        exits = exits a[2] ":" exit_count[key] ";"
      }
    }

    print node,total[node],failed[node]+0,states,exits
  }
}
' "$TASK_SUMMARY" | sort -k3,3nr -k2,2nr | tee "$NODE_SUMMARY"

echo
echo "Node summary written to:"
echo "  $NODE_SUMMARY"
echo

echo "============================================================"
echo "Failed / non-completed tasks"
echo "============================================================"

awk -F'\t' '
NR == 1 || $4 !~ /^COMPLETED/ {
  print
}
' "$TASK_SUMMARY" | tee "$FAILED_TASKS"

echo
echo "Failed task table written to:"
echo "  $FAILED_TASKS"
echo

if [[ -n "$LOG_DIR" ]]; then
  echo "============================================================"
  echo "Scanning logs in: $LOG_DIR"
  echo "============================================================"

  if [[ ! -d "$LOG_DIR" ]]; then
    echo "WARNING: LOG_DIR does not exist: $LOG_DIR"
  else
    {
      echo -e "LogFile\tIllegalInstruction\tZcatMention\tZeroInputReads\tNodeMention"

      find "$LOG_DIR" -type f \
        \( -name "*${JOB_ID}*" -o -name "*.out" -o -name "*.err" -o -name "*.log" \) \
        | sort \
        | while read -r log_file; do

            illegal="no"
            zcat_hit="no"
            zero_reads="no"
            node_hit=""

            if grep -qi "Illegal instruction" "$log_file"; then
              illegal="yes"
            fi

            if grep -qi "zcat" "$log_file"; then
              zcat_hit="yes"
            fi

            if grep -q "Number of input reads.*|[[:space:]]*0" "$log_file"; then
              zero_reads="yes"
            fi

            node_hit="$(grep -Eim1 'node|hostname|Host:' "$log_file" || true)"
            node_hit="${node_hit//$'\t'/ }"

            if [[ "$illegal" == "yes" || "$zcat_hit" == "yes" || "$zero_reads" == "yes" ]]; then
              echo -e "${log_file}\t${illegal}\t${zcat_hit}\t${zero_reads}\t${node_hit}"
            fi
          done
    } | tee "$LOG_HITS"

    echo
    echo "Log signature table written to:"
    echo "  $LOG_HITS"
  fi
else
  echo "No log directory provided."
  echo "To also scan logs, rerun like:"
  echo "  bash $0 $JOB_ID /path/to/slurm/logs"
  echo
fi

echo "============================================================"
echo "Suggested next commands"
echo "============================================================"
echo
echo "Show tasks sorted by node:"
echo "  column -t -s \$'\\t' $TASK_SUMMARY | sort -k6,6"
echo
echo "Show suspicious nodes first:"
echo "  column -t -s \$'\\t' $NODE_SUMMARY"
echo
echo "Show failed tasks:"
echo "  column -t -s \$'\\t' $FAILED_TASKS"
echo
echo "Done."