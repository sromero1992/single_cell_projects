#!/usr/bin/env bash
# =============================================================================
# scenic.sh  --  one-word driver for the pySCENIC run (wraps run_pyscenic2.sh)
# =============================================================================
#   ./scenic.sh status    # where am I? (running? + stage markers + done labels)
#   ./scenic.sh watch     # live-tail the log
#   ./scenic.sh run       # stop stragglers, launch detached, then tail
#   ./scenic.sh stop      # kill the running pyscenic job
#   ./scenic.sh clean     # delete outputs but KEEP input looms (dry-run + y/N)
#   ./scenic.sh fresh     # stop + clean + run (full rerun from scratch)
#
# Pass extra args to 12b after the subcommand, e.g.:
#   ./scenic.sh run --labels Stem_cells__Female Stem_cells__Male
#
# Portable: SCENIC_DIR / env name come from the SAME NR4A1_* env vars as config.R.
# Do NOT run 'fresh' or 'clean' while a job is in progress — it would wipe
# in-progress work. Use 'status' / 'watch' until it finishes.
# =============================================================================
set -uo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SCRIPT="${NR4A1_PYSCENIC_RUNNER:-$HERE/run_pyscenic2.sh}"
LOG="${NR4A1_SCENIC_LOG:-$HERE/log_run.txt}"
PYNAME="12b_run_pyscenic.py"

# Resolve SCENIC_DIR exactly like 12b does (for clean/status).
ROOT="${NR4A1_ROOT:-/home/ssromerogon/local_drive/optimus_drive/selim_working_dir/2026_nr4a1_ack/r_process}"
OUT="${NR4A1_OUTPUT:-$ROOT/seurat_output}"
SCENIC_DIR="${NR4A1_SCENIC_DIR:-$OUT/SCENIC}"

running_pids() { pgrep -f "$PYNAME" 2>/dev/null || true; }

cmd_status() {
  echo "SCENIC_DIR : $SCENIC_DIR"
  echo "runner     : $SCRIPT"
  local pids; pids="$(running_pids)"
  if [ -n "$pids" ]; then echo "RUNNING    : yes (pid $pids)"; else echo "RUNNING    : no"; fi
  if [ -f "$LOG" ]; then
    echo "--- stage markers (last 15 lines of interest) ---"
    grep -E '^#{3,}|=== pySCENIC done|\[grn\]|\[ctx\]|\[aucell\]|Nr4a1 regulon|\[FAIL\]|\[SKIP\]' "$LOG" 2>/dev/null | tail -15
  else
    echo "(no log yet at $LOG)"
  fi
  echo "--- finished labels (non-empty aucell.csv) ---"
  if [ -d "$SCENIC_DIR" ]; then
    local any=no
    for d in "$SCENIC_DIR"/*/; do
      [ -e "$d" ] || continue
      if [ -s "${d}aucell.csv" ]; then
        echo "  [done] $(basename "$d")  ($(wc -l < "${d}aucell.csv") lines)"; any=yes
      fi
    done
    [ "$any" = no ] && echo "  (none yet)"
  fi
}

cmd_watch() { if [ -f "$LOG" ]; then tail -f "$LOG"; else echo "no log at $LOG"; exit 1; fi; }

cmd_stop() {
  local pids; pids="$(running_pids)"
  if [ -n "$pids" ]; then
    echo "stopping pyscenic: $pids"
    # shellcheck disable=SC2086
    kill $pids 2>/dev/null || true; sleep 2
    # shellcheck disable=SC2046
    kill -9 $(running_pids) 2>/dev/null || true
    echo "stopped."
  else
    echo "nothing running."
  fi
}

cmd_run() {
  cmd_stop
  echo "launching detached -> $LOG"
  nohup bash "$SCRIPT" "$@" > "$LOG" 2>&1 &
  sleep 1
  echo "pid $(running_pids)"
  echo "(tailing; Ctrl-C stops watching, the job keeps running)"
  tail -f "$LOG"
}

cmd_clean() {
  if [ ! -d "$SCENIC_DIR" ]; then echo "nothing to clean ($SCENIC_DIR absent)."; return; fi
  echo "DRY RUN — would delete these (input *.loom / *_expr_cells_x_genes.csv.gz are KEPT):"
  find "$SCENIC_DIR" -type f ! -name '*.loom' ! -name '*_expr_cells_x_genes.csv.gz' -print | sed 's/^/  /'
  read -r -p "Delete the files above? [y/N] " ans
  case "${ans:-N}" in
    y|Y)
      find "$SCENIC_DIR" -type f ! -name '*.loom' ! -name '*_expr_cells_x_genes.csv.gz' -delete
      echo "cleaned (input looms kept)." ;;
    *) echo "aborted." ;;
  esac
}

cmd_fresh() { cmd_stop; cmd_clean; cmd_run "$@"; }

sub="${1:-}"; [ $# -gt 0 ] && shift || true
case "$sub" in
  status) cmd_status "$@" ;;
  watch)  cmd_watch  "$@" ;;
  run)    cmd_run    "$@" ;;
  stop)   cmd_stop   "$@" ;;
  clean)  cmd_clean  "$@" ;;
  fresh)  cmd_fresh  "$@" ;;
  *) echo "usage: ./scenic.sh {status|watch|run|stop|clean|fresh} [extra args passed to 12b]"; exit 1 ;;
esac
