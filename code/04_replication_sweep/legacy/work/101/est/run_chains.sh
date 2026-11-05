#!/bin/bash
# Segmented chains for 101 (m0, mDM interleaved): 16 segments x 25,000 iterations = 400,000 per model (burn-in 100,000 discarded at write), thin 250.
# The sweep ran each segment under a shared R lock (flock) so other jobs could interleave; the pipeline runs studies
# one at a time, so the lock is used only when flock exists and RLOCK is set (macOS has no flock by default).
# port: needs SWEEP_CODE / SWEEP_DATA / SWEEP_OUT / SWEEP_HIST (as set by run_one.R); RSCRIPT defaults to Rscript on PATH.
: "${SWEEP_CODE:?}" "${SWEEP_DATA:?}" "${SWEEP_OUT:?}" "${SWEEP_HIST:?}"; RSCRIPT=${RSCRIPT:-Rscript}
LOCK=(); if [ -n "${RLOCK:-}" ] && command -v flock >/dev/null 2>&1; then LOCK=(flock "$RLOCK"); fi
E=$SWEEP_OUT/work/101/est; mkdir -p "$E"
for k in $(seq 1 16); do
  for m in m0 mDM; do
    f=$(printf "$E/seg_%s_%03d.rds" $m $k)
    [ -f "$f" ] && continue
    echo "$(date +%T) start $m $k" >> "$E/runner.log"
    "${LOCK[@]}" "$RSCRIPT" "$SWEEP_CODE/work/101/rerun_est_101.R" seg $m $k 25000 250 >> "$E/runner.log" 2>&1 || { echo "FAILED $m $k" >> "$E/runner.log"; exit 1; }
  done
done
echo "$(date +%T) ALL DONE" >> "$E/runner.log"
# port: final collection step (run by hand in the sweep; see its log_est_101_write_100000.txt): burn-in 100,000 -> results_estaudit/101__smigoabi.csv
"${LOCK[@]}" "$RSCRIPT" "$SWEEP_CODE/work/101/rerun_est_101.R" write 100000 >> "$E/runner.log" 2>&1 || { echo "FAILED write" >> "$E/runner.log"; exit 1; }
