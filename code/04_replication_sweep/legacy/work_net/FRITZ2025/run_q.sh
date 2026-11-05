#!/bin/bash
# usage: run_q.sh m0 hist   (the final runs: Stepping MCMC-MLE, 2 PSOCK workers)
# Ported: paths from SWEEP_CODE / SWEEP_DATA / SWEEP_OUT; flock on the old sweep R.lock removed (the pipeline caller
# serialises R jobs); logs and outputs go under $SWEEP_OUT/work_net. Not called by the pipeline; legacy/_slow/<id>.R
# reproduces these calls.
: "${SWEEP_CODE:?}" "${SWEEP_DATA:?}" "${SWEEP_OUT:?}"
C=$SWEEP_CODE/work_net; B=$SWEEP_OUT/work_net; mkdir -p $B/FRITZ2025
for m in "$@"; do
  Rscript $C/FRITZ2025/fritz_run.R $m Stepping 2 > $B/FRITZ2025/log_${m}_Stepping.txt 2>&1
done
echo QDONE >> $B/FRITZ2025/qdone.txt
