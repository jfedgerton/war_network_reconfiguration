#!/bin/bash
# Original q2.sh: waited for FRITZ2025/qdone.txt, then FRITZ2025 MPLE points and CHEN2021 MPLE (current stack; superseded
# by the CHEN_OLD=1 MPLE runs in CHEN2021/run_mcmle.sh, which write_chen.R reads).
# Ported: paths from SWEEP_CODE / SWEEP_DATA / SWEEP_OUT; flock on the old sweep R.lock removed (the pipeline caller
# serialises R jobs); logs and outputs go under $SWEEP_OUT/work_net. Not called by the pipeline; legacy/_slow/<id>.R
# reproduces these calls.
: "${SWEEP_CODE:?}" "${SWEEP_DATA:?}" "${SWEEP_OUT:?}"
C=$SWEEP_CODE/work_net; B=$SWEEP_OUT/work_net
mkdir -p $B/FRITZ2025 $B/CHEN2021
while [ ! -f $B/FRITZ2025/qdone.txt ]; do sleep 20; done
Rscript $C/FRITZ2025/fritz_mple.R > $B/FRITZ2025/log_mple_points.txt 2>&1
for m in m0 hist m26 placebo; do
  Rscript $C/CHEN2021/chen_run.R $m MPLE 1 > $B/CHEN2021/log_${m}_MPLE.txt 2>&1
done
echo Q2DONE > $B/q2done.txt
