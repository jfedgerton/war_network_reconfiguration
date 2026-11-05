#!/bin/bash
# Ported: paths from SWEEP_CODE / SWEEP_DATA / SWEEP_OUT; flock on the old sweep R.lock removed (the pipeline caller
# serialises R jobs); logs and outputs go under $SWEEP_OUT/work_net. Not called by the pipeline; legacy/_slow/<id>.R
# reproduces these calls.
: "${SWEEP_CODE:?}" "${SWEEP_DATA:?}" "${SWEEP_OUT:?}"
: "${SWEEP_LIB_CHEN:?}"
export CHEN_OLD=1 R_LIBS=$SWEEP_LIB_CHEN OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
C=$SWEEP_CODE/work_net/CHEN2021; W=$SWEEP_OUT/work_net/CHEN2021; mkdir -p $W
# original first line waited for a running PID (while kill -0 1556 ...); dropped.
for m in m0 hist; do
  echo "$(date) start MCMLE $m" >> $W/mcmle_queue.log
  Rscript $C/chen_run.R $m MCMLE 2 2000 4000 8 > $W/log_${m}_MCMLE.txt 2>&1
  echo "$(date) end MCMLE $m rc=$?" >> $W/mcmle_queue.log
done
echo MCMLE_DONE >> $W/mcmle_queue.log
