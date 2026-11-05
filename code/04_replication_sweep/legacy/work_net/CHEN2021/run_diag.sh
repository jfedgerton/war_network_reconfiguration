#!/bin/bash
# Ported: paths from SWEEP_CODE / SWEEP_DATA / SWEEP_OUT; flock on the old sweep R.lock removed (the pipeline caller
# serialises R jobs); logs and outputs go under $SWEEP_OUT/work_net. Not called by the pipeline; legacy/_slow/<id>.R
# reproduces these calls.
: "${SWEEP_CODE:?}" "${SWEEP_DATA:?}" "${SWEEP_OUT:?}"
: "${SWEEP_LIB_CHEN:?}"
export CHEN_OLD=1 R_LIBS=$SWEEP_LIB_CHEN OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
C=$SWEEP_CODE/work_net/CHEN2021; W=$SWEEP_OUT/work_net/CHEN2021; mkdir -p $W
while ! grep -q MCMLE_DONE $W/mcmle_queue.log; do sleep 20; done
for m in m0 hist; do Rscript $C/chen_run.R $m DIAG 1 2000 4000 > $W/log_${m}_DIAG.txt 2>&1; done
echo DIAG_DONE >> $W/mcmle_queue.log
