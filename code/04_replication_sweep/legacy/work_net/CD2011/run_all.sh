#!/bin/bash
# Ported: paths from SWEEP_CODE / SWEEP_DATA / SWEEP_OUT; flock on the old sweep R.lock removed (the pipeline caller
# serialises R jobs); logs and outputs go under $SWEEP_OUT/work_net. Not called by the pipeline; legacy/_slow/<id>.R
# reproduces these calls.
: "${SWEEP_CODE:?}" "${SWEEP_DATA:?}" "${SWEEP_OUT:?}"
C=$SWEEP_CODE/work_net; B=$SWEEP_OUT/work_net; mkdir -p $B/CD2011
for m in m0 hist m26 ldv0 ldvhist pbase placebo; do
  Rscript $C/CD2011/est_cd2011.R $m 1000 >> $B/CD2011/log_cd2011.txt 2>&1
done
echo ALLDONE >> $B/CD2011/log_cd2011.txt
