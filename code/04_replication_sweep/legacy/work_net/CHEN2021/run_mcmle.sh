#!/bin/bash
# CHEN2021 MCMLE on the author's stack (ergm 3.10.1 / multilayer.ergm 0.1.6 / statnet.common 4.3.0 / network 1.15)
# (original: the MCMLE loop here returned within seconds; the MCMLE fits used by write_chen.R came from run_mcmle2.sh.
#  The MPLE loop below produced res_<m>_MPLE_old.rds.)
# Ported: paths from SWEEP_CODE / SWEEP_DATA / SWEEP_OUT; flock on the old sweep R.lock removed (the pipeline caller
# serialises R jobs); logs and outputs go under $SWEEP_OUT/work_net. Not called by the pipeline; legacy/_slow/<id>.R
# reproduces these calls.
: "${SWEEP_CODE:?}" "${SWEEP_DATA:?}" "${SWEEP_OUT:?}"
: "${SWEEP_LIB_CHEN:?}"
export CHEN_OLD=1 R_LIBS=$SWEEP_LIB_CHEN OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
C=$SWEEP_CODE/work_net/CHEN2021; W=$SWEEP_OUT/work_net/CHEN2021; mkdir -p $W
for m in m0 hist; do
  echo "$(date) start $m" >> $W/mcmle_queue.log
  Rscript $C/chen_run.R $m MCMLE 2 2000 4000 8 > $W/log_${m}_MCMLE.txt 2>&1
  echo "$(date) end $m rc=$?" >> $W/mcmle_queue.log
done
# MPLE point estimates (author stack) for all four variants; seconds each
for m in m0 hist m26 placebo; do Rscript $C/chen_run.R $m MPLE 1 > $W/log_${m}_MPLE.txt 2>&1; done
echo DONE >> $W/mcmle_queue.log
