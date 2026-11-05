#!/bin/sh
# resumable queue: skips a run whose log already has the final "elapsed" line
# Ported: paths from SWEEP_CODE / SWEEP_DATA / SWEEP_OUT; flock on the old sweep R.lock removed (the pipeline caller
# serialises R jobs); logs and outputs go under $SWEEP_OUT/work_net. Not called by the pipeline; legacy/_slow/<id>.R
# reproduces these calls.
: "${SWEEP_CODE:?}" "${SWEEP_DATA:?}" "${SWEEP_OUT:?}"
: "${SWEEP_LIB_AMEN:?}"   # amen 1.4 (pa2018_version) library for run_minhas.R
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
RS=Rscript
C=$SWEEP_CODE/work_net; B=$SWEEP_OUT/work_net; D=$SWEEP_DATA/work_net
run() { d=$1; shift; log=$1; shift; s=$1; shift; mkdir -p $B/$d; cd $B/$d
  if grep -q "^elapsed\|^\s*elapsed\|Total time" $log 2>/dev/null && [ -f $log.ok ]; then echo "$(date) skip $log" >> $B/queue6.log; return; fi
  echo "$(date) start $log" >> $B/queue6.log; $RS $C/$d/$s "$@" > $log 2>&1; rc=$?; [ $rc -eq 0 ] && touch $log.ok; echo "$(date) end $log rc=$rc" >> $B/queue6.log; }
run MINHAS2022 log_weeks_m0.txt run_minhas.R weeks m0 500 2500 10 1950
run MINHAS2022 log_weeks_hist.txt run_minhas.R weeks hist 500 2500 10 1950
# the two WSC runs below used an earlier cold-start run_wsc.R; they were discarded (killed; see queue6.log) and rerun by queue7.sh
run WSC2007 log_wsc_m26.txt run_wsc.R m26 25000
run WSC2007 log_wsc_placebo.txt run_wsc.R placebo 25000
run MINHAS2022 log_rs_m26.txt run_minhas.R rs m26 500 3000 10 1945
run MINHAS2022 log_rs_placebo.txt run_minhas.R rs placebo 500 3000 10 1945
run WH2011 log_wh_m26.txt probit_wh.R $D/WH2011/wh_m26.csv out_m26 12000 3000 12
run WH2011 log_wh_placebo.txt probit_wh.R $D/WH2011/wh_placebo.csv out_placebo 12000 3000 12
run MINHAS2022 log_weeks_m26.txt run_minhas.R weeks m26 500 2500 10 1950
run MINHAS2022 log_weeks_placebo.txt run_minhas.R weeks placebo 500 2500 10 1950
echo "$(date) Q6DONE" >> $B/queue6.log
