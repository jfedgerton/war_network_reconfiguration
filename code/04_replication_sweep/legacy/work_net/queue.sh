#!/bin/sh
# Ported: paths from SWEEP_CODE / SWEEP_DATA / SWEEP_OUT; flock on the old sweep R.lock removed (the pipeline caller
# serialises R jobs); logs and outputs go under $SWEEP_OUT/work_net. Not called by the pipeline; legacy/_slow/<id>.R
# reproduces these calls.
: "${SWEEP_CODE:?}" "${SWEEP_DATA:?}" "${SWEEP_OUT:?}"
: "${SWEEP_LIB_AMEN:?}"   # amen 1.4 (pa2018_version) library for run_minhas.R
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
RS=Rscript
C=$SWEEP_CODE/work_net; B=$SWEEP_OUT/work_net; D=$SWEEP_DATA/work_net
run() { d=$1; shift; log=$1; shift; s=$1; shift; mkdir -p $B/$d; cd $B/$d; echo "$(date) start $log" >> $B/queue.log; $RS $C/$d/$s "$@" > $log 2>&1; echo "$(date) end $log rc=$?" >> $B/queue.log; }
# the two WSC runs used an earlier cold-start run_wsc.R (outputs later renamed coldstart_*, discarded); final WSC chains: queue7.sh.
# rs m0/hist and wh m0/hist below are the final fits; this queue wrote no .ok markers (they were touched by hand before queue6/7).
run WSC2007 log_wsc_m0.txt run_wsc.R m0 25000
run WSC2007 log_wsc_hist.txt run_wsc.R hist 25000
run MINHAS2022 log_rs_m0.txt run_minhas.R rs m0 500 3000 10 1945
run MINHAS2022 log_rs_hist.txt run_minhas.R rs hist 500 3000 10 1945
run WH2011 log_wh_m0.txt probit_wh.R $D/WH2011/wh_m0.csv out_m0 12000 3000 12
run WH2011 log_wh_hist.txt probit_wh.R $D/WH2011/wh_hist.csv out_hist 12000 3000 12
run MINHAS2022 log_weeks_m0.txt run_minhas.R weeks m0 500 2500 10 1950
run MINHAS2022 log_weeks_hist.txt run_minhas.R weeks hist 500 2500 10 1950
echo QDONE >> $B/queue.log
