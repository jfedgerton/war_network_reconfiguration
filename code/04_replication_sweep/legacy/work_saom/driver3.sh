#!/bin/bash
# Ported driver: paths come from SWEEP_CODE / SWEEP_OUT (SWEEP_DATA must also be set for the R scripts);
# flock on the old sweep R.lock removed (the pipeline caller serialises R jobs). Run logs, frozen script copies and
# all fits go to $SWEEP_OUT/work_saom/<id>/. Not called by the pipeline; legacy/_slow/<id>.R reproduces these calls.
: "${SWEEP_CODE:?}" "${SWEEP_DATA:?}" "${SWEEP_OUT:?}"
RS=Rscript
C=$SWEEP_CODE/work_saom
W=$SWEEP_OUT/work_saom
# original first line waited for a running PID (while kill -0 9080 ...); dropped.
# Original X085 m0 part 1 (log renamed run_m0_part1.out afterwards). It was stopped by hand after ladder rep 6
# (n2start 800 accepted, switch to n3 = 5000 announced, ~1.65 h); driver4.sh then resumed from fit_m0_best.rds.
run() { d=$1; s=$2; m=$3; mkdir -p $W/$d && cp $C/$d/$s $W/$d/frozen_${m}.R && cd $W/$d && CAP_H=4 REP_MIN=40 $RS frozen_${m}.R $m > run_${m}.out 2>&1; }
run X085 rerun_X085.R m0
run X085 rerun_X085.R mH
echo ALLDONE > $W/driver.done
