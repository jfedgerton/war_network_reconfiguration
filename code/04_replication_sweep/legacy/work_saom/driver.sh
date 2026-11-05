#!/bin/bash
# Sequential SAOM runs (original: first attempt; X085 m0 of this driver was superseded by driver3/driver4).
# Each run executes a frozen copy of the script (R --file parses incrementally; never edit a running script).
# Ported driver: paths come from SWEEP_CODE / SWEEP_OUT (SWEEP_DATA must also be set for the R scripts);
# flock on the old sweep R.lock removed (the pipeline caller serialises R jobs). Run logs, frozen script copies and
# all fits go to $SWEEP_OUT/work_saom/<id>/. Not called by the pipeline; legacy/_slow/<id>.R reproduces these calls.
: "${SWEEP_CODE:?}" "${SWEEP_DATA:?}" "${SWEEP_OUT:?}"
RS=Rscript
C=$SWEEP_CODE/work_saom
W=$SWEEP_OUT/work_saom
run() { d=$1; s=$2; m=$3; mkdir -p $W/$d && cp $C/$d/$s $W/$d/frozen_${m}.R && cd $W/$d && CAP_H=5 $RS frozen_${m}.R $m > run_${m}.out 2>&1; }
run X060 rerun_X060.R m0
run X085 rerun_X085.R m0
run X060 rerun_X060.R mH
run X085 rerun_X085.R mH
echo ALLDONE > $W/driver.done
