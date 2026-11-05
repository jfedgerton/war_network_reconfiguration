#!/bin/bash
# Ported driver: paths come from SWEEP_CODE / SWEEP_OUT (SWEEP_DATA must also be set for the R scripts);
# flock on the old sweep R.lock removed (the pipeline caller serialises R jobs). Run logs, frozen script copies and
# all fits go to $SWEEP_OUT/work_saom/<id>/. Not called by the pipeline; legacy/_slow/<id>.R reproduces these calls.
: "${SWEEP_CODE:?}" "${SWEEP_DATA:?}" "${SWEEP_OUT:?}"
RS=Rscript
C=$SWEEP_CODE/work_saom
W=$SWEEP_OUT/work_saom
mkdir -p $W/X085
cd $W/X085 && cp $C/X085/rerun_X085.R frozen_m0b.R && RESUME=1 N2SW=1 N3SW=1 CAP_H=3 REP_MIN=90 $RS frozen_m0b.R m0 > run_m0_part2.out 2>&1
cd $W/X085 && cp $C/X085/rerun_X085.R frozen_mH.R && CAP_H=4.5 REP_MIN=90 $RS frozen_mH.R mH > run_mH.out 2>&1
echo ALLDONE > $W/driver.done
