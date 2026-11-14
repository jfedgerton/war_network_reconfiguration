#!/usr/bin/env bash
# [2026-11-14] Pipeline audit P1 fixes + reduced main pool (outcome-overlap studies -> appendix expanded set).
#  Reruns: P8_031 (t+1 lead anchor), 57 (printed control added), P9_029 (single affinity lag), P8_048 (Developed dummy),
#  X040 (independence working correlation). Kinne (11) sample==1 is a TEST only (prints, writes nothing).
cd "$HOME/Dropbox/ottoman_byzantine" || exit 1
LOG=logs/restructure_20261114_b.log
{
  set -x
  for ID in 57 P9_029 P8_048 X040; do   # P8_031 finished in the first pass
    if ! Rscript code/04_replication_sweep/run_one.R $ID; then
      REBASELINE_OK=$ID Rscript code/04_replication_sweep/92b_rebaseline_one.R $ID || { echo "STOP: rebaseline $ID"; exit 1; }
      Rscript code/04_replication_sweep/run_one.R $ID || { echo "STOP: run_one after rebaseline $ID"; exit 1; }
    fi
    if grep -q "^$ID," code/04_replication_sweep/sens_studies.csv; then ONLY_IDS=$ID Rscript code/04_replication_sweep/run_sensitivity.R || { echo "STOP: sensitivity $ID"; exit 1; }; fi
  done
  grep -n "NOT_REPRODUCED\|REPRODUCED" logs/replication_sweep/57.log | tail -3
  SWEEP_CODE=$PWD/code/04_replication_sweep/legacy SWEEP_DATA=$PWD/data/raw/replication_sweep SWEEP_OUT=$PWD/tmp/test11 \
    SWEEP_HIST=$PWD/data/processed/replication_sweep/hist/history_tables.rds \
    Rscript code/04_replication_sweep/legacy/work/11/test_sample1_11.R > logs/test11_20261114.log 2>&1 || echo "KINNE TEST FAILED (see logs/test11_20261114.log)"
  tail -12 logs/test11_20261114.log
  for s in 90_summarise 90b_broader_rule 90c_expanded_set 91_figures 93_sensitivity_summary 94_study_flow_table; do Rscript code/04_replication_sweep/$s.R || { echo "STOP at $s"; exit 1; }; done
  (cd code/04_replication_sweep && Rscript 95_study_table.R) || { echo "STOP at 95"; exit 1; }
  EXPECT_SKIP=1 Rscript code/03_figures/05_appendix_compact_figures.R > logs/05_restructure_20261114.log 2>&1 || { echo "STOP at 05 (see logs/05_restructure_20261114.log)"; exit 1; }
  grep -n "counts differ\|median" logs/05_restructure_20261114.log | head
  echo "ALL STEPS FINISHED"
} 2>&1 | tee "$LOG"
