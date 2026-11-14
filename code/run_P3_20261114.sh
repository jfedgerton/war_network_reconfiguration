#!/usr/bin/env bash
# [2026-11-14] Pipeline audit P3 fixes: P9_080 (missing ally power change -> NA), then summaries/tables/figures.
# P8_030 and 100 moved to outcome_overlap in the registry (no rerun needed).
cd "$HOME/Dropbox/ottoman_byzantine" || exit 1
LOG=logs/P3_20261114.log
{
  set -x
  for ID in P9_080; do
    if ! Rscript code/04_replication_sweep/run_one.R $ID; then
      REBASELINE_OK=$ID Rscript code/04_replication_sweep/92b_rebaseline_one.R $ID || { echo "STOP: rebaseline $ID"; exit 1; }
      Rscript code/04_replication_sweep/run_one.R $ID || { echo "STOP: run_one after rebaseline $ID"; exit 1; }
    fi
    if grep -q "^$ID," code/04_replication_sweep/sens_studies.csv; then ONLY_IDS=$ID Rscript code/04_replication_sweep/run_sensitivity.R || { echo "STOP: sensitivity $ID"; exit 1; }; fi
  done
  for s in 90_summarise 90b_broader_rule 90c_expanded_set 91_figures 93_sensitivity_summary 94_study_flow_table; do Rscript code/04_replication_sweep/$s.R || { echo "STOP at $s"; exit 1; }; done
  (cd code/04_replication_sweep && Rscript 95_study_table.R) || { echo "STOP at 95"; exit 1; }
  EXPECT_SKIP=1 Rscript code/03_figures/05_appendix_compact_figures.R > logs/05_P3_20261114.log 2>&1 || { echo "STOP at 05 (see logs/05_P3_20261114.log)"; exit 1; }
  grep -n "counts differ\|median" logs/05_P3_20261114.log | head
  echo "ALL STEPS FINISHED"
} 2>&1 | tee "$LOG"
