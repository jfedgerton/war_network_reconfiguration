#!/usr/bin/env bash
# [2026-11-14] Summaries/tables/figures after the restructure (studies already rerun; Powers 57 moved to not_reproduced).
cd "$HOME/Dropbox/ottoman_byzantine" || exit 1
LOG=logs/restructure_20261114_d.log
{
  set -x
  for s in 90_summarise 90b_broader_rule 90c_expanded_set 91_figures 93_sensitivity_summary 94_study_flow_table; do Rscript code/04_replication_sweep/$s.R || { echo "STOP at $s"; exit 1; }; done
  (cd code/04_replication_sweep && Rscript 95_study_table.R) || { echo "STOP at 95"; exit 1; }
  EXPECT_SKIP=1 Rscript code/03_figures/05_appendix_compact_figures.R > logs/05_restructure_20261114.log 2>&1 || { echo "STOP at 05 (see logs/05_restructure_20261114.log)"; exit 1; }
  grep -n "counts differ\|median" logs/05_restructure_20261114.log | head
  echo "ALL STEPS FINISHED"
} 2>&1 | tee "$LOG"
