#!/usr/bin/env bash
# [2026-11-13] Follow-up timing fixes found by the final spec audit (history anchored on the t+1 outcome year): X092, 12, 63, 11.
# Run AFTER run_final_20261113.sh has finished. Then rebuilds all summaries and tables.
cd "$HOME/Dropbox/ottoman_byzantine" || exit 1
LOG=logs/final3_20261113.log
{
  set -x
  for ID in 12 63 11; do
    Rscript code/04_replication_sweep/run_one.R $ID
    REBASELINE_OK=$ID Rscript code/04_replication_sweep/92b_rebaseline_one.R $ID || { echo "STOP: rebaseline $ID"; exit 1; }
    Rscript code/04_replication_sweep/run_one.R $ID || { echo "STOP: run_one after rebaseline $ID"; exit 1; }
    if grep -q "^$ID," code/04_replication_sweep/sens_studies.csv; then ONLY_IDS=$ID Rscript code/04_replication_sweep/run_sensitivity.R || { echo "STOP: sensitivity $ID"; exit 1; }; fi
  done
  for s in 90_summarise 90b_broader_rule 91_figures 93_sensitivity_summary 94_study_flow_table; do Rscript code/04_replication_sweep/$s.R || { echo "STOP at $s"; exit 1; }; done
  (cd code/04_replication_sweep && Rscript 95_study_table.R) || { echo "STOP at 95"; exit 1; }
  Rscript code/03_figures/05_appendix_compact_figures.R > logs/05_final3_20261113.log 2>&1 || { echo "STOP at 05"; exit 1; }
  echo "ALL STEPS FINISHED"
} 2>&1 | tee "$LOG"
