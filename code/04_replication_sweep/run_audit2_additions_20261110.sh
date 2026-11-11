#!/usr/bin/env bash
# Rerun after the 2026-11-10 two-audit fixes and the targeted additions. Run from anywhere:
#   bash ~/Dropbox/ottoman_byzantine/code/04_replication_sweep/run_audit2_additions_20261110.sh
cd "$HOME/Dropbox/ottoman_byzantine" || exit 1
LOG=logs/audit2_additions_20261110.log
{
  set -x
  for ID in ADD_T01 ADD_T03 ADD_T05 ADD_D10 ADD_D11 ADD_D01; do
    Rscript code/04_replication_sweep/run_one.R "$ID"
  done
  # P9_078: new GEE estimator -> first run fails by design, rebaseline, then must pass
  Rscript code/04_replication_sweep/run_one.R P9_078
  REBASELINE_OK=P9_078 Rscript code/04_replication_sweep/92b_rebaseline_one.R P9_078
  Rscript code/04_replication_sweep/run_one.R P9_078 || { echo "P9_078 did not pass after rebaseline"; exit 1; }
  ONLY_IDS=C2,P8_019,P8_031,X020,X022,X025,X077,X078,X092,X074,P9_078,ADD_T01,ADD_T05,ADD_D10 \
    Rscript code/04_replication_sweep/run_sensitivity.R || { echo "sensitivity failed"; exit 1; }
  Rscript code/04_replication_sweep/90_summarise.R || exit 1
  Rscript code/04_replication_sweep/91_figures.R || exit 1
  Rscript code/04_replication_sweep/93_sensitivity_summary.R || exit 1
  Rscript code/04_replication_sweep/94_study_flow_table.R || exit 1
  Rscript code/03_figures/05_appendix_compact_figures.R || exit 1
  Rscript code/03_figures/02_case_study_map.R || exit 1
  echo "ALL STEPS FINISHED"
} 2>&1 | tee "$LOG"
