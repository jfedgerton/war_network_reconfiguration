#!/usr/bin/env bash
# Finish the 2026-11-10 rerun: ADD_D11 (needed statmod) and the P9_078 sensitivity hook (N fix), then summaries and figures.
#   bash ~/Dropbox/ottoman_byzantine/code/04_replication_sweep/run_audit2_finish_20261111.sh
cd "$HOME/Dropbox/ottoman_byzantine" || exit 1
LOG=logs/audit2_finish_20261111.log
{
  set -x
  Rscript -e 'if (!requireNamespace("statmod", quietly = TRUE)) install.packages("statmod", repos = "https://cloud.r-project.org"); stopifnot(requireNamespace("statmod"))' || exit 1
  Rscript code/04_replication_sweep/run_one.R ADD_D11 || { echo "ADD_D11 failed"; exit 1; }
  ONLY_IDS=P9_078 Rscript code/04_replication_sweep/run_sensitivity.R || { echo "sensitivity failed"; exit 1; }
  Rscript code/04_replication_sweep/90_summarise.R || exit 1
  Rscript code/04_replication_sweep/91_figures.R || exit 1
  Rscript code/04_replication_sweep/93_sensitivity_summary.R || exit 1
  Rscript code/04_replication_sweep/94_study_flow_table.R || exit 1
  Rscript code/03_figures/05_appendix_compact_figures.R || exit 1
  Rscript code/03_figures/02_case_study_map.R || exit 1
  echo "ALL STEPS FINISHED"
} 2>&1 | tee "$LOG"
