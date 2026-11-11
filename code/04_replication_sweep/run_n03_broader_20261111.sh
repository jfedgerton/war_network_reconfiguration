#!/usr/bin/env bash
# N03 back in the pool (conflict-outcome term) + broader-rule appendix count; summaries and figures only (no re-estimation).
#   bash ~/Dropbox/ottoman_byzantine/code/04_replication_sweep/run_n03_broader_20261111.sh
cd "$HOME/Dropbox/ottoman_byzantine" || exit 1
{
  set -x
  Rscript code/04_replication_sweep/90_summarise.R || exit 1
  Rscript code/04_replication_sweep/90b_broader_rule.R || exit 1
  Rscript code/04_replication_sweep/91_figures.R || exit 1
  Rscript code/04_replication_sweep/94_study_flow_table.R || exit 1
  Rscript code/03_figures/05_appendix_compact_figures.R || exit 1
  echo "ALL STEPS FINISHED"
} 2>&1 | tee logs/n03_broader_20261111.log
