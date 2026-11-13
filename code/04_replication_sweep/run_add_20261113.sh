#!/usr/bin/env bash
# [2026-11-13] Add inv34 (Suhas 2023, author data) and inv8 (Long et al. 2007, ATOP v3.0 rebuild) to the pool (140 -> 142),
# reclassify the five remaining pending articles as no data, rerun the sweep summaries, figures and tables.
#   bash ~/Dropbox/ottoman_byzantine/code/04_replication_sweep/run_add_20261113.sh
cd "$HOME/Dropbox/ottoman_byzantine" || exit 1
mkdir -p logs
LOG=logs/add_20261113.log
{
  set -x
  for ID in inv34 inv8; do
    Rscript code/04_replication_sweep/run_one.R "$ID" || { echo "STOP: run_one $ID failed (see logs/replication_sweep/$ID.log)"; exit 1; }
  done
  ONLY_IDS=inv34,inv8 Rscript code/04_replication_sweep/run_sensitivity.R || { echo "STOP: sensitivity failed"; exit 1; }
  Rscript code/04_replication_sweep/90_summarise.R || exit 1
  Rscript code/04_replication_sweep/90b_broader_rule.R || exit 1
  Rscript code/04_replication_sweep/91_figures.R || exit 1
  Rscript code/04_replication_sweep/93_sensitivity_summary.R || exit 1
  Rscript code/04_replication_sweep/94_study_flow_table.R || exit 1
  Rscript code/04_replication_sweep/95_study_table.R || exit 1
  Rscript code/03_figures/05_appendix_compact_figures.R || exit 1
  # in-plot titles removed (B31): regenerate the multiverse figures and the permutation grid (draws are cached; no re-estimation)
  Rscript code/02_analysis/07c_multiverse_figures_bellig.R > logs/07c_titles_20261113.log 2>&1 || { echo "STOP: 07c failed"; exit 1; }
  MV_PERM_ALT=1 Rscript code/02_analysis/32_multiverse_permutation.R > logs/32_titles_20261113.log 2>&1 || { echo "STOP: 32 failed"; exit 1; }
  Rscript code/03_figures/06_appendix_grid_fig3_4.R > logs/06_titles_20261113.log 2>&1 || { echo "STOP: 06 failed"; exit 1; }
  echo "ALL STEPS FINISHED"
} 2>&1 | tee "$LOG"
