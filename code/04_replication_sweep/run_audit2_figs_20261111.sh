#!/usr/bin/env bash
# Last two steps of the 2026-11-11 finish run (05 expected counts updated), appended to the finish log.
#   bash ~/Dropbox/ottoman_byzantine/code/04_replication_sweep/run_audit2_figs_20261111.sh
cd "$HOME/Dropbox/ottoman_byzantine" || exit 1
{
  set -x
  Rscript code/03_figures/05_appendix_compact_figures.R || exit 1
  Rscript code/03_figures/02_case_study_map.R || exit 1
  echo "ALL STEPS FINISHED"
} 2>&1 | tee -a logs/audit2_finish_20261111.log
