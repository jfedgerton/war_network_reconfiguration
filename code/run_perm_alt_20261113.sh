#!/bin/bash
# [2026-11-13] Multiverse permutation (500 draws) with state draws on churn_alt only -> mv_perm_alt/, *_alt outputs.
# Old draws in data/processed/mv_perm/ are untouched. Resumable: rerun the same command if interrupted.
set -euo pipefail
cd ~/Dropbox/ottoman_byzantine
mkdir -p logs
MV_PERM_ALT=1 Rscript code/02_analysis/32_multiverse_permutation.R > logs/32_alt_20261113.log 2>&1
Rscript code/03_figures/06_appendix_grid_fig3_4.R > logs/06_alt_20261113.log 2>&1
cat results/appendix/multiverse_permutation_summary_alt.csv
echo DONE
