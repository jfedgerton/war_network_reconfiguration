#!/bin/zsh
# 2026-11-12: redraw appendix Figs 11-12 (fatal MIDs only) and rerun 19 for the 1x2 contamination figure.
# Contamination analyses stay on ALL MIDs (author decision 2026-11-12). Supersedes run_contam_fig_20261112.sh.
#   11 war samples, new state fatal-year cut points (Fig 11)  ~minutes
#   05 compact appendix figures (B2 Fig 6, B3 Fig 11, B5 Fig 12) ~minutes
#   19 contamination multiverse (appendix Fig 4)              ~2-4 h, run LAST
# Old outputs are copied to tmp/backup_20261112_fatal/outputs/ first. Then the changed scripts are pushed.
#   zsh ~/Dropbox/ottoman_byzantine/code/02_analysis/run_fatal_20261112.sh
set -u
P=~/Dropbox/ottoman_byzantine; R=$P/github_release
LOG=$P/logs/fatal_20261112.log
mkdir -p $P/logs $P/tmp/backup_20261112_fatal/outputs
cd $P || exit 1
for f in data/processed/contamination_*.csv data/processed/robustness_war_samples.csv \
         code/03_figures/output/manuscript/fig_ms_contamination_*.png code/03_figures/output/manuscript/fig_app_compact_*.png; do
  [[ -f $f ]] && cp -p $f $P/tmp/backup_20261112_fatal/outputs/
done
echo "START $(date)" | tee $LOG
run() { echo "== $1  $(date)" | tee -a $LOG; Rscript $1 >> $LOG 2>&1 || { echo "STOP: $1 failed (see $LOG)"; exit 1; }; }
run code/02_analysis/11_robustness_war_samples.R
run code/03_figures/05_appendix_compact_figures.R
run code/02_analysis/19_contamination_multiverse.R
echo "ALL STEPS FINISHED $(date)" | tee -a $LOG

cd $R || exit 1
git pull --ff-only || { echo "STOP: pull failed"; exit 1; }
FILES=(code/02_analysis/19_contamination_multiverse.R code/02_analysis/11_robustness_war_samples.R code/03_figures/05_appendix_compact_figures.R
       code/02_analysis/run_fatal_20261112.sh)   # outputs (csv, png) are gitignored in the release repo
for f in $FILES; do [[ -f $P/$f ]] || { echo "STOP: missing $P/$f"; exit 1; }; mkdir -p $(dirname $f); cp -p $P/$f $f; done
git add $FILES
git commit -F - <<'MSG'
Appendix figures: single-row fatal war samples, 2x3 IGO-share grid, 1x2 contamination multiverse

- 19: contamination multiverse figure one column x two rows
- 11: state-level samples for years with >= 3/5/10 fatal MIDs system-wide
- 05: B3 single row, log(1 + fatal MIDs) only; new B5 IGO-share grid (fatal only)

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01Dof3mMRXGfgtxqVmTQMu55
MSG
git push || { echo "STOP: push failed"; exit 1; }
echo "PUSHED"
open $P/code/03_figures/output/manuscript
