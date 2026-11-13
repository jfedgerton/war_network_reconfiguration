#!/bin/zsh
# [2026-11-13] Rerun 07c (Fig. 1 heatmap now cividis) and push today's code changes to the release repo.
#   zsh ~/Dropbox/ottoman_byzantine/code/run_push_cividis_20261113.sh
set -u
P=~/Dropbox/ottoman_byzantine; R=$P/github_release
cd $P || exit 1
mkdir -p logs
Rscript code/02_analysis/07c_multiverse_figures_bellig.R > logs/07c_cividis_20261113.log 2>&1 || { echo "STOP: 07c failed, see logs/07c_cividis_20261113.log"; tail -20 logs/07c_cividis_20261113.log; exit 1; }
ls -l code/03_figures/output/fig_multiverse_hypotheses_bellig.png
echo "07c OK"

cd $R || exit 1
git pull --ff-only || { echo "STOP: pull failed"; exit 1; }
FILES=(code/02_analysis/07c_multiverse_figures_bellig.R
       code/02_analysis/01_system_level_analysis.R
       code/02_analysis/02_state_level_analysis.R
       code/02_analysis/19_contamination_multiverse.R
       code/02_analysis/32_multiverse_permutation.R
       code/run_perm_alt_20261113.sh
       code/run_push_cividis_20261113.sh)
for f in $FILES; do [[ -f $P/$f ]] || { echo "STOP: missing $P/$f"; exit 1; }; mkdir -p $(dirname $f); cp -p $P/$f $f; done
git add $FILES
git diff --cached --quiet && { echo "Nothing new to commit"; exit 0; }
git commit -F - <<'MSG'
Fig. 1 heatmap in cividis; table notes, signed contamination prediction, permutation fix

- 07c: multiverse heatmap uses a sequential cividis scale (0-100) instead of red-white-blue diverging at 50; white labels on dark cells
- 01, 02: shortened table notes to match the manuscript
- 19: contamination multiverse scores the prediction on the signed shift (prediction_held); old label kept as label_held
- 32: bellig_rule column added when absent before the figure step (fixes semi_join error)

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01Dof3mMRXGfgtxqVmTQMu55
MSG
git push || { echo "STOP: push failed"; exit 1; }
echo "PUSHED $(git rev-parse --short HEAD)"
