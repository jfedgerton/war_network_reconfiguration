#!/bin/zsh
# Push the 2026-11-13 final audit changes (pool 126, drops, term switches, history-timing fixes, labels) to the release repo.
#   zsh ~/Dropbox/ottoman_byzantine/code/push_20261113_final.sh
set -u
P=~/Dropbox/ottoman_byzantine; R=$P/github_release
cd $R || exit 1
git pull --ff-only || { echo "STOP: pull failed"; exit 1; }
FILES=(code/02_analysis/07c_multiverse_figures_bellig.R
       code/02_analysis/32_multiverse_permutation.R
       code/03_figures/05_appendix_compact_figures.R
       code/04_replication_sweep/91_figures.R
       code/04_replication_sweep/93_sensitivity_summary.R
       code/04_replication_sweep/95_study_table.R
       code/04_replication_sweep/drop_20261113.py
       code/04_replication_sweep/registry.csv
       code/04_replication_sweep/extract_map.csv
       code/04_replication_sweep/sens_studies.csv
       code/04_replication_sweep/study_table_labels.csv
       code/04_replication_sweep/study_flow.csv
       code/04_replication_sweep/sweep_studies.bib
       code/04_replication_sweep/expected_est.rds
       code/04_replication_sweep/legacy/work_add20261113/inv34/rerun_inv34.R
       code/04_replication_sweep/legacy/work_add20261113/inv8/rebuild_inv8.R
       code/04_replication_sweep/legacy/work/103/rerun_103.R
       code/04_replication_sweep/legacy/work/52/rebuild_52.R
       code/04_replication_sweep/legacy/work_p4/X040/rerun_verify_X040.R
       code/04_replication_sweep/legacy/work_demcap/D2/rerun_D2.R
       code/04_replication_sweep/legacy/work_p8/P8_003/rerun_P8_003.R
       code/04_replication_sweep/legacy/work_p9/P9_028/rebuild_P9_028.R
       code/04_replication_sweep/legacy/work_p9/P9_029/rebuild_P9_029.R
       code/04_replication_sweep/legacy/work_p9/P9_061/rebuild_P9_061.R
       code/04_replication_sweep/legacy/work_p9/P9_089/rebuild_P9_089.R
       code/04_replication_sweep/legacy/work_p4/X092/rerun_est_X092.R
       code/04_replication_sweep/legacy/work/12/rerun_12.R
       code/04_replication_sweep/legacy/work/63/rerun_est_63.R
       code/04_replication_sweep/legacy/work/11/rerun_est_11.R
       code/04_replication_sweep/run_add_20261113.sh
       code/run_final_20261113.sh code/run_final2_20261113.sh code/run_final3_20261113.sh code/run_final4_20261113.sh code/run_final5_20261113.sh
       code/push_20261113_final.sh)
for f in $FILES; do [[ -f $P/$f ]] || { echo "STOP: missing $P/$f"; exit 1; }; mkdir -p $(dirname $f); cp -p $P/$f $f; done
git add $FILES
git diff --cached --quiet && { echo "Nothing new to commit"; exit 0; }
git commit -F - <<'MSG'
Replication sweep final audit: pool 126 (80 tie / 36 dem / 10 cap)

- Drop out-of-rule studies (ICB-only, IMI, resistance/reciprocation, compellent threats, composite threat, methods/measurement papers, Benson 2004)
- Switch to headline terms: Bak 2018, Bayer & Bernhard 2010 (dichotomous), Daxecker 2007 (H2 interaction)
- Anchor conflict history on the outcome year where the DV is a t+1 lead (verified against COW dyadic MID data)
- Table 9 labels, units, outcomes and published values corrected against the PDFs; rebuilt studies flagged

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01Dof3mmrxgfgtxqvmtqmu55
MSG
git push || { echo "STOP: push failed"; exit 1; }
echo "PUSHED"
