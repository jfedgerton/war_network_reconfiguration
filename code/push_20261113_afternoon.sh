#!/bin/zsh
# Push the 2026-11-13 afternoon changes (two added replication studies, figure titles, stars) to the release repo.
#   zsh ~/Dropbox/ottoman_byzantine/code/push_20261113_afternoon.sh
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
       code/04_replication_sweep/registry.csv
       code/04_replication_sweep/extract_map.csv
       code/04_replication_sweep/sens_studies.csv
       code/04_replication_sweep/study_table_labels.csv
       code/04_replication_sweep/study_flow.csv
       code/04_replication_sweep/sweep_studies.bib
       code/04_replication_sweep/expected_est.rds
       code/04_replication_sweep/legacy/work_add20261113/inv34/rerun_inv34.R
       code/04_replication_sweep/legacy/work_add20261113/inv8/rebuild_inv8.R
       code/04_replication_sweep/run_add_20261113.sh
       code/push_20261113_afternoon.sh)
for f in $FILES; do [[ -f $P/$f ]] || { echo "STOP: missing $P/$f"; exit 1; }; mkdir -p $(dirname $f); cp -p $P/$f $f; done
git add $FILES
git diff --cached --quiet && { echo "Nothing new to commit"; exit 0; }
git commit -F - <<'MSG'
Replication sweep: add Suhas 2023 (author data) and Long et al. 2007 (ATOP v3.0 rebuild); pool 140 -> 142

- inv34 Suhas 2023 ISQ: author's ISQ supplemental data; Table 5.0 model 1 reproduced exactly; common-sample exception 0.75
- inv8 Long, Nordstrom & Baek 2007 JOP: protocol-B rebuild from ATOP v3.0 + GML MIDs (failure = onset of a mutual-force MID); REBUILT_MATCH on non-aggression
- remaining five pending articles reclassified as no replication data
- 95: pooled count assertion 142; 91/93/32/07c: in-plot titles removed; 05: Table stars match main text (* .05, ** .01, *** .001)
- Raw inputs (not in repo): data/raw/replication_sweep/_incoming/Suhas2023 (OUP sqad070 supplemental zip), _incoming/Long2007 (ATOP v3.0 from atopdata.org)

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01Dof3mMRXGfgtxqVmTQMu55
MSG
git push || { echo "STOP: push failed"; exit 1; }
echo "PUSHED $(git rev-parse --short HEAD)"
