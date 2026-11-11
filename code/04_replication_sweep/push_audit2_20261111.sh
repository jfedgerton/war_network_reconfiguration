#!/bin/zsh
# Sync the 2026-11-10/11 replication-sweep audit fixes and additions to the GitHub release clone, commit and push.
# Run only after run_audit2_finish_20261111.sh printed "ALL STEPS FINISHED".
#   zsh ~/Dropbox/ottoman_byzantine/code/04_replication_sweep/push_audit2_20261111.sh
set -u
P=~/Dropbox/ottoman_byzantine; R=$P/github_release
grep -q "ALL STEPS FINISHED" $P/logs/audit2_finish_20261111.log 2>/dev/null || { echo "STOP: rerun did not finish cleanly"; exit 1; }
cd $R || exit 1
git pull --ff-only || { echo "STOP: pull failed"; exit 1; }
S=code/04_replication_sweep
FILES=(
  $S/registry.csv $S/study_flow.csv $S/sens_studies.csv $S/extract_map.csv $S/expected_est.rds
  $S/apply_audit2_20261110.py $S/run_audit2_additions_20261110.sh $S/run_audit2_finish_20261111.sh $S/push_audit2_20261111.sh $S/run_audit2_figs_20261111.sh code/03_figures/05_appendix_compact_figures.R
  $S/legacy/work_p8/P8_019/rerun_P8_019.R $S/legacy/work_p4/X020/rerun_X020.R $S/legacy/work_p4/X025/rerun_est_X025.R
  $S/legacy/work_demcap/C2/rerun_C2.R $S/legacy/work_p8/P8_031/rerun_P8_031.R $S/legacy/work_p4/X022/rerun_X022.R
  $S/legacy/work_p4/X077/rerun_est_X077.R $S/legacy/work_p4/X078/rerun_X078.R $S/legacy/work_p4/X092/rerun_est_X092.R
  $S/legacy/work_p4/X074/rebuild_X074.R $S/legacy/work_p9/P9_078/rebuild_P9_078.R
  $S/legacy/work_add20261110/ADD_T01/rerun_ADD_T01.R $S/legacy/work_add20261110/ADD_T03/rerun_ADD_T03.R
  $S/legacy/work_add20261110/ADD_T05/rerun_ADD_T05.R $S/legacy/work_add20261110/ADD_D01/rerun_ADD_D01.R
  $S/legacy/work_add20261110/ADD_D10/rerun_ADD_D10.R $S/legacy/work_add20261110/ADD_D11/rerun_ADD_D11.R
)
for f in $FILES; do
  [[ -f $P/$f ]] || { echo "STOP: missing $P/$f"; exit 1; }
  mkdir -p $(dirname $f); cp -p $P/$f $f
done
git add $FILES
git status --short
git commit -F - <<'MSG'
Replication sweep: two-audit adjudication, sensitivity hooks, targeted additions

- Adjudicated two independent replication audits: headline terms corrected (P8_019, X020, X025, C2, P8_042,
  X022, P8_031, N08), family corrections (X085, X043), out-of-scope and duplicate exclusions, P9_078 GEE AR(1).
- Sensitivity hooks moved to the current headline terms; X077 split-only.
- Six targeted additions (ADD_T01, ADD_T03, ADD_T05, ADD_D01, ADD_D10, ADD_D11) plus study-flow rows
  for the candidates without data. Pool: 139 articles (88 tie, 40 democracy, 11 capitalism).

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01Dof3mMRXGfgtxqVmTQMu55
MSG
git push || { echo "STOP: push failed"; exit 1; }
echo "PUSHED"
