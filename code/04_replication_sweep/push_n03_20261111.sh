#!/bin/zsh
# Push the N03 headline switch and the broader-rule appendix count to the GitHub release clone.
#   zsh ~/Dropbox/ottoman_byzantine/code/04_replication_sweep/push_n03_20261111.sh
set -u
P=~/Dropbox/ottoman_byzantine; R=$P/github_release
grep -q "ALL STEPS FINISHED" $P/logs/n03_broader_20261111.log 2>/dev/null || { echo "STOP: rerun did not finish cleanly"; exit 1; }
cd $R || exit 1
git pull --ff-only || { echo "STOP: pull failed"; exit 1; }
S=code/04_replication_sweep
FILES=($S/registry.csv $S/study_flow.csv $S/90b_broader_rule.R $S/run_n03_broader_20261111.sh $S/push_n03_20261111.sh)
for f in $FILES; do [[ -f $P/$f ]] || { echo "STOP: missing $P/$f"; exit 1; }; mkdir -p $(dirname $f); cp -p $P/$f $f; done
git add $FILES
git commit -F - <<'MSG'
Replication sweep: N03 pooled on its conflict-outcome term; broader-rule appendix count (90b)

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01Dof3mMRXGfgtxqVmTQMu55
MSG
git push || { echo "STOP: push failed"; exit 1; }
echo "PUSHED"
