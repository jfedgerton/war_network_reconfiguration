#!/bin/zsh
# 2026-11-14: (1) regenerate the appendix permutation table with R and check it matches the version uploaded to Overleaf,
#             (2) push all Nov 14 replication-sweep changes (88-study pool) and the permutation-table script to the release repo.
#   zsh ~/Dropbox/ottoman_byzantine/code/push_20261114_final.sh
set -u
P=~/Dropbox/ottoman_byzantine; R=$P/github_release
cd $P || exit 1

# ---- 0. expanded-set flip: summaries, tables and figures from the saved estimates (no re-estimation; a few minutes) ----
LOG=logs/flip_20261114.log
{ for s in 90_summarise 90b_broader_rule 90c_expanded_set 91_figures 93_sensitivity_summary 94_study_flow_table; do
    Rscript code/04_replication_sweep/$s.R || { echo "STOP at $s"; exit 1; }; done
  (cd code/04_replication_sweep && Rscript 95_study_table.R) || { echo "STOP at 95"; exit 1; }
  Rscript code/03_figures/05_appendix_compact_figures.R || { echo "STOP at 05"; exit 1; }
  Rscript code/03_figures/07_table_contam_spec.R || { echo "STOP at 07"; exit 1; }
  Rscript code/03_figures/09_tables_falsification_war.R || { echo "STOP at 09"; exit 1; }
} > $LOG 2>&1 || { echo "STOP: summaries failed (see $LOG)"; tail -5 $LOG; exit 1; }
for f in tableA_sweep_subsets.tex tableA_sens_subsets.tex tableA_sens_split.tex tableA_contam_falsification.tex tableA_war_samples.tex tableA_contam_spec.tex tableA_sweep_families.tex tableA_igo_share.tex tableA_sweep_overlap.tex; do
  cmp -s results/appendix/$f tmp/flip_upload/$f && echo "OK: $f identical to the Overleaf upload" \
    || { echo "STOP: $f differs from the Overleaf upload"; diff tmp/flip_upload/$f results/appendix/$f; exit 1; }
done

# ---- 1. permutation table (reads saved CSVs; seconds, no re-permutation) ----
Rscript code/03_figures/08_table_permutation.R > logs/08_table_permutation_20261114.log 2>&1 || { echo "STOP: 08_table_permutation.R failed (see logs/08_table_permutation_20261114.log)"; exit 1; }
cmp -s results/appendix/tableA_permutation.tex tmp/tableA_permutation_uploaded_20261114.tex \
  && echo "OK: R table identical to the Overleaf upload" \
  || { echo "STOP: R table differs from the Overleaf upload:"; diff tmp/tableA_permutation_uploaded_20261114.tex results/appendix/tableA_permutation.tex; exit 1; }

# ---- 2. push ----
cd $R || exit 1
git pull --ff-only || { echo "STOP: pull failed"; exit 1; }
FILES=(code/02_analysis/04_robustness_placebos.R
       code/02_analysis/07_robustness_multiverse.R
       code/02_analysis/11_robustness_war_samples.R
       code/02_analysis/20_contamination_falsification.R
       code/03_figures/05_appendix_compact_figures.R
       code/03_figures/07_table_contam_spec.R
       code/03_figures/08_table_permutation.R
       code/03_figures/09_tables_falsification_war.R
       code/04_replication_sweep/90_summarise.R
       code/04_replication_sweep/90c_expanded_set.R
       code/04_replication_sweep/91_figures.R
       code/04_replication_sweep/93_sensitivity_summary.R
       code/04_replication_sweep/94_study_flow_table.R
       code/04_replication_sweep/95_study_table.R
       code/04_replication_sweep/chk_outcome_overlap.py
       code/04_replication_sweep/DECISIONS_20261114.md
       code/04_replication_sweep/expected_est.rds
       code/04_replication_sweep/registry.csv
       code/04_replication_sweep/study_flow.csv
       code/04_replication_sweep/study_table_labels.csv
       code/04_replication_sweep/legacy/work/57/rerun_57.R
       code/04_replication_sweep/legacy/work_p4/X040/rerun_verify_X040.R
       code/04_replication_sweep/legacy/work_p8/P8_031/rerun_P8_031.R
       code/04_replication_sweep/legacy/work_p8/P8_048/rebuild_P8_048.R
       code/04_replication_sweep/legacy/work_p9/P9_029/rebuild_P9_029.R
       code/04_replication_sweep/legacy/work_p9/P9_080/rebuild_P9_080.R
       code/run_restructure_20261114.sh code/run_restructure_summaries_20261114.sh
       code/run_P3_20261114.sh code/run_P4_20261114.sh
       code/push_20261114_final.sh)
for f in $FILES; do [[ -f $P/$f ]] || { echo "STOP: missing $P/$f"; exit 1; }; mkdir -p $(dirname $f); cp -p $P/$f $f; done
# run_all.R: the release copy differs on purpose (RUN_STAGE5 defaults to FALSE); the 08 step was added to it directly
grep -q "08_table_permutation.R" run_all.R || { echo "STOP: release run_all.R lacks the 08 step"; exit 1; }
git add $FILES run_all.R
git diff --cached --quiet && { echo "Nothing new to commit"; exit 0; }
git commit -F - <<'MSG'
Replication sweep: expanded set (122 studies; 77 tie) as headline; 88-study pooled subset; permutation table

- Outcome-overlap rule: 34 studies whose outcome codes every year of a dispute (or is conditional on an
  ongoing dispute) move to pool_status outcome_overlap, reported as the appendix expanded set
  (90c_expanded_set.R, Table 10; classification check in chk_outcome_overlap.py against COW MID 4.02)
- Dropped: Leeds & Mattes 2007, McDonald 2015, Powers 2006, Park & Colaresi 2014 (do not reproduce / outside rule)
- Script fixes: P8_031 t+1 lead, 57 control, P9_029 double lag, P8_048 developed dummy, X040 working
  correlation, P9_080 NA handling
- Headline = expanded set (pooled + outcome_overlap): 90_summarise.R expanded block (expanded_headline.csv,
  tableA_sweep_subsets.tex), 91_figures.R main tie figure with 77 studies (open circles = overlap),
  93_sensitivity_summary.R on 84 studies plus subset table (tableA_sens_subsets.tex) and a check that each
  sensitivity estimate reproduces the main sweep
- Figures: no study labels in the main-text figure; overlap group in the all-attempts figure
- New code/03_figures/08_table_permutation.R and 09_tables_falsification_war.R replace the appendix permutation,
  falsification and fatal-MID-sample figures with tables; subset table gains a "What the subset shows" column
- Full decision log: code/04_replication_sweep/DECISIONS_20261114.md

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01Dof3mmrxgfgtxqvmtqmu55
MSG
git push || { echo "STOP: push failed"; exit 1; }
echo "PUSHED"
