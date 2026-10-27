###############################################################################
# run_all.R
# Master script: runs entire pipeline with renv dependency management,
# sink-based logging, and timing.
#
# Tremors in the System: Conflict Alignment in the International System
#
# Pipeline:
#   Stage 0: Clean cache (FRESH_RERUN = TRUE wipes multiverse checkpoints)
#   Stage 1: Data cleaning        code/01_data_cleaning/01-03
#   Stage 2: Main-text analysis   code/02_analysis/01-06
#   Stage 3: Appendix robustness  code/02_analysis/07-30
#   Stage 4: Figures              code/03_figures/01-03
#   Stage 5: Replication sweep    code/04_replication_sweep/run_sweep.R (RUN_STAGE5)
#
# Usage:
#   Rscript run_all.R               # full pipeline, fresh rerun
#   Rscript -e 'source("run_all.R")'
#
# To preserve checkpoints (faster partial rerun): set FRESH_RERUN <- FALSE
# inside the Stage 0 block below.
#
# ── CHANGE LOG ───────────────────────────────────────────────────────────────
# 2026-07-14  Agent-based model removed from the manuscript AND the appendix.
#             code/03_abm/ archived to archive/abm_removed_20260714/.
#             Former Stage 3 (ABM) deleted.
#             code/04_figures/ renamed to code/03_figures/.
#             02_analysis scripts renumbered sequentially (01-26); see the
#             mapping table in README.md.
#             Archived: 05_svar.R (orphan -- no appendix figure consumed it),
#             plus five patch/diagnostic scripts. See
#             archive/code_archive_20260714/.
###############################################################################

# ── 0. Setup ─────────────────────────────────────────────────────────────────

setwd(here::here())

# ── 0a. renv: reproducible dependency management ─────────────────────────────
if (!requireNamespace("renv", quietly = TRUE)) {
  install.packages("renv")
}

if (!file.exists("renv.lock")) {
  cat("=== Initialising renv (first run) ===\n")
  renv::init(bare = TRUE)
  renv::install(c(
    "tidyverse", "fixest", "peacesciencer", "data.table", "igraph",
    "vars", "svars", "strucchange", "tseries", "sandwich",
    "CausalImpact", "tsDyn", "MSwM", "quantreg",
    "patchwork", "broom", "readxl", "here", "texreg", "scales"
  ))
  renv::snapshot()
  cat("=== renv snapshot saved ===\n")
} else {
  cat("=== Restoring renv library ===\n")
  renv::restore(prompt = FALSE)
}

# ── 0b. Logging ──────────────────────────────────────────────────────────────

LOG_DIR <- here::here("logs")
dir.create(LOG_DIR, showWarnings = FALSE, recursive = TRUE)

timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
log_file  <- file.path(LOG_DIR, paste0("run_all_", timestamp, ".log"))

log_con <- file(log_file, open = "wt")
sink(log_con, type = "output")
sink(log_con, type = "message")

cat(sprintf("=== Pipeline started: %s ===\n", Sys.time()))
cat(sprintf("R version: %s\n", R.version.string))
cat(sprintf("Working directory: %s\n", getwd()))
cat(sprintf("Log file: %s\n\n", log_file))

# ── 0c. Helpers ──────────────────────────────────────────────────────────────

run_script <- function(path, label) {
  cat(sprintf("\n%s\n", paste(rep("=", 72), collapse = "")))
  cat(sprintf("  START: %s\n", label))
  cat(sprintf("  File:  %s\n", path))
  cat(sprintf("  Time:  %s\n", Sys.time()))
  cat(sprintf("%s\n\n", paste(rep("=", 72), collapse = "")))

  t0 <- proc.time()

  result <- tryCatch({
    source(path, local = new.env(parent = globalenv()), echo = FALSE)
    "SUCCESS"
  }, error = function(e) {
    cat(sprintf("\n*** ERROR in %s ***\n%s\n\n", label, conditionMessage(e)))
    paste("FAILED:", conditionMessage(e))
  })

  elapsed <- (proc.time() - t0)["elapsed"]
  cat(sprintf("\n  END:     %s\n", label))
  cat(sprintf("  Status:  %s\n", result))
  cat(sprintf("  Elapsed: %.1f seconds\n\n", elapsed))

  list(script = label, status = result, seconds = elapsed)
}

# Convenience wrapper: gc() then run.
step <- function(dir, file, label) {
  gc()
  results[[length(results) + 1]] <<- run_script(
    here::here("code", dir, file), label
  )
}

results <- list()

# ═══════════════════════════════════════════════════════════════════════════════
# STAGE 0: CLEAN CACHE (fresh-rerun mode)
# ═══════════════════════════════════════════════════════════════════════════════
# Wipes intermediate caches so every stage genuinely reruns rather than loading
# stale artifacts. Set FRESH_RERUN to FALSE for incremental partial reruns.
FRESH_RERUN <- TRUE

# ── PARALLEL WORKERS ─────────────────────────────────────────────────────────
# 07_robustness_multiverse.R parallelises with furrr/multisession. EACH worker
# serialises its own copy of the globals, which include dyad_data (~1.9 GB), so
# peak memory is roughly  N_WORKERS x 2 GB.
#
# The script's own default is availableCores() - 14, which on a many-core box
# spawns enough workers to exhaust RAM -- this is what killed previous runs.
# Cap it here. Raise only if you are confident about headroom.
#
#   4 workers ~  8-10 GB
#   6 workers ~ 12-14 GB
#   8 workers ~ 16-18 GB
N_WORKERS <- 4L
Sys.setenv(ROBUSTNESS_WORKERS = as.character(N_WORKERS))
cat(sprintf("\nROBUSTNESS_WORKERS set to %d (peak RAM ~%d GB in stage 3)\n",
            N_WORKERS, N_WORKERS * 2L))

if (FRESH_RERUN) {
  cat("\n\n###  STAGE 0: CLEAN CACHE (FRESH_RERUN = TRUE)  ###\n\n")

  # 07_robustness_multiverse.R checkpoint dir. Without clearing this, the
  # multiverse loop loads cached per-(transform, method) batches from disk
  # and skips actual recomputation.
  ckpt_dir <- here::here("data", "processed", ".robustness_checkpoints")
  if (dir.exists(ckpt_dir)) {
    n_ck <- length(list.files(ckpt_dir, pattern = "\\.rds$"))
    unlink(list.files(ckpt_dir, full.names = TRUE))
    cat(sprintf("  Cleared %d checkpoint files in %s\n", n_ck, ckpt_dir))
  }
  Sys.setenv(ROBUSTNESS_RESET = "1")  # 07_robustness_multiverse honors this

  cat("  Analysis CSVs are overwritten by their producing scripts on rerun\n")
  cat("  (no preemptive delete).\n")
  cat("\nFresh-rerun cache cleared. Proceeding to Stage 1.\n")
}

# ═══════════════════════════════════════════════════════════════════════════════
# STAGE 1: DATA CLEANING
# ═══════════════════════════════════════════════════════════════════════════════

cat("\n\n###  STAGE 1: DATA CLEANING  ###\n\n")

step("01_data_cleaning", "01_build_network_layers.R",
     "01 Build network layers (COW, 1816-2014)")
step("01_data_cleaning", "02_compute_network_measures.R",
     "02 Compute network measures (KL, JS, churn, Leiden)")
# [2026-10-12] 03_compute_spillover_canonical.R skipped: its only consumer (17_spillover_nonfighters.R) is
# archived with H3. The script is left in place; uncomment to rebuild the spillover RDS.
# step("01_data_cleaning", "03_compute_spillover_canonical.R",
#      "03 Spillover canonical: inv_dist_to_war + 2-hop (predetermined network)")

# ═══════════════════════════════════════════════════════════════════════════════
# STAGE 2: MAIN-TEXT ANALYSIS  (02_analysis/01-06)
# ═══════════════════════════════════════════════════════════════════════════════
# Everything the manuscript's tables and figures depend on.

cat("\n\n###  STAGE 2: MAIN-TEXT ANALYSIS  ###\n\n")

step("02_analysis", "01_system_level_analysis.R",
     "01 System level (OLS, VAR, Granger, FEVD)  -> Tables 1-3")
step("02_analysis", "02_state_level_analysis.R",
     "02 State level (TWFE, panel LP)            -> Tables 1, 4; Figure 4")
step("02_analysis", "31_table1_sumstats.R",
     "31 Table 1 summary stats (all + peace dyads) -> Table 1")
step("02_analysis", "03_contamination_demo.R",
     "03 Contamination ladder + tie equation     -> Figure 8")
step("02_analysis", "04_robustness_placebos.R",
     "04 Placebo / permutation tests             -> Figure 6 (+ appendix)")
step("02_analysis", "05_robustness_alt_mechanisms.R",
     "05 Alternative mechanisms (mediation)      -> Figure 7 (+ appendix)")
# [2026-10-12] 06 moved to Stage 3, right after 07: it reads robustness_specifications.csv, which 07 writes.

# ═══════════════════════════════════════════════════════════════════════════════
# STAGE 3: APPENDIX ROBUSTNESS  (02_analysis/07-26)
# ═══════════════════════════════════════════════════════════════════════════════
# Except 06 (sensitivity figures, which needs 07's multiverse CSV), nothing in this stage feeds the main text.
# The multiverse (07) is by far the longest-running script here.

cat("\n\n###  STAGE 3: APPENDIX ROBUSTNESS  ###\n\n")

step("02_analysis", "07_robustness_multiverse.R",
     "07 Multiverse: spec curve, meta-regression, LOO")
step("02_analysis", "06_sensitivity_plots.R",
     "06 Renders sensitivity figures (needs 04, 05 and 07 outputs) -> Figures 6, 7")
step("02_analysis", "08_validate_edge_weights.R",
     "08 Edge-weight validation (S-score, AUC, FE)")
step("02_analysis", "09_mediation_power_delta.R",
     "09 Mediation by power change")
step("02_analysis", "10_robustness_year_drops.R",
     "10 Random + sliding-window year drops")
step("02_analysis", "11_robustness_war_samples.R",
     "11 War-sample restrictions (fatal, interstate-only)")
step("02_analysis", "12_robustness_hostlev.R",
     "12 Hostility-level (CoW MID HostLev) IV")
step("02_analysis", "13_robustness_functional_form.R",
     "13 Functional form / binned treatment / long difference")
step("02_analysis", "14_igo_membership_sliding.R",
     "14 IGO membership sliding-window shares")
step("02_analysis", "15_igo_density_break.R",
     "15 IGO density structural break")
step("02_analysis", "16_robustness_igo_cadence.R",
     "16 IGO measurement cadence, leave-IGO-out, 1966+ alt-spec")
# [2026-10-12] 17 (spillover to non-belligerents, H3) and 18 (churn excluding belligerent dyads, built on the
# partner-war-exposure control) are archived in code/02_analysis/_archive/. H3 and partner war exposure were
# dropped from the paper; the peace-dyad panels of Tables 1 and 3 replace the excluded-belligerent check.

# ── Contamination robustness (appendix) ──
# All three are siblings of 03_contamination_demo.R (main text, Figure 8).
# 19 is SLOW (~2-4 h). 20 and 21 are ~30-45 min each.
step("02_analysis", "19_contamination_multiverse.R",
     "19 Contamination multiverse (24 layers x weighting x transform x window x monadic)")
step("02_analysis", "20_contamination_falsification.R",
     "20 Contamination falsification: non-conflict placebo, permutation, future conflict")
step("02_analysis", "21_contamination_specification.R",
     "21 Contamination specification: sample (at-risk / all / PRD), gapped history window")

step("02_analysis", "22_threshold_regime.R",
     "22 Markov-switching + quantile regression")
step("02_analysis", "23_wild_cluster_bootstrap.R",
     "23 Wild cluster bootstrap inference")
step("02_analysis", "24_robustness_cinc_nonmobilization.R",
     "24 CINC without mobilization components / pre-post war windows (R3)")

# Restored 2026-10-01 (from code/02_analysis/_archive and archive/code_archive_20260714/orphans):
# producers of appendix exhibits that were displayed but no longer regenerated.
step("02_analysis", "25_robustness_strength_control.R",
     "25 Strength control (fig_ms_strength_control)")
step("02_analysis", "26_robustness_gdp_control.R",
     "26 GDP control (fig_ms_gdp_robustness, tableA_gdp_robustness)")
step("02_analysis", "27_robustness_nonparticipant_kl.R",
     "27 Non-participant KL (WWII) robustness")
step("02_analysis", "28_bayesian_structural_ts.R",
     "28 Bayesian structural time series (CausalImpact)")
step("02_analysis", "29_oos_prediction_period.R",
     "29 Out-of-sample prediction (period holdouts)")
step("02_analysis", "30_winner_loser_analysis.R",
     "30 Winner/loser analysis (centrality, war outcome TWFE)")

# Added 2026-10-07: full regression tables behind the main-text figures (appendix section
# "Full Regression Tables for Main-Text Figures"). Was a standalone script, so its tables went stale.
step("02_analysis", "export_figure_tables_appendix.R",
     "Full regression tables for main-text figures (tableA_fig4/6/7/8_*.tex)")

# NOTE: former 15_igo_membership_status.R was ARCHIVED 2026-07-14 -- the file on
# disk was truncated (unclosed labs() at line 618) and the robustness it tested
# changed nothing (alternative IGO membership codings). Scripts 16-24 shifted
# down by one to close the gap.
#
# NOTE: 26_robustness_gdp_control.R has been FOLDED INTO
# 07_robustness_multiverse.R as the `extra_control` dimension (base / gdp /
# strength, applied to the FE_full model class). Same for the archived
# 25_robustness_strength_control.R. Neither runs standalone any more.

# ── NOT IN THE PIPELINE ──────────────────────────────────────────────────────
# archive/code_archive_20260714/orphans/05_svar.R
#   Structural VAR. Produced fig_svar_fevd.png, which no appendix section
#   included. The appendix's "Structural VAR" subsection actually displays
#   fig_fevd_horizon.png + fig_granger_forward_vs_reverse.png, both from
#   07_robustness_multiverse.R.
#
# code/02_analysis/_archive/25_robustness_strength_control.R
#   log(strength_t1) control. Still produces fig_ms_strength_control.png,
#   which the appendix DOES display -- the PNG is committed, but the figure
#   is no longer reproducible from the pipeline. Re-enable if a reviewer
#   asks for it.
#
# archive/abm_removed_20260714/code/03_abm/
#   Agent-based model + ABC. Removed from the manuscript and the appendix.

# ═══════════════════════════════════════════════════════════════════════════════
# STAGE 4: FIGURES  (03_figures/01-02)
# ═══════════════════════════════════════════════════════════════════════════════
# Main-text figures 6, 7, and 8 are written directly by their producing
# scripts in Stage 2. Only these two are standalone.

cat("\n\n###  STAGE 4: FIGURES  ###\n\n")

step("03_figures", "01_multiverse_summary_plots.R",
     "01 Multiverse summary plots (appendix)")
step("03_figures", "02_case_study_map.R",
     "02 Ottoman-Byzantine four-panel map        -> Figure 5")
step("03_figures", "03_mediation_table.R",
     "03 Mediation LaTeX table (appendix companion to Figure 7)")

# ═══════════════════════════════════════════════════════════════════════════════
# STAGE 5: REPLICATION SWEEP  (04_replication_sweep)
# ═══════════════════════════════════════════════════════════════════════════════
# Re-estimates published headline models with and without conflict-history
# controls. Needs the authors' replication data in data/raw/replication_sweep/
# (code/04_replication_sweep/01_fetch_data.R). Writes est_<id>.rds keyed by study
# ID only, pooled statistics, and the contamination figures; see that folder's
# README. Slow network models load cached estimates unless RERUN_SLOW=TRUE.
RUN_STAGE5 <- as.logical(Sys.getenv("RUN_STAGE5", "FALSE"))   # RUN_STAGE5=FALSE skips the sweep and keeps its existing results
stopifnot(!is.na(RUN_STAGE5))

if (RUN_STAGE5) {
  cat("\n\n###  STAGE 5: REPLICATION SWEEP  ###\n\n")
  step("04_replication_sweep", "run_sweep.R",
       "Replication sweep -> results/main/fig_contamination_ties + appendix")
}

# ═══════════════════════════════════════════════════════════════════════════════
# SUMMARY
# ═══════════════════════════════════════════════════════════════════════════════

cat("\n\n")
cat(paste(rep("=", 72), collapse = ""))
cat("\n  PIPELINE SUMMARY\n")
cat(paste(rep("=", 72), collapse = ""))
cat("\n\n")

summary_df <- do.call(rbind, lapply(results, function(r) {
  data.frame(
    script  = r$script,
    status  = r$status,
    seconds = round(r$seconds, 1),
    stringsAsFactors = FALSE
  )
}))

print(summary_df, row.names = FALSE)

n_ok       <- sum(summary_df$status == "SUCCESS")
n_total    <- nrow(summary_df)
total_time <- sum(summary_df$seconds)

cat(sprintf("\n  Completed: %d / %d scripts\n", n_ok, n_total))
cat(sprintf("  Total time: %.1f seconds (%.1f minutes)\n",
            total_time, total_time / 60))
cat(sprintf("  Finished: %s\n", Sys.time()))

if (n_ok < n_total) {
  cat("\n  *** FAILURES ***\n")
  failed <- summary_df[summary_df$status != "SUCCESS", ]
  for (i in seq_len(nrow(failed))) {
    cat(sprintf("    - %s: %s\n", failed$script[i], failed$status[i]))
  }
}

cat(sprintf("\n  Full log: %s\n", log_file))

# ── Close log ─────────────────────────────────────────────────────────────────

sink(type = "message")
sink(type = "output")
close(log_con)

# Steps run to the end so one failure does not waste a long run, but the pipeline must not look successful.
if (n_ok < n_total) stop(n_total - n_ok, " of ", n_total, " pipeline steps failed; see ", log_file)

cat(sprintf("\nPipeline complete. %d/%d succeeded. Log: %s\n",
            n_ok, n_total, log_file))
