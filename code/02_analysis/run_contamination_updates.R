###############################################################################
# run_contamination_updates.R
#
# Re-source the updated contamination scripts (+ the belligerent-drop
# robustness) in one batch. Each script is self-contained -- it reads from
# data/processed and writes to results/ -- so they are independent and can be
# run in any order. Fast scripts are run first so an error surfaces early; the
# two slow ones (19 multiverse, 20 permutation) run last.
#
# MEMORY: each script is sourced in its OWN environment and gc() is called
# between them, so objects from one script (e.g. the ~2 GB panel in 19) do not
# accumulate across the batch. The multiverse is additionally forced onto a
# SINGLE worker (CONTAM_WORKERS = "1"), read at the top of script 19, so it
# cannot be OOM-killed.
#
# A failing script does NOT stop the batch: the error is caught, reported, and
# the run continues. A PASS/FAIL summary is printed at the end.
#
# Usage:
#   Rscript code/02_analysis/run_contamination_updates.R
#   # or, from the project root in an R session:
#   source("code/02_analysis/run_contamination_updates.R")
###############################################################################

library(here)
set.seed(123)

ANALYSIS_DIR <- here::here("code", "02_analysis")

# Force the multiverse (script 19) onto a single worker. Must be set BEFORE
# sourcing 19, which reads Sys.getenv("CONTAM_WORKERS", "4") at the top.
Sys.setenv(CONTAM_WORKERS = "1")

# Run order: fast first, the two long-running scripts last.
SCRIPTS <- c(
  "03_contamination_demo.R",         # main demo / figure
  "21_contamination_specification.R",# sample/window checks + no-controls ladder
  "18_robustness_excl_belligerent.R",# belligerent-drop robustness
  "19_contamination_multiverse.R",   # SLOW; 1 worker (see above)
  "20_contamination_falsification.R" # SLOW; permutation null (~80 min)
)

results <- data.frame(script = character(), status = character(),
                      minutes = numeric(), stringsAsFactors = FALSE)

batch_t0 <- proc.time()

for (s in SCRIPTS) {
  path <- file.path(ANALYSIS_DIR, s)

  cat("\n\n############################################################\n")
  cat(sprintf("# RUNNING: %s\n", s))
  cat(sprintf("# start:   %s\n", format(Sys.time())))
  cat("############################################################\n")

  if (!file.exists(path)) {
    cat(sprintf("\n!!! NOT FOUND: %s\n", path))
    results <- rbind(results, data.frame(script = s, status = "NOT FOUND",
                                         minutes = NA_real_))
    next
  }

  t0 <- proc.time()
  ok <- tryCatch({
    # Fresh environment per script so memory is released between runs.
    source(path, local = new.env(parent = globalenv()), echo = FALSE)
    TRUE
  }, error = function(e) {
    cat(sprintf("\n!!! ERROR in %s:\n%s\n", s, conditionMessage(e)))
    FALSE
  })
  mins <- unname((proc.time() - t0)["elapsed"] / 60)

  results <- rbind(results, data.frame(
    script = s, status = if (ok) "OK" else "FAILED",
    minutes = round(mins, 1)))

  gc(verbose = FALSE)
  cat(sprintf("\n# %-34s -> %s  (%.1f min)\n", s, if (ok) "OK" else "FAILED", mins))
}

total_min <- unname((proc.time() - batch_t0)["elapsed"] / 60)

cat("\n\n=====================================================================\n")
cat("BATCH SUMMARY\n")
cat("=====================================================================\n")
print(results, row.names = FALSE)
cat(sprintf("\nTotal elapsed: %.1f min\n", total_min))

failed <- results$script[results$status != "OK"]
if (length(failed)) {
  cat(sprintf("\n%d script(s) did NOT complete: %s\n",
              length(failed), paste(failed, collapse = ", ")))
} else {
  cat("\nAll scripts completed successfully.\n")
}
