###############################################################################
# run_sensitivity.R -- appendix sensitivity analyses for the replication sweep (separate from run_sweep.R).
#
# For every study listed in sens_studies.csv it re-runs the study's ported script with SWEEP_SENS=1, which activates the
# hook at the end of that script (legacy/hist/sens_helper.R): the tie coefficient with dyadic history only, monadic
# history only and both; the strength of history as a confounder next to the two benchmarks; a matched placebo; and,
# for studies with dial = Y, the lambda dial. Studies with scope = split (headline term is an interaction or product
# term) get the three-way split only. See sens_helper.R for definitions.
#
# Run AFTER run_sweep.R: every study's hook must reproduce that study's verified baseline and with-history estimates
# (est_<id>.rds) before its results are kept.
#
# Outputs (numbers keyed by study id only; no file pairs study names with results):
#   results/replication_sweep/sens/sens_<id>.rds
#   results/replication_sweep/sens_summary.rds, results/appendix/fig_sweep_*.{pdf,png}, results/appendix/tableA_sens_*.tex
#     (written by 93_sensitivity_summary.R)
# Flags (environment): SENS_WORKERS (2) studies run in parallel; SENS_DRAWS (5) placebo draws per specification;
#   ONLY_IDS comma list (subset; the pooled summary is skipped); RESUME (FALSE) keep studies already verified.
###############################################################################
PIPE <- Sys.getenv("PIPE_ROOT", unset = here::here())
SW   <- file.path(PIPE, "code/04_replication_sweep")
EST  <- file.path(PIPE, "results/replication_sweep")
SENS <- file.path(EST, "sens"); BUILD <- file.path(EST, "_build_sens")
LOG  <- file.path(PIPE, "logs/replication_sweep")
HIST <- file.path(PIPE, "data/processed/replication_sweep/hist/history_tables.rds")
RS   <- file.path(R.home("bin"), "Rscript")
WORKERS  <- as.integer(Sys.getenv("SENS_WORKERS", "2")); DRAWS <- as.integer(Sys.getenv("SENS_DRAWS", "5"))
ONLY_IDS <- setdiff(strsplit(Sys.getenv("ONLY_IDS", ""), ",")[[1]], "")
RESUME   <- as.logical(Sys.getenv("RESUME", "FALSE"))
stopifnot(!is.na(WORKERS), WORKERS >= 1, !is.na(DRAWS), DRAWS >= 1, !is.na(RESUME), file.exists(HIST))
for (d in c(SENS, LOG)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
set.seed(123)

sens_main <- function() {
unlink(BUILD, recursive = TRUE); on.exit(unlink(BUILD, recursive = TRUE), add = TRUE)   # intermediate CSVs carry citations
reg <- read.csv(file.path(SW, "registry.csv"), colClasses = "character")
ss  <- read.csv(Sys.getenv("SENS_STUDIES", unset = file.path(SW, "sens_studies.csv")), colClasses = "character")   # SENS_STUDIES: development only
ex  <- readRDS(file.path(SW, "expected_est.rds")); ex$row_id <- as.character(ex$row_id)
stopifnot(!anyDuplicated(ss$id), all(ss$id %in% reg$id), all(ss$dial %in% c("Y", "N")), all(ss$scope %in% c("full", "split")), all(ss$dial[ss$scope == "split"] == "N"), all(reg$slow[reg$id %in% ss$id] != "Y"))
ids <- ss$id; if (length(ONLY_IDS)) { stopifnot(all(ONLY_IDS %in% ids)); ids <- intersect(ids, ONLY_IDS) }
head_row <- function(id) { r <- reg[reg$id == id & reg$study_headline_term == "Y", ]; stopifnot(nrow(r) == 1); r }

# a study's result is kept only if its hook reproduces the verified estimates of the study's headline row
verify <- function(id) {
  f <- file.path(SENS, paste0("sens_", id, ".rds")); if (!file.exists(f)) return("NO OUTPUT")
  o <- readRDS(f); x <- ex[ex$row_id == head_row(id)$row_id, ]; stopifnot(nrow(x) == 1, o$id == id)
  rel <- function(a, b) abs(a - b) / max(abs(b), 1e-12)
  TOL <- max(x$tol, 1e-5)   # some study scripts store their estimates rounded to 6 significant digits; the hook returns them unrounded
  if (rel(o$split["none", "b"], x$b0) >= TOL) return("BASELINE MISMATCH")
  if (rel(o$split["both", "b"], x$b_hist) >= TOL) return("HISTORY MISMATCH")
  PT <- function(a, b, tol) abs(a - b) <= tol * max(abs(b), 1e-12) + 1e-12
  if (!PT(o$split["none", "p"], x$p0, max(x$ptol, 1e-4))) return("BASELINE P MISMATCH")
  if (!PT(o$split["both", "p"], x$p_hist, max(x$ptol, 1e-4))) return("HISTORY P MISMATCH")
  if (!all(is.finite(o$split))) return("NONFINITE SPLIT")
  if (!identical(isTRUE(o$split_only), ss$scope[ss$id == id] == "split")) return("SCOPE MISMATCH")
  if (isTRUE(o$split_only)) return("PASS")                     # split-only studies have no placebo to check
  if (is.finite(x$n) && o$n != x$n) return("N MISMATCH")       # full scope: the helper's sample must be the model's sample
  if (!isTRUE(o$check_ok)) return("LAMBDA1 MISMATCH")
  if (o$draws != DRAWS) return("STALE DRAWS")
  "PASS"
}

one <- function(id) {
  if (RESUME && identical(tryCatch(verify(id), error = function(e) "VERIFY ERROR"), "PASS")) return("PASS")
  unlink(file.path(SENS, paste0("sens_", id, ".rds")))
  r <- head_row(id); out <- file.path(BUILD, id); dir.create(out, recursive = TRUE, showWarnings = FALSE)
  logf <- file.path(LOG, paste0("sens_", id, ".log")); t0 <- Sys.time()
  rc <- system2(RS, shQuote(file.path(SW, "legacy", r$entry_script)), stdout = logf, stderr = logf,
                env = c(paste0("SWEEP_CODE=", file.path(SW, "legacy")), paste0("SWEEP_DATA=", file.path(PIPE, "data/raw/replication_sweep")),
                        paste0("SWEEP_OUT=", out), paste0("SWEEP_HIST=", HIST), "SWEEP_SENS=1", paste0("SWEEP_SENS_OUT=", SENS),
                        paste0("SENS_DRAWS=", DRAWS), paste0("SENS_DIAL=", as.integer(ss$dial[ss$id == id] == "Y")),
                        "OMP_NUM_THREADS=1", "OPENBLAS_NUM_THREADS=1", "MKL_NUM_THREADS=1"))
  unlink(out, recursive = TRUE)
  st <- if (!identical(as.integer(rc), 0L)) "ERROR" else tryCatch(verify(id), error = function(e) "VERIFY ERROR")
  if (st != "PASS") unlink(file.path(SENS, paste0("sens_", id, ".rds")))          # only verified results are kept
  cat(sprintf("%-8s %-18s %5.1f min\n", id, st, as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  st
}
cat(sprintf("sensitivity: %d studies, %d worker(s), %d placebo draws per specification\n", length(ids), WORKERS, DRAWS))
status <- unlist(parallel::mclapply(ids, one, mc.cores = WORKERS, mc.preschedule = FALSE)); names(status) <- ids
print(table(status))
if (any(status != "PASS")) stop("sensitivity: ", sum(status != "PASS"), " studies failed (", paste(ids[status != "PASS"], collapse = ", "),
                                "); see logs/replication_sweep/sens_<id>.log")
if (length(ONLY_IDS)) { cat("ONLY_IDS set: pooled summary skipped\n"); return(invisible(status)) }
# Shared study scripts (e.g. McDonald 2015, P9_086/P9_087) can write sens files for unregistered sibling ids. Move any file
# whose id is not in sens_studies.csv to sens_excluded/ so 93's guard does not stop the run [2026-11-02]
extra <- setdiff(sub("^sens_(.*)\\.rds$", "\\1", list.files(SENS, pattern = "^sens_.*\\.rds$")), ss$id)
if (length(extra)) {
  EXCL <- file.path(EST, "sens_excluded"); dir.create(EXCL, showWarnings = FALSE)
  for (x in extra) stopifnot(file.rename(file.path(SENS, paste0("sens_", x, ".rds")), file.path(EXCL, paste0("sens_", x, ".rds"))))
  cat("Moved sens files for ids not in sens_studies.csv to sens_excluded/:", paste(extra, collapse = ", "), "\n")
}
stopifnot(system2(RS, shQuote(file.path(SW, "93_sensitivity_summary.R"))) == 0)
cat("Sensitivity analyses complete\n"); invisible(status)
}
sens_main()
