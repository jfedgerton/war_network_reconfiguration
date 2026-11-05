###############################################################################
# run_sweep.R -- Stage 5: replication sweep (conflict history as an omitted confounder)
#
# Re-estimates each published study's headline model with and without conflict-history
# controls, using the article's own estimator. Inputs: registry.csv (studies and their
# modeling approach) + authors' replication data in data/raw/replication_sweep/.
#
# Outputs (no file pairs study names with results):
#   data/processed/replication_sweep/hist/history_tables.rds  conflict-history lookup tables (generated)
#   results/replication_sweep/est_<id>.rds       estimates keyed by neutral study id only
#   results/replication_sweep/summary_stats.rds  pooled statistics (no study rows)
#   results/main/fig_contamination_ties.pdf, results/appendix/fig_contamination_{democracy,capitalism}.pdf
# Anyone who wants the study-by-study table must join est_<id>.rds to registry.csv themselves.
#
# Flags (environment variables)
#   RUN_APPENDIX_FAMILIES  TRUE  also run democratic- and capitalist-peace studies (appendix)
#   RERUN_SLOW             FALSE re-estimate the 12 slow studies (TERGM/SAOM/AME/latent space and one
#                                Bayesian model; ~35 h in total). FALSE checks and loads their cached estimates.
#   ONLY_IDS               ""    comma list: run a subset (debugging; pooled statistics and figures are skipped)
#   RESUME                 FALSE keep est_<id>.rds of studies whose check log ends in "DONE <id> ... PASS"
#   CHECK_MD5              TRUE  verify the md5 of every input file before estimating (FALSE: existence and size only)
###############################################################################
RUN_APPENDIX_FAMILIES <- as.logical(Sys.getenv("RUN_APPENDIX_FAMILIES", "TRUE"))
RERUN_SLOW            <- as.logical(Sys.getenv("RERUN_SLOW", "FALSE"))
ONLY_IDS              <- setdiff(strsplit(Sys.getenv("ONLY_IDS", ""), ",")[[1]], "")
RESUME                <- as.logical(Sys.getenv("RESUME", "FALSE"))
CHECK_MD5             <- as.logical(Sys.getenv("CHECK_MD5", "TRUE"))
stopifnot(!is.na(RUN_APPENDIX_FAMILIES), !is.na(RERUN_SLOW), !is.na(RESUME), !is.na(CHECK_MD5))

PIPE <- Sys.getenv("PIPE_ROOT", unset = here::here())
SW   <- file.path(PIPE, "code/04_replication_sweep")
DATA <- file.path(PIPE, "data/raw/replication_sweep")
HIST <- file.path(PIPE, "data/processed/replication_sweep/hist/history_tables.rds")
EST  <- file.path(PIPE, "results/replication_sweep")
LOG  <- file.path(PIPE, "logs/replication_sweep")
RS   <- file.path(R.home("bin"), "Rscript")
dir.create(LOG, recursive = TRUE, showWarnings = FALSE)
Sys.setenv(PIPE_ROOT = PIPE, SWEEP_HIST = HIST)
set.seed(123)

BUILD <- file.path(EST, "_build")
unlink(BUILD, recursive = TRUE)                                  # leftovers of an interrupted run
# run_all.R source()s this file, so the body runs inside a function: on.exit cleanup works and nothing calls quit().
sweep_main <- function() {
on.exit(unlink(BUILD, recursive = TRUE), add = TRUE)

reg <- read.csv(file.path(SW, "registry.csv"), colClasses = "character")
stopifnot(!anyDuplicated(reg$row_id), all(reg$family %in% c("TIE", "DEMOCRACY", "CAPITALISM")))

# ---- 0. history tables (generated data, kept out of data/raw) ----
if (!file.exists(HIST)) stopifnot(system2(RS, shQuote(file.path(SW, "00_build_history.R"))) == 0)
stopifnot(!is.null(readRDS(HIST)$exh))

# ---- 0b. peacesciencer's external data (cow_trade_ddy.rds etc.) are not shipped with the package; several ports read them ----
if (!nzchar(system.file("extdata", "cow_trade_ddy.rds", package = "peacesciencer"))) peacesciencer::download_extdata()
stopifnot(nzchar(system.file("extdata", "cow_trade_ddy.rds", package = "peacesciencer")))

# ---- 1. which studies ----
ids <- unique(reg$id[reg$family == "TIE" | RUN_APPENDIX_FAMILIES])
if (length(ONLY_IDS)) { stopifnot(all(ONLY_IDS %in% reg$id)); ids <- intersect(ids, ONLY_IDS) }
slow <- unique(reg$id[reg$slow == "Y"])
cat(sprintf("studies: %d (%d slow, %s)\n", length(ids), sum(ids %in% slow), if (RERUN_SLOW) "re-estimated" else "cached"))

# ---- 2. preflight: every input file of the selected studies is present (and matches its md5) ----
src <- read.csv(file.path(SW, "data_sources.csv"), colClasses = "character")
need <- src[vapply(strsplit(src$id, ";"), function(x) any(x %in% setdiff(ids, if (RERUN_SLOW) character(0) else slow)), TRUE), ]
f <- file.path(DATA, need$rel_path)
miss <- !file.exists(f)
bad_size <- !miss & file.size(f) != as.numeric(need$bytes)
bad_md5 <- rep(FALSE, nrow(need))
if (CHECK_MD5) bad_md5[!miss & !bad_size] <- unname(tools::md5sum(f[!miss & !bad_size])) != need$md5[!miss & !bad_size]
if (any(miss | bad_size | bad_md5)) {
  cat("Input files missing or not matching data_sources.csv:\n")
  print(data.frame(id = need$id, file = need$rel_path, problem = ifelse(miss, "missing", ifelse(bad_size, "size", "md5")))[miss | bad_size | bad_md5, ],
        row.names = FALSE)
  stop("run code/04_replication_sweep/01_fetch_data.R (and place the manual PRIO files) first")
}
cat(sprintf("preflight: %d input files present%s\n", nrow(need), if (CHECK_MD5) " with matching md5" else ""))

# ---- 3. estimate (or check cached estimates) one study at a time, each in a fresh R process ----
if (!RESUME && !length(ONLY_IDS)) unlink(list.files(EST, pattern = "^(failed_)?est_.*\\.rds$", full.names = TRUE))
status <- setNames(character(length(ids)), ids)
for (id in ids) {
  chk <- file.path(LOG, paste0(id, "_check.txt"))
  if (RESUME && file.exists(file.path(EST, paste0("est_", id, ".rds"))) && file.exists(chk) &&
      any(grepl(paste0("^DONE ", id, " .* PASS$"), readLines(chk)))) { status[id] <- "PASS"; next }
  unlink(file.path(EST, paste0("est_", id, ".rds")))
  t0 <- Sys.time()
  a <- c(shQuote(file.path(SW, "run_one.R")), id, if (id %in% slow && !RERUN_SLOW) "--cached")
  rc <- system2(RS, a, stdout = chk, stderr = chk)
  status[id] <- c("0" = "PASS", "2" = "MISMATCH")[as.character(rc)]
  if (is.na(status[id])) status[id] <- "ERROR"
  if (status[id] == "PASS" && id %in% slow && !RERUN_SLOW) status[id] <- "CACHED"
  cat(sprintf("%-8s %-8s %5.1f min\n", id, status[id], as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  unlink(BUILD, recursive = TRUE)
}
print(table(status))
if (any(!status %in% c("PASS", "CACHED"))) stop("replication sweep: ", sum(!status %in% c("PASS", "CACHED")),
                                                 " studies failed; see logs/replication_sweep/<id>_check.txt and <id>.log")

# ---- 4. pooled statistics and figures (computed in memory from est_<id>.rds + registry) ----
if (length(ONLY_IDS)) { cat("ONLY_IDS set: pooled statistics and figures skipped\n"); return(invisible(status)) }
for (s in c("90_summarise.R", "91_figures.R", "94_study_flow_table.R")) stopifnot(system2(RS, shQuote(file.path(SW, s))) == 0)
cat("Stage 5 complete\n")
invisible(status)
}
sweep_main()
