# run_one.R -- run one study's ported script(s), extract ID-keyed estimates, check against expected values.
# Usage: Rscript code/04_replication_sweep/run_one.R <study_id> [--cached]
#   --cached  do not estimate: take results/replication_sweep/cache_slow/est_<id>.rds (slow studies) and run the same check
# Writes (only when every row passes the check):
#         results/replication_sweep/est_<id>.rds   (numbers keyed by row_id/id/key_term; no citations)
# Always: stdout = the check lines, ending with "DONE <id> <n> rows PASS" or "FAILED <id>" (run_sweep.R saves it as
#         logs/replication_sweep/<id>_check.txt); script stdout/stderr go to logs/replication_sweep/<id>.log.
# Intermediate result CSVs (which carry citations) live in results/replication_sweep/_build/ and are deleted on exit.
# The body runs inside main() so that on.exit() (the _build cleanup) fires: at top level under Rscript it never runs.
# Errors (stop) unwind main() and run the cleanup before Rscript exits with status 1; a failed check returns 2 and quits after.
main <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  stopifnot(length(args) >= 1)
  ID <- args[1]; CACHED <- "--cached" %in% args

  PIPE <- Sys.getenv("PIPE_ROOT", unset = here::here())
  SW   <- file.path(PIPE, "code/04_replication_sweep")
  CODE <- file.path(SW, "legacy")
  DATA <- file.path(PIPE, "data/raw/replication_sweep")
  HIST <- Sys.getenv("SWEEP_HIST", unset = file.path(PIPE, "data/processed/replication_sweep/hist/history_tables.rds"))
  EST  <- file.path(PIPE, "results/replication_sweep")
  OUT  <- file.path(EST, "_build")
  LOG  <- file.path(PIPE, "logs/replication_sweep")
  for (d in c(EST, LOG)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
  stopifnot(file.exists(HIST))

  reg <- read.csv(file.path(SW, "registry.csv"), colClasses = "character")
  em  <- read.csv(file.path(SW, "extract_map.csv"), colClasses = "character")
  rows <- reg[reg$id == ID, ]
  if (nrow(rows) == 0) stop("unknown study id: ", ID)
  stopifnot(all(rows$row_id %in% em$row_id), all(em$id[match(rows$row_id, em$row_id)] == ID))
  COLS <- c("b0", "p0", "b_hist", "p_hist", "b_26", "p_26", "placebo_pct", "jointdem_pct", "n")

  if (CACHED) {
    # ---- 1'. slow study: cached estimates (built from the sweep's logged chains; see SLOW_MODELS.md) ----
    stopifnot(all(rows$slow == "Y"))
    est <- readRDS(file.path(EST, "cache_slow", paste0("est_", ID, ".rds")))
    est$row_id <- as.character(est$row_id)
    stopifnot(setequal(est$row_id, rows$row_id), all(est$id == ID), all(COLS %in% names(est)))
    est <- est[match(rows$row_id, est$row_id), c("row_id", "id", "key_term", COLS)]
  } else {
    # ---- 1. run the ported estimating script(s) in a clean _build ----
    unlink(OUT, recursive = TRUE); dir.create(OUT, recursive = TRUE)
    on.exit(unlink(OUT, recursive = TRUE), add = TRUE)            # named intermediate CSVs never persist, even on error
    logf <- file.path(LOG, paste0(ID, ".log"))
    cat("", file = logf)
    for (s in unique(rows$entry_script)) {
      f <- file.path(CODE, s); stopifnot(file.exists(f))
      cat(sprintf("[%s] running %s\n", format(Sys.time()), s), file = logf, append = TRUE)
      t0 <- Sys.time()
      rc <- system2(file.path(R.home("bin"), "Rscript"), shQuote(f), stdout = logf, stderr = logf,
                    env = c(paste0("SWEEP_CODE=", CODE), paste0("SWEEP_DATA=", DATA), paste0("SWEEP_OUT=", OUT),
                            paste0("SWEEP_HIST=", HIST),
                            "OMP_NUM_THREADS=1", "OPENBLAS_NUM_THREADS=1", "MKL_NUM_THREADS=1"))   # single-threaded BLAS: deterministic sums
      cat(sprintf("[%s] exit %s after %.1f min\n", format(Sys.time()), rc, as.numeric(difftime(Sys.time(), t0, units = "mins"))),
          file = logf, append = TRUE)
      if (!identical(as.integer(rc), 0L)) { cat("FAILED", ID, "\n"); stop("script failed for ", ID, " (", s, "); see ", logf) }
    }

    # ---- 2. extract numbers (the named result CSVs stay inside _build) ----
    num <- function(x) suppressWarnings(as.numeric(x))
    est <- do.call(rbind, lapply(rows$row_id, function(rid) {
      m <- em[em$row_id == rid, ]; stopifnot(nrow(m) == 1)
      f <- file.path(OUT, m$out_file)
      if (!file.exists(f)) stop("expected output missing for row ", rid, ": ", m$out_file)
      x <- read.csv(f, check.names = FALSE, colClasses = "character")
      if (nzchar(m$row_filter)) x <- x[with(x, eval(parse(text = m$row_filter))), , drop = FALSE]
      if (nrow(x) != 1) stop("row ", rid, ": filter matched ", nrow(x), " rows in ", m$out_file)
      kt <- rows$key_term[rows$row_id == rid]
      # the output row must be the registry's term (where the output names its term)
      # (compared after dropping punctuation and case, and one may contain the other: outputs write interactions as a:b and
      #  some registry labels add a qualifier, e.g. jally_new2 vs jally_new2_participant)
      nk <- function(z) tolower(gsub("[^[:alnum:]]", "", z))
      if ("key_term" %in% names(x) && !grepl(nk(x$key_term), nk(kt), fixed = TRUE) && !grepl(nk(kt), nk(x$key_term), fixed = TRUE))
        stop("row ", rid, ": output key_term '", x$key_term, "' does not match registry '", kt, "'")
      v <- sapply(COLS, function(k) { col <- m[[paste0("col_", k)]]; if (nzchar(col)) num(x[[col]]) else NA_real_ })
      data.frame(row_id = rid, id = ID, key_term = kt, t(v), check.names = FALSE)
    }))
  }

  # ---- 3. check against expected values ----
  # Per-row tolerances in expected_est.rds: tol (relative, b0/b_hist/b_26; default 1e-6, 1e-4 where reference values were rounded
  # copies), ptol (absolute, p values; default 1e-5), btol (absolute, benchmark shifts in percentage points; default 1e-2).
  # Rows estimated with lme4::glmer get tol 1e-3, ptol 1e-3, btol 2: adaptive quadrature results move at that scale across lme4 /
  # Matrix versions and CPU architectures (seen 2026-10-01 for X036). n: exact. Missing on both sides counts as a match; missing on one side is a failure.
  ex <- readRDS(file.path(SW, "expected_est.rds")); ex$row_id <- as.character(ex$row_id)
  ex <- ex[match(est$row_id, ex$row_id), ]
  stopifnot(!anyNA(ex$row_id), all(ex$id == ID))
  rel <- function(a, b) ifelse(is.na(a) & is.na(b), 0, abs(a - b) / pmax(abs(b), 1e-12))
  dif <- function(a, b) ifelse(is.na(a) & is.na(b), 0, abs(a - b))
  ptol <- if (is.null(ex$ptol)) rep(1e-5, nrow(ex)) else ex$ptol      # per-row p tolerance (looser for glmer rows)
  btol <- if (is.null(ex$btol)) rep(1e-2, nrow(ex)) else ex$btol      # per-row benchmark tolerance, percentage points
  chk <- cbind(b0 = rel(est$b0, ex$b0) < ex$tol, b_hist = rel(est$b_hist, ex$b_hist) < ex$tol, b_26 = rel(est$b_26, ex$b_26) < ex$tol,
               p0 = dif(est$p0, ex$p0) < ptol, p_hist = dif(est$p_hist, ex$p_hist) < ptol, p_26 = dif(est$p_26, ex$p_26) < ptol,
               placebo = dif(est$placebo_pct, ex$placebo_pct) < btol, jointdem = dif(est$jointdem_pct, ex$jointdem_bench_pct) < btol,
               n = if (is.null(ex$n)) rep(NA, nrow(est)) else dif(est$n, ex$n) < 0.5)
  chk[is.na(chk)] <- FALSE
  # SWEEP_CHECK_COLS (comma list, default all) restricts the gate; used only when rebuilding expected_est.rds after a deliberate
  # change to the history controls (then only the history-independent columns b0,p0,placebo,jointdem are checked).
  CC <- strsplit(Sys.getenv("SWEEP_CHECK_COLS", paste(colnames(chk), collapse = ",")), ",")[[1]]
  stopifnot(all(CC %in% colnames(chk)))
  if (length(CC) < ncol(chk)) cat("NOTE: checking only", paste(CC, collapse = ","), "\n")
  chk <- chk[, CC, drop = FALSE]
  ok <- apply(chk, 1, all)
  for (i in seq_len(nrow(est)))
    cat(sprintf("%s row %s %s: b0 %.6g (exp %.6g)  b_hist %.6g (exp %.6g)%s\n", ID, est$row_id[i], if (ok[i]) "PASS" else "FAIL",
                est$b0[i], ex$b0[i], est$b_hist[i], ex$b_hist[i],
                if (ok[i]) "" else paste0("  failed: ", paste(colnames(chk)[!chk[i, ]], collapse = ","))))
  if (!all(ok)) {
    saveRDS(est, file.path(EST, paste0("failed_est_", ID, ".rds")))   # for inspection only; aggregation never reads it
    cat("FAILED", ID, "\n"); return(2L)   # exit code set after main() returns, so on.exit cleanup runs first
  }
  unlink(file.path(EST, paste0("failed_est_", ID, ".rds")))
  saveRDS(est, file.path(EST, paste0("est_", ID, ".rds")))
  cat(sprintf("DONE %s %d rows PASS\n", ID, nrow(est)))
  0L
}
rc <- main()
if (!identical(rc, 0L)) quit(status = rc)
