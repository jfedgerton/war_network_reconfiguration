# 92_rebuild_expected.R -- rebuild the history-dependent reference values in expected_est.rds.
# Use ONLY after a deliberate change to the conflict-history controls (00_build_history.R / legacy/hist/add_history.R).
# Steps:
#   1. SWEEP_CHECK_COLS=b0,p0,placebo,jointdem Rscript code/04_replication_sweep/run_sweep.R
#      (every study must reproduce its history-independent columns; est_<id>.rds is written only then)
#   2. Rscript code/04_replication_sweep/92_rebuild_expected.R
#      (copies b_hist, p_hist, b_26, p_26 from est_<id>.rds into expected_est.rds; the old file is kept with a date suffix)
#   3. Rscript code/04_replication_sweep/run_sweep.R   (full check, all columns)
# Nothing else in expected_est.rds is touched: b0, p0, the benchmark shifts, n and the tolerances stay as they were.
PIPE <- Sys.getenv("PIPE_ROOT", unset = here::here())
SW <- file.path(PIPE, "code/04_replication_sweep"); EST <- file.path(PIPE, "results/replication_sweep")
reg <- read.csv(file.path(SW, "registry.csv"), colClasses = "character")
ex  <- readRDS(file.path(SW, "expected_est.rds")); ex$row_id <- as.character(ex$row_id)
stopifnot(setequal(ex$row_id, reg$row_id), !anyDuplicated(ex$row_id))

ids <- unique(reg$id)
f <- file.path(EST, paste0("est_", ids, ".rds"))
if (!all(file.exists(f))) stop("no est file for: ", paste(ids[!file.exists(f)], collapse = ", "), " -- run step 1 first")
est <- do.call(rbind, lapply(f, function(x) { e <- readRDS(x); e$row_id <- as.character(e$row_id); e }))
stopifnot(setequal(est$row_id, ex$row_id), !anyDuplicated(est$row_id))
est <- est[match(ex$row_id, est$row_id), ]
stopifnot(all(est$id == ex$id))

# the history-independent columns must agree with the reference before anything is overwritten
rel <- function(a, b) ifelse(is.na(a) & is.na(b), 0, abs(a - b) / pmax(abs(b), 1e-12))
dif <- function(a, b) ifelse(is.na(a) & is.na(b), 0, abs(a - b))
ptol <- if (is.null(ex$ptol)) rep(1e-5, nrow(ex)) else ex$ptol
ok <- rel(est$b0, ex$b0) < ex$tol & dif(est$p0, ex$p0) < ptol
ok[is.na(ok)] <- FALSE
if (!all(ok)) stop("b0/p0 differ from the reference for rows: ", paste(ex$row_id[!ok], collapse = ", "))

HIST_COLS <- c("b_hist", "p_hist", "b_26", "p_26")
changed <- sapply(HIST_COLS, function(k) sum(dif(est[[k]], ex[[k]]) > 0, na.rm = TRUE))
na_flip <- sapply(HIST_COLS, function(k) sum(is.na(est[[k]]) != is.na(ex[[k]])))
if (any(na_flip > 0)) stop("a history column is missing on one side only: ", paste(names(na_flip)[na_flip > 0], collapse = ", "))
# rows whose history columns differ from the current reference by more than the check tolerance (the rest already agree)
beyond <- rel(est$b_hist, ex$b_hist) >= ex$tol | rel(est$b_26, ex$b_26) >= ex$tol | dif(est$p_hist, ex$p_hist) >= ptol | dif(est$p_26, ex$p_26) >= ptol
cat("studies whose history columns differ from the current reference beyond tolerance:", length(unique(ex$id[beyond])), "of", length(ids), "\n")
if (any(beyond)) cat(" ", paste(sort(unique(ex$id[beyond])), collapse = " "), "\n")
old <- file.path(SW, paste0("expected_est_before_rebuild_", format(Sys.Date(), "%Y%m%d"), ".rds"))
stopifnot(file.copy(file.path(SW, "expected_est.rds"), old, overwrite = TRUE))
for (k in HIST_COLS) ex[[k]] <- est[[k]]
saveRDS(ex, file.path(SW, "expected_est.rds"))
cat("rows:", nrow(ex), " studies:", length(ids), "\nvalues changed per column:\n"); print(changed)
cat("previous reference kept as", basename(old), "\nnow rerun run_sweep.R without SWEEP_CHECK_COLS\n")
