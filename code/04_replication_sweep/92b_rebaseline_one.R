# 92b_rebaseline_one.R -- replace ONE study's reference values in expected_est.rds after a deliberate change to that study's
# estimator (first use: X092, calendar-spaced AR(1), audit 2026-11-08).
# Steps:
#   1. Rscript code/04_replication_sweep/run_one.R <ID>      -> fails the check by design, writes results/replication_sweep/failed_est_<ID>.rds
#   2. REBASELINE_OK=<ID> Rscript code/04_replication_sweep/92b_rebaseline_one.R <ID>
#      (prints old vs new for every row; the old file is kept as expected_est_before_<ID>_<date>.rds)
#   3. Rscript code/04_replication_sweep/run_one.R <ID>      -> must now PASS and write est_<ID>.rds
# Only that study's rows change; tolerances (tol, ptol, btol) are kept. Every other row is checked to be byte-identical.
PIPE <- Sys.getenv("PIPE_ROOT", unset = here::here())
SW <- file.path(PIPE, "code/04_replication_sweep"); EST <- file.path(PIPE, "results/replication_sweep")
ID <- commandArgs(trailingOnly = TRUE)[1]
stopifnot(!is.na(ID), nzchar(ID), identical(Sys.getenv("REBASELINE_OK"), ID))   # explicit confirmation, same id
f <- file.path(EST, paste0("failed_est_", ID, ".rds"))
if (!file.exists(f)) stop("no ", basename(f), " -- run step 1 first (run_one.R ", ID, ")")
est <- readRDS(f); est$row_id <- as.character(est$row_id)
ex <- readRDS(file.path(SW, "expected_est.rds")); ex$row_id <- as.character(ex$row_id); ex0 <- ex
stopifnot(all(est$id == ID), setequal(est$row_id, ex$row_id[ex$id == ID]), !anyDuplicated(est$row_id))
i <- match(est$row_id, ex$row_id)
MAP <- c(b0 = "b0", p0 = "p0", b_hist = "b_hist", p_hist = "p_hist", b_26 = "b_26", p_26 = "p_26",
         placebo_pct = "placebo_pct", jointdem_bench_pct = "jointdem_pct", n = "n")   # expected column = est column
stopifnot(all(names(MAP) %in% names(ex)), all(MAP %in% names(est)))
cmp <- do.call(rbind, lapply(names(MAP), function(k) data.frame(row_id = est$row_id, col = k, old = ex[[k]][i], new = est[[MAP[[k]]]])))
print(cmp, digits = 6, row.names = FALSE)
stopifnot(all(is.finite(est$b0)), all(is.finite(est$b_hist)), !anyNA(est$p0), !anyNA(est$p_hist))
for (k in names(MAP)) ex[[k]][i] <- est[[MAP[[k]]]]
stopifnot(identical(ex[-i, ], ex0[-i, ]), identical(names(ex), names(ex0)))
old <- file.path(SW, sprintf("expected_est_before_%s_%s.rds", ID, format(Sys.Date(), "%Y%m%d")))
stopifnot(file.copy(file.path(SW, "expected_est.rds"), old, overwrite = FALSE))
saveRDS(ex, file.path(SW, "expected_est.rds"))
cat("updated", length(i), "rows for", ID, "; previous reference kept as", basename(old), "\nnow rerun: Rscript code/04_replication_sweep/run_one.R", ID, "\n")
