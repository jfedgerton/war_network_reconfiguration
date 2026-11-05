# X057 full chain (TERGM, btergm bootstrapped MPLE over 6 community-embedding files; m0/mH R = 500, m26/mP R = 100;
# seed 1912 per fit; 2 cores). As run (no driver script; manual call, log work_saom/X057/run_R500.out):
# rerun_X057.R 500 -> res_X057_R500.rds (the script default R = 2000 run was aborted; run_R2000.out).
# Logged runtime (run_R500.out): 24 fits, ~74 min total.
source(file.path(Sys.getenv("SWEEP_CODE"), "_slow/_common.R"))
ID <- "X057"; W <- file.path(OUT, "work_saom/X057")
if (!SKIP_EST) run_r("work_saom/X057/rerun_X057.R", "500", log = file.path(W, "run_R500.out"), wd = W)
# Extraction. No saved writer in the sweep: results_saom/X057__weak_comms.csv = mean over the 6 embedding fits (as the script's
# own summary), b rounded to 5 decimals, p to 4. pct_bC uses the sweep-wide signed convention 100 * (b_placebo - b0) / |b0|
# (rounded to 1 decimal); the original CSV used 100 * (|b0| - |b_placebo|) / |b0|, changed 2026-10-01 after the Codex audit.
# m26 is estimated but was not carried into the master table, so b_26 stays blank.
r <- readRDS(file.path(W, "res_X057_R500.rds"))
avg <- function(m) { x <- r[r$model == m & r$term == "edgecov.weak_comms", ]; stopifnot(nrow(x) == 6, length(unique(x$file)) == 6)
  c(b = mean(x$b), p = mean(x$p_boot)) }
a <- avg("m0"); h <- avg("mH"); pl <- avg("mP")
d <- data.frame(id = ID, key_term = kt("242", ID), b0 = round(a[["b"]], 5), p0 = round(a[["p"]], 4), b_dm = round(h[["b"]], 5),
                p_dm = round(h[["p"]], 4), b_26 = NA, p_26 = NA, pct_bC = round(100 * (pl[["b"]] - a[["b"]]) / abs(a[["b"]]), 1),
                pct_bD = NA, n_used = NA)
write_slow(ID, d)
