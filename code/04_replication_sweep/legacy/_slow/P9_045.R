# P9_045 full chain (heckprob ML via work_p8/heckprob_b3.R, cluster-robust SE on ddyadid; six steps m0/mD/mDM/m26/bC/bD plus the
# full-sample fit on 882-907k directed dyad-years; year reconstructed from row order inside the script; seed 123 via common_b3.R).
# rerun_P9_045.R writes results_p8/P9_045__Ldepkg1.csv and moves it to results_p9/.
# Runtime: ~46 min on 1 BLAS thread (7 heckprob fits, ~6.5 min each) -> slow = Y, cached by default.
source(file.path(Sys.getenv("SWEEP_CODE"), "_slow/_common.R"))
ID <- "P9_045"; W <- file.path(OUT, "work_p9/P9_045")
if (!SKIP_EST) run_r("work_p9/P9_045/rerun_P9_045.R", log = file.path(W, "run.out"), wd = W)
x <- read_res("results_p9/P9_045__Ldepkg1.csv"); stopifnot(x$pid == ID, x$key_term == "Ldepkg1")
d <- data.frame(id = ID, key_term = kt("300", ID), b0 = num(x$b0), p0 = num(x$p0), b_dm = num(x$b_dm), p_dm = num(x$p_dm),
                b_26 = num(x$b_26), p_26 = num(x$p_26), pct_bC = num(x$pct_bC), pct_bD = num(x$pct_bD), n_used = num(x$n_used))
write_slow(ID, d)
