# P9_057 full chain (heckprob ML via work_p8/heckprob_b3.R, OIM SE; six steps m0/mD/mDM/m26/bC/bD plus the full-sample fit on
# 585k directed challenger-target dyad-years; seed 123 via work_p8/common_b3.R).
# rerun_P9_057.R writes results_p8/P9_057__ally_cost.csv and moves it to results_p9/.
# Runtime: ~11 min under run_one.R on 1 BLAS thread (7 heckprob fits; 2026-10-10) -> slow = Y (> 10 min), cached by default.
source(file.path(Sys.getenv("SWEEP_CODE"), "_slow/_common.R"))
ID <- "P9_057"; W <- file.path(OUT, "work_p9/P9_057")
if (!SKIP_EST) run_r("work_p9/P9_057/rerun_P9_057.R", log = file.path(W, "run.out"), wd = W)
x <- read_res("results_p9/P9_057__ally_cost.csv"); stopifnot(x$pid == ID, x$key_term == "ally_cost")
d <- data.frame(id = ID, key_term = kt("314", ID), b0 = num(x$b0), p0 = num(x$p0), b_dm = num(x$b_dm), p_dm = num(x$p_dm),
                b_26 = num(x$b_26), p_26 = num(x$p_26), pct_bC = num(x$pct_bC), pct_bD = num(x$pct_bD), n_used = num(x$n_used))
write_slow(ID, d)
