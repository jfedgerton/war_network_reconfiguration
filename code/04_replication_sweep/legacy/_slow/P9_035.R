# P9_035 full chain (protocol-B rebuild; censored (Heckman) probit, heckprob ML via work_p8/heckprob_b3.R, OIM SE; six steps
# m0/mD/mDM/m26/bC/bD plus the full-sample fit on 690-740k directed dyad-years; seed 123 via work_p8/common_b3.R).
# As run (run_one.R, 2026-10-10): rebuild_P9_035.R writes results_p8/P9_035__out_allypc.csv and moves it to results_p9/.
# Runtime: ~25 min on 1 BLAS thread (see SLOW_MODELS convention: > 10 min -> slow = Y, cached by default).
source(file.path(Sys.getenv("SWEEP_CODE"), "_slow/_common.R"))
ID <- "P9_035"; W <- file.path(OUT, "work_p9/P9_035")
if (!SKIP_EST) run_r("work_p9/P9_035/rebuild_P9_035.R", log = file.path(W, "run.out"), wd = W)
# Extraction: the write_p8 row file (same columns run_one.R maps for non-slow studies).
x <- read_res("results_p9/P9_035__out_allypc.csv"); stopifnot(x$pid == ID, x$key_term == "out_allypc")
d <- data.frame(id = ID, key_term = kt("291", ID), b0 = num(x$b0), p0 = num(x$p0), b_dm = num(x$b_dm), p_dm = num(x$p_dm),
                b_26 = num(x$b_26), p_26 = num(x$p_26), pct_bC = num(x$pct_bC), pct_bD = num(x$pct_bD), n_used = num(x$n_used))
write_slow(ID, d)
