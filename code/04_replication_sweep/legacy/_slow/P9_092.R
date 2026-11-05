# P9_092 full chain (port of Millerd 2023, author's Dataverse archive doi:10.7910/DVN/I9ZSPR; random-effects logit via
# lme4::glmer, random intercept by dyad, nAGQ = 12; six steps m0/mD/mDM/m26/bC/bD plus the full-sample fit on 50,216 PR dyad-years).
# Common sample uses the named 0.75 guard exception in work_p8/common_b3.R (PI decision 2026-10-11; 2012-2014 outside the history window).
# As run (2026-10-11): rerun_P9_092.R writes results_p8/P9_092__gendense_lag.csv and moves it to results_p9/.
# Runtime: ~20 min per fit, ~2.5 h on 1 BLAS thread (see SLOW_MODELS.md: > 15 min -> slow = Y, cached by default).
source(file.path(Sys.getenv("SWEEP_CODE"), "_slow/_common.R"))
ID <- "P9_092"; W <- file.path(OUT, "work_p9/P9_092")
if (!SKIP_EST) run_r("work_p9/P9_092/rerun_P9_092.R", log = file.path(W, "run.out"), wd = W)
x <- read_res("results_p9/P9_092__gendense_lag.csv"); stopifnot(x$pid == ID, x$key_term == "gendense_lag")
d <- data.frame(id = ID, key_term = kt("335", ID), b0 = num(x$b0), p0 = num(x$p0), b_dm = num(x$b_dm), p_dm = num(x$p_dm),
                b_26 = num(x$b_26), p_26 = num(x$p_26), pct_bC = num(x$pct_bC), pct_bD = num(x$pct_bD), n_used = num(x$n_used))
write_slow(ID, d)
