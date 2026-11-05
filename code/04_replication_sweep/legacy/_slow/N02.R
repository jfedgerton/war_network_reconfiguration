# CHEN2021 full chain (multilayer ERGM on the author's stack: ergm 3.10.1, multilayer.ergm 0.1.6, statnet.common 4.3.0,
# network 1.15 from SWEEP_LIB_CHEN; MCMLE seed 30748, burn-in 1e5, interval 2000, samplesize 4000, maxit 8, 2 parallel chains;
# m0 63 min, hist 70 min wall logged). As run: CHEN2021/run_mcmle2.sh (MCMLE m0, hist), CHEN2021/run_mcmle.sh's MPLE loop
# (author stack, m0 hist m26 placebo; seconds), CHEN2021/run_diag.sh (post-hoc DIAG m0, hist), then CHEN2021/write_chen.R.
source(file.path(Sys.getenv("SWEEP_CODE"), "_slow/_common.R"))
ID <- "N02"; W <- file.path(OUT, "work_net/CHEN2021")
if (!SKIP_EST) {
  lib <- Sys.getenv("SWEEP_LIB_CHEN"); stopifnot(nzchar(lib), file.exists(file.path(lib, "ergm/DESCRIPTION")))
  ev <- c(CHEN_OLD = "1", R_LIBS = lib, SWEEP_LIB_CHEN = lib, OPENBLAS_NUM_THREADS = "1", OMP_NUM_THREADS = "1")
  for (m in c("m0", "hist")) run_r("work_net/CHEN2021/chen_run.R", c(m, "MCMLE", "2", "2000", "4000", "8"), log = file.path(W, paste0("log_", m, "_MCMLE.txt")), wd = W, env = ev)
  for (m in c("m0", "hist", "m26", "placebo")) run_r("work_net/CHEN2021/chen_run.R", c(m, "MPLE", "1"), log = file.path(W, paste0("log_", m, "_MPLE.txt")), wd = W, env = ev)
  for (m in c("m0", "hist")) run_r("work_net/CHEN2021/chen_run.R", c(m, "DIAG", "1", "2000", "4000"), log = file.path(W, paste0("log_", m, "_DIAG.txt")), wd = W, env = ev)
}
run_r("work_net/CHEN2021/write_chen.R", log = file.path(W, "log_write_chen.txt"), wd = W)
row <- function(row_id, fn) { x <- read_res(paste0("results_net/CHEN2021__", fn, ".csv")); stopifnot(x$key_effect == kt(row_id, ID))
  data.frame(id = ID, key_term = kt(row_id, ID), b0 = num(x$b0), p0 = num(x$p0), b_dm = num(x$b_hist), p_dm = num(x$p_hist),
             b_26 = num(x$b_26), p_26 = NA, pct_bC = NA, pct_bD = NA, n_used = NA) }   # b_26 = "NOT_RUN" -> NA
d <- rbind(row("256", "ThreeWay"), row("263", "MIDxCivil_kstar"))
write_slow(ID, d)
