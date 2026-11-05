# MINHAS2022_WEEKS full chain (AME, amen 1.4 pa2018_version ame_repL, bin, R = 2, authors' seed 6886, warm start from the
# authors' final MCMC state; ~30-33 min per model, 4 models). Needs SWEEP_LIB_AMEN (library holding amen 1.4 built from s7minhas/amen@pa2018_version).
# As run (work_net/queue.sh, queue6.sh, queue7.sh): run_minhas.R weeks <m> 500 2500 10 1950 for m0, hist, m26, placebo (each touches
# log_weeks_<m>.txt.ok on success), then work_net/agg_mine.R (all tags at once; here only MINHAS2022_WEEKS).
source(file.path(Sys.getenv("SWEEP_CODE"), "_slow/_common.R"))
ID <- "N05"; W <- file.path(OUT, "work_net/MINHAS2022")
if (!SKIP_EST) {
  lib <- Sys.getenv("SWEEP_LIB_AMEN"); stopifnot(nzchar(lib), file.exists(file.path(lib, "amen/DESCRIPTION")))
  for (m in c("m0", "hist", "m26", "placebo"))
    run_r("work_net/MINHAS2022/run_minhas.R", c("weeks", m, "500", "2500", "10", "1950"), log = file.path(W, paste0("log_weeks_", m, ".txt")), wd = W,
          env = c(SWEEP_LIB_AMEN = lib, OPENBLAS_NUM_THREADS = "1", OMP_NUM_THREADS = "1"), ok_marker = TRUE)
}
run_r("work_net/agg_mine.R", "MINHAS2022_WEEKS", log = file.path(OUT, "work_net/log_agg_mine_MINHAS2022_WEEKS.txt"), wd = file.path(OUT, "work_net"))
row <- function(row_id, key) { x <- read_res(paste0("results_net/MINHAS2022_WEEKS__", key, ".csv")); stopifnot(x$key_effect == kt(row_id, ID))
  data.frame(id = ID, key_term = kt(row_id, ID), b0 = num(x$b0), p0 = num(x$p0), b_dm = num(x$b_hist), p_dm = num(x$p_hist),
             b_26 = num(x$b_26), p_26 = NA, pct_bC = num(x$pct_placebo), pct_bD = NA, n_used = NA) }
d <- rbind(row("253", "machine"), row("251", "junta"), row("246", "boss"), row("250", "strongman"))
write_slow(ID, d)
