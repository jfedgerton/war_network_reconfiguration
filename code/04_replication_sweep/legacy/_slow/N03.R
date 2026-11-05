# FRITZ2025 full chain (signed TERGM, ergm.sign 1.0 tsergm, Stepping MCMC-MLE, 2 PSOCK workers; ~89 min logged).
# As run: FRITZ2025/run_q.sh m0 hist (fritz_run.R <m> Stepping 2), then q2.sh's first step (fritz_mple.R, seconds),
# then parse_fritz.py (ported: work_net/parse_fritz.R). m26 / placebo were not estimated by MCMC-MLE (MPLE points only).
source(file.path(Sys.getenv("SWEEP_CODE"), "_slow/_common.R"))
ID <- "N03"; W <- file.path(OUT, "work_net/FRITZ2025")
if (!SKIP_EST) {
  for (m in c("m0", "hist")) run_r("work_net/FRITZ2025/fritz_run.R", c(m, "Stepping", "2"), log = file.path(W, paste0("log_", m, "_Stepping.txt")), wd = W)
  cat("QDONE\n", file = file.path(W, "qdone.txt"), append = TRUE)
  run_r("work_net/FRITZ2025/fritz_mple.R", log = file.path(W, "log_mple_points.txt"), wd = W)
}
run_r("work_net/parse_fritz.R", log = file.path(W, "log_parse_fritz.txt"), wd = W)
row <- function(row_id, fn) { x <- read_res(paste0("results_net/FRITZ2025__", fn, ".csv"))
  data.frame(id = ID, key_term = kt(row_id, ID), b0 = num(x$b0), p0 = num(x$p0), b_dm = num(x$b_hist), p_dm = num(x$p_hist),
             b_26 = num(x$b_26), p_26 = NA, pct_bC = NA, pct_bD = NA, n_used = NA) }   # b_26 = "NOT_RUN" -> NA
d <- rbind(row("245", "GWESE_pos"), row("252", "GWESF_pos"), row("260", "GWESE_neg"), row("266", "GWESF_neg"))
write_slow(ID, d)
