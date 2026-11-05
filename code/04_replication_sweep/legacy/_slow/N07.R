# WH2011 full chain (longitudinal social-relations probit, authors' LSNProbitMID.R; 12,000 iterations after 3,000 burn-in,
# thin 12; seed 123; ~36-44 min per model, 4 models).
# As run (work_net/queue.sh for m0/hist, queue7.sh for m26/placebo): probit_wh.R wh_<m>.csv out_<m> 12000 3000 12 (the .ok
# markers for m0/hist were touched by hand), then work_net/agg_mine.R (all tags at once; here only WH2011).
# Inputs wh_<m>.csv are the staged outputs of work_net/WH2011/build_wh.R + final_wh.R (data construction; not re-run here).
source(file.path(Sys.getenv("SWEEP_CODE"), "_slow/_common.R"))
ID <- "N07"; W <- file.path(OUT, "work_net/WH2011")
if (!SKIP_EST) for (m in c("m0", "hist", "m26", "placebo")) {
  dat <- file.path(ROOT, "work_net/WH2011", paste0("wh_", m, ".csv")); stopifnot(file.exists(dat))
  run_r("work_net/WH2011/probit_wh.R", c(dat, file.path(W, paste0("out_", m)), "12000", "3000", "12"),
        log = file.path(W, paste0("log_wh_", m, ".txt")), wd = W, env = c(OPENBLAS_NUM_THREADS = "1", OMP_NUM_THREADS = "1"), ok_marker = TRUE)
}
run_r("work_net/agg_mine.R", "WH2011", log = file.path(OUT, "work_net/log_agg_mine_WH2011.txt"), wd = file.path(OUT, "work_net"))
row <- function(row_id, key) { x <- read_res(paste0("results_net/WH2011__", key, ".csv")); stopifnot(x$key_effect == kt(row_id, ID))
  data.frame(id = ID, key_term = kt(row_id, ID), b0 = num(x$b0), p0 = num(x$p0), b_dm = num(x$b_hist), p_dm = num(x$p_hist),
             b_26 = num(x$b_26), p_26 = NA, pct_bC = num(x$pct_placebo), pct_bD = NA, n_used = NA) }
d <- rbind(row("262", "alliance_level"), row("248", "trade_log_exports_i_to_j"), row("255", "trade_log_imports_j_to_i"), row("258", "joint_IGOs"))
write_slow(ID, d)
