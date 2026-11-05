# WSC2007 full chain (latent space, authors' gbme binomial k = 3 directed, year 2000; 25,000 scans from a warm start at the
# last draw of the authors' archived chain; seed 123; ~30-34 min per model, 4 models).
# As run (work_net/queue7.sh): run_wsc.R <m> 25000 for m0, hist, m26, placebo (each touches log_wsc_<m>.txt.ok on success),
# then work_net/agg_mine.R (all tags at once; here only WSC2007).
source(file.path(Sys.getenv("SWEEP_CODE"), "_slow/_common.R"))
ID <- "N08"; W <- file.path(OUT, "work_net/WSC2007")
if (!SKIP_EST) for (m in c("m0", "hist", "m26", "placebo"))
  run_r("work_net/WSC2007/run_wsc.R", c(m, "25000"), log = file.path(W, paste0("log_wsc_", m, ".txt")), wd = W,
        env = c(OPENBLAS_NUM_THREADS = "1", OMP_NUM_THREADS = "1"), ok_marker = TRUE)
run_r("work_net/agg_mine.R", "WSC2007", log = file.path(OUT, "work_net/log_agg_mine_WSC2007.txt"), wd = file.path(OUT, "work_net"))
row <- function(row_id, key) { x <- read_res(paste0("results_net/WSC2007__", key, ".csv")); stopifnot(x$key_effect == kt(row_id, ID))
  data.frame(id = ID, key_term = kt(row_id, ID), b0 = num(x$b0), p0 = num(x$p0), b_dm = num(x$b_hist), p_dm = num(x$p_hist),
             b_26 = num(x$b_26), p_26 = NA, pct_bC = num(x$pct_placebo), pct_bD = NA, n_used = NA) }
d <- rbind(row("264", "joint_democracy_polity_i_x_polity_j"), row("257", "trade_import_ij"), row("247", "joint_IGO"))
write_slow(ID, d)
