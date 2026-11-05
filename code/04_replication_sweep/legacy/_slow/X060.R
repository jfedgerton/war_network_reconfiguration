# X060 full chain (slow: ~11 h on 2 cores). SAOM (RSiena 1.6.6), m0 then mH.
# As run (work_saom/driver.sh, driver2.sh): rerun_X060.R m0, then rerun_X060.R mH (warm start from res_m0.rds), CAP_H=5
# in the environment (not read by rerun_X060.R). Writes SWEEP_OUT/_slow/X060.csv.
source(file.path(Sys.getenv("SWEEP_CODE"), "_slow/_common.R"))
ID <- "X060"; W <- file.path(OUT, "work_saom/X060")
if (!SKIP_EST) {
  run_r("work_saom/X060/rerun_X060.R", "m0", log = file.path(W, "run_m0.out"), wd = W, env = c(CAP_H = "5"))
  run_r("work_saom/X060/rerun_X060.R", "mH", log = file.path(W, "run_mH.out"), wd = W, env = c(CAP_H = "5"))
}
# Extraction. The sweep had no saved writer: results_saom/X060__*.csv were written by hand from res_m0.rds / res_mH.rds
# with b rounded to 4 decimals and p to 5 (this reproduces every master value except row 243 p0, typed as 0.00433 vs 0.00434).
r0 <- readRDS(file.path(W, "res_m0.rds")); rH <- readRDS(file.path(W, "res_mH.rds"))
stopifnot(!r0$probe, !rH$probe)
pick <- function(r, eff) { k <- which(r$res$dv == "threat.net" & r$res$effect == eff & r$res$type == "eval"); stopifnot(length(k) == 1); r$res[k, ] }
row <- function(row_id, eff) { a <- pick(r0, eff); h <- pick(rH, eff)
  data.frame(id = ID, key_term = kt(row_id, ID), b0 = round(a$b, 4), p0 = round(a$p, 5), b_dm = round(h$b, 4), p_dm = round(h$p, 5),
             b_26 = NA, p_26 = NA, pct_bC = NA, pct_bD = NA, n_used = NA) }
d <- rbind(row("238", "threat.net: support_mi.net alter"),
           row("243", "threat.net: support_mi.net alter x support_ipcloseness.net alter"))
write_slow(ID, d)
