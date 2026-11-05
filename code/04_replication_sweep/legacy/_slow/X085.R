# X085 full chain (slow: ~3.2 h logged on 2 cores, plus up to 4 h of capped ladder time). SAOM (RSiena 1.6.6), m0 then mH.
# As run: m0 in two parts, then mH (work_saom/driver3.sh, driver4.sh):
#   part 1  CAP_H=4 REP_MIN=40 rerun_X085.R m0            (stopped BY HAND after ladder rep 6, i.e. once n2start=800 and the
#                                                          switch to n3=5000 were reached; 1.65 h)
#   part 2  RESUME=1 N2SW=1 N3SW=1 CAP_H=3 REP_MIN=90 m0  (resumes from fit_m0_best.rds; converged in 1 rep, 55.5 min)
#   mH      CAP_H=4.5 REP_MIN=90 rerun_X085.R mH          (warm start from res_m0.rds; converged after 3 reps, 2.07 h)
# The manual stop is reproduced with a wall-clock cap on part 1: X085_PART1_CAP_H (default 1.5 h: rep 6 started at 1.29 h and
# rep 7 would have started at 1.65 h in the original run). The cap and the REP_MIN timeouts are machine-speed dependent, so the
# script asserts that part 1 ended in the same ladder state (n3 switch announced, not converged) and stops otherwise.
source(file.path(Sys.getenv("SWEEP_CODE"), "_slow/_common.R"))
ID <- "X085"; W <- file.path(OUT, "work_saom/X085")
if (!SKIP_EST) {
  p1_cap <- Sys.getenv("X085_PART1_CAP_H", "1.5")
  run_r("work_saom/X085/rerun_X085.R", "m0", log = file.path(W, "run_m0_part1.out"), wd = W, env = c(CAP_H = p1_cap, REP_MIN = "40"))
  l1 <- readLines(file.path(W, "run_m0_part1.out"))
  if (!any(grepl("Switching to n3 = 5000", l1, fixed = TRUE)) || any(grepl("Estimation converged", l1, fixed = TRUE)) ||
      !file.exists(file.path(W, "fit_m0_best.rds")))
    stop("X085 m0 part 1 did not end in the recorded state (n2/n3 switched, not converged); adjust X085_PART1_CAP_H")
  run_r("work_saom/X085/rerun_X085.R", "m0", log = file.path(W, "run_m0_part2.out"), wd = W,
        env = c(RESUME = "1", N2SW = "1", N3SW = "1", CAP_H = "3", REP_MIN = "90"))
  run_r("work_saom/X085/rerun_X085.R", "mH", log = file.path(W, "run_mH.out"), wd = W, env = c(CAP_H = "4.5", REP_MIN = "90"))
  for (f in c("run_m0_part2.out", "run_mH.out")) if (!any(grepl("Estimation converged", readLines(file.path(W, f)), fixed = TRUE)))
    warning("X085: ", f, " has no 'Estimation converged' line (ladder stopped at its cap)")
}
# Extraction. No saved writer in the sweep: results_saom/X085__joint_democracy_in_SAOM.csv was written by hand from
# res_m0.rds / res_mH.rds with b rounded to 4 decimals and p to 5 (reproduces the master values).
r0 <- readRDS(file.path(W, "res_m0.rds")); rH <- readRDS(file.path(W, "res_mH.rds"))
stopifnot(!r0$probe, !rH$probe)
pick <- function(r, eff) { k <- which(r$res$dv == "con_est" & r$res$effect == eff & r$res$type == "eval"); stopifnot(length(k) == 1); r$res[k, ] }
a <- pick(r0, "con_est: same dem_est x dem_est alter"); h <- pick(rH, "con_est: same dem_est x dem_est alter")
d <- data.frame(id = ID, key_term = kt("240", ID), b0 = round(a$b, 4), p0 = round(a$p, 5), b_dm = round(h$b, 4), p_dm = round(h$p, 5),
                b_26 = NA, p_26 = NA, pct_bC = NA, pct_bD = NA, n_used = NA)
write_slow(ID, d)
