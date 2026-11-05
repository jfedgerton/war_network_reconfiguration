# CD2011 full chain (TERGM, btergm MPLE + 1000 network-year bootstrap reps, fastglm, seed 123; ~10-11 min per model,
# 7 models, 74 min logged). As run (work_net/CD2011/run_all.sh): est_cd2011.R <m> 1000 for m0 hist m26 ldv0 ldvhist pbase
# placebo, appended to log_cd2011.txt; then parse_cd.py (ported: work_net/parse_cd.R).
# Input cd2011_nets.rds is the staged output of work_net/CD2011/build_cd2011.R (data construction; not re-run here).
source(file.path(Sys.getenv("SWEEP_CODE"), "_slow/_common.R"))
ID <- "N01"; W <- file.path(OUT, "work_net/CD2011"); LOG <- file.path(W, "log_cd2011.txt")
if (!SKIP_EST) {
  dir.create(W, recursive = TRUE, showWarnings = FALSE); if (file.exists(LOG)) file.remove(LOG)   # run_all.sh appends
  for (m in c("m0", "hist", "m26", "ldv0", "ldvhist", "pbase", "placebo"))
    run_r("work_net/CD2011/est_cd2011.R", c(m, "1000"), log = LOG, wd = W, append = TRUE)
  cat("ALLDONE\n", file = LOG, append = TRUE)
}
run_r("work_net/parse_cd.R", log = file.path(W, "log_parse_cd.txt"), wd = W)
row <- function(row_id, fn, has26) { x <- read_res(paste0("results_net/CD2011__", fn, ".csv")); stopifnot(x$key_effect == kt(row_id, ID))
  data.frame(id = ID, key_term = kt(row_id, ID), b0 = num(x$b0), p0 = num(x$p0), b_dm = num(x$b_hist), p_dm = num(x$p_hist),
             b_26 = if (has26) num(x$b_26) else NA, p_26 = NA, pct_bC = if (has26) num(x$pct_placebo) else NA, pct_bD = NA, n_used = NA) }
# kstar2 / triangle: no results CSV in the sweep; the master rows were read off log_cd2011.txt (m0, hist, m26, pbase, placebo).
# Reproduced conventions of those two hand-made rows: a bootstrap p of 0 is reported as 0.5 / R = 0.0005, and
# pct_bC now uses the sweep-wide signed convention 100 * (b_placebo - b_pbase) / |b_pbase| (not rounded); the hand-made rows
# used |b_placebo| - |b_pbase|, changed 2026-10-01 after the Codex audit.
M <- readRDS(file.path(W, "log_cd2011_parsed.rds"))
pfl <- function(p) ifelse(p == 0, 0.5 / 1000, p)
row_log <- function(row_id, term) { g <- function(m) { v <- M[[m]][[term]]; stopifnot(!is.null(v)); v }
  data.frame(id = ID, key_term = kt(row_id, ID), b0 = g("m0")$b, p0 = pfl(g("m0")$p), b_dm = g("hist")$b, p_dm = pfl(g("hist")$p),
             b_26 = g("m26")$b, p_26 = pfl(g("m26")$p), pct_bC = 100 * (g("placebo")$b - g("pbase")$b) / abs(g("pbase")$b),
             pct_bD = NA, n_used = NA) }
d <- rbind(row("259", "SEq", TRUE), row("265", "MinRegime", TRUE), row("254", "SEq_LDV", FALSE),
           row_log("267", "kstar2"), row_log("268", "triangle"))
write_slow(ID, d)
