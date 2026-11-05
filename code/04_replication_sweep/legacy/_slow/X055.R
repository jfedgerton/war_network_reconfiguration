# X055 full chain (TERGM, btergm bootstrapped MPLE; m0/mH R = 1000, m26/mP R = 100; seed 123; 2 cores).
# As run (no driver script; manual call, log work_saom/X055/run.out): rebuild_X055.R 1000 -> res_X055_R1000.rds.
# Logged runtime (run.out): m0 1129 s, mH 1273 s, m26 138 s, mP 150 s (~45 min total).
source(file.path(Sys.getenv("SWEEP_CODE"), "_slow/_common.R"))
ID <- "X055"; W <- file.path(OUT, "work_saom/X055")
if (!SKIP_EST) run_r("work_saom/X055/rebuild_X055.R", "1000", log = file.path(W, "run.out"), wd = W)
# Extraction. No saved writer in the sweep (work_saom/X055/extract.R only printed the rows); results_saom/X055__*.csv hold
# b rounded to 4 decimals and the bootstrap p as is. pct_bC uses the sweep-wide signed convention 100 * (b_placebo - b0) / |b0|
# (rounded to 1 decimal); the original CSVs used 100 * (|b0| - |b_placebo|) / |b0|, changed 2026-10-01 after the Codex audit.
# m26 is estimated but was not carried into the master table, so b_26 stays blank.
r <- readRDS(file.path(W, "res_X055_R1000.rds"))
pick <- function(m, term) { k <- which(r$model == m & r$term == term); stopifnot(length(k) == 1); r[k, ] }
row <- function(row_id, term) { a <- pick("m0", term); h <- pick("mH", term); pl <- pick("mP", term)
  stopifnot(a$R == 1000, h$R == 1000)
  data.frame(id = ID, key_term = kt(row_id, ID), b0 = round(a$b, 4), p0 = a$p_boot, b_dm = round(h$b, 4), p_dm = h$p_boot,
             b_26 = NA, p_26 = NA, pct_bC = round(100 * (pl$b - a$b) / abs(a$b), 1), pct_bD = NA, n_used = NA) }
d <- rbind(row("241", "edgecov.ca1"), row("244", "edgecov.ca2"), row("239", "edgecov.ca3"))
write_slow(ID, d)
