# Rerun 103: Bak 2017 International Interactions
# Headline: Table 1 col 1, probit COW_cwinit (MID initiation vs protege), robust cluster(dyadid), directed dyad-years.
# Headline tie: relcincally (defensive-ally capability relative, weighted by proximity).
ROOT <- Sys.getenv("SWEEP_DATA"); OUT <- Sys.getenv("SWEEP_OUT"); CODE <- Sys.getenv("SWEEP_CODE")
stopifnot(nzchar(ROOT), nzchar(OUT), nzchar(CODE))
source(file.path(CODE, "work/rerun_helpers_b118.R"))
library(haven)
WO <- file.path(OUT, "work/103"); dir.create(WO, recursive = TRUE, showWarnings = FALSE)  # port: outputs
W <- file.path(ROOT, "work/103"); logf <- file.path(WO, "log_103.txt"); cat("", file = logf)
sink(logf, split = TRUE, append = TRUE)

# ---- 1. Load ----
control <- c("pchaloff","pchalneu","defensiveallydum","COW_contigdum","rival","pastconflict","COW_sscore",
             "COW_parity","jointdem","peace","peace2","peace3")
base_rhs <- c("cowcontigdum","relcincally","fpsimsrswvas", control)
dv <- "COW_cwinit"; tie <- "relcincally"
# read only needed columns (memory: machine shared)
d <- as.data.frame(read_dta(file.path(W, "data.dta"), col_select = all_of(c("ccode1","ccode2","year","dyadid", dv, base_rhs))))
n_orig_file <- nrow(d); stopifnot(n_orig_file == 585467)

# ---- 2. Reproduce ----
m_rep <- fit_glm(base_rhs, dv, d, "probit", "dyadid")
cat("\n==== REPRODUCTION (full sample) ====\n"); print(coeftable(m_rep)); cat("N =", nobs(m_rep), "\n")
stopifnot(nobs(m_rep) == 582126)
n_orig <- nobs(m_rep)
b_pub <- -0.478  # table1.txt col 1 (se 0.073), 3-decimal outreg2 rounding
b_repro <- unname(coef(m_rep)[tie]); p_repro <- coeftable(m_rep)[tie, 4]
rel <- abs(b_repro - b_pub) / abs(b_pub)
reproduced <- ifelse(rel <= 0.10 & p_repro < .05, "REPRODUCED", "NOT_REPRODUCED")
cat(sprintf("b_pub=%g b_repro=%g se=%g rel.diff=%.4f p=%.3g -> %s\n", b_pub, b_repro, coeftable(m_rep)[tie, 2], rel, p_repro, reproduced))

# ---- 3. History (directed dyad; add_history handles ordering) ----
d <- add_history(d, "ccode1", "ccode2", "year"); d$year_h <- d$year
# ---- 4. Battery ----
ds <- battery_sample(d, c(dv, base_rhs, "dyadid"))
fitter <- function(rhs, data) fit_glm(rhs, dv, data, "probit", "dyadid")
R <- run_battery(ds, base_rhs, fitter, function(m) tie_stats(m, tie), logf)
# ---- sensitivity hook: runs only under run_sensitivity.R (SWEEP_SENS=1), never in the main sweep ----
if (Sys.getenv("SWEEP_SENS") == "1") {
  source(file.path(CODE, "hist/sens_helper.R"))
  sens_run(id = "103", dat = ds, y = dv, tie = tie, ctrl = setdiff(base_rhs, tie), fe = NULL,
           fit = function(dat, add) tie_stats(fitter(c(base_rhs, add), dat), tie))
}

n_yr_drop <- sum(complete.cases(d[, c(dv, base_rhs)]) & d$year < 1822)
notes <- sprintf("Headline relcincally Tbl1 col1 (pub -0.478 rounded; rel.diff %.3f). Authors already control rival + pastconflict + cubic peace years. Dropped %d obs pre-1822, rest of drop from benchmark missingness. Directed dyads (history symmetric).",
                 rel, n_yr_drop)
write_result(103, "Bak 2017", reproduced, b_pub, b_repro, "probit, cluster dyadid", tie, n_orig, nrow(ds), R, "Y", notes)
sink()
