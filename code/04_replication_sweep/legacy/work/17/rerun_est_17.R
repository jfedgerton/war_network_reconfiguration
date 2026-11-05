# ESTAUDIT 17: Chen 2021 JOP, Table 1 Model 4. Copy of rerun_17.R (same data, spec, common sample) with the estimator fixed:
#   article: Stata relogit (Tomz-King-Zeng) midsidea ... _Iyear_1952-_Iyear_2010, cluster(ddyad)  [ReplicationCodes.do line 31; OA: "rare events logistic regressions"]
#   -> King & Zeng (2001) finite-sample bias correction applied in ALL six steps (no tau / no prior correction: full directed-dyad
#   population, not a case-control sample), variance (n/(n+k))^2 * cluster-robust V (relogit convention). Implemented directly (ea2/kz.R,
#   Zelig formula, self-tested against lm.influence reference). Year dummies as fixest FE (identical ML); bias computed with full year-dummy design.
ROOT <- Sys.getenv("SWEEP_DATA"); OUT <- Sys.getenv("SWEEP_OUT"); CODE <- Sys.getenv("SWEEP_CODE")
stopifnot(nzchar(ROOT), nzchar(OUT), nzchar(CODE))
source(file.path(CODE, "hist/add_history.R")); source(file.path(CODE, "ea2/kz.R")); source(file.path(CODE, "ea2/est_helpers.R"))
suppressPackageStartupMessages({library(haven); library(fixest)})
set.seed(123)
WO <- file.path(OUT, "work/17"); dir.create(WO, recursive = TRUE, showWarnings = FALSE)  # port: outputs
W <- file.path(ROOT, "work/17"); sink(file.path(WO, "log_est_17.txt"), split = TRUE)
kz_selftest()
TIE <- "extdeporig"
X <- c("extdeporig","dyaddep0","dpcapshare","comdpno","dyaddp","c2majpow","c1dem","c2dem","cont","dist",
       "pcyrmidsidea1","pcyrmidsidea2","pcyrmidsidea3")
HVv <- c("hist_dyad_w15","hist_mono_w15","hist_dyad_w26","hist_mono_w26","bm_cinc_w15","bm_democ_w15")
ADD <- list(m0 = character(0), mD = "hist_dyad_w15", mDM = c("hist_dyad_w15","hist_mono_w15"), m26 = c("hist_dyad_w26","hist_mono_w26"),
            bC = "bm_cinc_w15", bD = "bm_democ_w15")
SSC <- ssc(adj = FALSE, cluster.adj = TRUE)

d <- read_dta(file.path(W, "replication.dta"), col_select = any_of(c("midsidea", X, "ddyad","ccode1","ccode2","year"))) |> zap_labels() |> as.data.frame()
d <- d[complete.cases(d[, c("midsidea", X, "ddyad", "year")]), ]
n_orig <- nrow(d); cat("n_orig =", n_orig, " (rerun_17: 1378861)\n"); stopifnot(n_orig == 1378861)

relogit_kz <- function(dat, rhs, tag) {
  f <- as.formula(paste("midsidea ~", paste(rhs, collapse = "+"), "| year"))
  m <- feglm(f, family = binomial, data = dat, mem.clean = TRUE, cluster = ~ddyad, ssc = SSC, nthreads = 1)
  stopifnot(m$convStatus)
  ct <- coeftable(m); used <- obs(m); pi <- unname(fitted(m, type = "response")); n <- nobs(m)
  stopifnot(length(pi) == length(used), n == length(used))
  yrs <- sort(unique(dat$year[used])); k <- length(rhs) + length(yrs)
  Xm <- as.matrix(dat[used, rhs])
  yr <- dat$year[used]
  Xfun <- function(r) { D <- outer(yr[r], yrs, "==") * 1; colnames(D) <- paste0("y", yrs); cbind(Xm[r, , drop = FALSE], D) }
  bias <- kz_bias(Xfun, n, pi)
  b_hat <- ct[rhs, 1]; se_hat <- ct[rhs, 2]
  b_t <- b_hat - bias[rhs]; se_t <- se_hat * n / (n + k); p_t <- 2 * pnorm(-abs(b_t / se_t))
  out <- cbind(b_kz = b_t, se_kz = se_t, p_kz = p_t, b_logit = b_hat, se_logit = se_hat, p_logit = ct[rhs, 4], bias = bias[rhs])
  cat("\n====", tag, " N =", n, " k =", k, " removed obs (FE all-zero) =", nrow(dat) - n, "\n"); print(signif(out, 5))
  rm(m, Xm); invisible(gc())
  list(tab = out, n = n)
}

rep <- relogit_kz(d, X, "REPRO full sample (Table 1 M4)")
cat("repro relogit b =", rep$tab[TIE, "b_kz"], " (pub |b| .0093, SE .004)\n")

d <- add_history(d, "ccode1", "ccode2", "year")
d <- d[d$year >= 1822 & d$year <= 2011, ]
d <- d[complete.cases(d[, HVv]), ]
n_used <- nrow(d); cat("n_used =", n_used, " (rerun_17: 1252992)\n"); stopifnot(n_used == 1252992)

R <- list(); RL <- list()
for (nm in names(ADD)) {
  r <- relogit_kz(d, c(X, ADD[[nm]]), nm); stopifnot(r$n == n_used)
  R[[nm]] <- c(b = r$tab[TIE, "b_kz"], se = r$tab[TIE, "se_kz"], p = r$tab[TIE, "p_kz"])
  RL[[nm]] <- c(b = r$tab[TIE, "b_logit"], se = r$tab[TIE, "se_logit"], p = r$tab[TIE, "p_logit"])
}
print(do.call(rbind, R)); print(do.call(rbind, RL))

# ---- sensitivity hook: runs only under run_sensitivity.R (SWEEP_SENS=1), never in the main sweep ----
if (Sys.getenv("SWEEP_SENS") == "1") {
  source(file.path(CODE, "hist/sens_helper.R"))
  sens_run(id = "17", dat = d, y = "midsidea", tie = TIE, ctrl = setdiff(X, TIE), fe = "year",
           fit = function(dat, add) { r <- relogit_kz(dat, c(X, add), "sens"); c(b = unname(r$tab[TIE, "b_kz"]), se = unname(r$tab[TIE, "se_kz"]), p = unname(r$tab[TIE, "p_kz"])) })
}

pc <- function(x, b0) 100 * (x - b0) / abs(b0)
write_est("17", "Chen 2021", TIE,
  "rare-events logit (Stata relogit: King-Zeng bias correction, no tau) + year dummies, SE clustered on directed dyad (ddyad); cubic peace years",
  "code (ReplicationCodes.do Table 1 Model 4) + Online Appendix ('rare events logistic regressions'; SEs clustered by dyads)",
  "relogit: logit + year FE (fixest) + King-Zeng bias correction, cluster ddyad, V scaled (n/(n+k))^2",
  "N", "all six steps used ordinary logit (relogit approximated by logit; no King-Zeng bias correction)", "Y", R, n_used,
  sprintf(paste0("Corrected (KZ) b0=%.6g; uncorrected logit b0=%.6g (p=%.3g) [= previous results/17.csv]. Uncorrected mDM b=%.6g (p=%.3g, pct_dm %.1f). ",
    "Full-sample relogit repro b=%.6g (pub |b| .0093, sign negative per text). KZ correction is tiny at N=1.25M (bias ~%.2g). Common sample N=%d (drops %d with missing benchmarks)."),
    R$m0["b"], RL$m0["b"], RL$m0["p"], RL$mDM["b"], RL$mDM["p"], pc(RL$mDM["b"], RL$m0["b"]), rep$tab[TIE, "b_kz"], rep$tab[TIE, "bias"], n_used, n_orig - n_used))
sink()
