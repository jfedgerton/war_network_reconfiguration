# ESTAUDIT 53 Kim & Rousseau 2005 JPR -- copy of rerun_53.R with the estimator changed to the paper's two-stage probit.
# Author do-file K&R_TABLEI.do / log: Table I reports (col 1) single-equation probit AND (col 2) the simultaneous-equation
#   two-stage probit least squares (Keshk's cdsimeq: OLS depend_m <-> probit uof, Amemiya/Maddala corrected VCE). The paper's
#   argument (trade and conflict are endogenous; "half right") rests on col 2 (key audit: headline = 2SPROBIT). rerun_53.R ran only col 1.
# cdsimeq is ported line-for-line from cdsimeq.ado v2.0 (work/58/cdsimeq.ado; same port validated exactly in rerun_58.R):
#   stage 1: OLS cont_dep ~ all exog; probit dich_dep ~ all exog (V_o = OIM vcov); stage 2: probit dich_dep ~ I_cont + dich_iv;
#   corrected VCE alpha2 = dd + MA_d * dd G'V0^-1 (X'X)^-1 V0^-1 G dd.
# Author log target (full sample): I_depend_m = -.0122649 (corr. SE .0180213), N = 5476.
# History/benchmark terms enter the conflict (probit) equation, hence also the stage-1 exogenous list (as in rerun_58.R).
ROOT <- Sys.getenv("SWEEP_DATA"); OUT <- Sys.getenv("SWEEP_OUT"); CODE <- Sys.getenv("SWEEP_CODE")
stopifnot(nzchar(ROOT), nzchar(OUT), nzchar(CODE))
source(file.path(CODE, "work/rerun_common_k6.R"))
source(file.path(CODE, "work_p4/est_common.R"))
W <- file.path(ROOT, "work/53"); WO <- file.path(OUT, "work/53"); dir.create(WO, recursive = TRUE, showWarnings = FALSE); sink(file.path(WO, "log_est_53.txt"), split = TRUE)
oim_vcov <- function(X, y, b) {
  xb <- drop(X %*% b); q <- 2 * y - 1; lam <- q * dnorm(q * xb) / pnorm(q * xb); w <- lam * (lam + xb)
  V <- solve(crossprod(X * sqrt(w))); dimnames(V) <- list(colnames(X), colnames(X)); V
}
probit_nr <- function(X, y) {
  b <- rep(0, ncol(X)); b[ncol(X)] <- qnorm(mean(y))
  for (it in 1:100) {
    xb <- drop(X %*% b); q <- 2 * y - 1; lam <- q * dnorm(q * xb) / pnorm(q * xb)
    step <- drop(solve(crossprod(X * sqrt(lam * (lam + xb))), crossprod(X, lam))); b <- b + step
    if (max(abs(step)) < 1e-10) break
  }
  stopifnot(max(abs(step)) < 1e-8); names(b) <- colnames(X); b
}
cdsimeq_probit <- function(d, cont_dep, cont_iv, dich_dep, dich_iv, tie) {
  exog <- unique(c(cont_iv, dich_iv))
  d <- d[complete.cases(d[, c(cont_dep, dich_dep, exog)]), ]
  X <- cbind(as.matrix(d[, exog]), `_cons` = 1); y1 <- d[[cont_dep]]; y2 <- d[[dich_dep]]
  XtX <- crossprod(X); stopifnot(rcond(XtX) > 1e-14)
  pie_1 <- drop(solve(XtX, crossprod(X, y1))); names(pie_1) <- colnames(X); v1 <- y1 - drop(X %*% pie_1)
  sigma_11 <- sum(v1^2) / (nrow(X) - ncol(X))
  pie_2 <- probit_nr(X, y2); xb1 <- drop(X %*% pie_2); V_o <- oim_vcov(X, y2, pie_2)
  sigma_12 <- mean(y2 * v1 / dnorm(xb1))
  Z <- cbind(I_cont = drop(X %*% pie_1), as.matrix(d[, dich_iv]), `_cons` = 1)
  alpha_2 <- probit_nr(Z, y2); gamma_2 <- alpha_2[["I_cont"]]
  MA_d <- gamma_2^2 * sigma_11 - 2 * gamma_2 * sigma_12
  J2 <- matrix(0, ncol(X), ncol(Z) - 1, dimnames = list(colnames(X), colnames(Z)[-1])); for (v in colnames(Z)[-1]) J2[v, v] <- 1
  G <- cbind(pie_1, J2); aa <- solve(V_o); dd <- solve(t(G) %*% aa %*% G)
  vce <- dd + MA_d * dd %*% (t(G) %*% aa %*% solve(XtX) %*% aa %*% G) %*% dd
  se <- sqrt(diag(vce)); names(se) <- names(alpha_2); z <- alpha_2 / se
  tab <- cbind(b = alpha_2, se_corr = se, z = z, p = 2 * pnorm(-abs(z)), se_uncorr = sqrt(diag(oim_vcov(Z, y2, alpha_2))))
  cat("\n---- cdsimeq probit eq:", dich_dep, "~ I_", cont_dep, "+", paste(dich_iv, collapse = " + "), "| N =", nrow(X), "\n")
  cat("sigma_11", sigma_11, "sigma_12", sigma_12, "gamma_2", gamma_2, "MA_d", MA_d, "\n"); print(round(tab, 6))
  list(b = unname(alpha_2[tie]), se = unname(se[tie]), p = unname(2 * pnorm(-abs(z[tie]))), n = nrow(X))
}
d <- as.data.frame(zap_labels(read_dta(file.path(W, "KIM & ROUSSEAU (JPR).dta"))))
stopifnot(nrow(d) == 5652); names(d) <- sub("^_", "", names(d))
cont_dep <- "depend_m"; cont_iv <- c("dem_a","dem_b","dem_ab","allies","noncomm","dist","gdp","pop","pta","fcolon","oecd","region")
dich_dep <- "uof"; dich_iv <- c("dem_a","dem_b","dem_ab","bof","bof2","allies","satis","contig","dist","major","civil_s","iisbd","peaceyrs","spline1","spline2","spline3")
rp <- cdsimeq_probit(d, cont_dep, cont_iv, dich_dep, dich_iv, "I_cont")
cat("\nREPRO vs author log: I_depend_m", rp$b, "(log -.0122649) se", rp$se, "(log .0180213) N", rp$n, "(log 5476)\n")
stopifnot(rp$n == 5476, abs(rp$b - (-.0122649)) < 1e-4, abs(rp$se - .0180213) < 1e-4)
d <- add_history(d, "state_a", "state_b", "year")
d <- restrict_sample(d, unique(c(cont_dep, cont_iv, dich_dep, dich_iv)))
M <- list(m0 = cdsimeq_probit(d, cont_dep, cont_iv, dich_dep, dich_iv, "I_cont"),
  mD  = cdsimeq_probit(d, cont_dep, cont_iv, dich_dep, c(dich_iv, "hist_dyad_w15"), "I_cont"),
  mDM = cdsimeq_probit(d, cont_dep, cont_iv, dich_dep, c(dich_iv, "hist_dyad_w15", "hist_mono_w15"), "I_cont"),
  m26 = cdsimeq_probit(d, cont_dep, cont_iv, dich_dep, c(dich_iv, "hist_dyad_w26", "hist_mono_w26"), "I_cont"),
  bC  = cdsimeq_probit(d, cont_dep, cont_iv, dich_dep, c(dich_iv, "bm_cinc_w15"), "I_cont"),
  bD  = cdsimeq_probit(d, cont_dep, cont_iv, dich_dep, c(dich_iv, "bm_democ_w15"), "I_cont"))
stopifnot(length(unique(sapply(M, `[[`, "n"))) == 1)
# port: prior sweep result (results_key/53.csv, from rerun_53.R + write_key_mine.py) is only quoted in the notes string; not produced
# by this pipeline step, so fall back to NA when absent. No effect on estimates.
okf <- file.path(OUT, "results_key", "53.csv")
ok <- if (file.exists(okf)) read.csv(okf) else data.frame(b0 = NA_real_, p0 = NA_real_, b_dm = NA_real_, pct_dm = NA_real_)
est_write("53", "Kim & Rousseau 2005", "depend_m",
  "Table I: single-eq probit (col 1) and simultaneous-equation 2SPROBIT (Keshk cdsimeq; OLS trade eq <-> probit conflict eq, Amemiya-corrected SE) (col 2); headline = 2SPROBIT",
  "code (K&R_TABLEI.do + author log); article PDF not on file",
  "cdsimeq 2SPROBIT (R port of cdsimeq.ado), probit-eq key term = instrumented I_depend_m, corrected SE", "N",
  "previous rerun_53.R ran only the single-equation probit (col 1) in all six steps; the headline 2SPROBIT (col 2) was not estimated",
  "Y", M,
  sprintf("Fixed in rerun_est_53.R. Repro I_depend_m %.6f (se %.5f, N %d) vs author log -.0122649 (.0180213). Numbers above are 2SPROBIT (key = I_depend_m). Single-eq probit (col 1, previously run, matched): b0=%.5g p0=%.3g b_dm=%.5g pct_dm=%.1f. History/benchmarks enter the conflict equation (and stage-1 exogenous list).",
          rp$b, rp$se, rp$n, ok$b0, ok$p0, ok$b_dm, ok$pct_dm))
sink()
