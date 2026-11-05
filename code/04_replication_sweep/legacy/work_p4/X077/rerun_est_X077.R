# ESTAUDIT X077 Hwang 2010 II, Model 1.1: Stata "probit cwmid ... , robust" (hwang.do line 16). Copy of rerun_X077.R (same data/spec/common sample);
# SEs recomputed Stata-style (observed-Hessian bread, n/(n-1)) instead of fixest's expected-information "hetero" SE, in all six steps.
ROOT <- Sys.getenv("SWEEP_DATA"); OUT <- Sys.getenv("SWEEP_OUT"); CODE <- Sys.getenv("SWEEP_CODE")
stopifnot(nzchar(ROOT), nzchar(OUT), nzchar(CODE))
source(file.path(CODE, "hist/add_history.R")); source(file.path(CODE, "ea2/est_helpers.R")); source(file.path(CODE, "ea2/probit_oim.R")); set.seed(123)
WO <- file.path(OUT, "work_p4/X077"); dir.create(WO, recursive = TRUE, showWarnings = FALSE)  # port: outputs
W <- file.path(ROOT, "work_p4/X077"); sink(file.path(WO, "log_est_X077.txt"), split = TRUE)
d <- as.data.frame(haven::zap_labels(haven::read_dta(file.path(W, "hwang.tab"))))
names(d) <- sub("^_spline", "spline", names(d))
f0 <- cwmid ~ lnpower_1 + I_1 + Ipow + democracy + contig + peaceyrs + spline1 + spline2 + spline3
K <- c("Ipow", "democracy")
full <- probit_stata(f0, d, K, "robust"); cat("FULL (pub Ipow -0.129, N 502,432): N =", attr(full, "n"), "\n"); print(full)
d <- d[complete.cases(d[, all.vars(f0)]), ]; d <- add_history(d, "ccode1", "ccode2", "year"); d <- d[d$year >= 1822 & d$year <= 2011, ]
HV <- c("hist_dyad_w15","hist_mono_w15","hist_dyad_w26","hist_mono_w26","bm_cinc_w15","bm_democ_w15")
d <- d[complete.cases(d[, HV]) & is.finite(d$bm_cinc_w15) & is.finite(d$bm_democ_w15), ]
cat("common N:", nrow(d), "(results: 477029)\n"); stopifnot(nrow(d) == 477029)
ADD <- list(m0 = NULL, mD = "hist_dyad_w15", mDM = c("hist_dyad_w15","hist_mono_w15"), m26 = c("hist_dyad_w26","hist_mono_w26"), bC = "bm_cinc_w15", bD = "bm_democ_w15")
R <- lapply(ADD, function(a) { f <- if (is.null(a)) f0 else update(f0, as.formula(paste(". ~ . +", paste(a, collapse = "+")))); probit_stata(f, d, K, "robust") })
for (k in K) {
  tab <- t(sapply(R, function(r) r[, k])); print(k); print(tab); ratio <- tab[, "se"] / tab[, "se_fixest"]; cat("SE ratio:", round(ratio, 4), "\n")
  if (max(abs(ratio - 1)) > 0.02)
    write_est("X077", "Hwang 2010 II", k, "probit, robust SE (Stata: observed-Hessian bread, n/(n-1)), peace years + 3 splines (Table 1 Model 1.1)",
      "pdf (Hwang2010II) + code (hwang.do line 16: probit ..., robust)", "probit (fixest point estimates) with Stata-equivalent robust SE (observed-Hessian bread)", "N",
      sprintf("SE computation: fixest expected-information bread instead of Stata observed Hessian (SE ratio %.3f-%.3f)", min(ratio), max(ratio)), "Y",
      lapply(R, function(r) r[c("b", "se", "p"), k]), nrow(d), paste0(if (k == "democracy") "CONTROL row. " else "", sprintf("Previous fixest-SE p0 %.3g, p_dm %.3g.", tab["m0", "p_fixest"], tab["mDM", "p_fixest"])))
  else cat(k, ": SEs agree within 2% -> MATCH row kept\n")
}
sink()
