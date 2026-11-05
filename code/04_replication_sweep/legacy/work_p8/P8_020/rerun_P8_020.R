# P8_020 Reiter & Tillman 2002 JOP "Public, Legislative, and Executive Constraints on the Democratic Initiation of Conflict"
# Data + log: Dataverse doi:10.7910/DVN/IEH0V2 (jop02rep.dta, jop02rep.smcl). MONADIC state-year, 37 democracies 1919-1992.
# Headline = Table 1 Model 2 (log): xtgee midindm vanpart eff_ptie majdum coaldum presdum prempres assind treatrat majpow autbor
#   if year>1918 & year<1993, fam(bin) link(probit) corr(ar1) t(year) force robust    [i = ccode]
# KEY (a priori, abstract): vanpart -- "higher political participation levels decrease the likelihood of initiating".
# Estimator all steps: GEE probit AR(1), Stata `force` (lags by within-panel order), semi-robust SE by ccode (geepack).
CODE <- Sys.getenv("SWEEP_CODE"); stopifnot(nzchar(CODE))
source(file.path(CODE, "work_p8/common_b2.R"))
suppressPackageStartupMessages({library(haven); library(geepack)})
WD <- file.path(ROOT, "work_p8/P8_020"); W <- file.path(OUT, "work_p8/P8_020"); dir.create(W, recursive = TRUE, showWarnings = FALSE); sink(file.path(W, "log_P8_020.txt"), split = TRUE)
d <- zap_labels(read_dta(file.path(WD, "jop02rep.dta"))) %>% as.data.frame()
print(names(d)); stopifnot(anyDuplicated(d[, c("ccode", "year")]) == 0)
d <- d[d$year > 1918 & d$year < 1993, ]
rhs <- c("vanpart","eff_ptie","majdum","coaldum","presdum","prempres","assind","treatrat","majpow","autbor")
terms <- c("vanpart")
# Stata xtgee AR(1) replicated with work/gee_ar1_ea1.R (moment rho on Pearson residuals, lags by within-panel order = force);
# robust SE x sqrt(G/(G-1)) (Stata small-sample cluster factor). Verified: reproduces the log exactly (vanpart -.0126103, assind -1.126466).
source(file.path(CODE, "work/gee_ar1_ea1.R"))
fam_probit <- list(name = "probit", linkfun = qnorm, linkinv = function(e) pnorm(e), mu.eta = dnorm,
                   variance = function(mu) mu * (1 - mu), valideta = function(e) all(is.finite(e)))
fit_gee <- function(r, data) {
  f <- as.formula(paste("midindm ~", paste(r, collapse = " + ")))
  dd <- data[complete.cases(data[, c(all.vars(f), "ccode", "year")]), ]; dd <- dd[order(dd$ccode, dd$year), ]
  dd$seq <- ave(dd$year, dd$ccode, FUN = seq_along)
  X <- model.matrix(f, dd); b_init <- coef(glm(f, binomial("probit"), dd))
  g <- gee_ar1(X, dd$midindm, dd$ccode, dd$seq, fam_probit, b_init)
  se <- g$se * sqrt(g$G / (g$G - 1)); p <- 2 * pnorm(-abs(g$coef / se))
  cat("\n----", deparse(f, width.cutoff = 500), "| N =", g$n, "| G =", g$G, "| rho =", round(g$rho, 4), "\n")
  print(round(cbind(b = g$coef, se = se, p = p), 6))
  out <- lapply(setNames(terms, terms), function(t) c(b = unname(g$coef[t]), se = unname(se[t]), p = unname(p[t]))); attr(out, "n") <- g$n; out
}
r0 <- fit_gee(rhs, d)   # full-sample reproduction (pub: vanpart -.0126103 se .0053703 N 1363)
ds <- common_sample(d, c("midindm", rhs, "ccode", "year"), "ccode", NA, "year", monadic = TRUE)
M <- run_six(rhs, ds, fit_gee, monadic = TRUE)
dc <- d[complete.cases(d[, c("midindm", rhs)]), ]; lost <- setdiff(unique(dc$ccode), unique(ds$ccode)); cat("ccodes lost in common sample:", lost, "\n")
v_full <- verdict_fn(r0$vanpart["b"], r0$vanpart["p"], -0.0126103, TRUE, FALSE)
meta <- list(pid = "P8_020", cite = "Reiter & Tillman 2002", journal = "JOP", year = 2002, doi = "10.1111/0022-3816.00147", family = "DEMOCRACY")
write_p8(meta, "vanpart", "\"higher political participation levels decrease the likelihood of initiating an international dispute\"",
  "Table 1 Model 2 (authors' log jop02rep.smcl)", "GEE probit, AR(1) corr, Stata force, semi-robust SE clustered on ccode (xtgee)",
  "GEE probit AR(1) (Stata moment rho, force/order lags; work/gee_ar1_ea1.R), robust SE by ccode x sqrt(G/(G-1)) -- all steps",
  -0.0126103, 0.0053703, 1363, TRUE, FALSE, M, "N", verdict = v_full, notes =
  sprintf("MONADIC state-year (democracies only): mD not estimable (b_d NA); mDM = + own MID-onset history t-1..t-5; m26 = + t-2..t-6; bC/bD = own lagged CINC / own lagged democracy share. Full-sample repro b=%.5f (N=%d). Published = authors\' log Table 1 Model 2; estimator reproduces the log exactly on the full sample. Verdict from full-sample reproduction; common-sample m0 (N 1303, 2 panels lost: ccodes %s, missing lagged CINC/polity benchmarks) has p=%.3f.",
          r0$vanpart["b"], attr(r0, "n"), paste(lost, collapse = "/"), M$m0$vanpart["p"]))
sink()
