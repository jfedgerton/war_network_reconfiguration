# ESTAUDIT 11: Kinne 2012 JOP "Multilateral trade and militarized conflict" -- MONADIC state-year.
# Article headline (Kinne_JOP.do, Table 2 baseline):
#   quietly nbreg ftmids base if sample==1 ; local a = e(alpha)
#   xtgee ftmids base, family(nbinomial `a') link(nbinomial) corr(ar 1) vce(robust)      [tsset ccode year]
# The earlier rerun (rerun_11.R) did NOT estimate this model: it ran a SUPPLEMENTARY logit (Table 2 logit column, binary mid +
# peace years/splines) and left all six headline steps empty. This script runs the article's estimator in all six steps:
#   GEE, negative-binomial variance mu + a*mu^2 with the NB canonical link ln(mu/(mu+1/a)), AR(1) working correlation over
#   calendar years within ccode, sandwich (robust) SE; a re-estimated by NB2 ML (MASS::glm.nb, = Stata nbreg alpha) for each
#   step's covariate set, exactly as the do-file does for every model. Solver: work/gee_ar1_ea1.R (validated in
#   ea1/test_gee.R vs geepack fixed-correlation GEE (coef + robust SE identical) and ea1/test_gee2.R vs direct NB ML).
# Monadic history (as X032): own-state MID onsets t-1..t-5 (hist_mono_w15) / t-2..t-6 (w26); bC = own CINC mean t-1..t-5,
#   bD = own polity2>=6 share t-1..t-5 (add_history_monadic). No dyadic history exists -> mDM = mD.
ROOT <- Sys.getenv("SWEEP_DATA"); OUT <- Sys.getenv("SWEEP_OUT"); CODE <- Sys.getenv("SWEEP_CODE")
stopifnot(nzchar(ROOT), nzchar(OUT), nzchar(CODE))
source(file.path(CODE, "work/gee_ar1_ea1.R")); source(file.path(CODE, "work/est_write_ea1.R"))
source(file.path(CODE, "work_p4/common_p4_x8.R"))          # add_history_monadic()
suppressPackageStartupMessages({library(haven); library(MASS)})
set.seed(123)
WO <- file.path(OUT, "work/11"); dir.create(WO, recursive = TRUE, showWarnings = FALSE)  # port: outputs
W <- file.path(ROOT, "work/11"); sink(file.path(WO, "log_est_11.txt"), split = TRUE)
d <- as.data.frame(zap_labels(read_dta(file.path(W, "Kinne_JOP.dta"))))
stopifnot(nrow(d) == 7074, !anyDuplicated(d[, c("ccode", "year")]))
base <- c("cent50","lnopen","polity","lngdppc","lncinc","contigs","concap")
d <- d[complete.cases(d[, c("ftmids", base, "ccode", "year")]), ]
cat("GEE sample N =", nrow(d), " years", range(d$year), " states", length(unique(d$ccode)), "\n"); stopifnot(nrow(d) == 6800)

fit <- function(dat, add, tag) {
  rhs <- c(base, add)
  f <- reformulate(rhs, "ftmids")
  nb <- glm.nb(f, data = dat, control = glm.control(maxit = 100)); a <- 1 / nb$theta
  X <- model.matrix(f, dat); fam <- fam_nb_canon(a)
  b_start <- c(fam$linkfun(mean(dat$ftmids)), rep(0, ncol(X) - 1))
  g <- gee_ar1(X, dat$ftmids, dat$ccode, dat$year, fam, b_start)
  cat("\n====", tag, " N =", g$n, " panels =", g$G, " alpha(nbreg) =", signif(a, 5), " rho =", signif(g$rho, 4), " iter =", g$iter, "\n")
  print(signif(cbind(b = g$coef, se = g$se, z = g$coef / g$se, p = g$p), 5))
  list(r = c(b = unname(g$coef["cent50"]), se = unname(g$se["cent50"]), p = unname(g$p["cent50"])), n = g$n, alpha = a, rho = g$rho)
}
rep <- fit(d, character(0), "REPRO Table 2 baseline (full GEE sample)")

# [2026-11-13 audit] DV is initiations in t+1 (Kinne 2012 p.317, fn.19) -> anchor history on the outcome year
d$yr_out <- d$year + 1
d <- add_history_monadic(d, "ccode", "yr_out")
s <- d[d$year >= 1822 & d$year <= 2011, ]
s <- s[complete.cases(s[, c("hist_mono_w15","hist_mono_w26","bm_cinc_w15","bm_democ_w15")]) & is.finite(s$bm_cinc_w15) & is.finite(s$bm_democ_w15), ]
cat("common sample N =", nrow(s), " dropped =", nrow(d) - nrow(s), "\n"); stopifnot(nrow(s) > 0.9 * nrow(d))
ADD <- list(m0 = character(0), mD = "hist_mono_w15", mDM = "hist_mono_w15", m26 = "hist_mono_w26", bC = "bm_cinc_w15", bD = "bm_democ_w15")
F <- list(); for (nm in c("m0","mD","m26","bC","bD")) F[[nm]] <- fit(s, ADD[[nm]], nm)
F$mDM <- F$mD
R <- lapply(F, `[[`, "r"); stopifnot(all(sapply(F, `[[`, "n") == nrow(s)))
saveRDS(list(F = F, rep = rep), file.path(WO, "est_11.rds"))
write_est_ea1("11", "Kinne 2012", "cent50",
  "GEE negative binomial (fixed alpha from nbreg, NB canonical link), AR(1) working correlation within state, robust SE; state-year 1950-2000",
  "code (Kinne_JOP.do Table 2 baseline: xtgee ftmids ..., family(nbinomial `a') link(nbinomial) corr(ar 1) vce(robust)); audit_final headline; no PDF",
  c("GEE-NB AR(1) robust (work/gee_ar1_ea1.R), alpha from glm.nb", "GEE-NB AR(1) robust + own-state MID history w15",
    "= mD (monadic: no dyadic history)", "GEE-NB AR(1) robust + own-state MID history w26", "GEE-NB AR(1) robust + own CINC w15", "GEE-NB AR(1) robust + own polity>=6 w15"),
  "N", "previous rerun (rerun_11.R) did not estimate the headline GEE-NB AR(1) at all (MONADIC, all six steps NA); only a supplementary pooled logit of binary mid with peace-year splines, cluster ccode",
  "Y", R, nrow(s),
  sprintf(paste0("Full-sample repro cent50 b=%.4g (SE %.3g, p=%.3g), N=%d, alpha=%.3g, rho=%.3f; published value not available (no PDF/log) -> reproduction UNVERIFIED. ",
    "MONADIC: mDM = mD (own-state onsets only); bC/bD are own-state benchmarks. Supplementary logit (old) b=-2.891 -> -2.578 with own history. ",
    "rho estimated by lag-1 moment of Pearson residuals (Stata's exact AR estimator may differ slightly); robust SE without small-sample factor."),
    rep$r["b"], rep$r["se"], rep$r["p"], rep$n, rep$alpha, rep$rho))
sink()
