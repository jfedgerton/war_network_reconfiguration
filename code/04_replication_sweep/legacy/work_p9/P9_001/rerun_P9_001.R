# P9_001 Rahman 2026 JCR "Institutionalized Defense Cooperation and Regime Security" (10.1177/00220027261432196)
# Data + do-file: Harvard Dataverse doi:10.7910/DVN/8EOIJY (JCR_replication_data.dta, tables_main_text_do.do).
# HEADLINE: Table 1 Model 2 (do-file "Table 1 - Mixed effects logit models", m2):
#   xtmelogit any_coup_dum dca_dum defense prev_acoup polity_iv lnGDPpc || ccode: , intpoints(10)
#   MONADIC country-year 1980-2010, full sample. Published (article Table 1): dca_dum -0.734 (0.214), N = 4,524.
# Non-headline row: Table 1 Model 1 (same model with dcacount): published -0.550 (0.167), N = 4,524.
# KEY TERM (a priori, abstract): countries with "one or more DCAs exhibit significantly lower odds of coup attempts" -> dca_dum.
# Estimator: random-intercept logit by ccode, 10-point adaptive Gauss-Hermite quadrature -> lme4::glmer(..., nAGQ = 10), all six steps;
#   history terms enter the fixed part. No randomness (seed 123 set in common_b3.R).
CODE <- Sys.getenv("SWEEP_CODE"); stopifnot(nzchar(CODE))
source(file.path(CODE, "work_p8/common_b3.R"))
suppressPackageStartupMessages(library(lme4))
WD <- file.path(ROOT, "work_p9/P9_001"); W <- file.path(OUT, "work_p9/P9_001"); dir.create(W, recursive = TRUE, showWarnings = FALSE); sink(file.path(W, "log_P9_001.txt"), split = TRUE)
d <- as.data.frame(haven::zap_labels(haven::read_dta(file.path(WD, "JCR_replication_data.dta"))))
stopifnot(!anyDuplicated(d[, c("ccode", "year")]), all(d$year >= 1980 & d$year <= 2010))
ctrl <- c("defense", "prev_acoup", "polity_iv", "lnGDPpc")
fit_me <- function(term) function(rhs, data) {
  f <- as.formula(paste("any_coup_dum ~", paste(rhs, collapse = " + "), "+ (1 | ccode)"))
  m <- glmer(f, data = data, family = binomial("logit"), nAGQ = 10,
             control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)))
  ct <- summary(m)$coefficients
  cat("\n---- glmer", deparse(f, width.cutoff = 500), "| nAGQ = 10 | N =", nobs(m), "| groups =", ngrps(m), "| logLik =", round(as.numeric(logLik(m)), 3), "\n")
  print(round(ct, 6)); print(VarCorr(m))
  msg <- m@optinfo$conv$lme4$messages; if (length(msg)) cat("CONVERGENCE MSG:", msg, "\n")
  stopifnot(term %in% rownames(ct))
  out <- setNames(list(c(b = ct[term, 1], se = ct[term, 2], p = ct[term, 4])), term)
  attr(out, "n") <- nobs(m); out
}
pub <- list(dca_dum = c(b = -0.734, se = 0.214, n = 4524), dcacount = c(b = -0.550, se = 0.167, n = 4524))
R <- list()
for (term in c("dca_dum", "dcacount")) {
  rhs0 <- c(term, ctrl); vars <- c("any_coup_dum", rhs0, "ccode", "year"); fitter <- fit_me(term)
  full <- d[complete.cases(d[, vars]), ]
  cat("\n==== term", term, "complete-case N:", nrow(full), "\n")
  r_full <- fitter(rhs0, full); n_orig <- attr(r_full, "n")
  cat("REPRO vs article Table 1:", term, "b", r_full[[term]][["b"]], "(pub", pub[[term]][["b"]], ") se", r_full[[term]][["se"]],
      "(pub", pub[[term]][["se"]], ") N", n_orig, "(pub", pub[[term]][["n"]], ")\n")
  stopifnot(n_orig == pub[[term]][["n"]], abs(r_full[[term]][["b"]] - pub[[term]][["b"]]) / abs(pub[[term]][["b"]]) <= 0.01)
  ds <- common_sample(d, vars, "ccode", NULL, "year", monadic = TRUE)
  M <- run_six(rhs0, ds, fitter, monadic = TRUE)
  R[[term]] <- list(M = M, r_full = r_full, n_orig = n_orig, n_used = nrow(ds))
}
saveRDS(R, file.path(W, "M_P9_001.rds"))
# ---- result-CSV write step (write_p8 writes to results_p8; the file is moved to results_p9 unchanged) ----
meta <- list(pid = "P9_001", cite = "Rahman 2026", journal = "JCR", year = 2026, doi = "10.1177/00220027261432196", family = "TIE")
ev <- "abstract: 'one or more DCAs exhibit significantly lower odds of coup attempts'"
art <- "mixed-effects logit (xtmelogit), random intercept by country, 10 integration points; OIM SE"
est <- "lme4::glmer random intercept (1 | ccode), adaptive Gauss-Hermite nAGQ = 10 (= Stata intpoints(10)), all six steps"
ahc <- "Y (past coup dummy prev_acoup)"
mv <- function(term) {
  src <- file.path(OUT, "results_p8", paste0("P9_001__", term, ".csv")); dst <- file.path(OUT, "results_p9", paste0("P9_001__", term, ".csv"))
  dir.create(dirname(dst), recursive = TRUE, showWarnings = FALSE); stopifnot(file.exists(src), file.rename(src, dst))
}
nt <- function(term) sprintf(paste0("MONADIC country-year 1980-2010 (history = own-state MID onsets w15/w26; mD NA, mDM = +own w15, m26 = +own w26; ",
  "bC = own CINC, bD = own polity>=6). Authors' Dataverse data + do-file; published values from article Table 1. Full-sample repro b=%.4f se=%.4f N=%d. ",
  "Common 1822-2011 sample N=%d."), R[[term]]$r_full[[term]][["b"]], R[[term]]$r_full[[term]][["se"]], R[[term]]$n_orig, R[[term]]$n_used)
write_p8(meta, "dca_dum", ev, "Table 1 Model 2 (mixed-effects logit, coup attempt, full sample 1980-2010)", art, est,
  -0.734, 0.214, 4524, TRUE, FALSE, R$dca_dum$M, ahc, nt("dca_dum")); mv("dca_dum")
write_p8(meta, "dcacount", ev, "Table 1 Model 1 (mixed-effects logit, DCA count, full sample 1980-2010)", art, est,
  -0.550, 0.167, 4524, TRUE, FALSE, R$dcacount$M, ahc, nt("dcacount")); mv("dcacount")
sink()
