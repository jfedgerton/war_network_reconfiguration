# P8_036 McCallister 2016 II "Beyond Dyads: Regional Democratic Strength's Influence on Dyadic Conflict" (10.1080/03050629.2015.1069292)
# Data + do-file: Harvard Dataverse doi:10.7910/DVN/KT4XKZ (regdem_mccallister-1.tab original .dta; regdem_mccallister_II.do).
# HEADLINE: Table 1 Model 5 (full model, do-file "TABLE 1: Model 5"):
#   xtmelogit MID land_contiguity water_contiguity rivals icowsal JointAutoc dissimilarity interact_region_diss || regionalcode: regional_demstrength, cov(unstruct)
#   i.e. multilevel logit, random intercept + random slope on regional democratic strength by Lemke region (unstructured), intra-regional dyads.
# KEY TERM (a priori, abstract): "when the regional power balance favors democracies, it influences conflict between dyads, especially mixed
#   dyads" -> product term interact_region_diss (regional democratic strength x regime dissimilarity).
# Estimator: lme4::glmer binomial logit with (1 + regional_demstrength | regionalcode), unstructured; Laplace approximation (nAGQ = 1)
#   -- Stata xtmelogit default is 7-point adaptive quadrature, which lme4 cannot do with a random slope -> closest fast equivalent (flagged).
#   Same estimator in all six steps; history terms enter the fixed part.
CODE <- Sys.getenv("SWEEP_CODE"); stopifnot(nzchar(CODE))
source(file.path(CODE, "work_p8/common_b3.R"))
suppressPackageStartupMessages(library(lme4))
WD <- file.path(ROOT, "work_p8/P8_036"); W <- file.path(OUT, "work_p8/P8_036"); dir.create(W, recursive = TRUE, showWarnings = FALSE); sink(file.path(W, "log_P8_036.txt"), split = TRUE)
fx <- c("land_contiguity", "water_contiguity", "rivals", "icowsal", "JointAutoc", "dissimilarity", "interact_region_diss")
d <- haven::read_dta(file.path(WD, "regdem.dta"), col_select = all_of(c("ccode1", "ccode2", "year", "MID", fx, "regional_demstrength", "regionalcode")))
d <- as.data.frame(haven::zap_labels(d)); cat("rows", nrow(d), "\n")
# do-file data steps applied before Table 1 (duplicates with identical dyad-year-MID dropped; regionalcode 0/21 -> missing)
d$regionalcode[d$regionalcode %in% c(0, 21)] <- NA
d <- d[!is.na(d$regionalcode), ]
d <- d[!duplicated(d[, c("ccode1", "ccode2", "year", "MID")]), ]
fit_me <- function(rhs, data) {
  f <- as.formula(paste("MID ~", paste(rhs, collapse = " + "), "+ (1 + regional_demstrength | regionalcode)"))
  m <- glmer(f, data = data, family = binomial("logit"), nAGQ = 1,
             control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)))
  ct <- summary(m)$coefficients
  cat("\n---- glmer", deparse(f, width.cutoff = 500), "| N =", nobs(m), "| logLik =", round(as.numeric(logLik(m)), 3), "\n"); print(round(ct, 6))
  print(VarCorr(m)); msg <- m@optinfo$conv$lme4$messages; if (length(msg)) cat("CONVERGENCE MSG:", msg, "\n")
  out <- list(interact_region_diss = c(b = ct["interact_region_diss", 1], se = ct["interact_region_diss", 2], p = ct["interact_region_diss", 4]))
  attr(out, "n") <- nobs(m); out
}
vars <- c("MID", fx, "regional_demstrength", "regionalcode")
full <- d[complete.cases(d[, vars]), ]
r_full <- fit_me(fx, full); n_orig <- attr(r_full, "n")
ds <- common_sample(d, vars, "ccode1", "ccode2", "year")
M <- run_six(fx, ds, fit_me)
saveRDS(list(M = M, r_full = r_full, n_orig = n_orig, n_used = nrow(ds)), file.path(W, "M_P8_036.rds"))
# ---- port: result-CSV write step for this study, copied verbatim from work_p8/write_b3.R (paths: saved model list read from OUT) ----
NOPDF <- "Published PDF not accessible (paywall; the device that stages sweep PDFs was disconnected), and the archive has no log with the table numbers -> UNVERIFIED."
R <- readRDS(file.path(W, "M_P8_036.rds"))
meta <- list(pid = "P8_036", cite = "McCallister 2016", journal = "II", year = 2016, doi = "10.1080/03050629.2015.1069292", family = "DEMOCRACY")
write_p8(meta, "interact_region_diss", "abstract: 'when the regional power balance favors democracies, it influences conflict between dyads, especially mixed dyads'",
  "Table 1 Model 5 (do-file: xtmelogit MID ... interact_region_diss || regionalcode: regional_demstrength, cov(unstruct))",
  "multilevel logit (xtmelogit), random intercept + random slope on regional democratic strength by region, unstructured; OIM SE",
  "lme4::glmer, same random structure, Laplace (nAGQ=1) instead of Stata's 7-point adaptive quadrature (lme4 cannot do AGQ with random slopes) -- flagged, all six steps",
  NA, NA, NA, FALSE, FALSE, R$M, "N (rivalry dummy only; no MID history, no peace years)",
  sprintf("%s Authors' data (Dataverse KT4XKZ regdem file) + do-file. Intra-regional dyads (regionalcode non-missing, codes 0/21 set missing as in do-file). Full N=%d, b=%.4f p=%.3g; common 1822-2011 N=%d.",
    NOPDF, R$n_orig, R$r_full$interact_region_diss[["b"]], R$r_full$interact_region_diss[["p"]], R$n_used), verdict = "UNVERIFIED")
sink()
