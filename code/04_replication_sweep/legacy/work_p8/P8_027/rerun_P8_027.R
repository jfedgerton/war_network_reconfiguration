# P8_027 Lai & Slater 2006 AJPS "Institutions of the Offensive" (10.1111/j.1540-5907.2006.00173.x)
# Data + log: Harvard Dataverse doi:10.7910/DVN/XTUHQG (ajpsreplicate.dta, ajpsreplication.smcl, laislaterart2.doc = tables).
# HEADLINE: Table 1 Model 1 (= "replicate model 1" in the authors' log): nbreg initone initonelag total cap open totalally party military,
#   cluster(ccode). MONADIC state-year 1950-92, nondemocracies + democracies (democracy = reference category).
# KEY TERM (a priori, abstract): "military regimes are more likely than single-party regimes to initiate military disputes" ->
#   military (published .516, robust SE .224, p<.05, N=4960). The military-vs-party contrast is logged as a secondary quantity.
# Estimator: negative binomial (NB2) ML with country-clustered SE -> fixest::fenegbin, vcov ~ccode, all six steps.
CODE <- Sys.getenv("SWEEP_CODE"); stopifnot(nzchar(CODE))
source(file.path(CODE, "work_p8/common_b3.R"))
WD <- file.path(ROOT, "work_p8/P8_027"); W <- file.path(OUT, "work_p8/P8_027"); dir.create(W, recursive = TRUE, showWarnings = FALSE); sink(file.path(W, "log_P8_027.txt"), split = TRUE)
d <- as.data.frame(haven::zap_labels(haven::read_dta(file.path(WD, "ajpsreplicate.dta"))))
rhs0 <- c("initonelag", "total", "cap", "open", "totalally", "party", "military")
fit_nb <- function(rhs, data) {
  m <- fenegbin(as.formula(paste("initone ~", paste(rhs, collapse = " + "))), data = data, vcov = ~ccode, notes = FALSE)
  stopifnot(m$convStatus)
  ct <- coeftable(m); cat("\n---- nbreg initone ~", paste(rhs, collapse = " + "), "| N =", nobs(m), "| theta =", round(m$theta, 4), "\n"); print(round(ct, 6))
  V <- vcov(m); con <- unname(coef(m)["military"] - coef(m)["party"])
  cse <- sqrt(V["military", "military"] + V["party", "party"] - 2 * V["military", "party"])
  cat("contrast military - party:", round(con, 4), "se", round(cse, 4), "p", signif(2 * pnorm(-abs(con / cse)), 3), "\n")
  out <- list(military = c(b = ct["military", 1], se = ct["military", 2], p = ct["military", 4]))
  attr(out, "n") <- nobs(m); out
}
vars <- c("initone", rhs0, "ccode")
full <- d[complete.cases(d[, vars]), ]
r_full <- fit_nb(rhs0, full); n_orig <- attr(r_full, "n")
cat("\nREPRO vs log: military b", r_full$military[["b"]], "(log .516266) se", r_full$military[["se"]], "(log .224477) N", n_orig, "(log 4960)\n")
ds <- common_sample(d, vars, "ccode", NULL, "year", monadic = TRUE)
M <- run_six(rhs0, ds, fit_nb, monadic = TRUE)
saveRDS(list(M = M, r_full = r_full, n_orig = n_orig, n_used = nrow(ds)), file.path(W, "M_P8_027.rds"))
meta <- list(pid = "P8_027", cite = "Lai & Slater 2006", journal = "AJPS", year = 2006, doi = "10.1111/j.1540-5907.2006.00173.x", family = "DEMOCRACY")
write_p8(meta, "military", "abstract: 'military regimes are more likely than single-party regimes to initiate military disputes'",
  "Table 1 Model 1 (nbreg, any MID initiated, 1950-92)", "negative binomial (NB2), SE clustered on state; lagged DV",
  "fenegbin (NB2 ML), vcov clustered ccode, all six steps", 0.516, 0.224, 4960, TRUE, FALSE, M, "Y (lagged DV initonelag)",
  sprintf(paste0("MONADIC state-year (history = own-state MID onsets w15/w26; mD NA, mDM = +own w15, m26 = +own w26; bC = own CINC, bD = own polity>=6). ",
  "Full-sample repro b=%.4f se=%.4f N=%d vs author log .516266 (.224477) N=4960. Common 1822-2011 sample N=%d. Reference category = democracy; ",
  "military-vs-party contrast printed in log. Authors' model already contains lagged DV (count of initiations t-1)."),
  r_full$military[["b"]], r_full$military[["se"]], n_orig, nrow(ds)))
sink()
