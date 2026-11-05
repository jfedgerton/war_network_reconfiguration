# P9_046 Peterson 2014 JCR 58(4) "Dyadic Trade, Exit Costs, and Conflict" (10.1177/0022002713478794)
# Data + do-file: author page timothypeterson.org/research ("JCR-exit-costs-replication.zip": JCR_rep_data.dta, JCR_rep.do,
#   JCR_EC_appendix.pdf). No log in the archive. Erratum (author page) corrects an equation typo only (no estimates).
# Directed dyads 1984-2000. DV fincthc = count of high-conflict events A -> B in t+1 (do-file: "to lag IVs ... I created future
#   variables (t+1) for conflict and peace years"). zinb, inflate(lndist), cluster(dyadid).
# HEADLINE: Table 1 Model 1 (all commodities, share-based exit costs, high-conflict events):
#   "zinb fincthc lnsg1 lnsg2 lnsg_int lnval mingdp lndist caprat polity1_a polity2_a polity_i s_wt_glo fpy fpy2 fpy3, inflate(lndist) cluster(dyadid)"
#   KEY TERM lnsg_int = lnsg1 x lnsg2 (abstract: "jointly high exit costs are pacifying").
#   Stata abbreviations polity1_a / polity2_a / polity_i resolve to polity1_adj / polity2_adj / polity_int (only matches in the data).
#   PI DECISION 2026-10-10: HEADLINE is Table 1 Model 3 (strategic commodities, high-conflict events), do-file line 12:
#   "zinb fincthc lnsg_strat1 lnsg_strat2 lnsg_strat_int lnval_strat mingdp lndist caprat polity1_a polity2_a polity_i s_wt_glo fpy fpy2 fpy3, inflate(lndist) cluster(dyadid)"
#   KEY TERM lnsg_strat_int. Model 1 (lnsg_int) kept as a non-headline second row. Published values from the article PDF (Table 1).
# Estimator: ZINB ML (NB2 count, logit inflation on lndist + constant), cluster-robust SE on dyadid (Stata G/(G-1)),
#   custom ML in zinb_P9_046.R (validated against pscl::zeroinfl coefficients). History terms enter the count equation only.
# History controls and benchmarks are attached at the outcome year (year + 1), since the DV is a lead.
CODE <- Sys.getenv("SWEEP_CODE"); stopifnot(nzchar(CODE))
source(file.path(CODE, "work_p8/common_b3.R"))
source(file.path(CODE, "work_p9/P9_046/zinb_P9_046.R"))
WD <- file.path(ROOT, "work_p9/P9_046"); W <- file.path(OUT, "work_p9/P9_046"); dir.create(W, recursive = TRUE, showWarnings = FALSE)
sink(file.path(W, "log_P9_046.txt"), split = TRUE)
d <- as.data.frame(haven::zap_labels(haven::read_dta(file.path(WD, "JCR_rep_data.dta"))))
stopifnot(nrow(d) == 165006, !anyDuplicated(d[, c("ccode1", "ccode2", "year")]))
stopifnot(max(abs(d$lnsg_int - d$lnsg1 * d$lnsg2)) < 1e-5)
stopifnot(max(abs(d$lnsg_strat_int - d$lnsg_strat1 * d$lnsg_strat2), na.rm = TRUE) < 1e-5)
cat("year range:", range(d$year), "\n")
ctrl <- c("mingdp", "lndist", "caprat", "polity1_adj", "polity2_adj", "polity_int", "s_wt_glo", "fpy", "fpy2", "fpy3")
specs <- list(
  M3 = list(term = "lnsg_strat_int", rhs = c("lnsg_strat1", "lnsg_strat2", "lnsg_strat_int", "lnval_strat", ctrl), pub = c(-0.437, 0.069, 165006),
            hm = "Table 1 Model 3 (share-based exit costs, STRATEGIC commodities, high-conflict events; do-file line 12)",
            ev = "abstract: jointly high exit costs pacifying, 'a pattern most prominent for trade in strategic commodities'",
            note = "HEADLINE (PI decision 2026-10-10). Published Table 1 Model 3 -0.437*** (0.069), N 165,006; the model fully supporting H3, used for Fig 3 and the US-China example."),
  M1 = list(term = "lnsg_int", rhs = c("lnsg1", "lnsg2", "lnsg_int", "lnval", ctrl), pub = c(-0.138, 0.034, 165006),
            hm = "Table 1 Model 1 (share-based exit costs, all commodities, high-conflict events; do-file m1)",
            ev = "abstract: 'jointly high exit costs are pacifying'",
            note = "Non-headline second row (article's 'baseline model'). Published Table 1 Model 1 -0.138*** (0.034), N 165,006."))
mk_fitter <- function(term) function(r, data) {
  X <- cbind(`(Intercept)` = 1, as.matrix(data[, r])); Z <- cbind(`(Intercept)` = 1, lndist = data$lndist)
  f <- fit_zinb(X, Z, data$fincthc, cl = data$dyadid)
  se <- sqrt(diag(f$V)); z <- f$coef / se
  cat("\n---- zinb fincthc ~", paste(r, collapse = " + "), "| inflate(lndist) | N =", f$n, "| clusters:", f$G, "| ll =", f$ll, "\n")
  print(round(cbind(b = f$coef, se = se, p = 2 * pnorm(-abs(z))), 6))
  out <- list(c(b = unname(f$coef[term]), se = unname(se[term]), p = unname(2 * pnorm(-abs(z[term])))))
  names(out) <- term; attr(out, "n") <- f$n; out
}
d$outcome_year <- d$year + 1
meta <- list(pid = "P9_046", cite = "Peterson 2014", journal = "JCR", year = 2014, doi = "10.1177/0022002713478794", family = "TIE")
art <- "zero-inflated negative binomial (inflate: lndist), SE clustered on dyad, future peace years + square + cube"
used <- "custom ZINB ML (NB2 + logit inflation on lndist), cluster-robust SE dyadid, all six steps; history in count equation"
ahc <- "N (peace years + square + cube only)"
dir.create(file.path(OUT, "results_p9"), recursive = TRUE, showWarnings = FALSE)
res <- list()
for (s in names(specs)) {
  sp <- specs[[s]]; fitter <- mk_fitter(sp$term); vars <- c("fincthc", sp$rhs, "dyadid")
  cat("\n==================== spec", s, "key term", sp$term, "\n")
  d0 <- d[complete.cases(d[, vars]), ]
  full <- fitter(sp$rhs, d0)
  cat("\nFULL-SAMPLE:", sp$term, "b", full[[sp$term]][["b"]], "se", full[[sp$term]][["se"]], "N", attr(full, "n"), "| published", sp$pub, "\n")
  ds <- common_sample(d, vars, "ccode1", "ccode2", "outcome_year")
  M <- run_six(sp$rhs, ds, fitter)
  res[[s]] <- list(M = M, full = full, n_used = nrow(ds))
  write_p8(meta, sp$term, sp$ev, sp$hm, art, used, sp$pub[1], sp$pub[2], as.integer(sp$pub[3]), TRUE, FALSE, M, ahc,
    sprintf(paste0("Directed dyads 1984-2000 (author's data). %s Full-sample b=%.6f se=%.6f N=%d. Common 1822-2011 sample N=%d ",
      "(history/benchmarks attached at outcome year t+1). Interaction: both constituent terms kept in all steps."),
      sp$note, full[[sp$term]][["b"]], full[[sp$term]][["se"]], attr(full, "n"), nrow(ds)))
  src <- file.path(OUT, "results_p8", paste0("P9_046__", sp$term, ".csv")); dst <- file.path(OUT, "results_p9", basename(src))
  stopifnot(file.exists(src), file.rename(src, dst))
  rm(d0, ds, M, full); gc()
}
saveRDS(res, file.path(W, "M_P9_046.rds"))
sink()
