# P9_008 Salehyan 2008 AJPS 52(4) "The Externalities of Civil Strife: Refugees as a Source of International Conflict" (10.1111/j.1540-5907.2008.00343.x)
# Data + do-file: Harvard Dataverse doi:10.7910/DVN/BGFB87 ("RefugeesWar_directed.dta", "refugee_regressions2.do");
#   staged as salehyan_directed.dta / salehyan.do. No log in the archive.
# Directed dyads 1955-2000 (year >= 1955), DV mzinit_lead = cc1 initiates a MID vs cc2 in t+1. Probit, cluster(dyad) (undirected dyad id),
#   peace years + 3 splines (lpcyrs lpcyrs1-3 = Stata lpcyrs*). civwar1/2 regenerated from uppcivcon1/2 as in the do-file
#   (the stored civwar1/2 columns differ and are overwritten, as Stata's `capture drop; gen` does).
# HEADLINE: model 1 (stock model; "model 1 provides evidence in support of both hypotheses").
#   KEY TERM logref1 = log(1 + refugees hosted in cc1 from cc2, '000s) (receiving state initiates; first proposition in the abstract).
#   SECOND ROW logref2 = refugees from cc1 hosted in cc2 (sending state initiates).
#   Published (author's working-paper copy, Table 2 Model 1): logref1 .057 (.019), logref2 .056 (.013), N = 86,497.
# Estimator: fixest feglm probit, vcov ~dyad, all six steps. Both rows share one model and one common sample.
# History controls and benchmarks are attached at the outcome year (year + 1), since the DV is a lead (fix 2026-10-10).
CODE <- Sys.getenv("SWEEP_CODE"); stopifnot(nzchar(CODE))
source(file.path(CODE, "work_p8/common_b3.R"))
WD <- file.path(ROOT, "work_p9/P9_008"); W <- file.path(OUT, "work_p9/P9_008"); dir.create(W, recursive = TRUE, showWarnings = FALSE)
sink(file.path(W, "log_P9_008.txt"), split = TRUE)
d <- as.data.frame(haven::zap_labels(haven::read_dta(file.path(WD, "salehyan_directed.dta"))))
stopifnot(nrow(d) == 116352, !anyDuplicated(d[, c("ccode1", "ccode2", "year")]))
d$civwar1 <- d$uppcivcon1; d$civwar2 <- d$uppcivcon2          # do-file: gen civwar1=uppcivcon1, civwar2=uppcivcon2
d <- d[d$year >= 1955, ]
cat("year range:", range(d$year), "\n")
terms <- c("logref1", "logref2")
rhs0 <- c("logref1", "logref2", "civwar1", "civwar2", "dem1", "dem2", "demdem", "trans1", "trans2", "transtrans", "contig", "colcont",
          "capshare", "s_wt_glo", "depend1", "depend2", "igos", "lpcyrs", "lpcyrs1", "lpcyrs2", "lpcyrs3")
vars <- c("mzinit_lead", rhs0, "dyad")
pub <- list(logref1 = c(b = 0.057, se = 0.019), logref2 = c(b = 0.056, se = 0.013)); pub_n <- 86497
fitter <- function(r, data) fit_bin("mzinit_lead", r, data, "probit", "~dyad", terms = terms)
d0 <- d[complete.cases(d[, vars]), ]
full <- fitter(rhs0, d0)
for (t in terms) cat("\nREPRO vs published:", t, "b", full[[t]][["b"]], "(pub", pub[[t]][["b"]], ") se", full[[t]][["se"]], "(pub", pub[[t]][["se"]], ")")
cat("\nN", attr(full, "n"), "(pub", pub_n, ")\n")
# published values are rounded to 3 decimals: require |b - pub| <= 0.0005 (rounding) and the exact N
stopifnot(attr(full, "n") == pub_n, all(sapply(terms, function(t) abs(full[[t]][["b"]] - pub[[t]][["b"]]) <= 0.0005)))
# History timing (fix 2026-10-10): the DV is MID initiation in t+1, so history (t-1..t-5 / t-2..t-6) and benchmark windows must
# end the year before the OUTCOME year. Attach them at outcome_year = year + 1 (sample-year filter 1822-2011 then applies to t+1).
d$outcome_year <- d$year + 1
# assert the lead: on consecutive directed-dyad records, mzinit_lead[t] == mzinit[t+1] (NA matching NA)
stopifnot("mzinit" %in% names(d))
o <- d[order(d$ccode1, d$ccode2, d$year), c("ccode1", "ccode2", "year", "mzinit", "mzinit_lead")]; no <- nrow(o)
cs <- o$ccode1[-1] == o$ccode1[-no] & o$ccode2[-1] == o$ccode2[-no] & o$year[-1] == o$year[-no] + 1
la <- o$mzinit_lead[-no][cs]; lb <- o$mzinit[-1][cs]
stopifnot(sum(cs) > 1000, identical(is.na(la), is.na(lb)), all(la[!is.na(la)] == lb[!is.na(lb)]))
cat("lead check: mzinit_lead[t] == mzinit[t+1] on", sum(cs), "consecutive directed-dyad pairs\n"); rm(o, cs, la, lb, no)
ds <- common_sample(d, vars, "ccode1", "ccode2", "outcome_year")
M <- run_six(rhs0, ds, fitter)
out <- list(M = M, full = full, n_used = nrow(ds))
saveRDS(out, file.path(W, "M_P9_008.rds"))

# ---- sensitivity hook: runs only under run_sensitivity.R (SWEEP_SENS=1), never in the main sweep ----
if (Sys.getenv("SWEEP_SENS") == "1") {
  source(file.path(CODE, "hist/sens_helper.R"))
  # headline: logref1 enters additively -> full sensitivity (logref2 stays among the controls)
  sens_run(id = "P9_008", dat = ds, y = "mzinit_lead", tie = "logref1", ctrl = setdiff(rhs0, "logref1"), fe = NULL,
           fit = function(dat, add) fitter(c(rhs0, add), dat)[["logref1"]])
}

meta <- list(pid = "P9_008", cite = "Salehyan 2008", journal = "AJPS", year = 2008, doi = "10.1111/j.1540-5907.2008.00343.x", family = "TIE")
art <- "probit, SE clustered on (undirected) dyad, peace years + 3 splines"
used <- "feglm probit, vcov clustered dyad, all six steps"
ahc <- "N (peace years + 3 splines only)"
nt <- function(t, extra) sprintf(paste0("Directed dyads 1955-2000 (author's data). Published values from the author's working-paper copy ",
  "(Table 2 Model 1, 3 decimals; AJPS table number may differ). Full-sample repro b=%.6f se=%.6f N=%d vs pub %s (%s) N=%d. ",
  "Common 1822-2011 sample N=%d (history/benchmarks attached at outcome year t+1). %s"), full[[t]][["b"]], full[[t]][["se"]], attr(full, "n"), format(pub[[t]][["b"]]),
  format(pub[[t]][["se"]]), pub_n, nrow(ds), extra)
write_p8(meta, "logref1",
  "abstract: 'Refugee-receiving states are more likely to initiate MIDs'",
  "Table 2 Model 1 (refugee stocks; working-paper numbering)", art, used,
  0.057, 0.019, pub_n, TRUE, FALSE, M, ahc, nt("logref1", "Headline (receiving-state term). logref2 kept in all steps."))
write_p8(meta, "logref2",
  "abstract: 'refugee-sending states initiate MIDs as they violate borders'",
  "Table 2 Model 1 (refugee stocks; working-paper numbering)", art, used,
  0.056, 0.013, pub_n, TRUE, FALSE, M, ahc, nt("logref2", "Non-headline second row (sending-state term) from the same model."))
dir.create(file.path(OUT, "results_p9"), recursive = TRUE, showWarnings = FALSE)
for (tm in terms) {
  src <- file.path(OUT, "results_p8", paste0("P9_008__", tm, ".csv")); dst <- file.path(OUT, "results_p9", paste0("P9_008__", tm, ".csv"))
  stopifnot(file.exists(src), file.rename(src, dst))
}
sink()
