# P9_054 Peterson & Venteicher 2013 FPA 9(2) "Trade Relationships and Asymmetric Crisis Perception" (10.1111/j.1743-8594.2012.00197.x)
# Data + do-file: author's research page (Dropbox fpa-asymmetry-replication.zip: fpa_rep.do, fpa_rep.dta). No log. Article closed access
#   and no open copy found -> published values PENDING (verdict NO_PUBLISHED).
# Directed dyads (ccode1 = crisis initiator / first perceiver, ccode2 = target), 1918-2001, ICB crises.
# HEADLINE (abstract's first claim: unilateral perception persists when the initiator does not depend on trade with the target):
#   Table 1 Model 2 = streg ldep1_OR ldep2_OR lcapratio lndist ljdem ls grav_dum trigtdum, dist(lognormal) cluster(dyadid)
#   on the stored stset (_st == 1: 594 crisis dyads; _t = days until mutual perception, _d = failure). Model 2 is the column the
#   do-file interprets (margins / stcurve after m2). KEY TERM ldep1_OR (initiator's Oneal-Russett trade dependence on the target).
# SECOND ROW (abstract's second claim: a state first perceives a threat as a target when its dependence is high):
#   Table 2 Model 5 = probit crisis ldep1_OR ldep2_OR lcapratio lndist ljdem ls cyear cy2 cy3, cluster(dyadid) (cyear abbreviates
#   cyears). KEY TERM ldep2_OR (do-file margins after m5 are over ldep2_OR). All directed dyad-years.
# Estimators: lognormal AFT by survival::survreg with Stata cluster-robust V (sandwich x G/(G-1)); fixest feglm probit vcov ~dyadid.
#   Each spec has its own common sample. DVs are contemporaneous with the row year (no lead).
CODE <- Sys.getenv("SWEEP_CODE"); stopifnot(nzchar(CODE))
source(file.path(CODE, "work_p8/common_b3.R"))
suppressPackageStartupMessages(library(survival))
WD <- file.path(ROOT, "work_p9/P9_054"); W <- file.path(OUT, "work_p9/P9_054"); dir.create(W, recursive = TRUE, showWarnings = FALSE)
sink(file.path(W, "log_P9_054.txt"), split = TRUE)
d <- as.data.frame(haven::zap_labels(haven::read_dta(file.path(WD, "fpa_rep.dta"))))
stopifnot(nrow(d) == 1216222, sum(d$`_st` == 1) == 594, all(d$`_t0`[d$`_st` == 1] == 0))
names(d)[names(d) == "_t"] <- "st_t"; names(d)[names(d) == "_d"] <- "st_d"; names(d)[names(d) == "_st"] <- "st_st"
cat("year range:", range(d$year), "\n")

fit_aft <- function(rhs, data, terms) {
  f <- as.formula(paste("Surv(st_t, st_d) ~", paste(rhs, collapse = " + ")))
  m <- survreg(f, data = data, dist = "lognormal", robust = TRUE, cluster = dyadid)
  stopifnot(m$iter < 30)
  G <- length(unique(data$dyadid[complete.cases(data[, c("st_t", "st_d", rhs)])])); V <- m$var * G / (G - 1)
  b <- coef(m); se <- setNames(sqrt(diag(V))[seq_along(b)], names(b)); z <- b / se
  cat("\n---- streg lognormal", deparse(f, width.cutoff = 500), "| N =", length(m$linear.predictors), "| clusters", G, "| log sigma", round(log(m$scale), 6), "\n")
  print(round(cbind(b = b, se = se, z = z), 6))
  out <- lapply(setNames(terms, terms), function(t) c(b = b[[t]], se = se[[t]], p = 2 * pnorm(-abs(z[[t]]))))
  attr(out, "n") <- length(m$linear.predictors); out
}

ctl1 <- c("lcapratio", "lndist", "ljdem", "ls")
specs <- list(
  T1M2 = list(term = "ldep1_OR", dv = "st_t", rhs = c("ldep1_OR", "ldep2_OR", ctl1, "grav_dum", "trigtdum"),
              sub = function(x) x[x$st_st == 1, ], fitter = function(r, data) fit_aft(r, data, "ldep1_OR"), extra = c("st_d")),
  T2M5 = list(term = "ldep2_OR", dv = "crisis", rhs = c("ldep1_OR", "ldep2_OR", ctl1, "cyears", "cy2", "cy3"),
              sub = function(x) x, fitter = function(r, data) fit_bin("crisis", r, data, "probit", "~dyadid", terms = "ldep2_OR"), extra = character(0)))
out <- list()
for (s in names(specs)) {
  sp <- specs[[s]]; rhs0 <- sp$rhs; vars <- c(sp$dv, sp$extra, rhs0, "dyadid")
  cat("\n==== spec", s, "\n")
  ds0 <- sp$sub(d)
  full <- sp$fitter(rhs0, ds0[complete.cases(ds0[, vars]), ])
  cat("\nFULL-SAMPLE:", sp$term, "b", full[[sp$term]][["b"]], "se", full[[sp$term]][["se"]], "N", attr(full, "n"), "(published values pending)\n")
  ds <- common_sample(ds0, vars, "ccode1", "ccode2", "year")
  out[[s]] <- list(M = run_six(rhs0, ds, sp$fitter), full = full, n_used = nrow(ds))
  rm(ds0, ds); invisible(gc())
}
saveRDS(out, file.path(W, "M_P9_054.rds"))

# ---- sensitivity hook: runs only under run_sensitivity.R (SWEEP_SENS=1), never in the main sweep ----
if (Sys.getenv("SWEEP_SENS") == "1") {
  source(file.path(CODE, "hist/sens_helper.R"))
  sp <- specs$T1M2; dsh <- common_sample(sp$sub(d), c(sp$dv, sp$extra, sp$rhs, "dyadid"), "ccode1", "ccode2", "year")
  sens_run(id = "P9_054", dat = dsh, y = "st_t", tie = "ldep1_OR", ctrl = setdiff(sp$rhs, "ldep1_OR"), fe = NULL,
           fit = function(dat, add) sp$fitter(c(sp$rhs, add), dat)[["ldep1_OR"]])
}

meta <- list(pid = "P9_054", cite = "Peterson 2013", journal = "FPA", year = 2013, doi = "10.1111/j.1743-8594.2012.00197.x", family = "TIE")
ahc <- "N (none in duration model; crisis-years cubic in probit)"
nt <- function(s, extra) { o <- out[[s]]; sp <- specs[[s]]
  sprintf(paste0("Directed dyads 1918-2001 (author's data). Published values from the article PDF (Tables 1-2). Full-sample b and N exact; the probit (ldep2_OR) clustered SE is 9%% below the published SE (0.0997 vs 0.110), within the 1-18%% SE gap also seen in P9_047; kept as is (PI decision 2026-10-10). Full-sample b=%.6f se=%.6f N=%d. ",
    "Common 1822-2011 sample N=%d. %s"), o$full[[sp$term]][["b"]], o$full[[sp$term]][["se"]], attr(o$full, "n"), o$n_used, extra) }
write_p8(meta, "ldep1_OR",
  "abstract: 'Unilateral crisis perception is more likely to persist when the initiator does not depend on trade'",
  "Table 1 Model 2 (lognormal AFT, time to mutual crisis perception; do-file)", "streg lognormal AFT, SE clustered on dyadid",
  "survreg lognormal AFT, Stata cluster-robust V (dyadid), all six steps",
  -4.886, 2.558, 425, FALSE, FALSE, out$T1M2$M, ahc, nt("T1M2", "Headline (Model 2 is the Fig 1 source). Published Table 1 Model 2 -4.886 (2.558), N 425, p<.10 only (two-tailed). AFT: negative b = shorter unilateral perception."))
write_p8(meta, "ldep2_OR",
  "abstract: 'a state is more likely first to perceive a threat ... when its trade dependence ... is high'",
  "Table 2 Model 5 (probit crisis onset, all directed dyad-years; do-file)", "probit, SE clustered on dyadid, crisis-years cubic",
  "feglm probit, vcov clustered dyadid, all six steps",
  0.242, 0.110, 1082087, TRUE, FALSE, out$T2M5$M, ahc, nt("T2M5", "Non-headline second row (target-side dependence, crisis onset). Published Table 2 Model 5 0.242* (0.110), N 1,082,087."))
dir.create(file.path(OUT, "results_p9"), recursive = TRUE, showWarnings = FALSE)
for (tm in c("ldep1_OR", "ldep2_OR")) {
  src <- file.path(OUT, "results_p8", paste0("P9_054__", tm, ".csv")); dst <- file.path(OUT, "results_p9", paste0("P9_054__", tm, ".csv"))
  stopifnot(file.exists(src), file.rename(src, dst))
}
sink()
