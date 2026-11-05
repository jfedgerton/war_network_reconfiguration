# P9_018 Allen, Campbell, Hernandez & Shepherd 2023 FPA 19(1) "US Military Deployments and the Risk of Coup d'Etat" (10.1093/fpa/orac027)
# Data + log: Harvard Dataverse doi:10.7910/DVN/IXD3DA ("primary_coup_replication.dta", "model design - RR.do",
#   "final_replication_log.smcl"); staged as adh_coup.dta / adh_model_design.do / adh_log.smcl.
# Country-years (host state ccode1), 1946-2019; tsset ccode1 year; all regressors lagged one year (Stata L.).
# HEADLINE (abstract: "the presence of US troops decreases the likelihood of a coup"): Table 1 Model 6 ("Cubic Polynomial"; the
#   full model the authors use for the margins plot) =
#   "logit coup l.lntroops100 l.dependence l.kappavv l.threat_env l.excluded l.ln_terrorist_attack l.civil_war l.lnmilex l.lnmilper
#    l.marxist l.coldwar l.time l.time2 l.time3 i.regime_case_fe, robust".
#   KEY TERM L.lntroops100. Author's log: -.1472121 (robust SE .062233), N = 3,646.
# Estimator: unconditional logit with regime-case dummies -> fixest feglm logit, fe = regime_case_fe (fixest drops groups that
#   predict failure/success perfectly, as Stata does), vcov hetero, all six steps.
# History linkage: dyad (US = 2, host ccode1) at the row year (the tie is US troops in the host state).
CODE <- Sys.getenv("SWEEP_CODE"); stopifnot(nzchar(CODE))
source(file.path(CODE, "work_p8/common_b3.R"))
WD <- file.path(ROOT, "work_p9/P9_018"); W <- file.path(OUT, "work_p9/P9_018"); dir.create(W, recursive = TRUE, showWarnings = FALSE)
sink(file.path(W, "log_P9_018.txt"), split = TRUE)
d <- as.data.frame(haven::zap_labels(haven::read_dta(file.path(WD, "adh_coup.dta"))))
stopifnot(nrow(d) == 11026, !anyDuplicated(d[, c("ccode1", "year")]), !any(d$ccode1 == 2))
cat("year range:", range(d$year), "\n")
lagv <- c("lntroops100", "dependence", "kappavv", "threat_env", "excluded", "ln_terrorist_attack", "civil_war", "lnmilex", "lnmilper",
          "marxist", "coldwar", "time", "time2", "time3")
# Stata L. under tsset ccode1 year: value of the same ccode1 at year - 1 (missing if that year is absent)
lg <- d[, c("ccode1", "year", lagv)]; lg$year <- lg$year + 1; names(lg)[-(1:2)] <- paste0("L_", lagv)
d <- left_join(d, lg, by = c("ccode1", "year")); stopifnot(nrow(d) == 11026)
d$ccode_us <- 2
term <- "L_lntroops100"
rhs0 <- paste0("L_", lagv)
vars <- c("coup", rhs0, "regime_case_fe")
fitter <- function(r, data) fit_bin("coup", r, data, "logit", "hetero", fe = "regime_case_fe", terms = term)
d0 <- d[complete.cases(d[, vars]), ]
full <- fitter(rhs0, d0)
cat("\nREPRO vs log:", term, "b", full[[term]][["b"]], "(log -.1472121) se", full[[term]][["se"]], "(log .062233) N", attr(full, "n"), "(log 3646)\n")
stopifnot(attr(full, "n") == 3646, abs(full[[term]][["b"]] - (-0.1472121)) / 0.1472121 <= 0.03)
# Gap check (verification only, not used in the six steps): Stata's log reports "950005.regime_case_fe omitted because of
# collinearity" (a 2-row regime case), which folds that case into the base category (lowest level in the sample, 41000).
# The full-rank model above has a higher log-likelihood; reproducing Stata's omission must return the logged b exactly.
m_full <- feglm(as.formula(paste("coup ~", paste(rhs0, collapse = " + "), "| regime_case_fe")), data = d0, family = binomial(), notes = FALSE)
dv <- d0[obs(m_full), ]; stopifnot(min(dv$regime_case_fe) == 41000, sum(dv$regime_case_fe == 950005) == 2)
dv$fe_stata <- ifelse(dv$regime_case_fe == 950005, 41000, dv$regime_case_fe)
m_st <- feglm(as.formula(paste("coup ~", paste(rhs0, collapse = " + "), "| fe_stata")), data = dv, family = binomial(), notes = FALSE)
cat("Stata-omission check: b", coef(m_st)[[term]], "(log -.1472121) logLik", logLik(m_st), "(log -917.76193); full-rank logLik",
    logLik(m_full), "\n")
stopifnot(abs(coef(m_st)[[term]] - (-0.1472121)) < 1e-6, abs(logLik(m_st) - (-917.76193)) < 1e-3, logLik(m_full) > logLik(m_st))
rm(m_full, m_st, dv); invisible(gc())
ds <- common_sample(d, vars, "ccode_us", "ccode1", "year")
M <- run_six(rhs0, ds, fitter)
saveRDS(list(M = M, full = full, n_used = nrow(ds)), file.path(W, "M_P9_018.rds"))

# ---- sensitivity hook: runs only under run_sensitivity.R (SWEEP_SENS=1), never in the main sweep ----
if (Sys.getenv("SWEEP_SENS") == "1") {
  source(file.path(CODE, "hist/sens_helper.R"))
  # L_lntroops100 enters the FE logit additively -> full sensitivity; FE logit drops all-0/all-1 groups -> fe_sep
  sens_run(id = "P9_018", dat = ds, y = "coup", tie = term, ctrl = setdiff(rhs0, term), fe = "regime_case_fe", fe_sep = TRUE,
           fit = function(dat, add) fitter(c(rhs0, add), dat)[[term]])
}

meta <- list(pid = "P9_018", cite = "Allen et al. 2023", journal = "FPA", year = 2023, doi = "10.1093/fpa/orac027", family = "TIE")
write_p8(meta, term,
  "abstract: 'the presence of US troops decreases the likelihood of a coup'",
  "Table 1 Model 6 (cubic polynomial, full controls; author's log)",
  "logit with regime-case dummies, robust SE, lagged regressors, cubic time",
  "feglm logit, fe regime_case_fe, vcov hetero, all six steps", -0.1472121, 0.062233, 3646, TRUE, FALSE, M,
  "N (time since last coup cubic only)",
  sprintf(paste0("Host country-years (Allen et al. data, %d-%d). Published values from the author's Stata log. Full-sample repro b=%.6f ",
    "se=%.6f N=%d vs log -.1472121 (.062233) N=3646. Common 1822-2011 sample N=%d (fixest count after dropping perfectly ",
    "predicted FE groups). History linked to dyad (2, host). Gap vs log (-1.9%%) is Stata's collinearity omission of the 2-row ",
    "regime dummy 950005 (folded into base 41000): reproducing that omission returns the logged b exactly (checked in script); ",
    "the six steps use the full-rank specification as written in the do-file."), min(d$year), max(d$year), full[[term]][["b"]], full[[term]][["se"]],
    attr(full, "n"), nrow(ds)))
dir.create(file.path(OUT, "results_p9"), recursive = TRUE, showWarnings = FALSE)
src <- file.path(OUT, "results_p8", paste0("P9_018__", term, ".csv")); dst <- file.path(OUT, "results_p9", paste0("P9_018__", term, ".csv"))
stopifnot(file.exists(src), file.rename(src, dst))
sink()
