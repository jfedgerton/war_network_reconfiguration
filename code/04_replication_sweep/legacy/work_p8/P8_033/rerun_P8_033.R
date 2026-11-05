# P8_033 Colgan & Weeks 2015 IO "Revolution, Personalist Dictatorships, and International Conflict" (10.1017/s0020818314000307)
# Data + do-file: Cambridge supplementary zip S0020818314000307sup001.zip (sup003.dta = ColganWeeksDyadic-public, sup002.do).
# HEADLINE: Table 1, model 1.1 (do-file): logit mzinit on revpers10_lag1 rev10_lag1 pers10_lag1 + controls + pcyrsmzinit + 3 splines,
#   robust cluster(dirdyadid). Directed dyad-years.
# KEY TERM (fixed a priori from abstract): revpers10_lag1 -- "revolutions that result in personalist dictatorships are significantly
#   more likely to lead to international conflict".
CODE <- Sys.getenv("SWEEP_CODE"); stopifnot(nzchar(CODE))
source(file.path(CODE, "work_p8/common_b3.R"))
WD <- file.path(ROOT, "work_p8/P8_033"); W <- file.path(OUT, "work_p8/P8_033"); dir.create(W, recursive = TRUE, showWarnings = FALSE); sink(file.path(W, "log_P8_033.txt"), split = TRUE)
meta <- list(pid = "P8_033", cite = "Colgan & Weeks 2015", journal = "IO", year = 2015, doi = "10.1017/s0020818314000307", family = "DEMOCRACY")
rhs0 <- c("revpers10_lag1", "rev10_lag1", "pers10_lag1", "cap_1_lag", "cap_2_lag", "initshare_lag", "dependlow_lag", "civilwar_1_lag",
          "civilwar_2_lag", "majmaj_lag", "minmaj_lag", "majmin_lag", "contigdum_lag", "logdist_lag", "s_wt_glo_lag", "s_lead_1_lag",
          "s_lead_2_lag", "pcyrsmzinit", "pcyrsmzinits1", "pcyrsmzinits2", "pcyrsmzinits3")
d <- haven::read_dta(file.path(WD, "S0020818314000307sup003.dta"), col_select = all_of(c("ccode1", "ccode2", "year", "dirdyadid", "mzinit", rhs0)))
d <- as.data.frame(haven::zap_labels(d)); cat("rows", nrow(d), "years", range(d$year), "\n")
vars <- c("mzinit", rhs0, "dirdyadid")
fitter <- function(rhs, data) fit_bin("mzinit", rhs, data, "logit", "~dirdyadid", terms = "revpers10_lag1")
# 1. reproduce on authors' full estimation sample
full <- d[complete.cases(d[, vars]), ]
r_full <- fitter(rhs0, full); n_orig <- attr(r_full, "n")
# 2. common 1822-2011 sample + history
ds <- common_sample(d, vars, "ccode1", "ccode2", "year"); rm(d, full); invisible(gc())
M <- run_six(rhs0, ds, fitter)
# ---- sensitivity hook: runs only under run_sensitivity.R (SWEEP_SENS=1), never in the main sweep ----
if (Sys.getenv("SWEEP_SENS") == "1") {
  source(file.path(CODE, "hist/sens_helper.R"))
  sens_run(id = "P8_033", dat = ds, y = "mzinit", tie = "revpers10_lag1", ctrl = setdiff(rhs0, "revpers10_lag1"), fe = NULL,
           fit = function(dat, add) fitter(c(rhs0, add), dat)[["revpers10_lag1"]])
}
saveRDS(list(M = M, r_full = r_full, n_orig = n_orig, n_used = nrow(ds)), file.path(W, "M_P8_033.rds"))
# ---- port: result-CSV write step for this study, copied verbatim from work_p8/write_b3.R (paths: saved model list read from OUT) ----
NOPDF <- "Published PDF not accessible (paywall; the device that stages sweep PDFs was disconnected), and the archive has no log with the table numbers -> UNVERIFIED."
R <- readRDS(file.path(W, "M_P8_033.rds"))
meta <- list(pid = "P8_033", cite = "Colgan & Weeks 2015", journal = "IO", year = 2015, doi = "10.1017/s0020818314000307", family = "DEMOCRACY")
write_p8(meta, "revpers10_lag1", "abstract: 'revolutions that result in personalist dictatorships are significantly more likely to lead to international conflict'",
  "Table 1 Model 1.1 (do-file: logit mzinit, all MIDs, all revolutions)", "logit, robust SE clustered on directed dyad, peace years + 3 splines; directed dyad-years",
  "logit (fixest), SE clustered dirdyadid, all six steps", NA, NA, NA, FALSE, FALSE, R$M, "N (peace years + splines only)",
  sprintf("%s Authors' data + do-file (Cambridge supplementary zip). Full-sample N=%d, b=%.4f p=%.3g. Common 1822-2011 sample N=%d.",
    NOPDF, R$n_orig, R$r_full$revpers10_lag1[["b"]], R$r_full$revpers10_lag1[["p"]], R$n_used), verdict = "UNVERIFIED")
sink()
