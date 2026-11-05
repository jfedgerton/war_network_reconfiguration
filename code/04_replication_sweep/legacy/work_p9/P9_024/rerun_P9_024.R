# P9_024 Thyne 2010 JPR 47(4) "Supporter of stability or agent of agitation? The effect of US foreign policy on coups in
#   Latin America, 1960-99" (10.1177/0022343310368350)
# Data + do-file: PRIO JPR replication archive, Thyne_JPR_coups_data_appendix 47(4).zip (thyne_JPR_coups.dta, thyne_JPR_coups.do,
#   thyne_JPR_appendix.pdf); members extracted into DATA/work_p9/P9_024/. No Stata log in the archive.
# MONADIC country-MONTH panel (Latin American states, ccode x year x month). DV coup = coup attempt in the month.
# HEADLINE (do-file "Table II", first model): logit coup US_signals US_aid US_mids democracy military_reg instability civil_war
#   gdppc chgdp peace _spline1 _spline2 _spline3, cluster(ccode).  KEY TERM US_signals (net US->state COPDAB/WEIS signal score).
# Published values: none in the archive (no log); article PDF not accessible -> pub_b NA (verdict NO_PUBLISHED) until supplied.
# Estimator all steps: fixest feglm logit, vcov clustered on ccode. History: monadic convention (ADD_MON), joined on ccode-year,
#   so all months of a year carry that year's history and benchmark values.
CODE <- Sys.getenv("SWEEP_CODE"); stopifnot(nzchar(CODE))
source(file.path(CODE, "work_p8/common_b3.R"))
WD <- file.path(ROOT, "work_p9/P9_024"); W <- file.path(OUT, "work_p9/P9_024"); dir.create(W, recursive = TRUE, showWarnings = FALSE)
sink(file.path(W, "log_P9_024.txt"), split = TRUE)
d <- as.data.frame(haven::zap_labels(haven::read_dta(file.path(WD, "thyne_JPR_coups.dta"))))
stopifnot(nrow(d) == 14364, !anyDuplicated(d[, c("ccode", "year", "month")]), all(d$month %in% 1:12))
names(d) <- sub("^_", "", names(d))   # Stata _spline1.. -> spline1.. (R formula-safe)
stopifnot(all(c("spline1", "spline2", "spline3") %in% names(d)))
rhs0 <- c("US_signals", "US_aid", "US_mids", "democracy", "military_reg", "instability", "civil_war", "gdppc", "chgdp",
          "peace", "spline1", "spline2", "spline3")
vars <- c("coup", rhs0, "ccode")
fitter <- function(r, data) fit_bin("coup", r, data, "logit", "~ccode", terms = "US_signals")
d0 <- d[complete.cases(d[, vars]), ]
cat("full complete-case N:", nrow(d0), " years:", range(d0$year), " countries:", length(unique(d0$ccode)), " coups:", sum(d0$coup), "\n")
full <- fitter(rhs0, d0)
ds <- common_sample(d, vars, "ccode", NULL, "year", monadic = TRUE)
M <- run_six(rhs0, ds, fitter, monadic = TRUE)
saveRDS(list(M = M, full = full, n_used = nrow(ds)), file.path(W, "M_P9_024.rds"))
# sensitivity hook: NOT wired. sens_run() requires a dyadic history column (hist_dyad_w15), which the monadic convention does not
#   create; see tmp/p9/P9_024_proposal.md.

meta <- list(pid = "P9_024", cite = "Thyne 2010", journal = "JPR", year = 2010, doi = "10.1177/0022343310368350", family = "TIE")
write_p8(meta, "US_signals",
  "title: effect of US foreign policy signals on coups in Latin America",
  "Table II Model 1 (do-file order; net US signals)", "logit, SE clustered on country, peace months + 3 splines",
  "feglm logit, vcov clustered ccode, all six steps", NA, NA, NA, NA, FALSE, M, "N (peace + 3 splines; US_mids = MIDs with the US)",
  sprintf(paste0("MONADIC country-month (Latin America); mD not estimable. Published values pending (no log in archive; PDF needed). ",
    "Full-sample b=%.6f se=%.6f N=%d (years %d-%d). Common 1822-2011 sample N=%d; history joined on ccode-year."),
    full$US_signals[["b"]], full$US_signals[["se"]], attr(full, "n"), min(d0$year), max(d0$year), nrow(ds)))
dir.create(file.path(OUT, "results_p9"), recursive = TRUE, showWarnings = FALSE)
src <- file.path(OUT, "results_p8", "P9_024__US_signals.csv"); dst <- file.path(OUT, "results_p9", "P9_024__US_signals.csv")
stopifnot(file.exists(src), file.rename(src, dst))
sink()
