# P9_051 Rahman 2026 JPR "Reliability as resolve: reputational transference in international relations" (10.1093/jopres/xjag019)
#   (not the Rahman 2026 JCR coup paper, P9_001)
# Data + do-file: Harvard Dataverse doi:10.7910/DVN/C4LZBG ("jpr_replication_data.dta", "dofile_models.do", "Online Appendix.pdf").
#   No log in the archive.
# Directed dyads (challenger -> target) 1816-2000. DV dispute = challenger initiates a MID against the target in year t
#   (outcome = the target being targeted in interstate conflict; abstract: poor reliability "increases the likelihood of being targeted").
# HEADLINE: Table 1 Model 2 (all dyads, full controls):
#   "logit dispute L.poor_rep_t_def ptargdef pchaloff pchalneu capprop majorpower jdem s_un_glo ln_distance cold_war peaceyrs
#    peaceyrs2 peaceyrs3, vce(robust)"  with xtset id year.
#   KEY TERM L.poor_rep_t_def = target's reputation as an unreliable (defense) ally, lagged one year (gap-aware lag on id-year,
#   asserted identical to the stored l1_poor_rep_t_def). Column name here: l1_poor_rep_t_def.
#   PUBLISHED VALUES PENDING: article page 403 from the shell; the online appendix does not repeat Table 1; no author log.
#   Appendix N for full-control lagged models = 1,045,074 (Tables A3-A5), matched exactly below.
# Estimator: fixest feglm logit, vcov "hetero" (Stata vce(robust)), all six steps. DV contemporaneous -> history at year.
CODE <- Sys.getenv("SWEEP_CODE"); stopifnot(nzchar(CODE))
source(file.path(CODE, "work_p8/common_b3.R"))
WD <- file.path(ROOT, "work_p9/P9_051"); W <- file.path(OUT, "work_p9/P9_051"); dir.create(W, recursive = TRUE, showWarnings = FALSE)
sink(file.path(W, "log_P9_051.txt"), split = TRUE)
d <- as.data.frame(haven::zap_labels(haven::read_dta(file.path(WD, "jpr_replication_data.dta"))))
stopifnot(nrow(d) == 1077992, !anyDuplicated(d[, c("id", "year")]), !anyDuplicated(d[, c("challenger", "target", "year")]))
# Stata L. under xtset id year: value at (id, year - 1), NA if that record is absent
k <- d[, c("id", "year", "poor_rep_t_def")]; k$year <- k$year + 1; names(k)[3] <- "lag_chk"
d <- left_join(d, k, by = c("id", "year")); stopifnot(nrow(d) == 1077992)
stopifnot(identical(is.na(d$lag_chk), is.na(d$l1_poor_rep_t_def)), all(d$lag_chk == d$l1_poor_rep_t_def, na.rm = TRUE))
cat("year range:", range(d$year), "\n")
term <- "l1_poor_rep_t_def"
rhs0 <- c("l1_poor_rep_t_def", "ptargdef", "pchaloff", "pchalneu", "capprop", "majorpower", "jdem", "s_un_glo", "ln_distance", "cold_war",
          "peaceyrs", "peaceyrs2", "peaceyrs3")
vars <- c("dispute", rhs0)
pub_n <- 1045074
fitter <- function(r, data) fit_bin("dispute", r, data, "logit", "hetero", terms = term)
d0 <- d[complete.cases(d[, vars]), ]
full <- fitter(rhs0, d0); rm(d0); invisible(gc())
cat("\nFULL-SAMPLE:", term, "b", full[[term]][["b"]], "se", full[[term]][["se"]], "N", attr(full, "n"), "(appendix N", pub_n, "; b/SE pending)\n")
stopifnot(attr(full, "n") == pub_n)
ds <- common_sample(d, vars, "challenger", "target", "year")
M <- run_six(rhs0, ds, fitter)
out <- list(M = M, full = full, n_used = nrow(ds))
saveRDS(out, file.path(W, "M_P9_051.rds"))

meta <- list(pid = "P9_051", cite = "Rahman 2026", journal = "JPR", year = 2026, doi = "10.1093/jopres/xjag019", family = "TIE")
art <- "logit, robust (Huber-White) SE, peace years + square + cube; directed dyads 1816-2000"
used <- "feglm logit, vcov hetero, all six steps"
ahc <- "N (peace years + square + cube only)"
write_p8(meta, term,
  "abstract: poor reliability 'significantly increases the likelihood of being targeted'",
  "Table 1 Model 2 (lagged reputation, full controls, all directed dyads)", art, used,
  0.645, 0.132, pub_n, TRUE, FALSE, M, ahc,
  sprintf(paste0("Directed dyads 1816-2000 (author's data). Published Table 1 Model 2 0.645*** (0.132), N 1,045,074 (article PDF). ",
    "Full-sample b=%.6f se=%.6f N=%d. Common 1822-2011 sample N=%d."),
    full[[term]][["b"]], full[[term]][["se"]], attr(full, "n"), nrow(ds)))
dir.create(file.path(OUT, "results_p9"), recursive = TRUE, showWarnings = FALSE)
src <- file.path(OUT, "results_p8", paste0("P9_051__", term, ".csv")); dst <- file.path(OUT, "results_p9", paste0("P9_051__", term, ".csv"))
stopifnot(file.exists(src), file.rename(src, dst))
sink()
