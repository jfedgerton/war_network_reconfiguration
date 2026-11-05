# P9_092 Millerd 2023, International Interactions 49(6) "Peaceful conflict resolution through densely gender-equal international organizations"
#   (10.1080/03050629.2023.2277959). TIE (shared membership in densely gender-equal IGOs). PORT of the author's archive.
# Archive: Harvard Dataverse doi:10.7910/DVN/I9ZSPR (densely-ge-replication.do; revthreshold50_io_gender_mid.dta, 204 MB).
# Found in the 2026-10-11 round-2 full-text re-screen: the abstract frames the outcome as non-violent dispute resolution, but the dyadic
#   models' DV is MID onset (COW MID 5.0, wars set missing) among politically relevant dyads 1989-2014.
# HEADLINE (fixed before estimating): Table 2 Model 3 (H2a, all PR dyads), Number of shared densely gender-equal IOs -0.007 (0.003) **, N 50,216.
#   Do-file: xtset dyad_id; xtlogit mid gendense_lag numio_lag libdem1 libdem2 inclusive1 inclusive2 defense_lag cont_lag kmdist_lag trade_lag
#   cinc_lag cumulative_lag time time2 time3 if poli_relevant==1   (random-effects logit; Stata default 12-point adaptive Gauss-Hermite).
#   Stars: * p<.05 ** p<.01 *** p<.001 -> alpha_pub 0.05.
# Estimator here: lme4::glmer binomial logit, random intercept by dyad_id, nAGQ = 12 (adaptive GH, as Stata), OIM SE; all six steps.
# Common sample: 79.1% of complete cases (2012-2014 outside the history window); named guard exception 0.75 (PI 2026-10-11).
# Runtime: ~20 min per glmer fit, ~2.5 h for the chain -> slow = Y, run via _slow/P9_092.R, cached.
# Data note: the slim .rds (22 columns of the 204 MB .dta) is a manifested derived input under SWEEP_DATA (md5 checked, never rewritten);
#   without it the script rebuilds the extract from the .dta into SWEEP_OUT (Codex fix 2026-10-11).
ROOT <- Sys.getenv("SWEEP_DATA"); OUT <- Sys.getenv("SWEEP_OUT"); CODE <- Sys.getenv("SWEEP_CODE")
stopifnot(nzchar(ROOT), nzchar(OUT), nzchar(CODE))
source(file.path(CODE, "work_p8/common_b3.R"))
suppressPackageStartupMessages(library(lme4)); cat("lme4", as.character(packageVersion("lme4")), "\n")
W <- file.path(OUT, "work_p9/P9_092"); dir.create(W, recursive = TRUE, showWarnings = FALSE)
sink(file.path(W, "log_P9_092.txt"), split = TRUE)
V <- c("ccode1", "ccode2", "year", "dyad_id", "mid", "gendense_lag", "numio_lag", "libdem1", "libdem2", "inclusive1", "inclusive2", "defense_lag",
       "cont_lag", "kmdist_lag", "trade_lag", "cinc_lag", "cumulative_lag", "time", "time2", "time3", "poli_relevant", "genderweak")
fd <- file.path(ROOT, "work_p9/P9_092/revthreshold50_io_gender_mid.dta"); fs <- file.path(ROOT, "work_p9/P9_092/revthreshold50_P9_092_slim.rds")
# SWEEP_DATA is read-only here (Codex fix 2026-10-11). The slim extract is a manifested derived input (data_sources.csv, md5 below).
# If it is absent, it is rebuilt from the author's .dta into SWEEP_OUT (never into SWEEP_DATA) and checked for identical content.
if (file.exists(fs)) {
  stopifnot(unname(tools::md5sum(fs)) == "21489ac5b98c0e36baef3806daf5420c")   # immutable derived input
  d <- readRDS(fs)
} else {
  stopifnot(file.exists(fd), unname(tools::md5sum(fd)) == "b34d96e1376f8f976cbb850af9d16d68")
  d <- as.data.frame(haven::zap_labels(haven::read_dta(fd, col_select = all_of(V))))
  fo <- file.path(W, "revthreshold50_P9_092_slim.rds"); saveRDS(d, fo); cat("slim extract rebuilt into SWEEP_OUT:", fo, "\n")
}
stopifnot(nrow(d) == 458713, all(V %in% names(d)))
d <- d[d$poli_relevant %in% 1, ]
TERM <- "gendense_lag"
rhs0 <- c("gendense_lag", "numio_lag", "libdem1", "libdem2", "inclusive1", "inclusive2", "defense_lag", "cont_lag", "kmdist_lag", "trade_lag",
          "cinc_lag", "cumulative_lag", "time", "time2", "time3")
vars <- c("mid", rhs0, "dyad_id")
fit_re <- function(rhs, data) {
  f <- as.formula(paste("mid ~", paste(rhs, collapse = " + "), "+ (1 | dyad_id)"))
  # autoscale: lme4 centres/scales the fixed-effect columns internally and reports coefficients on the original scale
  m <- glmer(f, data = data, family = binomial("logit"), nAGQ = 12,
             control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5), autoscale = TRUE))
  ct <- summary(m)$coefficients
  cat("\n---- glmer", deparse(f, width.cutoff = 500), "| N =", nobs(m), "| logLik =", round(as.numeric(logLik(m)), 3), "\n"); print(round(ct, 6))
  print(VarCorr(m)); msg <- m@optinfo$conv$lme4$messages; if (length(msg)) cat("CONVERGENCE MSG:", msg, "\n")
  out <- list(gendense_lag = c(b = ct[TERM, 1], se = ct[TERM, 2], p = ct[TERM, 4])); attr(out, "n") <- nobs(m); out
}
d0 <- d[complete.cases(d[, vars]), ]
cat("complete-case N:", nrow(d0), "(published 50,216)\n")
full <- fit_re(rhs0, d0)
cat("\nPORTED vs published: gendense_lag b", full[[TERM]][["b"]], "(pub -0.007) se", full[[TERM]][["se"]], "(pub 0.003) N", attr(full, "n"), "(pub 50,216)\n")
ds <- common_sample(d, vars, "ccode1", "ccode2", "year", id = "P9_092")   # named guard exception 0.75 (PI 2026-10-11; see common_b3.R)
M <- run_six(rhs0, ds, fit_re)
if (Sys.getenv("SWEEP_SENS") == "1") {
  source(file.path(CODE, "hist/sens_helper.R"))
  sens_run(id = "P9_092", dat = ds, y = "mid", tie = TERM, ctrl = setdiff(rhs0, TERM), fe = NULL,
           fit = function(dat, add) fit_re(c(rhs0, add), dat)[[TERM]])
}
saveRDS(list(M = M, full = full), file.path(W, "M_P9_092.rds"))
meta <- list(pid = "P9_092", cite = "Millerd 2023", journal = "II", year = 2023, doi = "10.1080/03050629.2023.2277959", family = "TIE")
nt <- sprintf(paste0("PORT: author's Dataverse archive (doi:10.7910/DVN/I9ZSPR) and do-file; Table 2 Model 3 (xtlogit RE by dyad, PR dyads 1989-2014). ",
  "glmer random intercept, nAGQ 12. Full-sample b=%.6f se=%.6f N=%d vs published -0.007 (0.003) N 50,216. Common sample N=%d."),
  full[[TERM]][["b"]], full[[TERM]][["se"]], attr(full, "n"), nrow(ds))
write_p8(meta, TERM, "Table 2 Model 3 / H2a: dyads sharing more densely gender-equal IOs are less likely to engage in MIDs",
  "Table 2 Model 3 (random-effects logit, MID onset, PR dyads 1989-2014)", "xtlogit, re (12-point adaptive quadrature), cumulative MIDs + cubic time",
  "lme4::glmer random intercept by dyad, nAGQ 12, all six steps", -0.007, 0.003, 50216, TRUE, FALSE, M, "Y (cumulative MIDs + cubic time)", nt)
dir.create(file.path(OUT, "results_p9"), recursive = TRUE, showWarnings = FALSE)
src <- file.path(OUT, "results_p8", "P9_092__gendense_lag.csv"); stopifnot(file.exists(src), file.rename(src, file.path(OUT, "results_p9", basename(src))))
sink()
