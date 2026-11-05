# P9_017 Anderson, Mitchell & Schilling 2016 II 42(4) "Kantian Dynamics Revisited: Time-Varying Analyses of Dyadic IGO-Conflict
#   Relationships" (10.1080/03050629.2016.1185714)
# Data + code: Harvard Dataverse doi:10.7910/DVN/XN7EFD ("table1logitmodeldataset.dta" = Dorussen & Ward 2008 data; "table1dofile.do";
#   TVEM SAS programs). Staged as amsm_table1logit.dta / amsm_table1.do / amsm_*.sas. No author log.
# The article's main models are TVEM varying-coefficient logits (Penn State SAS macro, P-splines in year); not portable faithfully.
# PORTED = the constant-coefficient baseline the paper reports, Table 1 "Static Logit Model of MID Occurrence, 1950-2000":
#   "logit force_y1 i2_numigo low_depend_ksg2 low_polity2 noncontiguity force_pceyrs force_spline1 force_spline2 force_spline3"
#   (default, non-robust SE). KEY TERM i2_numigo (number of shared IGO memberships).
# Published Table 1 (static logit): i2_numigo b .0322, SE .0026, N 469,517; the full-sample re-run matches. No log in the archive.
# Undirected dyad-years; dyadid = ccode1 * 1000 + ccode2 (checked against country names). DV force_y1 as in the author's data.
# Estimator: fixest feglm logit, iid vcov, all six steps.
CODE <- Sys.getenv("SWEEP_CODE"); stopifnot(nzchar(CODE))
source(file.path(CODE, "work_p8/common_b3.R"))
WD <- file.path(ROOT, "work_p9/P9_017"); W <- file.path(OUT, "work_p9/P9_017"); dir.create(W, recursive = TRUE, showWarnings = FALSE)
sink(file.path(W, "log_P9_017.txt"), split = TRUE)
d <- as.data.frame(haven::zap_labels(haven::read_dta(file.path(WD, "amsm_table1logit.dta"))))
stopifnot(nrow(d) == 604946)
cat("rows with missing dyadid:", sum(is.na(d$dyadid)), "\n")
d$ccode1 <- d$dyadid %/% 1000; d$ccode2 <- d$dyadid %% 1000
k <- !is.na(d$dyadid)
stopifnot(!anyDuplicated(d[k, c("dyadid", "year")]), all(d$ccode1[k] < d$ccode2[k]),
          all(d$ccode1[k & d$country1 == "USA"] == 2), all(d$ccode2[k & d$country2 == "Canada"] == 20))
cat("year range:", range(d$year), "\n")
term <- "i2_numigo"
rhs0 <- c(term, "low_depend_ksg2", "low_polity2", "noncontiguity", "force_pceyrs", "force_spline1", "force_spline2", "force_spline3")
vars <- c("force_y1", rhs0)
fitter <- function(r, data) fit_bin("force_y1", r, data, "logit", "iid", terms = term)
d0 <- d[complete.cases(d[, vars]), ]
stopifnot(!anyNA(d0$dyadid))
full <- fitter(rhs0, d0)
cat("\nFULL-SAMPLE (published Table 1: b .0322, SE .0026, N 469,517):", term, "b", full[[term]][["b"]], "se", full[[term]][["se"]], "N", attr(full, "n"), "\n")
stopifnot(attr(full, "n") == 469517)
ds <- common_sample(d, vars, "ccode1", "ccode2", "year")
M <- run_six(rhs0, ds, fitter)
saveRDS(list(M = M, full = full, n_used = nrow(ds)), file.path(W, "M_P9_017.rds"))

# ---- sensitivity hook: runs only under run_sensitivity.R (SWEEP_SENS=1), never in the main sweep ----
if (Sys.getenv("SWEEP_SENS") == "1") {
  source(file.path(CODE, "hist/sens_helper.R"))
  # i2_numigo enters the static logit additively -> full sensitivity
  sens_run(id = "P9_017", dat = ds, y = "force_y1", tie = term, ctrl = setdiff(rhs0, term), fe = NULL,
           fit = function(dat, add) fitter(c(rhs0, add), dat)[[term]])
}

meta <- list(pid = "P9_017", cite = "Anderson, Mitchell & Schilling 2016", journal = "II", year = 2016,
             doi = "10.1080/03050629.2016.1185714", family = "TIE")
write_p8(meta, term,
  "title/abstract: time-varying analyses of dyadic IGO-conflict relationships (shared IGO memberships)",
  "Table 1 static logit of MID occurrence 1950-2000 (constant-coefficient baseline; main TVEM models not ported)",
  "logit, default SE, peace years + 3 splines (Table 1); main models TVEM varying-coefficient logit (SAS)",
  "feglm logit, iid vcov, all six steps", 0.0322, 0.0026, 469517, TRUE, FALSE, M, "N (peace years + 3 splines only)",
  sprintf(paste0("Undirected dyad-years 1950-2000 (Dorussen & Ward 2008 data). Published values from article Table 1 static logit (b .0322, SE .0026, N 469,517; no log). ",
    "Full-sample b=%.6f se=%.6f N=%d. Common 1822-2011 sample N=%d. Constant-coefficient baseline only; the article's ",
    "headline TVEM (time-varying) estimates are not reproduced."), full[[term]][["b"]], full[[term]][["se"]], attr(full, "n"), nrow(ds)))
dir.create(file.path(OUT, "results_p9"), recursive = TRUE, showWarnings = FALSE)
src <- file.path(OUT, "results_p8", paste0("P9_017__", term, ".csv")); dst <- file.path(OUT, "results_p9", paste0("P9_017__", term, ".csv"))
stopifnot(file.exists(src), file.rename(src, dst))
sink()
