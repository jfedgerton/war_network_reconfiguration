# P9_048 Bareis 2024 JCR 68(7-8) "Interstate Resource Conflicts: A Network-Centric Resource Access Security Perspective"
#   (10.1177/00220027231181051). TIE family (state A's position in the oil-import supply network).
# Data + do-file supplied by the PI (SAGE replication archive; staged read-only): replication_code.do, replication_data_dyadic.dta.
#   Table 1 needs only the dyadic file (the *_all, no_ongo and only_orig variants are used only in appendix robustness tables). No log.
# Directed politically relevant dyad-years 1962-2010. DV Conflict_initiation (state A initiates a MID in t; contemporaneous; all IVs
#   lagged one year by construction) -> history attached at Year.
# HEADLINE: Table 1 Model 2 (do-file line 24; the column variable FE is a copy of Conflict_initiation, asserted below):
#   logit FE c.Outside_oil_dependence##c.Degree_centrality Net_oil_imports $control i.ccode1 i.ccode2, vce(cluster group)
#   KEY TERM oil_dep_x_degree = Outside_oil_dependence x Degree_centrality (abstract: diversified position in resource supply
#   networks mitigates conflict propensity of outside-dependent states). Published -3.891*** (1.495), N 75,868.
# Estimator: logit with state-A and state-B dummies = fixest feglm logit with ccode1 + ccode2 fixed effects (Stata drops
#   observations whose dummy perfectly predicts the outcome; fixest drops the same all-zero FE groups), vcov clustered on group
#   (= group(ccode1 ccode2)), all six steps.
CODE <- Sys.getenv("SWEEP_CODE"); stopifnot(nzchar(CODE))
source(file.path(CODE, "work_p8/common_b3.R"))
WD <- file.path(ROOT, "work_p9/P9_048"); W <- file.path(OUT, "work_p9/P9_048"); dir.create(W, recursive = TRUE, showWarnings = FALSE)
sink(file.path(W, "log_P9_048.txt"), split = TRUE)
ctrl <- c("Territorial_contiguity", "Distance_log", "Colonial_contiguity", "Democracy_both", "Capability_ratio_log", "Alliance",
          "Trade_dependence", "Shared_IGO_memberships", "Civil_conflict_one", "Civil_conflict_both", "Peace_years", "Spline_1", "Spline_2", "Spline_3")
d <- as.data.frame(haven::zap_labels(haven::read_dta(file.path(WD, "replication_data_dyadic.dta"),
  col_select = all_of(c("ccode1", "ccode2", "group", "Year", "Conflict_initiation", "FE", "Outside_oil_dependence", "Degree_centrality",
                        "Net_oil_imports", ctrl)))))
cat("rows", nrow(d), "years", range(d$Year), "\n")
stopifnot(nrow(d) == 128325, !anyDuplicated(d[, c("ccode1", "ccode2", "Year")]), identical(d$FE, d$Conflict_initiation))
stopifnot(all(tapply(d$group, paste(d$ccode1, d$ccode2), function(x) length(unique(x))) == 1))
TERM <- "oil_dep_x_degree"
d[[TERM]] <- d$Outside_oil_dependence * d$Degree_centrality          # c.Outside_oil_dependence#c.Degree_centrality
rhs0 <- c("Outside_oil_dependence", "Degree_centrality", TERM, "Net_oil_imports", ctrl)
vars <- c("Conflict_initiation", rhs0, "group", "ccode1", "ccode2")
fitter <- function(r, data) fit_bin("Conflict_initiation", r, data, "logit", "~group", fe = "ccode1 + ccode2", terms = TERM)
full <- fitter(rhs0, d[complete.cases(d[, vars]), ])
pub <- c(b = -3.891, se = 1.495, n = 75868)
cat("\nFULL-SAMPLE:", TERM, "b", full[[TERM]][["b"]], "se", full[[TERM]][["se"]], "N", attr(full, "n"), "| published", pub, "\n")
ds <- common_sample(d, vars, "ccode1", "ccode2", "Year")
M <- run_six(rhs0, ds, fitter)
saveRDS(list(M = M, full = full, n_used = nrow(ds)), file.path(W, "M_P9_048.rds"))

# ---- sensitivity hook: runs only under run_sensitivity.R (SWEEP_SENS=1), never in the main sweep (split: product term) ----
if (Sys.getenv("SWEEP_SENS") == "1") {
  source(file.path(CODE, "hist/sens_helper.R"))
  sens_split(id = "P9_048", dat = ds, fit = function(dat, add) fitter(c(rhs0, add), dat)[[TERM]])   # product term -> split only
}

meta <- list(pid = "P9_048", cite = "Bareis 2024", journal = "JCR", year = 2024, doi = "10.1177/00220027231181051", family = "TIE")
write_p8(meta, TERM, "abstract: 'a diversified position in resource supply networks mitigates the conflict propensity of states characterized by outside dependence'",
  "Table 1 Model 2 (FE logit; outside oil dependence x degree centrality; do-file)",
  "logit with state-A and state-B dummies, SE clustered on dyad, peace years + cubic splines",
  "feglm logit, ccode1 + ccode2 FE, vcov clustered group, all six steps", pub[["b"]], pub[["se"]], pub[["n"]], TRUE, FALSE, M,
  "N (peace years + cubic splines only)",
  sprintf(paste0("Directed politically relevant dyad-years 1962-2010 (author's data + do-file; IVs lagged one year by construction). ",
    "Published Table 1 Model 2 -3.891*** (1.495), N 75,868. Full-sample b=%.6f se=%.6f N=%d. Common 1822-2011 sample N=%d (history at Year). ",
    "Interaction: constituent terms kept in all steps."), full[[TERM]][["b"]], full[[TERM]][["se"]], attr(full, "n"), nrow(ds)))
dir.create(file.path(OUT, "results_p9"), recursive = TRUE, showWarnings = FALSE)
src <- file.path(OUT, "results_p8", paste0("P9_048__", TERM, ".csv")); stopifnot(file.exists(src), file.rename(src, file.path(OUT, "results_p9", basename(src))))
sink()
