# P9_064 Anderson & Souva 2010 JCR 54(4) "The Accountability Effects of Political Institutions and Capitalism on Interstate Conflict"
#   (10.1177/0022002710370134). DEMOCRACY family (winning-coalition size, selectorate theory).
# Data + do-files supplied by the PI (staged read-only): "Anderson - Souva JCR Replication Data.dta" -> anderson_souva.dta,
#   "Anderson - Souva JCR replication do file.do" -> anderson_souva.do (DV mzinitt1 = Initiate, the DV the article discusses; matches
#   Tables 1-2), "...for revision DV.do" (DV mzrevis1t1, revisionist initiation; robustness only, not used). No log.
# Directed dyad-years. DV mzinitt1 = "mzinit var advanced one year" (a LEAD) -> history attached at outcome_year = year + 1.
# HEADLINE: Table 1 Model 1 (do-file line 17):
#   logit mzinitt1 w1 powerratio c.w1#c.powerratio cappublic1 c.cappublic1#c.powerratio contig150 lndist allies sglo peaceyrs
#     _spline1 _spline2 _spline3, cluster(ddyadid)
#   KEY TERM w1power = w1 x powerratio (H1b; abstract: "larger winning coalitions are more selective about using military force").
#   Published 1.775*** (0.542), N 497,457 (Table 1; two-tailed, * 10% ** 5% *** 1%).
# SECOND (non-headline) row: Table 2 Model 3 (do-file line 118), key term w1w2 = w1 x w2 (H2, initiator x target coalition size),
#   sample year >= 1970: published -1.299** (0.532), N 796,793.
# Estimator: fixest feglm logit, vcov clustered ddyadid, all six steps.
CODE <- Sys.getenv("SWEEP_CODE"); stopifnot(nzchar(CODE))
source(file.path(CODE, "work_p8/common_b3.R"))
WD <- file.path(ROOT, "work_p9/P9_064"); W <- file.path(OUT, "work_p9/P9_064"); dir.create(W, recursive = TRUE, showWarnings = FALSE)
sink(file.path(W, "log_P9_064.txt"), split = TRUE)
keep <- c("ccode1", "ccode2", "year", "ddyadid", "mzinitt1", "w1", "w2", "w1w2", "powerratio", "cappublic1", "contig150", "lndist",
          "allies", "sglo", "peaceyrs", "_spline1", "_spline2", "_spline3", "w1power", "cap1power")
d <- as.data.frame(haven::zap_labels(haven::read_dta(file.path(WD, "anderson_souva.dta"), col_select = all_of(keep))))
names(d) <- sub("^_spline", "spline", names(d))
cat("rows", nrow(d), "years", range(d$year), "\n")
stopifnot(!anyDuplicated(d[, c("ccode1", "ccode2", "year")]), all(d$ddyadid == 1000 * d$ccode1 + d$ccode2))
# product terms: stored columns equal the Stata factor-variable products
stopifnot(max(abs(d$w1power - d$w1 * d$powerratio), na.rm = TRUE) < 1e-5,
          max(abs(d$cap1power - d$cappublic1 * d$powerratio), na.rm = TRUE) < 1e-5,
          max(abs(d$w1w2 - d$w1 * d$w2), na.rm = TRUE) < 1e-5)
# lead diagnostic: mzinitt1[t] should equal the initiation recorded in t+1; peace years reset at t+1 after an event
o <- d[order(d$ddyadid, d$year), ]; no <- nrow(o)
cs <- o$ddyadid[-1] == o$ddyadid[-no] & o$year[-1] == o$year[-no] + 1
ev <- which(cs & o$mzinitt1[-no] %in% 1)
cat("P(peaceyrs[t+1] == 0 | mzinitt1[t] == 1):", mean(o$peaceyrs[ev + 1] == 0, na.rm = TRUE), " on", length(ev), "events\n"); rm(o, cs, ev, no)
ctrl <- c("contig150", "lndist", "allies", "sglo", "peaceyrs", "spline1", "spline2", "spline3")
specs <- list(
  M1 = list(term = "w1power", rhs = c("w1", "powerratio", "w1power", "cappublic1", "cap1power", ctrl), yr_min = -Inf,
            pub = c(1.775, 0.542, 497457), hm = "Table 1 Model 1 (initiation; W_A x power ratio, capitalism x power ratio; do-file)",
            ev = "abstract: 'robust support for selectorate theory's contention that larger winning coalitions are more selective about using military force'",
            note = "HEADLINE. Published Table 1 Model 1 1.775*** (0.542), N 497,457."),
  M3 = list(term = "w1w2", rhs = c("w1", "powerratio", "w1power", "w2", "w1w2", ctrl), yr_min = 1970,
            pub = c(-1.299, 0.532, 796793), hm = "Table 2 Model 3 (initiation and targeting, year >= 1970; W_A x W_B; do-file)",
            ev = "text: 'the coefficient on the interaction variable for state A and state B's winning coalition sizes will be negative and significant'",
            note = "Non-headline second row (H2, dyadic coalition-size interaction). Published Table 2 Model 3 -1.299** (0.532), N 796,793."))
d$outcome_year <- d$year + 1
meta <- list(pid = "P9_064", cite = "Anderson & Souva 2010", journal = "JCR", year = 2010, doi = "10.1177/0022002710370134", family = "DEMOCRACY")
dir.create(file.path(OUT, "results_p9"), recursive = TRUE, showWarnings = FALSE)
res <- list()
for (s in names(specs)) {
  sp <- specs[[s]]; vars <- c("mzinitt1", sp$rhs, "ddyadid")
  fitter <- function(r, data) fit_bin("mzinitt1", r, data, "logit", "~ddyadid", terms = sp$term)
  dd <- d[d$year >= sp$yr_min, ]
  cat("\n==================== spec", s, "key term", sp$term, "\n")
  full <- fitter(sp$rhs, dd[complete.cases(dd[, vars]), ])
  cat("\nFULL-SAMPLE:", sp$term, "b", full[[sp$term]][["b"]], "se", full[[sp$term]][["se"]], "N", attr(full, "n"), "| published", sp$pub, "\n")
  ds <- common_sample(dd, vars, "ccode1", "ccode2", "outcome_year"); rm(dd); invisible(gc())
  M <- run_six(sp$rhs, ds, fitter)
  res[[s]] <- list(M = M, full = full, n_used = nrow(ds))
  if (s == "M1" && Sys.getenv("SWEEP_SENS") == "1") {   # sensitivity hook (split: product term); runs only under run_sensitivity.R
    source(file.path(CODE, "hist/sens_helper.R"))
    sens_split(id = "P9_064", dat = ds, fit = function(dat, add) fitter(c(sp$rhs, add), dat)[[sp$term]])   # product term -> split only
  }
  write_p8(meta, sp$term, sp$ev, sp$hm, "logit, SE clustered on directed dyad, peace years + cubic splines",
    "feglm logit, vcov clustered ddyadid, all six steps", sp$pub[1], sp$pub[2], sp$pub[3], TRUE, FALSE, M,
    "N (peace years + cubic splines only)",
    sprintf(paste0("Directed dyad-years (author's data + do-file; DV MID initiation advanced one year). %s Full-sample b=%.6f se=%.6f N=%d. ",
      "Common 1822-2011 sample N=%d (history attached at outcome year t+1). Interaction: constituent terms kept in all steps."),
      sp$note, full[[sp$term]][["b"]], full[[sp$term]][["se"]], attr(full, "n"), nrow(ds)))
  src <- file.path(OUT, "results_p8", paste0("P9_064__", sp$term, ".csv")); stopifnot(file.exists(src), file.rename(src, file.path(OUT, "results_p9", basename(src))))
  rm(ds, M, full); invisible(gc())
}
saveRDS(res, file.path(W, "M_P9_064.rds"))
sink()
