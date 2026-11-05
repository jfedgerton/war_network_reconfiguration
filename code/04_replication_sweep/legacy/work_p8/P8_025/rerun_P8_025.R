# P8_025 Souva 2004 II "Institutional Similarity and Interstate Conflict" (10.1080/03050620490492213)
# Data: author's replication zip (myweb.fsu.edu/msouva/ii 2004 inst sim and conflict.zip -> "Souva ii data stata 7 final version.dta").
# Paper (author PDF, myweb.fsu.edu/msouva): Table 2 Model B (all MIDs, 1982-1996), N = 61,338.
#   Estimator (p. 272): GEE, logit link, AR(1) working correlation, robust SE; IVs lagged one year (DV = mid3t1 = MID onset in t+1).
#   Published: econ sim -4.63 (1.67), pol sim -5.17 (2.02), econ*pol 6.56 (2.53).
# KEY TERMS (a priori, H1-H3 p. 267): H1 political similarity (polsim), H2 economic similarity (icrgsim), H3 joint similarity (policrg product).
# History keyed to the outcome year (year + 1), since the DV is the t+1 onset.
CODE <- Sys.getenv("SWEEP_CODE"); stopifnot(nzchar(CODE))
source(file.path(CODE, "work_p8/p8_common.R"))
suppressMessages(library(haven))
WD <- file.path(ROOT, "work_p8/P8_025"); W <- file.path(OUT, "work_p8/P8_025"); dir.create(W, recursive = TRUE, showWarnings = FALSE); sink(file.path(W, "log_P8_025.txt"), split = TRUE)
d <- read_dta(file.path(WD, "Souva ii data stata 7 final version.dta")) %>% zap_labels() %>% as.data.frame()
stopifnot(!anyDuplicated(d[, c("ccode1", "ccode2", "year")]))
d$yr_dv <- d$year + 1
rhs <- c("icrgsim", "polsim", "policrg", "demweak", "devweak", "demdev", "depweak", "lnrelgdp", "contig150", "lndist", "majdyad")
vars <- c("mid3t1", rhs, "dyadid", "year")
cat("Complete-case N (paper Model B N = 61,338):", sum(complete.cases(d[, vars])), "\n")
terms <- c("icrgsim", "polsim", "policrg")
fitter <- function(r, data) fit_gee_logit_ar1("mid3t1", r, data, "dyadid", "year", terms)

# full-sample reproduction check
d0 <- d[complete.cases(d[, vars]), ]
full <- fitter(rhs, d0)

ds <- common_sample(d, vars, "ccode1", "ccode2", "yr_dv")
M <- run_six(rhs, ds, fitter)
# ---- sensitivity hook (split only): runs only under run_sensitivity.R (SWEEP_SENS=1), never in the main sweep ----
if (Sys.getenv("SWEEP_SENS") == "1") {
  source(file.path(CODE, "hist/sens_helper.R"))
  sens_split(id = "P8_025", dat = ds,
             fit = function(dat, add) { r <- fitter(c(rhs, add), dat)$polsim; c(b = unname(r[["b"]]), se = unname(r[["se"]]), p = unname(r[["p"]])) })
}
meta <- list(pid = "P8_025", short_cite = "Souva 2004", journal = "II", year = 2004, doi = "10.1080/03050620490492213",
  family = "DEMOCRACY", headline_table_model = "Table 2, Model B (all MIDs 1982-96)",
  article_estimator = "GEE logit, AR(1) working correlation, robust SE (p. 272)",
  est_used_all_steps = "GEE logit AR(1) robust (own estimating-equation solver, actual-year AR spacing = Stata force) in all six steps",
  author_has_history_control = "N")
ev <- c(icrgsim = "H2: 'Dyads with similar economic institutions are less likely to experience militarized conflict'",
        polsim = "H1: 'Dyads with similar political institutions are less likely to experience militarized conflict'",
        policrg = "H3: 'similar political and economic institutions are less likely to experience militarized conflict'")
pub <- list(icrgsim = c(-4.63, 1.67), polsim = c(-5.17, 2.02), policrg = c(6.56, 2.53))
for (t in terms) {
  fb <- full[[t]]
  write_p8(meta, t, ev[[t]], pub[[t]][1], pub[[t]][2], 61338, NA, FALSE, M,
    paste0("Author data + author spec. Full-sample (no history restriction) b = ", round(fb["b"], 3), " (SE ", round(fb["se"], 3),
           "). Paper reports one-tailed tests; verdict uses two-sided .05 from pub b/SE. Stata xtgee AR(1) rho estimator may differ slightly from the moment estimator used here. ",
           "History keyed to outcome year t+1. H3 term is the product; constituents interpreted at 0 of the other."))
}
sink()
