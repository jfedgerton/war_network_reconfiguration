# P8_019 Peceny, Beer & Sanchez-Terry 2002 APSR "Dictatorial Peace?" (10.1017/S0003055402004203)
# Data: Peceny et al.'s own dyad-year file as passed to Reiter & Stam (Dataverse doi:10.7910/DVN/FLJMRM, doublereduce.dta = the Peceny
#   data "doubled"); the Reiter-Stam log (dictpeacelog6-26) states that 'logit ... if dummy==1' is "a successful replication of Peceny et al,
#   p. 24, Table 3, Model 4" (N = 376,771). dummy==1 selects the original (undirected) row of each doubled pair.
# Headline: Table 3, Model 4: logit dispute on persdem personal military single democ contig majpow ally loglsrat advanced dispyrs + 3 splines.
# KEY TERMS (a priori, abstract/title: homogeneous autocratic dyads and the dictatorial peace): single (joint single-party), military (joint
#   military), personal (joint personalist), persdem (personalist-democracy mixed dyad), democ (joint democracy).
CODE <- Sys.getenv("SWEEP_CODE"); stopifnot(nzchar(CODE))
source(file.path(CODE, "work_p8/p8_common.R"))
suppressMessages(library(haven))
WD <- file.path(ROOT, "work_p8/P8_019"); W <- file.path(OUT, "work_p8/P8_019"); dir.create(W, recursive = TRUE, showWarnings = FALSE); sink(file.path(W, "log_P8_019.txt"), split = TRUE)
rhs <- c("persdem", "personal", "military", "single", "democ", "contig", "majpow", "ally", "loglsrat", "advanced", "dispyrs",
         "dspline1", "dspline2", "dspline3")
d <- read_dta(file.path(WD, "doublereduce.dta"), col_select = all_of(c("idyr", "statea", "stateb", "year", "dispute", rhs))) %>%
  zap_labels() %>% as.data.frame()
# Reiter-Stam log: 'dummy==1' (alternating rows after sort idyr; not reproducible in R because Stata breaks idyr ties arbitrarily and
#   idyr is not always paired) and 'if statea>stateb' (N 376,728) both reproduce Peceny et al. (R&S: "OK, that worked"). We use statea>stateb.
d <- d[!is.na(d$statea) & !is.na(d$stateb) & d$statea > d$stateb, ]
cat("rows:", nrow(d), " years:", range(d$year, na.rm = TRUE), " NA years:", sum(is.na(d$year)), "\n")
d <- d[!is.na(d$year), ]
cat("duplicate statea-stateb-year rows:", sum(duplicated(d[, c("statea", "stateb", "year")])), "\n")
terms <- c("single", "military", "personal", "persdem", "democ")
fitter <- function(r, data) fit_bin("dispute", r, data, "logit", "iid", terms)
vars <- c("dispute", rhs)
d0 <- d[complete.cases(d[, vars]), ]; cat("complete-case N:", nrow(d0), "(Reiter-Stam log 376,771)\n")
full <- fitter(rhs, d0)
ds <- common_sample(d, vars, "statea", "stateb", "year")
M <- run_six(rhs, ds, fitter)

# ---- sensitivity hook: runs only under run_sensitivity.R (SWEEP_SENS=1), never in the main sweep ----
if (Sys.getenv("SWEEP_SENS") == "1") {
  source(file.path(CODE, "hist/sens_helper.R"))
  # [2026-11-10 audit2] hook moved to the new headline term (single)
  sens_run(id = "P8_019", dat = ds, y = "dispute", tie = "single", ctrl = setdiff(rhs, "single"), fe = NULL,
           fit = function(dat, add) fitter(c(rhs, add), dat)[["single"]])
}

saveRDS(list(M = M, full = full), file.path(W, "M_P8_019.rds"))
# ---- port: result-CSV write step for this study, copied verbatim from work_p8/write_pending.R (paths: saved model list read from OUT) ----
{

  o <- readRDS(file.path(W, "M_P8_019.rds"))
  meta <- list(pid = "P8_019", short_cite = "Peceny et al. 2002", journal = "APSR", year = 2002, doi = "10.1017/S0003055402004203",
    family = "DEMOCRACY", headline_table_model = "Table 3, Model 4 (per Reiter & Stam replication log)",
    article_estimator = "logit, peace years + 3 splines (Stata default SE per Reiter & Stam replication of Table 3 Model 4)",
    est_used_all_steps = "logit iid SE (fixest), all six steps", author_has_history_control = "N (peace years/splines only)")
  # Reiter & Stam log (dummy==1 run, N 376,771) -- described by R&S as reproducing Peceny et al. Table 3 Model 4; article PDF not accessible
  rs <- list(single = c(-0.8406913, 0.1303307), military = c(-0.7936331, 0.7214446), personal = c(0.2219903, 0.2342579),
             persdem = c(0.613281, 0.0994073), democ = c(-1.052696, 0.1285853))
  ev <- c(single = "Abstract: 'single-party states are more peaceful toward one another than are mixed dyads'",
          military = "Abstract: military-regime dyads 'were not less likely to engage in militarized interstate disputes'",
          personal = "Abstract: personalist dyads 'were not less likely to engage in militarized interstate disputes'",
          persdem = "Autocratic subtype x democracy mixed dyad (personalist-democracy) in the headline model",
          democ = "Joint democracy benchmark in the dictatorial-peace model")
  for (t in names(rs)) write_p8(meta, t, ev[[t]], rs[[t]][1], rs[[t]][2], 376771, NA, FALSE, o$M,
    "Peceny et al.'s own data (via Reiter & Stam Dataverse FLJMRM). 'Published' b/SE are the Reiter-Stam log values for the run they state reproduces Peceny et al. Table 3 Model 4 (article PDF not accessible here) -- treat verdict as verified against that log, not the article. Sample: statea>stateb (R&S alternative that also reproduces).",
    full = o$full)
}
sink()
