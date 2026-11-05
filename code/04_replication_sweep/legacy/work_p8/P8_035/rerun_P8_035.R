# P8_035 Bak, Kenwick & Palmer 2016 EJIR "Who's careful: Regime type and target selection" (10.1177/1354066115611479)
# Data + do + log: Dataverse doi:10.7910/DVN/YLKAJH (demselect.dta, bkp_ejir_rep.do/.smcl). Directed dyad-years.
# Headline = primary models, Table (Models 1-4): logit initiate ..., robust cluster(dyadid).
#   Model 1 (Polity): relcap_targ demchal6 relcap_targ_demchal6 cap_1 cap_2 s_wt_glo contigdum dependlow cwpceyrs cwpceyrs2 cwpceyrs3
#   Model 2 (CGV):    relcap_targ cgv_dem relcap_targ_cgv_dem cgv_mil relcap_targ_cgv_mil cgv_civ relcap_targ_cgv_civ + same controls
# KEY TERMS (a priori, abstract; interaction claims -> product terms):
#   relcap_targ_demchal6 (Model 1): "democracies are not more selective than most forms of autocracy"
#   relcap_targ_cgv_mil  (Model 2): "Only military regimes demonstrate unique patterns of target selection"
# Published (authors' log): Model 1 relcap_targ_demchal6 -.2498175 (.2454769), N 918,421; Model 2 relcap_targ_cgv_mil 1.175286 (.4935685), N 1,012,769.
# Estimator all steps: logit, SE clustered by dyadid (Stata-style G/(G-1)).
CODE <- Sys.getenv("SWEEP_CODE"); stopifnot(nzchar(CODE))
source(file.path(CODE, "work_p8/common_b2.R"))
suppressPackageStartupMessages(library(haven))
WD <- file.path(ROOT, "work_p8/P8_035"); W <- file.path(OUT, "work_p8/P8_035"); dir.create(W, recursive = TRUE, showWarnings = FALSE); sink(file.path(W, "log_P8_035.txt"), split = TRUE)
nm <- names(read_dta(file.path(WD, "demselect.dta"), n_max = 1)); print(nm)
cc <- intersect(c("ccode1","ccode2","statea","stateb","ccodea","ccodeb","state1","state2"), nm); cat("ccode cols:", cc, "\n"); stopifnot(length(cc) == 2)
ctl <- c("cap_1","cap_2","s_wt_glo","contigdum","dependlow","cwpceyrs","cwpceyrs2","cwpceyrs3")
m1 <- c("relcap_targ","demchal6","relcap_targ_demchal6", ctl)
m2 <- c("relcap_targ","cgv_dem","relcap_targ_cgv_dem","cgv_mil","relcap_targ_cgv_mil","cgv_civ","relcap_targ_cgv_civ", ctl)
d <- zap_labels(read_dta(file.path(WD, "demselect.dta"), col_select = all_of(unique(c("initiate","dyadid","year", cc, m1, m2))))) %>% as.data.frame()
cat("rows", nrow(d), "years", range(d$year), "\n"); stopifnot(anyDuplicated(d[, c(cc, "year")]) == 0)
meta <- list(pid = "P8_035", cite = "Bak, Kenwick & Palmer 2016", journal = "EJIR", year = 2016, doi = "10.1177/1354066115611479", family = "DEMOCRACY")
ae <- "logit, robust SE clustered by directed dyad (dyadid), cubic peace-years"
eu <- "logit (fixest feglm), SE cluster(dyadid) with G/(G-1), all steps"
for (spec in list(list(rhs = m1, term = "relcap_targ_demchal6", pb = -0.2498175, ps = 0.2454769, pn = 918421, sig = FALSE, hm = "Model 1 (Polity)",
                       ev = "\"democracies are not more selective than most forms of autocracy\""),
                  list(rhs = m2, term = "relcap_targ_cgv_mil", pb = 1.175286, ps = 0.4935685, pn = 1012769, sig = TRUE, hm = "Model 2 (CGV)",
                       ev = "\"Only military regimes demonstrate unique patterns of target selection\""))) {
  fitter <- function(r, data) fit_bin("initiate", r, data, "logit", "~dyadid", terms = spec$term)
  r0 <- fitter(spec$rhs, d)
  ds <- common_sample(d, c("initiate", spec$rhs, "dyadid"), cc[1], cc[2], "year")
  M <- run_six(spec$rhs, ds, fitter)
  # ---- sensitivity hook (split only): runs only under run_sensitivity.R (SWEEP_SENS=1), never in the main sweep ----
  if (Sys.getenv("SWEEP_SENS") == "1" && spec$term == "relcap_targ_cgv_mil") {
    source(file.path(CODE, "hist/sens_helper.R"))
    sens_split(id = "P8_035", dat = ds,
               fit = function(dat, add) { r <- fitter(c(spec$rhs, add), dat)[[spec$term]]; c(b = unname(r[["b"]]), se = unname(r[["se"]]), p = unname(r[["p"]])) })
  }
  write_p8(meta, spec$term, spec$ev, paste(spec$hm, "(authors' log bkp_ejir_rep.smcl; primary models)"), ae, eu, spec$pb, spec$ps, spec$pn, spec$sig, FALSE, M,
    "N (peace-years cubic only)",
    sprintf("Full-sample repro N=%d: b=%.4f (se %.4f). Interaction (regime x target relative capability) product term per KEY_AUDIT rule. Common sample N=%d.",
            attr(r0, "n"), r0[[spec$term]]["b"], r0[[spec$term]]["se"], nrow(ds)),
    verdict = verdict_fn(r0[[spec$term]]["b"], r0[[spec$term]]["p"], spec$pb, spec$sig, FALSE))
  rm(ds, M); gc()
}
sink()
