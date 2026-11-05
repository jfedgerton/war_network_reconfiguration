# P8_050 Mousseau 2019 PEPS "Four Ways We Know the Democratic Peace Correlation Does Not Exist..." (10.1515/peps-2019-0023)
# Data: Dataverse doi:10.7910/DVN/YJD6CA (Replication Data.dta). PDF: author copy (UCF site), Table 3.
# Headline = Table 3 Model 3 (Model 2: Both Contractualist perfectly predicts -> dropped, N 442,840): logit fMID on
#   Democracy (dml), Both contractualist or axial (baxecon), One or both axial (obaxe), capability ratio (lncprt), major power (mjpw),
#   contiguity (cntg), distance (dist), system size (sysze), cubic spline temporal corrections (fpceyrs + fspl1-3), SE clustered by dyad.
#   Published: dml -0.01 (0.02); baxecon -3.28 (0.85); N 536,962; clusters 14,035. 1921-2010, non-directed dyads, IVs lagged 1 yr.
# KEY TERMS (a priori): baxecon -- "Both Contractualist or Axial (-3.28) is negative and highly significant"; dml -- abstract:
#   "democracy without contractualist economy has a near-zero correlation with peace".
CODE <- Sys.getenv("SWEEP_CODE"); stopifnot(nzchar(CODE))
source(file.path(CODE, "work_p8/common_b2.R"))
suppressPackageStartupMessages(library(haven))
WD <- file.path(ROOT, "work_p8/P8_050"); W <- file.path(OUT, "work_p8/P8_050"); dir.create(W, recursive = TRUE, showWarnings = FALSE); sink(file.path(W, "log_P8_050.txt"), split = TRUE)
d <- zap_labels(read_dta(file.path(WD, "rep.dta"))) %>% as.data.frame()
stopifnot(anyDuplicated(d[, c("dcode", "year")]) == 0); cat("years", range(d$year), " ccode1<ccode2:", mean(d$ccode1 < d$ccode2), "\n")
rhs <- c("dml","baxecon","obaxe","lncprt","mjpw","cntg","dist","sysze","fpceyrs","fspl1","fspl2","fspl3")
terms <- c("baxecon", "dml")
fitter <- function(r, data) fit_bin("fMID", r, data, "logit", "~dcode", terms = terms)
r0 <- fitter(rhs, d); n_orig <- attr(r0, "n")
ds <- common_sample(d, c("fMID", rhs, "dcode"), "ccode1", "ccode2", "year")
M <- run_six(rhs, ds, fitter)
# ---- sensitivity hook: runs only under run_sensitivity.R (SWEEP_SENS=1), never in the main sweep ----
if (Sys.getenv("SWEEP_SENS") == "1") {
  source(file.path(CODE, "hist/sens_helper.R"))
  sens_run(id = "P8_050", dat = ds, y = "fMID", tie = "baxecon", ctrl = setdiff(rhs, "baxecon"), fe = NULL,
           fit = function(dat, add) fitter(c(rhs, add), dat)[["baxecon"]])
}
meta <- list(pid = "P8_050", cite = "Mousseau 2019", journal = "PEPS", year = 2019, doi = "10.1515/peps-2019-0023", family = "CAPITALISM")
est <- "logit, SE clustered by non-directed dyad (dcode), cubic-spline peace-years; all steps (fixest feglm, Stata-style G/(G-1))"
nt <- sprintf("Full-sample repro (N=%d; pub 536,962): baxecon %.3f (%.3f), dml %.4f (%.4f). Common 1822-2011 sample N=%d. Rare-events logit reported by author as identical -> ordinary logit headline.",
  n_orig, r0$baxecon["b"], r0$baxecon["se"], r0$dml["b"], r0$dml["se"], nrow(ds))
write_p8(meta, "baxecon", "\"Both Contractualist or Axial (-3.28) is negative and highly significant\"", "Table 3 Model 3",
  "logit, dyad-clustered SE, cubic splines", est, -3.28, 0.85, 536962, TRUE, FALSE, M, "N (peace-years splines only)",
  paste(nt, "Key CIE term = both contractualist-or-axial (Model 2's both-contractualist term is perfectly predicting)."),
  verdict = verdict_fn(r0$baxecon["b"], r0$baxecon["p"], -3.28, TRUE, FALSE))
write_p8(meta, "dml", "\"democracy without contractualist economy has a near-zero correlation with peace\"", "Table 3 Model 3",
  "logit, dyad-clustered SE, cubic splines", est, -0.01, 0.02, 536962, FALSE, FALSE, M, "N (peace-years splines only)",
  paste(nt, "Democracy = lower Polity2 in dyad; paper claims null. Published b rounded to 2 dp (-0.01) so % comparison is coarse."),
  verdict = verdict_fn(r0$dml["b"], r0$dml["p"], -0.01, FALSE, FALSE))
sink()
