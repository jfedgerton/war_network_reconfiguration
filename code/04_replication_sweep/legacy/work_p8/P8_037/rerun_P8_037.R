# P8_037 Goldsmith, Semenovich, Sowmya & Grgic 2017 WP "Political Competition and the Initiation of International Conflict" (10.1017/S0043887116000307)
# Data: authors' Dataverse (doi:10.7910/DVN/J06MRI) Main_Analysis/BGetal_WP2017_dataset_main.dta (slimmed to model columns, 1950-2005).
# Headline (authors' Table_2.do, model 3 "components"): directed dyads 1950-2005, DV = ICB crisis initiation t+1 (ficbinitdum), probit,
#   dyad-clustered SE, cubic peace years; RHS = CompetitionA, CompetitionA x RegimeTypeB, ConstraintsA (+x), ParticipationA (+x), RegimeTypeB, controls.
# KEY TERM (a priori, title/abstract: 'institutionalized political competition' of initiator conditioned on target regime): the product
#   BGcomp1pXpolity22p (CompetitionA x RegimeTypeB); constituent BGcomp1p reported in the log only.
# DV is ICB crisis initiation in t+1 -> history keyed to year + 1 (GML MID history as elsewhere).
CODE <- Sys.getenv("SWEEP_CODE"); stopifnot(nzchar(CODE))
source(file.path(CODE, "work_p8/p8_common.R"))
WD <- file.path(ROOT, "work_p8/P8_037"); W <- file.path(OUT, "work_p8/P8_037"); dir.create(W, recursive = TRUE, showWarnings = FALSE); sink(file.path(W, "log_P8_037.txt"), split = TRUE)
d <- readRDS(file.path(WD, "goldsmith_slim.rds"))
stopifnot(!anyDuplicated(d[, c("ccode1", "ccode2", "year")]))
d$yr_dv <- d$year + 1
rhs <- c("BGcomp1p", "BGcomp1pXpolity22p", "BGexconst1p", "BGexconst1pXpolity22p", "BGpart1p", "BGpart1pXpolity22p", "polity22p",
         "lnGDPshare_1", "IntraUCDP_1", "IntraUCDP_2", "ally_dum_w", "NUMIGOfill", "cap_1", "cap_2", "contig24", "lndistance",
         "icbinitpcyrs", "icbinitpcyrs2", "icbinitpcyrs3")
vars <- c("ficbinitdum", rhs)
fitter <- function(r, data) fit_bin("ficbinitdum", r, data, "probit", "~dyadid", c("BGcomp1pXpolity22p", "BGcomp1p"))
d0 <- d[complete.cases(d[, vars]), ]; cat("complete-case N:", nrow(d0), "\n")
full <- fitter(rhs, d0)
ds <- common_sample(d, vars, "ccode1", "ccode2", "yr_dv")
M <- run_six(rhs, ds, fitter)
# ---- sensitivity hook (split only): runs only under run_sensitivity.R (SWEEP_SENS=1), never in the main sweep ----
if (Sys.getenv("SWEEP_SENS") == "1") {
  source(file.path(CODE, "hist/sens_helper.R"))
  sens_split(id = "P8_037", dat = ds,
             fit = function(dat, add) { r <- fitter(c(rhs, add), dat)$BGcomp1pXpolity22p; c(b = unname(r[["b"]]), se = unname(r[["se"]]), p = unname(r[["p"]])) })
}
saveRDS(list(M = M, full = full), file.path(W, "M_P8_037.rds"))
# ---- port: result-CSV write step for this study, copied verbatim from work_p8/write_pending.R (paths: saved model list read from OUT) ----
{

  o <- readRDS(file.path(W, "M_P8_037.rds"))
  pub <- as.numeric(Sys.getenv("PUB037_B", NA)); pse <- as.numeric(Sys.getenv("PUB037_SE", NA)); pn <- as.numeric(Sys.getenv("PUB037_N", NA))
  meta <- list(pid = "P8_037", short_cite = "Goldsmith et al. 2017", journal = "WP", year = 2017, doi = "10.1017/S0043887116000307",
    family = "DEMOCRACY", headline_table_model = "Table 2 (authors' Table_2.do model 3: regime components with controls)",
    article_estimator = "probit, dyad-clustered SE, cubic peace years (authors' Table_2.do)",
    est_used_all_steps = "probit, dyad-clustered SE (fixest), all six steps", author_has_history_control = "N (peace years only)")
  write_p8(meta, "BGcomp1pXpolity22p", "Title/abstract: institutionalized political competition of initiator, conditioned on target regime (CompetitionA x RegimeTypeB)",
    pub, pse, pn, NA, FALSE, o$M,
    paste0("Authors' data and do-file spec reproduced exactly as coded (N ", attr(o$full, "n"), "). DV = ICB crisis initiation t+1 (not MIDs); history keyed to t+1. ",
           if (is.na(pub)) "Published coefficient not available (article PDF not accessible) -> NO_PUBLISHED. " else "",
           "Constituent CompetitionA (BGcomp1p) m0 b = ", signif(o$M$m0$BGcomp1p["b"], 3), ", mDM b = ", signif(o$M$mDM$BGcomp1p["b"], 3), "."), full = o$full)
}
sink()
