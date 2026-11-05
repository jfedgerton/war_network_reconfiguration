# P8_042 Boldt 2026 II "Democratic transitions and international conflict" (10.1080/03050629.2025.2609092)
# Data + code: Harvard Dataverse doi:10.7910/DVN/NSHORX (dem-transitions-replication-data.rds, pooled-mids.R, log; PanelMatch 3.1.1).
# HEADLINE (abstract: "Pooling across three- and five-year periods after transition, democratic transitions reduce conflict, including fatal
#   conflicts"; Figure 6 = pooled-mids.R): PanelMatch TSCS matching + DiD, treatment tr = 1[v2x_regime >= 2], lag = 5, CBPS weighting,
#   exact match on UN region, covariates lag 1:5 of the outcome, ln_cinc, ln_pecpp, terr_dispute; ATT pooled over leads 1:3;
#   bootstrap SE (10,000 iterations, seed 1). MONADIC country-years (COW codes), 1900-2014.
# KEY TERMS (a priori): pooled 3-year ATT on MID participation (att_mid_pool3) and on fatal MIDs (att_fatal_pool3).
# Estimator: PanelMatch 3.1.3 (CRAN; author used 3.1.1) with the authors' exact call in all six steps; history/benchmarks are added to
#   covs.formula as time-t values (they summarize t-1..t-5 / t-2..t-6, i.e. pre-treatment). MONADIC convention: mD = NA; mDM = +own w15;
#   m26 = +own w26; bC = +own CINC mean; bD = +own polity>=6 mean. Common sample: years <= 2011 (GML coverage), all steps.
# p-values: two-sided normal approximation from the bootstrap SE (authors report 90% percentile CIs).
CODE <- Sys.getenv("SWEEP_CODE"); stopifnot(nzchar(CODE))
source(file.path(CODE, "work_p8/common_b3.R"))
suppressPackageStartupMessages(library(PanelMatch))
args <- commandArgs(TRUE); NB <- if (length(args)) as.integer(args[1]) else 10000L
WD <- file.path(ROOT, "work_p8/P8_042"); W <- file.path(OUT, "work_p8/P8_042"); dir.create(W, recursive = TRUE, showWarnings = FALSE); sink(file.path(W, paste0("log_P8_042_nb", NB, ".txt")), split = TRUE)
m <- as.data.frame(readRDS(file.path(WD, "dem-transitions-replication-data.rds")))
m$country <- as.integer(m$country); m$year <- as.integer(m$year)
m$tr <- as.numeric(ifelse(m$v2x_regime >= 2, 1, 0))
for (v in c("mid", "fatal", "ln_cinc", "ln_pecpp", "terr_dispute")) m[[v]] <- as.numeric(m[[v]])
m <- distinct(m, country, year, .keep_all = TRUE)
m <- add_history_monadic(m, "country", "year")
cat("rows", nrow(m), " NA bm_cinc:", sum(!is.finite(m$bm_cinc_w15)), " NA bm_democ:", sum(!is.finite(m$bm_democ_w15)), "\n")
m$bm_cinc_w15[!is.finite(m$bm_cinc_w15)] <- NA; m$bm_democ_w15[!is.finite(m$bm_democ_w15)] <- NA

pm_att <- function(dat, outcome, extra, nb, tag) {
  f <- paste0("~ I(lag(", outcome, ", 1:5)) + I(lag(ln_cinc, 1:5)) + I(lag(ln_pecpp, 1:5)) + I(lag(terr_dispute, 1:5))",
              if (length(extra)) paste0(" + ", paste(extra, collapse = " + ")) else "")
  set.seed(1)
  panel <- PanelData(panel.data = dat, unit.id = "country", time.id = "year", treatment = "tr", outcome = outcome)
  pm <- PanelMatch(panel.data = panel, lag = 5, refinement.method = "CBPS.weight", match.missing = TRUE, exact.match.variables = c("un_region"),
                   covs.formula = as.formula(f), qoi = "att", lead = 1:3, forbid.treatment.reversal = FALSE, use.diagonal.variance.matrix = TRUE)
  pe <- PanelEstimate(panel.data = panel, sets = pm, se.method = "bootstrap", number.iterations = nb, confidence.level = .90, pooled = TRUE)
  est <- unname(pe[["estimate"]]); se <- unname(pe[["standard.error"]]); q <- quantile(pe[["bootstrapped.estimates"]], c(.05, .95))
  cat(sprintf("---- %s | %s | covs %s | ATT(pooled 1:3) = %.5f  boot se %.5f  90%% CI [%.4f, %.4f]  p(normal) %.4g | treated sets %d\n",
              tag, outcome, f, est, se, q[1], q[2], 2 * pnorm(-abs(est / se)), length(pm[["att"]])))
  c(b = est, se = se, p = 2 * pnorm(-abs(est / se)))
}
ADDM <- list(m0 = character(0), mD = NULL, mDM = "hist_mono_w15", m26 = "hist_mono_w26", bC = "bm_cinc_w15", bD = "bm_democ_w15")
t0 <- Sys.time()
full <- list(mid = pm_att(m, "mid", character(0), NB, "FULL 1900-2014"), fatal = pm_att(m, "fatal", character(0), NB, "FULL 1900-2014"))
cat("time for full-sample pair:", format(Sys.time() - t0), "\n")
s <- m[m$year <= 2011, ]
res <- list()
for (k in names(ADDM)) {
  if (is.null(ADDM[[k]])) { res[[k]] <- NULL; next }
  res[[k]] <- list(mid = pm_att(s, "mid", ADDM[[k]], NB, k), fatal = pm_att(s, "fatal", ADDM[[k]], NB, k))
  saveRDS(list(full = full, res = res), file.path(W, paste0("M_P8_042_nb", NB, ".rds")))
}
cat("total time:", format(Sys.time() - t0), "\n")
# ---- port: result-CSV write step for this study, copied verbatim from work_p8/P8_042/write_P8_042.R (paths: saved model list read from OUT) ----
R <- readRDS(file.path(W, "M_P8_042_nb10000.rds")); stopifnot(NB == 10000L); stopifnot(all(c("m0", "mDM", "m26", "bC", "bD") %in% names(R$res)))
m <- as.data.frame(readRDS(file.path(WD, "dem-transitions-replication-data.rds"))); n_used <- sum(!duplicated(m[m$year <= 2011, c("country", "year")]))
meta <- list(pid = "P8_042", cite = "Boldt 2026", journal = "II", year = 2026, doi = "10.1080/03050629.2025.2609092", family = "DEMOCRACY")
for (o in c("mid", "fatal")) {
  term <- paste0("att_", o, "_pool3")
  M <- lapply(c(m0 = "m0", mD = "mD", mDM = "mDM", m26 = "m26", bC = "bC", bD = "bD"), function(k) {
    if (is.null(R$res[[k]])) return(NULL); x <- setNames(list(R$res[[k]][[o]]), term); attr(x, "n") <- n_used; x })
  f <- R$full[[o]]
  write_p8(meta, term, "abstract: 'Pooling across three- and five-year periods after transition, democratic transitions reduce conflict, including fatal conflicts'",
    paste0("Figure 6 (pooled-mids.R): PanelMatch ATT pooled over leads 1-3, outcome = ", o), "PanelMatch TSCS matching + DiD (lag 5, CBPS weights, exact UN region), bootstrap SE (10,000)",
    "PanelMatch 3.1.3 with the authors' call (10,000 bootstraps, seed 1), all steps; history/benchmarks added to covs.formula", NA, NA, NA, FALSE, FALSE, M,
    "Y (matching on 5 lags of the outcome)",
    sprintf(paste("MONADIC country-years (mD NA; mDM = +own MID-onset history w15; m26 = +own w26; bC = own CINC mean; bD = own polity>=6 mean). Published result is a",
      "figure only (no numeric table accessible; paywalled PDF; author log has no estimates) -> UNVERIFIED. Full-sample (1900-2014) rerun ATT = %.4f (boot se %.4f).",
      "Common sample = years <= 2011, %d country-years; p from normal approx. to bootstrap SE (authors use 90%% percentile CIs)."), f[["b"]], f[["se"]], n_used),
    verdict = "UNVERIFIED")
}
sink()
