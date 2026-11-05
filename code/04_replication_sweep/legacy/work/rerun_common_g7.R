# Shared helpers for rerun batch cids 43 62 15 104 32 100 91.
# fit_glm_tie(): one glm, cluster/robust SEs Stata-style (HC0 meat, G/(G-1) adj), return tie b/se/p/n.
# run_six():  m0, mD, mDM, m26, bC, bD on one fixed sample (fitter(add_string, data)).
# restrict_sample(): years 1822-2011 + complete cases on model vars and all history/benchmark vars.
# write_row(): results/<cid>.csv per RERUN_PROTOCOL step 5.
ROOT <- Sys.getenv("SWEEP_DATA"); OUT <- Sys.getenv("SWEEP_OUT"); CODE <- Sys.getenv("SWEEP_CODE")
stopifnot(nzchar(ROOT), nzchar(OUT), nzchar(CODE))
suppressPackageStartupMessages({library(sandwich); library(lmtest); library(dplyr)})
source(file.path(CODE, "hist/add_history.R"))
set.seed(123)

HIST_VARS <- c("hist_dyad_w15","hist_mono_w15","hist_dyad_w26","hist_mono_w26","bm_cinc_w15","bm_democ_w15")

# cluster = column name; "ROBUST" = observation-level robust (Stata ", robust": HC0 * n/(n-1)); NULL = model SEs
fit_glm_tie <- function(f, data, family, tie, cluster = NULL, quiet = FALSE) {
  vv <- c(all.vars(f), if (!is.null(cluster) && !identical(cluster, "ROBUST")) cluster)
  data <- data[complete.cases(as.data.frame(data)[, vv]), ]   # listwise deletion up front so cluster vector aligns
  m <- glm(f, data = data, family = family)
  stopifnot(nobs(m) == nrow(data))
  stopifnot(m$converged)
  V <- if (is.null(cluster)) vcov(m) else if (identical(cluster, "ROBUST")) vcovCL(m, cluster = seq_len(nobs(m))) else vcovCL(m, cluster = data[[cluster]])
  ct <- coeftest(m, vcov. = V)
  stopifnot(tie %in% rownames(ct))
  if (!quiet) { cat("\n---- ", paste(deparse(f, width.cutoff = 500), collapse = ""), " | N =", nobs(m), "\n"); print(signif(ct, 5)) }
  list(b = unname(ct[tie, 1]), se = unname(ct[tie, 2]), p = unname(ct[tie, 4]), n = nobs(m))
}

restrict_sample <- function(d, model_vars, yearvar = "year") {
  n_in <- nrow(d)
  d <- d[d[[yearvar]] >= 1822 & d[[yearvar]] <= 2011, ]
  n_yr <- nrow(d)
  n_mv <- sum(complete.cases(as.data.frame(d)[, model_vars]))
  keep <- complete.cases(as.data.frame(d)[, c(model_vars, HIST_VARS)]) &
          is.finite(d$bm_cinc_w15) & is.finite(d$bm_democ_w15)
  d <- d[keep, ]
  cat(sprintf("\nrestrict_sample: in=%d, after 1822-2011=%d, complete on model vars=%d, after complete cases (model+hist+bm)=%d\n", n_in, n_yr, n_mv, nrow(d)))
  RS_INFO <<- sprintf("1822-2011 & model-var complete N=%d; bm/hist missingness drops %d more", n_mv, n_mv - nrow(d))
  stopifnot(nrow(d) > 0)
  d
}

run_six <- function(fitter, d) {
  cat("\n==== m0 (baseline on restricted sample) ====\n"); m0  <- fitter("", d)
  cat("\n==== mD  ====\n"); mD  <- fitter(" + hist_dyad_w15", d)
  cat("\n==== mDM ====\n"); mDM <- fitter(" + hist_dyad_w15 + hist_mono_w15", d)
  cat("\n==== m26 ====\n"); m26 <- fitter(" + hist_dyad_w26 + hist_mono_w26", d)
  cat("\n==== bC  ====\n"); bC  <- fitter(" + bm_cinc_w15", d)
  cat("\n==== bD  ====\n"); bD  <- fitter(" + bm_democ_w15", d)
  R <- list(m0 = m0, mD = mD, mDM = mDM, m26 = m26, bC = bC, bD = bD)
  ns <- sapply(R, `[[`, "n"); cat("\nN by model:", ns, "\n")
  stopifnot(length(unique(ns)) == 1)   # same estimation sample in all six
  R
}

pct <- function(b, b0) 100 * (b - b0) / abs(b0)

write_row <- function(cid, short_cite, reproduced, b_pub, b_repro, estimator, tie_var, n_orig, n_used,
                      R = NULL, ahc, notes) {
  g <- function(m, k) if (is.null(m)) NA else m[[k]]
  m0 <- R$m0; mD <- R$mD; mDM <- R$mDM; m26 <- R$m26; bC <- R$bC; bD <- R$bD
  b0 <- g(m0, "b")
  row <- data.frame(cid = cid, short_cite = short_cite, reproduced = reproduced,
    b_pub = b_pub, b_repro = b_repro, estimator = estimator, tie_var = tie_var,
    n_orig = n_orig, n_used = n_used,
    b0 = b0, se0 = g(m0, "se"), p0 = g(m0, "p"),
    b_dm = g(mDM, "b"), se_dm = g(mDM, "se"), p_dm = g(mDM, "p"), pct_dm = pct(g(mDM, "b"), b0),
    b_d = g(mD, "b"), pct_d = pct(g(mD, "b"), b0),
    b_26 = g(m26, "b"), se_26 = g(m26, "se"), p_26 = g(m26, "p"), pct_26 = pct(g(m26, "b"), b0),
    pct_bC = pct(g(bC, "b"), b0), pct_bD = pct(g(bD, "b"), b0),
    author_has_history_control = ahc, notes = notes, stringsAsFactors = FALSE)
  dir.create(file.path(OUT, "results"), recursive = TRUE, showWarnings = FALSE)
  write.csv(row, file.path(OUT, "results", paste0(cid, ".csv")), row.names = FALSE)
  print(t(row[, names(row) != "short_cite", drop = FALSE]))
  invisible(row)
}

# Sensitivity: benchmark vars (bm_*) are NaN for states with no NMC/Polity data in t-1..t-5 (new states), so the
# protocol's common sample drops them from ALL models. Re-fit m0 and mDM on the sample requiring only model + hist vars
# (hist vars are never missing) and return a note string. Only run when that sample is larger.
sens_hist_only <- function(fitter, d_all, model_vars, yearvar = "year") {
  d2 <- d_all[d_all[[yearvar]] >= 1822 & d_all[[yearvar]] <= 2011, ]
  d2 <- d2[complete.cases(as.data.frame(d2)[, c(model_vars, HIST_VARS[1:4])]), ]
  cat("\n==== Sensitivity: sample requiring model + hist vars only (N =", nrow(d2), ") ====\n")
  s0 <- fitter("", d2); sdm <- fitter(" + hist_dyad_w15 + hist_mono_w15", d2)
  sprintf("hist-only sample N=%d: b0=%.3f (p=%.3g), b_dm=%.3f (p=%.3g), pct_dm=%.1f", s0$n, s0$b, s0$p, sdm$b, sdm$p, pct(sdm$b, s0$b))
}

# Memory-lean alternative for very large data: fixest::feglm logit, cluster SE with Stata ML adjustment (G/(G-1) only).
fit_feglm_tie <- function(f, data, tie, cluster, quiet = FALSE) {
  suppressPackageStartupMessages(library(fixest))
  vc <- if (identical(cluster, "ROBUST")) "hetero" else as.formula(paste0("~", cluster))
  m <- feglm(f, data = data, family = binomial(), vcov = vc,
             ssc = ssc(adj = FALSE, cluster.adj = TRUE), notes = FALSE, lean = TRUE, mem.clean = TRUE)
  stopifnot(m$convStatus)
  ct <- coeftable(m)
  stopifnot(tie %in% rownames(ct))
  if (!quiet) { cat("\n---- ", paste(deparse(f, width.cutoff = 500), collapse = ""), " | N =", nobs(m), "\n"); print(signif(ct, 5)) }
  out <- list(b = unname(ct[tie, 1]), se = unname(ct[tie, 2]), p = unname(ct[tie, 4]), n = nobs(m), all = ct)
  rm(m); invisible(gc())
  out
}
