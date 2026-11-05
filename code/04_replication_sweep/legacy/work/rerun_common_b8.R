# Shared helpers for rerun batch cids 5 12 14 25 122 106 86.
# fit_tie(): fit one glm, return b/se/p/n of the tie variable (cluster-robust via sandwich, Stata-style G/(G-1)).
# run_six(): m0, mD, mDM, m26, bC, bD on one fixed sample.
# write_row(): results/<cid>.csv per RERUN_PROTOCOL step 5.
ROOT <- Sys.getenv("SWEEP_DATA"); OUT <- Sys.getenv("SWEEP_OUT"); CODE <- Sys.getenv("SWEEP_CODE")
stopifnot(nzchar(ROOT), nzchar(OUT), nzchar(CODE))
suppressPackageStartupMessages({library(sandwich); library(lmtest); library(dplyr)})
set.seed(123)

HIST_VARS <- c("hist_dyad_w15","hist_mono_w15","hist_dyad_w26","hist_mono_w26","bm_cinc_w15","bm_democ_w15")

fit_tie <- function(f, data, family, tie, cluster = NULL, quiet = FALSE) {
  m <- glm(f, data = data, family = family)
  stopifnot(m$converged)
  # sandwich default for glm: HC0 meat + G/(G-1) cluster adjustment == Stata ML cluster SEs
  V <- if (is.null(cluster)) vcov(m) else vcovCL(m, cluster = data[[cluster]])
  ct <- coeftest(m, vcov. = V)
  stopifnot(tie %in% rownames(ct))
  if (!quiet) { cat("\n---- ", paste(deparse(f, width.cutoff = 500), collapse = ""), " | N =", nobs(m), "\n"); print(signif(ct, 5)) }
  list(b = ct[tie, 1], se = ct[tie, 2], p = ct[tie, 4], n = nobs(m))
}

# Faster equivalent for large data: fixest::feglm, Stata-style ML cluster SE (G/(G-1) only, no (N-1)/(N-K)); z-based p.
fit_tie_fe <- function(f, data, family, tie, cluster, weights = NULL) {
  m <- fixest::feglm(f, data = data, family = family, weights = weights, cluster = cluster,
                     ssc = fixest::ssc(adj = FALSE, cluster.adj = TRUE), notes = FALSE, mem.clean = TRUE)
  stopifnot(m$convStatus)
  ct <- fixest::coeftable(m)
  stopifnot(tie %in% rownames(ct))
  cat("\n---- ", paste(deparse(f, width.cutoff = 500), collapse = ""), " | N =", nobs(m), "\n"); print(signif(ct, 5))
  out <- list(b = ct[tie, 1], se = ct[tie, 2], p = ct[tie, 4], n = nobs(m)); rm(m); gc(verbose = FALSE)
  out
}

# fitter: function(formula_string_addition, data) -> list(b,se,p,n)
run_six <- function(fitter, d) {
  cat("\n==== m0 (baseline on restricted sample) ====\n"); m0  <- fitter("", d)
  cat("\n==== mD  ====\n"); mD  <- fitter(" + hist_dyad_w15", d)
  cat("\n==== mDM ====\n"); mDM <- fitter(" + hist_dyad_w15 + hist_mono_w15", d)
  cat("\n==== m26 ====\n"); m26 <- fitter(" + hist_dyad_w26 + hist_mono_w26", d)
  cat("\n==== bC  ====\n"); bC  <- fitter(" + bm_cinc_w15", d)
  cat("\n==== bD  ====\n"); bD  <- fitter(" + bm_democ_w15", d)
  list(m0 = m0, mD = mD, mDM = mDM, m26 = m26, bC = bC, bD = bD)
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
