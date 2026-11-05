# Shared helpers for rerun batch cids 102 64 89 54 53 13.
# fit_tie():        one glm (logit/probit, optional weights), Stata-style cluster SEs (HC0 meat * G/(G-1)); returns tie b/se/p/n.
# restrict_sample(): years 1822-2011 + complete cases on model vars and all history/benchmark vars.
# run_six():        m0, mD, mDM, m26, bC, bD on one fixed sample via fitter(add_string, data).
# write_row():      results/<cid>.csv per RERUN_PROTOCOL step 5.
ROOT <- Sys.getenv("SWEEP_DATA"); OUT <- Sys.getenv("SWEEP_OUT"); CODE <- Sys.getenv("SWEEP_CODE")
stopifnot(nzchar(ROOT), nzchar(OUT), nzchar(CODE))
suppressPackageStartupMessages({library(sandwich); library(lmtest); library(dplyr); library(haven)})
source(file.path(CODE, "hist/add_history.R"))
set.seed(123)

HIST_VARS <- c("hist_dyad_w15","hist_mono_w15","hist_dyad_w26","hist_mono_w26","bm_cinc_w15","bm_democ_w15")

# cluster: NULL = model-based SEs; column name = Stata cluster(); weights: column name or NULL
fit_tie <- function(f, data, link = "logit", tie, cluster = NULL, weights = NULL, quiet = FALSE) {
  fam <- binomial(link = link)
  data <- as.data.frame(data); rownames(data) <- NULL
  if (is.null(weights)) {
    m <- glm(f, data = data, family = fam, control = glm.control(maxit = 100))
  } else {
    data$.w <- data[[weights]]
    m <- suppressWarnings(glm(f, data = data, family = fam, weights = .w, control = glm.control(maxit = 100)))
  }
  stopifnot(m$converged)
  if (is.null(cluster)) {
    V <- vcov(m)
  } else {
    cl <- data[[cluster]][as.integer(rownames(model.frame(m)))]
    stopifnot(length(cl) == nrow(model.frame(m)), !anyNA(cl))
    V <- vcovCL(m, cluster = cl, type = "HC0", cadjust = TRUE)
  }
  ct <- coeftest(m, vcov. = V)
  stopifnot(tie %in% rownames(ct))
  n <- if (is.null(weights)) nobs(m) else sum(m$prior.weights > 0)
  if (!quiet) { cat("\n---- ", paste(deparse(f, width.cutoff = 500), collapse = ""), " | link =", link, "| N =", n, "\n"); print(signif(ct, 5)) }
  list(b = unname(ct[tie, 1]), se = unname(ct[tie, 2]), p = unname(ct[tie, 4]), n = n)
}

restrict_sample <- function(d, model_vars, yearvar = "year") {
  n_in <- nrow(d)
  d <- d[d[[yearvar]] >= 1822 & d[[yearvar]] <= 2011, ]
  n_yr <- nrow(d)
  keep <- complete.cases(as.data.frame(d)[, c(model_vars, HIST_VARS)]) &
          is.finite(d$bm_cinc_w15) & is.finite(d$bm_democ_w15)
  d <- d[keep, ]
  cat(sprintf("\nrestrict_sample: in=%d, after 1822-2011=%d, after complete cases (model+hist+bm)=%d\n", n_in, n_yr, nrow(d)))
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

repro_status <- function(b_pub, p_pub, b_rep, p_rep) {
  if (is.na(b_pub)) return("UNVERIFIED")
  ok <- abs(b_rep - b_pub) / abs(b_pub) <= 0.10 && ((p_pub < .05) == (p_rep < .05))
  if (ok) "REPRODUCED" else "NOT_REPRODUCED"
}

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
