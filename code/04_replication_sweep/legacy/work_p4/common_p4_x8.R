# Phase-4 helpers for batch x8 (X008 X016 ... X088). Sequential, fail-loud.
ROOT <- Sys.getenv("SWEEP_DATA"); OUT <- Sys.getenv("SWEEP_OUT"); CODE <- Sys.getenv("SWEEP_CODE")
stopifnot(nzchar(ROOT), nzchar(OUT), nzchar(CODE))
suppressPackageStartupMessages({library(dplyr)})
source(file.path(CODE, "hist/add_history.R"))
set.seed(123)

pct <- function(b, b0) 100 * (b - b0) / abs(b0)

# ct: coefficient matrix (rows = terms; cols Estimate, SE, ..., p in col 4). Returns named list per term.
grab <- function(ct, terms, n, label) {
  cat("\n----", label, "| N =", n, "\n"); print(round(ct, 6))
  stopifnot(all(terms %in% rownames(ct)))
  setNames(lapply(terms, function(t) list(b = ct[t, 1], se = ct[t, 2], p = ct[t, 4], n = n)), terms)
}

# fixest logit/probit/OLS with clustered or hetero-robust SEs
fit_fx <- function(f, data, terms, family = "logit", vcov = "hetero", label = "") {
  if (family == "ols") m <- fixest::feols(f, data = data, vcov = vcov, notes = FALSE, lean = TRUE)
  else m <- fixest::feglm(f, data = data, family = family, vcov = vcov, notes = FALSE, lean = TRUE,
                          ssc = fixest::ssc(adj = FALSE, cluster.adj = TRUE))
  if (family != "ols") stopifnot(m$convStatus)
  ct <- fixest::coeftable(m); n <- nobs(m); rm(m); gc()
  grab(ct, terms, n, paste(label, deparse(f, width.cutoff = 500)))
}

# one CSV per term: results_p4/<xid>__<term>.csv
write_p4 <- function(meta, term, family, role, evidence, reproduced, b_pub, b_repro, estimator, n_orig,
                     M = NULL, ahc, notes) {
  g <- function(k, f) if (is.null(M) || is.null(M[[k]])) NA else M[[k]][[term]][[f]]
  b0 <- g("m0", "b")
  row <- data.frame(xid = meta$xid, short_cite = meta$short_cite, journal = meta$journal, year = meta$year,
    doi = meta$doi, key_family = family, role = role, key_term = term, key_evidence = evidence,
    reproduced = reproduced, b_pub = b_pub, b_repro = b_repro, estimator = estimator,
    n_orig = n_orig, n_used = g("m0", "n"),
    b0 = b0, se0 = g("m0", "se"), p0 = g("m0", "p"),
    b_d = g("mD", "b"), pct_d = pct(g("mD", "b"), b0),
    b_dm = g("mDM", "b"), se_dm = g("mDM", "se"), p_dm = g("mDM", "p"), pct_dm = pct(g("mDM", "b"), b0),
    b_26 = g("m26", "b"), p_26 = g("m26", "p"), pct_26 = pct(g("m26", "b"), b0),
    pct_bC = pct(g("bC", "b"), b0), pct_bD = pct(g("bD", "b"), b0),
    author_has_history_control = ahc, notes = notes, stringsAsFactors = FALSE)
  dir.create(file.path(OUT, "results_p4"), recursive = TRUE, showWarnings = FALSE)
  fn <- file.path(OUT, "results_p4", paste0(meta$xid, "__", gsub("[^A-Za-z0-9_]", "_", term), ".csv"))
  write.csv(row, fn, row.names = FALSE); cat("wrote", fn, "\n"); print(t(row[, intersect(c("xid", "id", "key_term", "n_used", "b0", "se0", "p0", "b_dm", "se_dm", "p_dm", "b_26", "p_26", "pct_bC", "pct_bD"), names(row)), drop = FALSE]))  # port: id only, no citation printed
  invisible(row)
}

# Run the six models on one common sample. fitfun(formula_rhs_extra, data) -> grab() list.
run_six <- function(fitfun, s) {
  list(m0 = fitfun("", s), mD = fitfun("+ hist_dyad_w15", s), mDM = fitfun("+ hist_dyad_w15 + hist_mono_w15", s),
       m26 = fitfun("+ hist_dyad_w26 + hist_mono_w26", s), bC = fitfun("+ bm_cinc_w15", s), bD = fitfun("+ bm_democ_w15", s))
}
HISTV <- c("hist_dyad_w15", "hist_mono_w15", "hist_dyad_w26", "hist_mono_w26", "bm_cinc_w15", "bm_democ_w15")

# Monadic analogue: own-state MID-onset history (w15 = t-1..t-5, w26 = t-2..t-6) + own CINC / own democracy benchmarks.
add_history_monadic <- function(df, cc, year) {
  H <- .H
  k <- data.frame(ccode = as.numeric(df[[cc]]), year = as.numeric(df[[year]])) %>%
    left_join(H$sth, by = c("ccode", "year")) %>% mutate(across(c(s15, s26), ~coalesce(.x, 0)))
  stopifnot(nrow(k) == nrow(df))
  df$hist_mono_w15 <- log1p(k$s15); df$hist_mono_w26 <- log1p(k$s26)
  bm <- lapply(1:5, function(L) {
    x <- data.frame(ccode = k$ccode, yl = k$year - L) %>%
      left_join(H$nmc, by = c("ccode", "yl" = "year")) %>% left_join(H$pol, by = c("ccode", "yl" = "year"))
    data.frame(c = x$cinc, d = as.numeric(x$polity2 >= 6))
  })
  df$bm_cinc_w15  <- rowMeans(sapply(bm, `[[`, "c"), na.rm = TRUE)
  df$bm_democ_w15 <- rowMeans(sapply(bm, `[[`, "d"), na.rm = TRUE)
  df
}
