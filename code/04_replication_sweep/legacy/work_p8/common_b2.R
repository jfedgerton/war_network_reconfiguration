# Phase 8, batch 2 shared helpers. Sequential, fail-loud. Seed 123 unless authors set one.
ROOT <- Sys.getenv("SWEEP_DATA"); OUT <- Sys.getenv("SWEEP_OUT"); CODE <- Sys.getenv("SWEEP_CODE")
stopifnot(nzchar(ROOT), nzchar(OUT), nzchar(CODE))
suppressPackageStartupMessages({library(dplyr); library(fixest)})
source(file.path(CODE, "hist/add_history.R"))   # add_history(), .H
source(file.path(CODE, "work_p8/round2_rule.R"))   # rebuild_verdict(): round-2 rule with registry alpha_pub (2026-10-11)
set.seed(123)
HV <- c("hist_dyad_w15", "hist_mono_w15", "hist_dyad_w26", "hist_mono_w26", "bm_cinc_w15", "bm_democ_w15")
ADD  <- list(m0 = character(0), mD = "hist_dyad_w15", mDM = c("hist_dyad_w15", "hist_mono_w15"),
             m26 = c("hist_dyad_w26", "hist_mono_w26"), bC = "bm_cinc_w15", bD = "bm_democ_w15")
# monadic convention (as work_p4 X052): no dyadic term -> mD NA; mDM = + own history w15; m26 = + own history w26
ADD_MON <- list(m0 = character(0), mD = NULL, mDM = "hist_mono_w15", m26 = "hist_mono_w26", bC = "bm_cinc_w15", bD = "bm_democ_w15")

add_history_monadic <- function(df, cc, year) {
  k <- data.frame(ccode = as.numeric(df[[cc]]), year = as.numeric(df[[year]])) %>%
    left_join(.H$sth, by = c("ccode", "year")) %>% mutate(across(c(s15, s26), ~coalesce(.x, 0)))
  stopifnot(nrow(k) == nrow(df))
  df$hist_mono_w15 <- log1p(k$s15); df$hist_mono_w26 <- log1p(k$s26)
  bm <- lapply(1:5, function(L) {
    x <- data.frame(ccode = k$ccode, yl = k$year - L) %>%
      left_join(.H$nmc, by = c("ccode", "yl" = "year")) %>% left_join(.H$pol, by = c("ccode", "yl" = "year"))
    data.frame(c = x$cinc, d = as.numeric(x$polity2 >= 6))
  })
  df$bm_cinc_w15  <- rowMeans(sapply(bm, `[[`, "c"), na.rm = TRUE)
  df$bm_democ_w15 <- rowMeans(sapply(bm, `[[`, "d"), na.rm = TRUE)
  df
}

# common sample 1822-2011, complete on model vars + history vars
common_sample <- function(d, vars, c1, c2, yr, monadic = FALSE) {
  n_cc <- sum(complete.cases(d[, vars]))
  d <- if (monadic) add_history_monadic(d, c1, yr) else add_history(d, c1, c2, yr)
  hv <- if (monadic) c("hist_mono_w15", "hist_mono_w26", "bm_cinc_w15", "bm_democ_w15") else HV
  d <- d[d[[yr]] >= 1822 & d[[yr]] <= 2011, ]
  d <- d[complete.cases(d[, c(vars, hv)]) & is.finite(d$bm_cinc_w15) & is.finite(d$bm_democ_w15), ]
  cat("complete-case N:", n_cc, " common N:", nrow(d), " dropped:", n_cc - nrow(d), "\n")
  stopifnot(nrow(d) > 0.8 * n_cc, sd(d$hist_mono_w15) > 0)
  if (!monadic) stopifnot(sd(d$hist_dyad_w15) > 0)
  attr(d, "n_cc") <- n_cc
  d
}

# run six models; fitter(rhs_vector, data) -> named list term -> c(b,se,p), attr n
run_six <- function(rhs0, ds, fitter, monadic = FALSE) {
  A <- if (monadic) ADD_MON else ADD
  M <- lapply(names(ADD), function(k) if (is.null(A[[k]]) && k == "mD" && monadic) NULL else fitter(c(rhs0, A[[k]]), ds))
  names(M) <- names(ADD)
  ns <- sapply(M[!sapply(M, is.null)], function(m) attr(m, "n")); print(ns)
  stopifnot(length(unique(ns)) == 1)
  M
}

pct <- function(b, b0) 100 * (b - b0) / abs(b0)
verdict_fn <- function(b0, p0, pub_b, pub_sig, rebuilt) {
  if (is.na(pub_b)) return("NO_PUBLISHED")
  same_sig <- (p0 < .05) == pub_sig; same_sign <- sign(b0) == sign(pub_b)
  if (rebuilt) stop("rebuild verdicts go through rebuild_verdict() (round-2 rule, work_p8/round2_rule.R)")
  if (abs(b0 - pub_b) / abs(pub_b) <= .10 && same_sig && same_sign) "REPRODUCED" else if (same_sign && same_sig) "CLOSE" else "NOT_REPRODUCED"
}

write_p8 <- function(meta, term, evidence, hm, art_est, est_used, pub_b, pub_se, pub_n, pub_sig, rebuilt, M, ahc, notes,
                     verdict = NULL) {
  g <- function(m, x) if (is.null(M[[m]])) NA_real_ else unname(M[[m]][[term]][[x]])
  b0 <- g("m0", "b"); p0 <- g("m0", "p")
  v <- if (!is.null(verdict)) verdict else if (rebuilt) rebuild_verdict(meta$pid, term, b0, p0, pub_b, pub_sig, pub_se) else verdict_fn(b0, p0, pub_b, pub_sig, rebuilt)
  row <- data.frame(pid = meta$pid, short_cite = meta$cite, journal = meta$journal, year = meta$year, doi = meta$doi,
    family = meta$family, key_term = term, key_evidence = evidence, headline_table_model = hm, article_estimator = art_est,
    est_used_all_steps = est_used, pub_b = pub_b, pub_se = pub_se, pub_n = pub_n, verdict = v,
    b0 = b0, se0 = g("m0", "se"), p0 = p0, n_used = attr(M$m0, "n"),
    pct_vs_pub = if (is.na(pub_b)) NA else 100 * (b0 - pub_b) / abs(pub_b),
    b_d = g("mD", "b"), pct_d = pct(g("mD", "b"), b0),
    b_dm = g("mDM", "b"), se_dm = g("mDM", "se"), p_dm = g("mDM", "p"), pct_dm = pct(g("mDM", "b"), b0),
    b_26 = g("m26", "b"), p_26 = g("m26", "p"), pct_26 = pct(g("m26", "b"), b0),
    pct_bC = pct(g("bC", "b"), b0), pct_bD = pct(g("bD", "b"), b0),
    author_has_history_control = ahc, notes = notes, stringsAsFactors = FALSE)
  dir.create(file.path(OUT, "results_p8"), recursive = TRUE, showWarnings = FALSE)
  write.csv(row, file.path(OUT, "results_p8", paste0(meta$pid, "__", term, ".csv")), row.names = FALSE)
  print(t(row[, c("pid", "key_term", "verdict", "n_used", "b0", "se0", "p0", "b_dm", "p_dm", "b_26", "p_26", "pct_bC", "pct_bD")])); invisible(row)  # port: id only, no citation
}

# row for studies without estimates (NO_DATA / DEFERRED)
write_p8_empty <- function(meta, term, evidence, hm, art_est, pub_b, pub_se, pub_n, verdict, ahc, notes) {
  row <- data.frame(pid = meta$pid, short_cite = meta$cite, journal = meta$journal, year = meta$year, doi = meta$doi,
    family = meta$family, key_term = term, key_evidence = evidence, headline_table_model = hm, article_estimator = art_est,
    est_used_all_steps = NA, pub_b = pub_b, pub_se = pub_se, pub_n = pub_n, verdict = verdict,
    b0 = NA, se0 = NA, p0 = NA, n_used = NA, pct_vs_pub = NA, b_d = NA, pct_d = NA, b_dm = NA, se_dm = NA, p_dm = NA, pct_dm = NA,
    b_26 = NA, p_26 = NA, pct_26 = NA, pct_bC = NA, pct_bD = NA, author_has_history_control = ahc, notes = notes, stringsAsFactors = FALSE)
  dir.create(file.path(OUT, "results_p8"), recursive = TRUE, showWarnings = FALSE)
  write.csv(row, file.path(OUT, "results_p8", paste0(meta$pid, "__", term, ".csv")), row.names = FALSE)
  print(t(row[, c("pid", "key_term", "verdict", "n_used", "b0", "se0", "p0", "b_dm", "p_dm", "b_26", "p_26", "pct_bC", "pct_bD")])); invisible(row)  # port: id only, no citation
}

# fixest binary GLM returning term list
fit_bin <- function(dv, rhs, data, link = "logit", vc = "iid", fe = NULL, terms) {
  f <- paste(dv, "~", paste(rhs, collapse = " + "))
  if (!is.null(fe)) f <- paste(f, "|", fe)
  v <- if (grepl("^~", vc)) as.formula(vc) else vc
  m <- feglm(as.formula(f), data = data, family = binomial(link = link), vcov = v, ssc = ssc(adj = FALSE, cluster.adj = TRUE),
             notes = FALSE, lean = TRUE, nthreads = 1)   # Stata ML cluster SE: G/(G-1) only
  stopifnot(m$convStatus)
  ct <- coeftable(m); stopifnot(all(terms %in% rownames(ct)))
  cat("\n----", f, "| N =", nobs(m), "| vcov:", vc, "| link:", link, "\n"); print(round(ct, 6))
  out <- lapply(setNames(terms, terms), function(t) c(b = ct[t, 1], se = ct[t, 2], p = ct[t, 4]))
  attr(out, "n") <- nobs(m); rm(m); invisible(gc()); out
}
