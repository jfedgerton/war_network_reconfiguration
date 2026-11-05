# Phase 8 shared helpers (batch 1). Sequential, fail-loud. Seed 123.
ROOT <- Sys.getenv("SWEEP_DATA"); OUT <- Sys.getenv("SWEEP_OUT"); CODE <- Sys.getenv("SWEEP_CODE")
stopifnot(nzchar(ROOT), nzchar(OUT), nzchar(CODE))
suppressPackageStartupMessages({library(dplyr); library(fixest); library(sandwich)})
source(file.path(CODE, "hist/add_history.R"))
source(file.path(CODE, "work_p8/round2_rule.R"))   # rebuild_verdict(): round-2 rule with registry alpha_pub (2026-10-11)
set.seed(123)
HV <- c("hist_dyad_w15", "hist_mono_w15", "hist_dyad_w26", "hist_mono_w26", "bm_cinc_w15", "bm_democ_w15")
ADD <- list(m0 = character(0), mD = "hist_dyad_w15", mDM = c("hist_dyad_w15", "hist_mono_w15"),
            m26 = c("hist_dyad_w26", "hist_mono_w26"), bC = "bm_cinc_w15", bD = "bm_democ_w15")

# common sample: history keyed to the OUTCOME year (hy); 1822-2011; complete on model vars + HV
common_sample <- function(d, vars, c1, c2, hy, minfrac = 0.8) {
  n_cc <- sum(complete.cases(d[, vars]))
  d <- add_history(d, c1, c2, hy)
  d <- d[d[[hy]] >= 1822 & d[[hy]] <= 2011, ]
  d <- d[complete.cases(d[, c(vars, HV)]) & is.finite(d$bm_cinc_w15) & is.finite(d$bm_democ_w15), ]
  cat("complete-case N:", n_cc, " common N:", nrow(d), " dropped:", n_cc - nrow(d), "\n")
  stopifnot(nrow(d) > minfrac * n_cc, sd(d$hist_dyad_w15) > 0)
  d
}

# fitter(rhs, data) must return named list term -> c(b, se, p) with attr "n"
run_six <- function(rhs0, ds, fitter) {
  M <- lapply(ADD, function(a) fitter(c(rhs0, a), ds))
  ns <- sapply(M, function(m) attr(m, "n")); print(ns)
  stopifnot(length(unique(ns)) == 1)
  M
}

# fixest binary GLM (logit/probit), clustered or iid
fit_bin <- function(dv, rhs, data, link = "logit", vc = "iid", terms, fe = NULL) {
  f <- paste(dv, "~", paste(rhs, collapse = " + ")); if (!is.null(fe)) f <- paste(f, "|", fe)
  v <- if (grepl("^~", vc)) as.formula(vc) else vc
  m <- feglm(as.formula(f), data = data, family = binomial(link = link), vcov = v, notes = FALSE, lean = TRUE)
  stopifnot(m$convStatus)
  ct <- coeftable(m); stopifnot(all(terms %in% rownames(ct)))
  cat("\n----", f, "| N =", nobs(m), "| vcov:", vc, "| link:", link, "\n"); print(round(ct, 6))
  out <- lapply(setNames(terms, terms), function(t) c(b = ct[t, 1], se = ct[t, 2], p = ct[t, 4]))
  attr(out, "n") <- nobs(m); rm(m); invisible(gc()); out
}

pct <- function(b, b0) 100 * (b - b0) / abs(b0)
verdict <- function(b0, p0, pub_b, pub_sig, rebuilt) {
  if (is.na(pub_b)) return("NO_PUBLISHED")
  same_sig <- (p0 < .05) == pub_sig; same_sign <- sign(b0) == sign(pub_b)
  if (rebuilt) stop("rebuild verdicts go through rebuild_verdict() (round-2 rule, work_p8/round2_rule.R)")
  if (abs(b0 - pub_b) / abs(pub_b) <= .10 && same_sig && same_sign) "REPRODUCED" else if (same_sign && same_sig) "CLOSE" else "NOT_REPRODUCED"
}

# meta: list(pid, short_cite, journal, year, doi, family, headline_table_model, article_estimator, est_used_all_steps, author_has_history_control)
write_p8 <- function(meta, term, key_evidence, pub_b, pub_se, pub_n, pub_sig, rebuilt, M, notes, verdict_override = NULL, full = NULL) {
  g <- function(m, x) unname(M[[m]][[term]][[x]])
  b0 <- g("m0", "b"); p0 <- g("m0", "p")
  if (is.na(pub_sig) && !is.na(pub_b) && !is.na(pub_se)) pub_sig <- 2 * pnorm(-abs(pub_b / pub_se)) < .05
  vfun <- function(b, p) if (rebuilt) rebuild_verdict(meta$pid, term, b, p, pub_b, pub_sig, pub_se) else verdict(b, p, pub_b, pub_sig, rebuilt)
  v <- if (!is.null(verdict_override)) verdict_override else vfun(b0, p0)
  if (is.null(verdict_override) && !is.null(full) && !is.na(pub_b)) {   # reproduction judged on the authors' full sample
    vf <- vfun(unname(full[[term]]["b"]), unname(full[[term]]["p"]))
    notes <- paste0(notes, " | VERDICT BASIS: authors' full sample (b = ", signif(full[[term]]["b"], 4), ", p = ", signif(full[[term]]["p"], 3),
                    ", N = ", attr(full, "n"), "); verdict of b0 on the 1822-2011 common sample would be ", v, ".")
    v <- vf
  }
  row <- data.frame(pid = meta$pid, short_cite = meta$short_cite, journal = meta$journal, year = meta$year, doi = meta$doi,
    family = meta$family, key_term = term, key_evidence = key_evidence, headline_table_model = meta$headline_table_model,
    article_estimator = meta$article_estimator, est_used_all_steps = meta$est_used_all_steps,
    pub_b = pub_b, pub_se = pub_se, pub_n = pub_n, verdict = v,
    b0 = b0, se0 = g("m0", "se"), p0 = p0, n_used = attr(M$m0, "n"),
    pct_vs_pub = if (is.na(pub_b)) NA else 100 * (b0 - pub_b) / abs(pub_b),
    b_d = g("mD", "b"), pct_d = pct(g("mD", "b"), b0),
    b_dm = g("mDM", "b"), se_dm = g("mDM", "se"), p_dm = g("mDM", "p"), pct_dm = pct(g("mDM", "b"), b0),
    b_26 = g("m26", "b"), p_26 = g("m26", "p"), pct_26 = pct(g("m26", "b"), b0),
    pct_bC = pct(g("bC", "b"), b0), pct_bD = pct(g("bD", "b"), b0),
    author_has_history_control = meta$author_has_history_control, notes = notes, stringsAsFactors = FALSE)
  dir.create(file.path(OUT, "results_p8"), recursive = TRUE, showWarnings = FALSE)
  write.csv(row, file.path(OUT, "results_p8", paste0(meta$pid, "__", term, ".csv")), row.names = FALSE)
  print(t(row[, c("pid", "key_term", "verdict", "n_used", "b0", "se0", "p0", "b_dm", "p_dm", "b_26", "p_26", "pct_bC", "pct_bD")])); invisible(row)  # port: id only, no citation
}

# row for studies that could not be estimated (NO_DATA / DEFERRED)
write_p8_empty <- function(meta, term, key_evidence, pub_b, pub_se, pub_n, verdict, notes) {
  row <- data.frame(pid = meta$pid, short_cite = meta$short_cite, journal = meta$journal, year = meta$year, doi = meta$doi,
    family = meta$family, key_term = term, key_evidence = key_evidence, headline_table_model = meta$headline_table_model,
    article_estimator = meta$article_estimator, est_used_all_steps = meta$est_used_all_steps,
    pub_b = pub_b, pub_se = pub_se, pub_n = pub_n, verdict = verdict,
    b0 = NA, se0 = NA, p0 = NA, n_used = NA, pct_vs_pub = NA, b_d = NA, pct_d = NA, b_dm = NA, se_dm = NA, p_dm = NA,
    pct_dm = NA, b_26 = NA, p_26 = NA, pct_26 = NA, pct_bC = NA, pct_bD = NA,
    author_has_history_control = meta$author_has_history_control, notes = notes, stringsAsFactors = FALSE)
  write.csv(row, file.path(OUT, "results_p8", paste0(meta$pid, "__", term, ".csv")), row.names = FALSE)
  invisible(row)
}

# GEE AR(1) (Stata xtgee ... corr(ar 1) force vce(robust)) for binary logit via the generic solver
source(file.path(CODE, "work/gee_ar1_ea1.R"))
fam_logit <- list(name = "binomial logit", linkfun = qlogis, linkinv = plogis,
                  mu.eta = function(eta) { p <- plogis(eta); p * (1 - p) }, variance = function(mu) mu * (1 - mu),
                  valideta = function(eta) all(is.finite(eta)))
fit_gee_logit_ar1 <- function(dv, rhs, data, id, tt, terms) {
  X <- model.matrix(as.formula(paste("~", paste(rhs, collapse = " + "))), data)
  y <- data[[dv]]; stopifnot(nrow(X) == length(y))
  # column rescaling (no centering) for numerical conditioning; coefficients/SEs back-transformed exactly
  s <- apply(X, 2, sd); s[colnames(X) == "(Intercept)"] <- 1; stopifnot(all(s > 0))
  Xs <- sweep(X, 2, s, "/")
  b0 <- coef(glm.fit(Xs, y, family = binomial()))
  r <- gee_ar1(Xs, y, data[[id]], data[[tt]], fam_logit, beta0 = b0)
  r$coef <- r$coef / s; r$se <- r$se / s
  cat("\n---- GEE AR1 logit robust:", dv, "~", paste(rhs, collapse = " + "), "| N =", r$n, "| G =", r$G, "| rho =", round(r$rho, 4), "\n")
  print(round(cbind(b = r$coef, se = r$se, p = r$p), 6))
  out <- lapply(setNames(terms, terms), function(t) c(b = unname(r$coef[t]), se = unname(r$se[t]), p = unname(r$p[t])))
  attr(out, "n") <- r$n; out
}
