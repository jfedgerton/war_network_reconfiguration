# Verify-phase helpers (batch 4). Sequential, fail-loud. Seed 123.
ROOT <- Sys.getenv("SWEEP_DATA"); OUT <- Sys.getenv("SWEEP_OUT"); CODE <- Sys.getenv("SWEEP_CODE")
stopifnot(nzchar(ROOT), nzchar(OUT), nzchar(CODE))
suppressPackageStartupMessages({library(dplyr); library(fixest)})
source(file.path(CODE, "hist/add_history.R"))
source(file.path(CODE, "work_p8/round2_rule.R"))   # rebuild_verdict(): round-2 rule with registry alpha_pub (2026-10-11)
set.seed(123)
HV <- c("hist_dyad_w15", "hist_mono_w15", "hist_dyad_w26", "hist_mono_w26", "bm_cinc_w15", "bm_democ_w15")
ADD <- list(m0 = character(0), mD = "hist_dyad_w15", mDM = c("hist_dyad_w15", "hist_mono_w15"),
            m26 = c("hist_dyad_w26", "hist_mono_w26"), bC = "bm_cinc_w15", bD = "bm_democ_w15")

# fixest binary GLM. vc = "iid" | "hetero" | "~dyad"-style string; fe = optional FE var name
fit_bin <- function(dv, rhs, data, link = "logit", vc = "iid", fe = NULL, terms) {
  f <- paste(dv, "~", paste(rhs, collapse = " + "))
  if (!is.null(fe)) f <- paste(f, "|", fe)
  v <- if (grepl("^~", vc)) as.formula(vc) else vc
  m <- feglm(as.formula(f), data = data, family = binomial(link = link), vcov = v, notes = FALSE, lean = TRUE)
  stopifnot(m$convStatus)
  ct <- coeftable(m); stopifnot(all(terms %in% rownames(ct)))
  cat("\n----", f, "| N =", nobs(m), "| vcov:", vc, "| link:", link, "\n"); print(round(ct, 6))
  out <- lapply(setNames(terms, terms), function(t) c(b = ct[t, 1], se = ct[t, 2], p = ct[t, 4]))
  attr(out, "n") <- nobs(m); rm(m); invisible(gc()); out
}

# six models on one common sample; fitter(rhs, data) returns fit_bin-style list
run_six <- function(rhs0, ds, fitter) {
  M <- lapply(ADD, function(a) fitter(c(rhs0, a), ds))
  ns <- sapply(M, function(m) attr(m, "n")); print(ns)
  stopifnot(length(unique(ns)) == 1)
  M
}

# common sample: 1822-2011 on hist year, complete on model vars + HV
common_sample <- function(d, vars, c1, c2, yr) {
  n_cc <- sum(complete.cases(d[, vars]))
  d <- add_history(d, c1, c2, yr)
  d <- d[d[[yr]] >= 1822 & d[[yr]] <= 2011, ]
  d <- d[complete.cases(d[, c(vars, HV)]) & is.finite(d$bm_cinc_w15) & is.finite(d$bm_democ_w15), ]
  cat("complete-case N:", n_cc, " common N:", nrow(d), " dropped:", n_cc - nrow(d), "\n")
  stopifnot(nrow(d) > 0.8 * n_cc, sd(d$hist_dyad_w15) > 0)
  d
}

pct <- function(b, b0) 100 * (b - b0) / abs(b0)
verdict <- function(b0, p0, pub_b, pub_sig, rebuilt) {
  if (is.na(pub_b)) return("NO_PUBLISHED")
  same_sig <- (p0 < .05) == pub_sig; same_sign <- sign(b0) == sign(pub_b)
  if (rebuilt) stop("rebuild verdicts go through rebuild_verdict() (round-2 rule, work_p8/round2_rule.R)")
  if (abs(b0 - pub_b) / abs(pub_b) <= .10 && same_sig) "REPRODUCED" else if (same_sign && same_sig) "CLOSE" else "NOT_REPRODUCED"
}

write_verify <- function(id, cite, fam, term, hm, pub_b, pub_se, pub_n, pub_est, pub_sig, mism, fixed, unfix, rebuilt, M, notes) {
  g <- function(m, x) unname(M[[m]][[term]][[x]])
  b0 <- g("m0", "b"); p0 <- g("m0", "p")
  v <- if (rebuilt) rebuild_verdict(id, term, b0, p0, pub_b, pub_sig, pub_se) else verdict(b0, p0, pub_b, pub_sig, rebuilt)
  row <- data.frame(id = id, short_cite = cite, key_family = fam, key_term = term, headline_table_model = hm,
    pub_b = pub_b, pub_se = pub_se, pub_n = pub_n, pub_estimator = pub_est, spec_mismatches_before = mism,
    fixed = fixed, unfixable = unfix, verdict = v, b0 = b0, se0 = g("m0", "se"), p0 = p0, n_used = attr(M$m0, "n"),
    pct_vs_pub = if (is.na(pub_b)) NA else 100 * (b0 - pub_b) / abs(pub_b),
    b_dm = g("mDM", "b"), se_dm = g("mDM", "se"), p_dm = g("mDM", "p"), pct_dm = pct(g("mDM", "b"), b0),
    b_26 = g("m26", "b"), p_26 = g("m26", "p"), pct_26 = pct(g("m26", "b"), b0),
    pct_bC = pct(g("bC", "b"), b0), pct_bD = pct(g("bD", "b"), b0), notes = notes, stringsAsFactors = FALSE)
  dir.create(file.path(OUT, "results_verify"), recursive = TRUE, showWarnings = FALSE)
  write.csv(row, file.path(OUT, "results_verify", paste0(id, "__", term, ".csv")), row.names = FALSE)
  print(t(row[, intersect(c("xid", "id", "key_term", "n_used", "b0", "se0", "p0", "b_dm", "se_dm", "p_dm", "b_26", "p_26", "pct_bC", "pct_bD"), names(row)), drop = FALSE])); invisible(row)  # port: id only, no citation printed
}
