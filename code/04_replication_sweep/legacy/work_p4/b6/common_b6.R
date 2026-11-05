# Phase 4 helpers (batch b6: X006 X014 X022 X030 X038 X046 X054 X062 X070 X078 X086).
# Adapted from work_demcap/common_c3467.R: adds probit/iid options and the phase-4 output schema.
ROOT <- Sys.getenv("SWEEP_DATA"); OUT <- Sys.getenv("SWEEP_OUT"); CODE <- Sys.getenv("SWEEP_CODE")
stopifnot(nzchar(ROOT), nzchar(OUT), nzchar(CODE))
suppressPackageStartupMessages({library(dplyr); library(fixest)})
source(file.path(CODE, "hist/add_history.R"))
set.seed(123)

# One GLM; returns b/se/p for each requested term. vc = "~cluster" string or "iid". Cached (resume-safe).
fit_keys <- function(f, data, terms, vc, link = "logit", cache = NULL) {
  if (!is.null(cache) && file.exists(cache)) { o <- readRDS(cache); cat("\n---- [cached]", o$txt, sep = "\n"); return(o) }
  m <- feglm(f, data = data, family = binomial(link = link),
             vcov = if (vc == "iid") "iid" else as.formula(vc), notes = FALSE, mem.clean = TRUE, lean = TRUE)
  stopifnot(m$convStatus)
  ct <- coeftable(m); stopifnot(all(terms %in% rownames(ct)))
  txt <- c(paste("----", deparse(f, width.cutoff = 500), "| N =", nobs(m), "| vcov:", vc, "| link:", link),
           capture.output(print(round(ct, 6))))
  cat("\n", txt, sep = "\n")
  out <- list(k = lapply(setNames(terms, terms), function(t) c(b = ct[t, 1], se = ct[t, 2], p = ct[t, 4])), n = nobs(m), txt = txt)
  rm(m); gc()
  if (!is.null(cache)) saveRDS(out, cache)
  out
}

# Common sample: complete on model vars (+cluster), add history anchored on `yr` (the DV year),
# restrict to 1822-2011, complete on all history + benchmark vars.
common_sample <- function(d, f0, extra, c1, c2, yr) {
  d <- d[complete.cases(d[, unique(c(all.vars(f0), extra))]), ]
  n_cc <- nrow(d)
  d <- add_history(d, c1, c2, yr)
  d <- d[d[[yr]] >= 1822 & d[[yr]] <= 2011, ]
  hv <- c("hist_dyad_w15", "hist_mono_w15", "hist_dyad_w26", "hist_mono_w26", "bm_cinc_w15", "bm_democ_w15")
  d <- d[complete.cases(d[, hv]) & is.finite(d$bm_cinc_w15) & is.finite(d$bm_democ_w15), ]
  cat("complete-case N:", n_cc, " common sample N:", nrow(d), " dropped:", n_cc - nrow(d), "\n")
  stopifnot(nrow(d) > 0.8 * n_cc)
  attr(d, "n_cc") <- n_cc
  d
}

run_six <- function(f0, ds, terms, vc, link, W) {
  C <- function(k) file.path(W, paste0("cache_", k, ".rds"))
  list(m0  = fit_keys(f0, ds, terms, vc, link, C("m0")),
       mD  = fit_keys(update(f0, . ~ . + hist_dyad_w15), ds, terms, vc, link, C("mD")),
       mDM = fit_keys(update(f0, . ~ . + hist_dyad_w15 + hist_mono_w15), ds, terms, vc, link, C("mDM")),
       m26 = fit_keys(update(f0, . ~ . + hist_dyad_w26 + hist_mono_w26), ds, terms, vc, link, C("m26")),
       bC  = fit_keys(update(f0, . ~ . + bm_cinc_w15), ds, terms, vc, link, C("bC")),
       bD  = fit_keys(update(f0, . ~ . + bm_democ_w15), ds, terms, vc, link, C("bD")))
}

pct <- function(b, b0) 100 * (b - b0) / abs(b0)
# REPRODUCED if |b_repro-b_pub|/|b_pub| <= .10 and p<.05 agreement (sig_pub = published significance at .05)
repro_status <- function(b_pub, sig_pub, rep) {
  if (is.na(b_pub)) return("UNVERIFIED")
  if (abs(rep["b"] - b_pub) / abs(b_pub) <= 0.10 && sig_pub == (rep["p"] < .05)) "REPRODUCED" else "NOT_REPRODUCED"
}

META <- list()  # xid-level metadata filled by each script
write_p4 <- function(m, fam, role, term, evidence, reproduced, b_pub, b_repro, estimator, n_orig, M, hist_ctrl, notes) {
  g <- function(k, x) M[[k]]$k[[term]][[x]]
  b0 <- g("m0", "b")
  stopifnot(all(sapply(M, function(z) z$n) == M$m0$n))  # one common sample
  row <- data.frame(xid = m$xid, short_cite = m$short_cite, journal = m$journal, year = m$year, doi = m$doi,
    key_family = fam, role = role, key_term = term, key_evidence = evidence, reproduced = reproduced,
    b_pub = b_pub, b_repro = b_repro, estimator = estimator, n_orig = n_orig, n_used = M$m0$n,
    b0 = b0, se0 = g("m0", "se"), p0 = g("m0", "p"),
    b_d = g("mD", "b"), pct_d = pct(g("mD", "b"), b0),
    b_dm = g("mDM", "b"), se_dm = g("mDM", "se"), p_dm = g("mDM", "p"), pct_dm = pct(g("mDM", "b"), b0),
    b_26 = g("m26", "b"), p_26 = g("m26", "p"), pct_26 = pct(g("m26", "b"), b0),
    pct_bC = pct(g("bC", "b"), b0), pct_bD = pct(g("bD", "b"), b0),
    author_has_history_control = hist_ctrl, notes = notes, stringsAsFactors = FALSE)
  dir.create(file.path(OUT, "results_p4"), recursive = TRUE, showWarnings = FALSE)
  write.csv(row, file.path(OUT, "results_p4", paste0(m$xid, "__", term, ".csv")), row.names = FALSE)
  print(t(row[, intersect(c("xid", "id", "key_term", "n_used", "b0", "se0", "p0", "b_dm", "se_dm", "p_dm", "b_26", "p_26", "pct_bC", "pct_bD"), names(row)), drop = FALSE])); invisible(row)  # port: id only, no citation printed
}

# Rows for studies that are not estimated (EXCLUDED / NO_DATA / INFEASIBLE)
write_stub <- function(m, fam, term, evidence, status, estimator, hist_ctrl, notes, suffix = status) {
  cols <- c("b_pub","b_repro","n_orig","n_used","b0","se0","p0","b_d","pct_d","b_dm","se_dm","p_dm","pct_dm","b_26","p_26","pct_26","pct_bC","pct_bD")
  row <- data.frame(xid = m$xid, short_cite = m$short_cite, journal = m$journal, year = m$year, doi = m$doi,
    key_family = fam, role = "KEY", key_term = term, key_evidence = evidence, reproduced = status, stringsAsFactors = FALSE)
  for (cc in cols[1:2]) row[[cc]] <- NA
  row$estimator <- estimator
  for (cc in cols[3:18]) row[[cc]] <- NA
  row$author_has_history_control <- hist_ctrl; row$notes <- notes
  dir.create(file.path(OUT, "results_p4"), recursive = TRUE, showWarnings = FALSE); write.csv(row, file.path(OUT, "results_p4", paste0(m$xid, "__", suffix, ".csv")), row.names = FALSE)  # port: OUT dir
  invisible(row)
}
