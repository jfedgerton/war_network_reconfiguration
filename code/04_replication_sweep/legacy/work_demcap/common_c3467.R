# Shared helpers for demcap reruns C3 C4 C6 C7 (DEMCAP_PROTOCOL.md + RERUN_PROTOCOL.md).
ROOT <- Sys.getenv("SWEEP_DATA"); OUT <- Sys.getenv("SWEEP_OUT"); CODE <- Sys.getenv("SWEEP_CODE")
stopifnot(nzchar(ROOT), nzchar(OUT), nzchar(CODE))
suppressPackageStartupMessages({library(dplyr); library(fixest)})
source(file.path(CODE, "hist/add_history.R"))
set.seed(123)
LG <- binomial(link = "logit")

# Fit one logit (dyad-clustered SE), return b/se/p for each key term. Cached (resume-safe on shared box).
fit_keys <- function(f, data, terms, cluster, cache = NULL) {
  if (!is.null(cache) && file.exists(cache)) { o <- readRDS(cache); cat("\n---- [cached]", o$txt, sep = "\n"); return(o) }
  m <- feglm(f, data = data, family = LG, vcov = as.formula(paste0("~", cluster)), notes = FALSE, mem.clean = TRUE, lean = TRUE)
  stopifnot(m$convStatus)
  ct <- coeftable(m); stopifnot(all(terms %in% rownames(ct)))
  txt <- c(paste("----", deparse(f, width.cutoff = 500), "| N =", nobs(m), "| cluster:", cluster), capture.output(print(round(ct, 5))))
  cat("\n", txt, sep = "\n")
  out <- list(k = lapply(setNames(terms, terms), function(t) c(b = ct[t, 1], se = ct[t, 2], p = ct[t, 4])), n = nobs(m), txt = txt)
  rm(m); gc()
  if (!is.null(cache)) saveRDS(out, cache)
  out
}

# Build the common sample: complete on model vars + cluster, add history, restrict 1822-2011,
# complete on all history AND benchmark vars (one common sample for m0/mD/mDM/m26/bC/bD).
common_sample <- function(d, f0, cl, c1, c2, yr) {
  d <- d[complete.cases(d[, c(all.vars(f0), cl)]), ]
  n_cc <- nrow(d)
  d <- add_history(d, c1, c2, yr)
  d <- d[d[[yr]] >= 1822 & d[[yr]] <= 2011, ]
  hv <- c("hist_dyad_w15", "hist_mono_w15", "hist_dyad_w26", "hist_mono_w26", "bm_cinc_w15", "bm_democ_w15")
  d <- d[complete.cases(d[, hv]) & is.finite(d$bm_cinc_w15) & is.finite(d$bm_democ_w15), ]
  cat("complete-case N:", n_cc, " common sample N:", nrow(d), "\n")
  stopifnot(nrow(d) > 0.8 * n_cc)
  d
}

run_six <- function(f0, ds, terms, cl, W) {
  C <- function(k) file.path(W, paste0("cache_", k, ".rds"))
  list(m0  = fit_keys(f0, ds, terms, cl, C("m0")),
       mD  = fit_keys(update(f0, . ~ . + hist_dyad_w15), ds, terms, cl, C("mD")),
       mDM = fit_keys(update(f0, . ~ . + hist_dyad_w15 + hist_mono_w15), ds, terms, cl, C("mDM")),
       m26 = fit_keys(update(f0, . ~ . + hist_dyad_w26 + hist_mono_w26), ds, terms, cl, C("m26")),
       bC  = fit_keys(update(f0, . ~ . + bm_cinc_w15), ds, terms, cl, C("bC")),
       bD  = fit_keys(update(f0, . ~ . + bm_democ_w15), ds, terms, cl, C("bD")))
}

pct <- function(b, b0) 100 * (b - b0) / abs(b0)
repro_status <- function(b_pub, rep) {
  if (is.na(b_pub)) return("UNVERIFIED")
  ok <- abs(rep["b"] - b_pub) / abs(b_pub) <= 0.10
  # significance agreement is judged against the published/log p (passed as attr)
  if (ok && isTRUE(attr(b_pub, "sig") == (rep["p"] < .05))) "REPRODUCED" else "NOT_REPRODUCED"
}

write_dc <- function(study_id, short_cite, fam, term, evidence, reproduced, b_pub, b_repro, estimator, n_orig, M, notes) {
  g <- function(m, x) M[[m]]$k[[term]][[x]]
  b0 <- g("m0", "b")
  row <- data.frame(study_id = study_id, short_cite = short_cite, key_family = fam, key_term = term, key_evidence = evidence,
    reproduced = reproduced, b_pub = as.numeric(b_pub), b_repro = b_repro, estimator = estimator, n_orig = n_orig, n_used = M$m0$n,
    b0 = b0, se0 = g("m0", "se"), p0 = g("m0", "p"),
    b_d = g("mD", "b"), pct_d = pct(g("mD", "b"), b0),
    b_dm = g("mDM", "b"), se_dm = g("mDM", "se"), p_dm = g("mDM", "p"), pct_dm = pct(g("mDM", "b"), b0),
    b_26 = g("m26", "b"), p_26 = g("m26", "p"), pct_26 = pct(g("m26", "b"), b0),
    pct_bC = pct(g("bC", "b"), b0), pct_bD = pct(g("bD", "b"), b0), notes = notes, stringsAsFactors = FALSE)
  stopifnot(all(sapply(M, function(m) m$n) == M$m0$n))  # one common sample
  dir.create(file.path(OUT, "results_demcap"), recursive = TRUE, showWarnings = FALSE)
  write.csv(row, file.path(OUT, "results_demcap", paste0(study_id, "__", term, ".csv")), row.names = FALSE)
  print(t(row[, c("study_id", "key_term", "n_used", "b0", "se0", "p0", "b_dm", "p_dm", "b_26", "p_26", "pct_bC", "pct_bD")])); invisible(row)  # port: id only, no citation
}
