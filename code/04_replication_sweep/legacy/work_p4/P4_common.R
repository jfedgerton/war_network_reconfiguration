# Phase-4 helpers (batch X005..X093 step 8). Adapted from work_demcap/D_common.R.
ROOT <- Sys.getenv("SWEEP_DATA"); OUT <- Sys.getenv("SWEEP_OUT"); CODE <- Sys.getenv("SWEEP_CODE")
stopifnot(nzchar(ROOT), nzchar(OUT), nzchar(CODE))
suppressPackageStartupMessages({library(dplyr); library(fixest)})
source(file.path(CODE, "hist/add_history.R"))
set.seed(123)
META <- read.csv(file.path(ROOT, "phase4_studies.csv"), stringsAsFactors = FALSE)

# vc: "iid", "hetero", or a cluster variable name. fam: binomial link.
fit_keys <- function(f, data, terms, vc, cache = NULL, fam = binomial("logit")) {
  if (!is.null(cache) && file.exists(cache)) { o <- readRDS(cache); cat("\n---- [cached]", o$txt, sep = "\n"); return(o) }
  V <- if (vc %in% c("iid", "hetero")) vc else as.formula(paste0("~", vc))
  m <- feglm(f, data = data, family = fam, vcov = V, notes = FALSE, mem.clean = TRUE, lean = TRUE)
  stopifnot(m$convStatus)
  ct <- coeftable(m); stopifnot(all(terms %in% rownames(ct)))
  txt <- c(paste("----", paste(deparse(f, width.cutoff = 500), collapse = " "), "| N =", nobs(m), "| vcov:", vc),
           capture.output(print(round(ct, 6))))
  cat("\n", txt, sep = "\n")
  out <- list(k = lapply(setNames(terms, terms), function(t) c(b = ct[t, 1], se = ct[t, 2], p = ct[t, 4])), n = nobs(m), txt = txt)
  rm(m); gc()
  if (!is.null(cache)) saveRDS(out, cache)
  out
}

# complete on model vars (+cluster), add history, restrict 1822-2011, complete on history+benchmarks
common_sample <- function(d, f0, extra, c1, c2, yr, min_share = 0.8) {
  d <- d[complete.cases(d[, unique(c(all.vars(f0), extra))]), ]
  n_cc <- nrow(d)
  d <- add_history(d, c1, c2, yr)
  d <- d[d[[yr]] >= 1822 & d[[yr]] <= 2011, ]
  hv <- c("hist_dyad_w15", "hist_mono_w15", "hist_dyad_w26", "hist_mono_w26", "bm_cinc_w15", "bm_democ_w15")
  d <- d[complete.cases(d[, hv]) & is.finite(d$bm_cinc_w15) & is.finite(d$bm_democ_w15), ]
  cat("complete-case N:", n_cc, " common sample N:", nrow(d), "\n")
  stopifnot(nrow(d) > min_share * n_cc)
  d
}

run_six <- function(f0, ds, terms, vc, W, tag = "", fam = binomial("logit")) {
  C <- function(k) file.path(W, paste0("cache_", tag, k, ".rds"))
  list(m0  = fit_keys(f0, ds, terms, vc, C("m0"), fam),
       mD  = fit_keys(update(f0, . ~ . + hist_dyad_w15), ds, terms, vc, C("mD"), fam),
       mDM = fit_keys(update(f0, . ~ . + hist_dyad_w15 + hist_mono_w15), ds, terms, vc, C("mDM"), fam),
       m26 = fit_keys(update(f0, . ~ . + hist_dyad_w26 + hist_mono_w26), ds, terms, vc, C("m26"), fam),
       bC  = fit_keys(update(f0, . ~ . + bm_cinc_w15), ds, terms, vc, C("bC"), fam),
       bD  = fit_keys(update(f0, . ~ . + bm_democ_w15), ds, terms, vc, C("bD"), fam))
}

pct <- function(b, b0) 100 * (b - b0) / abs(b0)
repro_status <- function(b_pub, sig_pub, rep) {
  if (is.na(b_pub)) return("UNVERIFIED")
  if (abs(rep[["b"]] - b_pub) / abs(b_pub) <= 0.10 && sig_pub == (rep[["p"]] < .05)) "REPRODUCED" else "NOT_REPRODUCED"
}

COLS <- c("xid","short_cite","journal","year","doi","key_family","role","key_term","key_evidence","reproduced","b_pub","b_repro",
  "estimator","n_orig","n_used","b0","se0","p0","b_d","pct_d","b_dm","se_dm","p_dm","pct_dm","b_26","p_26","pct_26","pct_bC","pct_bD",
  "author_has_history_control","notes")

write_p4 <- function(xid, short_cite, fam, role, term, evidence, reproduced, b_pub, b_repro, estimator, n_orig, M, hist_ctrl, notes) {
  mt <- META[META$xid == xid, ]; stopifnot(nrow(mt) == 1)
  g <- function(m, x) M[[m]]$k[[term]][[x]]
  stopifnot(all(sapply(M, function(m) m$n) == M$m0$n))  # one common sample
  b0 <- g("m0", "b")
  row <- data.frame(xid = xid, short_cite = short_cite, journal = mt$journal, year = mt$year, doi = mt$doi, key_family = fam, role = role,
    key_term = term, key_evidence = evidence, reproduced = reproduced, b_pub = as.numeric(b_pub), b_repro = b_repro, estimator = estimator,
    n_orig = n_orig, n_used = M$m0$n, b0 = b0, se0 = g("m0", "se"), p0 = g("m0", "p"),
    b_d = g("mD", "b"), pct_d = pct(g("mD", "b"), b0),
    b_dm = g("mDM", "b"), se_dm = g("mDM", "se"), p_dm = g("mDM", "p"), pct_dm = pct(g("mDM", "b"), b0),
    b_26 = g("m26", "b"), p_26 = g("m26", "p"), pct_26 = pct(g("m26", "b"), b0),
    pct_bC = pct(g("bC", "b"), b0), pct_bD = pct(g("bD", "b"), b0),
    author_has_history_control = hist_ctrl, notes = notes, stringsAsFactors = FALSE)
  stopifnot(identical(names(row), COLS))
  dir.create(file.path(OUT, "results_p4"), recursive = TRUE, showWarnings = FALSE); write.csv(row, file.path(OUT, "results_p4", paste0(xid, "__", term, ".csv")), row.names = FALSE)  # port: OUT dir
  print(t(row[, intersect(c("xid", "id", "key_term", "n_used", "b0", "se0", "p0", "b_dm", "se_dm", "p_dm", "b_26", "p_26", "pct_bC", "pct_bD"), names(row)), drop = FALSE])); invisible(row)  # port: id only, no citation printed
}
