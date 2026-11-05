# Phase-4 helpers (FROZEN COPY of common_p4_b04.R for batch X002..X090 step 8). Generic fitter + common sample + six-model battery + writer.
ROOT <- Sys.getenv("SWEEP_DATA"); OUT <- Sys.getenv("SWEEP_OUT"); CODE <- Sys.getenv("SWEEP_CODE")
stopifnot(nzchar(ROOT), nzchar(OUT), nzchar(CODE))
suppressPackageStartupMessages({library(dplyr); library(fixest); library(sandwich); library(lmtest)})
source(file.path(CODE, "hist/add_history.R"))
set.seed(123)
HV <- c("hist_dyad_w15", "hist_mono_w15", "hist_dyad_w26", "hist_mono_w26", "bm_cinc_w15", "bm_democ_w15")

# fitfun(f, data) must return a coefficient matrix with columns b, se, p (rownames = terms)
fit_keys <- function(f, data, terms, fitfun, cache = NULL) {
  if (!is.null(cache) && file.exists(cache)) { o <- readRDS(cache); cat("\n---- [cached]", o$txt, sep = "\n"); return(o) }
  ct <- fitfun(f, data)
  stopifnot(all(terms %in% rownames(ct)), all(is.finite(ct[terms, 1])))
  txt <- c(paste("----", paste(deparse(f, width.cutoff = 500), collapse = " "), "| N =", attr(ct, "n")), capture.output(print(round(ct, 5))))
  cat("\n", txt, sep = "\n")
  out <- list(k = lapply(setNames(terms, terms), function(t) c(b = ct[t, 1], se = ct[t, 2], p = ct[t, 3])), n = attr(ct, "n"), txt = txt)
  gc()
  if (!is.null(cache)) saveRDS(out, cache)
  out
}

# common sample: complete on model vars (+extra), add history, restrict 1822-2011, complete on history+benchmarks
common_sample <- function(d, f0, extra, c1, c2, yr) {
  d <- d[complete.cases(d[, unique(c(all.vars(f0), extra))]), ]
  n_cc <- nrow(d)
  d <- add_history(d, c1, c2, yr)
  d <- d[d[[yr]] >= 1822 & d[[yr]] <= 2011, ]
  d <- d[complete.cases(d[, HV]) & is.finite(d$bm_cinc_w15) & is.finite(d$bm_democ_w15), ]
  cat("complete-case N:", n_cc, " common sample N:", nrow(d), " dropped:", n_cc - nrow(d), "\n")
  stopifnot(nrow(d) > 0.8 * n_cc)
  d
}

run_six <- function(f0, ds, terms, fitfun, W) {
  C <- function(k) file.path(W, paste0("cache_", k, ".rds"))
  list(m0  = fit_keys(f0, ds, terms, fitfun, C("m0")),
       mD  = fit_keys(update(f0, . ~ . + hist_dyad_w15), ds, terms, fitfun, C("mD")),
       mDM = fit_keys(update(f0, . ~ . + hist_dyad_w15 + hist_mono_w15), ds, terms, fitfun, C("mDM")),
       m26 = fit_keys(update(f0, . ~ . + hist_dyad_w26 + hist_mono_w26), ds, terms, fitfun, C("m26")),
       bC  = fit_keys(update(f0, . ~ . + bm_cinc_w15), ds, terms, fitfun, C("bC")),
       bD  = fit_keys(update(f0, . ~ . + bm_democ_w15), ds, terms, fitfun, C("bD")))
}

# ---- fitters ----
ff_logit_cl <- function(cl) function(f, data) {           # logit, cluster-robust (fixest, Stata-like small-sample adj.)
  m <- feglm(f, data = data, family = binomial("logit"), vcov = as.formula(paste0("~", cl)), notes = FALSE, mem.clean = TRUE)
  stopifnot(m$convStatus); ct <- coeftable(m)[, c(1, 2, 4), drop = FALSE]; attr(ct, "n") <- nobs(m); ct
}
ff_logit_hc1 <- function(f, data) {                        # logit, Stata ", robust" (HC1)
  m <- feglm(f, data = data, family = binomial("logit"), vcov = "hetero", notes = FALSE, mem.clean = TRUE)
  stopifnot(m$convStatus); ct <- coeftable(m)[, c(1, 2, 4), drop = FALSE]; attr(ct, "n") <- nobs(m); ct
}
ff_logit_iid <- function(f, data) {                        # plain logit (ML SEs)
  m <- glm(f, data = data, family = binomial("logit")); stopifnot(m$converged)
  ct <- summary(m)$coefficients[, c(1, 2, 4), drop = FALSE]; attr(ct, "n") <- nobs(m); ct
}
ff_ols_iid <- function(f, data) {
  m <- feols(f, data = data, vcov = "iid", notes = FALSE); ct <- coeftable(m)[, c(1, 2, 4), drop = FALSE]; attr(ct, "n") <- nobs(m); ct
}

pct <- function(b, b0) 100 * (b - b0) / abs(b0)

# ---- writer: one row per term ----
P4COLS <- c("xid","short_cite","journal","year","doi","key_family","role","key_term","key_evidence","reproduced","b_pub","b_repro",
  "estimator","n_orig","n_used","b0","se0","p0","b_d","pct_d","b_dm","se_dm","p_dm","pct_dm","b_26","p_26","pct_26","pct_bC","pct_bD",
  "author_has_history_control","notes")
write_p4 <- function(meta, term, fam, role, evidence, reproduced, b_pub, b_repro, estimator, n_orig, M, hist_ctrl, notes, fname = term) {
  g <- function(m, x) if (is.null(M[[m]])) NA_real_ else unname(M[[m]]$k[[term]][[x]])
  b0 <- g("m0", "b")
  row <- data.frame(xid = meta$xid, short_cite = meta$short_cite, journal = meta$journal, year = meta$year, doi = meta$doi,
    key_family = fam, role = role, key_term = term, key_evidence = evidence, reproduced = reproduced,
    b_pub = as.numeric(b_pub), b_repro = as.numeric(b_repro), estimator = estimator, n_orig = n_orig, n_used = M$m0$n,
    b0 = b0, se0 = g("m0", "se"), p0 = g("m0", "p"),
    b_d = g("mD", "b"), pct_d = pct(g("mD", "b"), b0),
    b_dm = g("mDM", "b"), se_dm = g("mDM", "se"), p_dm = g("mDM", "p"), pct_dm = pct(g("mDM", "b"), b0),
    b_26 = g("m26", "b"), p_26 = g("m26", "p"), pct_26 = pct(g("m26", "b"), b0),
    pct_bC = pct(g("bC", "b"), b0), pct_bD = pct(g("bD", "b"), b0),
    author_has_history_control = hist_ctrl, notes = notes, stringsAsFactors = FALSE)
  ns <- sapply(Filter(Negate(is.null), M), function(m) m$n); stopifnot(all(ns == M$m0$n))   # one common sample
  stopifnot(identical(names(row), P4COLS))
  dir.create(file.path(OUT, "results_p4"), recursive = TRUE, showWarnings = FALSE); write.csv(row, file.path(OUT, "results_p4", paste0(meta$xid, "__", gsub("[^A-Za-z0-9_.]", "_", fname), ".csv")), row.names = FALSE)  # port: OUT dir
  print(t(row[, intersect(c("xid", "id", "key_term", "n_used", "b0", "se0", "p0", "b_dm", "se_dm", "p_dm", "b_26", "p_26", "pct_bC", "pct_bD"), names(row)), drop = FALSE])); invisible(row)  # port: id only, no citation printed
}

# stub row for EXCLUDED / NO_DATA / INFEASIBLE
write_stub <- function(meta, status, fam, term, evidence, notes, fname = status) {
  row <- as.data.frame(setNames(as.list(rep(NA, length(P4COLS))), P4COLS))
  row$xid <- meta$xid; row$short_cite <- meta$short_cite; row$journal <- meta$journal; row$year <- meta$year; row$doi <- meta$doi
  row$key_family <- fam; row$role <- "KEY"; row$key_term <- term; row$key_evidence <- evidence; row$reproduced <- status; row$notes <- notes
  dir.create(file.path(OUT, "results_p4"), recursive = TRUE, showWarnings = FALSE); write.csv(row, file.path(OUT, "results_p4", paste0(meta$xid, "__", fname, ".csv")), row.names = FALSE)  # port: OUT dir
  invisible(row)
}

# ---- extra fitters for batch X002..X090 ----
ff_logit_cl_fe <- function(cl, fe) function(f, data) {    # logit with absorbed FE (e.g. year FE), cluster-robust
  f2 <- as.formula(paste(paste(deparse(f, width.cutoff = 500), collapse = " "), "|", fe))
  m <- feglm(f2, data = data, family = binomial("logit"), vcov = as.formula(paste0("~", cl)), notes = FALSE, mem.clean = TRUE)
  stopifnot(m$convStatus); ct <- coeftable(m)[, c(1, 2, 4), drop = FALSE]; attr(ct, "n") <- nobs(m); ct
}
ff_felogit_iid <- function(fe) function(f, data) {        # unconditional FE logit (dyad FE), ML SEs (approximates Stata clogit)
  f2 <- as.formula(paste(paste(deparse(f, width.cutoff = 500), collapse = " "), "|", fe))
  m <- feglm(f2, data = data, family = binomial("logit"), vcov = "iid", notes = FALSE, mem.clean = TRUE)
  stopifnot(m$convStatus); ct <- coeftable(m)[, c(1, 2, 4), drop = FALSE]; attr(ct, "n") <- nobs(m); ct
}
ff_gee_gauss_exch <- function(id) function(f, data) {     # Stata: xtgee y x, robust (gaussian/identity/exchangeable), semi-robust SE
  suppressPackageStartupMessages(library(geepack))
  data <- data[complete.cases(data[, c(all.vars(f), id)]), ]
  data <- data[order(data[[id]], data$.ord), ]
  data$.id <- data[[id]]
  m <- geeglm(f, data = data, id = .id, family = gaussian, corstr = "exchangeable", std.err = "san.se")
  ct <- summary(m)$coefficients[, c(1, 2, 4), drop = FALSE]; attr(ct, "n") <- nrow(data)
  cat("   [gee alpha =", round(unlist(m$geese$alpha), 4), "]\n"); rm(m); gc(); ct
}
