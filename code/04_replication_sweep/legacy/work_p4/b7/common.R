# Phase-4 helpers for batch b7 (X007 X015 X023 X031 X039 X047 X055 X063 X071 X079 X087)
suppressPackageStartupMessages({library(dplyr); library(fixest)})
ROOT <- Sys.getenv("SWEEP_DATA"); OUT <- Sys.getenv("SWEEP_OUT"); CODE <- Sys.getenv("SWEEP_CODE")
stopifnot(nzchar(ROOT), nzchar(OUT), nzchar(CODE))
source(file.path(CODE, "hist/add_history.R"))
set.seed(123)
HV <- c("hist_dyad_w15","hist_mono_w15","hist_dyad_w26","hist_mono_w26","bm_cinc_w15","bm_democ_w15")
ADD <- list(m0 = character(0), mD = "hist_dyad_w15", mDM = c("hist_dyad_w15","hist_mono_w15"),
            m26 = c("hist_dyad_w26","hist_mono_w26"), bC = "bm_cinc_w15", bD = "bm_democ_w15")

# Binary GLM, Stata-style SEs: iid (vc=NULL) or clustered (G/(G-1) only)
fit_bin <- function(rhs, dv, data, link = "logit", cluster = NULL) {
  f <- as.formula(paste(dv, "~", paste(rhs, collapse = " + ")))
  vc <- if (is.null(cluster)) "iid" else as.formula(paste0("~", cluster))
  feglm(f, data = data, family = binomial(link = link), vcov = vc,
        ssc = ssc(adj = FALSE, cluster.adj = TRUE), lean = TRUE, notes = FALSE, mem.clean = TRUE)
}
kstat <- function(ct, term) { stopifnot(term %in% rownames(ct)); c(b = unname(ct[term, 1]), se = unname(ct[term, 2]), p = unname(ct[term, 4])) }

# Common sample: complete on model vars, add history, restrict 1822-2011, complete on history+benchmarks
common_sample <- function(d, vars, c1, c2, yr) {
  d <- d[complete.cases(as.data.frame(d)[, vars]), ]
  n_cc <- nrow(d)
  d <- add_history(d, c1, c2, yr)
  d <- d[d[[yr]] >= 1822 & d[[yr]] <= 2011, ]
  d <- d[complete.cases(as.data.frame(d)[, HV]) & is.finite(d$bm_cinc_w15) & is.finite(d$bm_democ_w15), ]
  cat(sprintf("complete-case N: %d  common sample N: %d (dropped %d)\n", n_cc, nrow(d), n_cc - nrow(d)))
  stopifnot(nrow(d) > 0.5 * n_cc)
  attr(d, "n_cc") <- n_cc
  d
}

# fitter(rhs_vector, data) -> coeftable; runs six models, returns list of key stats per term
run_six <- function(d, base_rhs, fitter, terms, logf) {
  out <- list(); ns <- c()
  for (nm in names(ADD)) {
    r <- fitter(c(base_rhs, ADD[[nm]]), d)
    ct <- r$ct; ns[nm] <- r$n
    cat("\n==========", nm, " N =", r$n, "==========\n", file = logf, append = TRUE)
    capture.output(print(round(ct, 5)), file = logf, append = TRUE)
    out[[nm]] <- lapply(setNames(terms, terms), function(t) kstat(ct, t))
    cat(sprintf("[%s] N=%d %s\n", nm, r$n, paste(sprintf("%s b=%.4g p=%.3g", terms, sapply(out[[nm]], `[`, "b"), sapply(out[[nm]], `[`, "p")), collapse = "; ")))
    invisible(gc())
  }
  stopifnot(length(unique(ns)) == 1)
  attr(out, "n") <- unname(ns[1])
  out
}

P4COLS <- c("xid","short_cite","journal","year","doi","key_family","role","key_term","key_evidence","reproduced","b_pub","b_repro",
            "estimator","n_orig","n_used","b0","se0","p0","b_d","pct_d","b_dm","se_dm","p_dm","pct_dm","b_26","p_26","pct_26",
            "pct_bC","pct_bD","author_has_history_control","notes")
pct <- function(b, b0) 100 * (b - b0) / abs(b0)

write_p4 <- function(meta, term, fam, role, evidence, reproduced, b_pub, b_repro, estimator, n_orig, R, ahc, notes, fname_term = term) {
  g <- function(m, x) unname(R[[m]][[term]][x])
  b0 <- g("m0", "b")
  row <- data.frame(xid = meta$xid, short_cite = meta$cite, journal = meta$journal, year = meta$year, doi = meta$doi,
    key_family = fam, role = role, key_term = term, key_evidence = evidence, reproduced = reproduced,
    b_pub = b_pub, b_repro = b_repro, estimator = estimator, n_orig = n_orig, n_used = attr(R, "n"),
    b0 = b0, se0 = g("m0", "se"), p0 = g("m0", "p"),
    b_d = g("mD", "b"), pct_d = pct(g("mD", "b"), b0),
    b_dm = g("mDM", "b"), se_dm = g("mDM", "se"), p_dm = g("mDM", "p"), pct_dm = pct(g("mDM", "b"), b0),
    b_26 = g("m26", "b"), p_26 = g("m26", "p"), pct_26 = pct(g("m26", "b"), b0),
    pct_bC = pct(g("bC", "b"), b0), pct_bD = pct(g("bD", "b"), b0),
    author_has_history_control = ahc, notes = notes, stringsAsFactors = FALSE)
  stopifnot(identical(names(row), P4COLS))
  dir.create(file.path(OUT, "results_p4"), recursive = TRUE, showWarnings = FALSE); write.csv(row, file.path(OUT, "results_p4", paste0(meta$xid, "__", fname_term, ".csv")), row.names = FALSE)  # port: OUT dir
  print(t(row[, intersect(c("xid", "id", "key_term", "n_used", "b0", "se0", "p0", "b_dm", "se_dm", "p_dm", "b_26", "p_26", "pct_bC", "pct_bD"), names(row)), drop = FALSE])); invisible(row)  # port: id only, no citation printed
}

write_p4_stub <- function(meta, term, fam, role, evidence, status, estimator, ahc, notes, fname_term = term, b_pub = NA) {
  row <- as.data.frame(setNames(as.list(rep(NA, length(P4COLS))), P4COLS))
  row$xid <- meta$xid; row$short_cite <- meta$cite; row$journal <- meta$journal; row$year <- meta$year; row$doi <- meta$doi
  row$key_family <- fam; row$role <- role; row$key_term <- term; row$key_evidence <- evidence; row$reproduced <- status
  row$estimator <- estimator; row$author_has_history_control <- ahc; row$notes <- notes; row$b_pub <- b_pub
  dir.create(file.path(OUT, "results_p4"), recursive = TRUE, showWarnings = FALSE); write.csv(row, file.path(OUT, "results_p4", paste0(meta$xid, "__", fname_term, ".csv")), row.names = FALSE)  # port: OUT dir
  print(t(row[, intersect(c("xid", "id", "key_term", "n_used", "b0", "se0", "p0", "b_dm", "se_dm", "p_dm", "b_26", "p_26", "pct_bC", "pct_bD"), names(row)), drop = FALSE])); invisible(row)  # port: id only, no citation printed
}
