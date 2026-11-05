# Phase-4 helpers for batch xb1 (X001 X009 X017 X025 X033 ...). Sequential, fail-loud.
ROOT <- Sys.getenv("SWEEP_DATA"); OUT <- Sys.getenv("SWEEP_OUT"); CODE <- Sys.getenv("SWEEP_CODE")
stopifnot(nzchar(ROOT), nzchar(OUT), nzchar(CODE))
suppressPackageStartupMessages({library(dplyr); library(fixest)})
source(file.path(CODE, "hist/add_history.R"))
set.seed(123)
HV <- c("hist_dyad_w15", "hist_mono_w15", "hist_dyad_w26", "hist_mono_w26", "bm_cinc_w15", "bm_democ_w15")

# one GLM (logit/probit), clustered SE; returns b/se/p for each term; cached (resume-safe)
fit_keys <- function(f, data, terms, cluster, link = "logit", cache = NULL) {
  if (!is.null(cache) && file.exists(cache)) { o <- readRDS(cache); cat("\n---- [cached]", o$txt, sep = "\n"); return(o) }
  m <- feglm(f, data = data, family = binomial(link = link), vcov = as.formula(paste0("~", cluster)),
             notes = FALSE, mem.clean = TRUE, lean = TRUE)
  stopifnot(m$convStatus)
  ct <- coeftable(m); stopifnot(all(terms %in% rownames(ct)))
  txt <- c(paste("----", deparse(f, width.cutoff = 500), "| N =", nobs(m), "| cluster:", cluster, "| link:", link),
           capture.output(print(round(ct, 6))))
  cat("\n", txt, sep = "\n")
  out <- list(k = lapply(setNames(terms, terms), function(t) c(b = ct[t, 1], se = ct[t, 2], p = ct[t, 4])), n = nobs(m), txt = txt)
  rm(m); gc()
  if (!is.null(cache)) saveRDS(out, cache)
  out
}

# common sample: complete on model vars + cluster; add history (hist_year = DV year); restrict 1822-2011; complete on HV
common_sample <- function(d, f0, cl, c1, c2, hist_year) {
  d <- d[complete.cases(d[, unique(c(all.vars(f0), cl, c1, c2, hist_year))]), ]
  n_cc <- nrow(d)
  d <- add_history(d, c1, c2, hist_year)
  d <- d[d[[hist_year]] >= 1822 & d[[hist_year]] <= 2011, ]
  d <- d[complete.cases(d[, HV]) & is.finite(d$bm_cinc_w15) & is.finite(d$bm_democ_w15), ]
  cat("complete-case N:", n_cc, " common sample N:", nrow(d), " dropped:", n_cc - nrow(d), "\n")
  print(summary(d[, HV]))
  stopifnot(nrow(d) > 0.8 * n_cc, sd(d$hist_dyad_w15) > 0)
  attr(d, "n_cc") <- n_cc
  d
}

run_six <- function(f0, ds, terms, cl, W, link = "logit") {
  C <- function(k) file.path(W, paste0("cache_", k, ".rds"))
  list(m0  = fit_keys(f0, ds, terms, cl, link, C("m0")),
       mD  = fit_keys(update(f0, . ~ . + hist_dyad_w15), ds, terms, cl, link, C("mD")),
       mDM = fit_keys(update(f0, . ~ . + hist_dyad_w15 + hist_mono_w15), ds, terms, cl, link, C("mDM")),
       m26 = fit_keys(update(f0, . ~ . + hist_dyad_w26 + hist_mono_w26), ds, terms, cl, link, C("m26")),
       bC  = fit_keys(update(f0, . ~ . + bm_cinc_w15), ds, terms, cl, link, C("bC")),
       bD  = fit_keys(update(f0, . ~ . + bm_democ_w15), ds, terms, cl, link, C("bD")))
}

pct <- function(b, b0) 100 * (b - b0) / abs(b0)
# REPRODUCED if |b_repro-b_pub|/|b_pub| <= .10 and significance (p<.05) agrees
repro <- function(b_pub, p_pub_sig, b_rep, p_rep) {
  if (is.na(b_pub)) return("UNVERIFIED")
  if (abs(b_rep - b_pub) / abs(b_pub) <= 0.10 && (p_rep < .05) == p_pub_sig) "REPRODUCED" else "NOT_REPRODUCED"
}

META <- list(
  X001 = list(short_cite = "Reiter & Stam 2003 APSR", journal = "APSR", year = 2003, doi = "10.1017/s0003055403000704"),
  X009 = list(short_cite = "Gibler & Braithwaite 2013 BJPS", journal = "BJPS", year = 2012, doi = "10.1017/s000712341200052x"),
  X017 = list(short_cite = "Gartzke & Weisiger 2013 ISQ", journal = "ISQ", year = 2013, doi = "10.1111/isqu.12046"),
  X025 = list(short_cite = "Gelpi & Grieco 2001 JCR", journal = "JCR", year = 2001, doi = "10.1177/0022002701045006005"),
  X033 = list(short_cite = "Rasler & Thompson 2001 JPR", journal = "JPR", year = 2001, doi = "10.1177/0022343301038006002"))

write_p4 <- function(xid, fam, role, term, evidence, reproduced, b_pub, b_repro, estimator, n_orig, M, hist_ctrl, notes) {
  g <- function(m, x) M[[m]]$k[[term]][[x]]
  stopifnot(all(sapply(M, function(m) m$n) == M$m0$n))  # one common sample
  b0 <- g("m0", "b"); mt <- META[[xid]]
  row <- data.frame(xid = xid, short_cite = mt$short_cite, journal = mt$journal, year = mt$year, doi = mt$doi,
    key_family = fam, role = role, key_term = term, key_evidence = evidence, reproduced = reproduced,
    b_pub = as.numeric(b_pub), b_repro = as.numeric(b_repro), estimator = estimator, n_orig = n_orig, n_used = M$m0$n,
    b0 = b0, se0 = g("m0", "se"), p0 = g("m0", "p"),
    b_d = g("mD", "b"), pct_d = pct(g("mD", "b"), b0),
    b_dm = g("mDM", "b"), se_dm = g("mDM", "se"), p_dm = g("mDM", "p"), pct_dm = pct(g("mDM", "b"), b0),
    b_26 = g("m26", "b"), p_26 = g("m26", "p"), pct_26 = pct(g("m26", "b"), b0),
    pct_bC = pct(g("bC", "b"), b0), pct_bD = pct(g("bD", "b"), b0),
    author_has_history_control = hist_ctrl, notes = notes, stringsAsFactors = FALSE)
  dir.create(file.path(OUT, "results_p4"), recursive = TRUE, showWarnings = FALSE); write.csv(row, file.path(OUT, "results_p4", paste0(xid, "__", term, ".csv")), row.names = FALSE)  # port: OUT dir
  print(t(row[, intersect(c("xid", "id", "key_term", "n_used", "b0", "se0", "p0", "b_dm", "se_dm", "p_dm", "b_26", "p_26", "pct_bC", "pct_bD"), names(row)), drop = FALSE])); invisible(row)  # port: id only, no citation printed
}

# rows for studies with no estimates
write_stub <- function(xid, short_cite, journal, year, doi, fam, term, evidence, status, notes) {
  cols <- c("xid","short_cite","journal","year","doi","key_family","role","key_term","key_evidence","reproduced","b_pub","b_repro",
            "estimator","n_orig","n_used","b0","se0","p0","b_d","pct_d","b_dm","se_dm","p_dm","pct_dm","b_26","p_26","pct_26",
            "pct_bC","pct_bD","author_has_history_control","notes")
  row <- as.data.frame(setNames(as.list(rep(NA, length(cols))), cols))
  row$xid <- xid; row$short_cite <- short_cite; row$journal <- journal; row$year <- year; row$doi <- doi
  row$key_family <- fam; row$role <- "KEY"; row$key_term <- term; row$key_evidence <- evidence; row$reproduced <- status; row$notes <- notes
  fn <- if (status == "EXCLUDED") paste0(xid, "__EXCLUDED.csv") else paste0(xid, "__", term, ".csv")
  dir.create(file.path(OUT, "results_p4"), recursive = TRUE, showWarnings = FALSE); write.csv(row, file.path(OUT, "results_p4", fn), row.names = FALSE, na = "")  # port: OUT dir
}
