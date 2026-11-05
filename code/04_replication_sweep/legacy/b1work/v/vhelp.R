# verify-phase writer: results_verify/<id>__<term>.csv (VERIFY_PROTOCOL step 5)
ROOT <- Sys.getenv("SWEEP_DATA"); OUT <- Sys.getenv("SWEEP_OUT"); CODE <- Sys.getenv("SWEEP_CODE")
stopifnot(nzchar(ROOT), nzchar(OUT), nzchar(CODE))
source(file.path(CODE, "work_p8/round2_rule.R"))   # rebuild_verdict(): round-2 rule with registry alpha_pub (2026-10-11)
pctf <- function(b, b0) 100 * (b - b0) / abs(b0)
verdict_f <- function(b0, p0, pub_b, pub_sig, rebuilt = FALSE) {
  if (is.na(pub_b)) return("NO_PUBLISHED")
  same_sig <- (p0 < .05) == pub_sig; same_sign <- sign(b0) == sign(pub_b)
  if (rebuilt) stop("rebuild verdicts go through rebuild_verdict() (round-2 rule, work_p8/round2_rule.R)")
  if (!same_sign || !same_sig) return("NOT_REPRODUCED")
  if (abs(b0 - pub_b) / abs(pub_b) <= 0.10) "REPRODUCED" else "CLOSE"
}
vwrite <- function(id, cite, fam, term, tm, pub_b, pub_se, pub_n, pub_est, mism, fixed, unfix, R, pub_sig, rebuilt = FALSE, notes = "") {
  g <- function(m, k) if (is.null(R[[m]])) NA else unname(R[[m]][[k]])
  b0 <- g("m0","b"); p0 <- g("m0","p")
  v <- if (rebuilt) rebuild_verdict(id, term, b0, p0, pub_b, pub_sig, pub_se) else verdict_f(b0, p0, pub_b, pub_sig, rebuilt)
  row <- data.frame(id = id, short_cite = cite, key_family = fam, key_term = term, headline_table_model = tm,
    pub_b = pub_b, pub_se = pub_se, pub_n = pub_n, pub_estimator = pub_est, spec_mismatches_before = mism,
    fixed = fixed, unfixable = unfix, verdict = v, b0 = b0, se0 = g("m0","se"), p0 = p0, n_used = g("m0","n"),
    pct_vs_pub = if (is.na(pub_b)) NA else pctf(b0, pub_b),
    b_dm = g("mDM","b"), se_dm = g("mDM","se"), p_dm = g("mDM","p"), pct_dm = pctf(g("mDM","b"), b0),
    b_26 = g("m26","b"), p_26 = g("m26","p"), pct_26 = pctf(g("m26","b"), b0),
    pct_bC = pctf(g("bC","b"), b0), pct_bD = pctf(g("bD","b"), b0), notes = notes, stringsAsFactors = FALSE)
  dir.create(file.path(OUT, "results_verify"), recursive = TRUE, showWarnings = FALSE)
  fn <- file.path(OUT, "results_verify", paste0(id, "__", term, ".csv"))
  write.csv(row, fn, row.names = FALSE); cat("wrote", fn, "\n"); print(t(row[, c("verdict","b0","se0","p0","n_used","pct_vs_pub","b_dm","p_dm","pct_dm","b_26","p_26","pct_bC","pct_bD")]))
  invisible(row)
}
