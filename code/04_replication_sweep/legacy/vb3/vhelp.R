# Verify-phase helper (batch 3): write results_verify/<id>__<term>.csv per VERIFY_PROTOCOL step 5.
ROOT <- Sys.getenv("SWEEP_DATA"); OUT <- Sys.getenv("SWEEP_OUT"); CODE <- Sys.getenv("SWEEP_CODE")
stopifnot(nzchar(ROOT), nzchar(OUT), nzchar(CODE))
dir.create(file.path(OUT, "results_verify"), recursive = TRUE, showWarnings = FALSE)
source(file.path(CODE, "work_p8/round2_rule.R"))   # rebuild_verdict(): round-2 rule with registry alpha_pub (2026-10-11)
verdict_fn <- function(pub_b, pub_p_sig, b0, p0, rebuilt = FALSE) {
  # pub_p_sig: TRUE if published term significant at .05
  if (is.na(pub_b)) return("NO_PUBLISHED")
  same_sign <- sign(pub_b) == sign(b0); same_sig <- pub_p_sig == (p0 < 0.05)
  if (rebuilt) stop("rebuild verdicts go through rebuild_verdict() (round-2 rule, work_p8/round2_rule.R)")
  if (!(same_sign && same_sig)) return("NOT_REPRODUCED")
  if (abs(b0 - pub_b) / abs(pub_b) <= 0.10) "REPRODUCED" else "CLOSE"
}
write_verify <- function(id, short_cite, key_family, key_term, headline, pub_b, pub_se, pub_n, pub_est,
                         mism, fixed, unfix, pub_sig, R, n_used, notes, rebuilt = FALSE, verdict = NULL, fsuffix = key_term) {
  g <- function(m, k) if (is.null(R[[m]])) NA_real_ else unname(R[[m]][k])
  b0 <- g("m0", "b"); p0 <- g("m0", "p")
  pc <- function(x) 100 * (x - b0) / abs(b0)
  if (is.null(verdict)) verdict <- if (rebuilt) rebuild_verdict(id, key_term, b0, p0, pub_b, pub_sig, pub_se) else verdict_fn(pub_b, pub_sig, b0, p0, rebuilt)
  row <- data.frame(id = id, short_cite = short_cite, key_family = key_family, key_term = key_term,
    headline_table_model = headline, pub_b = pub_b, pub_se = pub_se, pub_n = pub_n, pub_estimator = pub_est,
    spec_mismatches_before = mism, fixed = fixed, unfixable = unfix, verdict = verdict,
    b0 = b0, se0 = g("m0", "se"), p0 = p0, n_used = n_used,
    pct_vs_pub = if (is.na(pub_b)) NA else 100 * (b0 - pub_b) / abs(pub_b),
    b_dm = g("mDM", "b"), se_dm = g("mDM", "se"), p_dm = g("mDM", "p"), pct_dm = pc(g("mDM", "b")),
    b_26 = g("m26", "b"), p_26 = g("m26", "p"), pct_26 = pc(g("m26", "b")),
    pct_bC = pc(g("bC", "b")), pct_bD = pc(g("bD", "b")), notes = notes, stringsAsFactors = FALSE)
  dir.create(file.path(OUT, "results_verify"), recursive = TRUE, showWarnings = FALSE); f <- file.path(OUT, "results_verify", paste0(id, "__", fsuffix, ".csv"))  # port: OUT dir
  write.csv(row, f, row.names = FALSE); print(t(row[, intersect(c("xid", "id", "key_term", "n_used", "b0", "se0", "p0", "b_dm", "se_dm", "p_dm", "b_26", "p_26", "pct_bC", "pct_bD"), names(row)), drop = FALSE])); invisible(row)  # port: id only, no citation printed
}
