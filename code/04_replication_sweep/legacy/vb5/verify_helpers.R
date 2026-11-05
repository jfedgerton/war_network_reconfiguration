# Batch-5 verify helpers: write results_verify/<id>__<term>.csv (VERIFY_PROTOCOL step 5)
ROOT <- Sys.getenv("SWEEP_DATA"); OUT <- Sys.getenv("SWEEP_OUT"); CODE <- Sys.getenv("SWEEP_CODE")
stopifnot(nzchar(ROOT), nzchar(OUT), nzchar(CODE))
source(file.path(CODE, "work_p8/round2_rule.R"))   # rebuild_verdict(): round-2 rule with registry alpha_pub (2026-10-11)
pctf <- function(b, b0) 100 * (b - b0) / abs(b0)
write_verify <- function(id, short_cite, key_family, key_term, headline, pub_b, pub_se, pub_n, pub_est,
                         mism, fixed, unfix, verdict, R, notes, term_file = key_term) {
  g <- function(m, k) if (is.null(m)) NA else m[[k]]
  b0 <- g(R$m0, "b")
  row <- data.frame(id = id, short_cite = short_cite, key_family = key_family, key_term = key_term,
    headline_table_model = headline, pub_b = pub_b, pub_se = pub_se, pub_n = pub_n, pub_estimator = pub_est,
    spec_mismatches_before = mism, fixed = fixed, unfixable = unfix, verdict = verdict,
    b0 = b0, se0 = g(R$m0, "se"), p0 = g(R$m0, "p"), n_used = g(R$m0, "n"),
    pct_vs_pub = if (is.na(pub_b)) NA else pctf(b0, pub_b),
    b_dm = g(R$mDM, "b"), se_dm = g(R$mDM, "se"), p_dm = g(R$mDM, "p"), pct_dm = pctf(g(R$mDM, "b"), b0),
    b_26 = g(R$m26, "b"), p_26 = g(R$m26, "p"), pct_26 = pctf(g(R$m26, "b"), b0),
    pct_bC = pctf(g(R$bC, "b"), b0), pct_bD = pctf(g(R$bD, "b"), b0), notes = notes, stringsAsFactors = FALSE)
  dir.create(file.path(OUT, "results_verify"), recursive = TRUE, showWarnings = FALSE); write.csv(row, file.path(OUT, "results_verify", paste0(id, "__", term_file, ".csv")), row.names = FALSE)  # port: OUT dir
  print(t(row[, intersect(c("xid", "id", "key_term", "n_used", "b0", "se0", "p0", "b_dm", "se_dm", "p_dm", "b_26", "p_26", "pct_bC", "pct_bD"), names(row)), drop = FALSE])); invisible(row)  # port: id only, no citation printed
}
# verdict for reproduced (non-rebuilt) studies
verdict_repro <- function(b0, p0, pub_b, pub_sig) {
  if (is.na(pub_b)) return("NO_PUBLISHED")
  same_sig <- (p0 < .05) == pub_sig
  if (sign(b0) != sign(pub_b) || !same_sig) return("NOT_REPRODUCED")
  if (abs(b0 - pub_b) / abs(pub_b) <= 0.10) "REPRODUCED" else "CLOSE"
}
verdict_rebuilt <- function(b0, p0, pub_b, pub_sig, id, term, pub_se = NA) rebuild_verdict(id, term, b0, p0, pub_b, pub_sig, pub_se)
