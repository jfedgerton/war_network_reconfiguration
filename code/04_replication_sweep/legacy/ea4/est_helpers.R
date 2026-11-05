# Estimator-audit (phase 6) output helper, batch 4. Writes results_estaudit/<id>__<term>.csv.
# R: list(m0, mD, mDM, m26, bC, bD), each c(b=, se=, p=) for the key term (mD may be NULL for monadic).
ROOT <- Sys.getenv("SWEEP_DATA"); OUT <- Sys.getenv("SWEEP_OUT"); CODE <- Sys.getenv("SWEEP_CODE")
stopifnot(nzchar(ROOT), nzchar(OUT), nzchar(CODE))
write_est <- function(id, cite, term, art_est, src, steps, match_before, mism, fixed, R, n_used, notes) {
  stopifnot(length(steps) == 6, all(c("m0", "mDM", "m26", "bC", "bD") %in% names(R)))
  g <- function(m, x) if (is.null(R[[m]])) NA else unname(R[[m]][[x]])
  b0 <- g("m0", "b"); pct <- function(b) 100 * (b - b0) / abs(b0)
  stopifnot(is.finite(b0))
  row <- data.frame(id = id, short_cite = cite, key_term = term, article_estimator = art_est, source_of_estimator = src,
    est_m0 = steps[1], est_mD = steps[2], est_mDM = steps[3], est_m26 = steps[4], est_bC = steps[5], est_bD = steps[6],
    match_before = match_before, mismatch_detail = mism, fixed = fixed,
    b0 = b0, se0 = g("m0", "se"), p0 = g("m0", "p"), n_used = n_used,
    b_dm = g("mDM", "b"), se_dm = g("mDM", "se"), p_dm = g("mDM", "p"), pct_dm = pct(g("mDM", "b")),
    b_26 = g("m26", "b"), p_26 = g("m26", "p"), pct_26 = pct(g("m26", "b")),
    pct_bC = pct(g("bC", "b")), pct_bD = pct(g("bD", "b")), notes = notes, stringsAsFactors = FALSE, row.names = NULL)
  dir.create(file.path(OUT, "results_estaudit"), recursive = TRUE, showWarnings = FALSE)
  write.csv(row, file.path(OUT, "results_estaudit", paste0(id, "__", term, ".csv")), row.names = FALSE)
  print(t(row[, intersect(c("xid", "id", "key_term", "n_used", "b0", "se0", "p0", "b_dm", "se_dm", "p_dm", "b_26", "p_26", "pct_bC", "pct_bD"), names(row)), drop = FALSE])); invisible(row)  # port: id only, no citation printed
}
