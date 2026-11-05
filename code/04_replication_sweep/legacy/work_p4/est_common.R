# Shared writer for ESTAUDIT (phase 6) rows: results_estaudit/<id>__<term>.csv
ROOT <- Sys.getenv("SWEEP_DATA"); OUT <- Sys.getenv("SWEEP_OUT"); CODE <- Sys.getenv("SWEEP_CODE")
stopifnot(nzchar(ROOT), nzchar(OUT), nzchar(CODE))
est_write <- function(id, cite, term, art_est, src, steps_est, match_before, mism, fixed, M, notes) {
  # M: list m0,mD,mDM,m26,bC,bD each list(b,se,p,n)
  stopifnot(all(c("m0","mD","mDM","m26","bC","bD") %in% names(M)))
  g <- function(m, x) unname(M[[m]][[x]])
  b0 <- g("m0", "b"); pc <- function(m) 100 * (g(m, "b") - b0) / abs(b0)
  if (length(steps_est) == 1) steps_est <- rep(steps_est, 6)
  row <- data.frame(id = id, short_cite = cite, key_term = term, article_estimator = art_est, source_of_estimator = src,
    est_m0 = steps_est[1], est_mD = steps_est[2], est_mDM = steps_est[3], est_m26 = steps_est[4], est_bC = steps_est[5], est_bD = steps_est[6],
    match_before = match_before, mismatch_detail = mism, fixed = fixed,
    b0 = b0, se0 = g("m0", "se"), p0 = g("m0", "p"), n_used = g("m0", "n"),
    b_dm = g("mDM", "b"), se_dm = g("mDM", "se"), p_dm = g("mDM", "p"), pct_dm = pc("mDM"),
    b_26 = g("m26", "b"), p_26 = g("m26", "p"), pct_26 = pc("m26"), pct_bC = pc("bC"), pct_bD = pc("bD"),
    notes = notes, stringsAsFactors = FALSE)
  dir.create(file.path(OUT, "results_estaudit"), recursive = TRUE, showWarnings = FALSE)
  write.csv(row, file.path(OUT, "results_estaudit", paste0(id, "__", term, ".csv")), row.names = FALSE)
  print(t(row[, intersect(c("xid", "id", "key_term", "n_used", "b0", "se0", "p0", "b_dm", "se_dm", "p_dm", "b_26", "p_26", "pct_bC", "pct_bD"), names(row)), drop = FALSE])); invisible(row)  # port: id only, no citation printed
}
