# Estimator-audit (phase 6) writer: results_estaudit/<id>__<term>.csv, ESTAUDIT_PROTOCOL step 5 columns.
ROOT <- Sys.getenv("SWEEP_DATA"); OUT <- Sys.getenv("SWEEP_OUT"); CODE <- Sys.getenv("SWEEP_CODE")
stopifnot(nzchar(ROOT), nzchar(OUT), nzchar(CODE))
EA_COLS <- c("id","short_cite","key_term","article_estimator","source_of_estimator","est_m0","est_mD","est_mDM","est_m26","est_bC","est_bD",
  "match_before","mismatch_detail","fixed","b0","se0","p0","n_used","b_dm","se_dm","p_dm","pct_dm","b_26","p_26","pct_26","pct_bC","pct_bD","notes")
# R: list(m0=, mD=, mDM=, m26=, bC=, bD=), each c(b=, se=, p=) or list with b/se/p/n; n_used given separately
write_est <- function(id, short_cite, key_term, art_est, src, est_steps, match_before, mism, fixed, R, n_used, notes, fname = key_term) {
  g <- function(m, k) { x <- R[[m]]; if (is.null(x)) NA_real_ else unname(as.numeric(x[[k]])) }
  b0 <- g("m0", "b"); pc <- function(x) 100 * (x - b0) / abs(b0)
  es <- if (length(est_steps) == 1) rep(est_steps, 6) else est_steps; stopifnot(length(es) == 6)
  row <- data.frame(id = id, short_cite = short_cite, key_term = key_term, article_estimator = art_est, source_of_estimator = src,
    est_m0 = es[1], est_mD = es[2], est_mDM = es[3], est_m26 = es[4], est_bC = es[5], est_bD = es[6],
    match_before = match_before, mismatch_detail = mism, fixed = fixed,
    b0 = b0, se0 = g("m0", "se"), p0 = g("m0", "p"), n_used = n_used,
    b_dm = g("mDM", "b"), se_dm = g("mDM", "se"), p_dm = g("mDM", "p"), pct_dm = pc(g("mDM", "b")),
    b_26 = g("m26", "b"), p_26 = g("m26", "p"), pct_26 = pc(g("m26", "b")),
    pct_bC = pc(g("bC", "b")), pct_bD = pc(g("bD", "b")), notes = notes, stringsAsFactors = FALSE)
  stopifnot(identical(names(row), EA_COLS), is.finite(b0))
  dir.create(file.path(OUT, "results_estaudit"), recursive = TRUE, showWarnings = FALSE); f <- file.path(OUT, "results_estaudit", paste0(id, "__", gsub("[^A-Za-z0-9_.]", "_", fname), ".csv"))  # port: OUT dir
  write.csv(row, f, row.names = FALSE); print(t(row[, intersect(c("xid", "id", "key_term", "n_used", "b0", "se0", "p0", "b_dm", "se_dm", "p_dm", "b_26", "p_26", "pct_bC", "pct_bD"), names(row)), drop = FALSE])); invisible(row)  # port: id only, no citation printed
}
