# Writer for results_estaudit/<id>__<term>.csv (ESTAUDIT_PROTOCOL step 5). R = list(m0,mD,mDM,m26,bC,bD), each c(b,se,p).
write_est_ea1 <- function(id, short_cite, key_term, article_estimator, source_of_estimator, est_steps, match_before,
                          mismatch_detail, fixed, R, n_used, notes, b0C = NULL, b0D = NULL) {
  ROOT <- Sys.getenv("SWEEP_DATA"); OUT <- Sys.getenv("SWEEP_OUT"); CODE <- Sys.getenv("SWEEP_CODE")
  stopifnot(nzchar(ROOT), nzchar(OUT), nzchar(CODE))
  stopifnot(all(c("m0","mD","mDM","m26","bC","bD") %in% names(R)))
  es <- if (length(est_steps) == 1) rep(est_steps, 6) else est_steps; stopifnot(length(es) == 6)
  pc <- function(x, b) 100 * (x - b) / abs(b)
  b0 <- unname(R$m0["b"]); rC <- if (is.null(b0C)) b0 else b0C; rD <- if (is.null(b0D)) b0 else b0D
  row <- data.frame(id = id, short_cite = short_cite, key_term = key_term, article_estimator = article_estimator,
    source_of_estimator = source_of_estimator, est_m0 = es[1], est_mD = es[2], est_mDM = es[3], est_m26 = es[4], est_bC = es[5], est_bD = es[6],
    match_before = match_before, mismatch_detail = mismatch_detail, fixed = fixed,
    b0 = b0, se0 = unname(R$m0["se"]), p0 = unname(R$m0["p"]), n_used = n_used,
    b_dm = unname(R$mDM["b"]), se_dm = unname(R$mDM["se"]), p_dm = unname(R$mDM["p"]), pct_dm = pc(unname(R$mDM["b"]), b0),
    b_26 = unname(R$m26["b"]), p_26 = unname(R$m26["p"]), pct_26 = pc(unname(R$m26["b"]), b0),
    pct_bC = pc(unname(R$bC["b"]), rC), pct_bD = pc(unname(R$bD["b"]), rD), notes = notes, stringsAsFactors = FALSE)
  dir.create(file.path(OUT, "results_estaudit"), recursive = TRUE, showWarnings = FALSE); f <- file.path(OUT, "results_estaudit", paste0(id, "__", gsub("[^A-Za-z0-9_]", "_", key_term), ".csv"))  # port: OUT dir
  write.csv(row, f, row.names = FALSE); cat("wrote", f, "\n"); print(t(row[, intersect(c("xid", "id", "key_term", "n_used", "b0", "se0", "p0", "b_dm", "se_dm", "p_dm", "b_26", "p_26", "pct_bC", "pct_bD"), names(row)), drop = FALSE]))  # port: id only, no citation printed
  invisible(row)
}
