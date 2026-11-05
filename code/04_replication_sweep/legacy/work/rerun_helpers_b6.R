# Shared helpers for cids 88 17 66 29 50 1 (sourced by each rerun_<cid>.R)
ROOT <- Sys.getenv("SWEEP_DATA"); OUT <- Sys.getenv("SWEEP_OUT"); CODE <- Sys.getenv("SWEEP_CODE")
stopifnot(nzchar(ROOT), nzchar(OUT), nzchar(CODE))
source(file.path(CODE, "hist/add_history.R"))
set.seed(123)

# extract coef/se/p for a term from a coeftest-like matrix
getc <- function(ct, term) {
  stopifnot(term %in% rownames(ct))
  c(b = unname(ct[term, 1]), se = unname(ct[term, 2]), p = unname(ct[term, 4]))
}
pct <- function(bx, b0) 100 * (bx - b0) / abs(b0)

RESULT_COLS <- c("cid","short_cite","reproduced","b_pub","b_repro","estimator","tie_var","n_orig","n_used",
  "b0","se0","p0","b_dm","se_dm","p_dm","pct_dm","b_d","pct_d","b_26","se_26","p_26","pct_26",
  "pct_bC","pct_bD","author_has_history_control","notes")

write_result <- function(row) {
  row <- row[RESULT_COLS]
  row[sapply(row, is.null)] <- NA
  # round non-integer numerics to 6 significant digits; keep counts (integer-valued) exact
  df <- as.data.frame(lapply(row, function(x) if (is.numeric(x) && !isTRUE(x == round(x))) signif(x, 6) else x),
                      stringsAsFactors = FALSE)
  names(df) <- RESULT_COLS
  dir.create(file.path(OUT, "results"), recursive = TRUE, showWarnings = FALSE)
  write.csv(df, file.path(OUT, "results", paste0(row$cid, ".csv")), row.names = FALSE)
  print(t(df[, names(df) != "short_cite", drop = FALSE]))
}

# fill the model-sequence numbers from a list of coef vectors (m0,mD,mDM,m26,bC,bD)
fill_seq <- function(row, r) {
  row$b0 <- r$m0["b"]; row$se0 <- r$m0["se"]; row$p0 <- r$m0["p"]
  row$b_d <- r$mD["b"]; row$pct_d <- pct(r$mD["b"], r$m0["b"])
  row$b_dm <- r$mDM["b"]; row$se_dm <- r$mDM["se"]; row$p_dm <- r$mDM["p"]; row$pct_dm <- pct(r$mDM["b"], r$m0["b"])
  row$b_26 <- r$m26["b"]; row$se_26 <- r$m26["se"]; row$p_26 <- r$m26["p"]; row$pct_26 <- pct(r$m26["b"], r$m0["b"])
  row$pct_bC <- pct(r$bC["b"], r$m0["b"]); row$pct_bD <- pct(r$bD["b"], r$m0["b"])
  lapply(row, function(x) if (is.numeric(x)) unname(x) else x)
}

ADD_SETS <- list(m0 = character(0), mD = "hist_dyad_w15", mDM = c("hist_dyad_w15","hist_mono_w15"),
                 m26 = c("hist_dyad_w26","hist_mono_w26"), bC = "bm_cinc_w15", bD = "bm_democ_w15")

# Memory-constrained box (shared cgroup, OOM kills): cache each expensive step so a killed run can resume.
cached <- function(W, key, expr) {
  f <- file.path(W, paste0("cache_", key, ".rds"))
  if (file.exists(f)) return(readRDS(f))
  val <- expr; saveRDS(val, f); invisible(gc()); val
}
