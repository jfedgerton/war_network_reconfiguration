# Shared helpers for recovery batch A (copy of rerun_common_b7.R; Stata-style ML cluster SE: G/(G-1) only).
# fit_tie(): fit one model, return b/se/p of the tie variable.
# write_row(): assemble + write results/<cid>.csv row per RERUN_PROTOCOL step 5.
ROOT <- Sys.getenv("SWEEP_DATA"); OUT <- Sys.getenv("SWEEP_OUT"); CODE <- Sys.getenv("SWEEP_CODE")
stopifnot(nzchar(ROOT), nzchar(OUT), nzchar(CODE))
suppressPackageStartupMessages({library(sandwich); library(lmtest); library(dplyr)})
source(file.path(CODE, "hist/add_history.R"))
set.seed(123)

# fixest::feglm (lean on memory). cluster = NULL -> iid ML SEs; else column name for clustered SEs
# (fixest small-sample adj ~ Stata's G/(G-1)). Formula may carry "| fe" fixed effects.
# cache: optional .rds path; if present the stored result is reused (shared box gets OOM-killed; resume-safe).
fit_tie <- function(f, data, family, tie, cluster = NULL, cache = NULL) {
  if (!is.null(cache) && file.exists(cache)) { o <- readRDS(cache); cat("\n---- [cached]", o$txt, sep = "\n"); return(o) }
  vc <- if (is.null(cluster)) "iid" else as.formula(paste0("~", cluster))
  m <- fixest::feglm(f, data = data, family = family, vcov = vc, ssc = fixest::ssc(adj = FALSE, cluster.adj = TRUE), notes = FALSE, mem.clean = TRUE, lean = TRUE)
  stopifnot(m$convStatus)
  ct <- fixest::coeftable(m)
  stopifnot(tie %in% rownames(ct))
  txt <- c(paste("----", deparse(f, width.cutoff = 500), "| N =", nobs(m), "| vcov:", if (is.null(cluster)) "iid" else cluster),
           capture.output(print(round(ct, 5))))
  cat("\n", txt, sep = "\n")
  out <- list(b = ct[tie, 1], se = ct[tie, 2], p = ct[tie, 4], n = nobs(m), txt = txt); rm(m); gc()
  if (!is.null(cache)) saveRDS(out, cache)
  out
}

pct <- function(b, b0) 100 * (b - b0) / abs(b0)

write_row <- function(cid, short_cite, reproduced, b_pub, b_repro, estimator, tie_var, n_orig, n_used,
                      m0 = NULL, mD = NULL, mDM = NULL, m26 = NULL, bC = NULL, bD = NULL,
                      ahc, notes) {
  g <- function(m, k) if (is.null(m)) NA else m[[k]]
  b0 <- g(m0, "b")
  row <- data.frame(cid = cid, short_cite = short_cite, reproduced = reproduced,
    b_pub = b_pub, b_repro = b_repro, estimator = estimator, tie_var = tie_var,
    n_orig = n_orig, n_used = n_used,
    b0 = b0, se0 = g(m0, "se"), p0 = g(m0, "p"),
    b_dm = g(mDM, "b"), se_dm = g(mDM, "se"), p_dm = g(mDM, "p"), pct_dm = pct(g(mDM, "b"), b0),
    b_d = g(mD, "b"), pct_d = pct(g(mD, "b"), b0),
    b_26 = g(m26, "b"), se_26 = g(m26, "se"), p_26 = g(m26, "p"), pct_26 = pct(g(m26, "b"), b0),
    pct_bC = pct(g(bC, "b"), b0), pct_bD = pct(g(bD, "b"), b0),
    author_has_history_control = ahc, notes = notes, stringsAsFactors = FALSE)
  dir.create(file.path(OUT, "results"), recursive = TRUE, showWarnings = FALSE)
  write.csv(row, file.path(OUT, "results", paste0(cid, ".csv")), row.names = FALSE)
  print(t(row[, names(row) != "short_cite", drop = FALSE]))
  invisible(row)
}
