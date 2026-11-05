# Shared helpers for rerun batch (cids 68 61 73 11 69 55 57).
# fit_glm(): glm logit/probit with optional cluster-robust SE (Stata-style: G/(G-1) only,
#            which is what Stata's logit/probit vce(cluster) uses), returns tie b/se/p/n.
# write_row(): assemble + write results/<cid>.csv row per RERUN_PROTOCOL step 5.
ROOT <- Sys.getenv("SWEEP_DATA"); OUT <- Sys.getenv("SWEEP_OUT"); CODE <- Sys.getenv("SWEEP_CODE")
stopifnot(nzchar(ROOT), nzchar(OUT), nzchar(CODE))
suppressPackageStartupMessages({library(sandwich); library(lmtest); library(dplyr)})
source(file.path(CODE, "hist/add_history.R"))
set.seed(123)

# cluster: NULL (ML SE), "robust" (HC0 * n/(n-1), Stata robust), or a column name
fit_glm <- function(f, data, tie, link = "logit", cluster = NULL, label = "") {
  cc <- complete.cases(data[, unique(c(all.vars(f), if (!is.null(cluster) && cluster != "robust") cluster))])
  if (!all(cc)) data <- data[cc, , drop = FALSE]   # avoid a copy when already complete (memory)
  if (isTRUE(getOption("use_fixest", FALSE))) {
    # memory-lean path (large N): same ML point estimates; SEs set to Stata conventions
    m <- fixest::feglm(f, data = data, family = binomial(link = link), glm.iter = 100,
                       vcov = if (is.null(cluster)) "iid" else if (cluster == "robust") "hetero" else as.formula(paste0("~", cluster)),
                       ssc = fixest::ssc(adj = FALSE, cluster.adj = TRUE), mem.clean = TRUE, nthreads = 1)
    stopifnot(m$convStatus, nobs(m) == nrow(data))
    n <- nobs(m); ct <- fixest::coeftable(m)
  } else {
  m <- glm(f, data = data, family = binomial(link = link), control = glm.control(maxit = 100))
  stopifnot(m$converged)
  n <- nobs(m)
  V <- if (is.null(cluster)) vcov(m) else if (cluster == "robust") {
    sandwich(m) * n / (n - 1)
  } else {
    vcovCL(m, cluster = data[[cluster]], type = "HC0", cadjust = TRUE)  # G/(G-1)
  }
  ct <- coeftest(m, vcov. = V)
  }
  stopifnot(tie %in% rownames(ct))
  cat("\n---- [", label, "] ", deparse(f, width.cutoff = 500), " | N =", n, "\n")
  print(round(ct, 5))
  out <- list(b = unname(ct[tie, 1]), se = unname(ct[tie, 2]), p = unname(ct[tie, 4]), n = n)
  rm(m); invisible(gc()); out
}

pct <- function(b, b0) 100 * (b - b0) / abs(b0)

write_row <- function(cid, short_cite, reproduced, b_pub, b_repro, estimator, tie_var, n_orig, n_used,
                      m0 = NULL, mD = NULL, mDM = NULL, m26 = NULL, bC = NULL, bD = NULL,
                      b0C = NULL, b0D = NULL, ahc, notes) {
  g <- function(m, k) if (is.null(m)) NA else m[[k]]
  b0 <- g(m0, "b")
  # benchmarks compared to baseline on the benchmark's own non-missing sample when supplied
  rC <- if (is.null(b0C)) b0 else b0C; rD <- if (is.null(b0D)) b0 else b0D
  row <- data.frame(cid = cid, short_cite = short_cite, reproduced = reproduced,
    b_pub = b_pub, b_repro = b_repro, estimator = estimator, tie_var = tie_var,
    n_orig = n_orig, n_used = n_used,
    b0 = b0, se0 = g(m0, "se"), p0 = g(m0, "p"),
    b_dm = g(mDM, "b"), se_dm = g(mDM, "se"), p_dm = g(mDM, "p"), pct_dm = pct(g(mDM, "b"), b0),
    b_d = g(mD, "b"), pct_d = pct(g(mD, "b"), b0),
    b_26 = g(m26, "b"), se_26 = g(m26, "se"), p_26 = g(m26, "p"), pct_26 = pct(g(m26, "b"), b0),
    pct_bC = pct(g(bC, "b"), rC), pct_bD = pct(g(bD, "b"), rD),
    author_has_history_control = ahc, notes = notes, stringsAsFactors = FALSE)
  dir.create(file.path(OUT, "results"), recursive = TRUE, showWarnings = FALSE)
  write.csv(row, file.path(OUT, "results", paste0(cid, ".csv")), row.names = FALSE)
  print(t(row[, names(row) != "short_cite", drop = FALSE]))
  invisible(row)
}

fail_row <- function(cid, short_cite, status, estimator, tie_var, ahc, notes, b_pub = NA) {
  write_row(cid, short_cite, status, b_pub, NA, estimator, tie_var, NA, NA, ahc = ahc, notes = notes)
}
