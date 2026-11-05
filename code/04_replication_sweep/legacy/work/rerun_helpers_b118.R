# Shared helpers for reruns of cids 118 65 52 119 103 4 (steps 3-5 of RERUN_PROTOCOL.md)
suppressPackageStartupMessages({library(dplyr); library(fixest)})
ROOT <- Sys.getenv("SWEEP_DATA"); OUT <- Sys.getenv("SWEEP_OUT"); CODE <- Sys.getenv("SWEEP_CODE")
stopifnot(nzchar(ROOT), nzchar(OUT), nzchar(CODE))
source(file.path(CODE, "hist/add_history.R"))
set.seed(123)

# Binary GLM with Stata-style clustered SEs (G/(G-1) only; Stata ML models do not apply (n-1)/(n-k))
fit_glm <- function(rhs, dv, data, link = "logit", cluster = NULL) {
  f <- as.formula(paste(dv, "~", paste(rhs, collapse = " + ")))
  vc <- if (is.null(cluster)) "iid" else as.formula(paste0("~", cluster))
  m <- feglm(f, data = data, family = binomial(link = link), vcov = vc,
             ssc = ssc(adj = FALSE, cluster.adj = TRUE), lean = TRUE, notes = FALSE)
  m
}

tie_stats <- function(m, tie) {
  ct <- coeftable(m)
  stopifnot(tie %in% rownames(ct))
  c(b = unname(ct[tie, 1]), se = unname(ct[tie, 2]), p = unname(ct[tie, 4]))
}

# Run m0, mD, mDM, m26, bC, bD on a common non-missing sample.
# fitter(rhs, data) must return a model; stat(m) returns c(b,se,p) for the tie.
run_battery <- function(d, base_rhs, fitter, stat, logf) {
  add <- list(m0 = character(0), mD = "hist_dyad_w15", mDM = c("hist_dyad_w15", "hist_mono_w15"),
              m26 = c("hist_dyad_w26", "hist_mono_w26"), bC = "bm_cinc_w15", bD = "bm_democ_w15")
  out <- list()
  for (nm in names(add)) {
    m <- fitter(c(base_rhs, add[[nm]]), d)
    out[[nm]] <- stat(m)
    cat("\n==========", nm, "==========\n", file = logf, append = TRUE)
    capture.output(print(summary_model(m)), file = logf, append = TRUE)
    cat(sprintf("[%s] b=%.6g se=%.6g p=%.4g\n", nm, out[[nm]]["b"], out[[nm]]["se"], out[[nm]]["p"]))
    rm(m); invisible(gc())
  }
  out
}
summary_model <- function(m) if (inherits(m, "fixest")) coeftable(m) else if (!is.null(m$coeftable)) m$coeftable else summary(m)

write_result <- function(cid, short_cite, reproduced, b_pub, b_repro, estimator, tie_var,
                         n_orig, n_used, R, ahc, notes) {
  b0 <- R$m0["b"]
  pct <- function(x) 100 * (x - b0) / abs(b0)
  row <- data.frame(cid = cid, short_cite = short_cite, reproduced = reproduced,
    b_pub = b_pub, b_repro = b_repro, estimator = estimator, tie_var = tie_var,
    n_orig = n_orig, n_used = n_used,
    b0 = b0, se0 = R$m0["se"], p0 = R$m0["p"],
    b_dm = R$mDM["b"], se_dm = R$mDM["se"], p_dm = R$mDM["p"], pct_dm = pct(R$mDM["b"]),
    b_d = R$mD["b"], pct_d = pct(R$mD["b"]),
    b_26 = R$m26["b"], se_26 = R$m26["se"], p_26 = R$m26["p"], pct_26 = pct(R$m26["b"]),
    pct_bC = pct(R$bC["b"]), pct_bD = pct(R$bD["b"]),
    author_has_history_control = ahc, notes = notes, row.names = NULL)
  dir.create(file.path(OUT, "results"), recursive = TRUE, showWarnings = FALSE)
  write.csv(row, file.path(OUT, "results", paste0(cid, ".csv")), row.names = FALSE)
  print(row[, names(row) != "short_cite", drop = FALSE])
  invisible(row)
}

write_fail <- function(cid, short_cite, status, estimator, tie_var, ahc, notes, b_pub = NA, b_repro = NA, n_orig = NA) {
  cols <- c("cid","short_cite","reproduced","b_pub","b_repro","estimator","tie_var","n_orig","n_used",
            "b0","se0","p0","b_dm","se_dm","p_dm","pct_dm","b_d","pct_d","b_26","se_26","p_26","pct_26",
            "pct_bC","pct_bD","author_has_history_control","notes")
  row <- as.data.frame(setNames(as.list(rep(NA, length(cols))), cols))
  row$cid <- cid; row$short_cite <- short_cite; row$reproduced <- status; row$estimator <- estimator
  row$tie_var <- tie_var; row$author_has_history_control <- ahc; row$notes <- notes
  row$b_pub <- b_pub; row$b_repro <- b_repro; row$n_orig <- n_orig
  dir.create(file.path(OUT, "results"), recursive = TRUE, showWarnings = FALSE)
  write.csv(row, file.path(OUT, "results", paste0(cid, ".csv")), row.names = FALSE)
  print(row[, names(row) != "short_cite", drop = FALSE])
}

# Restrict to GML w15 range and non-missing on all vars used in any battery model
battery_sample <- function(d, vars) {
  n_in <- nrow(d)
  d <- d[d$year_h >= 1822 & d$year_h <= 2011, ]
  n_yr <- nrow(d)
  hv <- c("hist_dyad_w15","hist_mono_w15","hist_dyad_w26","hist_mono_w26","bm_cinc_w15","bm_democ_w15")
  cc <- complete.cases(as.data.frame(d)[, c(vars, hv)]) & is.finite(d$bm_cinc_w15) & is.finite(d$bm_democ_w15)
  d <- d[cc, unique(c(vars, hv, "year_h"))]
  invisible(gc())
  cat(sprintf("battery sample: %d in -> %d after year restriction -> %d after non-missing\n", n_in, n_yr, nrow(d)))
  d
}
