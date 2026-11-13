# ESTAUDIT X092 Goldsmith 2006 EJIR -- copy of rerun_verify_X092.R with ONE estimator change.
# Author's do-file (models.do, Table 1 Model 2): xtgee ... , family(binomial) link(logit) corr(ar1) force robust (tsset dyadid year).
# Stata `force` fits AR(1) on unequally spaced panels by treating consecutive observations (ordered by year) as lag 1,
#   i.e. AR(1) on the within-panel ORDER, ignoring calendar gaps. rerun_verify_X092.R used geeglm(waves = year), which
#   spaces the AR(1) by calendar year (rho^|t-s|) -> different working correlation whenever a panel has gaps
#   (1,028 gap steps in 612 of 1,154 dyads in the model sample).
# Fix: waves = within-dyad sequence index (1..n_i) after ordering by year = Stata force. Stata xtgee with corr(ar1) also cannot use
#   single-observation panels; these are dropped (reported). Everything else identical. Seed 123.
# [2026-11-08] REVERTED to calendar-spaced AR(1) (waves = year) for the pipeline estimates. Audit: the order-based fit gives
#   smigoabi -0.0201 on the full sample; the calendar-spaced fit gives -0.0139 (smldmat -0.047), matching published -0.014
#   (-0.050). Singletons are still dropped. Both full-sample fits are still printed below for the record.
ROOT <- Sys.getenv("SWEEP_DATA"); OUT <- Sys.getenv("SWEEP_OUT"); CODE <- Sys.getenv("SWEEP_CODE")
stopifnot(nzchar(ROOT), nzchar(OUT), nzchar(CODE))
source(file.path(CODE, "vb4/vhelp.R"))
source(file.path(CODE, "work_p4/est_common.R"))
suppressPackageStartupMessages({library(haven); library(geepack)})
WO <- file.path(OUT, "work_p4/X092"); dir.create(WO, recursive = TRUE, showWarnings = FALSE)  # port: outputs
W <- file.path(ROOT, "work_p4/X092"); sink(file.path(WO, "log_est_X092.txt"), split = TRUE)
K <- c("smldmat", "smldep", "smigoabi")
X <- c("dispute1", K, "lcaprat2", "allies", "noncontg", "logdstab", "minrpwrs", "lsmgdppc", paste0("r", 1:5))
d <- zap_labels(read_dta(file.path(W, "data.dta"), col_select = all_of(c(X, "dyadid", "year", "statea", "stateb")))) %>% as.data.frame()
d <- d[d$year > 1919, ]; stopifnot(!anyDuplicated(d[, c("dyadid", "year")]))
rhs0 <- X[-1]
prep <- function(data, r) {
  data <- data[complete.cases(data[, c("dispute1", r, "dyadid", "year")]), ]
  data <- data[order(data$dyadid, data$year), ]
  tb <- table(data$dyadid); ns <- sum(tb == 1)
  data <- data[data$dyadid %in% as.numeric(names(tb)[tb > 1]), ]
  data$seq <- ave(data$year, data$dyadid, FUN = seq_along)
  attr(data, "n_single") <- ns; data
}
FG <- function(r, data, waves = "seq") {
  data <- prep(data, r)
  f <- as.formula(paste("dispute1 ~", paste(r, collapse = " + ")))
  data$w <- data[[waves]]
  m <- geeglm(f, data = data, id = dyadid, waves = w, family = binomial, corstr = "ar1")
  ct <- summary(m)$coefficients
  cat("\n----", deparse(f, width.cutoff = 500), "| N =", nrow(data), "| singletons dropped =", attr(data, "n_single"),
      "| waves =", waves, "| alpha =", round(unlist(m$geese$alpha), 4), "\n"); print(round(ct[K, ], 5))
  out <- lapply(setNames(K, K), function(t) c(b = ct[t, 1], se = ct[t, 2], p = ct[t, 4])); attr(out, "n") <- nrow(data)
  rm(m); invisible(gc()); out
}
# Full-sample check vs published Table 1 Model 2 (N 35,685): order-based (force) vs calendar-spaced AR(1)
rep_seq <- FG(rhs0, d, "seq"); rep_yr <- FG(rhs0, d, "year")
cat("\nPUBLISHED M2: smldmat -.050 (.010); smldep -32.460 (12.635); smigoabi -.014 (.005); N 35,685\n")
# [2026-11-13 audit] TRIANGLE dispute1 in row year t is MID at t+1 (checked vs COW dyadic MID) -> anchor history on the outcome year
d$yr_out <- d$year + 1
s <- common_sample(d, c("dispute1", rhs0, "dyadid"), "statea", "stateb", "yr_out")
M <- run_six(rhs0, s, function(r, dd) FG(r, dd, "year"))   # [2026-11-08] was "seq"
# ---- sensitivity hook: runs only under run_sensitivity.R (SWEEP_SENS=1), never in the main sweep ----
if (Sys.getenv("SWEEP_SENS") == "1") {
  source(file.path(CODE, "hist/sens_helper.R"))
  # dat = prep(s, rhs0): the rows FG() actually estimates on (FG calls prep(), which drops single-observation dyads from s)
  # [2026-11-10 audit2] hook moved to the new headline term (smldmat)
  sens_run(id = "X092", dat = prep(s, rhs0), y = "dispute1", tie = "smldmat", ctrl = setdiff(rhs0, "smldmat"), fe = NULL,
           fit = function(dat, add) { r <- FG(c(rhs0, add), dat, "year")$smldmat; c(b = unname(r[["b"]]), se = unname(r[["se"]]), p = unname(r[["p"]])) })
}
conv <- function(t) lapply(M, function(m) list(b = unname(m[[t]]["b"]), se = unname(m[[t]]["se"]), p = unname(m[[t]]["p"]), n = attr(m, "n")))
f_old <- file.path(OUT, "results_verify", "X092__smldep.csv")  # port: prior result (unused below); absent in pipeline
old <- if (file.exists(f_old)) read.csv(f_old) else NULL
for (k in K) {
  f_ov <- file.path(OUT, "results_verify", paste0("X092__", k, ".csv"))  # port: prior result used in notes only; absent in pipeline
  ov <- if (file.exists(f_ov)) read.csv(f_ov) else list(b0 = NA_real_, pct_dm = NA_real_)
  est_write("X092", "Goldsmith 2006", k, "GEE logit, AR(1) working corr. spaced by calendar year, semi-robust SE by dyad (author: xtgee ... corr(ar1) force robust)",
    "code (models.do Table 1 Model 2) + pdf p.541-542 (GEE logit, AR1, semi-robust)",
    "geeglm logit AR(1), waves = calendar year, singletons dropped, robust SE",
    "N", "calendar-spaced AR(1) instead of Stata force order-based lags (reproduces the published estimates; the order-based version does not); singleton panels dropped",
    "Y", conv(k),
    sprintf("Reverted 2026-11-08. Full-sample repro (calendar-spaced): b=%.5g N=%d (order-based force: b=%.5g); pub M2 N 35,685. Previous (verify) b0=%.5g pct_dm=%.1f. Common sample, %s.",
            rep_yr[[k]]["b"], attr(rep_yr, "n"), rep_seq[[k]]["b"], ov$b0, ov$pct_dm, "all six steps on the same calendar-spaced AR(1) GEE"))
}
sink()
