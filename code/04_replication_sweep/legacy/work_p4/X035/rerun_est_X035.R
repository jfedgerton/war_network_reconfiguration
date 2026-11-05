# ESTAUDIT X035 Goldsmith, Chalup & Quinlan 2008 JPR -- copy of rerun_verify_X035.R with ONE estimator change.
# Authors' do-file ("jpr final.do", Model 10): xtgee fcwfataldinit ..., family(binomial) link(logit) corr(ar1) force robust.
# Stata `force` = AR(1) on within-panel ORDER (lags from data ordered by year, calendar gaps ignored). rerun_verify_X035.R used
#   geeglm(waves = year) -> AR(1) spaced by calendar year; the model sample has 13,4xx gap steps (most dyad panels have gaps).
# Evidence (ea3/x035_waves.R, full sample N 28,431): order-based waves reproduce ALL published z's to 2 dp (jntdem -1.36, jntaut -1.52,
#   polsim3 -2.18); calendar waves give -1.37/-1.51/-2.20.
# Fix: waves = within-dyad sequence (1..n_i, ordered by year); singleton panels dropped as before. Everything else identical.
suppressPackageStartupMessages(library(geepack))
ROOT <- Sys.getenv("SWEEP_DATA"); OUT <- Sys.getenv("SWEEP_OUT"); CODE <- Sys.getenv("SWEEP_CODE")
stopifnot(nzchar(ROOT), nzchar(OUT), nzchar(CODE))
source(file.path(CODE, "work_p4/b3/common_b3.R"))
source(file.path(CODE, "work_p4/est_common.R"))
WO <- file.path(OUT, "work_p4/X035"); dir.create(WO, recursive = TRUE, showWarnings = FALSE)  # port: outputs
W <- file.path(ROOT, "work_p4/X035"); sink(file.path(WO, "log_est_X035.txt"), split = TRUE)
d <- haven::read_dta(file.path(W, "GCQjpr.dta"), col_select = c("fcwfataldinit","jntdem","jntaut","polsim3","lnlowdep","lnnumigofill","contigland","majmaj","lndistance","lnparity","year","rgroups","bgdyadid","ccode1","ccode2"))
d <- as.data.frame(lapply(d, as.numeric)); d <- d[d$year > 1885 & d$year < 2001 & d$rgroups < 6, ]
rhs <- c("jntdem","jntaut","polsim3","lnlowdep","lnnumigofill","contigland","majmaj","lndistance","lnparity")
f0 <- as.formula(paste("fcwfataldinit ~", paste(rhs, collapse = " + "))); terms <- c("jntdem","jntaut","polsim3")
drop_single <- function(dd) { tb <- table(dd$bgdyadid); k <- dd$bgdyadid %in% as.numeric(names(tb)[tb > 1]); cat("dropping singleton-panel obs:", sum(!k), "\n"); dd[k, ] }
fit_gee <- function(f, data) {
  dd <- data[complete.cases(data[, c(all.vars(f), "bgdyadid", "year")]), ]; dd <- dd[order(dd$bgdyadid, dd$year), ]
  dd$seq <- ave(dd$year, dd$bgdyadid, FUN = seq_along)
  g <- geeglm(f, family = binomial("logit"), data = dd, id = bgdyadid, waves = seq, corstr = "ar1"); ct <- summary(g)$coefficients
  cat("\n----", deparse(f, width.cutoff = 500), "| N =", nrow(dd), "| alpha =", round(unlist(g$geese$alpha), 4), "\n"); print(round(ct, 6))
  out <- list(k = lapply(setNames(terms, terms), function(t) c(b = ct[t, 1], se = ct[t, 2], p = ct[t, 4])), n = nrow(dd))
  rm(g); gc(); out
}
d1 <- drop_single(d[complete.cases(d[, c(all.vars(f0), "bgdyadid")]), ]); stopifnot(nrow(d1) == 28431)
rep <- fit_gee(f0, d1)
d$yh <- d$year + 1   # DV is a lead (as in verify)
ds <- common_sample(d, f0, "bgdyadid", "ccode1", "ccode2", "yh"); ds <- drop_single(ds)
M <- list(m0 = fit_gee(f0, ds), mD = fit_gee(update(f0, . ~ . + hist_dyad_w15), ds),
  mDM = fit_gee(update(f0, . ~ . + hist_dyad_w15 + hist_mono_w15), ds), m26 = fit_gee(update(f0, . ~ . + hist_dyad_w26 + hist_mono_w26), ds),
  bC = fit_gee(update(f0, . ~ . + bm_cinc_w15), ds), bD = fit_gee(update(f0, . ~ . + bm_democ_w15), ds))
stopifnot(length(unique(sapply(M, `[[`, "n"))) == 1)
# ---- sensitivity hook: runs only under run_sensitivity.R (SWEEP_SENS=1), never in the main sweep ----
if (Sys.getenv("SWEEP_SENS") == "1") {
  source(file.path(CODE, "hist/sens_helper.R"))
  sens_run(id = "X035", dat = ds, y = "fcwfataldinit", tie = "polsim3",
           ctrl = setdiff(rhs, "polsim3"), fe = NULL,
           fit = function(dat, add) { r <- fit_gee(if (length(add)) update(f0, as.formula(paste(". ~ . +", paste(add, collapse = " + ")))) else f0, dat)$k$polsim3
             c(b = unname(r[["b"]]), se = unname(r[["se"]]), p = unname(r[["p"]])) })
}
conv <- function(t) lapply(M, function(m) list(b = unname(m$k[[t]]["b"]), se = unname(m$k[[t]]["se"]), p = unname(m$k[[t]]["p"]), n = m$n))
PUBZ <- c(jntdem = "-0.84 (z -1.36)", jntaut = "-0.50 (z -1.52)", polsim3 = "-0.58 (z -2.18)")
for (t in terms) {
  f_ov <- file.path(OUT, "results_verify", paste0("X035__", t, ".csv"))  # port: prior result used in notes only; absent in pipeline
  ov <- if (file.exists(f_ov)) read.csv(f_ov) else list(b0 = NA_real_, pct_dm = NA_real_)
  est_write("X035", "Goldsmith, Chalup & Quinlan 2008", t,
    "GEE logit, AR(1) with Stata force (order-based lags), semi-robust SE by dyad (xtgee ... corr(ar1) force robust)",
    "code (jpr final.do Model 10) + pdf Table II (GEE logit AR1)",
    "geeglm logit AR(1), waves = within-dyad order (Stata force), singletons dropped, robust SE", "N",
    "previous geeglm used waves = calendar year (rho^|t-s| across gaps) instead of Stata force order-based AR(1)", "Y", conv(t),
    sprintf("Fixed in rerun_est_X035.R. Full-sample repro (force) b=%.4f z=%.2f N=%d vs pub %s. Previous (verify, calendar waves) b0=%.4f pct_dm=%.1f.",
            rep$k[[t]]["b"], rep$k[[t]]["b"] / rep$k[[t]]["se"], rep$n, PUBZ[[t]], ov$b0, ov$pct_dm))
}
sink()
