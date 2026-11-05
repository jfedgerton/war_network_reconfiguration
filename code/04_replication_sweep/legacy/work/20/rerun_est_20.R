# ESTAUDIT 20 Oneal & Russett 2001 IO -- copy of rerun_verify_20.R with ONE estimator change.
# Paper (p.472-3, fn 15, fn 35): GEE logit, AR(1), robust SE by dyad, estimated in Stata (xtgee). Stata xtgee corr(ar1) refuses
#   unequally spaced panels unless `force` is given, in which case lags are taken from the data ORDER (calendar gaps ignored).
#   The model sample has 1,186 gap steps in 553 of 884 dyads, so a Stata fit implies order-based AR(1). (Authors' code not available;
#   X035 test showed order-based geeglm waves reproduce a Stata `force` table to 2 dp.)
# rerun_verify_20.R used geeglm(waves = year) (calendar-spaced AR(1)). Fix: waves = within-dyad order; singleton panels dropped
#   (Stata cannot use them for AR(1)). Everything else identical (TRIANGLE PR dyads, >=20-yr dyads, history on outcome year).
ROOT <- Sys.getenv("SWEEP_DATA"); OUT <- Sys.getenv("SWEEP_OUT"); CODE <- Sys.getenv("SWEEP_CODE")
stopifnot(nzchar(ROOT), nzchar(OUT), nzchar(CODE))
source(file.path(CODE, "hist/add_history.R")); source(file.path(CODE, "work_p4/est_common.R"))
suppressPackageStartupMessages(library(geepack)); set.seed(123)
WO <- file.path(OUT, "work/20"); dir.create(WO, recursive = TRUE, showWarnings = FALSE)
sink(file.path(WO, "log_est_20.txt"), split = TRUE)
d <- read.delim(file.path(ROOT, "recA/or/TRIANGLE.tab"))
d$smldmat <- pmin(d$demauta, d$demautb); d$smldep <- pmin(d$dependa, d$dependb); d$dyadid <- 1000 * d$statea + d$stateb
d <- d[order(d$dyadid, d$year), ]; d$nyrs <- ave(d$year, d$dyadid, FUN = length); d <- d[d$nyrs >= 20, ]
d$dvyear <- d$year + 1
X <- "noncontg + lcaprat2 + allies + smldmat + smldep + logdstab + minrpwrs"; f0 <- as.formula(paste("dispute1 ~", X))
fit_gee <- function(f, data, tie = "smldep", waves = "seq") {
  data <- data[complete.cases(data[, all.vars(f)]), ]; data <- data[order(data$dyadid, data$year), ]
  tb <- table(data$dyadid); ns <- sum(tb == 1); data <- data[data$dyadid %in% as.numeric(names(tb)[tb > 1]), ]
  data$seq <- ave(data$year, data$dyadid, FUN = seq_along); data$w <- data[[waves]]
  m <- geeglm(f, data = data, id = dyadid, waves = w, family = binomial, corstr = "ar1")
  ct <- summary(m)$coefficients; cat("\n----", deparse(f, width.cutoff = 500), "| N =", nrow(data), "| singletons dropped", ns, "| waves", waves, "| alpha", round(unlist(m$geese$alpha), 4), "\n"); print(round(ct, 5))
  list(b = ct[tie, 1], se = ct[tie, 2], p = ct[tie, 4], n = nrow(data))
}
full <- fit_gee(f0, d)   # full-sample (all years) order-based fit, for reference vs published -50.4
d <- add_history(d, "statea", "stateb", "dvyear")
s <- subset(d, dvyear >= 1822 & dvyear <= 2011)
s <- s[complete.cases(s[, c(all.vars(f0), "hist_dyad_w15","hist_mono_w15","hist_dyad_w26","hist_mono_w26","bm_cinc_w15","bm_democ_w15")]), ]
R <- list(m0 = fit_gee(f0, s), mD = fit_gee(update(f0, . ~ . + hist_dyad_w15), s),
  mDM = fit_gee(update(f0, . ~ . + hist_dyad_w15 + hist_mono_w15), s), m26 = fit_gee(update(f0, . ~ . + hist_dyad_w26 + hist_mono_w26), s),
  bC = fit_gee(update(f0, . ~ . + bm_cinc_w15), s), bD = fit_gee(update(f0, . ~ . + bm_democ_w15), s))
stopifnot(length(unique(sapply(R, `[[`, "n"))) == 1)
# ---- sensitivity hook: runs only under run_sensitivity.R (SWEEP_SENS=1), never in the main sweep ----
# (fit_gee() drops single-observation dyads before fitting; s is already complete on all model columns, so the model sample is s without those dyads)
if (Sys.getenv("SWEEP_SENS") == "1") {
  source(file.path(CODE, "hist/sens_helper.R"))
  sens_run(id = "20", dat = s[ave(s$year, s$dyadid, FUN = length) > 1, ], y = "dispute1", tie = "smldep",
           ctrl = c("noncontg", "lcaprat2", "allies", "smldmat", "logdstab", "minrpwrs"), fe = NULL,
           fit = function(dat, add) { r <- fit_gee(if (length(add)) update(f0, as.formula(paste(". ~ . +", paste(add, collapse = " + ")))) else f0, dat)
             c(b = unname(r$b), se = unname(r$se), p = unname(r$p)) })
}
# port: prior sweep result (results_verify/20__smldep.csv, from rerun_verify_20.R) is only quoted in the notes string; not produced
# by this pipeline step, so fall back to NA when absent. No effect on estimates.
ovf <- file.path(OUT, "results_verify", "20__smldep.csv")
ov <- if (file.exists(ovf)) read.csv(ovf) else data.frame(b0 = NA_real_, pct_dm = NA_real_)
est_write("20", "Oneal & Russett 2001", "smldep",
  "GEE logit, AR(1), robust SE by dyad (Stata xtgee; order-based AR(1) implied by force on gapped panels)",
  "pdf p.472-473 (GEE logit AR1, robust SE clustered by dyad; Stata per fn 15/35)",
  "geeglm logit AR(1), waves = within-dyad order (Stata force), singletons dropped, robust SE", "N",
  "previous geeglm used waves = calendar year (calendar-spaced AR(1)); Stata xtgee needs force on gapped panels -> order-based lags",
  "Y", R,
  sprintf("Fixed in rerun_est_20.R (REBUILT TRIANGLE PR-dyad sample as in verify; all-dyad sample still unavailable). Full-sample order-based b=%.4g (pub -50.4, all dyads). Previous (verify, calendar waves) b0=%.4g pct_dm=%.1f. Authors' code not available: force/order-based AR(1) is inferred from Stata behavior, not documented.",
          full$b, ov$b0, ov$pct_dm))
sink()
