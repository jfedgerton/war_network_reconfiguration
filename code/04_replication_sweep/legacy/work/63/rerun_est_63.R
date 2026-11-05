# ESTAUDIT 63 Bussmann 2010 (Table I Model 1: Stata probit ..., cluster(dyadid)). Copy of rerun_verify_63.R (same data, spec, common sample).
# Check of the SE computation: Stata's probit cluster-robust V uses the OBSERVED Hessian as bread (Newton-Raphson ML);
# fixest feglm uses the expected-information (IRLS) bread. Point estimates are identical. Here both are computed for all six steps.
ROOT <- Sys.getenv("SWEEP_DATA"); OUT <- Sys.getenv("SWEEP_OUT"); CODE <- Sys.getenv("SWEEP_CODE")
stopifnot(nzchar(ROOT), nzchar(OUT), nzchar(CODE))
source(file.path(CODE, "hist/add_history.R")); source(file.path(CODE, "ea2/est_helpers.R"))
suppressPackageStartupMessages(library(fixest)); set.seed(123)
W <- file.path(ROOT, "work/63"); WO <- file.path(OUT, "work/63"); dir.create(WO, recursive = TRUE, showWarnings = FALSE); sink(file.path(WO, "log_est_63.txt"), split = TRUE)
d <- as.data.frame(haven::zap_labels(haven::read_dta(file.path(W, "BussmanndyadicFDIconflictdata.dta"))))
names(d) <- sub("^_spline", "spline", names(d)); stopifnot(nrow(d) == 304424)
f0 <- fatal1 ~ smlinflowsgdpim + smldem + lrgdem + smlrgdp96pcG + demdev + smldep + dircont + lncaprat + py + spline1 + spline2 + spline3
T <- "smlinflowsgdpim"
# Stata-equivalent probit cluster SE: V = H^-1 (sum_g s_g s_g') H^-1 * G/(G-1), H = observed Hessian
probit_oim_cl <- function(f, data, cl) {
  m <- feglm(f, data = data, family = binomial("probit"), vcov = as.formula(paste0("~", cl)), ssc = ssc(adj = FALSE, cluster.adj = TRUE), notes = FALSE)
  stopifnot(m$convStatus); used <- obs(m); dd <- data[used, ]
  X <- model.matrix(f, dd); y <- dd[[all.vars(f)[1]]]; b <- coef(m)[colnames(X)]; stopifnot(!anyNA(b))
  eta <- drop(X %*% b); q <- 2 * y - 1
  lam <- q * dnorm(q * eta) / pnorm(q * eta)
  H <- crossprod(X, X * (-lam * (lam + eta)))                  # observed Hessian (negative definite)
  S <- rowsum(X * lam, dd[[cl]]); G <- nrow(S)
  Hi <- solve(-H); V <- Hi %*% crossprod(S) %*% Hi * G / (G - 1)
  se_oim <- sqrt(diag(V))[T]; ct <- coeftable(m)
  c(b = unname(b[T]), se = unname(se_oim), p = unname(2 * pnorm(-abs(b[T] / se_oim))), se_fixest = unname(ct[T, 2]), p_fixest = unname(ct[T, 4]), n = nobs(m))
}
full <- probit_oim_cl(f0, d[complete.cases(d[, c(all.vars(f0), "dyadid")]), ], "dyadid")
cat("FULL SAMPLE (pub b -2.0411, SE .6197, N 207,089):\n"); print(full)
d <- add_history(d, "statea", "stateb", "year")
HV <- c("hist_dyad_w15","hist_mono_w15","hist_dyad_w26","hist_mono_w26","bm_cinc_w15","bm_democ_w15")
ds <- d[d$year >= 1822 & d$year <= 2011, ]
ds <- ds[complete.cases(ds[, c(all.vars(f0), "dyadid", HV)]) & is.finite(ds$bm_cinc_w15) & is.finite(ds$bm_democ_w15), ]
cat("common sample N:", nrow(ds), "\n"); stopifnot(nrow(ds) == 202241)
ADD <- list(m0 = NULL, mD = "hist_dyad_w15", mDM = c("hist_dyad_w15","hist_mono_w15"), m26 = c("hist_dyad_w26","hist_mono_w26"), bC = "bm_cinc_w15", bD = "bm_democ_w15")
R <- lapply(ADD, function(a) { f <- if (is.null(a)) f0 else update(f0, as.formula(paste(". ~ . +", paste(a, collapse = "+")))); probit_oim_cl(f, ds, "dyadid") })
tab <- do.call(rbind, R); print(tab)
ratio <- tab[, "se"] / tab[, "se_fixest"]; cat("SE ratio OIM/fixest by step:", round(ratio, 4), "\n")
saveRDS(list(full = full, tab = tab), file.path(WO, "est_63.rds"))

# ---- sensitivity hook: runs only under run_sensitivity.R (SWEEP_SENS=1), never in the main sweep ----
if (Sys.getenv("SWEEP_SENS") == "1") {
  source(file.path(CODE, "hist/sens_helper.R"))
  sens_run(id = "63", dat = ds, y = "fatal1", tie = "smlinflowsgdpim",
           ctrl = c("smldem", "lrgdem", "smlrgdp96pcG", "demdev", "smldep", "dircont", "lncaprat", "py", "spline1", "spline2", "spline3"), fe = NULL,
           fit = function(dat, add) probit_oim_cl(if (length(add)) update(f0, as.formula(paste(". ~ . +", paste(add, collapse = " + ")))) else f0, dat, "dyadid")[c("b", "se", "p")])
}

if (max(abs(ratio - 1)) > 0.02) {
  write_est("63", "Bussmann 2010", T, "probit, SE clustered by dyad (Stata: observed-Hessian bread), peace years + 3 splines (Table I Model 1)",
    "pdf (Table I 'probit') + code (Bussmann_testsJPR2010.do)", "probit (fixest point estimates) with Stata-equivalent cluster-robust SE (observed-Hessian bread, G/(G-1))", "N",
    sprintf("SE computation: fixest expected-information bread instead of Stata's observed Hessian (SE ratio OIM/fixest %.3f-%.3f across steps); point estimates identical", min(ratio), max(ratio)),
    "Y", lapply(R, function(r) r[c("b", "se", "p")]), nrow(ds),
    sprintf("Full-sample Stata-equivalent SE %.4f (pub .6197; fixest %.4f). Previous (fixest-SE) p0 %.3g, p_dm %.3g, p_26 %.3g.", full["se"], full["se_fixest"], tab["m0", "p_fixest"], tab["mDM", "p_fixest"], tab["m26", "p_fixest"]))
} else cat("OIM and fixest SEs agree within 2% -> MATCH row kept\n")
sink()
