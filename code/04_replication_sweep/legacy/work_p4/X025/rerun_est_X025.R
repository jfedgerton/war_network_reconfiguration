# ESTAUDIT X025 Gelpi & Grieco 2001 JCR -- copy of rerun_verify_X025.R (models A, B) and rerun_X025.R (full model, dfdemexp) with the
# estimator changed from ordinary logit to the paper's rare-events logit.
# PDF (OCR, p.~806-809): population ~1 million leader dyad-years, ~0.06% crises; case-control design: 100% of crisis dyad-years +
#   small random sample of non-crisis dyad-years (sample 2,418; N used 2,345). "Our logit estimates are corrected for this oversampling
#   on the dependent variable and the skewed distribution ... in the manner suggested by King and Zeng" -> relogit with prior
#   correction (tau stated: ~0.06% = 0.0006) + finite-sample bias correction. Table 1 note: Huber-White SE clustered on leader dyad.
# Implementation (direct R port of Zelig 5.1.7 relogit(), file Zelig/R/model-relogit.R; King & Zeng 2001 eq. 11):
#   b_hat = ML logit; pi = fitted; w1 = tau/ybar, w0 = (1-tau)/(1-ybar), w_i = w1*y + w0*(1-y); W = pi(1-pi)w;
#   Q_ii = x_i'(X'WX)^-1 x_i; xi_i = .5 Q_ii ((1+w1)pi_i - w1); bias = (X'WX)^-1 X'W xi; b_tilde = b_hat - bias;
#   prior correction: intercept -= log(((1-tau)/tau) * (ybar/(1-ybar))) (slopes unaffected).
#   VCE: dyad-clustered sandwich of the ML logit (HC0 * G/(G-1)), times (n/(n+k))^2 (King & Zeng 2001 / Stata relogit).
# Seed 123 (no randomness). Uncorrected (ordinary logit) b0 reported alongside.
ROOT <- Sys.getenv("SWEEP_DATA"); OUT <- Sys.getenv("SWEEP_OUT"); CODE <- Sys.getenv("SWEEP_CODE")
stopifnot(nzchar(ROOT), nzchar(OUT), nzchar(CODE))
source(file.path(CODE, "work_p4/xb1/p4c.R"))
source(file.path(CODE, "work_p4/est_common.R"))
suppressPackageStartupMessages({library(readstata13); library(sandwich)})
WO <- file.path(OUT, "work_p4/X025"); dir.create(WO, recursive = TRUE, showWarnings = FALSE)  # port: outputs
W <- file.path(ROOT, "work_p4/X025"); sink(file.path(WO, "log_est_X025.txt"), split = TRUE)
TAU <- 0.0006
relogit_fit <- function(f, data, terms, cluster = "dyadid", tau = TAU) {
  dd <- data[complete.cases(data[, c(all.vars(f), cluster)]), ]
  m <- glm(f, data = dd, family = binomial("logit"), control = glm.control(epsilon = 1e-12, maxit = 100)); stopifnot(m$converged)
  X <- model.matrix(m); y <- m$y; n <- nrow(X); k <- ncol(X); ybar <- mean(y); pihat <- fitted(m)
  if (is.null(tau)) { w1 <- 1; wi <- rep(1, n) } else { w1 <- tau / ybar; wi <- w1 * y + (1 - tau) / (1 - ybar) * (1 - y) }
  Wt <- pihat * (1 - pihat) * wi
  A <- solve(crossprod(X * sqrt(Wt)))
  Q <- rowSums((X %*% A) * X)
  xi <- if (is.null(tau)) 0.5 * Q * (2 * pihat - 1) else 0.5 * Q * ((1 + w1) * pihat - w1)
  bias <- drop(A %*% crossprod(X, Wt * xi))
  b <- coef(m) - bias
  if (!is.null(tau)) b["(Intercept)"] <- b["(Intercept)"] - log(((1 - tau) / tau) * (ybar / (1 - ybar)))
  V <- vcovCL(m, cluster = dd[[cluster]], type = "HC0", cadjust = TRUE) * (n / (n + k))^2
  se <- sqrt(diag(V)); z <- b / se
  ct <- cbind(b_relogit = b, se = se, p = 2 * pnorm(-abs(z)), b_logit = coef(m), se_logit = sqrt(diag(vcovCL(m, cluster = dd[[cluster]], type = "HC0", cadjust = TRUE))))
  cat("\n----", deparse(f, width.cutoff = 500), "| N =", n, "| ybar =", round(ybar, 4), "| tau =", tau, "\n"); print(round(ct, 6))
  list(k = lapply(setNames(terms, terms), function(t) c(b = ct[t, 1], se = ct[t, 2], p = ct[t, 3], b_logit = ct[t, 4], se_logit = ct[t, 5])), n = n)
}
six <- function(f0, ds, terms) list(m0 = relogit_fit(f0, ds, terms), mD = relogit_fit(update(f0, . ~ . + hist_dyad_w15), ds, terms),
  mDM = relogit_fit(update(f0, . ~ . + hist_dyad_w15 + hist_mono_w15), ds, terms), m26 = relogit_fit(update(f0, . ~ . + hist_dyad_w26 + hist_mono_w26), ds, terms),
  bC = relogit_fit(update(f0, . ~ . + bm_cinc_w15), ds, terms), bD = relogit_fit(update(f0, . ~ . + bm_democ_w15), ds, terms))
conv <- function(M, t) lapply(M, function(m) list(b = unname(m$k[[t]]["b"]), se = unname(m$k[[t]]["se"]), p = unname(m$k[[t]]["p"]), n = m$n))
d <- read.dta13(file.path(W, "attracting-trouble.dta"), convert.factors = FALSE); names(d) <- sub("^_", "x_", names(d))
CTL <- c("chalpol4","dymilbal","defmajp","chmajp","border","dist","chatrans","deftrans","dayint","peaceyrs","x_spline1","x_spline2","x_spline3")
fA <- reformulate(c("defpol4", CTL), "conflict"); fB <- reformulate(c("lndefexp","defpol4","lnchaexp", CTL), "conflict")
fF <- reformulate(c("lndefexp","defpol4","dfdemexp","lnchaexp", CTL), "conflict")
AE <- "rare-events logit (King-Zeng relogit: prior correction for case-control sample, tau ~0.0006, + finite-sample bias correction), Huber-White SE clustered on leader dyad"
SRC <- "pdf (OCR) research-design section + Table 1 note (p.~806-809)"
SE_ <- "relogit (R port of Zelig relogit: bias corr. + prior corr. tau=0.0006), dyad-clustered SE x (n/(n+k))^2"
MIS <- "previous rerun used ordinary logit (no King-Zeng bias correction, no prior correction) in all six steps"
PUB <- list(defpol4_notenure = "0.027 (0.011)", defpol4 = "0.016 (0.012)", lndefexp = "-0.274 (0.078)", dfdemexp = "not published (fn 16)")
note_row <- function(M, t, rep, old_file, tt) {
  ov <- if (file.exists(old_file)) read.csv(old_file) else list(pct_dm = NA_real_)  # port: prior result used in notes only; absent in pipeline
  sprintf("Fixed in rerun_est_X025.R. b0 relogit=%.5g vs uncorrected logit b0=%.5g (same sample). Full-sample repro relogit b=%.5g (logit %.5g), N=%d; pub %s. Previous (logit) pct_dm=%.1f. Slopes differ from logit only via the bias correction (prior correction moves the intercept only).",
          M$m0$k[[tt]]["b"], M$m0$k[[tt]]["b_logit"], rep$k[[tt]]["b"], rep$k[[tt]]["b_logit"], rep$n, PUB[[t]], ov$pct_dm)
}
# A: Democracy Model (Table 1 col 1)
rA <- relogit_fit(fA, d, "defpol4"); sA <- common_sample(d, fA, "dyadid", "chall", "def", "year"); MA <- six(fA, sA, "defpol4")
stopifnot(length(unique(sapply(MA, `[[`, "n"))) == 1)
est_write("X025", "Gelpi & Grieco 2001 JCR", "defpol4_notenure", AE, SRC, SE_, "N", MIS, "Y", conv(MA, "defpol4"),
  note_row(MA, "defpol4_notenure", rA, file.path(OUT, "results_verify/X025__defpol4_notenure.csv"), "defpol4"))
# B: Leadership Tenure Model (Table 1 col 2)
rB <- relogit_fit(fB, d, c("defpol4","lndefexp")); sB <- common_sample(d, fB, "dyadid", "chall", "def", "year"); MB <- six(fB, sB, c("defpol4","lndefexp"))
stopifnot(length(unique(sapply(MB, `[[`, "n"))) == 1)
# ---- sensitivity hook: runs only under run_sensitivity.R (SWEEP_SENS=1), never in the main sweep ----
# (headline row = lndefexp in model B, the Leadership Tenure Model)
if (Sys.getenv("SWEEP_SENS") == "1") {
  source(file.path(CODE, "hist/sens_helper.R"))
  sens_run(id = "X025", dat = sB, y = "conflict", tie = "lndefexp", ctrl = c("defpol4", "lnchaexp", CTL), fe = NULL,
           fit = function(dat, add) { r <- relogit_fit(if (length(add)) update(fB, as.formula(paste(". ~ . +", paste(add, collapse = " + ")))) else fB, dat, "lndefexp")$k[["lndefexp"]]
             c(b = unname(r["b"]), se = unname(r["se"]), p = unname(r["p"])) })
}
est_write("X025", "Gelpi & Grieco 2001 JCR", "defpol4", AE, SRC, SE_, "N", MIS, "Y", conv(MB, "defpol4"),
  note_row(MB, "defpol4", rB, file.path(OUT, "results_verify/X025__defpol4.csv"), "defpol4"))
est_write("X025", "Gelpi & Grieco 2001 JCR", "lndefexp", AE, SRC, SE_, "N", MIS, "Y", conv(MB, "lndefexp"),
  note_row(MB, "lndefexp", rB, file.path(OUT, "results_verify/X025__lndefexp.csv"), "lndefexp"))
# Full model with tenure x regime (dfdemexp; fn 16, not published)
rF <- relogit_fit(fF, d, "dfdemexp"); sF <- common_sample(d, fF, "dyadid", "chall", "def", "year"); MF <- six(fF, sF, "dfdemexp")
stopifnot(length(unique(sapply(MF, `[[`, "n"))) == 1)
est_write("X025", "Gelpi & Grieco 2001 JCR", "dfdemexp", AE, SRC, SE_, "N", MIS, "Y", conv(MF, "dfdemexp"),
  note_row(MF, "dfdemexp", rF, file.path(OUT, "results_p4/X025__dfdemexp.csv"), "dfdemexp"))
# Sensitivity: bias correction without prior correction (tau = NULL) for m0 / mDM of model B
cat("\nSENSITIVITY tau=NULL (bias correction only), model B:\n")
s0 <- relogit_fit(fB, sB, c("defpol4","lndefexp"), tau = NULL); s1 <- relogit_fit(update(fB, . ~ . + hist_dyad_w15 + hist_mono_w15), sB, c("defpol4","lndefexp"), tau = NULL)
sink()
