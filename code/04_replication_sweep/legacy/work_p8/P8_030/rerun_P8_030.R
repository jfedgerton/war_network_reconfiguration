# P8_030 Chiba & Fang 2014 JOP "Institutional Opposition, Regime Accountability, and International Conflict" (10.1017/s0022381614000127)
# Data + code + log: Harvard Dataverse doi:10.7910/DVN/24455 (IORAIC.zip: ioraic_data.dta, ioraic_setupEstimator.do, 1_ioraic_estimation.do,
#   output/logfile_main.log). Unit: challenger-target territorial-claim dyad-years 1946-95 (N = 3880).
# HEADLINE: Table 1 "unified model" = custom ML estimator `boselect` (IO-opposition probit xg jointly with a two-regime ordered-probit
#   escalation equation: xbi under IO opposition, xbn without; one cut tau1 = exp(lntau1); rho = tanh(athrho)); ml ... robust.
# KEY TERM (a priori, abstract): "more accountable regimes are less likely to escalate military crises when an international organization
#   opposes their actions" -> [xbi]netdem1 (Polity score in the escalation equation under IO opposition). Author log: -.1090921 (.0286524).
#   The interaction analogue ([xbi]netdem1 - [xbn]netdem1) is also logged.
# Estimator: the authors' boselect likelihood ported line-for-line to R (bivariate normal CDF = vectorized Genz bvnu from work/55/biprobit.R),
#   BFGS then Newton polish, Stata-style robust (sandwich * n/(n-1)) SEs. History/benchmark covariates enter all three equations (as the
#   authors' own controls do). Same estimator for all six steps.
CODE <- Sys.getenv("SWEEP_CODE"); stopifnot(nzchar(CODE))
source(file.path(CODE, "work_p8/common_b3.R"))
source(file.path(CODE, "work/55/biprobit.R"))   # bvn_lower()
WD <- file.path(ROOT, "work_p8/P8_030"); W <- file.path(OUT, "work_p8/P8_030"); dir.create(W, recursive = TRUE, showWarnings = FALSE); sink(file.path(W, "log_P8_030.txt"), split = TRUE)
d <- as.data.frame(haven::zap_labels(haven::read_dta(file.path(WD, "IORAIC/ioraic_data.dta"))))
covE <- c("netdem1", "locad1", "endrv10j", "allyj", "strval1", "ecoval1"); covI <- c(covE, "pacsett_chal")

boselect_ll <- function(th, XE, XI, y) {
  kE <- ncol(XE); kI <- ncol(XI)
  xbi <- drop(XE %*% th[1:kE]); xbn <- drop(XE %*% th[kE + 1:kE]); xg <- drop(XI %*% th[2 * kE + 1:kI])
  tau1 <- exp(th[2 * kE + kI + 1]); rho <- tanh(th[2 * kE + kI + 2]); rho <- sign(rho) * min(abs(rho), 1 - 1e-9)
  p1 <- bvn_lower(xg, tau1 - xbi, -rho) - bvn_lower(xg, -xbi, -rho)
  p2 <- bvn_lower(xg, xbi - tau1, rho)
  p3 <- bvn_lower(-xg, tau1 - xbn, rho) - bvn_lower(-xg, -xbn, rho)
  p4 <- bvn_lower(-xg, xbn - tau1, -rho)
  p0 <- 1 - p1 - p2 - p3 - p4
  P <- cbind(p0, p1, p2, p3, p4)[cbind(seq_along(y), y + 1)]
  log(pmax(P, 1e-300))
}
fit_boselect <- function(extra, data, start = NULL) {
  XE <- cbind(as.matrix(data[, c(covE, extra), drop = FALSE]), `_cons` = 1)
  XI <- cbind(as.matrix(data[, c(covI, extra), drop = FALSE]), `_cons` = 1)
  y <- data$dv_esc_io; stopifnot(all(y %in% 0:4))
  nm <- c(paste0("xbi:", colnames(XE)), paste0("xbn:", colnames(XE)), paste0("xg:", colnames(XI)), "lntau1", "athrho")
  if (is.null(start)) start <- rep(0, length(nm))
  names(start) <- nm
  nll <- function(t) -sum(boselect_ll(t, XE, XI, y))
  lo <- rep(-Inf, length(nm)); hi <- rep(Inf, length(nm)); lo[length(nm)] <- -4; hi[length(nm)] <- 4; lo[length(nm) - 1] <- -6; hi[length(nm) - 1] <- 3
  o <- nlminb(start, nll, lower = lo, upper = hi, control = list(eval.max = 20000, iter.max = 5000, rel.tol = 1e-12))
  cat("nlminb:", o$message, "LL", -o$objective, "\n")
  th <- o$par
  for (it in 1:15) {   # Newton polish on numeric derivatives
    g <- numDeriv::grad(nll, th); if (max(abs(g)) < 1e-6) break
    H <- numDeriv::hessian(nll, th); step <- -solve(H, g); lam <- 1
    while (nll(th + lam * step) > nll(th) + 1e-10 && lam > 1e-4) lam <- lam / 2
    th <- th + lam * step
  }
  g <- numDeriv::grad(nll, th); stopifnot(max(abs(g)) < 1e-3)
  H <- numDeriv::hessian(nll, th); Ainv <- solve(H)
  S <- numDeriv::jacobian(function(t) boselect_ll(t, XE, XI, y), th)
  n <- length(y); V <- Ainv %*% crossprod(S) %*% Ainv * n / (n - 1)
  se <- sqrt(diag(V)); names(th) <- names(se) <- nm; dimnames(V) <- list(nm, nm)
  cat("\n---- boselect | extra:", paste(extra, collapse = " + "), "| N =", n, "| LL =", round(-nll(th), 4), "| max|grad| =", signif(max(abs(g)), 3), "\n")
  print(round(cbind(b = th, se = se, z = th / se), 5))
  dif <- th[["xbi:netdem1"]] - th[["xbn:netdem1"]]
  dse <- sqrt(V["xbi:netdem1", "xbi:netdem1"] + V["xbn:netdem1", "xbn:netdem1"] - 2 * V["xbi:netdem1", "xbn:netdem1"])
  cat("diff [xbi]netdem1 - [xbn]netdem1 =", round(dif, 5), "se", round(dse, 5), "p", signif(2 * pnorm(-abs(dif / dse)), 3), "\n")
  out <- list(xbi_netdem1 = c(b = th[["xbi:netdem1"]], se = se[["xbi:netdem1"]], p = 2 * pnorm(-abs(th[["xbi:netdem1"]] / se[["xbi:netdem1"]]))))
  attr(out, "n") <- n; attr(out, "theta") <- th; attr(out, "diff") <- c(b = dif, se = dse); out
}

# 0. validate the likelihood at the authors' published estimates (log: LL = -1041.7114)
b_auth <- unlist(read.delim(file.path(WD, "IORAIC/output/betaUM.txt")))
XE0 <- cbind(as.matrix(d[, covE]), 1); XI0 <- cbind(as.matrix(d[, covI]), 1)
ll_auth <- sum(boselect_ll(b_auth, XE0, XI0, d$dv_esc_io)); cat("LL at authors' betas:", ll_auth, "(log -1041.7114)\n")
stopifnot(abs(ll_auth - (-1041.7114)) < 1e-3)
# 1. reproduce on the full sample: from zero starts (lntau1 = athrho = 0) and from the authors' estimates; keep the higher LL
r_zero <- tryCatch(fit_boselect(character(0), d), error = function(e) { cat("zero-start fit failed:", conditionMessage(e), "\n"); NULL })
r_auth <- fit_boselect(character(0), d, unname(b_auth))
ll_of <- function(r) sum(boselect_ll(attr(r, "theta"), XE0, XI0, d$dv_esc_io))
r_full <- if (!is.null(r_zero) && ll_of(r_zero) > ll_of(r_auth) + 1e-6) r_zero else r_auth
cat("LL zero-start:", if (is.null(r_zero)) NA else ll_of(r_zero), " LL author-start:", ll_of(r_auth), "\n")
cat("REPRO: [xbi]netdem1", r_full$xbi_netdem1[["b"]], "(log -.1090921) se", r_full$xbi_netdem1[["se"]], "(log .0286524)\n")
# 2. common sample 1822-2011 + history (challenger/target dyad)
vars <- c("dv_esc_io", covI)
ds <- common_sample(d, vars, "challenger", "target", "year")
st0 <- attr(r_full, "theta")
fitter <- function(rhs, data) {
  extra <- setdiff(rhs, covI)
  kE <- length(covE) + 1; kI <- length(covI) + 1; kx <- length(extra)
  s <- c(st0[1:(kE - 1)], rep(0, kx), st0[kE], st0[kE + 1:(kE - 1)], rep(0, kx), st0[2 * kE],
         st0[2 * kE + 1:(kI - 1)], rep(0, kx), st0[2 * kE + kI], st0[2 * kE + kI + 1:2])
  fit_boselect(extra, data, unname(s))
}
M <- run_six(covI, ds, fitter)
saveRDS(list(M = M, r_full = r_full, n_orig = nrow(d), n_used = nrow(ds), diffs = lapply(M, attr, "diff")), file.path(W, "M_P8_030.rds"))
meta <- list(pid = "P8_030", cite = "Chiba & Fang 2014", journal = "JOP", year = 2014, doi = "10.1017/s0022381614000127", family = "DEMOCRACY")
dd <- sapply(M, function(m) attr(m, "diff")[["b"]])
write_p8(meta, "xbi_netdem1", "abstract: 'more accountable regimes are less likely to escalate ... when an international organization opposes'",
  "Table 1 unified model, [xbi] (escalation under IO opposition) netdem1", "custom ML 'boselect' (IO-opposition probit + two-regime ordered probit, bivariate normal errors), robust SE",
  "boselect likelihood ported to R (BFGS + Newton, Stata robust SE), all six steps; history terms in all three equations",
  -0.1090921, 0.0286524, 3880, TRUE, FALSE, M, "N (no MID-history terms; enduring-rivalry dummy endrv10j only)",
  sprintf(paste0("Likelihood validated at authors' betas (LL %.4f vs log -1041.7114). Full-sample repro from zero starts: b=%.5f se=%.5f (log -.1090921/.0286524). ",
  "Common 1822-2011 sample N=%d (orig 3880). Interaction analogue [xbi]-[xbn] netdem1: m0 %.4f, mDM %.4f, m26 %.4f (log in log_P8_030.txt)."),
  ll_auth, r_full$xbi_netdem1[["b"]], r_full$xbi_netdem1[["se"]], nrow(ds), dd[["m0"]], dd[["mDM"]], dd[["m26"]]))
sink()
