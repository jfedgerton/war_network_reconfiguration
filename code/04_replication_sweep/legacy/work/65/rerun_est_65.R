# ESTAUDIT 65 (phase 6): copy of rerun_65.R. ONLY change: the conflict-equation variance is the cdsimeq
# Amemiya/Maddala (1983: 244-245) correction that the article reports (Table III note: "Standard errors from Maddala procedure"),
# implemented line-by-line from Keshk's cdsimeq.ado (SJ3-2 st0038, fetched to work/65/cdsimeq.ado), replacing Murphy-Topel.
#   V(alpha2) = D^-1 + d * D^-1 G' V0^-1 (X'X)^-1 V0^-1 G D^-1,  D = G' V0^-1 G,  G = [pi1, J2],
#   pi1 = stage-1 OLS coefs (lnrtrade on all exog), V0 = OIM covariance of stage-1 probit (fatal on all exog),
#   d = gamma2^2*sigma11 - 2*gamma2*sigma12, sigma11 = RSS/df, sigma12 = mean(fatal * v1 / phi(X pi2)).
# Point estimates, samples and the six steps are unchanged. Output -> results_estaudit/65__lnrtrade.csv (results/65.csv untouched).
# Rerun 65: Hegre, Oneal & Russett 2010 JPR "Trade does promote peace"
# Headline: Table 3 cols 3-4 (cdsimeq 2SPLS; state & year dummies in trade eq). Tie = I_lnrtrade in the fatal-MID probit eq.
# cdsimeq point estimates re-implemented: stage 1 OLS of lnrtrade on ALL exogenous vars (union of both eqs),
# stage 2 probit of fatal on fitted lnrtrade + conflict-eq exog. SEs: Murphy-Topel two-step correction
# (cdsimeq uses Amemiya/Maddala correction; asymptotically equivalent target, finite-sample differs).
ROOT <- Sys.getenv("SWEEP_DATA"); OUT <- Sys.getenv("SWEEP_OUT"); CODE <- Sys.getenv("SWEEP_CODE")
stopifnot(nzchar(ROOT), nzchar(OUT), nzchar(CODE))
source(file.path(CODE, "work/rerun_helpers_b118.R"))
source(file.path(CODE, "work/65/prep_65.R"))
WO <- file.path(OUT, "work/65"); dir.create(WO, recursive = TRUE, showWarnings = FALSE)  # port: outputs
W <- file.path(ROOT, "work/65"); logf <- file.path(WO, "log_est_65.txt"); cat("", file = logf)
sink(logf, split = TRUE, append = TRUE)

# ---- 1. Build authors' analysis data (HOR_JPR_table3.do) ----
d0 <- raw65()
# btscs py: reproduction check showed btscs skips rows with missing event (ongoing-dispute years coded missing):
# count runs over non-missing-fatal rows only (py missing where fatal missing). This matches the log far better
# (naive probit LL -2556.8 vs log -2557.24) than counting through missing rows (LL -2567.9).
ok <- !is.na(d0$fatal); first <- c(TRUE, d0$dyadid[-1] != d0$dyadid[-nrow(d0)])
py <- rep(NA_real_, nrow(d0)); lastev <- FALSE; cur <- -1
for (i in seq_len(nrow(d0))) { if (first[i]) { cur <- -1; lastev <- FALSE }
  if (ok[i]) { cur <- if (cur < 0 || lastev) 0 else cur + 1; py[i] <- cur; lastev <- d0$fatal[i] == 1 } }
d0$py <- py
# btscs nspline(3) knots: exact btscs.ado unavailable; basis form verified on btscs output elsewhere;
# knots set by rule (quartiles of py, Stata type-2), T = max(py). Not tuned to the log.
kn <- quantile(d0$py, c(.25, .5, .75), type = 2, na.rm = TRUE); Tm <- max(d0$py, na.rm = TRUE)
cat("spline knots:", kn, " T:", Tm, "\n")
d0 <- add_splines(d0, kn, Tm)
d <- build65(d0); rm(d0); invisible(gc())
stopifnot(nrow(d) == 527571)
names(d) <- sub("^_spline", "spl", names(d))
# Year dummies (Iyear_2..51 = 1951..2000) and the authors' state dummies (only the 164 codes listed in the do-file)
# are built inside cds() as SPARSE columns (memory: shared machine).
st_list <- scan(file.path(W, "dummies_do.txt"), quiet = TRUE); stopifnot(length(st_list) == 164)

conf_x <- c("smldem","lrgdem","dircont","logdstab","pwin_cap","allies","lnlrgcap","py","spl1","spl2","spl3","lnNt")
# ---- 2a. Check naive probit vs log (pre-drop sample, N=432524) ----
m_n <- feglm(as.formula(paste("fatal ~ lnrtrade +", paste(conf_x, collapse = "+"))), d, family = binomial("probit"), notes = FALSE)
cat("\nNaive probit (log: lnrtrade -.0573319, N 432524, LL -2557.2444):\n"); print(coeftable(m_n)[, 1:2]); cat("N", nobs(m_n), "LL", logLik(m_n), "\n")
naive_b <- unname(coef(m_n)["lnrtrade"]); rm(m_n)

# drop dyads involving the 36 states with no fatal disputes (as in do-file); state2 dummy dropped (baseline)
drop_st <- c(41,42,52,110,140,145,165,205,225,305,317,349,359,366,367,369,370,375,380,411,420,434,516,522,553,570,572,580,581,590,712,760,780,830,910,950)
d <- d[!(d$statea %in% drop_st | d$stateb %in% drop_st), ]
st_keep <- setdiff(st_list, c(2, drop_st)); stopifnot(length(st_keep) == 127)
trade_x <- c("lnsmlrgdp","lnlrgrgdp","lnsmlpop","lnlrgpop","logdstab","dircont","jtdem","allies","PTA","s_wt_glo",
             "py","spl1","spl2","spl3","lnNt")
m_n2 <- feglm(as.formula(paste("fatal ~ lnrtrade +", paste(conf_x, collapse = "+"))), d, family = binomial("probit"), notes = FALSE)
cat("Naive probit post-drop (log: -.0587106, N 279343):", coef(m_n2)["lnrtrade"], nobs(m_n2), "\n"); rm(m_n2)
suppressPackageStartupMessages(library(Matrix))

# ---- 2b. cdsimeq conflict equation ----
cds <- function(data, conf_rhs, trade_rhs, fe = TRUE) {
  exog <- unique(c(trade_rhs, conf_rhs))
  data <- data[complete.cases(data[, c("fatal", "lnrtrade", exog)]), ]
  n <- nrow(data)
  Zc <- as.matrix(data[, exog]); sc <- apply(abs(Zc), 2, max); sc[sc == 0] <- 1
  Zc <- sweep(Zc, 2, sc, "/")                      # column scaling for numerically stable normal equations
  Z <- cbind(Matrix(1, n, 1, sparse = TRUE), Matrix(Zc, sparse = TRUE))
  colnames(Z) <- c("(Intercept)", exog)
  if (fe) {
    iy <- which(data$year >= 1951)
    Dy <- sparseMatrix(i = iy, j = data$year[iy] - 1950, x = 1, dims = c(n, 50))
    ja <- match(data$statea, st_keep); jb <- match(data$stateb, st_keep)
    Ds <- sparseMatrix(i = c(which(!is.na(ja)), which(!is.na(jb))), j = c(ja[!is.na(ja)], jb[!is.na(jb)]), x = 1, dims = c(n, length(st_keep)))
    stopifnot(max(Ds) == 1)
    zn0 <- colnames(Z)
    Z <- cbind(Z, Dy, Ds)
    colnames(Z) <- c(zn0, paste0("Iyear_", 1951:2000), paste0("state_", st_keep))
  }
  ZtZ <- as.matrix(crossprod(Z)); keep <- diag(ZtZ) > 0
  Z <- Z[, keep]; ZtZ <- ZtZ[keep, keep]
  ch <- suppressWarnings(chol(ZtZ, pivot = TRUE)); rk <- attr(ch, "rank")
  if (rk < ncol(ZtZ)) { dropc <- attr(ch, "pivot")[(rk + 1):ncol(ZtZ)]; Z <- Z[, -dropc]; ZtZ <- ZtZ[-dropc, -dropc] }
  k1 <- ncol(Z); ZtZi <- chol2inv(chol(ZtZ))
  b1 <- ZtZi %*% as.vector(crossprod(Z, data$lnrtrade))
  data$I_lnrtrade <- as.vector(Z %*% b1); e <- data$lnrtrade - data$I_lnrtrade
  m2 <- feglm(as.formula(paste("fatal ~ I_lnrtrade +", paste(conf_rhs, collapse = "+"))), data,
              family = binomial("probit"), notes = FALSE)
  stopifnot(nobs(m2) == n)
  X2 <- cbind(1, as.matrix(data[, c("I_lnrtrade", conf_rhs)]))
  b2 <- coef(m2); xb <- drop(X2 %*% b2); P <- pnorm(xb); P <- pmin(pmax(P, 1e-12), 1 - 1e-12)
  lam <- dnorm(xb) * (data$fatal - P) / (P * (1 - P))
  V2 <- vcov(m2, vcov = "iid"); g <- b2["I_lnrtrade"]
  s2 <- sum(e^2) / (n - k1); V1 <- s2 * ZtZi
  C <- as.matrix(crossprod(X2 * (lam^2 * g), Z))   # sum score2 * dl2/dtheta1'
  R <- as.matrix(crossprod(X2 * (lam * e / s2), Z)) # sum score2 * score1'
  M <- C %*% V1 %*% t(C) - R %*% V1 %*% t(C) - C %*% V1 %*% t(R)
  V <- V2 + V2 %*% M %*% V2
  se_mt <- sqrt(diag(V)); names(se_mt) <- names(b2)
  # ---- cdsimeq Maddala correction (probit equation) ----
  # stage-1 probit reduced form: fatal on all exogenous Z (Newton-Raphson, observed information = Stata probit default OIM)
  y <- data$fatal; q <- 2 * y - 1
  pi2 <- numeric(ncol(Z)); pi2[1] <- qnorm(mean(y))
  for (it in 1:100) {
    eta <- as.vector(Z %*% pi2)
    lam <- q * exp(dnorm(q * eta, log = TRUE) - pnorm(q * eta, log.p = TRUE))
    g <- as.vector(crossprod(Z, lam))
    Hn <- as.matrix(crossprod(Z, Z * (lam * (lam + eta))))   # minus Hessian
    step <- solve(Hn, g); pi2 <- pi2 + step
    if (max(abs(step)) < 1e-10) break
  }
  stopifnot(it < 100, all(is.finite(pi2)))
  eta <- as.vector(Z %*% pi2)
  lam <- q * exp(dnorm(q * eta, log = TRUE) - pnorm(q * eta, log.p = TRUE))
  V0 <- solve(as.matrix(crossprod(Z, Z * (lam * (lam + eta)))))
  sigma11 <- s2
  sigma12 <- sum(y * e / dnorm(eta)) / n
  gamma2 <- unname(b2["I_lnrtrade"])
  MA_d <- gamma2^2 * sigma11 - 2 * gamma2 * sigma12
  # G = [pi1, J2]: J2 selects the conflict-eq exogenous columns of Z (+ constant, placed last as in Stata)
  cz <- colnames(Z); stopifnot(length(cz) == ncol(Z), !anyDuplicated(cz), all(c(conf_rhs, "(Intercept)") %in% cz))
  J2 <- sapply(c(conf_rhs, "(Intercept)"), function(v) as.numeric(cz == v))
  stopifnot(all(colSums(J2) == 1))
  G <- cbind(as.vector(b1), J2)
  V0i <- solve(V0)
  Dm <- solve(t(G) %*% V0i %*% G)
  Vm <- Dm + MA_d * Dm %*% (t(G) %*% V0i %*% ZtZi %*% V0i %*% G) %*% Dm
  se_mad <- sqrt(diag(Vm))[1]   # element 1 = I_lnrtrade (invariant to the column scaling of Z)
  cat(sprintf("  [maddala] probit RF iters %d; sigma11 %.5g sigma12 %.5g gamma2 %.5g d %.5g | se_Maddala %.6g se_MT %.6g se_naive %.6g\n",
              it, sigma11, sigma12, gamma2, MA_d, se_mad, se_mt["I_lnrtrade"], sqrt(V2["I_lnrtrade", "I_lnrtrade"])))
  se <- se_mt; se["I_lnrtrade"] <- se_mad   # only the tie SE is used downstream
  out <- list(b = b2, se = se, n = n, k1 = k1,
       coeftable = cbind(b = b2, se_MT = se_mt, se_naive = sqrt(diag(V2)), p_MT = 2 * pnorm(-abs(b2 / se_mt))))
  rm(Z, data, m2); invisible(gc()); out
}
r_rep <- cds(d, conf_x, trade_x)
cat("\n==== REPRODUCTION cdsimeq Table 3 cols 3-4 (log: I_lnrtrade -.0747785 se .0145783) ====\n")
print(r_rep$coeftable); cat("N =", r_rep$n, " stage-1 columns:", r_rep$k1, "\n")
b_pub <- -.0747785
b_repro <- unname(r_rep$b["I_lnrtrade"]); p_repro <- unname(r_rep$coeftable["I_lnrtrade", "p_MT"])
rel <- abs(b_repro - b_pub) / abs(b_pub)
reproduced <- ifelse(rel <= 0.10 & p_repro < .05, "REPRODUCED", "NOT_REPRODUCED")
cat(sprintf("b_pub=%g b_repro=%g rel.diff=%.4f p=%.3g -> %s\n", b_pub, b_repro, rel, p_repro, reproduced))
# also cols 1-2 as a check (no dummies): log -.0875165
r_c1 <- cds(d, conf_x, trade_x, fe = FALSE)
cat("Cols 1-2 check (log -.0875165 se .0163197):", r_c1$b["I_lnrtrade"], r_c1$se["I_lnrtrade"], "\n")
n_orig <- r_rep$n

# ---- 3. History ----
d <- add_history(d, "statea", "stateb", "year"); d$year_h <- d$year
# ---- 4. Battery (history vars enter the conflict eq and hence the stage-1 instrument set) ----
ds <- battery_sample(d, unique(c("fatal", "lnrtrade", conf_x, trade_x, "year", "statea", "stateb"))); rm(d); invisible(gc())
fitter <- function(rhs, data) { extra <- setdiff(rhs, conf_x); cds(data, rhs, c(trade_x, extra)) }
stat <- function(r) c(b = unname(r$b["I_lnrtrade"]), se = unname(r$se["I_lnrtrade"]),
                      p = unname(2 * pnorm(-abs(r$b["I_lnrtrade"] / r$se["I_lnrtrade"]))))
summary_model <- function(m) m$coeftable[!grepl("^Iyear_|^state", rownames(m$coeftable)), ]
R <- run_battery(ds, conf_x, fitter, stat, logf)

# single-equation probit robustness (not the headline): same battery with observed lnrtrade
fit_p <- function(rhs, data) fit_glm(c("lnrtrade", rhs), "fatal", data, "probit", NULL)
Rp <- lapply(list(m0 = character(0), mDM = c("hist_dyad_w15", "hist_mono_w15")),
             function(a) tie_stats(fit_p(c(conf_x, a), ds), "lnrtrade"))
cat("Single-eq probit lnrtrade: m0", Rp$m0, " mDM", Rp$mDM, "\n")

f_old <- file.path(OUT, "results/65.csv")  # port: prior rerun_65.R output (not ported); consistency check + notes only
old <- if (file.exists(f_old)) read.csv(f_old) else data.frame(b0 = NA, p0 = NA, p_dm = NA, p_26 = NA)
if (file.exists(f_old)) stopifnot(abs(R$m0["b"] - old$b0) < 1e-6)   # point estimates unchanged; only the SE correction differs
source(file.path(CODE, "ea4/est_helpers.R"))
st <- "2SPLS (cdsimeq re-impl.) probit conflict eq, Maddala/Amemiya SE"
notes <- sprintf(paste0("SE-only fix (Murphy-Topel -> cdsimeq Maddala). Repro Tbl3 cols3-4: b %.5f (log -.0747785), se_Maddala %.5f (log .0145783; MT was %.5f); ",
  "cols1-2 check b %.5f se_Maddala %.5f (log -.0875165/.0163197). Old MT results: p0 %.3g, p_dm %.3g, p_26 %.3g. ",
  "btscs py/splines rebuilt as before (knots %s). Hist vars in both stages. Single-eq probit lnrtrade (not headline): %.4f -> %.4f."),
  b_repro, r_rep$se["I_lnrtrade"], r_rep$coeftable["I_lnrtrade", "se_MT"], r_c1$b["I_lnrtrade"], r_c1$se["I_lnrtrade"],
  old$p0, old$p_dm, old$p_26, paste(kn, collapse = "/"), Rp$m0["b"], Rp$mDM["b"])
write_est(65, "Hegre, Oneal & Russett 2010", "lnrtrade",
  "two-stage probit least squares (Maddala 1983; Keshk cdsimeq), conflict eq = probit with instrumented ln trade; Maddala-corrected SE; peace years + 3 splines; year & state FE in trade eq",
  "pdf Table III (p.769) + note 'Standard errors from Maddala procedure'; authors' HOR_JPR_table3.do/log",
  rep(st, 6), "N",
  "all six steps used Murphy-Topel two-step SEs instead of the cdsimeq Maddala/Amemiya correction reported in the article; point estimates (2SPLS) already matched",
  "Y", R, nrow(ds), notes)
sink()
