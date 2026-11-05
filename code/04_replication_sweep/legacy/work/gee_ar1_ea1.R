# GEE with AR(1) working correlation for arbitrary (link, variance) pairs -- used for Stata
#   xtgee y X, family(nbinomial k) link(nbinomial) corr(ar 1) vce(robust)      (Kinne 2012 JOP, estimator audit ea1)
# geepack cannot do a negative-binomial variance with fixed k, so the estimating equations are solved directly:
#   sum_i D_i' V_i^-1 (y_i - mu_i) = 0,  D_i = diag(dmu/deta) X_i,  V_i = A_i^1/2 R_i(rho) A_i^1/2,  A = diag(V(mu)),
#   R_i[s,t] = rho^|year_s - year_t|  (uses actual years, so gaps are handled as Stata's tsset-based AR structure).
#   rho: moment estimator on Pearson residuals, rho = [sum_{t,t+1 adjacent} r_it r_i,t+1 / n_pairs] / [sum r^2 / N]
#   (scale cancels; Stata fixes phi = 1 for nbinomial -- irrelevant for the robust SE).
#   Robust (sandwich) SE: M^-1 (sum_i U_i U_i') M^-1, no small-sample factor (Stata xtgee vce(robust) default).
fam_nb_canon <- function(k) {      # Stata link(nbinomial): eta = ln(mu / (mu + 1/k)), variance mu + k mu^2
  list(name = paste0("nbinomial(k=", signif(k, 5), ") canonical link"),
       linkfun = function(mu) log(mu / (mu + 1 / k)),
       linkinv = function(eta) (1 / k) * exp(eta) / (1 - exp(eta)),
       mu.eta  = function(eta) { mu <- (1 / k) * exp(eta) / (1 - exp(eta)); mu * (1 + k * mu) },
       variance = function(mu) mu + k * mu^2,
       valideta = function(eta) all(is.finite(eta)) && all(eta < 0))
}
fam_pois_log <- list(name = "poisson log", linkfun = log, linkinv = exp, mu.eta = exp, variance = function(mu) mu,
                     valideta = function(eta) all(is.finite(eta)))

gee_ar1 <- function(X, y, id, tt, fam, beta0, fixed_rho = NULL, tol = 1e-9, maxit = 200, verbose = FALSE) {
  o <- order(id, tt); X <- X[o, , drop = FALSE]; y <- y[o]; id <- id[o]; tt <- tt[o]
  stopifnot(!anyDuplicated(data.frame(id, tt)), all(is.finite(X)), all(is.finite(y)))
  N <- length(y); p <- ncol(X)
  idx <- split(seq_len(N), id)
  adj <- which(c(id[-1] == id[-N] & diff(tt) == 1, FALSE))      # (t, t+1) adjacent pairs
  beta <- beta0; rho <- 0
  for (it in 1:maxit) {
    eta <- as.vector(X %*% beta); stopifnot(fam$valideta(eta))
    mu <- fam$linkinv(eta); v <- fam$variance(mu); de <- fam$mu.eta(eta)
    r <- (y - mu) / sqrt(v)
    rho <- if (!is.null(fixed_rho)) fixed_rho else (sum(r[adj] * r[adj + 1]) / length(adj)) / (sum(r^2) / N)
    stopifnot(abs(rho) < 1)
    M <- matrix(0, p, p); U <- numeric(p)
    for (g in idx) {
      Rg <- rho^abs(outer(tt[g], tt[g], "-"))
      sv <- sqrt(v[g]); Vg <- Rg * outer(sv, sv)
      Dg <- X[g, , drop = FALSE] * de[g]
      Wi <- solve(Vg, cbind(Dg, y[g] - mu[g]))
      M <- M + crossprod(Dg, Wi[, 1:p, drop = FALSE]); U <- U + crossprod(Dg, Wi[, p + 1])
    }
    step <- as.vector(solve(M, U)); lam <- 1
    while (!fam$valideta(as.vector(X %*% (beta + lam * step))) && lam > 1e-6) lam <- lam / 2
    beta <- beta + lam * step
    if (verbose) cat("it", it, "rho", round(rho, 5), "max|step|", signif(max(abs(lam * step)), 3), "\n")
    if (max(abs(lam * step)) < tol) break
  }
  stopifnot(it < maxit)
  # sandwich at the solution (rho held at its final value)
  eta <- as.vector(X %*% beta); mu <- fam$linkinv(eta); v <- fam$variance(mu); de <- fam$mu.eta(eta)
  M <- matrix(0, p, p); B <- matrix(0, p, p)
  for (g in idx) {
    Rg <- rho^abs(outer(tt[g], tt[g], "-")); sv <- sqrt(v[g]); Vg <- Rg * outer(sv, sv)
    Dg <- X[g, , drop = FALSE] * de[g]
    Wi <- solve(Vg, cbind(Dg, y[g] - mu[g]))
    M <- M + crossprod(Dg, Wi[, 1:p, drop = FALSE]); Ug <- crossprod(Dg, Wi[, p + 1]); B <- B + Ug %*% t(Ug)
  }
  Mi <- solve(M); Vr <- Mi %*% B %*% Mi
  names(beta) <- colnames(X); dimnames(Vr) <- list(colnames(X), colnames(X))
  se <- sqrt(diag(Vr)); z <- beta / se
  list(coef = beta, se = se, p = 2 * pnorm(-abs(z)), V = Vr, rho = rho, iter = it, n = N, G = length(idx))
}
