# Stata heckprob (probit with sample selection) by ML in R, analytic scores, OIM or cluster-robust SEs.
# Outcome y* = X b + e1 observed only if s = 1, selection s* = Z g + e2; corr(e1, e2) = rho = tanh(athrho).
#   s = 0: P = Phi(-Zg);  s = 1, y = q in {1,-1}: P = Phi2(q Xb, Zg, q rho)
# Bivariate normal CDF: bvn_lower() (vectorized Genz bvnu; work/55/biprobit.R).
source(file.path(CODE, "work/55/biprobit.R"))
hp_parts <- function(th, X, Z, y, s) {
  kx <- ncol(X); kz <- ncol(Z)
  b <- th[1:kx]; g <- th[kx + 1:kz]; rho <- tanh(th[kx + kz + 1]); rho <- sign(rho) * min(abs(rho), 1 - 1e-10)
  zg <- drop(Z %*% g); n <- length(s)
  ll <- numeric(n); Sx <- matrix(0, n, kx); Sz <- matrix(0, n, kz); Sr <- numeric(n)
  i0 <- s == 0
  P0 <- pmax(pnorm(-zg[i0]), 1e-300); ll[i0] <- log(P0); Sz[i0, ] <- Z[i0, , drop = FALSE] * (-dnorm(zg[i0]) / P0)
  for (q in c(1, -1)) {
    i <- s == 1 & (2 * y - 1) == q; if (!any(i)) next
    xb <- drop(X[i, , drop = FALSE] %*% b); w1 <- q * xb; w2 <- zg[i]; r <- q * rho; sr <- sqrt(1 - r^2)
    P <- pmax(bvn_lower(w1, w2, r), 1e-300); ll[i] <- log(P)
    g1 <- dnorm(w1) * pnorm((w2 - r * w1) / sr) / P
    g2 <- dnorm(w2) * pnorm((w1 - r * w2) / sr) / P
    phi2 <- exp(-(w1^2 - 2 * r * w1 * w2 + w2^2) / (2 * (1 - r^2))) / (2 * pi * sr)
    Sx[i, ] <- X[i, , drop = FALSE] * (g1 * q); Sz[i, ] <- Z[i, , drop = FALSE] * g2
    Sr[i] <- q * phi2 / P * (1 - rho^2)
  }
  list(ll = ll, S = cbind(Sx, Sz, Sr))
}
# X, Z: model matrices WITH constant columns. cl: optional cluster vector -> Stata-style cluster-robust V.
fit_heckprob <- function(X, Z, y, s, cl = NULL, start = NULL) {
  stopifnot(all(s %in% 0:1), all(is.na(y[s == 0]) | y[s == 0] %in% 0:1), all(y[s == 1] %in% 0:1))
  y[s == 0] <- 0
  sx <- apply(X, 2, sd); sx[sx == 0] <- 1; sz <- apply(Z, 2, sd); sz[sz == 0] <- 1
  Xs <- sweep(X, 2, sx, "/"); Zs <- sweep(Z, 2, sz, "/"); D <- c(1 / sx, 1 / sz, 1)
  if (is.null(start)) {
    g0 <- coef(glm(s ~ Zs - 1, family = binomial("probit")))
    b0 <- coef(glm(y[s == 1] ~ Xs[s == 1, ] - 1, family = binomial("probit")))
    start <- c(b0, g0, 0)
  } else start <- start / D
  cache <- new.env(); cache$th <- NULL
  getp <- function(th) { if (!identical(th, cache$th)) { cache$th <- th; cache$p <- hp_parts(th, Xs, Zs, y, s) }; cache$p }
  o <- nlminb(start, function(t) -sum(getp(t)$ll), function(t) -colSums(getp(t)$S), control = list(eval.max = 5000, iter.max = 3000, rel.tol = 1e-14))
  th <- o$par
  hfun <- function(t) { H <- numDeriv::jacobian(function(u) colSums(hp_parts(u, Xs, Zs, y, s)$S), t, method.args = list(r = 2)); (H + t(H)) / 2 }
  for (it in 1:20) {
    g <- colSums(getp(th)$S); if (max(abs(g)) < 1e-7) break
    step <- solve(-hfun(th), g); lam <- 1; ll0 <- sum(getp(th)$ll)
    while (sum(getp(th + lam * step)$ll) < ll0 - 1e-10 && lam > 1e-4) lam <- lam / 2
    th <- th + lam * step
  }
  P <- getp(th); gmax <- max(abs(colSums(P$S))); cat("heckprob nlminb:", o$message, " LL", sum(P$ll), " Newton", it, " max|score|", signif(gmax, 3), "\n")
  stopifnot(gmax < 1e-4)
  Ainv <- solve(-hfun(th))
  if (is.null(cl)) V <- Ainv else {
    Sc <- rowsum(P$S, cl); G <- nrow(Sc); V <- Ainv %*% crossprod(Sc) %*% Ainv * G / (G - 1)
  }
  coef <- th * D; V <- V * outer(D, D)
  nm <- c(paste0("out:", colnames(X)), paste0("sel:", colnames(Z)), "athrho"); names(coef) <- nm; dimnames(V) <- list(nm, nm)
  list(coef = coef, V = V, se = sqrt(diag(V)), ll = sum(P$ll), n = length(s), n_sel = sum(s), rho = tanh(coef[["athrho"]]))
}
