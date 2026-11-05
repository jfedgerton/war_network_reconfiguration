# Stata zinb (NB2 count + logit inflation) by ML in R, analytic scores, Stata-style cluster-robust V (G/(G-1)).
#   pi = invlogit(Z g); mu = exp(X b); alpha = exp(lnalpha); m = 1/alpha; f0 = (1 + alpha mu)^(-m)
#   y = 0: L = pi + (1 - pi) f0;   y > 0: L = (1 - pi) NB(y; mu, alpha)
# X, Z: model matrices WITH constant columns. Parameters ordered (b, g, lnalpha).
zinb_parts <- function(th, X, Z, y) {
  kx <- ncol(X); kz <- ncol(Z)
  b <- th[1:kx]; g <- th[kx + 1:kz]; la <- th[kx + kz + 1]; a <- exp(la); m <- 1 / a
  eta <- drop(X %*% b); mu <- exp(eta); zg <- drop(Z %*% g)
  pi <- plogis(zg); l1pa <- log1p(a * mu); q <- mu / (1 + a * mu)
  n <- length(y); ll <- numeric(n); Sb <- numeric(n); Sg <- numeric(n); Sa <- numeric(n)
  i0 <- y == 0; i1 <- !i0
  # y = 0
  lf0 <- -m * l1pa[i0]; f0 <- exp(lf0); p0 <- pi[i0]
  L0 <- p0 + (1 - p0) * f0; ll[i0] <- log(L0)
  w <- (1 - p0) * f0 / L0
  Sb[i0] <- w * (-q[i0])
  Sa[i0] <- w * (a * ((1 / a^2) * l1pa[i0] - q[i0] / a))
  Sg[i0] <- p0 * (1 - p0) * (1 - f0) / L0
  # y > 0
  yy <- y[i1]
  ll[i1] <- log1p(-pi[i1]) + lgamma(yy + m) - lgamma(m) - lgamma(yy + 1) - m * l1pa[i1] + yy * (log(a) + eta[i1]) - yy * l1pa[i1]
  Sb[i1] <- (yy - mu[i1]) / (1 + a * mu[i1])
  Sa[i1] <- a * ((digamma(yy + m) - digamma(m)) * (-1 / a^2) + (1 / a^2) * l1pa[i1] - m * q[i1] + yy / a - yy * q[i1])
  Sg[i1] <- -pi[i1]
  list(ll = ll, S = cbind(X * Sb, Z * Sg, Sa))
}

fit_zinb <- function(X, Z, y, cl = NULL, start = NULL) {
  stopifnot(all(y >= 0), all(y == round(y)), !anyNA(X), !anyNA(Z), !anyNA(y))
  sx <- apply(X, 2, sd); sx[sx == 0] <- 1; sz <- apply(Z, 2, sd); sz[sz == 0] <- 1
  Xs <- sweep(X, 2, sx, "/"); Zs <- sweep(Z, 2, sz, "/"); D <- c(1 / sx, 1 / sz, 1)
  kx <- ncol(X); kz <- ncol(Z)
  if (is.null(start)) {
    pg <- glm.fit(Xs, y, family = poisson())
    stopifnot(pg$converged)
    start <- c(pg$coefficients, rep(0, kz), 0)
  } else start <- start / D
  fn <- function(th) -sum(zinb_parts(th, Xs, Zs, y)$ll)
  gr <- function(th) -colSums(zinb_parts(th, Xs, Zs, y)$S)
  o <- optim(start, fn, gr, method = "BFGS", control = list(maxit = 5000, reltol = 1e-14))
  stopifnot(o$convergence == 0)
  # polish with Newton steps on the numerical Hessian of the analytic gradient
  for (it in 1:50) {
    H <- optimHess(o$par, fn, gr); gvec <- gr(o$par)
    step <- solve(H, gvec); o$par <- o$par - step
    if (max(abs(step)) < 1e-10) break
  }
  gfin <- gr(o$par); stopifnot(max(abs(gfin)) < 1e-4)
  H <- optimHess(o$par, fn, gr); Hi <- solve(H)
  S <- zinb_parts(o$par, Xs, Zs, y)$S
  if (is.null(cl)) V <- Hi else {
    Sg <- rowsum(S, cl); G <- nrow(Sg)
    V <- Hi %*% crossprod(Sg) %*% Hi * G / (G - 1)
  }
  th <- o$par * D; V <- V * outer(D, D)
  nm <- c(colnames(X), paste0("inflate:", colnames(Z)), "lnalpha")
  names(th) <- nm; dimnames(V) <- list(nm, nm)
  list(coef = th, V = V, ll = -fn(o$par), n = length(y), G = if (is.null(cl)) NA else length(unique(cl)))
}
