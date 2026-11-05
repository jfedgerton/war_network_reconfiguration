# Split-population (partial-observability) bivariate probit, Poirier 1980 / Xiang 2010 JOP:
#   P(y = 1) = Phi2(X b, Z g; rho)  (conflict equation X, relevance equation Z), y = 0 otherwise. Full ML, rho = 0.9 tanh(theta) (bounded inside the BVN branch).
# Bivariate normal CDF: Genz (2004) BVND algorithm, 20-point Gauss-Legendre branch (|rho| < 0.925 enforced), vectorised.
# Analytic gradient; OIM SEs from the numerical Jacobian of the analytic score (numDeriv). Fail-loud checks on synthetic data.
.gl20 <- local({ n <- 20; i <- 1:(n - 1); b <- i / sqrt(4 * i^2 - 1); J <- matrix(0, n, n); J[cbind(i, i + 1)] <- b; J[cbind(i + 1, i)] <- b
  e <- eigen(J, symmetric = TRUE); list(x = e$values, w = 2 * e$vectors[1, ]^2) })
pbvn <- function(a, b, r) {                       # P(X <= a, Y <= b), corr r (scalar), |r| < 0.925
  stopifnot(length(r) == 1, abs(r) < 0.925)
  h <- -a; k <- -b; hs <- (h * h + k * k) / 2; asr <- asin(r); s <- 0
  for (j in seq_along(.gl20$x)) { sn <- sin(asr * (1 + .gl20$x[j]) / 2); s <- s + .gl20$w[j] * exp((sn * h * k - hs) / (1 - sn^2)) }
  s * asr / (4 * pi) + pnorm(-h) * pnorm(-k)
}
dbvn <- function(a, b, r) exp(-(a^2 - 2 * r * a * b + b^2) / (2 * (1 - r^2))) / (2 * pi * sqrt(1 - r^2))
sp_ll <- function(th, X, Z, y, w = NULL) {
  kx <- ncol(X); kz <- ncol(Z); bx <- th[1:kx]; gz <- th[kx + 1:kz]; r <- 0.9 * tanh(th[kx + kz + 1])
  a <- as.vector(X %*% bx); b <- as.vector(Z %*% gz); P <- pmin(pmax(pbvn(a, b, r), 1e-300), 1 - 1e-15)
  ll <- ifelse(y == 1, log(P), log1p(-P)); sum(ll)
}
sp_score <- function(th, X, Z, y, by_obs = FALSE) {
  kx <- ncol(X); kz <- ncol(Z); bx <- th[1:kx]; gz <- th[kx + 1:kz]; r <- 0.9 * tanh(th[kx + kz + 1])
  a <- as.vector(X %*% bx); b <- as.vector(Z %*% gz); P <- pmin(pmax(pbvn(a, b, r), 1e-300), 1 - 1e-15); s1 <- sqrt(1 - r^2)
  dPa <- dnorm(a) * pnorm((b - r * a) / s1); dPb <- dnorm(b) * pnorm((a - r * b) / s1); dPr <- dbvn(a, b, r) * 0.9 * (1 - tanh(th[kx + kz + 1])^2)
  m <- ifelse(y == 1, 1 / P, -1 / (1 - P))
  S <- cbind(X * (m * dPa), Z * (m * dPb), m * dPr)
  if (by_obs) S else colSums(S)
}
fit_splitpop <- function(y, X, Z, start = NULL, cluster = NULL) {
  stopifnot(all(y %in% 0:1), nrow(X) == length(y), nrow(Z) == length(y), all(is.finite(X)), all(is.finite(Z)))
  if (is.null(start)) {
    p1 <- glm.fit(X, y, family = binomial("probit")); stopifnot(all(X[, 1] == 1), all(Z[, 1] == 1)); start <- c(p1$coefficients, 1, rep(0, ncol(Z) - 1), 0)   # relevance intercept 1, rho 0
    names(start) <- NULL
  }
  f <- function(t) -sp_ll(t, X, Z, y); g <- function(t) -sp_score(t, X, Z, y)
  o <- optim(start, f, g, method = "BFGS", control = list(maxit = 2000, reltol = 1e-12))
  o <- optim(o$par, f, g, method = "BFGS", control = list(maxit = 2000, reltol = 1e-14))
  stopifnot(o$convergence == 0)
  cat("split-pop ML: LL", -o$value, "par", paste(round(o$par, 4), collapse = " "), "\n")
  H <- numDeriv::jacobian(function(t) sp_score(t, X, Z, y), o$par)
  V <- tryCatch(solve(-H), error = function(e) stop("split-pop Hessian singular (boundary or unidentified parameter): ", conditionMessage(e)))
  if (!is.null(cluster)) { S <- rowsum(sp_score(o$par, X, Z, y, by_obs = TRUE), cluster); V <- V %*% crossprod(S) %*% V }
  nm <- c(paste0("out:", colnames(X)), paste0("rel:", colnames(Z)), "theta_rho")
  list(coef = setNames(o$par, nm), se = setNames(sqrt(diag(V)), nm), ll = -o$value, n = length(y), rho = 0.9 * tanh(o$par[length(o$par)]),
       grad_max = max(abs(sp_score(o$par, X, Z, y))))
}
splitpop_selftest <- function() {
  set.seed(123); n <- 40000
  x1 <- rnorm(n); x2 <- rnorm(n); z1 <- rnorm(n); e <- MASS::mvrnorm(n, c(0, 0), matrix(c(1, .3, .3, 1), 2))
  y <- as.integer((0.3 + 0.8 * x1 - 0.6 * x2 + e[, 1] > 0) & (0.5 + 0.7 * x1 - 0.9 * z1 + e[, 2] > 0))
  X <- cbind(`(Intercept)` = 1, x1 = x1, x2 = x2); Z <- cbind(`(Intercept)` = 1, x1 = x1, z1 = z1)
  th <- c(.2, .7, -.5, .4, .6, -.8, .2)
  gn <- numDeriv::grad(function(t) sp_ll(t, X, Z, y), th); ga <- sp_score(th, X, Z, y)
  stopifnot(max(abs(gn - ga)) < 1e-3 * max(1, max(abs(ga))))
  chk <- c(pbvn(0.3, -0.2, 0.5), mvtnorm::pmvnorm(upper = c(0.3, -0.2), corr = matrix(c(1, .5, .5, 1), 2))[1])
  stopifnot(abs(diff(chk)) < 1e-7)
  f <- fit_splitpop(y, X, Z)
  cat("splitpop self-test: true out (0.3, 0.8, -0.6) rel (0.5, 0.7, -0.9) rho 0.3 | est", round(f$coef[1:6], 3), "rho", round(f$rho, 3), "\n")
  stopifnot(max(abs(f$coef[1:6] - c(.3, .8, -.6, .5, .7, -.9))) < 0.2)
  invisible(TRUE)
}
