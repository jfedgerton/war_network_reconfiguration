# Newey (1987) two-step IV probit (minimum-chi-square / AGLS), as Stata's ivprobit, twostep - hand-coded (ivprobit is archived on CRAN).
#   y1 binary; y2 one continuous endogenous regressor; X1 included exogenous (no constant; added); X2 excluded instruments.
#   Step 1: OLS y2 on Z = [1, X1, X2] -> Pi, v. Step 2: probit y1 on [Z, v] -> alpha (on Z), lambda (on v), J_aa^-1.
#   Step 3: probit y1 on [y2, 1, X1, v] -> beta_y. Omega = Var(alpha) + (lambda - beta_y)^2 Var(Pi). D = [Pi, S] (S selects 1 and X1 in Z).
#   delta = (D' Omega^-1 D)^-1 D' Omega^-1 alpha; Var(delta) = (D' Omega^-1 D)^-1. delta is beta / sd(u1 | v) (Stata twostep scale). Var(alpha) and Var(Pi) cluster-robust when cluster given.
ivprobit_newey <- function(y1, y2, X1, X2, cluster = NULL) {
  stopifnot(all(y1 %in% 0:1), is.numeric(y2), nrow(X1) == length(y1), nrow(X2) == length(y1))
  Z <- cbind(`(Intercept)` = 1, X1, X2); kz <- ncol(Z)
  vc_cl <- function(X, sc) { if (is.null(cluster)) crossprod(sc) else crossprod(rowsum(sc, cluster)) }
  # step 1
  f1 <- lm.fit(Z, y2); Pi <- f1$coefficients; v <- f1$residuals
  ZtZi <- solve(crossprod(Z)); VPi <- ZtZi %*% vc_cl(Z, Z * v) %*% ZtZi
  # step 2 (reduced-form probit with v)
  R <- cbind(Z, v = v); p2 <- glm.fit(R, y1, family = binomial("probit")); stopifnot(p2$converged)
  eta <- as.vector(R %*% p2$coefficients); q <- 2 * y1 - 1; lam_g <- q * dnorm(q * eta) / pnorm(q * eta)   # generalized residual
  Hinv <- solve(crossprod(R, R * (lam_g * (lam_g + eta))))                                             # inverse of the (expected-form) information
  Va <- (Hinv %*% vc_cl(R, R * lam_g) %*% Hinv)[1:kz, 1:kz]
  alpha <- p2$coefficients[1:kz]; lambda <- p2$coefficients[kz + 1]
  # step 3 (2SIV beta_y)
  R3 <- cbind(y2 = y2, `(Intercept)` = 1, X1, v = v); p3 <- glm.fit(R3, y1, family = binomial("probit")); stopifnot(p3$converged)
  by <- p3$coefficients[1]
  Om <- Va + (lambda - by)^2 * VPi
  S <- diag(kz)[, 1:(1 + ncol(X1)), drop = FALSE]
  D <- cbind(Pi, S); Oi <- solve(Om)
  Vd <- solve(t(D) %*% Oi %*% D); delta <- as.vector(Vd %*% t(D) %*% Oi %*% alpha)
  nm <- c("y2", colnames(Z)[1:(1 + ncol(X1))])
  list(coef = setNames(delta, nm), se = setNames(sqrt(diag(Vd)), nm), n = length(y1), first_stage_F = NA)
}
ivprobit_selftest <- function() {
  set.seed(123); n <- 30000
  x1 <- rnorm(n); z1 <- rnorm(n); z2 <- rnorm(n); e <- MASS::mvrnorm(n, c(0, 0), matrix(c(1, .5, .5, 1), 2))
  y2 <- 0.5 * x1 + 0.8 * z1 - 0.6 * z2 + e[, 2]; y1 <- as.integer(-0.2 + 0.7 * y2 - 0.5 * x1 + e[, 1] > 0)
  f <- ivprobit_newey(y1, y2, cbind(x1 = x1), cbind(z1 = z1, z2 = z2))
  tru <- c(.7, -.2, -.5) / sqrt(1 - .5^2)   # Newey two-step identifies beta / sd(u1 | v) (as Stata's twostep)
  cat("ivprobit self-test: true/sd(u1|v)", round(tru, 3), "| est", round(f$coef, 3), "se", round(f$se, 3), "\n")
  stopifnot(max(abs(f$coef - tru)) < 0.06)
  invisible(TRUE)
}
