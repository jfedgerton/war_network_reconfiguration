# King & Zeng (2001) rare-events logit finite-sample bias correction (no prior correction; no tau), as in Stata relogit / Zelig relogit:
#   W = diag(pi(1-pi)); Q_ii = x_i' (X'WX)^-1 x_i; xi_i = 0.5 * Q_ii * (2 pi_i - 1); bias = (X'WX)^-1 X'W xi; b_tilde = b_hat - bias.
#   Variance: V(b_tilde) = (n/(n+k))^2 V(b_hat)  (V(b_hat) = the model's own vcov, here cluster-robust as in relogit, cluster()).
# kz_bias(Xfun, n, pi, idx_chunks): Xfun(rows) returns the dense design matrix (incl. intercept or full FE dummies) for those rows.
kz_bias <- function(Xfun, n, pi, chunk = 100000) {
  stopifnot(length(pi) == n, all(pi > 0 & pi < 1))
  ch <- split(seq_len(n), ceiling(seq_len(n) / chunk))
  A <- NULL
  for (r in ch) { Xc <- Xfun(r); w <- pi[r] * (1 - pi[r]); A <- if (is.null(A)) crossprod(Xc, Xc * w) else A + crossprod(Xc, Xc * w) }
  Ainv <- solve(A); g <- 0
  for (r in ch) {
    Xc <- Xfun(r); w <- pi[r] * (1 - pi[r])
    Q <- rowSums((Xc %*% Ainv) * Xc)
    xi <- 0.5 * Q * (2 * pi[r] - 1)
    g <- g + crossprod(Xc, w * xi)
  }
  b <- drop(Ainv %*% g); names(b) <- colnames(A); b
}
# Self-test vs the Zelig reference implementation (lm.influence hat values) on simulated rare-event data; stops loudly on mismatch.
kz_selftest <- function() {
  set.seed(123); n <- 3000; x1 <- rnorm(n); x2 <- rbinom(n, 1, .3); g <- factor(sample(1:5, n, TRUE))
  y <- rbinom(n, 1, plogis(-3.5 + .5 * x1 + .8 * x2 + as.numeric(g) * .1))
  X <- model.matrix(~ x1 + x2 + g)
  m <- glm(y ~ X - 1, family = binomial); pih <- fitted(m); W <- pih * (1 - pih)
  Qd <- lm.influence(lm(y ~ X - 1, weights = W))$hat / W          # Zelig relogit code
  xi <- 0.5 * Qd * (2 * pih - 1)
  bias_ref <- coef(lm(xi ~ X - 1, weights = W))
  bias_new <- kz_bias(function(r) X[r, , drop = FALSE], n, pih, chunk = 700)
  stopifnot(max(abs(bias_ref - bias_new)) < 1e-8)
  # reparametrization invariance of slope bias (intercept+dummies vs full dummies without intercept)
  X2 <- cbind(x1 = x1, x2 = x2, model.matrix(~ g - 1))
  b2 <- kz_bias(function(r) X2[r, , drop = FALSE], n, pih, chunk = 700)
  stopifnot(max(abs(b2[c("x1", "x2")] - bias_ref[c("Xx1", "Xx2")])) < 1e-8)
  cat("kz_selftest OK: max|diff| =", signif(max(abs(bias_ref - bias_new)), 3), "\n"); invisible(TRUE)
}
