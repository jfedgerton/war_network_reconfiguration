# Multinomial logit by Newton-Raphson with Stata-style cluster-robust SEs (G/(G-1)), base outcome = lowest category.
mlogit_fit <- function(y, X, cluster = NULL, maxit = 60, tol = 1e-9) {
  X <- as.matrix(as.data.frame(X)); nm <- colnames(X)
  s <- apply(X, 2, sd); s[s == 0] <- 1
  Xs <- cbind(`(Intercept)` = 1, sweep(X, 2, s, "/"))
  lev <- sort(unique(y)); stopifnot(length(lev) >= 2)
  J <- length(lev) - 1; n <- nrow(Xs); p <- ncol(Xs)
  Y <- sapply(lev[-1], function(k) as.numeric(y == k))
  B <- matrix(0, p, J)
  ll_fun <- function(B) { E <- cbind(0, Xs %*% B); mx <- apply(E, 1, max); lse <- mx + log(rowSums(exp(E - mx)))
    sum(rowSums(cbind(1 - rowSums(Y), Y) * E) - lse) }
  probs <- function(B) { E <- cbind(0, Xs %*% B); E <- exp(E - apply(E, 1, max)); P <- E / rowSums(E); P[, -1, drop = FALSE] }
  ll <- ll_fun(B)
  for (it in 1:maxit) {
    P <- probs(B); G <- crossprod(Xs, Y - P)
    H <- matrix(0, p * J, p * J)
    for (j in 1:J) for (l in j:J) {
      w <- if (j == l) P[, j] * (1 - P[, j]) else -P[, j] * P[, l]
      blk <- crossprod(Xs, Xs * w)
      H[(j-1)*p + 1:p, (l-1)*p + 1:p] <- blk; H[(l-1)*p + 1:p, (j-1)*p + 1:p] <- blk
    }
    step <- matrix(solve(H, as.vector(G)), p, J)
    t <- 1; repeat { Bn <- B + t * step; lln <- ll_fun(Bn); if (is.finite(lln) && lln >= ll - 1e-10 || t < 1e-4) break; t <- t / 2 }
    conv <- abs(lln - ll) < tol * (abs(ll) + 1); B <- Bn; ll <- lln
    if (conv) break
  }
  stopifnot(it < maxit)
  P <- probs(B); Hinv <- solve(H)   # H at previous iterate; recompute at optimum
  H <- matrix(0, p * J, p * J)
  for (j in 1:J) for (l in j:J) { w <- if (j == l) P[, j] * (1 - P[, j]) else -P[, j] * P[, l]
    blk <- crossprod(Xs, Xs * w); H[(j-1)*p + 1:p, (l-1)*p + 1:p] <- blk; H[(l-1)*p + 1:p, (j-1)*p + 1:p] <- blk }
  Hinv <- solve(H)
  S <- do.call(cbind, lapply(1:J, function(j) Xs * (Y[, j] - P[, j])))   # n x pJ scores
  if (!is.null(cluster)) { Sg <- rowsum(S, cluster); Gc <- nrow(Sg); meat <- crossprod(Sg) * Gc / (Gc - 1) }
  else meat <- crossprod(S)
  V <- Hinv %*% meat %*% Hinv
  D <- rep(c(1, 1 / s), J)                                       # back-transform to original scale
  b <- as.vector(B) * D; V <- V * outer(D, D)
  rn <- paste0(rep(paste0("y", lev[-1]), each = p), ":", rep(c("(Intercept)", nm), J))
  se <- sqrt(diag(V)); z <- b / se
  ct <- cbind(Estimate = b, `Std. Error` = se, `z value` = z, `Pr(>|z|)` = 2 * pnorm(-abs(z))); rownames(ct) <- rn
  list(coeftable = ct, loglik = ll, n = n, iter = it, levels = lev)
}
mlogit_stat <- function(m, var, outcome) {
  r <- paste0("y", outcome, ":", var); stopifnot(r %in% rownames(m$coeftable))
  c(b = unname(m$coeftable[r, 1]), se = unname(m$coeftable[r, 2]), p = unname(m$coeftable[r, 4]))
}
