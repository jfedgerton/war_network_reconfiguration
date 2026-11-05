# King & Zeng (2001) rare-events logit ("relogit") for the estimator audit (phase 6, agent ea1).
# Implemented directly in R (no Zelig on this box), following the Zelig/Stata relogit formulas:
#   1. ML logit  -> beta_hat, pi_hat.
#   2. W = diag(pi_hat (1 - pi_hat)); Q_ii = x_i' (X'WX)^-1 x_i;  xi_i = 0.5 * Q_ii * ((1 + w1) pi_hat_i - w1), w1 = 1
#      (no case-control weighting; with prior correction w1 is also 1 in Zelig/Stata).
#   3. bias = (X'WX)^-1 X'W xi  (= weighted LS of xi on X with weights W);  beta_tilde = beta_hat - bias.
#   4. V(beta_tilde) = (n / (n + k))^2 * V(beta_hat), applied to the SE type the article used (cluster-robust here).
#   5. Prior correction (only if tau supplied): intercept -= log(((1 - tau) / tau) * (ybar / (1 - ybar))); slopes unchanged.
# Fixed effects: fixest drops FE groups with all-0/all-1 outcomes (as Stata's logit/relogit drop perfect predictors); X for
#   the bias term = [covariates, full one-hot FE dummies] on the retained sample. The bias for slopes is invariant to the FE
#   parameterisation (depends only on the column space of X).
suppressPackageStartupMessages(library(fixest))

relogit_kz <- function(f, data, cluster = NULL, fe = NULL, tau = NULL, vcov_type = NULL) {
  ftxt <- paste(deparse(f, width.cutoff = 500), collapse = " ")
  f2 <- if (is.null(fe)) f else as.formula(paste(ftxt, "|", fe))
  vc <- if (!is.null(vcov_type)) vcov_type else if (is.null(cluster)) "iid" else as.formula(paste0("~", cluster))
  m <- feglm(f2, data = data, family = binomial("logit"), vcov = vc, ssc = ssc(adj = FALSE, cluster.adj = TRUE),
             notes = FALSE, nthreads = 1)
  stopifnot(m$convStatus)
  b_hat <- coef(m); V_hat <- vcov(m)
  keep <- fixest::obs(m)                           # row indices of data used in estimation
  dk <- data[keep, , drop = FALSE]
  X <- model.matrix(f, data = dk)                  # includes intercept when no FE
  if (!is.null(fe)) {
    X <- X[, colnames(X) != "(Intercept)", drop = FALSE]
    for (g in trimws(strsplit(fe, "\\+")[[1]])) {
      gg <- factor(dk[[g]])
      D <- Matrix::sparse.model.matrix(~ gg - 1)
      if (g != trimws(strsplit(fe, "\\+")[[1]])[1]) D <- D[, -1, drop = FALSE]   # avoid collinearity across multiple FEs
      X <- cbind(X, as.matrix(D))
    }
  }
  stopifnot(nrow(X) == nobs(m), all(names(b_hat) %in% colnames(X)))
  pi_hat <- as.numeric(fitted(m, type = "response")); stopifnot(length(pi_hat) == nrow(X))
  w <- pi_hat * (1 - pi_hat)
  # column scaling for numerical conditioning (cubic peace years ~1e6): Q_ii is invariant; bias_orig = bias_scaled / sc
  sc <- apply(X, 2, function(z) max(abs(z))); sc[sc == 0] <- 1
  X <- sweep(X, 2, sc, "/")
  ch <- split(seq_len(nrow(X)), ceiling(seq_len(nrow(X)) / 1e5))   # row chunks (memory)
  XtWX <- Reduce(`+`, lapply(ch, function(r) crossprod(X[r, , drop = FALSE], X[r, , drop = FALSE] * w[r])))
  A <- solve(XtWX)
  Qd <- unlist(lapply(ch, function(r) { Xr <- X[r, , drop = FALSE]; rowSums((Xr %*% A) * Xr) }), use.names = FALSE)
  xi <- 0.5 * Qd * (2 * pi_hat - 1)
  bias_all <- as.numeric(A %*% Reduce(`+`, lapply(ch, function(r) crossprod(X[r, , drop = FALSE], w[r] * xi[r])))) / sc
  names(bias_all) <- colnames(X)
  bias <- bias_all[names(b_hat)]
  n <- nrow(X); k <- ncol(X)
  b_til <- b_hat - bias
  se_til <- sqrt(diag(V_hat)) * (n / (n + k))
  p_til <- 2 * pnorm(-abs(b_til / se_til))
  ct <- cbind(b = b_til, se = se_til, p = p_til)
  unc <- cbind(b = b_hat, se = sqrt(diag(V_hat)), p = 2 * pnorm(-abs(b_hat / sqrt(diag(V_hat)))))
  rownames(unc) <- paste0(names(b_hat), "__unc")
  prior_int <- NA
  if (!is.null(tau) && "(Intercept)" %in% names(b_til)) {
    ybar <- mean(model.response(model.frame(f, dk)))
    prior_int <- b_til["(Intercept)"] - log(((1 - tau) / tau) * (ybar / (1 - ybar)))
  }
  out <- rbind(ct, unc)
  attr(out, "n") <- n; attr(out, "k") <- k; attr(out, "prior_intercept") <- prior_int
  rm(m, X, A); invisible(gc())
  out
}
