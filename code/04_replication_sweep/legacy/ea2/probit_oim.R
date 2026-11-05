# Stata-equivalent probit SEs from the OBSERVED Hessian (Stata probit is Newton-Raphson ML; vce(oim) default, robust/cluster sandwich with OIM bread).
# type: "oim" (model-based), "robust" (HC, n/(n-1)), or a cluster column name (G/(G-1)). Point estimates from fixest feglm (identical ML).
suppressPackageStartupMessages(library(fixest))
probit_stata <- function(f, data, terms, type = "oim") {
  vc <- if (type == "oim") "iid" else if (type == "robust") "hetero" else as.formula(paste0("~", type))
  m <- feglm(f, data = data, family = binomial("probit"), vcov = vc, notes = FALSE)
  stopifnot(m$convStatus); used <- obs(m); dd <- data[used, , drop = FALSE]
  X <- model.matrix(f, dd); y <- dd[[all.vars(f)[1]]]; b <- coef(m)[colnames(X)]; stopifnot(!anyNA(b), nrow(X) == nobs(m))
  eta <- drop(X %*% b); q <- 2 * y - 1; lam <- q * dnorm(q * eta) / pnorm(q * eta)
  Hi <- solve(crossprod(X, X * (lam * (lam + eta))))            # inverse of minus the observed Hessian
  V <- if (type == "oim") Hi else if (type == "robust") { n <- nrow(X); Hi %*% crossprod(X * lam) %*% Hi * n / (n - 1) } else {
    S <- rowsum(X * lam, dd[[type]]); G <- nrow(S); Hi %*% crossprod(S) %*% Hi * G / (G - 1) }
  se <- sqrt(diag(V)); ct <- coeftable(m)
  out <- sapply(terms, function(t) c(b = unname(b[t]), se = unname(se[t]), p = unname(2 * pnorm(-abs(b[t] / se[t]))), se_fixest = unname(ct[t, 2]), p_fixest = unname(ct[t, 4])))
  attr(out, "n") <- nobs(m); rm(m); invisible(gc()); out
}
