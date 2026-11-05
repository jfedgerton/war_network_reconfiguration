# Recursive bivariate probit (Stata biprobit) by ML, with analytic scores and Stata-style robust SE.
# Eq1: y1* = X1 b1 + e1 ; Eq2: y2* = X2 b2 + e2 (X2 may contain y1) ; corr(e1,e2) = rho = tanh(athrho)
# Bivariate normal CDF: vectorized port of Genz's bvnu (Drezner-Wesolowsky / Genz 2004), per-call scalar r.

bvn_lower <- function(h, k, r) {      # P(X < h, Y < k), corr r (scalar); h,k vectors
  stopifnot(length(r) == 1, abs(r) < 1)
  bvnu(-h, -k, r)
}
bvnu <- function(dh, dk, r) {         # P(X > dh, Y > dk)
  if (abs(r) < 0.3) {
    w <- c(0.1713244923791705, 0.3607615730481384, 0.4679139345726904)
    x <- c(0.9324695142031522, 0.6612093864662647, 0.2386191860831970)
  } else if (abs(r) < 0.75) {
    w <- c(.04717533638651177, 0.1069393259953183, 0.1600783285433464, 0.2031674267230659, 0.2334925365383547, 0.2491470458134029)
    x <- c(0.9815606342467191, 0.9041172563704750, 0.7699026741943050, 0.5873179542866171, 0.3678314989981802, 0.1252334085114692)
  } else {
    w <- c(.01761400713915212, .04060142980038694, .06267204833410906, .08327674157670475, 0.1019301198172404,
           0.1181945319615184, 0.1316886384491766, 0.1420961093183821, 0.1491729864726037, 0.1527533871307259)
    x <- c(0.9931285991850949, 0.9639719272779138, 0.9122344282513259, 0.8391169718222188, 0.7463319064601508,
           0.6360536807265150, 0.5108670019508271, 0.3737060887154196, 0.2277858511416451, 0.07652652113349733)
  }
  w <- c(w, w); x <- c(1 - x, 1 + x)
  h <- dh; k <- dk; hk <- h * k
  if (abs(r) < 0.925) {
    hs <- (h * h + k * k) / 2; asr <- asin(r) / 2
    sn <- sin(asr * x)                                   # length L
    E <- exp((outer(hk, sn) - hs) / matrix(1 - sn^2, length(h), length(sn), byrow = TRUE))
    bvn <- as.vector(E %*% w) * asr / (2 * pi) + pnorm(-h) * pnorm(-k)
  } else {
    if (r < 0) { k <- -k; hk <- -hk }
    as <- 1 - r^2; a <- sqrt(as); bs <- (h - k)^2
    asr <- -(bs / as + hk) / 2; c <- (4 - hk) / 8; d <- (12 - hk) / 80
    bvn <- ifelse(asr > -100, a * exp(asr) * (1 - c * (bs - as) * (1 - d * bs) / 3 + c * d * as^2), 0)
    b <- sqrt(bs); sp <- sqrt(2 * pi) * pnorm(-b / a)
    bvn <- ifelse(hk > -100, bvn - exp(-hk / 2) * sp * b * (1 - c * bs * (1 - d * bs) / 3), bvn)
    a2 <- a / 2; xs <- (a2 * x)^2                         # length L
    ASR <- -(outer(bs, 1 / xs) + hk) / 2                  # n x L
    SP <- 1 + outer(c, xs) * (1 + 5 * outer(d, xs))
    rs <- sqrt(1 - xs)
    EP <- exp(-outer(hk / 2, xs / (1 + rs)^2)) / matrix(rs, length(h), length(xs), byrow = TRUE)
    T <- ifelse(ASR > -100, exp(ASR) * (SP - EP), 0)
    bvn <- (a2 * as.vector(T %*% w) - bvn) / (2 * pi)
    if (r > 0) bvn <- bvn + pnorm(-pmax(h, k))
    else {
      L <- ifelse(h < 0, pnorm(k) - pnorm(h), pnorm(-h) - pnorm(-k))
      bvn <- ifelse(h >= k, -bvn, L - bvn)
    }
  }
  pmax(0, pmin(1, bvn))
}

# per-observation log-likelihood and analytic scores
bp_parts <- function(theta, X1, X2, y1, y2) {
  k1 <- ncol(X1); k2 <- ncol(X2)
  b1 <- theta[1:k1]; b2 <- theta[k1 + 1:k2]; ath <- theta[k1 + k2 + 1]; rho <- tanh(ath); rho <- sign(rho) * min(abs(rho), 1 - 1e-9)
  q1 <- 2 * y1 - 1; q2 <- 2 * y2 - 1
  w1 <- q1 * as.vector(X1 %*% b1); w2 <- q2 * as.vector(X2 %*% b2); rr <- q1 * q2 * rho
  P <- numeric(length(y1))
  same <- (q1 * q2) == 1                               # rr = +rho ; otherwise rr = -rho
  if (any(same))  P[same]  <- bvn_lower(w1[same],  w2[same],  rho)
  if (any(!same)) P[!same] <- bvn_lower(w1[!same], w2[!same], -rho)
  P <- pmax(P, 1e-300)
  s <- sqrt(1 - rr^2)
  g1 <- dnorm(w1) * pnorm((w2 - rr * w1) / s) / P
  g2 <- dnorm(w2) * pnorm((w1 - rr * w2) / s) / P
  phi2 <- exp(-(w1^2 - 2 * rr * w1 * w2 + w2^2) / (2 * (1 - rr^2))) / (2 * pi * s)
  gr <- phi2 / P * q1 * q2 * (1 - rho^2)                 # d/d athrho
  S <- cbind(X1 * (g1 * q1), X2 * (g2 * q2), gr)
  list(ll = log(P), S = S)
}

fit_biprobit <- function(X1, X2, y1, y2) {
  # scale non-constant columns by their sd for numerical conditioning (splines ~ 1e4); map back at the end
  sc1 <- apply(X1, 2, sd); sc1[sc1 == 0] <- 1; sc2 <- apply(X2, 2, sd); sc2[sc2 == 0] <- 1
  f <- fit_biprobit_core(sweep(X1, 2, sc1, "/"), sweep(X2, 2, sc2, "/"), y1, y2)
  D <- c(1 / sc1, 1 / sc2, 1)
  f$coef <- f$coef * D; f$V <- f$V * outer(D, D); f$Vml <- f$Vml * outer(D, D)
  nm <- c(paste0("eq1:", colnames(X1)), paste0("eq2:", colnames(X2)), "athrho")
  names(f$coef) <- nm; dimnames(f$V) <- dimnames(f$Vml) <- list(nm, nm)
  f
}

fit_biprobit_core <- function(X1, X2, y1, y2) {
  k <- ncol(X1) + ncol(X2) + 1
  s1 <- coef(glm(y1 ~ X1 - 1, family = binomial("probit")))
  s2 <- coef(glm(y2 ~ X2 - 1, family = binomial("probit")))
  start <- c(s1, s2, 0)
  cache <- new.env(); cache$th <- NULL
  getp <- function(th) { if (!identical(th, cache$th)) { cache$th <- th; cache$p <- bp_parts(th, X1, X2, y1, y2) }; cache$p }
  nll <- function(th) -sum(getp(th)$ll)
  ngr <- function(th) -colSums(getp(th)$S)
  o <- nlminb(start, nll, ngr, control = list(eval.max = 5000, iter.max = 2000))
  cat("nlminb:", o$message, " obj", o$objective, " iters", o$iterations, "\n")
  th <- o$par
  gfun <- function(t) colSums(getp(t)$S)
  hfun <- function(t) { H <- numDeriv::jacobian(function(u) colSums(bp_parts(u, X1, X2, y1, y2)$S), t, method.args = list(r = 2)); (H + t(H)) / 2 }
  # polish: BHHH steps, then Newton (numeric Hessian of analytic score) if needed
  for (it in 1:30) {
    P <- getp(th); g <- colSums(P$S); if (max(abs(g)) < 1e-7) break
    step <- solve(crossprod(P$S), g); lam <- 1; ll0 <- sum(P$ll)
    while (sum(getp(th + lam * step)$ll) < ll0 - 1e-10 && lam > 1e-4) lam <- lam / 2
    th <- th + lam * step
  }
  for (it2 in 1:10) {
    g <- gfun(th); if (max(abs(g)) < 1e-7) break
    step <- solve(-hfun(th), g); lam <- 1; ll0 <- -nll(th)
    while (-nll(th + lam * step) < ll0 - 1e-10 && lam > 1e-4) lam <- lam / 2
    th <- th + lam * step
  }
  cat("polish BHHH", it, "Newton", it2, " max|score|", max(abs(gfun(th))), "\n")
  stopifnot(max(abs(gfun(th))) < 1e-4)
  H <- hfun(th)
  Ainv <- solve(-H)
  S <- bp_parts(th, X1, X2, y1, y2)$S
  n <- length(y1)
  Vrob <- Ainv %*% crossprod(S) %*% Ainv * n / (n - 1)   # Stata robust
  list(coef = th, V = Vrob, Vml = Ainv, ll = -nll(th), n = n, grad_max = max(abs(colSums(S))), rho = tanh(th[k]))
}
