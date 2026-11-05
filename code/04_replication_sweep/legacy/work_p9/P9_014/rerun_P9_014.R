# P9_014 Findley & Teo 2006 JOP "Rethinking Third-Party Interventions into Civil Wars: An Actor-Centric Approach"
# (10.1111/j.1468-2508.2006.00473.x). Data + do-file: Harvard Dataverse doi:10.7910/DVN/SDTDQG (replication-data.zip:
# replicate.dta, wars-regression-replicate.do). No author log; article closed access -> published values pending.
# HEADLINE MODEL (do-file, both cureregr calls; Stata 8.2 cureregr = SSC cureregr8, Buxton):
#   stset time1, failure(<side>==1) origin(time time0) time0(time0) exit(time .) id(id)
#   cureregr thirdmajpow colhist coldwar sameregion contigdum logcapratio rival allydum ethnic ideology jointdem,
#     dist(lognormal) class(mixture) link(logistic)
#     scale(rival allydum allyongovtdum allyonoppdum rivalongovtdum rivalonoppdum logfatalities loglforced)       [OIM SE]
#   side = govtinter (intervention on the government side), then oppinter (opposition side).
#   Unit = civil-war phase record x potential intervener (id = 1000*conflict + intervener ccode), 1944-98, 49,395 records.
# Likelihood (transcribed from cureregr8's evaluator pcm000201_lf.ado, staged with the data; key in pcm_key.pdf):
#   pi = invlogit(x_cure b) = cure (never-intervene) fraction; lambda = exp(x_scale b); gamma = exp(shape _cons);
#   tt = (lambda t)^gamma; S(t) = pi + (1 - pi)(1 - Phi(ln tt)); per record (delayed entry at _t0):
#   event: ln(1 - pi) + ln(phi(ln tt) gamma / t) - ln S(t0);  censored: ln S(t) - ln S(t0)  (S(t0) = 1 when t0 = 0).
#   Positive scale coefficient = larger lambda = earlier intervention. ML by nlminb with the analytic score, polished by Newton
#   steps; OIM SE from the numerical Jacobian of the analytic score (Stata ml lf uses numerical derivatives). The analytic score is
#   checked against numerical derivatives of the transcribed likelihood. cureregr fits a constant-only model first and continues
#   from it; mirrored here in every step.
# VERIFICATION of the likelihood machinery (stage 0, fail-loud): the same generic evaluator code reproduces the cureregr package's
#   published example fit (cureregr_predict_eg.pdf: non-mixture Weibull, lml link, logLik -138.50247, all 8 coefficients/SEs).
# stset translation: origin(time time0) with id() = one origin per subject = the subject's first-record time0 (records sorted,
#   contiguous); _t0 = time0 - origin, _t = time1 - origin; failure per record (<= 1 failure per id, at its last record).
# KEY TERMS (a priori; tie terms enter only the scale equation): the abstract's actor-centric claim of "convergent and divergent
#   interests" between potential and actual interveners ->
#   headline gov_allyongovtdum = allyongovtdum, scale eq., government-side model (ally of a state already intervening for the
#   government); non-headline opp_rivalongovtdum = rivalongovtdum, scale eq., opposition-side model (rival of a state already
#   intervening for the government).
# History (precedent _slow/101.R -> work/101/rerun_est_101.R: custom multi-equation estimator, history added to the equation that
#   carries the key term): history and benchmark terms enter the scale equation; cure equation unchanged. Dyad = (civil-war state
#   ccode, intervener ccode = id %% 1000), keyed to the record year.
# No randomness (seed 123 set in common_b3.R).
CODE <- Sys.getenv("SWEEP_CODE"); stopifnot(nzchar(CODE))
source(file.path(CODE, "work_p8/common_b3.R"))
suppressPackageStartupMessages(library(numDeriv))
WD <- file.path(ROOT, "work_p9/P9_014"); W <- file.path(OUT, "work_p9/P9_014"); dir.create(W, recursive = TRUE, showWarnings = FALSE)
sink(file.path(W, "log_P9_014.txt"), split = TRUE)

# ---- likelihood: per-record lnf of cureregr8 evaluators (mixture lognormal logistic = pcm000201; non-mixture Weibull lml = pcm010102) ----
cure_ll <- function(par, X1, X2, X3, t, t0, d, class, dist, link) {
  k1 <- ncol(X1); k2 <- ncol(X2)
  th1 <- drop(X1 %*% par[1:k1]); th2 <- drop(X2 %*% par[k1 + 1:k2]); th3 <- drop(X3 %*% par[-(1:(k1 + k2))])
  lam <- exp(th2); gam <- exp(th3)
  tt <- (lam * t)^gam; tt0 <- (lam * t0)^gam
  pi <- if (link == "logistic") plogis(th1) else exp(-exp(th1))
  if (dist == "lognormal") { k <- pnorm(log(tt)); k0 <- pnorm(log(tt0)); dk <- dnorm(log(tt)) * gam / t }
  else { k <- 1 - exp(-tt); k0 <- 1 - exp(-tt0); dk <- exp(-tt) * gam * tt / t }
  k0[tt0 == 0] <- 0
  if (class == "mixture") {
    S <- 1 + (pi - 1) * k; S0 <- 1 + (pi - 1) * k0
    ifelse(d == 1, log(1 - pi) + log(dk) - log(S0), log(S) - log(S0))
  } else {
    lp <- log(pi)
    ifelse(d == 1, log(-lp) + k * lp + log(dk) - k0 * lp, k * lp - k0 * lp)
  }
}
# analytic score of the mixture / lognormal / logistic lnf (z = ln tt = gamma (ln lambda + ln t))
mix_score <- function(par, X1, X2, X3, t, t0, d) {
  k1 <- ncol(X1); k2 <- ncol(X2)
  th1 <- drop(X1 %*% par[1:k1]); th2 <- drop(X2 %*% par[k1 + 1:k2]); th3 <- drop(X3 %*% par[-(1:(k1 + k2))])
  pi <- plogis(th1); gam <- exp(th3)
  z <- gam * (th2 + log(t)); pos0 <- t0 > 0; z0 <- ifelse(pos0, gam * (th2 + log(ifelse(pos0, t0, 1))), 0)
  Phi <- pnorm(z); phi <- dnorm(z); Phi0 <- ifelse(pos0, pnorm(z0), 0); phi0 <- ifelse(pos0, dnorm(z0), 0); zphi0 <- ifelse(pos0, z0 * phi0, 0)
  S <- 1 - (1 - pi) * Phi; S0 <- 1 - (1 - pi) * Phi0; ev <- d == 1
  g1 <- ifelse(ev, -pi, pi * (1 - pi) * Phi / S) - pi * (1 - pi) * Phi0 / S0
  g2 <- ifelse(ev, -z * gam, -(1 - pi) * gam * phi / S) + (1 - pi) * gam * phi0 / S0
  g3 <- ifelse(ev, 1 - z * z, -(1 - pi) * z * phi / S) + (1 - pi) * zphi0 / S0
  c(crossprod(X1, g1), crossprod(X2, g2), crossprod(X3, g3))
}
newton_polish <- function(p, g, H_of) {
  for (i in 1:30) { st <- solve(H_of(p), g(p)); p <- p - st; if (max(abs(st)) < 1e-10) break }
  stopifnot(i < 30); p
}

# ---- stage 0: verify the evaluator machinery on the cureregr package example (published output, cureregr_predict_eg.pdf p.2) ----
e <- as.data.frame(haven::zap_labels(haven::read_dta(file.path(WD, "example_stata.dta")))); stopifnot(nrow(e) == 250)
E1 <- cbind(age = e$age, l500 = e$l500, m23 = e$m23, `_cons` = 1); E2 <- cbind(m23 = e$m23, `_cons` = 1)
fe <- function(p) -sum(cure_ll(p, E1, E2, E2, e$efsyrs, rep(0, nrow(e)), as.numeric(e$exit == 1), "nonmixture", "weibull", "lml"))
pe <- nlminb(c(0, 0, 0, -1, 0, 0, 0, 0), fe, control = list(eval.max = 5000, iter.max = 2000))$par
pe <- newton_polish(pe, function(p) grad(fe, p), function(p) hessian(fe, p)); see <- sqrt(diag(solve(hessian(fe, pe))))
pub_e <- c(.0676627, 1.35413, .5573143, -2.216853, -.2831445, .605446, .2480482, .0922628)
pub_se <- c(.0246766, .2706938, .2384032, .3448536, .2133382, .1500396, .1504704, .0929942)
cat("STAGE 0 cureregr example: logLik", -fe(pe), "(pub -138.50247); max |b - pub|", max(abs(pe - pub_e)), "; max |se - pub|", max(abs(see - pub_se)), "\n")
stopifnot(abs(-fe(pe) - (-138.50247)) < 1e-4, max(abs(pe - pub_e)) < 1e-5, max(abs(see - pub_se)) < 1e-5)

# ---- data and stset ----
unzip(file.path(WD, "ft_replication.zip"), files = "replication-data/replicate.dta", exdir = W, overwrite = TRUE)
d <- as.data.frame(haven::zap_labels(haven::read_dta(file.path(W, "replication-data/replicate.dta"))))
stopifnot(nrow(d) == 49395, length(unique(d$id)) == 20820)
d <- d[order(d$id, d$time0), ]
stopifnot(all(d$time1 > d$time0), all(unlist(tapply(seq_len(nrow(d)), d$id, function(i) d$time0[i][-1] == d$time1[i][-length(i)]))))  # contiguous
org <- ave(d$time0, d$id, FUN = min); d$t0 <- d$time0 - org; d$t <- d$time1 - org
stopifnot(all(tapply(d$govtinter, d$id, sum) <= 1), all(tapply(d$oppinter, d$id, sum) <= 1))
d$intervener <- d$id %% 1000
stopifnot(all(d$intervener != d$ccode), mean(d$intervener %in% .H$nmc$ccode) == 1)
cat("records", nrow(d), "subjects", length(unique(d$id)), "years", range(d$year), "gov failures", sum(d$govtinter), "opp failures", sum(d$oppinter), "\n")
V1 <- c("thirdmajpow", "colhist", "coldwar", "sameregion", "contigdum", "logcapratio", "rival", "allydum", "ethnic", "ideology", "jointdem")
rhs0 <- c("rival", "allydum", "allyongovtdum", "allyonoppdum", "rivalongovtdum", "rivalonoppdum", "logfatalities", "loglforced")

fit_cure_mix <- function(side, term, tie) function(rhs, data) {
  dd <- as.numeric(data[[side]] == 1)
  X1 <- cbind(as.matrix(data[, V1]), `_cons` = 1); X2 <- cbind(as.matrix(data[, rhs, drop = FALSE]), `_cons` = 1)
  X3 <- matrix(1, nrow(data), 1, dimnames = list(NULL, "_cons"))
  run <- function(A1, A2, start) {
    f <- function(p) -sum(cure_ll(p, A1, A2, X3, data$t, data$t0, dd, "mixture", "lognormal", "logistic"))
    g <- function(p) -mix_score(p, A1, A2, X3, data$t, data$t0, dd)
    p <- nlminb(start, f, g, control = list(eval.max = 5000, iter.max = 3000, rel.tol = 1e-14))$par
    p <- newton_polish(p, g, function(p) jacobian(g, p))
    H <- jacobian(g, p); H <- (H + t(H)) / 2
    stopifnot(max(abs(g(p))) < 1e-6, all(eigen(H, symmetric = TRUE, only.values = TRUE)$values > 0))
    list(b = p, se = sqrt(diag(solve(H))), ll = -f(p), f = f, g = g)
  }
  r0 <- run(X1[, "_cons", drop = FALSE], X2[, "_cons", drop = FALSE], c(5, -5, 0))           # constant-only (as cureregr)
  r <- run(X1, X2, c(rep(0, ncol(X1) - 1), r0$b[1], rep(0, ncol(X2) - 1), r0$b[2], r0$b[3]))  # full model continues from it
  gchk <- max(abs(grad(r$f, r$b) - r$g(r$b))); stopifnot(gchk < 1e-4)                         # analytic score vs numerical
  nm <- c(paste0("cure:", colnames(X1)), paste0("scale:", colnames(X2)), "shape:_cons")
  cat("\n---- cureregr mixture lognormal logistic | failure:", side, "| scale eq:", paste(rhs, collapse = " + "), "| records =", nrow(data),
      "| subjects =", length(unique(data$id)), "| failures =", sum(dd), "| logLik =", round(r$ll, 4), "(const-only", round(r0$ll, 4), ") | score check", signif(gchk, 2), "\n")
  print(data.frame(term = nm, b = signif(r$b, 7), se = signif(r$se, 7), z = round(r$b / r$se, 3)), row.names = FALSE)
  j <- which(nm == paste0("scale:", tie)); stopifnot(length(j) == 1)
  out <- setNames(list(c(b = r$b[j], se = r$se[j], p = 2 * pnorm(-abs(r$b[j] / r$se[j])))), term)
  attr(out, "n") <- nrow(data); out
}
SPEC <- list(gov_allyongovtdum = list(side = "govtinter", tie = "allyongovtdum"),
             opp_rivalongovtdum = list(side = "oppinter", tie = "rivalongovtdum"))
vars <- c("t", "t0", "govtinter", "oppinter", V1, rhs0, "id", "ccode", "intervener", "year")
full <- d[complete.cases(d[, vars]), ]; stopifnot(nrow(full) == 49395)
R <- list()
for (term in names(SPEC)) {
  fitter <- fit_cure_mix(SPEC[[term]]$side, term, SPEC[[term]]$tie)
  cat("\n==== term", term, "(", SPEC[[term]]$side, "model, scale eq.) ====\n")
  r_full <- fitter(rhs0, full); n_orig <- attr(r_full, "n")
  cat("FULL-SAMPLE FIT:", term, "b", r_full[[term]][["b"]], "se", r_full[[term]][["se"]], "N records", n_orig, "(published value pending: no log, closed PDF)\n")
  ds <- common_sample(d, vars, "ccode", "intervener", "year")
  M <- run_six(rhs0, ds, fitter)
  R[[term]] <- list(M = M, r_full = r_full, n_orig = n_orig, n_used = nrow(ds))
}
saveRDS(R, file.path(W, "M_P9_014.rds"))
# ---- result-CSV write step (write_p8 writes to results_p8; the file is moved to results_p9 unchanged) ----
meta <- list(pid = "P9_014", cite = "Findley & Teo 2006", journal = "JOP", year = 2006, doi = "10.1111/j.1468-2508.2006.00473.x", family = "TIE")
ev <- "abstract: 'emphasizing convergent and divergent interests and the connections between potential interveners, actual interveners'"
art <- "split-population (mixture cure) lognormal duration model, logistic cure link (cureregr), tie terms in the scale (timing) equation; OIM SE; separate models for government- and opposition-side intervention"
est <- "custom ML of cureregr8's mixture-lognormal-logistic likelihood (delayed entry), analytic score + Newton, OIM SE; history/benchmarks added to the scale equation, all six steps"
ahc <- "N"
mv <- function(term) {
  src <- file.path(OUT, "results_p8", paste0("P9_014__", term, ".csv")); dst <- file.path(OUT, "results_p9", paste0("P9_014__", term, ".csv"))
  dir.create(dirname(dst), recursive = TRUE, showWarnings = FALSE); stopifnot(file.exists(src), file.rename(src, dst))
}
nt <- function(term) sprintf(paste0("Scale-equation coefficient (%s), %s model; positive = earlier intervention. Published value pending (no author log; article closed access). ",
  "Likelihood transcribed from cureregr8 pcm000201_lf.ado; machinery verified on the package's published example. History/benchmarks in the scale equation; dyad = civil-war state x intervener (id %% 1000), record year. ",
  "Full-sample b=%.4f se=%.4f N=%d records. Common 1822-2011 sample N=%d."),
  SPEC[[term]]$tie, SPEC[[term]]$side, R[[term]]$r_full[[term]][["b"]], R[[term]]$r_full[[term]][["se"]], R[[term]]$n_orig, R[[term]]$n_used)
write_p8(meta, "gov_allyongovtdum", ev, "cureregr #1 (government-side intervention), scale equation: allyongovtdum", art, est,
  NA, NA, NA, NA, FALSE, R$gov_allyongovtdum$M, ahc, nt("gov_allyongovtdum")); mv("gov_allyongovtdum")
write_p8(meta, "opp_rivalongovtdum", ev, "cureregr #2 (opposition-side intervention), scale equation: rivalongovtdum", art, est,
  NA, NA, NA, NA, FALSE, R$opp_rivalongovtdum$M, ahc, nt("opp_rivalongovtdum")); mv("opp_rivalongovtdum")
sink()
