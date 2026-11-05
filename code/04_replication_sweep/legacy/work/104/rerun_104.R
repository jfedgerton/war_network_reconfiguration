# Rerun 104: Karreth 2018 II -- claims main model (SI Table A3 / Figure 1):
#   useforce ~ igo_lev3_count_use_bc + salint + saltan + terriss + jointpol7 + rivalry_th + absidealdiff + atopally
# Author: rstanarm stan_glm, priors N(0,10) on all coefs (diffuse) -> posterior mean ~ ML estimate; author also fits glm (m1).
# We use ML logit with model-based SEs (analogue of posterior SD). Published posterior mean for tie: -0.35 (SD 0.16), 168 claims.
ROOT <- Sys.getenv("SWEEP_DATA"); OUT <- Sys.getenv("SWEEP_OUT"); CODE <- Sys.getenv("SWEEP_CODE")
stopifnot(nzchar(ROOT), nzchar(OUT), nzchar(CODE))
source(file.path(CODE, "work/rerun_common_g7.R"))
W <- file.path(ROOT, "work/104")
WO <- file.path(OUT, "work/104"); dir.create(WO, recursive = TRUE, showWarnings = FALSE)  # port: outputs
sink(file.path(WO, "log_104.txt"), split = TRUE)

# 1. Data prep exactly as 11_claims_readdata.R (only lines relevant to model 1)
dat <- read.csv(file.path(W, "claims_full.csv"))
dat$useforce <- ifelse(dat$maxhost > 3, 1, 0)
dat <- dat[dat$yearbeg > 1945, ]
dat[dat$dyadid == 365705 & dat$year == 1991, ]$absidealdiff <- .472
dat[dat$dyadid == 373701 & dat$year == 1991, ]$absidealdiff <- .029
dat[dat$dyadid == 630701 & dat$year == 1991, ]$absidealdiff <- 1.232
dat[dat$dyadid == 640645 & dat$year == 1999, ]$absidealdiff <- 2.31
BASE <- "useforce ~ igo_lev3_count_use_bc + salint + saltan + terriss + jointpol7 + rivalry_th + absidealdiff + atopally"
MV <- c("useforce","igo_lev3_count_use_bc","salint","saltan","terriss","jointpol7","rivalry_th","absidealdiff","atopally")
fitter_ml <- function(add, dat) fit_glm_tie(as.formula(paste0(BASE, add)), dat, binomial(), "igo_lev3_count_use_bc", cluster = NULL)

# Bayesian logit matching the author's priors (N(0,10) on coefficients and intercept; no autoscaling -- rstanarm may
# autoscale by sd(x), immaterial at this diffuseness). Random-walk Metropolis, proposal = 2.38^2/k * ML vcov,
# 4 chains x (5000 burn + 25000 kept), seed 123. Returns posterior mean, posterior SD, two-sided tail prob 2*min(P>0,P<0).
fitter <- function(add, dat) {
  f <- as.formula(paste0(BASE, add))
  dat <- dat[complete.cases(dat[, all.vars(f)]), ]
  m <- glm(f, data = dat, family = binomial()); stopifnot(m$converged)
  X <- model.matrix(m); y <- m$y; k <- ncol(X)
  lp <- function(b) { eta <- drop(X %*% b); sum(y * eta - log1p(exp(eta))) + sum(dnorm(b, 0, 10, log = TRUE)) }
  L <- t(chol(vcov(m) * 2.38^2 / k))
  set.seed(123)
  draws <- do.call(rbind, lapply(1:4, function(ch) {
    b <- coef(m) + drop(L %*% rnorm(k)); cur <- lp(b); out <- matrix(NA, 25000, k); acc <- 0
    for (i in 1:30000) { prop <- b + drop(L %*% rnorm(k)); lpp <- lp(prop)
      if (log(runif(1)) < lpp - cur) { b <- prop; cur <- lpp; acc <- acc + 1 }
      if (i > 5000) out[i - 5000, ] <- b }
    cat(sprintf("  chain %d acceptance %.2f\n", ch, acc / 30000)); out }))
  colnames(draws) <- colnames(X)
  pm <- colMeans(draws); psd <- apply(draws, 2, sd); ptail <- apply(draws, 2, function(v) 2 * min(mean(v > 0), mean(v < 0)))
  cat("\n---- ", deparse(f, width.cutoff = 500), "| N =", nrow(X), "\n")
  print(signif(cbind(post_mean = pm, post_sd = psd, tail_p = ptail, ml = coef(m)), 4))
  tie <- "igo_lev3_count_use_bc"
  list(b = unname(pm[tie]), se = unname(psd[tie]), p = unname(ptail[tie]), n = nrow(X))
}

# 2. Reproduce
cat("\n==== Reproduction: claims model 1 (ML logit, author's glm m1) ====\n")
rp_ml <- fitter_ml("", dat)
cat("\n==== Reproduction: claims model 1 (Bayesian logit, MH) ====\n")
rp <- fitter("", dat)
b_pub <- -0.35
stopifnot(rp$n == 168)
# significance of the Bayesian estimate: 95% interval -0.35 +/- 1.96*0.16 excludes 0 -> treat as p<.05
repro <- if (abs(rp$b - b_pub) / abs(b_pub) <= 0.10 && rp$p < .05) "REPRODUCED" else "NOT_REPRODUCED"
cat("\nb_pub =", b_pub, " b_repro =", rp$b, " status:", repro, "\n")
n_orig <- rp$n

# 3. History merge (claim-year; ccode1/ccode2 are the claim dyad) and restriction (all post-1945 -> no year drop)
dat <- add_history(dat, "ccode1", "ccode2", "year")
d_all <- dat
dat <- restrict_sample(dat, MV)

# 4. Six models
R <- run_six(fitter, dat)

SENS <- sens_hist_only(fitter, d_all, MV)
cat("\n", SENS, "\n")

# 5. Results
write_row(104, "Karreth 2018", repro, b_pub, rp$b, "Bayesian logit, N(0,10) priors, random-walk Metropolis (4x25k draws); b=posterior mean, se=posterior SD, p=2*min tail prob", "igo_lev3_count_use_bc",
  n_orig, R$m0$n, R, ahc = "Y",
  notes = sprintf("Claims M1 (SI Table A3), 168 ICOW claims post-1945, DV = use of force within claim (conditional sample). b_pub = posterior mean (rounded to 2dp); ML glm on same data gives %.3f (p=%.3f). Author has Thompson rivalry control (Y). Small N=168; history vars are claim-dyad MID counts. Dropped %d obs vs orig (history/benchmark missing, bm_* NaN for new states e.g. post-Soviet 1991, Belize 1981); %s. Sensitivity %s.", rp_ml$b, rp_ml$p, n_orig - R$m0$n, RS_INFO, SENS))
sink()
