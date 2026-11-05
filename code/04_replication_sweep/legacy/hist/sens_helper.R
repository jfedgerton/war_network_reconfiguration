# sens_helper.R -- sensitivity analyses for one study (appendix). Sourced ONLY by the hook at the end of a study script,
# and only when run_sensitivity.R sets SWEEP_SENS=1. The main sweep (run_sweep.R) never runs this code.
#
# For the study's headline tie term it produces, with the article's own estimator:
#   1. SPLIT    the tie coefficient with no history, dyadic history only, monadic history only, and both (b, se, p).
#   2. STRENGTH how strongly conflict history predicts the tie variable and the outcome, net of the study's controls,
#               next to the two benchmark variables (lagged joint democracy, lagged capability ratio). These are the
#               two partial R2 of the omitted-variable-bias formula (Cinelli and Hazlett 2020), computed from linear
#               regressions on the study's estimation sample (a linear approximation for the nonlinear estimators).
#   3. PLACEBO  a matched placebo for each history specification. Each history variable h is split, net of the
#               controls X, into the part predicted by the tie variable D and a remainder u:
#                   h~ = beta * D~ + u            (h~, D~ = residuals of h and D on X)
#               The synthetic variable is  z(lambda) = beta * D~ + lambda * u + sqrt(1 - lambda^2) * u*, where u* is a
#               random permutation of u that is then residualised on D and X (and the absorbed fixed effects) and
#               rescaled to the length of u. Because u and u* are both orthogonal to D given X, z has exactly the same
#               partial association with the tie variable as h for every lambda and every draw. At lambda = 1 it is
#               h itself (given X); at lambda = 0 the remainder is a shuffled copy of u, so its link to the outcome is
#               broken by construction (rows are shuffled across the whole sample, including across fixed-effect groups). If adding history only moved the tie coefficient because history is
#               collinear with the tie variable, the lambda = 0 variable would move it just as much.
#   4. DIAL     (selected studies) lambda = 0.25, 0.5, 0.75 for the both-variables specification.
#
# Arguments of sens_run():
#   id    study id (registry.csv)
#   dat   the estimation sample the study's history models use (after add_history() and the sample restriction)
#   y     outcome column;  tie  the headline tie term's column
#   ctrl  the model's other right-hand-side terms, as R formula terms (character vector), WITHOUT the tie term
#   fe    names of absorbed fixed-effect variables (character vector) or NULL
#   fe_sep  TRUE for fixed-effect logit/probit studies: the estimator drops fixed-effect groups in which the outcome does
#         not vary, so dat is first cut to the rows fixest keeps (obs() of a fixed-effects-only binomial fit)
#   fit   function(dat, add): fits the study's model on `dat` with the columns named in `add` appended to the baseline
#         right-hand side and returns c(b = , se = , p = ) for the tie term
# Output: <SWEEP_SENS_OUT>/sens_<id>.rds (numbers keyed by study id; no citation).
# sens_split(id, dat, fit): the SPLIT only (step 1), for studies whose headline term is an interaction, a product term or
#   a constituent of one. Steps 2-4 need a single tie column whose association with history is defined net of the
#   controls; with an interaction the controls contain the constituent terms, so those steps are not computed.
# Settings (environment): SENS_DRAWS placebo draws per specification (default 5); SENS_DIAL = 1 adds the dial.
suppressPackageStartupMessages(library(fixest))

sens_run <- function(id, dat, y, tie, ctrl, fe = NULL, fit, fe_sep = FALSE) {
  OUTD <- Sys.getenv("SWEEP_SENS_OUT"); stopifnot(nzchar(OUTD), dir.exists(OUTD))
  DRAWS <- as.integer(Sys.getenv("SENS_DRAWS", "5")); DIAL <- Sys.getenv("SENS_DIAL", "0") == "1"
  HV <- c(dyad = "hist_dyad_w15", mono = "hist_mono_w15"); BM <- c(joint_democracy = "bm_democ_w15", capability_ratio = "bm_cinc_w15")
  stopifnot(is.data.frame(dat), length(y) == 1, length(tie) == 1, all(c(y, tie, HV, BM) %in% names(dat)), DRAWS >= 1,
            !anyNA(dat[, c(y, tie, HV, BM)]), is.function(fit))
  dat <- as.data.frame(dat)
  if (fe_sep) { stopifnot(length(fe) > 0, all(dat[[y]] %in% c(0, 1)))
    m <- feglm(as.formula(paste(y, "~ 1 |", paste(fe, collapse = " + "))), data = dat, family = binomial(), notes = FALSE)
    dat <- dat[obs(m), , drop = FALSE]; rownames(dat) <- NULL; rm(m) }
  n <- nrow(dat)
  set.seed(123)
  f3 <- function(add) { r <- fit(dat, add); stopifnot(is.numeric(r), all(c("b", "se", "p") %in% names(r))); r[c("b", "se", "p")] }

  # ---- 1. split: none / dyadic / monadic / both ----
  split <- rbind(none = f3(character(0)), dyad = f3(HV[["dyad"]]), mono = f3(HV[["mono"]]), both = f3(unname(HV)))
  stopifnot(all(is.finite(split[, "b"])))

  # ---- residuals net of the controls (and absorbed fixed effects) ----
  rhs <- function(extra = NULL) { t <- c(extra, ctrl); t <- if (length(t)) paste(t, collapse = " + ") else "1"
    if (length(fe)) paste(t, "|", paste(fe, collapse = " + ")) else t }
  res <- function(v, extra = NULL) {
    m <- feols(as.formula(paste(v, "~", rhs(extra))), data = dat, notes = FALSE, warn = FALSE, fixef.rm = "none")   # keep singleton groups (residual 0)
    if (nobs(m) != n) stop("residualising ", v, " used ", nobs(m), " of ", n, " rows (missing values in the controls?)")
    as.numeric(resid(m))
  }
  Dt <- res(tie)                       # tie | X
  Yt <- res(y, tie)                    # outcome | tie, X   (linear probability residual)
  if (sum(Dt^2) < 1e-10 * sum((dat[[tie]] - mean(dat[[tie]]))^2)) stop("tie variable is (almost) a linear function of the controls")
  # a benchmark that does not vary in this sample has no strength to report (NA), as in the main sweep
  resb <- function(v, extra = NULL) if (length(unique(dat[[v]])) < 2) rep(NA_real_, n) else res(v, extra)
  Hx <- sapply(HV, res); Bx <- sapply(BM, resb)                                 # history and benchmarks | X
  Hy <- sapply(HV, res, extra = tie); By <- sapply(BM, resb, extra = tie)       # ... | tie, X
  r2 <- function(target, V) { V <- as.matrix(V); if (anyNA(V)) return(NA_real_); 1 - sum(lm.fit(cbind(1, V), target)$residuals^2) / sum((target - mean(target))^2) }
  strength <- rbind(
    history_both     = c(r2_tie = r2(Dt, Hx),                r2_out = r2(Yt, Hy),                pcor_tie = NA),
    history_dyad     = c(r2_tie = r2(Dt, Hx[, "dyad"]),      r2_out = r2(Yt, Hy[, "dyad"]),      pcor_tie = cor(Dt, Hx[, "dyad"])),
    history_mono     = c(r2_tie = r2(Dt, Hx[, "mono"]),      r2_out = r2(Yt, Hy[, "mono"]),      pcor_tie = cor(Dt, Hx[, "mono"])),
    joint_democracy  = c(r2_tie = r2(Dt, Bx[, 1]),           r2_out = r2(Yt, By[, 1]),           pcor_tie = cor(Dt, Bx[, 1])),
    capability_ratio = c(r2_tie = r2(Dt, Bx[, 2]),           r2_out = r2(Yt, By[, 2]),           pcor_tie = cor(Dt, Bx[, 2])))

  # ---- 3. matched placebo (and 4. dial) ----
  beta <- colSums(Hx * Dt) / sum(Dt^2); U <- Hx - outer(Dt, beta)   # remainder of each history variable (orthogonal to D~)
  stopifnot(identical(colnames(Hx), names(HV)), max(abs(colSums(U * Dt))) < 1e-6 * sum(Dt^2))
  zfit <- function(vars, lam, perm) {
    z <- dat; nm <- paste0("sens_z_", vars)
    for (k in vars) {
      up <- 0
      if (lam < 1) { dat$sens_tmp_u <<- U[perm, k]; up <- res("sens_tmp_u", tie)      # shuffled remainder, made orthogonal to D and X again
        stopifnot(sum(up^2) > 0, abs(sum(up * Dt)) < 1e-6 * sqrt(sum(up^2) * sum(Dt^2))); up <- up * sqrt(sum(U[, k]^2) / sum(up^2)) }
      z[[paste0("sens_z_", k)]] <- beta[[k]] * Dt + lam * U[, k] + sqrt(1 - lam^2) * up
    }
    r <- fit(z, nm); r[["b"]]
  }
  SPECS <- list(dyad = "dyad", mono = "mono", both = c("dyad", "mono"))
  # lambda = 1 must reproduce the study's own history models (same column space given the controls)
  # (difference scaled by the larger of |b| and its standard error, so a history coefficient near zero does not break the check)
  check1 <- sapply(names(SPECS), function(s) (zfit(SPECS[[s]], 1, seq_len(n)) - split[s, "b"]) / max(abs(split[s, "b"]), split[s, "se"]))
  # one seed per fit, drawn up front: a study's fitting function may reset the random-number stream itself
  GRID <- rbind(data.frame(spec = names(SPECS), lambda = 0), if (DIAL) data.frame(spec = "both", lambda = c(.25, .5, .75)))
  GRID <- GRID[rep(seq_len(nrow(GRID)), each = DRAWS), ]; GRID$draw <- rep(seq_len(DRAWS), length.out = nrow(GRID))
  GRID$seed <- sample.int(.Machine$integer.max, nrow(GRID)); stopifnot(!anyDuplicated(GRID$seed))
  GRID$b <- sapply(seq_len(nrow(GRID)), function(i) { set.seed(GRID$seed[i]); zfit(SPECS[[GRID$spec[i]]], GRID$lambda[i], sample.int(n)) })
  placebo <- GRID[GRID$lambda == 0, c("spec", "lambda", "draw", "b")]; rownames(placebo) <- NULL

  out <- list(id = id, n = n, tie = tie, split = split, strength = strength, check_lambda1 = check1,
              check_ok = all(abs(check1) < 1e-4), placebo = `rownames<-`(GRID[, c("spec", "lambda", "draw", "b")], NULL), draws = DRAWS, dial = DIAL,
              split_only = FALSE)
  saveRDS(out, file.path(OUTD, paste0("sens_", id, ".rds")))
  cat("\n==== sensitivity", id, ": n", n, " b none/dyad/mono/both", signif(split[, "b"], 5), "\n     lambda=1 check (scaled difference)", signif(check1, 3),
      if (out$check_ok) "OK" else "FAILED", "\n     placebo/none (mean over draws):",
      paste(sapply(names(SPECS), function(s) sprintf("%s %.3f", s, mean(placebo$b[placebo$spec == s]) / split["none", "b"])), collapse = "  "), "\n")
  invisible(out)
}

sens_split <- function(id, dat, fit) {
  OUTD <- Sys.getenv("SWEEP_SENS_OUT"); stopifnot(nzchar(OUTD), dir.exists(OUTD))
  HV <- c(dyad = "hist_dyad_w15", mono = "hist_mono_w15")
  stopifnot(is.data.frame(dat), all(HV %in% names(dat)), !anyNA(dat[, HV]), is.function(fit))
  dat <- as.data.frame(dat)
  f3 <- function(add) { r <- fit(dat, add); stopifnot(is.numeric(r), all(c("b", "se", "p") %in% names(r))); r[c("b", "se", "p")] }
  split <- rbind(none = f3(character(0)), dyad = f3(HV[["dyad"]]), mono = f3(HV[["mono"]]), both = f3(unname(HV)))
  stopifnot(all(is.finite(split)))
  out <- list(id = id, n = nrow(dat), split = split, split_only = TRUE)
  saveRDS(out, file.path(OUTD, paste0("sens_", id, ".rds")))
  cat("\n==== sensitivity (split only)", id, ": n", nrow(dat), " b none/dyad/mono/both", signif(split[, "b"], 5), "\n")
  invisible(out)
}
