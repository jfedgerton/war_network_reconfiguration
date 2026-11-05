# P8_008 Nieman 2016 JPR "The return on social bonds: Social hierarchy and international conflict" (10.1177/0022343316642509)
# Data + do + log: PRIO JPR archive (files.prio.org/Journals/JPR/2016/53/5/Mark Nieman.zip).
# Headline = Table 3, two-stage strategic probit by statistical backwards induction (do-file):
#   Stage 1 (US punish, N 652): probit Punished_MID_Sanct_lag relUS_SecHierBA relUS_EconHierBA rus_ally_hierA PowerRatio100A_US
#       PowerRatio100A_US2 lndistance_US Democ_Chal ongoingMID_US Prev_Chal constantXB22, nocons  -> pr_pun
#   Stage 2 (challenge, directed dyads, N 549,570): probit Challenge FUS_SecHierA FUS_EconHierA F*(conflict utility x pr_pun)
#       FconstantA1 FconstantXA21, nocons, vce(bootstrap, reps(500)), seed 83.
# KEY TERMS (a priori, abstract): "a state's degree of subordination reduces its probability of conflict initiation against
#   other subordinates" -> challenger's subordination to US in the challenge (status-quo utility) equation: FUS_SecHierA (security
#   hierarchy, primary) and FUS_EconHierA (economic hierarchy). Both enter as -hierarchy; positive b = subordination lowers challenge.
# History: added as extra regressors in the stage-2 challenge equation only (stage 1 is US punishment of the challenger, N 652, and
#   is held fixed across steps, as the authors' bootstrap also holds pr_pun fixed). Same estimator (nocons probit, bootstrap 500) in all steps.
# Usage: Rscript rerun_P8_008.R prep | <model> | write
args <- commandArgs(TRUE); stage <- args[1]
# port: run_one.R calls the script without arguments -> run all stages in order, each in its own R process as in the sweep
if (length(args) == 0) {
  self <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)); stopifnot(length(self) == 1)
  for (s in c("prep", "m0", "mD", "mDM", "m26", "bC", "bD", "write")) stopifnot(system2(file.path(R.home("bin"), "Rscript"), c(shQuote(self), s)) == 0)
  quit(save = "no")
}
CODE <- Sys.getenv("SWEEP_CODE"); stopifnot(nzchar(CODE))
source(file.path(CODE, "work_p8/common_b2.R"))
suppressPackageStartupMessages(library(haven))
WD <- file.path(ROOT, "work_p8/P8_008"); W <- file.path(OUT, "work_p8/P8_008"); dir.create(W, recursive = TRUE, showWarnings = FALSE); lf <- file.path(W, paste0("log_P8_008_", stage, ".txt")); sink(lf, split = TRUE)
terms <- c("FUS_SecHierA", "FUS_EconHierA")
rhs <- c("FUS_SecHierA","FUS_EconHierA","FPowerRatioA","FPowerRatioA_2","FpowerchangeA","FPrev_Chal","Fcontig","FJointDem","Fally",
         "Fbitrade","FcivwarA","FconstantA1","FconstantXA21")
fit_probit <- function(r, data, start = NULL) {
  f <- as.formula(paste("Challenge ~ 0 +", paste(r, collapse = " + ")))
  m <- feglm(f, data = data, family = binomial("probit"), vcov = "iid", notes = FALSE, nthreads = 1, start = start, lean = TRUE)
  stopifnot(m$convStatus); m
}
if (stage == "prep") {
  keep <- c("ccode1","ccode2","year","Challenge","Punished_MID_Sanct_lag","relUS_SecHierBA","relUS_EconHierBA","rus_ally_hierA",
            "PowerRatio100A_US","PowerRatio100A_US2","lndistance_US","Democ_Chal","ongoingMID_US","Prev_Chal","constantXB22",
            "PowerRatioA","powerchange100","civwarA","lndistance","contig","JointDem","bitrade_gdp","ally","constantXA21",
            "US_SecHierA","US_EconHierA","constantXA1")
  # port: the .dta was unzipped from nieman.zip and later deleted in the sweep -> idempotent unzip into OUT
  zm <- unzip(file.path(WD, "nieman.zip"), list = TRUE)$Name; zm <- zm[basename(zm) == "BargHierGenData_JPR.dta"]; stopifnot(length(zm) == 1)
  if (!file.exists(file.path(W, "BargHierGenData_JPR.dta"))) { unzip(file.path(WD, "nieman.zip"), files = zm, exdir = W, junkpaths = TRUE) }
  d <- zap_labels(read_dta(file.path(W, "BargHierGenData_JPR.dta"), col_select = all_of(keep))) %>% as.data.frame()
  cat("rows", nrow(d), "years", range(d$year), "\n")
  p1 <- glm(Punished_MID_Sanct_lag ~ 0 + relUS_SecHierBA + relUS_EconHierBA + rus_ally_hierA + PowerRatio100A_US + PowerRatio100A_US2 +
            lndistance_US + Democ_Chal + ongoingMID_US + Prev_Chal + constantXB22, family = binomial("probit"), data = d)
  print(summary(p1)$coefficients); cat("stage-1 N", nobs(p1), "(pub 652; relUS_SecHierBA 1.031771)\n")
  stopifnot(nobs(p1) == 652, abs(coef(p1)["relUS_SecHierBA"] - 1.031771) < 1e-4)
  X1 <- model.matrix(~ 0 + relUS_SecHierBA + relUS_EconHierBA + rus_ally_hierA + PowerRatio100A_US + PowerRatio100A_US2 +
            lndistance_US + Democ_Chal + ongoingMID_US + Prev_Chal + constantXB22, model.frame(~ ., d[, c("relUS_SecHierBA","relUS_EconHierBA","rus_ally_hierA","PowerRatio100A_US","PowerRatio100A_US2","lndistance_US","Democ_Chal","ongoingMID_US","Prev_Chal","constantXB22")], na.action = na.pass))
  d$pr_pun <- as.vector(pnorm(X1 %*% coef(p1)))
  d <- d %>% mutate(FPowerRatioA = PowerRatioA * pr_pun, FPowerRatioA_2 = PowerRatioA^2 * pr_pun, FpowerchangeA = powerchange100 * pr_pun,
    FPrev_Chal = Prev_Chal * pr_pun, FcivwarA = civwarA * pr_pun, Fcontig = contig * pr_pun, FJointDem = JointDem * pr_pun,
    Fbitrade = bitrade_gdp * pr_pun, Fally = ally * pr_pun, FconstantXA21 = constantXA21 * (1 - pr_pun),
    FUS_SecHierA = -US_SecHierA, FUS_EconHierA = -US_EconHierA, FconstantA1 = -constantXA1)
  d <- d[, c("ccode1","ccode2","year","Challenge", rhs)]
  m <- fit_probit(rhs, d); print(coeftable(m)); cat("stage-2 full N", nobs(m), "(pub 549,570; FUS_SecHierA .1852821, FUS_EconHierA .0010267)\n")
  full <- list(b = coef(m), n = nobs(m))
  ds <- common_sample(d, c("Challenge", rhs), "ccode1", "ccode2", "year")
  ds <- ds[, c("ccode1","ccode2","year","Challenge", rhs, HV)]
  saveRDS(list(ds = ds, full = full), file.path(W, "prep.rds"))
} else if (stage %in% names(ADD)) {
  P <- readRDS(file.path(W, "prep.rds")); ds <- P$ds
  r <- c(rhs, ADD[[stage]]); m <- fit_probit(r, ds); b <- coef(m); print(coeftable(m))
  # Nonparametric bootstrap of stage 2 (resample obs with replacement) implemented as multinomial frequency weights + weighted
  # Fisher-scoring probit warm-started at b (identical estimator to refitting on the resampled rows; much faster).
  X <- as.matrix(ds[, r]); y <- ds$Challenge
  probit_w <- function(w, b) {
    for (it in 1:50) {
      eta <- as.vector(X %*% b); mu <- pmin(pmax(pnorm(eta), 1e-12), 1 - 1e-12); de <- dnorm(eta)
      ww <- w * de^2 / (mu * (1 - mu)); z <- w * de * (y - mu) / (mu * (1 - mu))
      step <- solve(crossprod(X, X * ww), crossprod(X, z))[, 1]; b <- b + step
      if (max(abs(step) / pmax(1, abs(b))) < 1e-7) return(b)
    }
    NULL
  }
  bchk <- probit_w(rep(1, length(y)), b * 0.9); print(rbind(b, bchk)); stopifnot(!is.null(bchk), max(abs(bchk - b) / pmax(1, abs(b))) < 1e-4)   # weights 1 reproduces the fixest fit
  set.seed(83); B <- 100; n <- nrow(ds); bs <- matrix(NA, B, length(b), dimnames = list(NULL, names(b)))
  for (i in 1:B) {
    w <- as.vector(rmultinom(1, n, rep(1 / n, n)))
    bb <- tryCatch(probit_w(w, b), error = function(e) NULL)
    if (!is.null(bb)) bs[i, ] <- bb
    if (i %% 50 == 0) cat("rep", i, format(Sys.time()), "\n")
  }
  ok <- complete.cases(bs); cat("bootstrap reps converged:", sum(ok), "\n"); stopifnot(sum(ok) >= 90)
  se <- apply(bs[ok, ], 2, sd); p <- 2 * pnorm(-abs(b / se))
  print(round(cbind(b, se, p), 6))
  out <- lapply(setNames(terms, terms), function(t) c(b = unname(b[t]), se = unname(se[t]), p = unname(p[t]))); attr(out, "n") <- nobs(m)
  saveRDS(out, file.path(W, paste0("fit_", stage, ".rds")))
} else if (stage == "write") {
  P <- readRDS(file.path(W, "prep.rds"))
  M <- lapply(setNames(names(ADD), names(ADD)), function(k) readRDS(file.path(W, paste0("fit_", k, ".rds"))))
  stopifnot(length(unique(sapply(M, attr, "n"))) == 1)
  meta <- list(pid = "P8_008", cite = "Nieman 2016", journal = "JPR", year = 2016, doi = "10.1177/0022343316642509", family = "TIE")
  ae <- "two-stage strategic probit (statistical backwards induction): stage-1 probit US punishment -> pr_pun; stage-2 nocons probit of challenge with bootstrap SE (500 reps, seed 83)"
  eu <- "same two-stage estimator in all steps: stage 1 fixed (reproduced exactly), stage-2 nocons probit (fixest) + nonparametric bootstrap SE with 100 reps (authors 500; reduced for runtime under the shared R lock, flagged), seed 83; history added to stage 2"
  nt <- sprintf("Full-sample stage-2 repro N=%d: FUS_SecHierA %.6f (pub .1852821), FUS_EconHierA %.6f (pub .0010267). Directed dyads 1950-2000; published values from authors' log (Table 3). Common sample N=%d. SE by bootstrap on common sample (R RNG differs from Stata's, so bootstrap SEs are not bit-identical to the log).",
    P$full$n, P$full$b["FUS_SecHierA"], P$full$b["FUS_EconHierA"], attr(M$m0, "n"))
  vf <- function(t, pb, sig) verdict_fn(P$full$b[t], if (sig) 0 else 1, pb, sig, FALSE)   # full-sample b; sig from full-sample bootstrap not rerun -> use common-sample p
  write_p8(meta, "FUS_SecHierA", "\"a state's degree of subordination reduces its probability of conflict initiation\"", "Table 3 challenge equation",
    ae, eu, 0.1852821, 0.0560249, 549570, TRUE, FALSE, M, "N (Prev_Chal = prior challenges, weighted by pr_pun)",
    paste(nt, "Coefficient is on -(security hierarchy): positive = subordination lowers initiation."),
    verdict = verdict_fn(M$m0$FUS_SecHierA["b"], M$m0$FUS_SecHierA["p"], 0.1852821, TRUE, FALSE))
  write_p8(meta, "FUS_EconHierA", "\"a state's degree of subordination reduces its probability of conflict initiation\"", "Table 3 challenge equation",
    ae, eu, 0.0010267, 0.0475494, 549570, FALSE, FALSE, M, "N (Prev_Chal = prior challenges, weighted by pr_pun)",
    paste(nt, "Economic hierarchy dimension; published insignificant (tiny b, so % comparisons unstable)."),
    verdict = verdict_fn(M$m0$FUS_EconHierA["b"], M$m0$FUS_EconHierA["p"], 0.0010267, FALSE, FALSE))
}
sink()
