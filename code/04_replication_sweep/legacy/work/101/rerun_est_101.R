# ESTAUDIT 101: Xiang 2016 II. Article headline = Table 1 Model 3: Bayesian split-population (Boolean) probit with dyad-specific
#   random effects in the relevance equation, estimated with the authors' own Gibbs sampler (dyadEffect package, replication.R /
#   regression.R: set.seed(0); burnin 100,000; mcmc 1,000,000; thin 500; relevance eq. G = smldep smigoabi contigkb logdstab majpower
#   (+ dyad effects alpha ~ N(mu, s2)); conflict eq. X = smldmat smldep smigoabi lcaprat2 allies py1 py2 py3 + constant; authors' start values).
# Previous reruns (rerun_101.R) used a pooled probit (Table 1 Model 1) in all six steps -> estimator mismatch.
# Here: the authors' sampler (src copied to est/xiang_mcmc.cpp; only includes for progress/omp removed), data arrays rebuilt from
#   kantdata.dta (validated cell-by-cell against the authors' regression_data/*.txt), history terms added to the CONFLICT equation
#   (history merged on outcome year = year+1, as in rerun_101.R), common sample = rerun_101's 114,896 dyad-years (other cells set to
#   missing, -999, exactly as the sampler treats missing Y).
# Usage: Rscript rerun_est_101.R time            -> short timing runs (m0, mDM)
#        Rscript rerun_est_101.R run <model> <mcmc> <burnin> <thin>
#        Rscript rerun_est_101.R seg <model> <segno> <iters> <thin>  -> one chain segment (continues from segment segno-1's final state)
#        Rscript rerun_est_101.R write <burnin_iters> -> collect saved segments, drop burn-in, write results_estaudit row
args <- commandArgs(TRUE)
ROOT <- Sys.getenv("SWEEP_DATA"); OUT <- Sys.getenv("SWEEP_OUT"); CODE <- Sys.getenv("SWEEP_CODE")
stopifnot(nzchar(ROOT), nzchar(OUT), nzchar(CODE))
W <- file.path(ROOT, "work/101"); P <- file.path(W, "dyadEffect/data")
E <- file.path(OUT, "work/101/est"); dir.create(E, recursive = TRUE, showWarnings = FALSE)   # port: logs, segments, chains, Rcpp cache
EC <- file.path(CODE, "work/101/est")                                                       # port: sampler source (code)
suppressPackageStartupMessages({library(haven); library(Rcpp)})
logf <- file.path(E, paste0("log_est_101_", paste(args, collapse = "_"), ".txt")); sink(logf, split = TRUE)
source(file.path(CODE, "hist/add_history.R"))

# ---- data arrays ----
d <- as.data.frame(zap_labels(read_dta(file.path(P, "kantdata.dta")))); stopifnot(nrow(d) == 116013)
dy <- unique(d[, c("statea","stateb")]); dy <- dy[order(dy$statea, dy$stateb), ]; stopifnot(nrow(dy) == 5116)
i <- match(paste(d$statea, d$stateb), paste(dy$statea, dy$stateb)); j <- d$year - 1949 + 1; stopifnot(all(j >= 1 & j <= 43), !anyNA(i))
mk <- function(v) { M <- matrix(-999, 5116, 43); M[cbind(i, j)] <- v; M }
Yfull <- mk(d$dispute1)
Aut <- as.matrix(read.table(file.path(P, "regression_data/dispute1.txt"))); stopifnot(all(Aut == Yfull))
known <- Yfull >= 0; stopifnot(sum(known) == 116013)
GV <- c("smldep","smigoabi","contigkb","logdstab","majpower")
XV <- c("smldmat","smldep","smigoabi","lcaprat2","allies","py1","py2","py3")
if (args[1] %in% c("time", "run") || (args[1] == "seg" && args[3] == "1")) for (v in unique(c(GV, XV))) {
  A <- as.matrix(read.table(file.path(P, paste0("regression_data/", v, ".txt")))); B <- mk(d[[v]])
  dev <- max(abs(A[known] - B[known]) / pmax(1, abs(B[known])))
  cat(sprintf("validate %-9s max rel dev vs authors' txt on known cells: %.2e\n", v, dev)); stopifnot(dev < 1e-4)
}
# history (merged on outcome year t+1, as rerun_101.R) and common sample
d$year_out <- d$year + 1
d <- add_history(d, "statea", "stateb", "year_out")
HVv <- c("hist_dyad_w15","hist_mono_w15","hist_dyad_w26","hist_mono_w26","bm_cinc_w15","bm_democ_w15")
cs <- d$year_out >= 1822 & d$year_out <= 2011 & complete.cases(d[, HVv]) & is.finite(d$bm_cinc_w15)
cat("common sample:", sum(cs), " (rerun_101: 114896)\n"); stopifnot(sum(cs) == 114896)
Ycs <- matrix(-999, 5116, 43); Ycs[cbind(i[cs], j[cs])] <- d$dispute1[cs]
ADD <- list(m0 = character(0), mD = "hist_dyad_w15", mDM = c("hist_dyad_w15","hist_mono_w15"), m26 = c("hist_dyad_w26","hist_mono_w26"),
            bC = "bm_cinc_w15", bD = "bm_democ_w15")
arr <- function(vars) { A <- array(0, c(5116, 43, length(vars))); for (k in seq_along(vars)) { M <- mk(d[[vars[k]]]); M[is.na(M)] <- 0; A[, , k] <- M }; A }
beta_start <- c(-.024, -5.681, .005, -.056, -.169, -.12, .002, -.00000882, -.725)   # authors' regression.R start values
gamma_start <- c(90.523, .044, 3.569, -.549, 2.175)

sourceCpp(file.path(EC, "xiang_mcmc.cpp"), cacheDir = file.path(E, "cppcache"))
run_model <- function(nm, mcmc, burnin, thin, verbose = 0, state = NULL, seed = 0) {
  xv <- c(XV, ADD[[nm]]); G <- arr(GV); X <- arr(xv)
  b0 <- c(beta_start[1:8], rep(0, length(ADD[[nm]])), beta_start[9]); g0 <- gamma_start; am <- 0; av <- 1; ai <- numeric(0)
  if (!is.null(state)) { b0 <- state$beta; g0 <- state$gamma; am <- state$alpha_mean; av <- state$alpha_var; ai <- state$alpha }
  set.seed(seed)                                                               # authors' seed 0 (segment k: 0 + k - 1)
  t0 <- proc.time()[["elapsed"]]
  out <- MCMCdyadeffect.matrix(Ycs, G, X, mcmc, burnin, thin, NA_integer_, verbose, NA_character_, am, av, g0, b0, 1, ai)
  el <- proc.time()[["elapsed"]] - t0
  pn <- c("alpha_mean","alpha_var", paste0("rel_", GV), paste0("conf_", c(xv, "const")))
  par <- out[, 5116 + seq_along(pn), drop = FALSE]; colnames(par) <- pn
  last <- nrow(out); np <- length(pn)
  st <- list(alpha = out[last, 1:5116], alpha_mean = par[last, 1], alpha_var = par[last, 2],
             gamma = unname(par[last, 2 + seq_along(GV)]), beta = unname(par[last, 2 + length(GV) + seq_len(length(xv) + 1)]))
  list(par = par, secs = el, iters = mcmc + burnin, mcmc = mcmc, burnin = burnin, thin = thin, state = st)
}

if (args[1] == "time") {
  for (nm in c("m0", "mDM")) {
    r <- run_model(nm, 400, 100, 1)
    cat(sprintf("TIMING %s: %d iterations in %.1f s -> %.4f s/iter; authors' 1.1M iters = %.1f h\n", nm, r$iters, r$secs, r$secs / r$iters, r$secs / r$iters * 1.1e6 / 3600))
    print(round(colMeans(r$par[-(1:50), ]), 4))
  }
} else if (args[1] == "run") {
  nm <- args[2]; mc <- as.integer(args[3]); bi <- as.integer(args[4]); th <- as.integer(args[5])
  r <- run_model(nm, mc, bi, th)
  saveRDS(r, file.path(E, paste0("chain_", nm, ".rds")))
  b <- r$par[, "conf_smigoabi"]
  cat(sprintf("%s: %d iters (burnin %d, thin %d) in %.1f min; smigoabi post mean %.5f sd %.5f 95%% CrI [%.5f, %.5f]\n", nm, r$iters, bi, th, r$secs / 60,
              mean(b), sd(b), quantile(b, .025), quantile(b, .975)))
  print(round(rbind(mean = colMeans(r$par), sd = apply(r$par, 2, sd)), 5))
}
if (args[1] == "seg") {
  nm <- args[2]; k <- as.integer(args[3]); it <- as.integer(args[4]); th <- as.integer(args[5])
  sf <- function(kk) file.path(E, sprintf("seg_%s_%03d.rds", nm, kk))
  stopifnot(!file.exists(sf(k)), k == 1 || file.exists(sf(k - 1)))
  state <- if (k == 1) NULL else readRDS(sf(k - 1))$state
  r <- run_model(nm, it, 0, th, state = state, seed = k - 1)
  r$par <- r$par[-1, , drop = FALSE]      # drop the iteration-0 draw: kept draws are every `th` iterations, contiguous across segments
  stopifnot(nrow(r$par) == it / th)
  saveRDS(r, sf(k))
  cat(sprintf("SEGMENT %s %d: %d iters in %.1f min; smigoabi mean this seg %.5f\n", nm, k, it, r$secs / 60, mean(r$par[, "conf_smigoabi"])))
}
if (args[1] == "write") {
  source(file.path(CODE, "ea2/est_helpers.R")); suppressPackageStartupMessages(library(coda))
  burn <- as.integer(args[2]); S <- list()
  nseg <- min(sapply(c("m0", "mDM"), function(nm) length(list.files(E, sprintf("^seg_%s_[0-9]+\\.rds$", nm)))))
  cat("using", nseg, "segments per model (equal chain lengths)\n"); stopifnot(nseg > 0)
  for (nm in c("m0", "mDM")) {
    fs <- sort(list.files(E, sprintf("^seg_%s_[0-9]+\\.rds$", nm), full.names = TRUE))[seq_len(nseg)]
    rs <- lapply(fs, readRDS); th <- rs[[1]]$thin; stopifnot(all(sapply(rs, `[[`, "thin") == th))
    par <- do.call(rbind, lapply(rs, `[[`, "par")); tot <- nrow(par) * th
    keep <- par[(burn / th + 1):nrow(par), , drop = FALSE]
    b <- keep[, "conf_smigoabi"]
    S[[nm]] <- list(b = mean(b), se = sd(b), p = 2 * min(mean(b > 0), mean(b < 0)), lo = unname(quantile(b, .025)), hi = unname(quantile(b, .975)),
                    ess = unname(effectiveSize(mcmc(b))), draws = length(b), total = tot, kept_iter = tot - burn,
                    gew = unname(geweke.diag(mcmc(b))$z), mins = sum(sapply(rs, `[[`, "secs")) / 60, all = colMeans(keep), dep = mean(keep[, "conf_smldep"]))
    cat(sprintf("%s: total iters %d, burn-in %d, kept draws %d (thin %d); smigoabi mean %.5f sd %.5f CrI [%.5f, %.5f] ESS %.0f Geweke z %.2f; CPU %.0f min\n",
                nm, tot, burn, length(b), th, S[[nm]]$b, S[[nm]]$se, S[[nm]]$lo, S[[nm]]$hi, S[[nm]]$ess, S[[nm]]$gew, S[[nm]]$mins))
    print(round(S[[nm]]$all, 5))
  }
  stopifnot(S$m0$total == S$mDM$total)
  R <- list(m0 = c(b = S$m0$b, se = S$m0$se, p = S$m0$p), mDM = c(b = S$mDM$b, se = S$mDM$se, p = S$mDM$p))
  AE <- "Bayesian split-population (Boolean) probit, dyad-specific random effects in relevance eq. (conflict eq.: smldmat smldep smigoabi lcaprat2 allies py1-3); Gibbs sampler, burn-in 100k + 1M iter, thin 500 (Table 1 Model 3)"
  ES <- c(sprintf("authors' Gibbs sampler (dyadEffect), %dk iter (burn-in %dk), thin %d", S$m0$total / 1000, burn / 1000, rs[[1]]$thin),
          "NOT RUN (runtime >4h CPU)", sprintf("authors' Gibbs sampler (dyadEffect), %dk iter (burn-in %dk), thin %d; history in conflict eq.", S$mDM$total / 1000, burn / 1000, rs[[1]]$thin),
          "NOT RUN (runtime >4h CPU)", "NOT RUN (runtime >4h CPU)", "NOT RUN (runtime >4h CPU)")
  write_est("101", "Xiang 2016", "smigoabi", AE, "pdf (Xiang2016II Table 1 Model 3) + authors' code (dyadEffect/regression.R, replication.R)", ES, "N",
    "all six steps used pooled probit (Table 1 Model 1) instead of the headline Bayesian split-population probit with dyad effects (Model 3); M3 conflict eq. excludes contiguity/distance/major power (relevance eq.)",
    "PARTIAL", R, 114896,
    sprintf(paste0("m0 and mDM only (authors' 1.1M iterations = %.1f h per model at 0.017 s/iter; all six = ~31 h CPU on a shared single R lock). ",
      "Chains: %d total iterations each (authors 1.1M), first %d discarded, %d kept draws (thin %d); se = posterior SD; p = two-sided posterior tail prob 2*min(P(b>0),P(b<0)). ",
      "m0 95%% CrI [%.4f, %.4f], ESS %.0f, Geweke z %.2f; mDM 95%% CrI [%.4f, %.4f], ESS %.0f, Geweke z %.2f. Published M3 smigoabi -0.008 [-0.014, -0.002] (full N 116,013). ",
      "Common sample N=114,896 (1,117 obs lacking benchmark vars treated as missing Y). History enters the conflict equation. Trade/GDP smldep (conflict eq.): m0 %.3f -> mDM %.3f. ",
      "Previous pooled-probit numbers (results_verify): b0 +0.0133, pct_dm -56%%. MC error: short chain relative to authors'; treat decimals beyond MC SE (sd/sqrt(ESS)) as noise."),
      1.1e6 * 0.017 / 3600, S$m0$total, burn, S$m0$draws, rs[[1]]$thin, S$m0$lo, S$m0$hi, S$m0$ess, S$m0$gew, S$mDM$lo, S$mDM$hi, S$mDM$ess, S$mDM$gew, S$m0$dep, S$mDM$dep))
}
sink()
