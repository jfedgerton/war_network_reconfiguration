# X085 Warren (2016) JPR -- SAOM (RSiena) coevolution of alliances, conflict, democracy, 1920-2000.
# Usage: Rscript rerun_X085.R <model> [probe]
#   model in m0 (authors' Table I baseline), mH (+hist w15 dyad X + mono egoPlusAltX on conflict net),
#            m26 (+hist w26), mP (+placebo lagged CINC-ratio dyadic X on conflict net)
#   probe  = single short siena07 call to time it
# Authors' settings (Baseline_1920_2000_full.R): 10-yr waves, con_lag 15, con_level 3, dem_level 5,
#   est_siena loop: n2start 20->800, n3 500->5000, diag 0.8, firstg 0.1, tconv<0.1, tcmax<0.25, <=100 reps.
# Deviations: 2 cores instead of 6; seed 123 in sienaModelCreate (authors set none); no out-of-sample prediction.
args <- commandArgs(TRUE); MODEL <- args[1]; PROBE <- length(args) > 1 && args[2] == "probe"
stopifnot(MODEL %in% c("m0","mH","m26","mP"))
ROOT <- Sys.getenv("SWEEP_DATA"); OUT <- Sys.getenv("SWEEP_OUT"); CODE <- Sys.getenv("SWEEP_CODE")
stopifnot(nzchar(ROOT), nzchar(OUT), nzchar(CODE))
WD <- file.path(OUT, "work_saom/X085"); dir.create(WD, recursive = TRUE, showWarnings = FALSE)
setwd(WD)   # fits, res_*.rds and RSiena project files are written/read here (pipeline output)
suppressPackageStartupMessages({library(RSiena); library(data.table); library(dplyr)})
cat("RSiena", as.character(packageVersion("RSiena")), R.version.string, "\n")

version_name <- "Baseline"; start_yr <- 1920; end_yr <- 2000; incr_yr <- 10
pred_sim <- FALSE; pred_yr <- 2010; est_min <- FALSE
con_lag <- 15; con_level <- 3; dem_level <- 5
num_cores <- 2; num_reps <- 100; num_diag <- 0.8; num_firstg <- 0.1
num_n2start <- c(20, 800); num_n3sims <- c(500, 5000); num_tconv <- 0.1; num_tcmax <- 0.25
source(file.path(CODE, "work_saom/X085/RSiena_helper_functions.R"))
# seed 123: override est_siena's sienaModelCreate by wrapping it
.smc <- RSiena::sienaModelCreate
# RSiena >=1.3 requires a named modelType; authors' scalar 3 applied to both networks (as in RSiena 1.1-289)
# seed: 123 for the first call, then 124, 125, ... (a fixed seed on every ladder rep makes prevAns reps identical,
# so a discarded rep would repeat forever; authors set no seed, so their reps differ)
.seedctr <- 122
sienaModelCreate <- function(..., modelType = NULL) { .seedctr <<- .seedctr + 1
  .smc(..., modelType = c(ally_est = 3, con_est = 3), seed = .seedctr) }
# RSiena 1.6 aborts if any initial value > thetaBound (default 50); 1940-50 alliance rate starts at ~59 (absent in 1.1-289)
siena07 <- function(...) RSiena::siena07(..., thetaBound = 1000)
source(file.path(CODE, "work_saom/X085/prep_warren.R"))   # authors' run_model.R data section verbatim -> est_data (baseline object)
stopifnot(exists("est_data"), length(est_data$nodeSets$States) > 100)
cat("States:", length(est_data$nodeSets$States), " waves:", length(est_years), "\n")

## ---- history / placebo covariates (conflict network only) ----
H <- readRDS(Sys.getenv("SWEEP_HIST"))
source(file.path(CODE, "hist/add_history.R"))
wave_yrs <- est_years[-length(est_years)]          # covariates for period t->t+1 measured at wave t
sts <- sort(unique(state_year$cow_ccode[state_year$year %in% wave_yrs]))
pairs <- expand.grid(cow_ccode1 = sts, cow_ccode2 = sts, year = wave_yrs)
pairs <- pairs[pairs$cow_ccode1 != pairs$cow_ccode2, ]
pairs <- add_history(pairs, "cow_ccode1", "cow_ccode2", "year")
stopifnot(all(pairs$year >= 1822), all(!is.na(pairs$hist_dyad_w15)))
hdyad15_data <- transmute(pairs, cow_ccode1, cow_ccode2, year, val = hist_dyad_w15)
hdyad26_data <- transmute(pairs, cow_ccode1, cow_ccode2, year, val = hist_dyad_w26)
bmcinc_data  <- transmute(pairs, cow_ccode1, cow_ccode2, year, val = bm_cinc_w15)
bmcinc_data$val[is.nan(bmcinc_data$val)] <- NA
sy <- expand.grid(cow_ccode = sts, year = wave_yrs) %>%
  left_join(H$sth, by = c("cow_ccode" = "ccode", "year")) %>%
  mutate(s15 = coalesce(s15, 0), s26 = coalesce(s26, 0))
hmono15_data <- transmute(sy, cow_ccode, year, val = log1p(s15))
hmono26_data <- transmute(sy, cow_ccode, year, val = log1p(s26))
cat("dyad hist w15 >0 share:", mean(pairs$hist_dyad_w15 > 0), " mono w15 mean:", mean(hmono15_data$val), "\n")

add_iv <- switch(MODEL, m0 = character(0), mH = "hmono15", m26 = "hmono26", mP = character(0))
add_dv <- switch(MODEL, m0 = character(0), mH = "hdyad15", m26 = "hdyad26", mP = "bmcinc")
if (MODEL == "m0") {
  dat <- est_data
} else {
  dat <- formatData(ee = entry_exit_data, yrs = est_years, ndv = ndv_list, bdv = bdv_list,
                    iv = c(iv_list, add_iv), dyadv = c(dyadv_list, add_dv), dyadc = dyadc_list,
                    ccodes = ccodes, ccode = ccode)
}
efx <- createFullEffects(data = dat, ndv = ndv_list, bdv = bdv_list, iv = iv_list,
                         dyadv = dyadv_list, dyadc = dyadc_list)   # authors' effects only
for (v in add_dv) efx <- includeEffects(efx, X, name = "con_est", interaction1 = paste0(v, "_est"))
for (v in add_iv) efx <- includeEffects(efx, egoPlusAltX, name = "con_est", interaction1 = paste0(v, "_est"))
print(efx)

t0 <- Sys.time()
if (PROBE) {
  mod <- sienaModelCreate(projname = paste0("probe_", MODEL), modelType = 3, truncation = 10, nsub = 1,
                          n2start = 20, n3 = 500, diagonalize = 0.8, firstg = 0.1)
  est <- siena07(mod, data = dat, effects = efx, useCluster = TRUE, nbrNodes = num_cores, batch = TRUE, verbose = FALSE)
} else {
  # Authors' est_siena ladder (RSiena_helper_functions.R) re-implemented with these documented changes:
  #  (1) restart rule: each rep starts from the last accepted theta via updateTheta with a fresh phase 1 (the paper:
  #      "repeated, using the coefficient estimates from the previous subphase as the starting values"); passing
  #      prevAns diverged on every retry under RSiena 1.6.6 (tcmax 85-97 after a rep-1 tcmax of 10.8);
  #  (2) seed 123, 124, ... per rep (authors set none; a fixed seed makes discarded reps repeat identically);
  #  (3) a rep that exceeds REP_MIN minutes is discarded like a diverged rep (one rep stalled >1 h in repeated phase-1
  #      restarts), every accepted rep is saved to disk, and the ladder stops at CAP_H hours.
  CAP_H <- as.numeric(Sys.getenv("CAP_H", "5")); REP_MIN <- as.numeric(Sys.getenv("REP_MIN", "40"))
  est_siena <- function(data, efx, nreps, n2start, n3sims, diag, firstg, ntconv, ntcmax, ncpu, stub) {
    st <- Sys.time(); n3 <- n3sims[1]; n2 <- n2start[1]; est <- NULL
    n2_switched <- FALSE; n3_switched <- FALSE; tcmax_last <- 20
    # RESUME=1: continue the ladder from the last accepted rep saved on disk (state flags restored from env)
    if (Sys.getenv("RESUME") == "1" && file.exists(sprintf("fit_%s_best.rds", MODEL))) {
      est <- readRDS(sprintf("fit_%s_best.rds", MODEL)); tcmax_last <- est$tconv.max
      n2_switched <- Sys.getenv("N2SW") == "1"; n3_switched <- Sys.getenv("N3SW") == "1"
      if (n2_switched) n2 <- n2start[2]; if (n3_switched) n3 <- n3sims[2]
      cat("RESUMED from saved rep: tconv.max", tcmax_last, " n2", n2, " n3", n3, "\n")
    }
    for (rep in 1:nreps) {
      if (as.numeric(difftime(Sys.time(), st, units = "hours")) > CAP_H) { cat("WALL-CLOCK CAP reached\n"); break }
      if (!is.null(est)) efx <- updateTheta(efx, est)
      model <- sienaModelCreate(projname = paste0("AllyConDem_EST_", stub), modelType = 3, truncation = 10, nsub = 1,
                                n2start = n2, n3 = n3, diagonalize = diag, firstg = firstg)
      t_rep <- Sys.time()
      new_est <- tryCatch({
        setTimeLimit(elapsed = REP_MIN * 60, transient = TRUE)
        siena07(model, data = data, effects = efx, useCluster = as.logical(ncpu > 1), nbrNodes = ncpu)
      }, error = function(e) { cat("rep", rep, "aborted:", conditionMessage(e), "\n"); NULL },
         finally = setTimeLimit(elapsed = Inf))
      if (is.null(new_est)) { closeAllConnections(); cat("Discarding estimate (timeout/error)\n\n"); next }
      tconv <- max(abs(new_est$tconv)); tcmax <- new_est$tconv.max
      cat("rep", rep, "n2start", n2, "n3", n3, " Max tconv:", tconv, " Max tconv ratio:", tcmax,
          " rep minutes:", round(as.numeric(difftime(Sys.time(), t_rep, units = "mins")), 1),
          " total hours:", round(as.numeric(difftime(Sys.time(), st, units = "hours")), 2), "\n")
      if (is.na(tcmax) || tcmax > tcmax_last * 1.5) { cat("Discarding estimate\n\n"); next }
      est <- new_est; tcmax_last <- tcmax
      saveRDS(est, sprintf("fit_%s_best.rds", MODEL))
      if (tconv < ntconv * 3 && tcmax < ntcmax * 3 && !n2_switched) { cat("Switching to n2start =", n2start[2], "\n"); n2 <- n2start[2]; n2_switched <- TRUE; next }
      if (tconv < ntconv * 2 && tcmax < ntcmax * 2 && !n3_switched) { cat("Switching to n3 =", n3sims[2], "\n"); n3 <- n3sims[2]; n3_switched <- TRUE; next }
      if (tconv < ntconv && tcmax < ntcmax) { cat("Estimation converged\n"); break }
    }
    est
  }
  if (MODEL != "m0" && file.exists("res_m0.rds")) {
    th0 <- readRDS("res_m0.rds")$res
    for (k in seq_len(nrow(th0))) {
      w <- which(efx$include & efx$effectName == th0$effect[k] & efx$name == th0$dv[k] & efx$type == th0$type[k])
      if (length(w) == 1) efx$initialValue[w] <- th0$b[k]
    }
    cat("Started from m0 estimates\n")
  }
  est <- est_siena(data = dat, efx = efx, nreps = num_reps, n2start = num_n2start, n3sims = num_n3sims,
                   diag = num_diag, firstg = num_firstg, ntconv = num_tconv, ntcmax = num_tcmax,
                   ncpu = num_cores, stub = paste0("X085_", MODEL))
}
el <- as.numeric(difftime(Sys.time(), t0, units = "hours"))
saveRDS(est, sprintf("fit_%s%s_raw.rds", MODEL, ifelse(PROBE, "_probe", "")))   # save before any post-processing
print(est)
cat("\nMODEL", MODEL, "elapsed hours", el, " max|t conv|", max(abs(est$tconv)), " overall max conv ratio", est$tconv.max, "\n")
re <- est$requestedEffects; stopifnot(nrow(re) == length(est$theta))
res <- data.frame(effect = re$effectName, dv = re$name, type = re$type,
                  b = est$theta, se = sqrt(diag(est$covtheta)), tconv = est$tconv)
res$p <- 2 * pnorm(-abs(res$b / res$se))
saveRDS(list(res = res, tconv.max = est$tconv.max, hours = el, probe = PROBE),
        sprintf("res_%s%s.rds", MODEL, ifelse(PROBE, "_probe", "")))
