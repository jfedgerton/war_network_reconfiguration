# X060 Beardsley (2024) JOP -- SAOM (RSiena) coevolution of support ties, threat ties and civil-war intensity,
# yearly waves 1980-2010, Model 1 (Fig 2a/3a; replication script 03Fig2a,3a,5&A2_JOPreplication.R).
# Usage: Rscript rerun_X060.R <m0|mH|m26|mP> [probe]
# Authors' algorithm: sienaAlgorithmCreate(useStdInits=F, nsub=4, n3=5000, maxlike=F, firstg=.03,
#   modelType=c(support.net=1, threat.net=1)); single siena07 call (no seed set; 'results may vary slightly').
# Deviations: seed 123; 2 cluster nodes; thetaBound raised if RSiena 1.6 requires it.
# History terms (threat network only): mH = X(hist_dyad_w15) + egoX/altX(log1p monadic MID onsets t-1..t-5);
#   m26 = same with t-2..t-6; mP = X(bm_cinc_w15) placebo. Covariates measured at wave year t for period t->t+1.
args <- commandArgs(TRUE); MODEL <- args[1]; PROBE <- length(args) > 1 && args[2] == "probe"
stopifnot(MODEL %in% c("m0","mH","m26","mP"))
ROOT <- Sys.getenv("SWEEP_DATA"); OUT <- Sys.getenv("SWEEP_OUT"); CODE <- Sys.getenv("SWEEP_CODE")
stopifnot(nzchar(ROOT), nzchar(OUT), nzchar(CODE))
WD <- file.path(OUT, "work_saom/X060"); dir.create(WD, recursive = TRUE, showWarnings = FALSE)
setwd(WD)   # all fits/res/RSiena project files are written here (pipeline output)
suppressPackageStartupMessages(library(dplyr))
source(file.path(CODE, "work_saom/X060/prep_data.R"))          # authors' data construction verbatim (reads the Dataverse .tab originals)
cat("RSiena", as.character(packageVersion("RSiena")), R.version.string, " N =", N, " waves =", nhi, "\n")
stopifnot(N >= 150, nhi == 31)

## ---- history covariates ----
H <- readRDS(Sys.getenv("SWEEP_HIST"))
source(file.path(CODE, "hist/add_history.R"))
ids <- as.numeric(rownames(emat))[order(id)]
stopifnot(all(ids == sort(id)))
dmat <- function(var) {
  arr <- array(NA_real_, dim = c(N, N, nlo))
  for (k in seq_along(yrlo)) {
    g <- expand.grid(a = seq_len(N), b = seq_len(N)); g$c1 <- ids[g$a]; g$c2 <- ids[g$b]; g$year <- yrlo[k]
    g <- add_history(g, "c1", "c2", "year"); v <- g[[var]]; v[is.nan(v)] <- NA
    m <- matrix(v, N, N); diag(m) <- 0; arr[, , k] <- m
  }
  arr
}
mmat <- function(var) {
  out <- matrix(NA_real_, N, nlo)
  for (k in seq_along(yrlo)) {
    s <- H$sth[H$sth$year == yrlo[k], ]
    out[, k] <- log1p(coalesce(s[[var]][match(ids, s$ccode)], 0))
  }
  out
}
if (MODEL == "mH")  { hd.net <- varDyadCovar(dmat("hist_dyad_w15")); hm.net <- varCovar(mmat("s15")) }
if (MODEL == "m26") { hd.net <- varDyadCovar(dmat("hist_dyad_w26")); hm.net <- varCovar(mmat("s26")) }
if (MODEL == "mP")  { bm <- dmat("bm_cinc_w15"); bm[is.na(bm)] <- mean(bm, na.rm = TRUE); hd.net <- varDyadCovar(bm) }
if (MODEL != "m0") {
  if (MODEL == "mP") {
    netdata <- sienaDataCreate(support.net, threat.net, intra_intens.net, mi1.net, energypc.net, mmp_change.net, dem.net,
      s_t_ipdistance.net, support_IGO.net, support_ipcloseness.net, support_mi.net, threat_mi.net, support_mi_pot.net,
      threat_density_nbrs.net, intra_intens_nbrs.net, contiguous.net, idealpointdistance.net, IGOcommon.net,
      support_i_j_ipdistance.net, hd.net, changes)
  } else {
    netdata <- sienaDataCreate(support.net, threat.net, intra_intens.net, mi1.net, energypc.net, mmp_change.net, dem.net,
      s_t_ipdistance.net, support_IGO.net, support_ipcloseness.net, support_mi.net, threat_mi.net, support_mi_pot.net,
      threat_density_nbrs.net, intra_intens_nbrs.net, contiguous.net, idealpointdistance.net, IGOcommon.net,
      support_i_j_ipdistance.net, hd.net, hm.net, changes)
  }
}
source(file.path(CODE, "work_saom/X060/prep_effects.R"))      # authors' effects verbatim -> eff2 (Model 1)
if (MODEL != "m0") {
  eff2 <- includeEffects(eff2, X, name = "threat.net", interaction1 = "hd.net")
  if (MODEL != "mP") eff2 <- includeEffects(eff2, egoX, altX, name = "threat.net", interaction1 = "hm.net")
}
if (MODEL != "m0" && file.exists("res_m0.rds")) {   # warm start from m0 (new effects start at 0)
  th0 <- readRDS("res_m0.rds")$res
  for (k in seq_len(nrow(th0))) {
    w <- which(eff2$include & eff2$effectName == th0$effect[k] & eff2$name == th0$dv[k] & eff2$type == th0$type[k])
    if (length(w) == 1) eff2$initialValue[w] <- th0$b[k]
  }
  cat("warm start from m0 estimates\n")
}
print(eff2)

model <- sienaAlgorithmCreate(useStdInits = F, projname = paste0("X060_", MODEL), nsub = if (PROBE) 1 else 4,
                              n3 = if (PROBE) 200 else 5000, maxlike = F, firstg = .03,
                              modelType = c(support.net = 1, threat.net = 1), seed = 123)
t0 <- Sys.time()
output <- siena07(model, data = netdata, effects = eff2, batch = TRUE, verbose = FALSE,
                  useCluster = TRUE, nbrNodes = 2, returnDeps = FALSE, thetaBound = 1000)
saveRDS(output, sprintf("fit_%s%s_raw.rds", MODEL, ifelse(PROBE, "_probe", "")))   # save before any post-processing
cat("first run: max|t conv|", max(abs(output$tconv)), " overall max conv ratio", output$tconv.max, "\n")
# not in authors' script: if not converged (authors report all t-ratios < 0.1), one continuation run with prevAns
if (!PROBE && (max(abs(output$tconv)) > 0.1 || output$tconv.max > 0.25)) {
  saveRDS(output, sprintf("fit_%s_run1.rds", MODEL))
  output <- siena07(model, data = netdata, effects = eff2, batch = TRUE, verbose = FALSE, prevAns = output,
                    useCluster = TRUE, nbrNodes = 2, returnDeps = FALSE, thetaBound = 1000)
  saveRDS(output, sprintf("fit_%s_raw.rds", MODEL))
  cat("continuation run: max|t conv|", max(abs(output$tconv)), " overall max conv ratio", output$tconv.max, "\n")
}
el <- as.numeric(difftime(Sys.time(), t0, units = "hours"))
print(output)
cat("\nMODEL", MODEL, "elapsed hours", el, " max|t conv|", max(abs(output$tconv)), " overall max conv ratio", output$tconv.max, "\n")
re <- output$requestedEffects; stopifnot(nrow(re) == length(output$theta))
res <- data.frame(effect = re$effectName, dv = re$name, type = re$type,
                  b = output$theta, se = output$se, tconv = output$tconv)
res$p <- 2 * pnorm(-abs(res$b / res$se))
saveRDS(list(res = res, tconv.max = output$tconv.max, hours = el, probe = PROBE, theta = output$theta),
        sprintf("res_%s%s.rds", MODEL, ifelse(PROBE, "_probe", "")))
saveRDS(output, sprintf("fit_%s%s.rds", MODEL, ifelse(PROBE, "_probe", "")))
