# MINHAS2022 rerun: authors' ame_repL (amen pa2018_version), R=2, probit; warm start from authors' final MCMC state (ameFit$startVals)
# usage: Rscript run_minhas.R <app: rs|weeks|gibler> <model: m0|hist> <burn> <nscan> <odens> <startyear>
args <- commandArgs(TRUE); app <- args[1]; mod <- args[2]; burn <- as.integer(args[3]); nscan <- as.integer(args[4]); odens <- as.integer(args[5]); y0 <- as.integer(args[6])
LIB_AMEN <- Sys.getenv("SWEEP_LIB_AMEN"); stopifnot(nzchar(LIB_AMEN), dir.exists(file.path(LIB_AMEN, "amen")))   # amen 1.4, s7minhas/amen@pa2018_version
.libPaths(c(LIB_AMEN, .libPaths())); library(amen)
stopifnot(packageVersion("amen") == "1.4", exists("rUV_rep_fc_cpp", asNamespace("amen")))
ROOT <- Sys.getenv("SWEEP_DATA"); OUT <- Sys.getenv("SWEEP_OUT"); CODE <- Sys.getenv("SWEEP_CODE")
stopifnot(nzchar(ROOT), nzchar(OUT), nzchar(CODE))
W <- file.path(OUT, "work_net/MINHAS2022"); dir.create(W, recursive = TRUE, showWarnings = FALSE)   # fit_<app>_<model>.rds -> OUT
P <- file.path(CODE, "work_net/MINHAS2022/x/replArchive")   # authors' helpers (code)
D <- file.path(ROOT, "work_net/MINHAS2022/x/replArchive/2_applications/application_data")   # authors' data (read-only)
source(file.path(P, "helpers/ame_repL.R"))
environment(ame_repL) <- asNamespace("amen")  # authors run it with amen attached; namespace gives access to internal helpers
cfg <- list(rs = list(dat = "reiter_stam/reiterStamData.rda", fit = "reiter_stam/ameFitReiterStam.rda", seed = 6886, sym = FALSE, nodes = FALSE),
            weeks = list(dat = "weeks/weeksData.rda", fit = "weeks/ameFitWeeks.rda", seed = 6886, sym = FALSE, nodes = TRUE),
            gibler = list(dat = "gibler/giblerData.rda", fit = "gibler/ameFitGibler_s119.rda", seed = 119, sym = TRUE, nodes = FALSE))[[app]]
if (app == "weeks") {   # authors' weeks_ameRun.R swaps in their rZ_bin_fc
  rZ_bin_fc2 <- dget(file.path(P, "helpers/rZBinfc.R")); assignInNamespace("rZ_bin_fc", rZ_bin_fc2, ns = "amen") }
load(file.path(D, cfg$dat)); load(file.path(D, cfg$fit)); sv <- ameFit$startVals; rm(ameFit); gc()
Tn <- length(yList); years <- y0 + seq_len(Tn) - 1
# GML history window t-1..t-5 needs onsets from 1816: keep years >= 1822 (same restricted sample for m0 and hist)
keep <- years >= 1822; yList <- yList[keep]; xDyadList <- xDyadList[keep]
if (exists("xNodeList.s")) { xNodeList.s <- xNodeList.s[keep]; xNodeList.r <- xNodeList.r[keep] }
if (any(!keep)) { stopifnot(dim(sv$Z)[3] == length(keep)); sv$Z <- sv$Z[, , keep, drop = FALSE] }
stopifnot(length(sort(unique(unlist(lapply(yList, rownames))))) == dim(sv$Z)[1])   # actor set unchanged by the restriction
cat("dropped", sum(!keep), "early years\n"); years <- years[keep]; Tn <- length(yList)
cat(app, mod, ": T =", Tn, "years", min(years), "-", max(years), "\n")
H <- readRDS(Sys.getenv("SWEEP_HIST"))
stopifnot(mod %in% c("m0", "hist", "m26", "placebo"))
cincr <- function(cc, yr) {  # mean over t-1..t-5 of stronger side's CINC share (as in add_history.R)
  n <- length(cc); acc <- array(NA, c(n, n, 5))
  for (L in 1:5) { cv <- H$nmc$cinc[H$nmc$year == yr - L][match(cc, H$nmc$ccode[H$nmc$year == yr - L])]
    acc[, , L] <- outer(cv, cv, pmax) / outer(cv, cv, "+") }
  apply(acc, c(1, 2), mean, na.rm = TRUE) }
if (mod != "m0") {
  dcol <- c(hist = "d15", m26 = "d26")[mod]; scol <- c(hist = "s15", m26 = "s26")[mod]
  labs <- switch(mod, hist = c("hist_dyad_w15", "hist_mono_w15"), m26 = c("hist_dyad_w26", "hist_mono_w26"), placebo = "bm_cinc_w15")
  nfill <- 0
  for (t in seq_len(Tn)) {
    cc <- as.numeric(rownames(yList[[t]])); n <- length(cc); yr <- years[t]
    X <- xDyadList[[t]]; mask <- is.na(X[, , 1])
    if (mod == "placebo") {
      cm <- cincr(cc, yr); nfill <- nfill + sum(is.na(cm[!mask]) & row(cm)[!mask] != col(cm)[!mask]); cm[is.na(cm)] <- mean(cm, na.rm = TRUE)
      cm[mask] <- NA; new <- list(cm)
    } else {
      dd <- matrix(0, n, n); h <- H$dyh[H$dyh$year == yr, ]
      if (nrow(h)) { key <- paste(h$c1, h$c2); A <- outer(cc, cc, pmin); B <- outer(cc, cc, pmax)
        v <- h[[dcol]][match(paste(A, B), key)]; dd[] <- ifelse(is.na(v), 0, v) }
      s <- H$sth[[scol]][H$sth$year == yr][match(cc, H$sth$ccode[H$sth$year == yr])]; s[is.na(s)] <- 0
      ee <- matrix(0, n, n); he <- H$exh[H$exh$year == yr, ]; ecol <- c(hist = "e15", m26 = "e26")[mod]
      if (nrow(he)) { v <- he[[ecol]][match(paste(outer(cc, cc, pmin), outer(cc, cc, pmax)), paste(he$c1, he$c2))]; ee[] <- ifelse(is.na(v), 0, v) }
      mono <- outer(s, s, "+") - ee; stopifnot(all(mono >= 0))   # participant-level opponent exclusion (2026-10-01)
      hd <- log1p(dd); hm <- log1p(mono); hd[mask] <- NA; hm[mask] <- NA; new <- list(hd, hm)
    }
    X2 <- array(NA, c(n, n, dim(X)[3] + length(new)), dimnames = list(rownames(X), colnames(X), c(dimnames(X)[[3]], labs)))
    X2[, , 1:dim(X)[3]] <- X; for (k in seq_along(new)) X2[, , dim(X)[3] + k] <- new[[k]]; xDyadList[[t]] <- X2
  }
  sv$beta <- c(sv$beta, rep(0, length(labs)))
  cat("added", labs, "; share of obs with first added term > 0:", mean(unlist(lapply(xDyadList, function(a) a[, , labs[1]] > 0)), na.rm = TRUE), "\n")
  if (mod == "placebo") cat("placebo: dyad-years with missing lagged CINC filled with the year mean:", nfill, "\n")
}
startVals0 <- sv
t0 <- Sys.time()
if (cfg$nodes) {
  fit <- ame_repL(Y = yList, Xdyad = xDyadList, Xrow = xNodeList.s, Xcol = xNodeList.r, model = "bin", symmetric = FALSE, R = 2,
                  nscan = nscan, seed = cfg$seed, burn = burn, odens = odens, plot = FALSE, print = FALSE, gof = FALSE, startVals = sv)
} else {
  fit <- ame_repL(Y = yList, Xdyad = xDyadList, Xrow = NULL, Xcol = NULL, model = "bin", symmetric = cfg$sym, intercept = TRUE, R = 2,
                  nscan = nscan, seed = cfg$seed, burn = burn, odens = odens, plot = FALSE, print = FALSE, gof = FALSE, startVals = sv)
}
cat("\nelapsed", format(Sys.time() - t0), "\n")
B <- fit$BETA; res <- data.frame(term = colnames(B), mean = colMeans(B), sd = apply(B, 2, sd),
  q025 = apply(B, 2, quantile, .025), q975 = apply(B, 2, quantile, .975), ess = coda::effectiveSize(coda::mcmc(B)),
  geweke_z = coda::geweke.diag(coda::mcmc(B))$z)
print(res, digits = 3)
saveRDS(list(res = res, BETA = B, VC = fit$VC, settings = list(burn = burn, nscan = nscan, odens = odens, seed = cfg$seed, years = range(years)), tryErrorChecks = fit$tryErrorChecks),
        file.path(W, paste0("fit_", app, "_", mod, ".rds")))
