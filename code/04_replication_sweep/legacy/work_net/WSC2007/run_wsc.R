# WSC2007 rerun: gbme (Hoff 2003) binomial, k=3, directed, year 2000; authors' uninformative priors (routine_uninform.R)
# usage: Rscript run_wsc.R <model: m0|hist|m26|placebo> <NS>
args <- commandArgs(TRUE); mod <- args[1]; NS <- as.integer(args[2])
ROOT <- Sys.getenv("SWEEP_DATA"); OUT <- Sys.getenv("SWEEP_OUT"); CODE <- Sys.getenv("SWEEP_CODE")
stopifnot(nzchar(ROOT), nzchar(OUT), nzchar(CODE))
W <- file.path(ROOT, "work_net/WSC2007")   # authors' data + archived chain (read-only)
WO <- file.path(OUT, "work_net/WSC2007"); dir.create(WO, recursive = TRUE, showWarnings = FALSE)   # chain output dirs
d <- dget(file.path(W, "uninform5k_2000/data2000"))
Yd <- d[[1]]; Xd <- d[[2]]; Xs <- d[[3]]
nm <- rownames(Xs); stopifnot(identical(nm, rownames(Yd)), identical(nm, dimnames(Xd)[[1]]))
# COW abbreviation -> ccode
cc <- read.csv(file.path(W, "cow_states.csv"))
cc <- unique(cc[, c("stateabb","ccode")])
map <- setNames(cc$ccode, cc$stateabb); map <- c(map, FJI = 950, RUM = 360, SER = 345)  # GW-style abbreviations
ccode <- unname(map[nm]); stopifnot(!any(is.na(ccode)))
H <- readRDS(Sys.getenv("SWEEP_HIST"))
yr <- 2000; n <- length(nm)
dyadmat <- function(col) {
  m <- matrix(0, n, n); h <- H$dyh[H$dyh$year == yr, ]
  for (i in 1:n) for (j in 1:n) if (i != j) {
    a <- min(ccode[i], ccode[j]); b <- max(ccode[i], ccode[j])
    v <- h[[col]][h$c1 == a & h$c2 == b]; if (length(v)) m[i, j] <- v }
  log1p(m) }
monovec <- function(col) { h <- H$sth[H$sth$year == yr, ]; v <- h[[col]][match(ccode, h$ccode)]; log1p(ifelse(is.na(v), 0, v)) }
cincr <- function() { # mean over t-1..t-5 of stronger side's CINC share (as in add_history.R)
  acc <- array(NA, c(n, n, 5))
  for (L in 1:5) { cv <- H$nmc$cinc[H$nmc$year == yr - L][match(ccode, H$nmc$ccode[H$nmc$year == yr - L])]
    acc[, , L] <- outer(cv, cv, pmax) / outer(cv, cv, "+") }
  m <- apply(acc, c(1, 2), mean, na.rm = TRUE); m[is.na(m)] <- mean(m, na.rm = TRUE); diag(m) <- 0; m }
addXd <- function(Xd, M, lab) { X2 <- array(0, c(n, n, dim(Xd)[3] + 1)); X2[, , 1:dim(Xd)[3]] <- Xd; X2[, , dim(Xd)[3] + 1] <- M
  dimnames(X2) <- list(nm, nm, c(dimnames(Xd)[[3]], lab)); X2 }
if (mod == "hist")    { Xd <- addXd(Xd, dyadmat("d15"), "hist_dyad_w15"); Xs <- cbind(Xs, hist_mono_w15 = monovec("s15")) }
if (mod == "m26")     { Xd <- addXd(Xd, dyadmat("d26"), "hist_dyad_w26"); Xs <- cbind(Xs, hist_mono_w26 = monovec("s26")) }
if (mod == "placebo") { Xd <- addXd(Xd, cincr(), "bm_cinc_w15") }
stopifnot(!any(is.na(Xd)), !any(is.na(Xs)))
cat(mod, ": dyadic layers", dimnames(Xd)[[3]], "| nodal", colnames(Xs), "\n")
if (mod %in% c("hist","m26")) { cat("share of dyads with history>0:", mean(Xd[, , 5] > 0), "; MIDs in 2000 on dyads with history:", sum(Yd[Xd[, , 5] > 0]), "of", sum(Yd), "\n") }
Xr <- Xs; Y <- Yd  # gbme helpers read global Y (as in authors script)
source(file.path(CODE, "work_net/WSC2007/gbme_ws.r"))  # authors' gbme.r + warm-start initialization of e, f, E
rd <- dim(Xd)[3]; rs <- ncol(Xs)
# authors' priors (routine_uninform.R), extended with the same diffuse settings for added terms
pi.Sab <- c(100, 10, 100, 4); pi.s2u <- c(2, 1); pi.s2v <- c(2, 1); pi.s2z <- matrix(nrow = 3, ncol = 2, c(rep(2, 3), rep(1, 3)))
pim.bd <- rep(0, rd); piS.bd <- diag(100, rd)
pim.b0sr <- c(-20, rep(0, 2 * rs)); piS.b0sr <- diag(c(400, rep(100, 2 * rs)))
priors <- list(pi.Sab = pi.Sab, pi.s2u = pi.s2u, pi.s2v = pi.s2v, pi.s2z = pi.s2z, pim.bd = pim.bd, piS.bd = piS.bd,
               pim.b0sr = pim.b0sr, piS.b0sr = piS.b0sr, pim.b0s = NULL, piS.b0s = NULL, pi.s2a = NULL)
# warm start: the LAST saved draw of the authors' uninformative-prior Y2K chain (uninform5k_2000, 497k scans):
# beta.d, beta.u, additive effects (A,B + sr.hat), bilinear e,f (U,V); dyadic residual E drawn from its prior at the authors' last (se, rho).
# (a cold start with E=0 needs ~500k scans for se to grow from ~0 to ~150-300; our first 25k-scan attempt had not converged: see coldstart_*)
set.seed(123)
ao <- read.table(file.path(W, "uninform5k_2000/OUT"), header = TRUE); last <- ao[nrow(ao), ]
A0 <- as.matrix(read.table(file.path(W, "uninform5k_2000/A"))); B0 <- as.matrix(read.table(file.path(W, "uninform5k_2000/B")))
stopifnot(ncol(A0) == n, nrow(A0) == nrow(ao), nrow(B0) == nrow(ao))
U0 <- as.matrix(read.table(file.path(W, "uninform5k_2000/U"))); V0 <- as.matrix(read.table(file.path(W, "uninform5k_2000/V")))
stopifnot(nrow(U0) == n * nrow(ao), ncol(U0) == 3, nrow(V0) == nrow(U0))
e0 <- U0[nrow(U0) - n + 1:n, ]; f0 <- V0[nrow(V0) - n + 1:n, ]; rm(U0, V0)
a0 <- A0[nrow(A0), ]; b0 <- B0[nrow(B0), ]
bd0 <- c(unlist(last[paste0("bd", 1:4)]), rep(0, rd - 4))
bu0 <- c(last$b0, unlist(last[paste0("bs", 1:3)]), rep(0, rs - 3), unlist(last[paste0("br", 1:3)]), rep(0, rs - 3))
Xu0 <- rbind(cbind(rep(.5, n), Xs, matrix(0, n, rs)), cbind(rep(.5, n), matrix(0, n, rs), Xr)); srhat <- Xu0 %*% bu0
se0 <- last$s2e; rho0 <- last$rho; su0 <- 2 * se0 * (1 + rho0); sv0 <- 2 * se0 * (1 - rho0)
ut <- upper.tri(matrix(0, n, n)); UU <- VV <- matrix(0, n, n)
uu <- rnorm(sum(ut), 0, sqrt(su0)); vv <- rnorm(sum(ut), 0, sqrt(sv0))
UU[ut] <- uu; UU <- UU + t(UU); VV[ut] <- vv; VV <- VV - t(VV); E0 <- (UU + VV) / 2; diag(E0) <- 0
cat("authors' last draw: scan", last$scan, "ll", last$ll, "se", se0, "rho", rho0, "\n")
startv <- list(beta.d = bd0, beta.u = bu0, s = a0 + srhat[1:n], r = b0 + srhat[n + 1:n], z = matrix(0, n, 3), e = e0, f = f0, E = E0)
cat("start beta.d:", round(bd0, 3), "\n")
od <- file.path(WO, paste0(Sys.getenv("WSC_OUTPREFIX", "out_"), mod)); dir.create(od, showWarnings = FALSE); setwd(od)
saveRDS(list(Xd_names = dimnames(Xd)[[3]], Xs_names = colnames(Xs)), "names.rds")
t0 <- Sys.time()
gbme(Y = Yd, Xd = Xd, Xs = Xs, Xr = Xr, k = 3, fam = "binomial", N = matrix(1, n, n), NS = NS, estart = F, directed = T,
     seed = 123, zwrite = FALSE, awrite = F, bwrite = F, odens = 20)
cat("elapsed", format(Sys.time() - t0), "\n")
