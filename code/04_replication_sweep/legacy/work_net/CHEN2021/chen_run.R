# CHEN2021: multilayer ERGM (multilayer.ergm), Table 3 cross-layer dependence model (+ history variants).
# Usage: Rscript chen_run.R <model: m0|hist|m26|placebo> <method: MPLE|MCMLE> [ncores]
ROOT <- Sys.getenv("SWEEP_DATA"); OUT <- Sys.getenv("SWEEP_OUT"); CODE <- Sys.getenv("SWEEP_CODE")
stopifnot(nzchar(ROOT), nzchar(OUT), nzchar(CODE))
W <- file.path(OUT, "work_net/CHEN2021"); dir.create(W, recursive = TRUE, showWarnings = FALSE)   # res/diag rds -> OUT
OLD <- Sys.getenv("CHEN_OLD") == "1"   # 1 = author's stack: ergm 3.10.1 + multilayer.ergm 0.1.6 (lib_chen)
if (OLD) { LIB_CHEN <- Sys.getenv("SWEEP_LIB_CHEN"); stopifnot(nzchar(LIB_CHEN), dir.exists(file.path(LIB_CHEN, "ergm"))); .libPaths(c(LIB_CHEN, .libPaths())) }
suppressPackageStartupMessages({library(multilayer.ergm); library(dplyr)})
if (OLD) stopifnot(packageVersion("ergm") == "3.10.1", packageVersion("multilayer.ergm") == "0.1.6", packageVersion("statnet.common") == "4.3.0")
cat("ergm", as.character(packageVersion("ergm")), " multilayer.ergm", as.character(packageVersion("multilayer.ergm")), "\n")
source(file.path(CODE, "hist/add_history.R"))
a <- commandArgs(TRUE); mod <- a[1]; meth <- a[2]; ncores <- if (length(a) > 2) as.integer(a[3]) else 1
load(file.path(ROOT, "work_net/CHEN2021/data"))   # Dataverse multilayer_data (RData): mid, intraconf, nsconf, nodeatts, demdem, contig, polrel, caprat, intoverlap
stopifnot(dim(mid) == c(189, 189), dim(nsconf) == c(544, 544), nrow(nodeatts) == 733)
set.seed(20842)
nw <- cbind(rbind(mid, t(intraconf)), rbind(intraconf, nsconf))
nw <- network(nw, directed = FALSE)
nw %v% "layer.mem" <- c(rep(1, 189), rep(2, 544))
for (v in c("dem","cinc","trans","interstate","intrastate","nonstate")) nw %v% v <- nodeatts[[v]]
absent <- matrix(0, ncol = 733, nrow = 733); absent[1:189, 1:189] <- (1 - intoverlap); diag(absent) <- 1
absent <- network(absent, directed = FALSE)
# ---- conflict history: MID onsets before the 1989-1995 period (t = 1989: onsets 1984-88; m26: 1983-87) ----
cc <- as.numeric(rownames(mid)); n <- length(cc); H <- readRDS(Sys.getenv("SWEEP_HIST"))
g <- expand.grid(i = 1:n, j = 1:n); g <- g[g$i < g$j, ]
dd <- add_history(data.frame(c1 = cc[g$i], c2 = cc[g$j], year = 1989), "c1", "c2", "year")
M189 <- function(v) { x <- matrix(0, n, n); x[cbind(g$i, g$j)] <- v; x[cbind(g$j, g$i)] <- v; x }  # layer-1 sized, like demdem/contig
hd15 <- M189(dd$hist_dyad_w15); hd26 <- M189(dd$hist_dyad_w26)
cr <- dd$bm_cinc_w15; cat("bm_cinc NA dyads (mean-imputed):", sum(is.na(cr)), "\n"); cr[is.na(cr)] <- mean(cr, na.rm = TRUE); cincr <- M189(cr)
s <- H$sth[H$sth$year == 1989, ]
nw %v% "mono15" <- c(log1p(coalesce(s$s15[match(cc, s$ccode)], 0)), rep(0, 544))
nw %v% "mono26" <- c(log1p(coalesce(s$s26[match(cc, s$ccode)], 0)), rep(0, 544))
cat("MID ties:", sum(mid) / 2, " with prior dyadic onsets (1984-88):", sum(mid[upper.tri(mid)] * (hd15[upper.tri(mid)] > 0)), "\n")
f0 <- "nw ~ edgecov_layer('edges', layer = 1) + degree_layer(0, layer = 1) + altkstar.fixed_layer(layer = 1) +
  nodefactor_layer('dem', layer = 1) + nodecov_layer('cinc', layer = 1) + edgecov_layer(demdem, layer = 1) + edgecov_layer(contig, layer = 1) +
  edgecov_layer(polrel, layer = 1) + edgecov_layer(caprat, layer = 1) + edgecov_layer('edges', layer = c(1, 2)) + altkstar.fixed_layer(layer = c(1,2)) +
  nodefactor_layer('dem', layer = c(1, 2)) + nodefactor_layer('trans', layer = c(1, 2)) + edgecov_layer('edges', layer = 2) + degree_layer(0, layer = 2) +
  kstar_layer(2, layer = 2) + altkstar.fixed_layer(layer = 2) + altkstar.fixed_crosslayer(lambda = 1, layers = list(1,c(1,2))) +
  altkstar.fixed_crosslayer(lambda = 1, layers = list(c(1,2),2)) + threetrail_crosslayer(layers = list(1,c(1,2),2), incident = 1)"
f <- switch(mod, m0 = f0,
  hist = paste(f0, "+ edgecov_layer(hd15, layer = 1) + nodecov_layer('mono15', layer = 1)"),
  m26 = paste(f0, "+ edgecov_layer(hd26, layer = 1) + nodecov_layer('mono26', layer = 1)"),
  placebo = paste(f0, "+ edgecov_layer(cincr, layer = 1)"))
# authors' MCMC settings (seed 30748, burnin 1e5, samplesize 2e4, interval 1e4; parallel 4 -> ncores here);
# ergm 3.10.1 was used by the author: we mimic its non-adaptive MCMLE (effectiveSize NULL, Hummel termination, maxit 20).
# Budget run (OLD stack): authors' seed and burn-in, but shorter chains (authors' 2e8 proposals/iteration = ~3.8 CPU-h/iteration here).
# Args 4-6: interval, samplesize, maxit (defaults 2000, 4000, 8). Everything else = ergm 3.10.1 defaults, as the author used.
iv <- if (length(a) > 3) as.integer(a[4]) else 2000L; ss <- if (length(a) > 4) as.integer(a[5]) else 4000L; mx <- if (length(a) > 5) as.integer(a[6]) else 8L
if (OLD && meth != "MPLE") cat("MCMLE settings: burnin 1e5, interval", iv, " samplesize", ss, " maxit", mx, " parallel", ncores, "\n")
ctrl <- if (meth == "MPLE") control.ergm(seed = 30748) else if (OLD)
  control.ergm(seed = 30748, MCMC.burnin = 100000, MCMC.samplesize = ss, MCMC.interval = iv, MCMLE.maxit = mx, parallel = ncores) else
  control.ergm(seed = 30748, MCMC.burnin = 100000, MCMC.samplesize = 20000, MCMC.interval = 10000, parallel = ncores,
               MCMLE.effectiveSize = NULL, MCMLE.termination = "Hummel", MCMLE.maxit = 20)
if (meth == "DIAG") {
  # Post-hoc MCMC diagnostics at the final MCMLE estimates (the MCMLE runs did not keep the sample: the as.mcmc.list call lacked coda::).
  # One chain, authors' burn-in 1e5, same interval as the fit (iv), ss draws; ESS + Geweke per statistic, and convergence t-ratio (mean sim - obs)/sd.
  r <- readRDS(file.path(W, paste0("res_", mod, "_MCMLE_old.rds"))); th <- setNames(r$b, r$term)
  t0 <- Sys.time()
  sim <- simulate(as.formula(f), coef = th, nsim = ss, statsonly = TRUE, constraints = ~fixedas(absent = absent), seed = 30748,
                  control = control.simulate.formula(MCMC.burnin = 100000, MCMC.interval = iv))
  obs <- summary(as.formula(f)); stopifnot(identical(names(obs), colnames(sim)), length(obs) == length(th))
  mc <- coda::mcmc(sim)
  dg <- data.frame(term = names(th), obs = obs, sim_mean = colMeans(sim), t_ratio = (colMeans(sim) - obs) / apply(sim, 2, sd),
                   ess = coda::effectiveSize(mc), geweke_z = coda::geweke.diag(mc)$z, row.names = NULL)
  cat("\n==== DIAG", mod, " draws", ss, " interval", iv, " minutes", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), "\n"); print(dg, digits = 3)
  saveRDS(dg, file.path(W, paste0("diag_", mod, ".rds"))); quit(save = "no")
}
t0 <- Sys.time()
m <- ergm(as.formula(f), estimate = if (meth == "MPLE") "MPLE" else "MLE", eval.loglik = FALSE, verbose = FALSE, control = ctrl, constraints = ~fixedas(absent = absent))
el <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
sm <- summary(m); cf <- if (!is.null(sm$coefficients)) sm$coefficients else sm$coefs
out <- data.frame(term = rownames(cf), b = cf[, 1], se = cf[, 2], p = cf[, ncol(cf)], row.names = NULL)
diag_ess <- tryCatch({ s <- m$sample; if (is.null(s)) NA else coda::effectiveSize(coda::as.mcmc.list(s)) }, error = function(e) NA)
gw <- tryCatch(coda::geweke.diag(coda::as.mcmc.list(m$sample)[[1]])$z, error = function(e) NA)
attr(out, "info") <- list(model = mod, method = meth, minutes = el, iterations = m$iterations, ess = diag_ess, geweke = gw)
cat("\n==== ", mod, meth, " minutes=", round(el, 1), " iterations=", m$iterations, "\n"); print(out, digits = 4)
if (meth != "MPLE") { cat("ESS:\n"); print(round(diag_ess)); cat("Geweke z:\n"); print(round(gw, 2)) }
saveRDS(out, file.path(W, paste0("res_", mod, "_", meth, if (OLD) "_old" else "", ".rds")))
