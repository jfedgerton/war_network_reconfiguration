# FRITZ2025: signed TERGM (ergm.sign::tsergm), Model 1 of Table 1 (+ history variants).
# Usage: Rscript fritz_run.R <model: m0|hist|m26|placebo> <method: MPLE|Stepping> [ncores]
ROOT <- Sys.getenv("SWEEP_DATA"); OUT <- Sys.getenv("SWEEP_OUT"); CODE <- Sys.getenv("SWEEP_CODE")
stopifnot(nzchar(ROOT), nzchar(OUT), nzchar(CODE))
W <- file.path(OUT, "work_net/FRITZ2025"); dir.create(W, recursive = TRUE, showWarnings = FALSE); setwd(W)   # outputs (res/fit/mple rds) -> OUT
suppressPackageStartupMessages({library(ergm.sign); library(dplyr); library(parallel)})
source(file.path(CODE, "hist/add_history.R"))
a <- commandArgs(TRUE); mod <- a[1]; meth <- a[2]; ncores <- if (length(a) > 2) as.integer(a[3]) else 1
tmp_data <- readRDS(file.path(ROOT, "work_net/FRITZ2025/src/application_dca/models/tmp_data.RDS"))
tmp_network <- tmp_data$target_matrix; yrs <- 2001:2010
stopifnot(length(tmp_network) == 10)
alpha <- matrix(log(2), 1, 1); alpha <- lapply(seq_along(tmp_network), function(X) alpha)
repetition_pos <- tmp_data$cov_repetition_pos; repetition_neg <- tmp_data$cov_repetition_neg
polity <- tmp_data$cov_polity; gdp <- tmp_data$cov_gdp; distance <- tmp_data$cov_distance; cinc <- tmp_data$cov_cinc
# ---- conflict history matrices (same node order as target matrices) ----
H <- readRDS(Sys.getenv("SWEEP_HIST"))
mk_hist <- function(k) {
  t <- yrs[k]; cc <- as.numeric(rownames(tmp_network[[k]])); n <- length(cc)
  g <- expand.grid(i = 1:n, j = 1:n); g <- g[g$i < g$j, ]
  dd <- data.frame(c1 = cc[g$i], c2 = cc[g$j], year = t); dd <- add_history(dd, "c1", "c2", "year")
  s <- H$sth[H$sth$year == t, ]; m15 <- log1p(coalesce(s$s15[match(cc, s$ccode)], 0)); m26 <- log1p(coalesce(s$s26[match(cc, s$ccode)], 0))
  M <- function(v) { x <- matrix(0, n, n, dimnames = dimnames(tmp_network[[k]])); x[cbind(g$i, g$j)] <- v; x[cbind(g$j, g$i)] <- v; x }
  cr <- dd$bm_cinc_w15; cr[is.na(cr)] <- mean(cr, na.rm = TRUE)
  list(d15 = M(dd$hist_dyad_w15), d26 = M(dd$hist_dyad_w26), m15 = M(m15[g$i] + m15[g$j]), m26 = M(m26[g$i] + m26[g$j]),
       cr = M(cr), n_cr_na = sum(is.na(dd$bm_cinc_w15)))
}
HL <- lapply(seq_along(yrs), mk_hist)
hd15 <- lapply(HL, `[[`, "d15"); hd26 <- lapply(HL, `[[`, "d26"); hm15 <- lapply(HL, `[[`, "m15"); hm26 <- lapply(HL, `[[`, "m26"); cincr <- lapply(HL, `[[`, "cr")
cat("bm_cinc NA dyads (mean-imputed):", sum(sapply(HL, `[[`, "n_cr_na")), "\n")
cat("MID ties with hist_dyad_w15>0:", sum(sapply(seq_along(yrs), function(k) sum(tmp_network[[k]] == -1 & hd15[[k]] > 0))) / 2,
    "of", sum(sapply(tmp_network, function(m) sum(m == -1))) / 2, "\n")
base <- "tmp_network ~ edges_pos + isolates_pos + isolates_neg + cov_dyad_pos(data = repetition_pos) + cov_dyad_neg(data = repetition_neg) +
  cov_dyad_pos(data = polity) + cov_dyad_neg(data = polity) + cov_dyad_pos(data = cinc) + cov_dyad_neg(data = cinc) + cov_dyad(data = distance) +
  cov_dyad_pos(data = gdp) + cov_dyad_neg(data = gdp) + gwese_pos(data = alpha) + gwese_neg(data = alpha) + gwesf_pos(data = alpha) +
  gwesf_neg(data = alpha) + gwdegree_pos(data = alpha) + gwdegree_neg(data = alpha) + edges_neg"
f <- switch(mod, m0 = base,
  hist = paste(base, "+ cov_dyad_neg(data = hd15) + cov_dyad_neg(data = hm15)"),
  m26 = paste(base, "+ cov_dyad_neg(data = hd26) + cov_dyad_neg(data = hm26)"),
  placebo = paste(base, "+ cov_dyad_neg(data = cincr)"))
clust <- NULL
if (ncores > 1) { clust <- makeCluster(ncores); clusterCall(clust, function() { library(ergm.sign); library(Rcpp); library(RcppArmadillo) }) }
# authors' control settings (application_dca.R); eval_likelihood=F (AIC path sampling only; does not affect estimates)
ctrl <- control.sergm(method_est = meth, tol = 0.2, method_var = "Fisher", max_it = 30, eval_likelihood = FALSE, cluster = clust,
  sampler_var = sampler.sergm(number_networks = 3000, init_empty = TRUE, n_proposals_burn_in = 100000, n_proposals = 10000, mh = FALSE, seed = 13),
  sampler_path_sampling = sampler.sergm(number_networks = 3000, init_empty = TRUE, n_proposals_burn_in = 100000, n_proposals = 1000, mh = FALSE, seed = 13),
  sampler_est = sampler.sergm(number_networks = 2000, init_empty = FALSE, n_proposals_burn_in = 1000 * 10, n_proposals = 1000, mh = FALSE, seed = 13),
  Stepping_number_grids = 1000)
set.seed(13)
t0 <- Sys.time()
m <- tsergm(formula = as.formula(f), control = ctrl)
el <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
if (!is.null(clust)) stopCluster(clust)
b <- m$coefficients; se <- sqrt(diag(m$var))
out <- data.frame(term = names(b), b = b, se = se, p = 2 * pnorm(-abs(b / se)), row.names = NULL)
ess <- tryCatch(coda::effectiveSize(m$mcmc_chain), error = function(e) rep(NA, length(b)))
gw <- tryCatch(coda::geweke.diag(m$mcmc_chain)$z, error = function(e) rep(NA, length(b)))
out$ess <- as.numeric(ess); out$geweke_z <- as.numeric(gw)
attr(out, "info") <- list(model = mod, method = meth, minutes = el, n_iter = length(m$coefficients_per_iteration))
cat("\n==== ", mod, meth, " minutes=", round(el, 1), " iterations=", length(m$coefficients_per_iteration), "\n"); print(out, digits = 4)
saveRDS(out, file.path(W, paste0("res_", mod, "_", meth, ".rds"))); saveRDS(m, file.path(W, paste0("fit_", mod, "_", meth, ".rds")))
