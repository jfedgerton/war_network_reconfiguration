# FRITZ2025: signed TERGM (ergm.sign::tsergm), Model 1 of Table 1 (+ history variants).
# MPLE point estimates (no variance) for all four specifications: fast cross-check / fallback for m26 & placebo
ROOT <- Sys.getenv("SWEEP_DATA"); OUT <- Sys.getenv("SWEEP_OUT"); CODE <- Sys.getenv("SWEEP_CODE")
stopifnot(nzchar(ROOT), nzchar(OUT), nzchar(CODE))
W <- file.path(OUT, "work_net/FRITZ2025"); dir.create(W, recursive = TRUE, showWarnings = FALSE); setwd(W)   # outputs (res/fit/mple rds) -> OUT
suppressPackageStartupMessages({library(ergm.sign); library(dplyr); library(parallel)})
source(file.path(CODE, "hist/add_history.R"))

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

res <- list()
for (mod in c("m0","hist","m26","placebo")) {
  f <- switch(mod, m0 = base,
    hist = paste(base, "+ cov_dyad_neg(data = hd15) + cov_dyad_neg(data = hm15)"),
    m26 = paste(base, "+ cov_dyad_neg(data = hd26) + cov_dyad_neg(data = hm26)"),
    placebo = paste(base, "+ cov_dyad_neg(data = cincr)"))
  pp <- ergm.sign:::t_formula_preprocess(as.formula(f))
  mm <- ergm.sign:::t_mple_est(networks = pp$network, n_actors = nrow(pp$network[[1]]), mod_return = TRUE, terms = pp$term_names,
                               data_lists = pp$data_list, type_list = pp$type_list, coef_names = pp$coef_names)
  b <- mm$coefficients; names(b) <- pp$coef_names; res[[mod]] <- b
  cat("\n== MPLE", mod, "\n"); print(round(b, 4))
}
saveRDS(res, file.path(W, "mple_points.rds"))
