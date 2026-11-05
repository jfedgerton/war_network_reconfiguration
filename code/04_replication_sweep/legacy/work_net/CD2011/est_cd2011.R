# CD2011 TERGM (MPLE + 1000 network-year bootstrap, btergm) — one model per call: Rscript est_cd2011.R <model> [R]
ROOT <- Sys.getenv("SWEEP_DATA"); OUT <- Sys.getenv("SWEEP_OUT"); CODE <- Sys.getenv("SWEEP_CODE")
stopifnot(nzchar(ROOT), nzchar(OUT), nzchar(CODE))
W <- file.path(OUT, "work_net/CD2011"); dir.create(W, recursive = TRUE, showWarnings = FALSE)   # res_<model>.rds -> OUT
suppressPackageStartupMessages({library(btergm); library(network)})
args <- commandArgs(TRUE); mod <- args[1]; R <- if (length(args) > 1) as.integer(args[2]) else 1000
B <- readRDS(file.path(ROOT, "work_net/CD2011/cd2011_nets.rds"))
plc <- mod %in% c("pbase", "placebo")
nets <- lapply(B, function(b) { g <- network(if (plc) b$yc else b$y, directed = FALSE)
  set.vertex.attribute(g, "mono15", b$mono15); set.vertex.attribute(g, "mono26", b$mono26); g })
X <- function(v) lapply(B, function(b) b$X[[v]])
minreg <- X("minreg302"); caprat <- X("caprat301"); dist <- X("distance"); seqv <- X("integcorrse")
lagmid <- X("lagmid"); hd15 <- X("hist_dyad_w15"); hd26 <- X("hist_dyad_w26"); cincr <- X("bm_cinc_w15")
base <- "nets ~ edges + edgecov(minreg) + edgecov(caprat) + edgecov(dist) + edgecov(seqv) + kstar(2) + triangle"
f <- switch(mod,
  logit   = "nets ~ edges + edgecov(minreg) + edgecov(caprat) + edgecov(dist) + edgecov(seqv)",
  m0      = base,
  hist    = paste(base, "+ edgecov(hd15) + nodecov('mono15')"),
  m26     = paste(base, "+ edgecov(hd26) + nodecov('mono26')"),
  pbase   = base,
  placebo = paste(base, "+ edgecov(cincr)"),
  ldv0    = paste(base, "+ edgecov(lagmid)"),
  ldvhist = paste(base, "+ edgecov(lagmid) + edgecov(hd15) + nodecov('mono15')"))
set.seed(123)
t0 <- Sys.time()
m <- btergm(as.formula(f), R = R, verbose = FALSE, usefastglm = TRUE)
el <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
bt <- m@boot$t; colnames(bt) <- names(coef(m))
ok <- apply(bt, 1, function(r) all(is.finite(r)))
ci <- t(apply(bt[ok, , drop = FALSE], 2, quantile, c(.025, .975)))
res <- data.frame(term = names(coef(m)), b = coef(m), boot_sd = apply(bt[ok, , drop = FALSE], 2, sd), lo = ci[, 1], hi = ci[, 2],
                  row.names = NULL)
res$p_boot <- sapply(seq_len(ncol(bt)), function(j) { x <- bt[ok, j]; min(1, 2 * min(mean(x <= 0), mean(x >= 0))) })
attr(res, "info") <- list(model = mod, R = R, R_ok = sum(ok), n_dyads = m@nobs, minutes = el, formula = f)
cat("\n==== ", mod, " R=", R, " ok=", sum(ok), " nobs=", m@nobs, " mins=", round(el, 1), "\n"); print(res, digits = 4)
saveRDS(res, file.path(W, paste0("res_", mod, ".rds")))
