# KEY rerun 100: Peterson & Quackenbush 2010 II. Key claim (abstract): "imposed settlements foster a pacifying effect of trade".
# Headline Model 1 (all settlement types) has no term for that claim; the authors test it in Table 1 Model 2 (imposed==1 spells):
#   stcox imposeedep imposerdep imp_interdep bothdem2 bothdem2xdurmo imp_relpow imp_relpowdif lndist strival if imposed==1, cluster(dyadid) nohr
# Key term pre-committed = imposeedep (imposee/recipient trade dependence; first trade term, Figure 1). imposerdep, imp_interdep and a
# joint Wald test of all three trade terms are logged/saved too. Same estimator/SE and six-model scheme as rerun_100.R.
ROOT <- Sys.getenv("SWEEP_DATA"); OUT <- Sys.getenv("SWEEP_OUT"); CODE <- Sys.getenv("SWEEP_CODE")
stopifnot(nzchar(ROOT), nzchar(OUT), nzchar(CODE))
source(file.path(CODE, "work/rerun_common_g7.R"))
suppressPackageStartupMessages({library(haven); library(survival)})
W <- file.path(ROOT, "work/100")
WO <- file.path(OUT, "work/100"); dir.create(WO, recursive = TRUE, showWarnings = FALSE)
sink(file.path(WO, "log_key_100.txt"), split = TRUE)
d <- as.data.frame(zap_labels(read_dta(file.path(W, "allmids.dta"))))
stopifnot(nrow(d) == 30765, all(d$`_st` == 1))
names(d)[names(d) == "_t0"] <- "t0"; names(d)[names(d) == "_t"] <- "t1"; names(d)[names(d) == "_d"] <- "dd"
d <- d[d$imposed == 1, ]
cat("imposed==1 rows:", nrow(d), "\n")
RHS <- "imposeedep + imposerdep + imp_interdep + bothdem2 + bothdem2xdurmo + imp_relpow + imp_relpowdif + lndist + strival"
MV <- c("t0","t1","dd","imposeedep","imposerdep","imp_interdep","bothdem2","bothdem2xdurmo","imp_relpow","imp_relpowdif","lndist","strival","dyadid")
TIE <- "imposeedep"; TR <- c("imposeedep","imposerdep","imp_interdep")
fitter <- function(add, dat) {
  dat <- dat[complete.cases(dat[, MV]), ]
  if (add != "") dat <- dat[complete.cases(dat[, trimws(strsplit(sub("^ \\+ ", "", add), "\\+")[[1]])]), ]
  f <- as.formula(paste0("Surv(t0, t1, dd) ~ ", RHS, add))
  m <- coxph(f, data = dat, ties = "breslow", cluster = dyadid)
  G <- length(unique(dat$dyadid)); V <- m$var * G / (G - 1)
  b <- coef(m); se <- sqrt(diag(V)); names(se) <- names(b); z <- b / se; p <- 2 * pnorm(-abs(z))
  cat("\n---- ", deparse(f, width.cutoff = 500), "| N =", m$n, " events =", m$nevent, " clusters =", G, "\n")
  print(signif(cbind(b = b, se = se, z = z, p = p), 5))
  k <- match(TR, names(b)); stopifnot(!anyNA(k))
  W3 <- as.numeric(t(b[k]) %*% solve(V[k, k]) %*% b[k]); pj <- pchisq(W3, 3, lower.tail = FALSE)
  cat(sprintf("joint Wald (3 trade terms): chi2=%.3f p=%.4g\n", W3, pj))
  list(b = unname(b[TIE]), se = unname(se[TIE]), p = unname(p[TIE]), n = m$n,
       b_er = unname(b["imposerdep"]), p_er = unname(p["imposerdep"]), b_int = unname(b["imp_interdep"]), p_int = unname(p["imp_interdep"]), p_joint = pj)
}
cat("\n==== Reproduction: Table 1 Model 2 ====\n")
rp <- fitter("", d); n_orig <- rp$n
d <- add_history(d, "ccode1", "ccode2", "year")
d <- restrict_sample(d, MV)
R <- run_six(fitter, d)
saveRDS(list(rp = rp, n_orig = n_orig, R = R, RS_INFO = RS_INFO), file.path(WO, "key_100.rds"))
print(sapply(R, function(x) unlist(x)))
# port (extraction only, no estimation): CSV rows for the two master terms from the fits above.
#  imposeedep: same numbers as key_100.rds (copied into results_key/100.csv by write_key_batch.py).
#  imp_interdep: master row was hand-copied by b1work/v/copyrows.py from log_key_100.txt (fitter prints signif(., 5));
#    the pipeline reports full precision instead (differences < 1e-5 relative).
pc <- function(b, b0) 100 * (b - b0) / abs(b0)
s5 <- function(k, v) R[[k]][[v]]
ex <- rbind(
  data.frame(id = "100", key_term = "imposeedep", b0 = R$m0$b, p0 = R$m0$p, b_dm = R$mDM$b, p_dm = R$mDM$p, b_26 = R$m26$b, p_26 = R$m26$p,
             pct_bC = pc(R$bC$b, R$m0$b), pct_bD = pc(R$bD$b, R$m0$b), n_used = R$m0$n),
  data.frame(id = "100", key_term = "imp_interdep", b0 = s5("m0", "b_int"), p0 = s5("m0", "p_int"), b_dm = s5("mDM", "b_int"), p_dm = s5("mDM", "p_int"),
             b_26 = s5("m26", "b_int"), p_26 = s5("m26", "p_int"), pct_bC = pc(s5("bC", "b_int"), s5("m0", "b_int")),
             pct_bD = pc(s5("bD", "b_int"), s5("m0", "b_int")), n_used = R$m0$n))
stopifnot(all(is.finite(unlist(ex[, -(1:2)]))))
write.csv(ex, file.path(WO, "key_100_extract.csv"), row.names = FALSE)
sink()
