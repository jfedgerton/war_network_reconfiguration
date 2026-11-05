# P8_003 Goldsmith 2013 DPE "International Trade and the Onset and Escalation of Interstate Conflict" (10.1080/10242694.2013.763637)
# PDF NOT ACCESSIBLE (T&F paywall; SSRN/ResearchGate 429; the user's device that stages sweep PDFs is disconnected) -> spec reconstructed.
# Data: the author's companion archive for Goldsmith 2013 IRAP "Different in Asia?" (Harvard Dataverse doi:10.7910/DVN/29717;
#   same author/year, same variables: fcwmidinit, fcwmid251, lnTradeVolumeGledd, GDPShareGledLow, controls). Its do-file estimates
#   heckprob (probit with selection): escalation fcwmid251 | onset fcwmidinit, excluding joiners (fcwjomid ~= 1), OIM SEs.
# Reconstructed DPE headline = the IRAP Model-1 equation set on the GLOBAL sample (no regional restriction), 1950+ (Gleditsch trade era;
#   IRAP Models 4/7 use year > 1949). Abstract: "two aspects of trade - volume and interdependence - ... two-stage process involving onset
#   and escalation"; El-Anis 2018 summarizes the DPE findings as interdependence inhibits onset, volume raises onset, interdependence has no
#   escalation effect, volume inhibits escalation.
# KEY TERMS (a priori): sel:GDPShareGledLow, sel:lnTradeVolumeGledd (onset), out:lnTradeVolumeGledd, out:GDPShareGledLow (escalation).
# Estimator: heckprob by ML (work_p8/heckprob_b3.R: analytic scores, Genz bivariate normal CDF), OIM SEs, all six steps; history/benchmark
#   terms enter both equations.
CODE <- Sys.getenv("SWEEP_CODE"); stopifnot(nzchar(CODE))
source(file.path(CODE, "work_p8/common_b3.R"))
source(file.path(CODE, "work_p8/heckprob_b3.R"))
WD <- file.path(ROOT, "work_p8/P8_003"); W <- file.path(OUT, "work_p8/P8_003"); dir.create(W, recursive = TRUE, showWarnings = FALSE); sink(file.path(W, "log_P8_003.txt"), split = TRUE)
outX <- c("lnTradeVolumeGledd", "GDPShareGledLow", "lncap_12", "jntautoc")
selZ <- c("lnTradeVolumeGledd", "GDPShareGledLow", "parity", "lncap_12", "s_wt_glo", "jntdemord", "jntautoc", "contig24", "lndistance",
          "cwpceyrs", "cwpceyrs2", "cwpceyrs3")
d <- haven::read_dta(file.path(WD, "irap.dta"), col_select = all_of(unique(c("ccode1", "ccode2", "year", "fcwmidinit", "fcwmid251", "fcwjomid", outX, selZ))))
d <- as.data.frame(haven::zap_labels(d)); cat("rows", nrow(d), "\n")
d <- d[(is.na(d$fcwjomid) | d$fcwjomid != 1) & d$year > 1949, ]   # Stata fcwjomid~=1 keeps missing fcwjomid
KEYS <- c("sel:GDPShareGledLow", "sel:lnTradeVolumeGledd", "out:lnTradeVolumeGledd", "out:GDPShareGledLow")
fit_hp <- function(extra, data, start = NULL) {
  X <- cbind(as.matrix(data[, c(outX, extra), drop = FALSE]), `_cons` = 1)
  Z <- cbind(as.matrix(data[, c(selZ, extra), drop = FALSE]), `_cons` = 1)
  s <- data$fcwmidinit; y <- ifelse(s == 1, data$fcwmid251, 0)
  f <- fit_heckprob(X, Z, y, s, start = start)
  z <- f$coef / f$se
  cat("---- heckprob | extra:", paste(extra, collapse = " + "), "| N =", f$n, "selected =", f$n_sel, "| LL =", round(f$ll, 3), "| rho =", round(f$rho, 4), "\n")
  print(round(cbind(b = f$coef, se = f$se, z = z), 5))
  out <- lapply(setNames(KEYS, KEYS), function(k) c(b = f$coef[[k]], se = f$se[[k]], p = 2 * pnorm(-abs(z[[k]]))))
  attr(out, "n") <- f$n; out
}
vars <- unique(c("fcwmidinit", outX, selZ))
# gradient check of the analytic heckprob scores on synthetic data (bvn_lower also checked vs mvtnorm in dbg.out; fail loudly if wrong)
set.seed(123); nn <- 4000; Xs0 <- cbind(a = rnorm(nn), c = 1); Zs0 <- cbind(a = rnorm(nn), b = rnorm(nn), c = 1)
ee <- MASS::mvrnorm(nn, c(0, 0), matrix(c(1, .4, .4, 1), 2)); ss <- as.numeric(Zs0 %*% c(.5, -.3, -.2) + ee[, 2] > 0)
yy <- as.numeric(Xs0 %*% c(.7, .1) + ee[, 1] > 0); yy[ss == 0] <- 0; th <- c(.6, .05, .45, -.25, -.15, .3)
ga <- colSums(hp_parts(th, Xs0, Zs0, yy, ss)$S); gn <- numDeriv::grad(function(t) sum(hp_parts(t, Xs0, Zs0, yy, ss)$ll), th)
cat("synthetic gradient check max abs diff:", max(abs(ga - gn)), "\n"); stopifnot(max(abs(ga - gn)) < 1e-5)
full <- d[complete.cases(d[, vars]) & (d$fcwmidinit == 0 | !is.na(d$fcwmid251)), ]
r_full <- fit_hp(character(0), full); n_orig <- attr(r_full, "n")
d$fcwmid251_chk <- ifelse(d$fcwmidinit == 1, d$fcwmid251, 0)
ds <- common_sample(d, c(vars, "fcwmid251_chk"), "ccode1", "ccode2", "year"); rm(d, full); invisible(gc())
fitter <- function(rhs, data) fit_hp(setdiff(rhs, vars), data)
M <- run_six(vars, ds, fitter)
saveRDS(list(M = M, r_full = r_full, n_orig = n_orig, n_used = nrow(ds)), file.path(W, "M_P8_003.rds"))
meta <- list(pid = "P8_003", cite = "Goldsmith 2013", journal = "DPE", year = 2013, doi = "10.1080/10242694.2013.763637", family = "TIE")
nt <- sprintf(paste0("PDF inaccessible -> headline spec reconstructed from the author's companion IRAP 2013 archive (same variables; heckprob escalation fcwmid251 | onset fcwmidinit; joiners excluded) ",
  "on the global sample 1950+; no published numbers available -> UNVERIFIED. Full-sample N=%d; common 1822-2011 N=%d. OIM SEs (IRAP do-file has no vce option). ",
  "History terms in both equations."), n_orig, nrow(ds))
for (k in KEYS) write_p8(meta, gsub(":", "_", k), "abstract: 'two aspects of trade - volume and interdependence - ... onset and escalation'",
  paste("reconstructed: heckprob (IRAP Model 1 equations, global sample),", k), "probit with sample selection (heckprob), onset = selection eq, escalation = outcome eq",
  "heckprob ML (R port, analytic scores), OIM SE, all six steps", NA, NA, NA, FALSE, FALSE,
  lapply(M, function(m) { o <- setNames(list(m[[k]]), gsub(":", "_", k)); attr(o, "n") <- attr(m, "n"); o }), "N (peace years + cubic)", nt, verdict = "UNVERIFIED")
sink()
