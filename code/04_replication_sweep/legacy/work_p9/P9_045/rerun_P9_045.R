# P9_045 Kim 2014 JCR 58(5) "Testing Two Explanations of the Liberal Peace: The Opportunity Cost and Signaling Arguments"
#   (10.1177/0022002713484280). TIE family (challenger's bilateral trade dependence on the target).
# Data + do-files supplied by the PI (UNL Digital Commons supplement; staged read-only): JCR_all.dta, JCR_onlyMID.dta,
#   "Table 1 & 2.do" -> kim_tables.do, "Figure 1s.do" -> kim_figures.do, "Supporting appendix.pdf" -> kim_appendix.pdf. No log.
# Directed dyad-years 1950-2001. HEADLINE: Table 1 Model 1 (kim_tables.do, first pass of the foreach loop):
#   heckprob recip Ldepkg1 $control1, select(initiation = Ldepkg1 $control2) cluster(ddyadid)
#   control1 = demodum1 demodum2 capshare1 majmin minmaj majmaj ally sglo lndist contig territory policy regimegovernment others
#   control2 = demodum1 demodum2 capshare1 majmin minmaj majmaj ally sglo lndist contig peaceyrs _spline1 _spline2 _spline3
#   KEY TERM out:Ldepkg1 (challenger's lagged bilateral trade / GDP) in the RECIPROCATION (outcome) equation.
#   Published -0.648* (0.343), N 907,454, 1,418 uncensored; significant at the 10% level only (article text + star legend).
# DV timing: initiation and reciprocation are in year t, economic variables are lagged (L prefix) -> history attached at year t.
# YEAR RECONSTRUCTION: JCR_all.dta carries no year column. Rows are stored in contiguous ccode2 blocks; within a block each year
#   holds each directed dyad at most once, so a new year starts when a dyad repeats. Block-local year indices are anchored to calendar
#   years by the block's last year = min(2001, the target's last COW system year) and validated against the dated MID rows in
#   JCR_onlyMID.dta (cwkeynum + dyad -> year). Fails loudly if < 97% of anchor rows agree.
# Stata marks out outcome-equation regressors only on selected rows (issue dummies exist only for MIDs); censored rows contribute
#   Phi(-Zg) only, so their X is unused and set to 0 (as P9_057).
# Estimator: heckprob by ML (work_p8/heckprob_b3.R), cluster-robust SE on ddyadid (Stata G/(G-1)), all six steps; history and
#   benchmark terms enter both equations (as P8_003, P9_057).
CODE <- Sys.getenv("SWEEP_CODE"); stopifnot(nzchar(CODE))
source(file.path(CODE, "work_p8/common_b3.R"))
source(file.path(CODE, "work_p8/heckprob_b3.R"))
WD <- file.path(ROOT, "work_p9/P9_045"); W <- file.path(OUT, "work_p9/P9_045"); dir.create(W, recursive = TRUE, showWarnings = FALSE)
sink(file.path(W, "log_P9_045.txt"), split = TRUE)
c1 <- c("demodum1", "demodum2", "capshare1", "majmin", "minmaj", "majmaj", "ally", "sglo", "lndist", "contig")
outX <- c("Ldepkg1", c1, "territory", "policy", "regimegovernment", "others")
selZ <- c("Ldepkg1", c1, "peaceyrs", "spline1", "spline2", "spline3")
d <- haven::read_dta(file.path(WD, "JCR_all.dta"), col_select = all_of(c("ccode1", "ccode2", "ddyadid", "cwkeynum", "initiation", "recip",
  setdiff(outX, character(0)), "peaceyrs", "_spline1", "_spline2", "_spline3")))
d <- as.data.frame(haven::zap_labels(d)); names(d) <- sub("^_spline", "spline", names(d))
for (v in names(d)) d[[v]] <- as.numeric(d[[v]])
cat("rows", nrow(d), "\n"); stopifnot(nrow(d) == 1088588, all(d$ddyadid == 1000 * d$ccode1 + d$ccode2), all(d$initiation %in% 0:1))

# ---- year reconstruction ----
n <- nrow(d); key <- d$ddyadid
r2 <- rle(d$ccode2); stopifnot(length(r2$lengths) == length(unique(d$ccode2)))   # contiguous target blocks
seg <- rep(seq_along(r2$lengths), r2$lengths)
blk <- integer(n); cur <- 0L; last_seg <- 0L; seen <- NULL
for (i in seq_len(n)) {
  if (seg[i] != last_seg) { cur <- 1L; last_seg <- seg[i]; seen <- new.env(hash = TRUE, size = 4096L) }
  k <- as.character(key[i]); if (!is.null(seen[[k]])) { cur <- cur + 1L; seen <- new.env(hash = TRUE, size = 4096L) }
  assign(k, TRUE, envir = seen); blk[i] <- cur
}
nblk <- tapply(blk, seg, max)
cs <- peacesciencer::cow_states
endy <- sapply(r2$values, function(cc) min(2001, max(cs$endyear[cs$ccode == cc])))
stopifnot(all(is.finite(endy)))
d$year <- (endy - nblk)[seg] + blk
stopifnot(min(d$year) >= 1950, max(d$year) <= 2001, !anyDuplicated(d[, c("ddyadid", "year")]))
m <- as.data.frame(haven::zap_labels(haven::read_dta(file.path(WD, "JCR_onlyMID.dta"), col_select = c("ccode1", "ccode2", "year", "cwkeynum"))))
an <- merge(m[m$cwkeynum > 0, ], d[d$cwkeynum > 0, c("ccode1", "ccode2", "cwkeynum", "year")], by = c("ccode1", "ccode2", "cwkeynum"),
            suffixes = c("_mid", "_rec"))
agree <- mean(an$year_mid == an$year_rec)
cat("year reconstruction: anchor rows", nrow(an), " agreement", round(agree, 4), "\n"); stopifnot(nrow(an) > 1000, agree >= 0.97)
o <- order(d$ddyadid, d$year); same <- diff(d$ddyadid[o]) == 0
cat("consecutive within-dyad years:", round(mean(diff(d$year[o])[same] == 1), 4),
    " peaceyrs +1 or reset:", round(mean((diff(d$peaceyrs[o]) == 1 | d$peaceyrs[o][-1] == 0)[same], na.rm = TRUE), 4), "\n"); rm(o, same, m, an, blk, seg)

# ---- heckprob setup ----
sel <- d$initiation == 1
stopifnot(all(d$recip[sel] %in% c(0, 1, NA)))
oc <- c("territory", "policy", "regimegovernment", "others")
cat("outcome-only vars missing on censored rows:", sapply(oc, function(v) mean(is.na(d[[v]][!sel]))), "\n")
d$recip_chk <- ifelse(sel, d$recip, 0)
for (v in oc) d[[v]] <- ifelse(sel, d[[v]], 0)
# Stata heckprob sample: selection regressors complete on all rows; outcome regressors + DV complete on selected rows
okrow <- complete.cases(d[, selZ]) & (!sel | complete.cases(d[, c("recip", outX)]))
d$use <- okrow
cat("estimation rows", sum(okrow), " selected", sum(okrow & sel), "\n")
KEY <- "out:Ldepkg1"
fit_hp <- function(extra, data) {
  X <- cbind(as.matrix(data[, c(outX, extra), drop = FALSE]), `_cons` = 1)
  Z <- cbind(as.matrix(data[, c(selZ, extra), drop = FALSE]), `_cons` = 1)
  f <- fit_heckprob(X, Z, data$recip_chk, data$initiation, cl = data$ddyadid)
  z <- f$coef / f$se
  cat("---- heckprob | extra:", paste(extra, collapse = " + "), "| N =", f$n, "selected =", f$n_sel, "| LL =", round(f$ll, 3), "| rho =", round(f$rho, 4), "\n")
  print(round(cbind(b = f$coef, se = f$se, z = z), 6))
  out <- list(Ldepkg1 = c(b = f$coef[[KEY]], se = f$se[[KEY]], p = 2 * pnorm(-abs(z[[KEY]]))))
  attr(out, "n") <- f$n; attr(out, "athrho") <- f$coef[["athrho"]]; out
}
# gradient check of the analytic heckprob scores on synthetic data (fail loudly; same check as P8_003 / P9_057)
set.seed(123); nn <- 4000; Xs0 <- cbind(a = rnorm(nn), c = 1); Zs0 <- cbind(a = rnorm(nn), b = rnorm(nn), c = 1)
ee <- MASS::mvrnorm(nn, c(0, 0), matrix(c(1, .4, .4, 1), 2)); ss <- as.numeric(Zs0 %*% c(.5, -.3, -.2) + ee[, 2] > 0)
yy <- as.numeric(Xs0 %*% c(.7, .1) + ee[, 1] > 0); yy[ss == 0] <- 0; th <- c(.6, .05, .45, -.25, -.15, .3)
ga <- colSums(hp_parts(th, Xs0, Zs0, yy, ss)$S); gn <- numDeriv::grad(function(t) sum(hp_parts(t, Xs0, Zs0, yy, ss)$ll), th)
cat("synthetic gradient check max abs diff:", max(abs(ga - gn)), "\n"); stopifnot(max(abs(ga - gn)) < 1e-5)
vars <- unique(c("initiation", "recip_chk", outX, selZ, "ddyadid"))
du <- d[d$use, ]; stopifnot(all(complete.cases(du[, vars])))
full <- fit_hp(character(0), du); rm(du); invisible(gc())
pub <- c(b = -0.648, se = 0.343); pub_n <- 907454
cat("\nFULL-SAMPLE:", KEY, "b", full$Ldepkg1[["b"]], "se", full$Ldepkg1[["se"]], "N", attr(full, "n"), "| published", pub, pub_n, "(1,418 uncensored)\n")
ds <- common_sample(d[d$use, ], vars, "ccode1", "ccode2", "year"); rm(d); invisible(gc())
cat("selected rows in common sample:", sum(ds$initiation), "\n")
fitter <- function(rhs, data) fit_hp(setdiff(rhs, vars), data)
M <- run_six(vars, ds, fitter)
saveRDS(list(M = M, full = full, n_used = nrow(ds)), file.path(W, "M_P9_045.rds"))

meta <- list(pid = "P9_045", cite = "Kim 2014", journal = "JCR", year = 2014, doi = "10.1177/0022002713484280", family = "TIE")
nt <- sprintf(paste0("Directed dyad-years 1950-2001 (author's data + do-file; year reconstructed from row order, see script header). ",
  "Published Table 1 Model 1 -0.648* (0.343), N 907,454 (1,418 uncensored); article treats p<.10 as significant. ",
  "Full-sample repro b=%.6f se=%.6f N=%d. Common 1822-2011 sample N=%d (%d selected). Cluster-robust SE on ddyadid. History terms in both equations."),
  full$Ldepkg1[["b"]], full$Ldepkg1[["se"]], attr(full, "n"), nrow(ds), sum(ds$initiation))
write_p8(meta, "Ldepkg1", "abstract: 'as challengers are more dependent on bilateral trade, targets are less likely to reciprocate disputes'",
  "Table 1 Model 1 (heckprob; outcome eq. reciprocation, selection eq. initiation; do-file)",
  "probit with sample selection (heckprob), SE clustered on directed dyad, peace years + splines in selection eq.",
  "heckprob ML (R port, analytic scores), cluster-robust SE ddyadid, all six steps",
  pub[["b"]], pub[["se"]], pub_n, FALSE, FALSE, M, "N (peace years + splines in selection eq. only)", nt)
dir.create(file.path(OUT, "results_p9"), recursive = TRUE, showWarnings = FALSE)
src <- file.path(OUT, "results_p8", "P9_045__Ldepkg1.csv"); stopifnot(file.exists(src), file.rename(src, file.path(OUT, "results_p9", basename(src))))
sink()
