# X079 Bell & Long 2016 II "Trade Interdependence and the Use of Force: Do Issues Matter?" (10.1080/03050629.2016.1145119)
# KEY (fixed from abstract before estimating): l.interdep (lagged trade interdependence), in the issue-specific IMI use-of-force models.
#   Evidence: "greater trade interdependence decreases the probability of military conflict over realpolitik issues like territory"
#   Two KEY rows: Model 6 (DV territorial; predicted -) and Model 2 (DV affect_policies; predicted +, 'regime, policies').
# Headline: author do-file: xtgee <dv> l.interdep l.dem_low l.relcap_alt l.alliance_agg l.majpow1 l.majpow2 distance,
#   family(binomial) link(probit) corr(ar1) vce(robust); xtset dyad year (directed dyads). R: geepack::geeglm, AR1 on id=dyad, waves=year.
# CONTROL (DEMOCRACY): l.dem_low (weak-link democracy) in the territorial model.
ROOT <- Sys.getenv("SWEEP_DATA"); OUT <- Sys.getenv("SWEEP_OUT"); CODE <- Sys.getenv("SWEEP_CODE")
stopifnot(nzchar(ROOT), nzchar(OUT), nzchar(CODE))
source(file.path(CODE, "work_p4/b7/common.R"))
suppressMessages(library(geepack))
WO <- file.path(OUT, "work_p4/X079"); dir.create(WO, recursive = TRUE, showWarnings = FALSE)  # port: outputs
W <- file.path(ROOT, "work_p4/X079"); logf <- file.path(WO, "log_X079.txt"); cat("", file = logf)
meta <- list(xid = "X079", cite = "Bell & Long 2016", journal = "II", year = 2016, doi = "10.1080/03050629.2016.1145119")
d <- readRDS(file.path(W, "bl_small.rds"))
rhs <- c("l_interdep","l_dem_low","l_relcap","l_ally","l_maj1","l_maj2","distance")
gfit <- function(r, dv, dat, tag) {
  cache <- file.path(WO, paste0("cache_", tag, ".rds"))
  if (file.exists(cache)) return(readRDS(cache))
  dat <- dat[order(dat$dyad, dat$year), ]
  f <- as.formula(paste(dv, "~", paste(r, collapse = " + ")))
  t0 <- Sys.time()
  m <- geeglm(f, family = binomial(link = "probit"), data = dat, id = dyad, waves = year, corstr = "ar1", std.err = "san.se")
  ct <- coef(summary(m)); ct <- cbind(ct[, 1:2], ct[, 3], 2 * pnorm(-abs(ct[, 1] / ct[, 2])))
  colnames(ct) <- c("Estimate","Std. Error","Wald","Pr(>|z|)")
  cat(tag, "alpha =", unlist(m$geese$alpha), " time:", format(Sys.time() - t0), "\n")
  out <- list(ct = ct, n = length(m$y)); rm(m); gc(); saveRDS(out, cache); out
}
res <- list()
for (dv in c("territorial","affect_policies")) {
  vars <- c(dv, rhs, "dyad", "year")
  d0 <- d[complete.cases(d[, vars]) & d$year <= 2011, ]   # year==9999 junk rows excluded (would also drop under Stata l. operator)
  rep <- gfit(rhs, dv, d0, paste0("repro_", dv)); cat("\n####", dv, "REPRO N =", rep$n, "\n"); print(round(rep$ct, 4))
  cat("\n####", dv, "REPRO\n", file = logf, append = TRUE); capture.output(print(rep$ct), file = logf, append = TRUE)
  s <- common_sample(d0, vars, "ccode1", "ccode2", "year")
  fitter <- function(r, dat) { tag <- paste0(dv, "_", paste(setdiff(r, rhs), collapse = "+")); gfit(r, dv, dat, tag) }
  cat("\n####", dv, "SIX MODELS\n", file = logf, append = TRUE)
  R <- run_six(s, rhs, fitter, c("l_interdep","l_dem_low"), logf)
  # ---- sensitivity hook: runs only under run_sensitivity.R (SWEEP_SENS=1), never in the main sweep ----
  # (headline row = DV territorial; gfit() caches by tag, so each sensitivity fit gets a throwaway tag and its cache file is removed)
  if (Sys.getenv("SWEEP_SENS") == "1" && dv == "territorial") {
    source(file.path(CODE, "hist/sens_helper.R"))
    sens_run(id = "X079", dat = s, y = "territorial", tie = "l_interdep", ctrl = setdiff(rhs, "l_interdep"), fe = NULL,
             fit = function(dat, add) { tg <- basename(tempfile("sens_")); r <- gfit(c(rhs, add), "territorial", dat, tg)
               unlink(file.path(WO, paste0("cache_", tg, ".rds"))); kstat(r$ct, "l_interdep") })
  }
  res[[dv]] <- list(rep = rep, R = R, n_orig = nrow(d0)); rm(s); gc()
}
est <- "GEE probit, AR(1) within directed dyad (waves=year), robust SEs (geepack geeglm)"
ahc <- "N (no peace-years/lagged-conflict term; AR1 working correlation only)"
nt <- function(x, lab) sprintf("%s. Paper tables not accessible (T&F paywall; no log in archive) -> UNVERIFIED. Lags built by dyad-year join (Stata l. semantics). repro N=%d; common N=%d. Outcome = IMI military intervention by stated issue, directed dyad-years 1946+.", lab, x$n_orig, attr(x$R, "n"))
T <- res$territorial; A <- res$affect_policies
write_p4(meta, "l_interdep", "TIE", "KEY", "abstract: 'trade interdependence decreases the probability of military conflict over...territory'",
         "UNVERIFIED", NA, unname(kstat(T$rep$ct, "l_interdep")["b"]), est, T$n_orig, T$R, ahc, nt(T, "Model 6, DV territorial"), fname_term = "l_interdep_territorial")
write_p4(meta, "l_interdep", "TIE", "KEY", "abstract: 'trade interdependence increases...use of force...to affect the regime or policy'",
         "UNVERIFIED", NA, unname(kstat(A$rep$ct, "l_interdep")["b"]), est, A$n_orig, A$R, ahc, nt(A, "Model 2, DV affect_policies"), fname_term = "l_interdep_affect_policies")
write_p4(meta, "l_dem_low", "DEMOCRACY", "CONTROL", "control in all models (weak-link democracy, lagged)",
         "UNVERIFIED", NA, unname(kstat(T$rep$ct, "l_dem_low")["b"]), est, T$n_orig, T$R, ahc, paste(nt(T, "Model 6, DV territorial"), "bD NOT a placebo for this term."))
