# Aggregate WSC2007, MINHAS2022 (RS, WEEKS, GIBLER), WH2011 into results_net/<tag>__<term>.csv (CD2011 column layout)
suppressPackageStartupMessages(library(coda))
OUT <- Sys.getenv("SWEEP_OUT"); stopifnot(nzchar(OUT))
B <- file.path(OUT, "work_net"); O <- file.path(OUT, "results_net"); dir.create(O, recursive = TRUE, showWarnings = FALSE)
# Port: optional args select which tags to aggregate (default all), so each study can be rebuilt on its own:
#   Rscript agg_mine.R [WSC2007] [MINHAS2022_RS] [MINHAS2022_WEEKS] [MINHAS2022_GIBLER] [WH2011]
SEL <- commandArgs(TRUE); if (!length(SEL)) SEL <- c("WSC2007", "MINHAS2022_RS", "MINHAS2022_WEEKS", "MINHAS2022_GIBLER", "WH2011")
stopifnot(all(SEL %in% c("WSC2007", "MINHAS2022_RS", "MINHAS2022_WEEKS", "MINHAS2022_GIBLER", "WH2011")))
summ <- function(x) { m <- mcmc(x); c(b = mean(x), se = sd(x), lo = unname(quantile(x, .025)), hi = unname(quantile(x, .975)),
  p = 2 * min(mean(x > 0), mean(x < 0)), ess = unname(effectiveSize(m)), gew = unname(geweke.diag(m)$z)) }
excl <- function(s) unname(s["lo"] > 0 | s["hi"] < 0)
pct <- function(s, s0) if (is.null(s)) NA else unname(100 * (s["b"] - s0["b"]) / abs(s0["b"]))
# REPRODUCED: within 10% (or rounding tolerance) and 95% CI significance agrees; CLOSE: significance agrees and same sign (or a null whose gap to pub < 1 posterior sd)
verdict <- function(s0, pub, pub_excl, tol_abs = 0) {
  if (is.na(pub)) return(NA)
  same_sig <- excl(s0) == pub_excl
  if (abs(s0["b"] - pub) <= max(0.10 * abs(pub), tol_abs) && same_sig) "REPRODUCED" else if (same_sig && (sign(s0["b"]) == sign(pub) || abs(s0["b"] - pub) < s0["se"])) "CLOSE" else "NOT_REPRODUCED" }
mkrow <- function(tag, study, key, term, pub, pub_ci, s0, sh, s26, spl, conv, verd, notes) {
  g <- function(s, k) if (is.null(s)) NA else unname(s[k])
  d <- data.frame(tag = tag, study = study, key_effect = key, model_term = term, pub_b = pub, pub_ci = pub_ci,
    b0 = g(s0, "b"), se0 = g(s0, "se"), p0 = g(s0, "p"), ci0_excl0 = if (is.null(s0)) NA else excl(s0),
    b_hist = g(sh, "b"), se_hist = g(sh, "se"), p_hist = g(sh, "p"), ci_hist_excl0 = if (is.null(sh)) NA else excl(sh),
    pct_hist = if (is.null(sh)) NA else pct(sh, s0), b_26 = g(s26, "b"), pct_26 = if (is.null(s26)) NA else pct(s26, s0),
    b_placebo = g(spl, "b"), pct_placebo = if (is.null(spl)) NA else pct(spl, s0), placebo_base_b0 = g(s0, "b"),
    convergence = conv, verdict = verd, notes = notes, row.names = NULL)
  write.csv(d, file.path(O, paste0(tag, "__", gsub("[^A-Za-z0-9_.]", "_", key), ".csv")), row.names = FALSE); d }
essfmt <- function(lst) paste(sprintf("%s ESS=%.0f Geweke z=%.2f", names(lst), sapply(lst, function(s) s["ess"]), sapply(lst, function(s) s["gew"])), collapse = "; ")
out <- list()

## ---------- WSC2007: gbme (Hoff 2003) binomial k=3 directed, 2000 cross-section (Table 2)
if ("WSC2007" %in% SEL) {
rd <- function(m) { f <- file.path(B, "WSC2007", paste0("out_", m), "OUT"); if (!file.exists(f) || !file.exists(file.path(B, "WSC2007", paste0("log_wsc_", m, ".txt.ok")))) return(NULL)
  o <- read.table(f, header = TRUE); stopifnot(max(o$scan) == 25000); o[o$scan > 5000, ] }
o0 <- rd("m0"); oh <- rd("hist"); o26 <- rd("m26"); opl <- rd("placebo"); stopifnot(!is.null(o0), !is.null(oh))
pub <- c(bd1 = -0.02, bd2 = -0.18, bd3 = 0.02); pubci <- c(bd1 = "[-0.08, 0.05]", bd2 = "[-0.40, 0.01]", bd3 = "[-0.22, 0.26]")
lab <- c(bd1 = "joint_democracy_polity_i_x_polity_j", bd2 = "trade_import_ij", bd3 = "joint_IGO")
auth <- c(bd1 = -0.017, bd2 = -0.167, bd3 = 0.017)   # authors' archived 2007 chain (Rep_2000_others/OUT.base, 500k scans, first 25% dropped)
xtra <- function(o, lbl) if (is.null(o)) "" else sprintf("%s coef %.2f [%.2f, %.2f]", lbl, mean(o$bd5), quantile(o$bd5, .025), quantile(o$bd5, .975))
for (k in names(pub)) {
  s0 <- summ(o0[[k]]); sh <- summ(oh[[k]]); s26 <- if (is.null(o26)) NULL else summ(o26[[k]]); spl <- if (is.null(opl)) NULL else summ(opl[[k]])
  cv <- essfmt(Filter(Negate(is.null), list(m0 = s0, hist = sh, m26 = s26, placebo = spl)))
  out[[length(out) + 1]] <- mkrow("WSC2007", "Ward, Siverson & Cao 2007 AJPS (gbme latent space, MIDs 2000)", lab[k], k, pub[k], pubci[k], s0, sh, s26, spl,
    paste0("gbme MCMC 25,000 scans (authors 500,000), first 5,000 dropped, odens 20 -> 1,000 draws; mean ll m0 ", round(mean(o0$ll), 1), " (authors' chain -29.2), mean se m0 ", round(mean(o0$s2e)), " (authors 174); ", cv),
    verdict(s0, pub[k], FALSE, tol_abs = 0.005),
    paste0("Headline = Table 2 (main text; year-2000 cross-section, 11 other years in Fig. 2). Authors' gbme.r, binomial, k=3, directed, authors' diffuse priors (routine_uninform.R), seed 123 (authors set none). ",
      "Warm start from the last saved draw of the authors' uninformative-prior Y2K chain (uninform5k_2000, scan 497,000): betas, additive and bilinear effects, and the dyadic residual matrix drawn from its prior at the authors' last (se, rho); gbme.r edited only to accept these starting values (gbme_ws.r). The 500k-scan setting (~11 CPU-h per model) is infeasible; a first attempt with the residual matrix started at 0 had not converged after 25k scans (ll about -45 vs authors' -29; se 0.3 vs ~175) and was discarded (coldstart_* files). Pub values are 2-decimal (tolerance +/-0.005); authors' own archived 2007 chain gives ", auth[k], ". ",
      "hist = dyadic log1p(MID onsets 1995-99) as Xd layer + state's onsets with others as sender and receiver covariates; m26 = same for 1994-98; placebo = mean CINC stronger-share 1995-99 as Xd layer. ",
      xtra(oh, "hist_dyad"), "; ", xtra(o26, "hist_dyad_w26"), "; ", xtra(opl, "bm_cinc"), ". Only 1.2% of dyads have history>0 but they carry 76 of 114 MIDs in 2000."))
}

}
## ---------- MINHAS2022: amen pa2018_version ame_repL, bin, R=2; warm start from authors' final MCMC state
ld <- function(app, m) { f <- file.path(B, "MINHAS2022", paste0("fit_", app, "_", m, ".rds")); ok <- file.path(B, "MINHAS2022", paste0("log_", app, "_", m, ".txt.ok"))
  if (!file.exists(f) || !file.exists(ok)) return(NULL); readRDS(f)$BETA }
apps <- list(
  rs = list(tag = "MINHAS2022_RS", study = "Minhas et al. 2022 PSRM, AME re-estimation of Reiter & Stam 2003 (directed dyads 1945-91)", tab = "B1",
    keys = list(c("pers_democ_directed_dyad", "pdemdtar.dyad", 0.252, 0.066), c("democ_pers_directed_dyad", "pdemdin.dyad", 0.114, 0.089)),
    set = "authors 100,000 scans (burn 50,000, odens 10, ~1 day); ours burn 500 + 3,000 scans, odens 10 -> 300 draws"),
  weeks = list(tag = "MINHAS2022_WEEKS", study = "Minhas et al. 2022 PSRM, AME re-estimation of Weeks 2012 (directed dyads 1950-99)", tab = "B2",
    keys = list(c("machine", "machinejlw_1.row", -0.003, 0.041), c("junta", "juntajlw_1.row", 0.031, 0.045), c("boss", "bossjlw_1.row", -0.040, 0.046), c("strongman", "strongmanjlw_1.row", 0.040, 0.044)),
    set = "authors 100,000 scans (burn 50,000, odens 25, ~1 day); ours burn 500 + 2,500 scans, odens 10 -> 250 draws"))
for (app in names(apps)) {
  if (!(apps[[app]]$tag %in% SEL)) next
  A <- apps[[app]]; b0 <- ld(app, "m0"); bh <- ld(app, "hist"); b26 <- ld(app, "m26"); bpl <- ld(app, "placebo")
  if (is.null(b0) || is.null(bh)) { cat("skip", app, "(m0/hist missing)\n"); next }
  hx <- function(bb, v) if (is.null(bb)) "" else { s <- summ(bb[, v]); sprintf("%s %.3f (sd %.3f, ESS %.0f)", sub(".dyad", "", v, fixed = TRUE), s["b"], s["se"], s["ess"]) }
  for (kk in A$keys) {
    term <- kk[2]; pb <- as.numeric(kk[3]); pse <- as.numeric(kk[4]); stopifnot(term %in% colnames(b0), term %in% colnames(bh))
    s0 <- summ(b0[, term]); sh <- summ(bh[, term]); s26 <- if (is.null(b26)) NULL else summ(b26[, term]); spl <- if (is.null(bpl)) NULL else summ(bpl[, term])
    cv <- essfmt(Filter(Negate(is.null), list(m0 = s0, hist = sh, m26 = s26, placebo = spl)))
    out[[length(out) + 1]] <- mkrow(A$tag, A$study, kk[1], term, pb, sprintf("[%.3f, %.3f] (+/-1.96 se)", pb - 1.96 * pse, pb + 1.96 * pse), s0, sh, s26, spl,
      paste0(A$set, "; ", cv), verdict(s0, pb, abs(pb / pse) > 1.96),
      paste0("AME is the paper's headline estimator (main-text Table 2 and Sec. 5); coefficients in Appendix Table ", A$tab, " (pub se ", pse, "). amen pa2018_version ame_repL, bin, R=2, authors' seed 6886, warm start from authors' final MCMC state (ameFit$startVals); short chain, so ESS is low. ",
        "Authors' model already has peace years + 3 cubic splines. hist adds dyadic log1p(dyad MID onsets t-1..t-5) and log1p(both states' onsets with others t-1..t-5) as dyadic covariates; m26 = t-2..t-6; placebo = mean CINC stronger-share t-1..t-5 (dyad-years with no lagged CINC filled with the year mean). ",
        hx(bh, "hist_dyad_w15.dyad"), "; ", hx(bh, "hist_mono_w15.dyad"), "; ", hx(b26, "hist_dyad_w26.dyad"), "; ", hx(bpl, "bm_cinc_w15.dyad"), "."))
  }
}
## Gibler: infeasible (memory)
if ("MINHAS2022_GIBLER" %in% SEL) for (kk in list(c("parity_at_entry_year", "firstparity", -0.221, 0.084), c("contemporary_parity", "parity4", 0.006, 0.085))) {
  pb <- as.numeric(kk[3]); pse <- as.numeric(kk[4])
  out[[length(out) + 1]] <- mkrow("MINHAS2022_GIBLER", "Minhas et al. 2022 PSRM, AME re-estimation of Gibler 2017 (undirected dyads 1816-2000)", kk[1], kk[2], pb,
    sprintf("[%.3f, %.3f] (+/-1.96 se)", pb - 1.96 * pse, pb + 1.96 * pse), NULL, NULL, NULL, NULL,
    "not estimated", "INFEASIBLE",
    paste0("AME headline (main text Table 2, Sec. 5.3, Fig. 9); pub = Appendix Table B3 AME averaged over the authors' 10 chains (se ", pse, "). Symmetric AME R=2 over 179 years (1822-2000 after GML restriction) x up to 185 states: m0 smoke test (10 scans) was OOM-killed at >5.6 GB on the ~7 GB shared machine; authors ran 10 parallel 50,000-scan chains on a 137 GB box (~4 days). Data and code downloaded; not estimated."))
}

## ---------- WH2011: longitudinal social-relations probit (authors' LSNProbitMID.R), year-specific betas 1991-2000
if ("WH2011" %in% SEL) {
bf <- function(m) { f <- file.path(B, "WH2011", paste0("out_", m), "Beta.txt"); if (!file.exists(f) || !file.exists(file.path(B, "WH2011", paste0("log_wh_", m, ".txt.ok")))) return(NULL); as.matrix(read.table(f)) }
w0 <- bf("m0"); wh <- bf("hist"); w26 <- bf("m26"); wpl <- bf("placebo"); stopifnot(!is.null(w0), !is.null(wh))
base <- c("int", "ally", "exp.ij", "exp.ji", "igo", "dist")
nmof <- list(m0 = base, hist = c(base, "hist_dyad_w15", "hist_mono_w15"), m26 = c(base, "hist_dyad_w26", "hist_mono_w26"), placebo = c(base, "bm_cinc_w15"))
avgv <- function(bb, nm, v) { K <- length(nm); stopifnot(ncol(bb) == 10 * K); rowMeans(bb[, (0:9) * K + match(v, nm), drop = FALSE]) }
nyr <- function(bb, nm, v) { K <- length(nm); sum(sapply(0:9, function(t) { x <- bb[, t * K + match(v, nm)]; quantile(x, .025) > 0 | quantile(x, .975) < 0 })) }
# pub: Fig. 5 medians read off the figure (approximate), averaged over 1991-2000
figpub <- list(ally = c(1.1, -0.3, -7, 1.6, 0, -2.7, -2, 0.2, -1.3, 0), exp.ij = c(-0.2, -0.2, -1, -0.3, -0.2, -0.9, -0.4, -0.4, -0.9, -7.3),
               exp.ji = c(-0.2, 0.1, -0.8, -0.4, 0.1, 0.6, 0.1, 0.2, 0.6, 3.2), igo = c(0, 0.02, 0.35, -0.1, -0.05, -0.02, 0.06, -0.15, 0.13, -0.25))
figsig <- c(ally = "1993 (-)", exp.ij = "1993, 1996, 1999, 2000 (-)", exp.ji = "1993 (-), 2000 (+)", igo = "1993 (+)")
wlab <- c(ally = "alliance_level", exp.ij = "trade_log_exports_i_to_j", exp.ji = "trade_log_imports_j_to_i", igo = "joint_IGOs")
peryr <- NULL
for (v in names(wlab)) {
  s0 <- summ(avgv(w0, nmof$m0, v)); sh <- summ(avgv(wh, nmof$hist, v))
  s26 <- if (is.null(w26)) NULL else summ(avgv(w26, nmof$m26, v)); spl <- if (is.null(wpl)) NULL else summ(avgv(wpl, nmof$placebo, v))
  for (t in 1:10) { x0 <- w0[, (t - 1) * 6 + match(v, base)]; xh <- wh[, (t - 1) * 8 + match(v, nmof$hist)]
    peryr <- rbind(peryr, data.frame(term = v, year = 1990 + t, pub_fig_median = figpub[[v]][t], b0 = mean(x0), lo0 = quantile(x0, .025), hi0 = quantile(x0, .975),
      b_hist = mean(xh), lo_hist = quantile(xh, .025), hi_hist = quantile(xh, .975), row.names = NULL)) }
  hdv <- avgv(wh, nmof$hist, "hist_dyad_w15")
  cv <- essfmt(Filter(Negate(is.null), list(m0 = s0, hist = sh, m26 = s26, placebo = spl)))
  out[[length(out) + 1]] <- mkrow("WH2011", "Westveld & Hoff 2011 AOAS, longitudinal SRM probit, MID initiation, Middle East + US 1991-2000", wlab[v], paste0(v, " (10-year average of year-specific coefs)"),
    round(mean(figpub[[v]]), 2), paste0("figure only; 95% CIs excl 0 in ", figsig[v]), s0, sh, s26, spl,
    paste0("15,000 iterations (authors 7,000,000), burn 3,000, thin 12 -> 1,000 draws; ", cv), "CLOSE",
    paste0("MID application is a co-equal main-text application (Sec. 6; the trade application in Sec. 5 comes first). Authors' LSNProbitMID.R, priors unchanged, seed 123. Pub = Fig. 5 medians read off the figure (approximate). Supplement MID.csv holds only 1991; 1992-2000 rebuilt from COW (directed MID ongoing, side A -> side B, matches released 1991 22/22; IGO r=.998; log real exports r=.955 vs released; alliance: authors' 1991 codes carried forward unless COW pattern changes). ",
      "Row = posterior of the 10-year average of the year-specific coefficient; years with 95% CI excl 0: m0 ", nyr(w0, nmof$m0, v), "/10, hist ", nyr(wh, nmof$hist, v), "/10 (per-year in work_net/WH2011/WH2011_peryear.csv). ",
      "Verdict CLOSE = qualitative 1991-99 pattern of Fig. 5 recovered on the rebuilt panel (e.g. 1993 alliance and IGO, negative distance), not a numeric check; the extreme 2000 values in Fig. 5 (exports -7.3, imports +3.2) are not recovered, which drives most of the gap in 10-year averages. The model already carries AR(1) sender/receiver effects (lagged dispute propensity). hist adds year-specific coefs on log1p dyad onsets t-1..t-5 and log1p both states' onsets with others; avg hist_dyad = ",
      sprintf("%.2f [%.2f, %.2f]", mean(hdv), quantile(hdv, .025), quantile(hdv, .975)), "."))
}
write.csv(peryr, file.path(B, "WH2011/WH2011_peryear.csv"), row.names = FALSE)
}
res <- do.call(rbind, out)
print(res[, c("tag", "key_effect", "pub_b", "b0", "se0", "ci0_excl0", "b_hist", "se_hist", "ci_hist_excl0", "pct_hist", "pct_26", "pct_placebo", "verdict")], digits = 3)
