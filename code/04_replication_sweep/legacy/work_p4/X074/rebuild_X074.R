# X074 Oneal & Russett 2005 CMPS "Rule of Three, Let It Be? When More Really Is Better" (10.1080/07388940500339209). REBUILD.
# No archive found (Dataverse title/author search; WebSearch). Paper (author PDF, saramitchell.org/russettoneal04.pdf) Table 1 col 1:
#   logit, all dyads 1885-2001, onset of any MID; lower & higher democracy, lower trade-to-GDP, allies, log capability ratio, contiguous,
#   log distance, major power, years of peace + 3 cubic splines, number of states; SEs adjusted for dyadic clusters (Table 3 note).
#   Published: lower trade-to-GDP -32.3 (se 9.0, p<.001 one-tailed); lower democracy -.0686 (se .0080, p<.001); N = 464,953.
# KEY terms (fixed from the paper BEFORE estimating): lowdep (TIE) and lowdem (DEMOCRACY, KEY: "We focus on the effects of the liberal
#   variables democracy and interdependence"). No capitalist terms in the model.
# Rebuild base: work_p4/X074/or_trade.rds (build_or_trade.R): Oneal's all-dyads 1816-2000 file (DOR 2013) + COW trade / SDP GDP.
ROOT <- Sys.getenv("SWEEP_DATA"); OUT <- Sys.getenv("SWEEP_OUT"); CODE <- Sys.getenv("SWEEP_CODE")
stopifnot(nzchar(ROOT), nzchar(OUT), nzchar(CODE))
source(file.path(CODE, "work_p4/common_p4_b02.R"))
source(file.path(CODE, "work_p8/round2_rule.R"))   # rebuild_verdict(): round-2 rule with registry alpha_pub (2026-10-11)
WO <- file.path(OUT, "work_p4/X074"); dir.create(WO, recursive = TRUE, showWarnings = FALSE)  # port: outputs
W <- file.path(ROOT, "work_p4/X074"); sink(file.path(WO, "log_X074.txt"), split = TRUE)
d <- readRDS(file.path(W, "or_trade.rds"))
d <- d[d$year >= 1885 & d$year <= 2000, ]
f0 <- mzmid1_orig ~ lowdem + highdem + lowdep + anyally + lcaprat + contig150 + ldist + majpow + py + spline1 + spline2 + spline3 + numstate
K <- c("lowdep", "lowdem")
FF <- ff_logit_cl("dyadid")
r0 <- FF(f0, d[complete.cases(d[, c(all.vars(f0), "dyadid")]), ]); print(round(r0, 5)); n_orig <- attr(r0, "n")
cat("rebuilt N:", n_orig, "(published 464,953)\n")
bp <- c(lowdep = -32.3, lowdem = -.0686)
st <- sapply(K, function(k) rebuild_verdict("X074", k, r0[k, 1], r0[k, 3], bp[[k]], TRUE)); print(st)   # round-2 rule (was: sign and p < .05)

# History timing (fix 2026-10-10): mzmid1_orig is MID onset in t+1 (DOR 2013 label "=mzmid[_n+1] if year[_n+1]==year+1"),
# so history and benchmark windows are keyed on the outcome year t+1. Years 1885-2000: the 1822-2011 filter is unaffected.
lab <- attr(haven::read_dta(file.path(ROOT, "work_demcap/D1/dor.dta"), col_select = "mzmid1_orig", n_max = 1)$mzmid1_orig, "label")
cat("mzmid1_orig label:", lab, "\n"); stopifnot(startsWith(lab, "=mzmid[_n+1] if year[_n+1]==year+1"))
# data check against COW MID 5.0 dyadic onsets (peacesciencer cow_mid_dirdisps): onset dyad-years in t+1 must mostly carry mzmid1_orig = 1
mk <- as.data.frame(peacesciencer::cow_mid_dirdisps); mk <- mk[mk$disponset == 1, ]
mk <- unique(paste(pmin(mk$ccode1, mk$ccode2), pmax(mk$ccode1, mk$ccode2), mk$year))
i1 <- paste(d$statea, d$stateb, d$year + 1) %in% mk & d$year < max(d$year); i0 <- paste(d$statea, d$stateb, d$year) %in% mk
s1 <- mean(d$mzmid1_orig[i1] == 1, na.rm = TRUE); s0 <- mean(d$mzmid1_orig[i0] == 1, na.rm = TRUE)
cat(sprintf("lead check: share of COW onset dyad-years with mzmid1_orig = 1: onset in t+1 %.4f (n %d), onset in t %.4f (n %d)\n", s1, sum(i1), s0, sum(i0)))
stopifnot(all(d$statea < d$stateb), sum(i1) > 500, s1 > 0.8, s1 > s0 + 0.3); rm(lab, mk, i1, i0, s1, s0)
d$outcome_year <- d$year + 1
ds <- common_sample(d, f0, "dyadid", "statea", "stateb", "outcome_year"); rm(d); gc()
M <- run_six(f0, ds, K, FF, WO)

# ---- sensitivity hook: runs only under run_sensitivity.R (SWEEP_SENS=1), never in the main sweep ----
if (Sys.getenv("SWEEP_SENS") == "1") {
  source(file.path(CODE, "hist/sens_helper.R"))
  sens_run(id = "X074", dat = ds, y = "mzmid1_orig", tie = "lowdep",   # [2026-11-10 audit2] hook moved to the headline row (lowdep)
           ctrl = c("highdem", "lowdem", "anyally", "lcaprat", "contig150", "ldist", "majpow", "py", "spline1", "spline2", "spline3", "numstate"), fe = NULL,
           fit = function(dat, add) { ct <- FF(if (length(add)) update(f0, as.formula(paste(". ~ . +", paste(add, collapse = " + ")))) else f0, dat)
             c(b = unname(ct["lowdep", 1]), se = unname(ct["lowdep", 2]), p = unname(ct["lowdep", 3])) })
}

meta <- list(xid = "X074", short_cite = "Oneal & Russett 2005", journal = "CMPS", year = 2005, doi = "10.1080/07388940500339209")
est <- "logit, dyad-clustered SE (fixest feglm)"
hc <- "Y (years of peace + 3 cubic splines only)"
nt <- paste0("REBUILT (no archive): Oneal's all-dyads file from Dafoe-Oneal-Russett 2013 (1816-2000; IVs at t, MID onset originators at t+1) restricted to 1885-2000; history/benchmarks keyed on outcome year t+1 ",
  "(paper 1885-2001) + lower trade/GDP from COW trade (smoothed, nominal) / SDP GDP (2011 USD) [paper: Gleditsch v4.1 + O&R pre-1950 data]. ",
  "allies = any COW alliance; contiguous = land or <=150 miles water; log cap ratio from CINC; btscs peace years (D1 rule) + 3 splines; number of states. ",
  "Missing trade not imputed. Rebuilt N ", n_orig, " vs published 464,953; common N ", M$m0$n, ". Sign+significance vs published = ")
write_p4(meta, "lowdep", "TIE", "KEY", "\"We focus on the effects of the liberal variables democracy and interdependence\"",
  st[["lowdep"]], bp[["lowdep"]], r0["lowdep", 1], est, n_orig, M, hc, paste0(nt, st[["lowdep"]], ". Lower trade-to-GDP ratio (weak link)."))
write_p4(meta, "lowdem", "DEMOCRACY", "KEY", "\"We focus on the effects of the liberal variables democracy and interdependence\"",
  st[["lowdem"]], bp[["lowdem"]], r0["lowdem", 1], est, n_orig, M, hc, paste0(nt, st[["lowdem"]], ". Lower Polity2 (weak link). bD NOT a placebo for this term."))
sink()
