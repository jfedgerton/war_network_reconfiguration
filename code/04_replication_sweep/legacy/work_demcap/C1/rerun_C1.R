# C1: Gartzke 2007 AJPS "The Capitalist Peace". Headline = Table 1 Model 5 (all liberal-economic terms, no region dummies):
#   logit maoznewl demlo demhi deplo capopenl rgdppclo gdpcontg sun2cati contig logdstab majpdyds alliesr lncaprt _spline1-3, cluster(dyadid)
#   (authors' do file capitalistpeace_ajps_0107.do line ~1241; b_pub from authors' log capitalistpeace_012007.log; N=166,140)
# Data: Gartzke's UCSD page capitalistpeace_012007.dta (the Choi 2011 Dataverse RESQGA deposit holds only raw inputs + same do/log).
# Key terms (abstract): "Economic development, free markets, and similar interstate interests all anticipate a lessening of militarized disputes"
#   capopenl (IMF Fin. Open. Low), rgdppclo (GDPPC Low; constituent of GDPPC x contig), sun2cati (interest similarity/affinity), demlo (joint democracy weak link).
# DV maoznewl = Maoz dyadic MID in t+1 -> history anchored on outcome year (year+1).
CODE <- Sys.getenv("SWEEP_CODE"); stopifnot(nzchar(CODE))
source(file.path(CODE, "work_demcap/common_c125.R"))
WD <- file.path(ROOT, "work_demcap/C1"); W <- file.path(OUT, "work_demcap/C1"); dir.create(W, recursive = TRUE, showWarnings = FALSE); sink(file.path(W, "log_C1.txt"), split = TRUE)
d <- haven::read_dta(file.path(WD, "cp.dta"), col_select = c(statea, stateb, year, dyadid, maoznewl, demlo, demhi, deplo, capopenl,
  rgdppclo, gdpcontg, sun2cati, contig, logdstab, majpdyds, alliesr, lncaprt, `_spline1`, `_spline2`, `_spline3`))
d <- as.data.frame(haven::zap_labels(d)); names(d) <- sub("^_spline", "spline", names(d))
stopifnot(all(d$maoznewl %in% c(0, 1)), !anyDuplicated(d[, c("statea", "stateb", "year")]))
d$yr_out <- d$year + 1
f <- maoznewl ~ demlo + demhi + deplo + capopenl + rgdppclo + gdpcontg + sun2cati + contig + logdstab + majpdyds + alliesr + lncaprt + spline1 + spline2 + spline3
K <- c("capopenl", "rgdppclo", "sun2cati", "demlo")
# published (authors' log): coef and p
PUB <- list(capopenl = c(-.2468004, 0.000), rgdppclo = c(.0002481, 0.000), sun2cati = c(-.9823773, 0.000), demlo = c(-.017071, 0.151))
r <- fit_keys(f, d, K, "dyadid", file.path(W, "cache_repro.rds"))
cat("repro N =", r$n, "(pub 166140)\n")
ds <- common_sample(d, f, "dyadid", "statea", "stateb", "yr_out"); rm(d); gc()
M <- run_six(f, ds, K, "dyadid", W)
# ---- sensitivity hook: runs only under run_sensitivity.R (SWEEP_SENS=1), never in the main sweep ----
if (Sys.getenv("SWEEP_SENS") == "1") {
  source(file.path(CODE, "hist/sens_helper.R"))
  # headline term capopenl (financial openness; additive). fit_keys without a cache file (run_six caches each model by file name)
  sens_run(id = "C1", dat = ds, y = "maoznewl", tie = "capopenl",
           ctrl = c("demlo", "demhi", "deplo", "rgdppclo", "gdpcontg", "sun2cati", "contig", "logdstab", "majpdyds", "alliesr", "lncaprt",
                    "spline1", "spline2", "spline3"), fe = NULL,
           fit = function(dat, add) { r <- fit_keys(if (length(add)) update(f, as.formula(paste(". ~ . +", paste(add, collapse = " + ")))) else f, dat, K, "dyadid")$k$capopenl
             c(b = unname(r[["b"]]), se = unname(r[["se"]]), p = unname(r[["p"]])) })
}
cat("years in common sample:", range(ds$year), "\n")
est <- "logit, dyad-clustered SE (Maoz MID onset t+1)"; sc <- "Gartzke 2007 AJPS"
ev <- "abstract: 'Economic development, free markets, and similar interstate interests all anticipate a lessening of militarized disputes'"
FAM <- c(capopenl = "CAPITALISM", rgdppclo = "CAPITALISM", sun2cati = "CAPITALISM", demlo = "DEMOCRACY")
for (k in K) {
  bp <- PUB[[k]][1]; attr(bp, "sig") <- PUB[[k]][2] < .05
  st <- repro_status(bp, r$k[[k]])
  nt <- paste0("Table 1 Model 5; b_pub from authors' log (p=", PUB[[k]][2], "). Authors control peace-year splines (_spline1-3, no peace-years count). ",
    "History anchored on DV year (t+1). Sample 1949-1992 all within GML range; ", r$n - M$m0$n, " obs dropped vs orig (benchmark/history missingness). ",
    if (k == "rgdppclo") "Constituent term of GDPPC x contiguity interaction (gdpcontg kept). " else "",
    if (k == "demlo") "Democracy key term: Gartzke's claim is that it is NOT significant once capitalist terms are included. bD (lagged joint democracy) is NOT a placebo here. " else "")
  write_dc("C1", sc, FAM[[k]], k, ev, st, bp, r$k[[k]][["b"]], est, r$n, M, nt)
}
sink()
