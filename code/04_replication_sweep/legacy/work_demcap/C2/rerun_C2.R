# C2: Choi 2011 ISQ "Re-Evaluating Capitalist and Democratic Peace Models". Preferred model = Table 2 Model 1
#   (MIDs, politically relevant dyads, "corrected for all four of the previously discussed errors"):
#   Gartzke T1M5 + peace years (maozyrs) + Political Distance (demhi - demlo) replacing Democracy (High), PRD sample.
# Data: Gartzke's capitalistpeace_012007.dta (Choi's Dataverse RESQGA deposit contains only Gartzke's raw inputs + do/log;
#   Choi's own code not deposited, so the corrections are rebuilt from the paper text). Paper: author's site via Wayback (user's machine).
# b_pub from Table 2 Model 1 (4 dp; significance stars are ONE-tailed; converted to two-tailed from b/se for the p<.05 check).
# Key terms: demlo (Democracy Low), capopenl (Fin. Open. Low), rgdppclo (GDPPC Low; constituent of interaction), sun2cati (INTERESTS).
# DV maoznewl = Maoz dyadic MID in t+1 -> history anchored on outcome year (year+1).
CODE <- Sys.getenv("SWEEP_CODE"); stopifnot(nzchar(CODE))
source(file.path(CODE, "work_demcap/common_c125.R"))
WD <- file.path(ROOT, "work_demcap/C2"); W <- file.path(OUT, "work_demcap/C2"); dir.create(W, recursive = TRUE, showWarnings = FALSE); sink(file.path(W, "log_C2.txt"), split = TRUE)
d <- haven::read_dta(file.path(ROOT, "work_demcap/C1/cp.dta"), col_select = c(statea, stateb, year, dyadid, maoznewl, demlo, demhi, deplo, capopenl,
  rgdppclo, gdpcontg, sun2cati, contig, logdstab, majpdyds, alliesr, lncaprt, maozyrs, `_spline1`, `_spline2`, `_spline3`))
d <- as.data.frame(haven::zap_labels(d)); names(d) <- sub("^_spline", "spline", names(d))
stopifnot(all(d$maoznewl %in% c(0, 1)), !anyDuplicated(d[, c("statea", "stateb", "year")]))
print(table(contig = d$contig, majp = d$majpdyds, useNA = "ifany"))
d$poldist <- d$demhi - d$demlo
d$yr_out <- d$year + 1
# sanity: Choi T1 Model 3 (all dyads, peace years + poldist): pub demlo -0.0370, N 166,140
f3 <- maoznewl ~ demlo + poldist + deplo + capopenl + rgdppclo + gdpcontg + sun2cati + contig + logdstab + majpdyds + alliesr + lncaprt + maozyrs + spline1 + spline2 + spline3
chk <- fit_keys(f3, d, c("demlo", "poldist", "maozyrs"), "dyadid", file.path(W, "cache_T1M3.rds"))
# politically relevant dyads: contiguous or includes a major power
d <- d[d$contig == 1 | d$majpdyds == 1, ]
f <- f3
K <- c("demlo", "capopenl", "rgdppclo", "sun2cati")
PUB <- list(demlo = c(-0.0210, 0.0124), capopenl = c(-0.1459, 0.0470), rgdppclo = c(0.0001, 0.0000), sun2cati = c(-0.6558, 0.1580))
r <- fit_keys(f, d, K, "dyadid", file.path(W, "cache_repro.rds"))
cat("repro N =", r$n, "(pub 14750)\n")
ds <- common_sample(d, f, "dyadid", "statea", "stateb", "yr_out"); rm(d); gc()
M <- run_six(f, ds, K, "dyadid", W)
# ---- sensitivity hook: runs only under run_sensitivity.R (SWEEP_SENS=1), never in the main sweep ----
if (Sys.getenv("SWEEP_SENS") == "1") {
  source(file.path(CODE, "hist/sens_helper.R"))
  sens_run(id = "C2", dat = ds, y = "maoznewl", tie = "demlo",   # [2026-11-10 audit2] hook moved to the new headline term (demlo)
           ctrl = c("sun2cati", "poldist", "deplo", "capopenl", "rgdppclo", "gdpcontg", "contig", "logdstab", "majpdyds", "alliesr", "lncaprt", "maozyrs", "spline1", "spline2", "spline3"), fe = NULL,
           fit = function(dat, add) { ff <- if (length(add)) update(f, as.formula(paste(". ~ . +", paste(add, collapse = " + ")))) else f
             fit_keys(ff, dat, "demlo", "dyadid")$k[["demlo"]] })
}

est <- "logit, dyad-clustered SE (Maoz MID onset t+1), politically relevant dyads"; sc <- "Choi 2011 ISQ"
ev <- "abstract: 'while capitalism does not emerge as a cause mitigating militarized disputes ... democracy does'"
FAM <- c(demlo = "DEMOCRACY", capopenl = "CAPITALISM", rgdppclo = "CAPITALISM", sun2cati = "CAPITALISM")
for (k in K) {
  bp <- PUB[[k]][1]
  p2 <- if (PUB[[k]][2] > 0) 2 * pnorm(-abs(bp / PUB[[k]][2])) else 0   # two-tailed p implied by published b/se
  attr(bp, "sig") <- p2 < .05
  st <- if (k == "rgdppclo") { # published only to 4 dp (0.0001): judge by rounding + significance
    if (round(r$k[[k]][["b"]], 4) == 0.0001 && r$k[[k]][["p"]] < .05) "REPRODUCED" else "NOT_REPRODUCED" } else repro_status(bp, r$k[[k]])
  nt <- paste0("Table 2 Model 1 (PRD, MIDs; Choi's fully corrected model). Choi's code not deposited: rebuilt as Gartzke T1M5 + maozyrs + poldist=demhi-demlo, PRD=contig|majpdyds. ",
    "Published stars one-tailed; implied two-tailed p=", signif(p2, 3), ". Authors control peace years + 3 splines. History anchored on DV year (t+1). ",
    "Repro N=", r$n, " vs pub 14,750; ", r$n - M$m0$n, " obs dropped for history/benchmark missingness. ",
    "Sanity T1M3 (all dyads) demlo=", signif(chk$k$demlo[["b"]], 4), " (pub -0.0370), N=", chk$n, ". ",
    if (k == "rgdppclo") "Published to 4 dp only (0.0001); reproduced judged by rounding + significance. Constituent of GDPPC x contig interaction. " else "",
    if (k == "demlo") "bD (lagged joint democracy) is NOT a placebo for this democracy term. " else "")
  write_dc("C2", sc, FAM[[k]], k, ev, st, bp, r$k[[k]][["b"]], est, r$n, M, nt)
}
sink()
