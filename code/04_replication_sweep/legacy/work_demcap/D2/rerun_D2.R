# D2: Dafoe 2011 AJPS, "Statistical Critiques of the Democratic Peace: Caveat Emptor"
# Headline (paper's preferred capitalist-controlled model, p.258-9: "Using Model G4, but with corrected temporal specification, no regional
#   dummies, and values of 0 imputed to FinancialOpennessLow for communist countries after 1965, a one s.d. increase in DemocracyLow is
#   estimated to reduce the probability of a MID by 10% (p < 0.001)"). Author's do-file ("10.04.20 Gartzke 2007 replication.do", l.1091):
#   logit maoznewl demlo demhi capopenlc1 deplo contig rgdppclo gdpcontg logdstab majpdyds alliesr lncaprt maozyrs _spline*, cluster(dyadid)
#   with the author's hand-copied coefficients in comments: demlo -.0727517, capopenlc1 -.1177215 (no log posted).
# Data: Gartzke 2007 capitalistpeace_012007.dta (same file the C1 rerun uses; hard-linked as cp.dta). maozyrs/_spline1-3 are Gartzke's btscs output.
# capopenlc1 = capopenl, set to 0 if missing & year>1965 & either side in Dafoe's communist list (kountry names -> COW codes below).
# Key terms: demlo (DemocracyLow; DEMOCRACY) and capopenlc1 (FinancialOpennessLow, communist-imputed; CAPITALISM).
CODE <- Sys.getenv("SWEEP_CODE"); stopifnot(nzchar(CODE))
source(file.path(CODE, "work_demcap/D_common.R"))
suppressPackageStartupMessages(library(haven))
WD <- file.path(ROOT, "work_demcap/D2"); W <- file.path(OUT, "work_demcap/D2"); dir.create(W, recursive = TRUE, showWarnings = FALSE); sink(file.path(W, "log_D2.txt"), split = TRUE)

d <- read_dta(file.path(WD, "cp.dta"), col_select = c(statea, stateb, year, maoznewl, demlo, demhi, capopenl, deplo, contig, rgdppclo,
                                                      gdpcontg, logdstab, majpdyds, alliesr, lncaprt, maozyrs, `_spline1`, `_spline2`, `_spline3`))
d <- as.data.frame(zap_labels(d)); names(d) <- sub("^_spline", "spline", names(d))
stopifnot(!anyDuplicated(d[, c("statea", "stateb", "year")]), all(d$maoznewl %in% c(0, 1, NA)))
d$dyadid <- d$statea * 1000 + d$stateb

# Dafoe's communist list (kountry standard names) as COW codes. 816 = Vietnam (DRV/unified); 817 (Rep. of Vietnam) not in list.
COMM <- c(540, 371, 373, 370, 346, 355, 811, 710, 40, 315, 366, 372, 265, 310, 705, 703, 367, 368, 541, 731, 290, 360, 365, 702, 701, 816, 369, 704)
d$capopenlc1 <- d$capopenl
imp <- is.na(d$capopenl) & d$year > 1965 & (d$statea %in% COMM | d$stateb %in% COMM)
d$capopenlc1[imp] <- 0
cat("capopenlc1 imputed to 0:", sum(imp), "rows\n")

f0 <- maoznewl ~ demlo + demhi + capopenlc1 + deplo + contig + rgdppclo + gdpcontg + logdstab + majpdyds + alliesr + lncaprt +
  maozyrs + spline1 + spline2 + spline3
K <- c("demlo", "capopenlc1")
r <- fit_keys(f0, d, K, "dyadid", file.path(W, "cache_repro.rds")); n_orig <- r$n
bp <- c(demlo = -.0727517, capopenlc1 = -.1177215)
st_k <- c(demlo = { b <- bp[["demlo"]]; attr(b, "sig") <- TRUE; repro_status(b, r$k$demlo) },                 # text: p < 0.001
          capopenlc1 = if (abs(r$k$capopenlc1[["b"]] - bp[["capopenlc1"]]) / abs(bp[["capopenlc1"]]) <= .10) "REPRODUCED" else "NOT_REPRODUCED")
print(st_k); cat("N repro:", n_orig, " years:", range(d$year[complete.cases(d[, all.vars(f0)])]), "\n")

d <- d[, c(all.vars(f0), "dyadid", "statea", "stateb", "year")]
ds <- common_sample(d, f0, "dyadid", "statea", "stateb", "year"); rm(d); gc()
M <- run_six(f0, ds, K, "dyadid", W)

# ---- sensitivity hook: runs only under run_sensitivity.R (SWEEP_SENS=1), never in the main sweep ----
if (Sys.getenv("SWEEP_SENS") == "1") {
  source(file.path(CODE, "hist/sens_helper.R"))
  sens_run(id = "D2", dat = ds, y = "maoznewl", tie = "demlo",
           ctrl = c("capopenlc1", "demhi", "deplo", "contig", "rgdppclo", "gdpcontg", "logdstab", "majpdyds", "alliesr", "lncaprt", "maozyrs", "spline1", "spline2", "spline3"), fe = NULL,
           fit = function(dat, add) { ff <- if (length(add)) update(f0, as.formula(paste(". ~ . +", paste(add, collapse = " + ")))) else f0
             fit_keys(ff, dat, "demlo", "dyadid")$k[["demlo"]] })
}

sc <- "Dafoe 2011 AJPS"; est <- "logit, dyad-clustered SE (fixest feglm)"
nt <- paste0("Preferred corrected Gartzke-G4 model (MID onset maoznewl, 1966-1992 via capopen; no regions, peace-years+3 splines, communist capopen=0). ",
  "b_pub = author's hand-copied coefficients in do-file (no log posted). ", n_orig - M$m0$n, " obs dropped vs orig (common sample incl. benchmarks): Iceland, Mauritius, Cambodia and 1991-92 post-Soviet dyads (no lagged CINC/polity2; many are communist-imputed capopen=0 rows), 15 of 882 MIDs; ",
  "this drop alone moves demlo from -0.074 (full) to b0=-0.045 and capopenlc1 from -0.115 to -0.183 (diag_drop.R). ")
write_dc("D2", sc, "DEMOCRACY", "demlo", "p.259: 'even when controlling for capitalist variables, democracy has a strong and significant association'",
  st_k[["demlo"]], bp[["demlo"]], r$k$demlo[["b"]], est, n_orig, M,
  paste0(nt, "Text: 1-SD demlo -> -10% P(MID), p<.001. bD (lagged joint democracy) NOT a placebo."))
write_dc("D2", sc, "CAPITALISM", "capopenlc1", "abstract: 'Gartzke's claim that capitalist dynamics explain away the democratic peace'",
  st_k[["capopenlc1"]], bp[["capopenlc1"]], r$k$capopenlc1[["b"]], est, n_orig, M,
  paste0(nt, "FinancialOpennessLow (Gartzke's IMF capital-openness, lower dyadic), communist-imputed; published p not reported (status on magnitude only)."))
sink()
