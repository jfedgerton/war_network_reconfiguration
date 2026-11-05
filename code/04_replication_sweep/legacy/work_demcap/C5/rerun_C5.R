# C5: Mousseau 2013 ISQ "The Democratic Peace Unraveled: It's the Economy". Headline = Table 1 Model 3:
#   logit mzfmidl CIEl dml PolDis lncprt mjpw cntg dist numstate fpceyrs fspl1 fspl2 fspl3, cl(ID)   (fatal MID onset t+1, 1961-2001, N=301,072)
# Data: Mousseau's original files as bundled in Dafoe, Oneal & Russett 2013 ISQ Dataverse doi:10.7910/DVN/Y2EMIO
#   ("Mous_ISQ DPUNR.do", intermediate DPUNRAV.dta, CINE 2011.2.15 ccode1/ccode2, Political Dissimilarity Werner).
#   None of Mousseau's own deposits (27254, 27248, N5RMUT) contain this paper's data. Deposited DPUNRAV.dta is the pre-merge
#   version, so the Table-1-relevant merges (CINE, PolDis) and CIEl/dml construction are rebuilt from the author's do file.
# b_pub: Table 1 Model 3 (2 dp): CIEl -0.85 (0.14)***, dml 0.03 (0.03) ns; DOR 2013 log (12-11-23_log.log) replicates identically.
# Key terms (abstract): "there is no correlation of democracy with peace, and contract-intensive economy is one of the most powerful..."
CODE <- Sys.getenv("SWEEP_CODE"); stopifnot(nzchar(CODE))
source(file.path(CODE, "work_demcap/common_c125.R"))
WD <- file.path(ROOT, "work_demcap/C5"); W <- file.path(OUT, "work_demcap/C5"); dir.create(W, recursive = TRUE, showWarnings = FALSE); sink(file.path(W, "log_C5.txt"), split = TRUE)
rd <- function(f) as.data.frame(haven::zap_labels(haven::read_dta(file.path(WD, f))))
d <- rd("DPUNRAV.dta")[, c("ccode1", "ccode2", "year", "ID", "mzfmidl", "polity21", "polity22", "lncprt", "mjpw", "cntg", "dist", "numstate", "fpceyrs", "fspl1", "fspl2", "fspl3")]
stopifnot(!anyDuplicated(d[, c("ccode1", "ccode2", "year")]), min(d$year) >= 1960)
c1 <- rd("cine1.dta")[, c("ccode1", "year", "CIE1")]; c2 <- rd("cine2.dta")[, c("ccode2", "year", "CIE2")]
stopifnot(!anyDuplicated(c1[, 1:2]), !anyDuplicated(c2[, 1:2]))
pdis <- rd("poldis1.dta")[, c("ccode1", "ccode2", "year", "PolDis")]; stopifnot(!anyDuplicated(pdis[, 1:3]))
n0 <- nrow(d)
d <- d |> left_join(c2, by = c("ccode2", "year")) |> left_join(c1, by = c("ccode1", "year")) |> left_join(pdis, by = c("ccode1", "ccode2", "year"))
stopifnot(nrow(d) == n0)
# author's construction: CIEl = lower of the two CIE; dml = lower polity2 (-99 already set missing in DPUNRAV)
d$CIEl <- ifelse(is.na(d$CIE1) | is.na(d$CIE2), NA, pmin(d$CIE1, d$CIE2))
d$dml  <- ifelse(is.na(d$polity21) | is.na(d$polity22), NA, pmin(d$polity21, d$polity22))
stopifnot(all(d$polity21 >= -10 & d$polity21 <= 10, na.rm = TRUE))
d <- d[!is.na(d$CIEl), ]
stopifnot(all(d$mzfmidl %in% c(0, 1, NA)))
d$yr_out <- d$year + 1
f <- mzfmidl ~ CIEl + dml + PolDis + lncprt + mjpw + cntg + dist + numstate + fpceyrs + fspl1 + fspl2 + fspl3
K <- c("CIEl", "dml")
PUB <- list(CIEl = c(-0.85, 0.14), dml = c(0.03, 0.03))
r <- fit_keys(f, d, K, "ID", file.path(W, "cache_repro.rds"))
cat("repro N =", r$n, "(pub 301072)\n")
# sanity vs DOR's 6-dp replication of Model 2 (no PolDis): CIEl -0.80234605, dml 0.0024645347, N 321568
chk <- fit_keys(mzfmidl ~ CIEl + dml + lncprt + mjpw + cntg + dist + numstate + fpceyrs + fspl1 + fspl2 + fspl3, d, K, "ID", file.path(W, "cache_M2.rds"))
ds <- common_sample(d, f, "ID", "ccode1", "ccode2", "yr_out"); rm(d); gc()
M <- run_six(f, ds, K, "ID", W)
cat("years in common sample:", range(ds$year), "\n")
# ---- sensitivity hook: runs only under run_sensitivity.R (SWEEP_SENS=1), never in the main sweep ----
if (Sys.getenv("SWEEP_SENS") == "1") {
  source(file.path(CODE, "hist/sens_helper.R"))
  sens_run(id = "C5", dat = ds, y = "mzfmidl", tie = "CIEl",
           ctrl = c("dml", "PolDis", "lncprt", "mjpw", "cntg", "dist", "numstate", "fpceyrs", "fspl1", "fspl2", "fspl3"), fe = NULL,
           fit = function(dat, add) { ff <- if (length(add)) update(f, as.formula(paste(". ~ . +", paste(add, collapse = " + ")))) else f
             fit_keys(ff, dat, "CIEl", "ID")$k[["CIEl"]] })
}

est <- "logit, dyad-clustered SE (fatal MID onset t+1)"; sc <- "Mousseau 2013 ISQ"
ev <- "abstract: 'there is no correlation of democracy with peace, and contract-intensive economy is one of the most powerful'"
FAM <- c(CIEl = "CAPITALISM", dml = "DEMOCRACY")
for (k in K) {
  bp <- PUB[[k]][1]; p2 <- 2 * pnorm(-abs(bp / PUB[[k]][2])); attr(bp, "sig") <- p2 < .05
  rb <- r$k[[k]][["b"]]
  st <- if (k == "dml") { # 2-dp publication: 0.03 compatible with [0.025,0.035); judge by rounding + significance
    if (round(rb, 2) == 0.03 && (r$k[[k]][["p"]] < .05) == FALSE) "REPRODUCED" else "NOT_REPRODUCED" } else repro_status(bp, r$k[[k]])
  nt <- paste0("Table 1 Model 3 (CIEl + dml + PolDis; fatal MIDs). b_pub 2 dp from paper (DOR 2013 log replicates identically). ",
    "Rebuilt CINE/PolDis merges + CIEl/dml from author's do file on deposited pre-merge DPUNRAV.dta. Repro N=", r$n, " vs pub 301,072. ",
    "Sanity Model 2: CIEl=", signif(chk$k$CIEl[["b"]], 6), " dml=", signif(chk$k$dml[["b"]], 4), " N=", chk$n, " (DOR 6dp: -0.802346, 0.002465, 321,568). ",
    "DV fatal MID; history vars are all-MID onsets (GML), anchored on DV year (t+1). Author controls fatal peace years + 3 splines. ",
    r$n - M$m0$n, " obs dropped vs repro for history/benchmark missingness. ",
    if (k == "dml") "Published 0.03 at 2 dp: reproduced judged by rounding + non-significance. Mousseau's claim is that dml is NOT significant once CIE included. bD (lagged joint democracy) is NOT a placebo here. " else "")
  write_dc("C5", sc, FAM[[k]], k, ev, st, bp, rb, est, r$n, M, nt)
}
sink()
