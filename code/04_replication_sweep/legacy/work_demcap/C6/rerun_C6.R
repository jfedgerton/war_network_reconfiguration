# C6: Gartzke & Hewitt 2010 II, "International Crises and the Capitalist Peace" (demcap key terms; tie-focused rerun in work/99)
# Headline: crisis-onset Model A5 (authors' log_01.txt; model used for Figure 2):
#   logit crisis1 demlo2 demhi2 deplo capopnl2 gdppcln gdplncon sun2cati contig logdstab majpdyds alliesr lncaprt nafmeast crisyrs _spline1-3, cluster(dyadid)
# Key terms: capopnl2 (financial openness, low), gdppcln (development, ln GDPpc low; constituent of gdppcln x contig), demlo2 (weak-link democracy).
# b_pub from authors' log (A5): capopnl2 -.1376699 (p .009), gdppcln .4852436 (p .003), demlo2 .0025094 (p .905); N = 254,018.
# Data: slim copy of authors' temp_data_01.dta (their btscs crisyrs/_spline*) saved by work/99/rerun_99.R.
CODE <- Sys.getenv("SWEEP_CODE"); stopifnot(nzchar(CODE))
source(file.path(CODE, "work_demcap/common_c3467.R"))
WD <- file.path(ROOT, "work_demcap/C6"); W <- file.path(OUT, "work_demcap/C6"); dir.create(W, recursive = TRUE, showWarnings = FALSE); sink(file.path(W, "log_C6.txt"), split = TRUE)
d <- readRDS(file.path(ROOT, "work/99/slim_99.rds"))
stopifnot(all(d$crisis1 %in% c(0, 1, NA)))
f0 <- crisis1 ~ demlo2 + demhi2 + deplo + capopnl2 + gdppcln + gdplncon + sun2cati + contig + logdstab + majpdyds +
  alliesr + lncaprt + nafmeast + crisyrs + spline1 + spline2 + spline3
K <- c("capopnl2", "gdppcln", "demlo2"); n_orig <- 254018
bp <- list(capopnl2 = structure(-0.1376699, sig = TRUE), gdppcln = structure(0.4852436, sig = TRUE), demlo2 = structure(0.0025094, sig = FALSE))
r <- fit_keys(f0, d, K, "dyadid", file.path(W, "cache_repro.rds")); stopifnot(r$n == n_orig)
ds <- common_sample(d, f0, "dyadid", "statea", "stateb", "year"); rm(d); gc()
M <- run_six(f0, ds, K, "dyadid", W)
# ---- sensitivity hook: runs only under run_sensitivity.R (SWEEP_SENS=1), never in the main sweep ----
if (Sys.getenv("SWEEP_SENS") == "1") {
  source(file.path(CODE, "hist/sens_helper.R"))
  sens_run(id = "C6", dat = ds, y = "crisis1", tie = "capopnl2",
           ctrl = c("demlo2", "demhi2", "deplo", "gdppcln", "gdplncon", "sun2cati", "contig", "logdstab", "majpdyds", "alliesr", "lncaprt", "nafmeast", "crisyrs", "spline1", "spline2", "spline3"), fe = NULL,
           fit = function(dat, add) { ff <- if (length(add)) update(f0, as.formula(paste(". ~ . +", paste(add, collapse = " + ")))) else f0
             fit_keys(ff, dat, "capopnl2", "dyadid")$k[["capopnl2"]] })
}

est <- "logit, dyad-clustered SE (ICB crisis onset)"; sc <- "Gartzke & Hewitt 2010 (II)"
nt <- paste0("Model A5, ICB crisis onset (kept as published; history vars are MID-based). Authors have crisyrs+3 splines. b_pub from authors' log_01. ",
  n_orig - M$m0$n, " obs dropped vs orig (common sample incl. CINC/democ benchmarks).")
ev <- "authors' A5 spec: capital openness, development (GDPpc x contig), weak-link democracy"
for (k in K) {
  fam <- if (k == "demlo2") "DEMOCRACY" else "CAPITALISM"
  extra <- switch(k, gdppcln = " gdppcln is the constituent term of gdppcln x contig (effect for non-contiguous dyads).",
                  demlo2 = " Null in published A5 (p=.905); pct changes off a near-zero b0 are unstable. bD NOT a placebo for this democracy term.", "")
  write_dc("C6", sc, fam, k, ev, repro_status(bp[[k]], r$k[[k]]), bp[[k]], r$k[[k]][["b"]], est, n_orig, M, paste0(nt, extra))
}
sink()
