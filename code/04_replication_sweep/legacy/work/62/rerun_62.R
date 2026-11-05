# Rerun 62: Dorussen & Ward 2010 JPR -- Table I Model 1:
#   logit force_y1 low_depend_ksg2 low_polity2 high_polity2 i2_numigo ln_lowgdpindex ln_highgdpindex contiguity ln_distance
#     minorpowers allies ln_capratio force_pceyrs force_spline1-3 lndyadtrade_index lntradelinks3_index maxflowtrade100_index, r
# Headline tie: lntradelinks3_index (trade-links network index); maxflowtrade100_index tracked as companion.
# b_pub: Table I point estimates are reproduced in the authors' robustness memo "Robustness 4" (identical spec, only SEs
#   clustered): lntradelinks3_index -0.20, maxflowtrade100_index -2.10 (2 dp). Significant at 5% even with clustered SEs.
ROOT <- Sys.getenv("SWEEP_DATA"); OUT <- Sys.getenv("SWEEP_OUT"); CODE <- Sys.getenv("SWEEP_CODE")
stopifnot(nzchar(ROOT), nzchar(OUT), nzchar(CODE))
source(file.path(CODE, "work/rerun_common_g7.R"))
suppressPackageStartupMessages(library(haven))
W <- file.path(ROOT, "work/62")
WO <- file.path(OUT, "work/62"); dir.create(WO, recursive = TRUE, showWarnings = FALSE)
sink(file.path(WO, "log_62.txt"), split = TRUE)

# 1. Authors' data (PRIO JPR replication zip, "Dorussen_Ward_Trade Network_JPR.dta")
RHSV <- c("low_depend_ksg2","low_polity2","high_polity2","i2_numigo","ln_lowgdpindex","ln_highgdpindex","contiguity","ln_distance",
          "minorpowers","allies","ln_capratio","force_pceyrs","force_spline1","force_spline2","force_spline3",
          "lndyadtrade_index","lntradelinks3_index","maxflowtrade100_index")
MV <- c("force_y1", RHSV)
d <- as.data.frame(zap_labels(read_dta(file.path(W, "dw_main.dta"), col_select = all_of(c(MV, "ccode1","ccode2","year","dyadid")))))
# lead check (fix 2026-10-10): force_y1[t] == force[t+1] on consecutive dyad records (force read only for this check)
lc <- as.data.frame(zap_labels(read_dta(file.path(W, "dw_main.dta"), col_select = all_of(c("ccode1","ccode2","year","force","force_y1")))))
lc <- lc[!is.na(lc$ccode1) & !is.na(lc$ccode2), ]; lc <- lc[order(lc$ccode1, lc$ccode2, lc$year), ]; nl <- nrow(lc)
cs <- which(lc$ccode1[-1] == lc$ccode1[-nl] & lc$ccode2[-1] == lc$ccode2[-nl] & lc$year[-1] == lc$year[-nl] + 1)
la <- lc$force_y1[cs]; lb <- lc$force[cs + 1]; ok <- !is.na(la) & !is.na(lb)
stopifnot(length(cs) > 1000, sum(ok) > 1000, all(la[ok] == lb[ok]))
cat("lead check: force_y1[t] == force[t+1] on", sum(ok), "consecutive dyad pairs with both non-missing\n"); rm(lc, cs, la, lb, ok, nl)
cat("rows:", nrow(d), " years:", range(d$year), "\n")
BASE <- paste("force_y1 ~", paste(RHSV, collapse = " + "))
fitter <- function(add, dat) fit_feglm_tie(as.formula(paste0(BASE, add)), dat, "lntradelinks3_index", cluster = "ROBUST")

# 2. Reproduce
cat("\n==== Reproduction: Table I Model 1 ====\n")
rp <- fitter("", d)
b_pub <- -0.20
cat("maxflowtrade100_index repro:", rp$all["maxflowtrade100_index", 1], " (pub -2.10)\n")
repro <- if (abs(rp$b - b_pub) / abs(b_pub) <= 0.10 && rp$p < .05) "REPRODUCED" else "NOT_REPRODUCED"
cat("\nb_pub =", b_pub, " b_repro =", rp$b, " p =", rp$p, " status:", repro, "\n")
n_orig <- rp$n

# 3. History merge (undirected dyad-year) and 1822-2011 restriction. DV is force at t+1, so history and benchmark windows are
#    keyed on the outcome year t+1 (fix 2026-10-10; previously keyed on predictor year t). Years 1948-2002: the 1822-2011
#    filter (on the row year) is unaffected by the shift.
d$outcome_year <- d$year + 1
d <- add_history(d, "ccode1", "ccode2", "outcome_year")
d <- d[complete.cases(d[, MV]), c(MV, "year", HIST_VARS)]    # keep only model-complete rows and needed columns (memory)
invisible(gc())
SENS <- sens_hist_only(fitter, d, MV)                          # before bm filter (see helper)
d <- restrict_sample(d, MV); invisible(gc())

# 4. Six models
R <- run_six(fitter, d)
mf0  <- fit_feglm_tie(as.formula(BASE), d, "maxflowtrade100_index", cluster = "ROBUST", quiet = TRUE)
mfdm <- fit_feglm_tie(as.formula(paste0(BASE, " + hist_dyad_w15 + hist_mono_w15")), d, "maxflowtrade100_index", cluster = "ROBUST", quiet = TRUE)

# ---- sensitivity hook: runs only under run_sensitivity.R (SWEEP_SENS=1), never in the main sweep ----
if (Sys.getenv("SWEEP_SENS") == "1") {
  source(file.path(CODE, "hist/sens_helper.R"))
  sens_run(id = "62", dat = d, y = "force_y1", tie = "lntradelinks3_index", ctrl = setdiff(RHSV, "lntradelinks3_index"), fe = NULL,
           fit = function(dat, add) { r <- fitter(if (length(add)) paste0(" + ", paste(add, collapse = " + ")) else "", dat); c(b = r$b, se = r$se, p = r$p) })
}

# 5. Results
write_row(62, "Dorussen & Ward 2010", repro, b_pub, rp$b, "logit, robust SE (fixest::feglm hetero)", "lntradelinks3_index",
  n_orig, R$m0$n, R, ahc = "Y",
  notes = sprintf("Table I M1; DV force MID at t+1 (history and benchmark windows keyed on outcome year t+1). b_pub (2dp) from authors' robustness memo R4 (same spec). Author controls force peace years + 3 splines (Y). Companion maxflowtrade100_index: repro %.3f (pub -2.10); b0=%.3f (p=%.3g) -> b_dm=%.3f (p=%.3g), pct=%.1f. Dropped %d obs vs orig; %s. Sensitivity %s.",
                  rp$all["maxflowtrade100_index", 1], mf0$b, mf0$p, mfdm$b, mfdm$p, pct(mfdm$b, mf0$b), n_orig - R$m0$n, RS_INFO, SENS))
sink()
