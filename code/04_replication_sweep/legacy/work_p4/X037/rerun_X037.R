# X037 Altman, Rojas-de-Galarreta & Urdinez 2021 JPR "An interactive model of democratic peace".
# Headline = Table 2 Model 1 (authors' script_ARU.do):
#   logit fatal_mid c.dem_mean##c.dem_spread trade_inter_grwth dependence cntgdumy logdist capabilities allydumy t t2 t3 if pol_rel==1, vce(cluster year)
#   dem_mean = 100*V-Dem polyarchy mean, dem_spread = 100*|diff|. NOTE t,t2,t3 are a cubic in group(year) (time trend), not peace years.
# Key terms fixed from abstract before estimating: dem_mean ("the higher a dyad's level of democracy ... the lower the probability of fatal MIDs")
#   and the product dem_mean x dem_spread (the paper's 'interactive model'). Estimation sample = authors' stored e(sample) _est_m1 == 1.
ROOT <- Sys.getenv("SWEEP_DATA"); OUT <- Sys.getenv("SWEEP_OUT"); CODE <- Sys.getenv("SWEEP_CODE")
stopifnot(nzchar(ROOT), nzchar(OUT), nzchar(CODE))
source(file.path(CODE, "work_p4/P4_common.R"))
WO <- file.path(OUT, "work_p4/X037"); dir.create(WO, recursive = TRUE, showWarnings = FALSE)  # port: outputs
W <- file.path(ROOT, "work_p4/X037"); sink(file.path(WO, "log_X037.txt"), split = TRUE)
V <- c("statea","stateb","year","dyadid","pol_rel","fatal_mid","mzfatald","dem_mean","dem_spread","v2x_polyarchy_mean","v2x_polyarchy_diff",
       "trade_inter_grwth","dependence","cntgdumy","logdist","capabilities","allydumy","t","t2","t3","_est_m1")
slim <- file.path(WO, "slim.rds")  # port: derived cache -> OUT
if (!file.exists(slim)) saveRDS(as.data.frame(haven::zap_labels(haven::read_dta(file.path(W, "Replication files (Stata)/data_ARU.dta"), col_select = tidyselect::all_of(V)))), slim)
d <- readRDS(slim); names(d)[names(d) == "_est_m1"] <- "est_m1"
dk <- duplicated(d[, c("statea","stateb","year")]) | duplicated(d[, c("statea","stateb","year")], fromLast = TRUE)
cat("rows in duplicated (statea,stateb,year) keys:", sum(dk), " exact dup rows:", sum(duplicated(d)), " in est_m1:", sum(dk & d$est_m1 %in% 1), " statea<stateb share:", mean(d$statea < d$stateb), "\n")
print(head(d[dk, c("statea","stateb","year","dyadid","pol_rel","fatal_mid","dem_mean","est_m1")], 6))   # kept as in authors' data
cat("dem_mean check:", max(abs(d$dem_mean - 100 * d$v2x_polyarchy_mean), na.rm = TRUE), " fatal_mid check:",
    mean(d$fatal_mid == as.numeric(d$mzfatald > 0), na.rm = TRUE), " est_m1 N:", sum(d$est_m1 == 1, na.rm = TRUE), "\n")
d$dmXds <- d$dem_mean * d$dem_spread
d <- d[d$pol_rel == 1, ]
f0 <- fatal_mid ~ dem_mean + dem_spread + dmXds + trade_inter_grwth + dependence + cntgdumy + logdist + capabilities + allydumy + t + t2 + t3
K <- c("dem_mean", "dmXds", "dem_spread")
r <- fit_keys(f0, d, K, "year", file.path(WO, "cache_repro.rds"))
cat("N repro", r$n, " vs authors e(sample)", sum(d$est_m1 == 1, na.rm = TRUE), "\n")
ds <- common_sample(d, f0, "year", "statea", "stateb", "year")
M <- run_six(f0, ds, K, "year", WO)
# ---- sensitivity hook (split only): runs only under run_sensitivity.R (SWEEP_SENS=1), never in the main sweep ----
if (Sys.getenv("SWEEP_SENS") == "1") {
  source(file.path(CODE, "hist/sens_helper.R"))
  # fit_keys with cache = NULL (what run_six runs, uncached): run_six would read the cache_m0/mD/mDM files
  sens_split(id = "X037", dat = ds,
             fit = function(dat, add) { r <- fit_keys(if (length(add)) update(f0, as.formula(paste(". ~ . +", paste(add, collapse = " + ")))) else f0, dat, K, "year")$k$dmXds
               c(b = unname(r[["b"]]), se = unname(r[["se"]]), p = unname(r[["p"]])) })
}
sc <- "Altman, Rojas-de-Galarreta & Urdinez 2021 JPR"; est <- "logit, year-clustered SE (fixest feglm)"
hc <- "N (cubic time trend in calendar year only; no peace years or lagged conflict)"
nt <- paste0("Table 2 Model 1 (fatal MIDs, politically relevant dyads). Paper table inaccessible (JPR paywall; repositorio.uc.cl returned HTML) -> UNVERIFIED; ",
  "online appendix Table A Model 5 (same spec, rare-events logit) reports ORs dem_mean .994 (n.s.), dem_spread 1.058***, product z=-4.92***. ",
  "N repro ", r$n, "; duplicated dyad-year keys in deposit kept as in authors' data; common sample drops ", r$n - M$m0$n, " obs. bD NOT a placebo for democracy terms.")
write_p4("X037", sc, "DEMOCRACY", "KEY", "dem_mean", "abstract: 'the higher a dyad's level of democracy is, the lower the probability of fatal' MIDs",
  "UNVERIFIED", NA, r$k$dem_mean[["b"]], est, r$n, M, hc, paste(nt, "dem_mean is the constituent of the product (effect at spread = 0)."))
write_p4("X037", sc, "DEMOCRACY", "KEY", "dmXds", "abstract: 'Interactions between the democratic level of the dyad ... and its democratic spread'",
  "UNVERIFIED", NA, r$k$dmXds[["b"]], est, r$n, M, hc, paste(nt, "dmXds = dem_mean x dem_spread (product term)."))
sink()
