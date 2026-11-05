# P9_070 Oneal 2003, JPR 40(6) "Measuring Interdependence and Its Pacific Benefits: A Reply to Gartzke & Li" (10.1177/00223433030406007).
#   TIE. REBUILD (protocol B; no archive).
# Article (staged PDF, Table I): logit, all dyads 1886-1992, onset of a MID; trade dependence = lower bilateral trade/GDP ratio in the dyad;
#   openness = total trade/GDP of the same (lower-dependence) state; "same controls as in Gartzke & Li (2003)" = the G&L 'Measure for
#   Measure' (JPR 2003) set, listed in the companion Table I of Gartzke & Li 2003b / Barbieri & Peters 2003: contiguity, joint democracy,
#   alliance, relative capabilities, peace years + 3 splines; one-tailed p-values.
# HEADLINE (fixed before estimating): Table I, all dyads, all disputes: trade dependence -41.6 (SE 11.7), p .001 one-tailed, N 282,814.
# Rebuild inputs: DOR/O-R all-dyad file + derived trade/GDP file (work_demcap/D1/dor.dta, work_p4/X074/or_trade.rds; mzmid1_orig is the
#   MID onset in t+1, so history is keyed on outcome_year = year + 1 as for X074); openness from COW national trade 3.0 totals over the same
#   GDP series (work_p8/P8_049). Logit with SEs clustered on dyad. Joint democracy = lower Polity (O/R weak link); peace years and
#   splines = or_trade py/spline1-3 (btscs on mzmid1_orig).
ROOT <- Sys.getenv("SWEEP_DATA"); OUT <- Sys.getenv("SWEEP_OUT"); CODE <- Sys.getenv("SWEEP_CODE")
stopifnot(nzchar(ROOT), nzchar(OUT), nzchar(CODE))
source(file.path(CODE, "work_p8/common_b3.R")); source(file.path(CODE, "work_p9/rebuild_common_p9r.R"))
set.seed(123)
W <- file.path(OUT, "work_p9/P9_070"); dir.create(W, recursive = TRUE, showWarnings = FALSE)
sink(file.path(W, "log_P9_070.txt"), split = TRUE)
d <- or_dyads(1885:1991)
nt_ <- cow_nat_trade(); n0 <- nrow(d)
d <- d %>% left_join(nt_ %>% rename(tota = tot), by = c("statea" = "ccode", "year")) %>% left_join(nt_ %>% rename(totb = tot), by = c("stateb" = "ccode", "year"))
stopifnot(nrow(d) == n0)
d <- d %>% mutate(opena = tota * 1e6 / gdp_a, openb = totb * 1e6 / gdp_b, alow = dep_a <= dep_b,
  open_low = ifelse(alow, opena, openb), outcome_year = year + 1)
d$open_low[!is.finite(d$open_low)] <- NA
rhs0 <- c("lowdep", "open_low", "contig150", "lowdem", "anyally", "lcaprat", "py", "spline1", "spline2", "spline3")
vars <- c("mzmid1_orig", rhs0, "dyadid")
cat("dyad-years 1885-1991 (outcome 1886-1992):", nrow(d), " missing:\n"); print(sapply(d[, vars], function(x) sum(is.na(x))))
d0 <- d[complete.cases(d[, vars]), ]
cat("complete-case N:", nrow(d0), "(published 282,814)  onsets", sum(d0$mzmid1_orig), "\n")
TERM <- "lowdep"
fitter <- function(r, data) fit_bin("mzmid1_orig", r, data, "logit", "~dyadid", terms = TERM)
full <- fitter(rhs0, d0)
cat("\nREBUILT vs published: lowdep b", full[[TERM]][["b"]], "(pub -41.6) se", full[[TERM]][["se"]], "(pub 11.7) N", attr(full, "n"), "(pub 282,814)\n")
ds <- common_sample(d, vars, "statea", "stateb", "outcome_year"); rm(d, d0); invisible(gc())
M <- run_six(rhs0, ds, fitter)
if (Sys.getenv("SWEEP_SENS") == "1") {
  source(file.path(CODE, "hist/sens_helper.R"))
  sens_run(id = "P9_070", dat = ds, y = "mzmid1_orig", tie = TERM, ctrl = setdiff(rhs0, TERM), fe = NULL,
           fit = function(dat, add) fitter(c(rhs0, add), dat)[[TERM]])
}
saveRDS(list(M = M, full = full), file.path(W, "M_P9_070.rds"))
meta <- list(pid = "P9_070", cite = "Oneal 2003", journal = "JPR", year = 2003, doi = "10.1177/00223433030406007", family = "TIE")
nt <- sprintf(paste0("REBUILT (no archive): DOR/O-R all dyads, IVs 1885-1991, DV mzmid1_orig = MID onset in t+1 (outcome years 1886-1992); ",
  "lowdep from or_trade; openness of the lower-dependence state = COW national trade 3.0 total / GDP; G&L JPR controls (contig150, lowdem, ",
  "anyally, lcaprat, py + 3 splines). Full-sample rebuilt b=%.6f se=%.6f N=%d vs published -41.6 (11.7) N 282,814. ",
  "Common sample N=%d (history at outcome year)."), full[[TERM]][["b"]], full[[TERM]][["se"]], attr(full, "n"), nrow(ds))
write_p8(meta, TERM, "abstract: high trade dependence (lower bilateral trade/GDP) reduces violent conflict; Table I all dyads, all disputes",
  "Table I, all dyads, all disputes (logit; MID onset 1886-1992)", "logit, G&L controls, one-tailed tests",
  "feglm logit, SE clustered on dyad, all six steps", -41.6, 11.7, 282814, TRUE, TRUE, M, "Y (peace years + splines)", nt)
dir.create(file.path(OUT, "results_p9"), recursive = TRUE, showWarnings = FALSE)
src <- file.path(OUT, "results_p8", "P9_070__lowdep.csv"); stopifnot(file.exists(src), file.rename(src, file.path(OUT, "results_p9", basename(src))))
sink()
