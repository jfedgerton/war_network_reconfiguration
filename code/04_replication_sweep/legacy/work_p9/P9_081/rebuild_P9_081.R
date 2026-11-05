# P9_081 Xiang, Xu & Keteku 2007, JCR 51(4) "Power: The Missing Link in the Trade Conflict Relationship" (10.1177/0022002707302802). TIE.
#   REBUILD (protocol B; no archive).
# Article (staged PDF pp. 652-656, Table 1): logit, robust SEs, dyad-years 1870-1992 (Barbieri's trade data), MID onset on t-1 covariates;
#   Salience_trade = sqrt[(trade_ij/total trade_i)(trade_ij/total trade_j)]; Model 2 adds Power_1 and Power_2 = each state's CINC score;
#   controls alliance, contiguity, joint democracy, relative capabilities, distance (log), peace years + 3 splines.
# HEADLINE (fixed before estimating): Table 1 Model 2, Salience_trade -3.176 (0.870) ***, N 137,314.
# Rebuild inputs: DOR/O-R all-dyad file + derived trade file (work_demcap/D1/dor.dta, work_p4/X074/or_trade.rds; COW dyadic trade; mzmid1_orig =
#   onset in t+1 against t covariates = the t-1 lag; history at outcome_year); COW national trade 3.0 totals (work_p8/P8_049); capital distance
#   ({peacesciencer}). Joint democracy = (polity2_a + 10) x (polity2_b + 10) (Barbieri's product form; assumption documented); relative
#   capabilities = log CINC ratio; contiguity = land; py + 3 splines from or_trade.
ROOT <- Sys.getenv("SWEEP_DATA"); OUT <- Sys.getenv("SWEEP_OUT"); CODE <- Sys.getenv("SWEEP_CODE")
stopifnot(nzchar(ROOT), nzchar(OUT), nzchar(CODE))
source(file.path(CODE, "work_p8/common_b3.R")); source(file.path(CODE, "work_p9/rebuild_common_p9r.R"))
set.seed(123)
W <- file.path(OUT, "work_p9/P9_081"); dir.create(W, recursive = TRUE, showWarnings = FALSE)
sink(file.path(W, "log_P9_081.txt"), split = TRUE)
d <- or_dyads(1869:1991)
cd <- suppressMessages(create_dyadyears(system = "cow", directed = FALSE, subset_years = 1869:1991) %>% add_capital_distance()) %>% as.data.frame()
nt_ <- cow_nat_trade(); n0 <- nrow(d)
d <- d %>% left_join(cd[, c("ccode1", "ccode2", "year", "capdist")], by = c("statea" = "ccode1", "stateb" = "ccode2", "year")) %>%
  left_join(nt_ %>% rename(tota = tot), by = c("statea" = "ccode", "year")) %>% left_join(nt_ %>% rename(totb = tot), by = c("stateb" = "ccode", "year"))
stopifnot(nrow(d) == n0)
d <- d %>% mutate(salience = sqrt((trade / tota) * (trade / totb)), jdem = (polity2a + 10) * (polity2b + 10), ldcap = log(capdist + 1),
  outcome_year = year + 1)
d$salience[!is.finite(d$salience)] <- NA
rhs0 <- c("salience", "cap_a", "cap_b", "anyally", "contigland", "jdem", "lcaprat", "ldcap", "py", "spline1", "spline2", "spline3")
vars <- c("mzmid1_orig", rhs0, "dyadid")
cat("dyad-years 1869-1991 (outcome 1870-1992):", nrow(d), " missing:\n"); print(sapply(d[, vars], function(x) sum(is.na(x))))
d0 <- d[complete.cases(d[, vars]), ]
cat("complete-case N:", nrow(d0), "(published 137,314)  onsets", sum(d0$mzmid1_orig), "\n")
TERM <- "salience"
fitter <- function(r, data) fit_bin("mzmid1_orig", r, data, "logit", "hetero", terms = TERM)
full <- fitter(rhs0, d0)
cat("\nREBUILT vs published: salience b", full[[TERM]][["b"]], "(pub -3.176) se", full[[TERM]][["se"]], "(pub 0.870) N", attr(full, "n"), "(pub 137,314)\n")
ds <- common_sample(d, vars, "statea", "stateb", "outcome_year"); rm(d, d0); invisible(gc())
M <- run_six(rhs0, ds, fitter)
if (Sys.getenv("SWEEP_SENS") == "1") {
  source(file.path(CODE, "hist/sens_helper.R"))
  sens_run(id = "P9_081", dat = ds, y = "mzmid1_orig", tie = TERM, ctrl = setdiff(rhs0, TERM), fe = NULL,
           fit = function(dat, add) fitter(c(rhs0, add), dat)[[TERM]])
}
saveRDS(list(M = M, full = full), file.path(W, "M_P9_081.rds"))
meta <- list(pid = "P9_081", cite = "Xiang, Xu & Keteku 2007", journal = "JCR", year = 2007, doi = "10.1177/0022002707302802", family = "TIE")
nt <- sprintf(paste0("REBUILT (no archive): DOR/O-R all dyads, IVs 1869-1991, DV mzmid1_orig (onset t+1); salience from COW dyadic and national trade; ",
  "CINC of each state; anyally, land contiguity, joint democracy (polity2+10 product), log CINC ratio, log capital distance, py + 3 splines; HC SE. ",
  "Full-sample rebuilt b=%.6f se=%.6f N=%d vs published -3.176 (0.870) N 137,314. Common sample N=%d."),
  full[[TERM]][["b"]], full[[TERM]][["se"]], attr(full, "n"), nrow(ds))
write_p8(meta, TERM, "abstract: with both states' power controlled, trade interdependence constrains conflict; Table 1 Model 2",
  "Table 1 Model 2 (logit, MID onset 1870-1992, with Power_1 and Power_2)", "logit, robust SE, peace years + 3 splines",
  "feglm logit, HC vcov, all six steps", -3.176, 0.870, 137314, TRUE, TRUE, M, "Y (peace years + splines)", nt)
dir.create(file.path(OUT, "results_p9"), recursive = TRUE, showWarnings = FALSE)
src <- file.path(OUT, "results_p8", "P9_081__salience.csv"); stopifnot(file.exists(src), file.rename(src, file.path(OUT, "results_p9", basename(src))))
sink()
