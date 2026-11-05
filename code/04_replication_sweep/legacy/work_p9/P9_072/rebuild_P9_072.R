# P9_072 Gartzke & Li 2003b, JPR 40(6) "All's Well That Ends Well: A Reply to Oneal, Barbieri & Peters" (10.1177/00223433030406008). TIE.
#   REBUILD (protocol B; no archive). Same base as P9_072 (Barbieri & Peters' Table I specification, which G&L's Table I re-estimates).
# Article (staged PDF, Table I): Model 3 adds major power and distance to Barbieri & Peters' model (all dyads 1949-92, logit, robust SEs);
#   salience = geometric mean of the two states' trade shares (bilateral trade / state's total trade; 'partner dependence' column);
#   dyadic openness = geometric mean of each state's total trade/GDP; contiguity; joint democracy; alliance; relative capabilities;
#   peace years + 3 splines; G&L data and controls.
# HEADLINE (fixed before estimating): Table I Model 3 (partner dependence, all controls): Salience_t-1 = -0.872 [1.044], n.s., N 118,517
#   (the reply's claim: salience loses significance once major power and distance enter).
# Rebuild inputs: DOR/O-R all-dyad file + derived trade/GDP file (mzmid1_orig = onset in t+1 against covariates at t, i.e. the t-1 lag;
#   history keyed on outcome_year = year + 1); state total trade from COW national trade 3.0; joint democracy = lower Polity (weak link);
#   relative capabilities = log capability ratio; py/splines from or_trade.
ROOT <- Sys.getenv("SWEEP_DATA"); OUT <- Sys.getenv("SWEEP_OUT"); CODE <- Sys.getenv("SWEEP_CODE")
stopifnot(nzchar(ROOT), nzchar(OUT), nzchar(CODE))
source(file.path(CODE, "work_p8/common_b3.R")); source(file.path(CODE, "work_p9/rebuild_common_p9r.R"))
set.seed(123)
W <- file.path(OUT, "work_p9/P9_072"); dir.create(W, recursive = TRUE, showWarnings = FALSE)
sink(file.path(W, "log_P9_072.txt"), split = TRUE)
d <- or_dyads(1948:1991)
cd <- suppressMessages(create_dyadyears(system = "cow", directed = FALSE, subset_years = 1948:1991) %>% add_capital_distance()) %>% as.data.frame()
d <- d %>% left_join(cd[, c("ccode1", "ccode2", "year", "capdist")], by = c("statea" = "ccode1", "stateb" = "ccode2", "year"))
nt_ <- cow_nat_trade(); n0 <- nrow(d)
d <- d %>% left_join(nt_ %>% rename(tota = tot), by = c("statea" = "ccode", "year")) %>% left_join(nt_ %>% rename(totb = tot), by = c("stateb" = "ccode", "year"))
stopifnot(nrow(d) == n0)
d <- d %>% mutate(sha = trade / tota, shb = trade / totb, salience = sqrt(sha * shb), ldcap = log(capdist + 1),
  openness = sqrt((tota * 1e6 / gdp_a) * (totb * 1e6 / gdp_b)), outcome_year = year + 1)
for (v in c("salience", "openness")) d[[v]][!is.finite(d[[v]])] <- NA
rhs0 <- c("salience", "openness", "contig150", "lowdem", "anyally", "lcaprat", "py", "spline1", "spline2", "spline3", "majpow", "ldcap")
vars <- c("mzmid1_orig", rhs0, "dyadid")
cat("dyad-years 1948-91 (outcome 1949-92):", nrow(d), " missing:\n"); print(sapply(d[, vars], function(x) sum(is.na(x))))
d0 <- d[complete.cases(d[, vars]), ]
cat("complete-case N:", nrow(d0), "(published 118,517)  onsets", sum(d0$mzmid1_orig), "\n")
TERM <- "salience"
fitter <- function(r, data) fit_bin("mzmid1_orig", r, data, "logit", "hetero", terms = TERM)
full <- fitter(rhs0, d0)
cat("\nREBUILT vs published: salience b", full[[TERM]][["b"]], "(pub -0.872) se", full[[TERM]][["se"]], "(pub 1.044) N", attr(full, "n"), "(pub 118,517)\n")
ds <- common_sample(d, vars, "statea", "stateb", "outcome_year"); rm(d, d0); invisible(gc())
M <- run_six(rhs0, ds, fitter)
if (Sys.getenv("SWEEP_SENS") == "1") {
  source(file.path(CODE, "hist/sens_helper.R"))
  sens_run(id = "P9_072", dat = ds, y = "mzmid1_orig", tie = TERM, ctrl = setdiff(rhs0, TERM), fe = NULL,
           fit = function(dat, add) fitter(c(rhs0, add), dat)[[TERM]])
}
saveRDS(list(M = M, full = full), file.path(W, "M_P9_072.rds"))
meta <- list(pid = "P9_072", cite = "Gartzke & Li 2003b", journal = "JPR", year = 2003, doi = "10.1177/00223433030406008", family = "TIE")
nt <- sprintf(paste0("REBUILT (no archive): DOR/O-R all dyads, IVs 1948-91, DV mzmid1_orig = MID onset in t+1 (1949-92); salience = sqrt(share_a share_b), ",
  "share = bilateral trade / COW national total trade; openness = geometric mean of total trade/GDP; contig150, lowdem (joint democracy), anyally, ",
  "lcaprat, py + 3 splines, majpow, log capital distance (peacesciencer); robust (HC) SE. Full-sample rebuilt b=%.6f se=%.6f N=%d vs published -0.872 (1.044) N 118,517. Common sample N=%d."),
  full[[TERM]][["b"]], full[[TERM]][["se"]], attr(full, "n"), nrow(ds))
write_p8(meta, TERM, "abstract: B&P's positive salience result is due to omitted major-power and distance controls; Table I Model 3",
  "Table I Model 3 (partner dependence, all controls; logit, MID onset, all dyads 1949-92)", "logit, robust SE, peace years + 3 splines",
  "feglm logit, HC vcov, all six steps", -0.872, 1.044, 118517, FALSE, TRUE, M, "Y (peace years + splines)", nt)
dir.create(file.path(OUT, "results_p9"), recursive = TRUE, showWarnings = FALSE)
src <- file.path(OUT, "results_p8", "P9_072__salience.csv"); stopifnot(file.exists(src), file.rename(src, file.path(OUT, "results_p9", basename(src))))
sink()
