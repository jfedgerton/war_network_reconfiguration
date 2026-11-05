# P9_069 Hegre 2009, CMPS 26(1) "Trade Dependence or Size Dependence? The Gravity Model of Trade and the Liberal Peace"
#   (10.1177/0738894208097665). TIE. REBUILD (protocol B; no archive).
# Article (staged PDF pp. 31-37, Table 3): logit, robust SEs, all dyads 1950-1992 from the Oneal-Russett (2005) data; DV = onset of a fatal
#   MID (Table 1 C-models and the log-likelihoods of Table 3); i = the larger economy (GDP) in the dyad; zero trade set to USD 200,000
#   (fn. 11). Model C3: ln(Tij/Gi) + ln Pi + ln Pj + ln Ii + ln Ij + Democracy_i + Democracy_j + Democracy_ij + Cij + ln Dij + Brevity of
#   peace + ln Nt. DV coded from the text ('onset of MIDs that lead to at least one battle-death'), all dyadic participants
#   (mzfatal1_all; no originator restriction is stated); brevity of peace counts years since any MID onset (mzmid1_all).
# HEADLINE (fixed before estimating): Table 3 column B (Model C3), ln(Tij/Gi) trade dependence = -0.243 (0.0282) ***, N 443,751.
# Rebuild inputs: Dafoe-Oneal-Russett all-dyad file (work_demcap/D1/dor.dta: mzfatal1_all, mzmid1_all) + derived O/R trade/GDP file
#   (work_p4/X074/or_trade.rds: trade $m, GDP $, polity2, numstate, contigland); population = COW NMC v4 tpop; capital distance from
#   {peacesciencer}. Brevity of peace = 2^(-y/5), y = years since the last MID onset (Raknerud-Hegre decay; half-life not stated in the
#   article - assumption documented). Democracy_ij = polity2_i x polity2_j.
ROOT <- Sys.getenv("SWEEP_DATA"); OUT <- Sys.getenv("SWEEP_OUT"); CODE <- Sys.getenv("SWEEP_CODE")
stopifnot(nzchar(ROOT), nzchar(OUT), nzchar(CODE))
source(file.path(CODE, "work_p8/common_b3.R")); source(file.path(CODE, "work_p9/rebuild_common_p9r.R"))
set.seed(123)
W <- file.path(OUT, "work_p9/P9_069"); dir.create(W, recursive = TRUE, showWarnings = FALSE)
sink(file.path(W, "log_P9_069.txt"), split = TRUE)
d <- or_dyads(1930:1992)
d <- d %>% mutate(dyadid = 1000 * statea + stateb)
d <- py_calc(d, "dyadid", "mzmid1_all"); d <- d[d$year >= 1950, ]
n <- read.csv(file.path(ROOT, "work_p8/P8_049/nmc_v4.csv"))[, c("ccode", "year", "tpop")]; n$tpop[n$tpop < 0] <- NA
cd <- suppressMessages(create_dyadyears(system = "cow", directed = FALSE, subset_years = 1950:1992) %>% add_capital_distance()) %>% as.data.frame()
n0 <- nrow(d)
d <- d %>% left_join(n %>% rename(popa = tpop), by = c("statea" = "ccode", "year")) %>% left_join(n %>% rename(popb = tpop), by = c("stateb" = "ccode", "year")) %>%
  left_join(cd[, c("ccode1", "ccode2", "year", "capdist")], by = c("statea" = "ccode1", "stateb" = "ccode2", "year"))
stopifnot(nrow(d) == n0)
d <- d %>% mutate(ai = gdp_a >= gdp_b,
  Gi = ifelse(ai, gdp_a, gdp_b), Pi = ifelse(ai, popa, popb), Pj = ifelse(ai, popb, popa), Gj = ifelse(ai, gdp_b, gdp_a),
  Di = ifelse(ai, polity2a, polity2b), Dj = ifelse(ai, polity2b, polity2a),
  T = ifelse(trade == 0, 0.2, trade) * 1e6,
  ltdep = log(T / Gi), lPi = log(Pi), lPj = log(Pj), lIi = log(Gi / Pi), lIj = log(Gj / Pj), Dij = Di * Dj,
  C = as.integer(contigland == 1), lD = log(capdist + 1), brev = 2^(-py / 5), lNt = log(numstate))
for (v in c("ltdep", "lPi", "lPj", "lIi", "lIj", "lD")) d[[v]][!is.finite(d[[v]])] <- NA
rhs0 <- c("ltdep", "lPi", "lPj", "lIi", "lIj", "Di", "Dj", "Dij", "C", "lD", "brev", "lNt")
vars <- c("mzfatal1_all", rhs0, "dyadid")
cat("dyad-years 1950-92:", nrow(d), " missing:\n"); print(sapply(d[, vars], function(x) sum(is.na(x))))
d0 <- d[complete.cases(d[, vars]), ]
cat("complete-case N:", nrow(d0), "(published 443,751)  fatal onsets", sum(d0$mzfatal1_all), "\n")
TERM <- "ltdep"
fitter <- function(r, data) fit_bin("mzfatal1_all", r, data, "logit", "hetero", terms = TERM)
full <- fitter(rhs0, d0)
cat("\nREBUILT vs published: ltdep b", full[[TERM]][["b"]], "(pub -0.243) se", full[[TERM]][["se"]], "(pub 0.0282) N", attr(full, "n"), "(pub 443,751)\n")
ds <- common_sample(d, vars, "statea", "stateb", "year"); rm(d, d0); invisible(gc())
M <- run_six(rhs0, ds, fitter)
if (Sys.getenv("SWEEP_SENS") == "1") {
  source(file.path(CODE, "hist/sens_helper.R"))
  sens_run(id = "P9_069", dat = ds, y = "mzfatal1_all", tie = TERM, ctrl = setdiff(rhs0, TERM), fe = NULL,
           fit = function(dat, add) fitter(c(rhs0, add), dat)[[TERM]])
}
saveRDS(list(M = M, full = full), file.path(W, "M_P9_069.rds"))
meta <- list(pid = "P9_069", cite = "Hegre 2009", journal = "CMPS", year = 2009, doi = "10.1177/0738894208097665", family = "TIE")
nt <- sprintf(paste0("REBUILT (no archive): DOR/O-R all dyads 1950-92 (dor.dta + derived or_trade.rds), fatal MID onset (mzfatal1_all); i = larger GDP; ",
  "zero trade = $0.2m; NMC population; peacesciencer capital distance; brevity of peace 2^(-y/5) (half-life assumed); polity2 product. ",
  "Full-sample rebuilt b=%.6f se=%.6f N=%d vs published -0.243 (0.0282) N 443,751. Common 1822-2011 sample N=%d."),
  full[[TERM]][["b"]], full[[TERM]][["se"]], attr(full, "n"), nrow(ds))
write_p8(meta, TERM, "abstract: correct specification (size controls) keeps a pacifying trade-dependence effect; Table 3 Model C3",
  "Table 3 column B (Model C3; fatal MID onset, all dyads 1950-92)", "logit, robust SE, brevity-of-peace decay, system size",
  "feglm logit, HC vcov, all six steps", -0.243, 0.0282, 443751, TRUE, TRUE, M, "Y (brevity of peace decay)", nt)
dir.create(file.path(OUT, "results_p9"), recursive = TRUE, showWarnings = FALSE)
src <- file.path(OUT, "results_p8", "P9_069__ltdep.csv"); stopifnot(file.exists(src), file.rename(src, file.path(OUT, "results_p9", basename(src))))
sink()
