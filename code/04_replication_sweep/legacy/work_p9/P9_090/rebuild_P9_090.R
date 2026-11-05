# P9_090 Goenner 2011, CMPS 28(5) "Simultaneity between Trade and Conflict: Endogenous Instruments of Mass Destruction"
#   (10.1177/0738894211418414). TIE. REBUILD (protocol B; no archive). Estimator: Newey (1987) two-step IV probit (work_p9/ivprobit_p9.R).
# Article (staged PDF pp. 470-473, Table 2): conflict equation, IV probit (Newey; Lee overidentification test), robust SEs; non-directed dyads
#   (Oneal 2003 / Oneal et al. 2003 data; KPR bilateral trade), n 119,508. Endogenous: LN bilateral trade; instruments: lagged LN bilateral trade
#   and the asymmetry of the dyad's land-to-population ratios (absolute difference). Controls: democracy score (L), growth rate (L), LN capability
#   ratio, allies, contiguity, LN distance, major power, peace years + 3 splines, LN real GDP (H).
# HEADLINE (fixed before estimating): Table 2 conflict equation, LN Bilateral Trade -0.0401 (0.0160) *, n 119,508. Legend *p<.05 **p<.01.
# Rebuild inputs: DOR/O-R all-dyad file + derived trade/GDP file (mzmid1_orig = MID onset in t+1 against t covariates; or_trade: trade $m,
#   gdp, lowdem, lowgrowth, lcaprat, anyally, contig150, majpow, py, spline1-3); LN trade = log(trade + 1) ($m; KPR real trade not available);
#   capital distance ({peacesciencer}); land area = WDI AG.LND.TOTL.K2 (state mean 1961-2000; work_p9/P9_090/wdi_land_P9_090.csv), population =
#   COW NMC tpop. Years 1951-1992 (lagged trade). SEs clustered on dyad.
ROOT <- Sys.getenv("SWEEP_DATA"); OUT <- Sys.getenv("SWEEP_OUT"); CODE <- Sys.getenv("SWEEP_CODE")
stopifnot(nzchar(ROOT), nzchar(OUT), nzchar(CODE))
source(file.path(CODE, "work_p8/common_b3.R")); source(file.path(CODE, "work_p9/rebuild_common_p9r.R")); source(file.path(CODE, "work_p9/ivprobit_p9.R"))
set.seed(123)
W <- file.path(OUT, "work_p9/P9_090"); dir.create(W, recursive = TRUE, showWarnings = FALSE)
sink(file.path(W, "log_P9_090.txt"), split = TRUE)
ivprobit_selftest()
d <- or_dyads(1950:1992)
lt <- d %>% transmute(statea, stateb, year = year + 1, ltrade_l = log(trade + 1))
cd <- suppressMessages(create_dyadyears(system = "cow", directed = FALSE, subset_years = 1950:1992) %>% add_capital_distance()) %>% as.data.frame()
land <- read.csv(file.path(ROOT, "work_p9/P9_090/wdi_land_P9_090.csv"))
land$ccode <- suppressWarnings(countrycode::countrycode(land$iso3c, "iso3c", "cown")); land <- land[!is.na(land$ccode), ] %>% group_by(ccode) %>% summarise(land = mean(land_km2), .groups = "drop")
n <- read.csv(file.path(ROOT, "work_p8/P8_049/nmc_v4.csv"))[, c("ccode", "year", "tpop")]; n$tpop[n$tpop <= 0] <- NA
n0 <- nrow(d)
d <- d %>% left_join(lt, by = c("statea", "stateb", "year")) %>% left_join(cd[, c("ccode1", "ccode2", "year", "capdist")], by = c("statea" = "ccode1", "stateb" = "ccode2", "year")) %>%
  left_join(land %>% rename(landa = land), by = c("statea" = "ccode")) %>% left_join(land %>% rename(landb = land), by = c("stateb" = "ccode")) %>%
  left_join(n %>% rename(popa = tpop), by = c("statea" = "ccode", "year")) %>% left_join(n %>% rename(popb = tpop), by = c("stateb" = "ccode", "year"))
stopifnot(nrow(d) == n0)
d <- d %>% filter(year >= 1951) %>% mutate(ltrade = log(trade + 1), lgdph = log(pmax(gdp_a, gdp_b)), ldcap = log(capdist + 1),
  lpasym = abs(landa / popa - landb / popb), outcome_year = year + 1)
X1n <- c("lowdem", "lowgrowth", "lcaprat", "anyally", "contig150", "ldcap", "majpow", "py", "spline1", "spline2", "spline3", "lgdph")
X2n <- c("ltrade_l", "lpasym")
vars <- c("mzmid1_orig", "ltrade", X1n, X2n, "dyadid")
cat("dyad-years 1951-92:", nrow(d), " missing:\n"); print(sapply(d[, vars], function(x) sum(is.na(x))))
TERM <- "ltrade"
fitter <- function(r, data) {
  extra <- setdiff(r, "BASE")
  f <- ivprobit_newey(data$mzmid1_orig, data$ltrade, as.matrix(data[, c(X1n, extra)]), as.matrix(data[, X2n]), cluster = data$dyadid)
  z <- f$coef / f$se
  cat("\n---- IV probit (Newey two-step) | extra:", paste(extra, collapse = " + "), "| N =", f$n, "\n"); print(round(cbind(b = f$coef, se = f$se, z = z), 6))
  out <- list(c(b = f$coef[["y2"]], se = f$se[["y2"]], p = 2 * pnorm(-abs(z[["y2"]])))); names(out) <- TERM; attr(out, "n") <- f$n; out
}
d0 <- d[complete.cases(d[, vars]), ]
cat("complete-case N:", nrow(d0), "(published 119,508)  onsets", sum(d0$mzmid1_orig), "\n")
full <- fitter("BASE", d0)
cat("\nREBUILT vs published: ltrade b", full[[TERM]][["b"]], "(pub -0.0401) se", full[[TERM]][["se"]], "(pub 0.0160) N", attr(full, "n"), "(pub 119,508)\n")
ds <- common_sample(d, vars, "statea", "stateb", "outcome_year"); rm(d, d0); invisible(gc())
M <- run_six("BASE", ds, fitter)
saveRDS(list(M = M, full = full), file.path(W, "M_P9_090.rds"))
b0 <- M$m0[[TERM]][["b"]]; p0 <- M$m0[[TERM]][["p"]]; r_ <- b0 / -0.0401
v <- if (sign(b0) == -1 && p0 < .05 && r_ >= 0.5 && r_ <= 2) "REBUILT_MATCH" else "REBUILT_MISMATCH"
meta <- list(pid = "P9_090", cite = "Goenner 2011", journal = "CMPS", year = 2011, doi = "10.1177/0738894211418414", family = "TIE")
nt <- sprintf(paste0("REBUILT (no archive): DOR/O-R all dyads, IVs 1951-91 (DV mzmid1_orig onset t+1); Newey two-step IV probit (work_p9/ivprobit_p9.R), ",
  "instruments lagged log trade and land/pop asymmetry (WDI land, NMC pop); log(trade+1) for KPR real trade; SEs clustered on dyad. Full-sample rebuilt ",
  "b=%.6f se=%.6f N=%d vs published -0.0401 (0.0160) N 119,508. Common N=%d. Round-2 rule."), full[[TERM]][["b"]], full[[TERM]][["se"]], attr(full, "n"), nrow(ds))
write_p8(meta, TERM, "abstract and Table 2: with relevant and exogenous instruments, trade reduces conflict",
  "Table 2 conflict equation (Newey IV probit; MID, non-directed dyads)", "IV probit (Newey 1987 two-step), robust SE",
  "hand-coded Newey two-step IV probit, cluster-robust components, all six steps", -0.0401, 0.0160, 119508, TRUE, TRUE, M, "Y (peace years + splines)", nt, verdict = v)
dir.create(file.path(OUT, "results_p9"), recursive = TRUE, showWarnings = FALSE)
src <- file.path(OUT, "results_p8", "P9_090__ltrade.csv"); stopifnot(file.exists(src), file.rename(src, file.path(OUT, "results_p9", basename(src))))
sink()
