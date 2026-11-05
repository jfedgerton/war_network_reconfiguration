# P8_005 Beyene 2015 Foreign Trade Review 50(2) "Does International Trade Reduce Political Disputes?" (10.1177/0015732515572058)
# REBUILD with peacesciencer (no replication archive). Spec from the staged PDF:
#   Table 3, Equation 8 (full model, pp.109-111): logit, DV = MID onset ("1 if a MID of any type occurred between dyads"),
#   8 East African states (Djibouti 522, Eritrea 531, Ethiopia 530, Kenya 501, Somalia 520, Sudan 625, Tanzania 510, Uganda 500),
#   "All Dyads, 1991-2011", N = 1175, all IVs lagged one year:
#   OpennessLow/High(Log) = log((X_i+M_i)/GDP_i) lower/higher in dyad; TradeShareLow/High = dyadic trade/(X_i+M_i);
#   DependenceLow/High = dyadic trade/GDP_i; DemocracyLow/High = polity2 lower/higher; Alliance (defence, entente or neutrality);
#   National Capability(Log) = log(CINC_stronger/CINC_weaker). No peace-years / temporal terms. SEs = SPSS model-based (iid).
#   Published Eq 8: DependenceHigh -11.272 (6.077), p = .064; DependenceLow -13.229 (49.885), ns.
# KEY TERM (fixed before estimating): trade_dependence = DependenceHigh(t-1). Abstract: "the pacification effect of ... economic
#   dependence on trade ... is significant"; the only dependence term the paper reports as significant (Eqs 6-8) is DependenceHigh.
# Unit: the paper's N = 1175 matches DIRECTED dyad-years (8 states, Eritrea from 1993, 1991-2011: 1148-1190 rows), so directed
#   dyad-years are used (flag: all regressors are symmetric, so each undirected dyad-year enters twice, as in the paper's N).
# Not replicable exactly: IMF DOTS trade -> COW dyadic trade (cow_trade_ddy, flow1+flow2) and COW national trade (cow_trade_sy);
#   GDP from peacesciencer cow_sdp_gdp (exp(wbgdp2011est), 2011 USD) vs current-USD trade (unit mismatch only rescales ratios);
#   Polity IV (Marshall 2013) -> polity2 from ccode_democracy (NOT add_democracy(), known NA bug); MID onset -> GML onset.
CODE <- Sys.getenv("SWEEP_CODE"); stopifnot(nzchar(CODE))
source(file.path(CODE, "work_p8/p8_common.R"))
suppressPackageStartupMessages(library(peacesciencer))
WD <- file.path(ROOT, "work_p8/P8_005"); W <- file.path(OUT, "work_p8/P8_005"); dir.create(W, recursive = TRUE, showWarnings = FALSE); sink(file.path(W, "log_P8_005.txt"), split = TRUE)
EA <- c(500, 501, 510, 520, 522, 530, 531, 625)
dy <- create_dyadyears(directed = TRUE, subset_years = 1991:2011) %>% filter(ccode1 %in% EA, ccode2 %in% EA) %>% add_gml_mids() %>%
  select(ccode1, ccode2, year, gmlmidonset) %>% as.data.frame()
cat("directed dyad-years 1991-2011:", nrow(dy), " onsets:", sum(dy$gmlmidonset), "\n"); stopifnot(nrow(dy) > 1000, nrow(dy) < 1250)
# state-year inputs (t-1 values joined on year - 1)
sy_tr <- cow_trade_sy %>% transmute(ccode, year, tot = ifelse(imports < 0 | exports < 0, NA, imports + exports))
# fallback where COW national trade is missing (Eritrea entirely): sum of COW dyadic flows with ALL partners (flagged in notes)
ftr <- system.file("extdata", "cow_trade_ddy.rds", package = "peacesciencer"); stopifnot(nzchar(ftr), file.exists(ftr))
dsum <- readRDS(ftr) %>% filter(ccode1 %in% EA) %>%
  mutate(f = ifelse(flow1 < 0 | flow2 < 0, NA, flow1 + flow2)) %>% group_by(ccode = ccode1, year) %>%
  summarise(tot_dsum = if (all(is.na(f))) NA_real_ else sum(f, na.rm = TRUE), .groups = "drop")
sy_tr <- sy_tr %>% full_join(dsum, by = c("ccode", "year")) %>% mutate(tot_fallback = is.na(tot) & !is.na(tot_dsum), tot = coalesce(tot, tot_dsum)) %>%
  select(ccode, year, tot, tot_fallback)
print(sy_tr %>% filter(ccode %in% EA, year >= 1990, year <= 2010) %>% group_by(ccode) %>% summarise(n = n(), fallback = sum(tot_fallback), na = sum(is.na(tot))))
gdp <- cow_sdp_gdp %>% transmute(ccode, year, gdp = exp(wbgdp2011est) / 1e6)      # millions 2011 USD (trade in millions current USD)
pol <- ccode_democracy %>% select(ccode, year, polity2) %>% as.data.frame(); stopifnot(!anyDuplicated(pol[, 1:2]))
nmc <- cow_nmc %>% select(ccode, year, cinc)
stp <- sy_tr %>% full_join(gdp, by = c("ccode", "year")) %>% full_join(pol, by = c("ccode", "year")) %>% full_join(nmc, by = c("ccode", "year")) %>%
  filter(ccode %in% EA) %>% as.data.frame()
stopifnot(!anyDuplicated(stp[, c("ccode", "year")]))
tr <- readRDS(ftr) %>% filter(ccode1 %in% EA, ccode2 %in% EA) %>%
  mutate(flow1 = ifelse(flow1 < 0, NA, flow1), flow2 = ifelse(flow2 < 0, NA, flow2), dtrade = flow1 + flow2) %>%
  select(ccode1, ccode2, year, dtrade) %>% as.data.frame()
stopifnot(!anyDuplicated(tr[, 1:3]))
ca <- peacesciencer::cow_alliance %>% mutate(ally = as.numeric(coalesce(cow_defense, 0) + coalesce(cow_neutral, 0) + coalesce(cow_entente, 0) > 0)) %>%
  select(ccode1, ccode2, year, ally)
ca2 <- rbind(ca, ca %>% rename(ccode1 = ccode2, ccode2 = ccode1)) %>% distinct(ccode1, ccode2, year, .keep_all = TRUE)   # both orders
d <- dy %>% mutate(yl = year - 1) %>%
  left_join(stp %>% rename_with(~paste0(.x, "_1"), -c(ccode, year)), by = c("ccode1" = "ccode", "yl" = "year")) %>%
  left_join(stp %>% rename_with(~paste0(.x, "_2"), -c(ccode, year)), by = c("ccode2" = "ccode", "yl" = "year")) %>%
  left_join(tr, by = c("ccode1", "ccode2", "yl" = "year")) %>%
  left_join(ca2, by = c("ccode1", "ccode2", "yl" = "year")) %>%
  mutate(ally = coalesce(ally, 0),
         open1 = log(tot_1 / gdp_1), open2 = log(tot_2 / gdp_2), OpenLow = pmin(open1, open2), OpenHigh = pmax(open1, open2),
         sh1 = dtrade / tot_1, sh2 = dtrade / tot_2, ShareLow = pmin(sh1, sh2), ShareHigh = pmax(sh1, sh2),
         dep1 = dtrade / gdp_1, dep2 = dtrade / gdp_2, DepLow = pmin(dep1, dep2), DepHigh = pmax(dep1, dep2),
         DemLow = pmin(polity2_1, polity2_2), DemHigh = pmax(polity2_1, polity2_2),
         lncap = log(pmax(cinc_1, cinc_2) / pmin(cinc_1, cinc_2)), trade_dependence = DepHigh) %>% as.data.frame()
stopifnot(nrow(d) == nrow(dy))
rhs0 <- c("OpenLow", "OpenHigh", "ShareLow", "ShareHigh", "DepLow", "trade_dependence", "DemLow", "DemHigh", "ally", "lncap")
vars <- c("gmlmidonset", rhs0)
for (v in rhs0) cat(sprintf("%-18s NA %4d  mean %9.4g  sd %9.4g\n", v, sum(!is.finite(d[[v]])), mean(d[[v]][is.finite(d[[v]])]), sd(d[[v]][is.finite(d[[v]])])))
for (v in rhs0) d[[v]][!is.finite(d[[v]])] <- NA          # log(0) etc. -> missing (listwise, as SPSS)
fitter <- function(rhs, data) fit_bin("gmlmidonset", rhs, data, "logit", "iid", terms = c("trade_dependence", "DepLow"))
full <- d[complete.cases(d[, vars]), ]
cat("rebuilt complete-case N =", nrow(full), "(published 1175)\n")
r_full <- fitter(rhs0, full); n_orig <- attr(r_full, "n")
ds <- common_sample(d, vars, "ccode1", "ccode2", "year")
M <- run_six(rhs0, ds, fitter)
# ---- sensitivity hook: runs only under run_sensitivity.R (SWEEP_SENS=1), never in the main sweep ----
if (Sys.getenv("SWEEP_SENS") == "1") {
  source(file.path(CODE, "hist/sens_helper.R"))
  sens_run(id = "P8_005", dat = ds, y = "gmlmidonset", tie = "trade_dependence", ctrl = setdiff(rhs0, "trade_dependence"), fe = NULL,
           fit = function(dat, add) { r <- fitter(c(rhs0, add), dat)$trade_dependence; c(b = unname(r[["b"]]), se = unname(r[["se"]]), p = unname(r[["p"]])) })
}
saveRDS(list(M = M, r_full = r_full, n_orig = n_orig, n_used = nrow(ds)), file.path(W, "M_P8_005.rds"))
meta <- list(pid = "P8_005", short_cite = "Beyene 2015", journal = "FTR", year = 2015, doi = "10.1177/0015732515572058", family = "TIE",
  headline_table_model = "Table 3 Equation 8 (full model; DependenceHigh t-1)",
  article_estimator = "logit (SPSS), model-based SE, IVs lagged 1 yr, no temporal terms; 8 East African states 1991-2011, N=1175",
  est_used_all_steps = "logit (fixest feglm), iid SE, all six steps; directed dyad-years (matches published N)",
  author_has_history_control = "N (no peace years or MID history)")
write_p8(meta, "trade_dependence", "abstract: 'the pacification effect of economic openness and economic dependence on trade ... is significant'",
  -11.272, 6.077, 1175, FALSE, TRUE, M,
  paste0("REBUILT (peacesciencer; no archive). Key = DependenceHigh(t-1), the only dependence term reported significant (Eq 6 ***, Eq 8 * at .10). ",
         "Published p = .064 (not significant at .05). Rebuilt full-sample N=", n_orig, " (pub 1175); full-sample DependenceHigh b=",
         signif(r_full$trade_dependence["b"], 4), " p=", signif(r_full$trade_dependence["p"], 3), "; DependenceLow b=", signif(r_full$DepLow["b"], 4),
         " p=", signif(r_full$DepLow["p"], 3), " (pub -13.229 (49.885), ns). Substitutions: IMF DOTS -> COW trade (dyadic + national), ",
         "SDP GDP (2011 USD), polity2 via ccode_democracy, GML MID onset (GML ends 2010, so DV year 2011 is lost). Eritrea has no COW national trade: its total trade = sum of COW dyadic flows with all partners (fallback).  Directed dyad-years: symmetric regressors, each dyad twice (as the published N implies)."),
  full = r_full)
sink()
