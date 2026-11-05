# P8_004 Aaronson, Abouharb & Wang 2015 Global Economy Journal 15(4) "The Liberal Illusion Is Not a Complete Delusion: The WTO
# Helps Member States Keep the Peace Only When It Increases Trade" (10.1515/gej-2014-0049)
# REBUILD (no replication archive). Spec from the staged PDF (pp.469-475; Table 2 read from the rendered page image, p.475):
#   Table 2 Model 4 (MIDs): logit, all non-directed dyad-years 1950-2001, DV = MID onset in year t; regressors: GATT/WTO x Trade,
#   Joint GATT/WTO membership, DemocracyL, DemocracyH (Polity IV), Annual change in GDP per capita (lower-growth state), Trade
#   ("sum of country i's imports from and exports to country j, as a percentage of country i's GDP"), GDPH, GDPL (PWT), Alliance
#   (COW defence/neutrality/entente), ln distance, Major power (at least one), Free Trade Agreements (count of signed dyadic FTAs,
#   Mansfield-Pevehouse updated); peace years + "natural spline function with three knots"; robust SEs; N = 540,991.
#   Published Model 4: GATT/WTO x Trade = -0.01* (SE printed "-0.004", i.e. .004); Joint GATT/WTO 1.08* (.46); Trade 0.004 (0.003).
# KEY TERM (fixed in phase 8; abstract: "when states are both members of the GATT/WTO and benefit from increased trade, they are
#   less likely to engage in militarized interstate disputes"): wto_membership_x_trade = JointGATTWTO x Trade.
# Rebuild choices (paper silent or data unavailable; none tuned):
#   - GATT/WTO membership years: Wikipedia "Member states of the WTO" table (WTO accession + GATT membership dates; downloaded
#     wiki_wto.html) plus hard-coded GATT contracting parties that are not WTO members: Czechoslovakia (315) 1948-1992 and
#     Yugoslavia (345) 1966-1992; Germany's 1951 date applied to West Germany (260) and to Germany (255). Member in year t if joined
#     in or before t. GATT withdrawals (ROC 1950, Lebanon/Syria 1951, Liberia 1953) are not coded.
#   - Trade = 100 x dyadic trade / GDP of the LOWER-dependence state (weak link; "country i" is undefined in non-directed data);
#     dyadic trade and nominal GDP from Gleditsch v4.1 (expanded trade, pwt_cow; IMF DOTS/PWT-based, 1950-2000).
#   - GDPH/GDPL = real GDP (Gleditsch rgdp96pc x pop, millions 1996 USD); growth = % change in real GDP per capita, lower of the two.
#   - IVs contemporaneous (the paper does not say they are lagged); Gleditsch data end in 2000, so DV-year 2001 is lost.
#   - FTA count = number of DESTA base PTAs signed on or before t with both states as members (substitute for Mansfield-Pevehouse).
#   - Peace years + cubic (py, py^2, py^3) instead of the unreported natural-spline knots. SEs clustered by dyad ("robust").
#   - Polity: polity2 via ccode_democracy (NOT add_democracy()). MID onset: GML (paper: COW MID).
CODE <- Sys.getenv("SWEEP_CODE"); stopifnot(nzchar(CODE))
source(file.path(CODE, "work_p8/p8_common.R"))
suppressPackageStartupMessages({library(peacesciencer); library(countrycode)})
WD <- file.path(ROOT, "work_p8/P8_004"); W <- file.path(OUT, "work_p8/P8_004"); dir.create(W, recursive = TRUE, showWarnings = FALSE); G <- file.path(ROOT, "work_p8/P8_012"); sink(file.path(W, "log_P8_004.txt"), split = TRUE)
# ---- membership ----
wg <- read.csv(file.path(WD, "wto_gatt_dates.csv"), stringsAsFactors = FALSE)
wg$ccode <- suppressWarnings(countrycode(wg$name, "country.name", "cown"))
wg$ccode[wg$name == "Serbia"] <- NA
print(wg[is.na(wg$ccode), "name"])
wg$join <- pmin(ifelse(is.na(wg$gatt_year), Inf, wg$gatt_year), ifelse(is.na(wg$wto_year), Inf, wg$wto_year))
mem <- wg %>% filter(!is.na(ccode), is.finite(join)) %>% select(ccode, join)
mem <- bind_rows(mem, data.frame(ccode = 260, join = 1951)) %>% group_by(ccode) %>% summarise(join = min(join), .groups = "drop")
memyr <- expand.grid(ccode = mem$ccode, year = 1950:2001) %>% left_join(mem, by = "ccode") %>% mutate(m = as.numeric(year >= join)) %>%
  select(ccode, year, m) %>% bind_rows(data.frame(ccode = 315, year = 1950:1992, m = 1), data.frame(ccode = 345, year = 1966:1992, m = 1)) %>%
  group_by(ccode, year) %>% summarise(m = max(m), .groups = "drop") %>% as.data.frame()
cat("GATT/WTO members 1950/1970/1990/2001:", sapply(c(1950, 1970, 1990, 2001), function(y) sum(memyr$m[memyr$year == y])), "\n")
# ---- FTA counts from DESTA (base treaties) ----
de <- read.csv(file.path(ROOT, "work_p8/P8_015/desta_list_of_treaties_02_02.txt"), stringsAsFactors = FALSE)
de <- de[de$entry_type == "base_treaty", ]
cc <- grep("^c[0-9]+$", names(de), value = TRUE)
fta <- do.call(rbind, lapply(seq_len(nrow(de)), function(i) {
  m <- unique(na.omit(as.numeric(unlist(de[i, cc])))); m <- m[m > 0]
  if (length(m) < 2) return(NULL)
  cw <- suppressWarnings(countrycode(m, "iso3n", "cown")); cw <- unique(na.omit(cw)); if (length(cw) < 2) return(NULL)
  p <- t(combn(sort(cw), 2)); data.frame(ccode1 = p[, 1], ccode2 = p[, 2], sy = de$year[i])
}))
cat("DESTA base treaties:", nrow(de), " dyad-treaty rows:", nrow(fta), "\n")
# ---- Gleditsch trade/GDP ----
gg <- read.table(file.path(G, "pwt_cow.asc"), header = TRUE) %>% transmute(ccode = statenum, year, gdp = gdppc * pop / 1000, rgdp = rgdp96pc * pop / 1000, rgdppc = rgdp96pc) %>%
  arrange(ccode, year) %>% group_by(ccode) %>% mutate(growth = 100 * (rgdppc / lag(rgdppc) - 1), growth = ifelse(year - lag(year) == 1, growth, NA)) %>% ungroup() %>% as.data.frame()
stopifnot(!anyDuplicated(gg[, 1:2]))
gt <- read.table(file.path(G, "uddtrade_cow.asc"), header = TRUE) %>% filter(year >= 1950, year <= 2001) %>%
  transmute(ccode1 = pmin(numa, numb), ccode2 = pmax(numa, numb), year, tr = ifelse(expab < 0 | impab < 0, NA, expab + impab)) %>%
  distinct(ccode1, ccode2, year, .keep_all = TRUE)
# ---- dyad-years ----
nd <- create_dyadyears(directed = FALSE, subset_years = 1816:2001) %>% add_gml_mids() %>%
  select(ccode1, ccode2, year, gmlmidonset) %>% arrange(ccode1, ccode2, year) %>%
  group_by(ccode1, ccode2) %>% mutate(spell = cumsum(lag(gmlmidonset, default = 0) == 1)) %>%
  group_by(ccode1, ccode2, spell) %>% mutate(py = row_number() - 1) %>% ungroup() %>% filter(year >= 1950, year <= 2001) %>% as.data.frame()
pol <- ccode_democracy %>% select(ccode, year, polity2) %>% as.data.frame(); stopifnot(!anyDuplicated(pol[, 1:2]))
ca <- peacesciencer::cow_alliance %>% mutate(ally = as.numeric(coalesce(cow_defense, 0) + coalesce(cow_neutral, 0) + coalesce(cow_entente, 0) > 0)) %>%
  select(ccode1, ccode2, year, ally)
ctl <- create_dyadyears(directed = FALSE, subset_years = 1950:2001) %>% add_capital_distance() %>% add_cow_majors() %>%
  select(ccode1, ccode2, year, capdist, cowmaj1, cowmaj2) %>% as.data.frame()
d <- nd %>% left_join(ctl, by = c("ccode1", "ccode2", "year")) %>%
  left_join(pol %>% rename(p1 = polity2), by = c("ccode1" = "ccode", "year")) %>% left_join(pol %>% rename(p2 = polity2), by = c("ccode2" = "ccode", "year")) %>%
  left_join(gg %>% rename_with(~paste0(.x, "1"), -c(ccode, year)), by = c("ccode1" = "ccode", "year")) %>%
  left_join(gg %>% rename_with(~paste0(.x, "2"), -c(ccode, year)), by = c("ccode2" = "ccode", "year")) %>%
  left_join(gt, by = c("ccode1", "ccode2", "year")) %>% left_join(ca, by = c("ccode1", "ccode2", "year")) %>%
  left_join(memyr %>% rename(m1 = m), by = c("ccode1" = "ccode", "year")) %>% left_join(memyr %>% rename(m2 = m), by = c("ccode2" = "ccode", "year"))
stopifnot(nrow(d) == nrow(nd))
fy <- d %>% select(ccode1, ccode2, year) %>% inner_join(fta, by = c("ccode1", "ccode2"), relationship = "many-to-many") %>% filter(sy <= year) %>%
  count(ccode1, ccode2, year, name = "fta")
d <- d %>% left_join(fy, by = c("ccode1", "ccode2", "year")) %>%
  mutate(fta = coalesce(fta, 0L), ally = coalesce(ally, 0), jgatt = as.numeric(coalesce(m1, 0) == 1 & coalesce(m2, 0) == 1),
         DemL = pmin(p1, p2), DemH = pmax(p1, p2), GrowthL = pmin(growth1, growth2),
         TradePct = 100 * pmin(tr / gdp1, tr / gdp2), GDPH = pmax(rgdp1, rgdp2), GDPL = pmin(rgdp1, rgdp2),
         lndist = log(capdist + 1), major = as.numeric(cowmaj1 == 1 | cowmaj2 == 1),
         wto_membership_x_trade = jgatt * TradePct, py2 = py^2, py3 = py^3, dyad = ccode1 * 1000 + ccode2) %>% as.data.frame()
rm(nd, ctl, gt, gg, fy, fta); invisible(gc())
rhs0 <- c("wto_membership_x_trade", "jgatt", "DemL", "DemH", "GrowthL", "TradePct", "GDPH", "GDPL", "ally", "lndist", "major", "fta", "py", "py2", "py3")
vars <- c("gmlmidonset", rhs0, "dyad")
for (v in rhs0) cat(sprintf("%-24s NA %7d  mean %10.4g  sd %10.4g\n", v, sum(is.na(d[[v]])), mean(d[[v]], na.rm = TRUE), sd(d[[v]], na.rm = TRUE)))
fitter <- function(rhs, data) fit_bin("gmlmidonset", rhs, data, "logit", "~dyad", terms = c("wto_membership_x_trade", "jgatt", "TradePct"))
full <- d[complete.cases(d[, vars]), ]
cat("rebuilt complete-case N =", nrow(full), "(published 540,991) years", range(full$year), "\n")
r_full <- fitter(rhs0, full); n_orig <- attr(r_full, "n"); rm(full); invisible(gc())
ds <- common_sample(d, vars, "ccode1", "ccode2", "year"); rm(d); invisible(gc())
M <- run_six(rhs0, ds, fitter)
# ---- sensitivity hook (split only): runs only under run_sensitivity.R (SWEEP_SENS=1), never in the main sweep ----
if (Sys.getenv("SWEEP_SENS") == "1") {
  source(file.path(CODE, "hist/sens_helper.R"))
  sens_split(id = "P8_004", dat = ds,
             fit = function(dat, add) { r <- fitter(c(rhs0, add), dat)$wto_membership_x_trade; c(b = unname(r[["b"]]), se = unname(r[["se"]]), p = unname(r[["p"]])) })
}
saveRDS(list(M = M, r_full = r_full, n_orig = n_orig, n_used = nrow(ds)), file.path(W, "M_P8_004.rds"))
meta <- list(pid = "P8_004", short_cite = "Aaronson et al. 2015", journal = "GEJ", year = 2015, doi = "10.1515/gej-2014-0049", family = "TIE",
  headline_table_model = "Table 2 Model 4 (MIDs; GATT/WTO x Trade)",
  article_estimator = "logit, robust SE, peace years + natural spline (3 knots); all non-directed dyads 1950-2001; N=540,991",
  est_used_all_steps = "logit (fixest feglm), SE clustered by dyad, peace years cubic, all six steps",
  author_has_history_control = "N (peace years + splines only)")
write_p8(meta, "wto_membership_x_trade", "abstract: 'when states are both members of the GATT/WTO and benefit from increased trade, they are less likely to engage in militarized interstate disputes'",
  -0.01, 0.004, 540991, TRUE, TRUE, M,
  paste0("REBUILT (no archive). Product JointGATT/WTO x Trade (Model 4; published -0.01*, SE .004, rounded to 2 decimals). Full-sample rebuilt N=", n_orig,
         " (pub 540,991): product b=", signif(r_full$wto_membership_x_trade["b"], 4), " p=", signif(r_full$wto_membership_x_trade["p"], 3),
         "; joint GATT/WTO b=", signif(r_full$jgatt["b"], 4), " p=", signif(r_full$jgatt["p"], 3), " (pub 1.08*); Trade b=", signif(r_full$TradePct["b"], 4),
         " p=", signif(r_full$TradePct["p"], 3), " (pub 0.004 ns). Substitutions: membership from Wikipedia WTO/GATT dates (+CZE/YUG hard-coded), ",
         "Trade = 100 x lower dyadic trade/GDP (Gleditsch v4.1), real GDP/growth from Gleditsch PWT (1950-2000, so 2001 lost), FTA count from DESTA base treaties, ",
         "cubic peace years for the unreported spline knots, polity2 via ccode_democracy, GML onsets, SE clustered by dyad."),
  full = r_full)
sink()
