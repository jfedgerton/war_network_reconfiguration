# P8_043 Francis 2009 DPE "The Human Capital Peace" (10.1080/10242690802682109) -- REBUILD (no replication data found).
# Paper (author copy, Francis Google site): Table 4 (E): logit of conflict (dyad involved in a MID with use of force or war) on Human capital L/H,
#   politically relevant dyads (contiguous incl. <=400 mi water, or major power), 1960-2000; all IVs lagged 1 yr and standardized (except
#   contiguity/alliance); country fixed effects (delta_i + delta_j), year dummies, peace-spell dummies; robust SE.
#   Published: Human capital L = -0.962 (0.293), Human capital H = -0.894 (0.262), N = 13,653.
#   Human capital = one-factor score of fertility, primary/secondary/tertiary gross enrollment (WDI, linearly interpolated).
# KEY TERMS (a priori, abstract: 'human capital may significantly decrease the likelihood of militarized conflict'): hc_lo, hc_hi.
# Rebuild choices (documented, not tuned): WDI current vintage (not WDI 2002); GML MIDs (ongoing, hostlev >= 4) instead of MID 3.10;
#   COW trade / larger SDP GDP (not Gleditsch 2002 / PWT); SDP GDPpc (not PWT 6.1); polity2 from ccode_democracy; alliance = any COW alliance
#   or both allied with US; capability ratio from COW NMC; peace spells from GML conflict series since 1816 (dummies = factor(py));
#   country FE as multi-membership dummies (1 if country is either side); countries whose dummy perfectly predicts 0 dropped (as Stata).
#   Population/UN-membership screen of the paper not applied.
CODE <- Sys.getenv("SWEEP_CODE"); stopifnot(nzchar(CODE))
source(file.path(CODE, "work_p8/p8_common.R"))
suppressMessages({library(peacesciencer); library(countrycode); library(tidyr)})
WD <- file.path(ROOT, "work_p8/P8_043"); W <- file.path(OUT, "work_p8/P8_043"); dir.create(W, recursive = TRUE, showWarnings = FALSE); sink(file.path(W, "log_P8_043.txt"), split = TRUE)

# 1. human capital factor score (country-year)
rd <- function(ind) {
  f <- unzip(file.path(WD, paste0(ind, ".zip")), list = TRUE)$Name; f <- f[grepl("^API_", f)]
  x <- read.csv(unz(file.path(WD, paste0(ind, ".zip")), f), skip = 4, check.names = FALSE)
  x <- x[, c("Country Code", as.character(1950:2005)[as.character(1950:2005) %in% names(x)])]
  x <- pivot_longer(x, -`Country Code`, names_to = "year", values_to = "v"); x$year <- as.integer(x$year); names(x)[1] <- "iso3"; x
}
hc <- rd("SP.DYN.TFRT.IN") %>% rename(fert = v) %>% full_join(rd("SE.PRM.ENRR") %>% rename(prim = v), by = c("iso3", "year")) %>%
  full_join(rd("SE.SEC.ENRR") %>% rename(sec = v), by = c("iso3", "year")) %>% full_join(rd("SE.TER.ENRR") %>% rename(ter = v), by = c("iso3", "year"))
hc$ccode <- suppressWarnings(countrycode(hc$iso3, "iso3c", "cown"))
hc <- hc[!is.na(hc$ccode), ]
ip <- function(y, v) if (sum(!is.na(v)) >= 2) approx(y, v, xout = y, rule = 1)$y else v     # linear interpolation, no extrapolation
hc <- hc %>% arrange(ccode, year) %>% group_by(ccode, iso3) %>%
  mutate(fert = ip(year, fert), prim = ip(year, prim), sec = ip(year, sec), ter = ip(year, ter)) %>% ungroup()
hc <- hc %>% group_by(ccode, year) %>% summarise(across(c(fert, prim, sec, ter), ~ mean(.x, na.rm = TRUE)), .groups = "drop") %>%
  mutate(across(c(fert, prim, sec, ter), ~ ifelse(is.nan(.x), NA, .x)))
cc <- complete.cases(hc[, c("fert", "prim", "sec", "ter")])
fa <- factanal(hc[cc, c("fert", "prim", "sec", "ter")], factors = 1, scores = "regression")
print(fa$loadings)
hc$hcap <- NA; hc$hcap[cc] <- as.numeric(fa$scores)
if (cor(hc$hcap, hc$sec, use = "complete.obs") < 0) hc$hcap <- -hc$hcap          # orient: higher = more human capital
cat("hcap country-years:", sum(!is.na(hc$hcap)), " years:", range(hc$year[!is.na(hc$hcap)]), "\n")

# 2. dyad-years 1816-2000 for DV / peace spells; PRD 1959-2000 for covariates
nd <- create_dyadyears(directed = FALSE, subset_years = 1816:2000) %>% add_gml_mids(keep = c("hostlev")) %>% as.data.frame()
nd$conflict <- as.numeric(nd$gmlmidongoing == 1 & !is.na(nd$hostlev) & nd$hostlev >= 4)
nd <- nd %>% arrange(ccode1, ccode2, year) %>% group_by(ccode1, ccode2) %>%
  mutate(spell = cumsum(lag(conflict, default = 0) == 1)) %>% group_by(ccode1, ccode2, spell) %>% mutate(py = row_number() - 1) %>%
  ungroup() %>% select(ccode1, ccode2, year, conflict, py) %>% as.data.frame()

pol <- ccode_democracy %>% select(ccode, year, polity2) %>% as.data.frame()
ftr <- system.file("extdata", "cow_trade_ddy.rds", package = "peacesciencer"); stopifnot(nzchar(ftr), file.exists(ftr))
tr <- readRDS(ftr) %>%
  mutate(flow1 = ifelse(flow1 < 0, NA, flow1), flow2 = ifelse(flow2 < 0, NA, flow2), lo = pmin(ccode1, ccode2), hi = pmax(ccode1, ccode2)) %>%
  group_by(lo, hi, year) %>% summarise(trade = mean(flow1 + flow2), .groups = "drop")
ca <- peacesciencer::cow_alliance %>% mutate(any = as.numeric(coalesce(cow_defense, 0) + coalesce(cow_neutral, 0) + coalesce(cow_nonagg, 0) + coalesce(cow_entente, 0) > 0))
usal <- ca %>% filter(ccode1 == 2 | ccode2 == 2, any == 1) %>% transmute(ccode = ifelse(ccode1 == 2, ccode2, ccode1), year, us = 1) %>% distinct()
ca2 <- ca %>% transmute(lo = pmin(ccode1, ccode2), hi = pmax(ccode1, ccode2), year, any) %>% group_by(lo, hi, year) %>% summarise(ally_f = max(any), .groups = "drop")
ctl <- create_dyadyears(directed = FALSE, subset_years = 1959:2000) %>% add_contiguity() %>% add_cow_majors() %>% filter_prd() %>%
  add_capital_distance() %>% add_nmc() %>% {suppressWarnings(add_sdp_gdp(.))} %>% as.data.frame()
cat("PRD dyad-years 1959-2000:", nrow(ctl), "\n")
ctl <- ctl %>% left_join(pol %>% rename(p1 = polity2), by = c("ccode1" = "ccode", "year")) %>% left_join(pol %>% rename(p2 = polity2), by = c("ccode2" = "ccode", "year")) %>%
  left_join(tr, by = c("ccode1" = "lo", "ccode2" = "hi", "year")) %>% left_join(ca2, by = c("ccode1" = "lo", "ccode2" = "hi", "year")) %>%
  left_join(usal %>% rename(us1 = us), by = c("ccode1" = "ccode", "year")) %>% left_join(usal %>% rename(us2 = us), by = c("ccode2" = "ccode", "year")) %>%
  left_join(hc %>% select(ccode, year, h1 = hcap), by = c("ccode1" = "ccode", "year")) %>% left_join(hc %>% select(ccode, year, h2 = hcap), by = c("ccode2" = "ccode", "year")) %>%
  mutate(contig = as.numeric(conttype %in% 1:5), dist = capdist,
         dem_lo = pmin(p1, p2), dem_hi = pmax(p1, p2),
         alliance = as.numeric(coalesce(ally_f, 0) == 1 | (coalesce(us1, 0) == 1 & coalesce(us2, 0) == 1) | (ccode1 == 2 & coalesce(us2, 0) == 1) | (ccode2 == 2 & coalesce(us1, 0) == 1)),
         btrade = coalesce(trade, 0) * 1e6 / pmax(exp(wbgdp2011est1), exp(wbgdp2011est2)),
         capratio = log(pmax(cinc1, cinc2) / pmin(cinc1, cinc2)),
         gdppc_lo = pmin(wbgdppc2011est1, wbgdppc2011est2), gdppc_hi = pmax(wbgdppc2011est1, wbgdppc2011est2),
         hc_lo = pmin(h1, h2), hc_hi = pmax(h1, h2)) %>%
  select(ccode1, ccode2, year, contig, dist, dem_lo, dem_hi, alliance, btrade, capratio, gdppc_lo, gdppc_hi, hc_lo, hc_hi)
# 3. lag all IVs one year; DV and peace spell at t
ivs <- c("contig", "dist", "dem_lo", "dem_hi", "alliance", "btrade", "capratio", "gdppc_lo", "gdppc_hi", "hc_lo", "hc_hi")
lagd <- ctl %>% mutate(year = year + 1)
d <- nd %>% filter(year >= 1960, year <= 2000) %>% inner_join(lagd, by = c("ccode1", "ccode2", "year"))
d <- d[is.finite(d$capratio) | is.na(d$capratio), ]
stz <- c("dist", "dem_lo", "dem_hi", "btrade", "capratio", "gdppc_lo", "gdppc_hi", "hc_lo", "hc_hi")
cc0 <- complete.cases(d[, c("conflict", ivs)])
for (v in stz) { m <- mean(d[[v]][cc0]); s <- sd(d[[v]][cc0]); d[[v]] <- (d[[v]] - m) / s }   # standardized on the complete-case sample
cat("complete-case dyad-years:", sum(cc0), " conflicts:", sum(d$conflict[cc0]), "(published N 13,653 after FE drops)\n")
# 4. multi-membership country dummies; iteratively drop countries whose dummy perfectly predicts failure (Stata logit behaviour)
d <- d[cc0, ]
repeat {
  cs <- sort(unique(c(d$ccode1, d$ccode2)))
  anyc <- sapply(cs, function(c) any(d$conflict[d$ccode1 == c | d$ccode2 == c] == 1))
  if (all(anyc)) break
  bad <- cs[!anyc]; d <- d[!(d$ccode1 %in% bad | d$ccode2 %in% bad), ]
}
cs <- sort(unique(c(d$ccode1, d$ccode2)))
for (c in cs[-1]) d[[paste0("cD", c)]] <- as.numeric(d$ccode1 == c | d$ccode2 == c)     # first country = reference
cdum <- paste0("cD", cs[-1])
cat("after country perfect-prediction drops: N =", nrow(d), " countries:", length(cs), "\n")
rhs <- c("hc_lo", "hc_hi", "contig", "dist", "dem_lo", "dem_hi", "alliance", "btrade", "capratio", "gdppc_lo", "gdppc_hi", cdum)
terms <- c("hc_lo", "hc_hi")
fitter <- function(r, data) {
  o <- fit_bin("conflict", r, data, "logit", "hetero", terms, fe = "year + py")
  o
}
# print only key rows of fit (country dummies are many) -> wrap
fit_quiet <- function(r, data) { out <- NULL; capture.output(out <- fitter(r, data)); for (t in terms) cat(t, round(out[[t]], 4), " N", attr(out, "n"), "\n"); out }
full <- fit_quiet(rhs, d)
ds <- common_sample(d, c("conflict", rhs, "py"), "ccode1", "ccode2", "year")
M <- run_six(rhs, ds, fit_quiet)
meta <- list(pid = "P8_043", short_cite = "Francis 2009", journal = "DPE", year = 2009, doi = "10.1080/10242690802682109",
  family = "CAPITALISM", headline_table_model = "Table 4, column E (logit, country FE, year and peace-spell dummies)",
  article_estimator = "logit with country fixed effects (delta_i + delta_j), year dummies, peace-spell dummies, robust SE; IVs lagged, standardized",
  est_used_all_steps = "logit; multi-membership country dummies + year FE + peace-year FE (fixest); heteroskedasticity-robust SE; all six steps",
  author_has_history_control = "N (peace-spell dummies only)")
nt <- "REBUILT (no replication data): see script header for every substitution (WDI current vintage factor score, GML use-of-force MID involvement, COW trade/SDP GDP, polity2, COW alliance + US-ally rule). "
write_p8(meta, "hc_lo", "Abstract: 'human capital may significantly decrease the likelihood of militarized conflict' (Table 4 E, HC L)",
  -0.962, 0.293, 13653, TRUE, TRUE, M, nt, full = full)
write_p8(meta, "hc_hi", "Abstract: 'human capital may significantly decrease the likelihood of militarized conflict' (Table 4 E, HC H)",
  -0.894, 0.262, 13653, TRUE, TRUE, M, nt, full = full)
sink()
