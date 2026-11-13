# P9_089 Reuveny & Li 2003, ISQ 47(3) "The Joint Democracy-Dyadic Conflict Nexus: A Simultaneous Equations Model" (10.1111/1468-2478.4703002).
#   DEMOCRACY. REBUILD (protocol B; no archive). Estimator: Maddala-type two-stage (2SPLS) system, hand-coded.
# Article (staged PDF pp. 333-341, Tables 1 and 4): politically relevant dyads 1950-92 (Oneal-Russett 1999 data); three equations:
#   DEML, DEMH (OLS; Polity 98 democ - autoc of the less / more democratic state) and MIDAB (probit; Maoz dyadic MID involvement).
#   Stage 1: each endogenous variable on all exogenous variables; stage 2: the endogenous right-hand-side variables replaced by first-stage
#   predictions (linear predictor for MIDAB); SEs = second-stage robust SEs clustered on dyads (no Maddala correction; fn. 25).
#   Exogenous, all lagged one year: democracy equations - prior democracy, year, inflation, inflation x year, GDP pc (PWT 5.6, log),
#   GDP pc x year, growth, trade (openness), diffusion (democracy in the state's politically relevant environment), third-party MID dummy
#   (Maoz); MID equation - contiguity (land or <= 150 miles), alliance, capability ratio (log), trade dependence (lower), growth (lower).
# HEADLINE (fixed before estimating): Table 4 Model 1 (MID equation), DemL -0.0207 (0.0067) ***, N 11,644 (two-tailed; * 10%, ** 5%, *** 1%).
# Rebuild inputs: TRIANGLE.tab (PR dyads; dispute1 = MID involvement, demaut, depend, allies, lcaprat2, noncontg); Polity IV 2014 polity2
#   (work_p8/P8_049); WDI CPI inflation (work_p9/P9_089/wdi_inflation_P9_089.csv; starts 1960, as the article's N implies); {peacesciencer}
#   simulated GDP pc / PWT GDP (stand-in for PWT 5.6), COW national trade 3.0 for openness, COW contiguity (diffusion = mean polity2 of
#   states contiguous within 150 miles), GML MIDs (third-party MID = the state has an ongoing MID with a state other than the dyad partner).
#   History terms (six-step design) enter both stages as exogenous variables.
ROOT <- Sys.getenv("SWEEP_DATA"); OUT <- Sys.getenv("SWEEP_OUT"); CODE <- Sys.getenv("SWEEP_CODE")
stopifnot(nzchar(ROOT), nzchar(OUT), nzchar(CODE))
source(file.path(CODE, "work_p8/common_b3.R")); source(file.path(CODE, "work_p9/rebuild_common_p9r.R"))
set.seed(123)
W <- file.path(OUT, "work_p9/P9_089"); dir.create(W, recursive = TRUE, showWarnings = FALSE)
sink(file.path(W, "log_P9_089.txt"), split = TRUE)
# ---- state-year exogenous variables (1946-1992) ----
yrs <- 1946:1992
s <- suppressMessages(create_stateyears(system = "cow", subset_years = yrs) %>% add_sim_gdp_pop()) %>% as.data.frame() %>% select(ccode, year, mrgdppc, pwtrgdp)
p <- polity_cow("polity2"); nt_ <- cow_nat_trade()
inf <- read.csv(file.path(ROOT, "work_p9/P9_089/wdi_inflation_P9_089.csv"))
inf$ccode <- suppressWarnings(countrycode::countrycode(inf$iso3c, "iso3c", "cown")); inf <- inf[!is.na(inf$ccode), ] %>% group_by(ccode, year) %>% summarise(infl = max(inflation), .groups = "drop")
s <- s %>% left_join(p, by = c("ccode", "year")) %>% left_join(nt_, by = c("ccode", "year")) %>% left_join(inf, by = c("ccode", "year")) %>%
  arrange(ccode, year) %>% group_by(ccode) %>%
  mutate(lgdppc = log(mrgdppc), gr1 = 100 * (mrgdppc / dplyr::lag(mrgdppc) - 1), growth = rowMeans(cbind(gr1, dplyr::lag(gr1), dplyr::lag(gr1, 2)), na.rm = TRUE),
         open = tot / pwtrgdp) %>% ungroup() %>% as.data.frame()
s$growth[is.nan(s$growth)] <- NA; s$open[!is.finite(s$open)] <- NA
dd <- ps_dyads(yrs, directed = TRUE, trade = FALSE) %>% select(ccode1, ccode2, year, conttype, gmlmidongoing)
dif <- dd %>% filter(conttype %in% 1:4) %>% left_join(p %>% rename(pn = polity2), by = c("ccode2" = "ccode", "year")) %>% group_by(ccode = ccode1, year) %>%
  summarise(diffusion = mean(pn, na.rm = TRUE), .groups = "drop") %>% mutate(diffusion = ifelse(is.nan(diffusion), NA, diffusion))
nmid <- dd %>% group_by(ccode = ccode1, year) %>% summarise(nmid = sum(gmlmidongoing == 1, na.rm = TRUE), .groups = "drop")
pair <- dd %>% transmute(ccode1, ccode2, year, pm = as.integer(gmlmidongoing %in% 1)); rm(dd)
s <- s %>% left_join(dif, by = c("ccode", "year")) %>% left_join(nmid, by = c("ccode", "year")) %>% mutate(nmid = coalesce(nmid, 0L))
# ---- TRIANGLE PR dyads, L/H assignment at t, exogenous lagged one year ----
t <- read.delim(file.path(ROOT, "recA/or/TRIANGLE.tab")); stopifnot(nrow(t) == 39996)
t <- t %>% filter(year >= 1950, year <= 1992) %>% mutate(dyadid = 1000 * statea + stateb)
# [2026-11-13 audit] dispute1 in row year t is MID at t+1 (checked vs COW dyadic MID). The article pairs MID_t with DEM_t and exogenous at t-1,
# so the DV is the same dyad's dispute1 from row t-1.
t <- t %>% left_join(read.delim(file.path(ROOT, "recA/or/TRIANGLE.tab")) %>% transmute(statea, stateb, year = year + 1, dispute_t = dispute1), by = c("statea", "stateb", "year"))
n0 <- nrow(t)
lag1 <- s %>% mutate(year = year + 1) %>% select(ccode, year, polity2l = polity2, infl, lgdppc, growth, open, diffusion, nmid)
t <- t %>% left_join(p %>% rename(pa = polity2), by = c("statea" = "ccode", "year")) %>% left_join(p %>% rename(pb = polity2), by = c("stateb" = "ccode", "year")) %>%
  left_join(lag1 %>% rename_with(~ paste0(.x, "_a"), -c(ccode, year)), by = c("statea" = "ccode", "year")) %>%
  left_join(lag1 %>% rename_with(~ paste0(.x, "_b"), -c(ccode, year)), by = c("stateb" = "ccode", "year")) %>%
  left_join(pair %>% mutate(year = year + 1) %>% rename(pm_l = pm), by = c("statea" = "ccode1", "stateb" = "ccode2", "year"))
stopifnot(nrow(t) == n0)
t <- t %>% mutate(aL = pa <= pb, DEML = ifelse(aL, pa, pb), DEMH = ifelse(aL, pb, pa), pm_l = coalesce(pm_l, 0L))
for (v in c("polity2l", "infl", "lgdppc", "growth", "open", "diffusion", "nmid")) {
  t[[paste0(v, "L")]] <- ifelse(t$aL, t[[paste0(v, "_a")]], t[[paste0(v, "_b")]]); t[[paste0(v, "H")]] <- ifelse(t$aL, t[[paste0(v, "_b")]], t[[paste0(v, "_a")]]) }
# MID-equation exogenous (lagged one year from TRIANGLE's own previous row)
tl <- read.delim(file.path(ROOT, "recA/or/TRIANGLE.tab")) %>% transmute(statea, stateb, year = year + 1, dep_l = pmin(dependa, dependb), ally_l = allies,
        caprat_l = lcaprat2, contig_l = -noncontg)
t <- t %>% left_join(tl, by = c("statea", "stateb", "year")); stopifnot(nrow(t) == n0)
t <- t %>% mutate(yr = year - 1950, inflyL = inflL * yr, inflyH = inflH * yr, gdpyL = lgdppcL * yr, gdpyH = lgdppcH * yr,
  thirdL = as.integer(nmidL - pm_l > 0), thirdH = as.integer(nmidH - pm_l > 0), growth_low = pmin(growthL, growthH)) %>% as.data.frame()
demX <- c("polity2lL", "polity2lH", "yr", "inflL", "inflH", "inflyL", "inflyH", "lgdppcL", "lgdppcH", "gdpyL", "gdpyH", "growthL", "growthH", "openL", "openH",
          "diffusionL", "diffusionH", "thirdL", "thirdH")
midX <- c("contig_l", "ally_l", "caprat_l", "dep_l", "growth_low")
allX <- unique(c(demX, midX))
vars <- c("dispute_t", "DEML", "DEMH", allX, "dyadid")
cat("PR dyad-years 1950-92:", nrow(t), " missing:\n"); print(sapply(t[, vars], function(x) sum(is.na(x))))
TERM <- "DEML_hat"
fit2s <- function(data, extra = character(0)) {
  X1 <- paste(c(allX, extra), collapse = " + ")
  data$DEML_hat <- fitted(lm(as.formula(paste("DEML ~", X1)), data = data))
  data$DEMH_hat <- fitted(lm(as.formula(paste("DEMH ~", X1)), data = data))
  fit_bin("dispute_t", c("DEML_hat", "DEMH_hat", midX, extra), data, "probit", "~dyadid", terms = TERM)
}
d0 <- t[complete.cases(t[, vars]), ]
cat("complete-case N:", nrow(d0), "(published 11,644)  MID dyad-years", sum(d0$dispute_t), "\n")
full <- fit2s(d0)
cat("\nREBUILT vs published: DEML_hat b", full[[TERM]][["b"]], "(pub -0.0207) se", full[[TERM]][["se"]], "(pub 0.0067) N", attr(full, "n"), "(pub 11,644)\n")
ds <- common_sample(t, vars, "statea", "stateb", "year")
fitter <- function(r, data) fit2s(data, setdiff(r, "BASE"))
M <- run_six("BASE", ds, fitter)
saveRDS(list(M = M, full = full), file.path(W, "M_P9_089.rds"))
b0 <- M$m0[[TERM]][["b"]]; p0 <- M$m0[[TERM]][["p"]]; r_ <- b0 / -0.0207
v <- rebuild_verdict("P9_089", TERM, b0, p0, -0.0207, NA, 0.0067)   # round-2 rule, alpha_pub from registry.csv
meta <- list(pid = "P9_089", cite = "Reuveny & Li 2003", journal = "ISQ", year = 2003, doi = "10.1111/1468-2478.4703002", family = "DEMOCRACY")
nt <- sprintf(paste0("REBUILT (no archive): TRIANGLE PR dyads 1950-92; two-stage system (OLS for DEML/DEMH on all exogenous variables, probit of MID involvement on fitted ",
  "DEML/DEMH + MID-equation exogenous variables), second-stage SEs clustered on dyad as in the article; WDI inflation, peacesciencer GDP pc, COW trade openness, ",
  "contiguity-based diffusion, GML third-party MIDs. Full-sample rebuilt b=%.6f se=%.6f N=%d vs published -0.0207 (0.0067) N 11,644. Common N=%d. Round-2 rule."),
  full[[TERM]][["b"]], full[[TERM]][["se"]], attr(full, "n"), nrow(ds))
write_p8(meta, TERM, "abstract and Table 4 Model 1: joint democracy reduces dyadic conflict once reciprocity is modelled",
  "Table 4 Model 1 (MID equation of the 2SPLS system; PR dyads 1950-92)", "two-stage (Maddala) system: OLS democracy equations, probit MID equation, robust SE clustered on dyad",
  "hand-coded two-stage system, probit stage 2 clustered on dyad, all six steps", -0.0207, 0.0067, 11644, TRUE, TRUE, M, "N (no temporal terms)", nt, verdict = v)
dir.create(file.path(OUT, "results_p9"), recursive = TRUE, showWarnings = FALSE)
src <- file.path(OUT, "results_p8", "P9_089__DEML_hat.csv"); stopifnot(file.exists(src), file.rename(src, file.path(OUT, "results_p9", basename(src))))
sink()
