# P8_012 Choi 2022/2023 RIPE 30(5) "When does liberal peace fail? Trade and nationalism" (10.1080/09692290.2022.2130959)
# REBUILD (no replication archive). Spec from the staged PDF, pp.1916-1921:
#   Table 1 Model 4 (All MIDs, multiplicative model): logit, all interstate dyads (non-directed), 1900-2001, DV = MID onset
#   (ongoing dyad-years dropped), all IVs lagged one year; Trade (lower dyadic trade/GDP; Gleditsch 2002 post-WWII, Oneal-Russett
#   before), Nationalism (lower of the two states' V-Dem v2exl_legitideolcr_0, 2019 V-Dem), Trade x Nationalism, Democracy (lower
#   Polity score), Geo distance (ln great-circle capital distance), Relative capabilities (ln CINC stronger/weaker), Alliance (COW
#   defence/neutrality/entente), Contiguity ("share a border"), Major power (at least one), System size; t, t^2, t^3; robust SEs
#   (author's other dyadic do-files cluster on dyad -> cluster by dyad). N = 455,963.
#   Published Model 4: Trade x Nationalism = +56.976 (18.055) (text: "significant and positive predictor for (fatal) MIDs");
#   Trade -63.737 (17.145); Nationalism .384 (.170). (Minus signs are lost in the PDF text layer; signs taken from the text.)
# KEY TERM (fixed in phase 8, interaction claim -> product term): trade_x_nationalist_leader = TradeLow x NationalismLow.
# Data used (substitutions flagged): trade/GDP for IV-years 1950-2000 = Gleditsch v4.1 expanded trade (uddtrade_cow, expab+impab,
#   current USD m) / (gdppc x pop/1000, current USD m) (the paper's own post-WWII source). Pre-1950 O&R trade/GDP not available ->
#   rebuilt sample is DV years 1951-2001 only (paper: 1900-2001).
#   Nationalism = V-Dem v16 (vdemdata GitHub) instead of v9 (values differ slightly, e.g. Thatcher .714 vs .667 in the paper).
#   Polity: polity2 via ccode_democracy (NOT add_democracy()). MID onset: GML (Gibler-Miller-Little) instead of COW MID.
#   System size = ln(number of COW system members in year t-1) (the paper does not state the transformation; a raw count would
#   imply a coefficient ~100x smaller than the published .449). Contiguity = direct land contiguity (conttype 1).
CODE <- Sys.getenv("SWEEP_CODE"); stopifnot(nzchar(CODE))
source(file.path(CODE, "work_p8/p8_common.R"))
suppressPackageStartupMessages(library(peacesciencer))
WD <- file.path(ROOT, "work_p8/P8_012"); W <- file.path(OUT, "work_p8/P8_012"); dir.create(W, recursive = TRUE, showWarnings = FALSE); sink(file.path(W, "log_P8_012.txt"), split = TRUE)
# ---- base dyad-years, GML onset, peace years (onset-based, ongoing years dropped) ----
nd <- create_dyadyears(directed = FALSE, subset_years = 1816:2001) %>% add_gml_mids() %>%
  select(ccode1, ccode2, year, gmlmidonset, gmlmidongoing) %>% arrange(ccode1, ccode2, year) %>%
  group_by(ccode1, ccode2) %>% mutate(spell = cumsum(lag(gmlmidonset, default = 0) == 1)) %>%
  group_by(ccode1, ccode2, spell) %>% mutate(py = row_number() - 1) %>% ungroup() %>%
  filter(year >= 1900, year <= 2001) %>% as.data.frame()
stopifnot(!anyDuplicated(nd[, c("ccode1", "ccode2", "year")]))
nd <- nd[!(nd$gmlmidongoing == 1 & nd$gmlmidonset == 0), ]        # drop ongoing, non-onset dyad-years
cat("dyad-years 1900-2001 after dropping ongoing:", nrow(nd), " onsets:", sum(nd$gmlmidonset), "\n")
# ---- lagged dyadic controls (year t-1) ----
pol <- ccode_democracy %>% select(ccode, year, polity2) %>% as.data.frame(); stopifnot(!anyDuplicated(pol[, 1:2]))
vd <- readRDS(file.path(WD, "vdem_nat.rds")) %>% filter(!is.na(COWcode), !is.na(v2exl_legitideolcr_0)) %>%
  group_by(ccode = COWcode, year) %>% summarise(nat = mean(v2exl_legitideolcr_0), .groups = "drop") %>% as.data.frame()
ssz <- create_stateyears(system = "cow", subset_years = 1899:2000) %>% count(year, name = "nstates") %>% mutate(lnsys = log(nstates))
ca <- peacesciencer::cow_alliance %>% mutate(ally = as.numeric(coalesce(cow_defense, 0) + coalesce(cow_neutral, 0) + coalesce(cow_entente, 0) > 0)) %>%
  select(ccode1, ccode2, year, ally)
ctl <- create_dyadyears(directed = FALSE, subset_years = 1899:2000) %>% add_contiguity() %>% add_capital_distance() %>% add_nmc() %>%
  add_cow_majors() %>% select(ccode1, ccode2, year, conttype, capdist, cinc1, cinc2, cowmaj1, cowmaj2) %>%
  left_join(pol %>% rename(pol1 = polity2), by = c("ccode1" = "ccode", "year")) %>% left_join(pol %>% rename(pol2 = polity2), by = c("ccode2" = "ccode", "year")) %>%
  left_join(vd %>% rename(nat1 = nat), by = c("ccode1" = "ccode", "year")) %>% left_join(vd %>% rename(nat2 = nat), by = c("ccode2" = "ccode", "year")) %>%
  left_join(ca, by = c("ccode1", "ccode2", "year")) %>% left_join(ssz, by = "year") %>%
  mutate(ally = coalesce(ally, 0), DemLow = pmin(pol1, pol2), NatLow = pmin(nat1, nat2), lndist = log(capdist + 1),
         lncap = log(pmax(cinc1, cinc2) / pmin(cinc1, cinc2)), contig = as.numeric(conttype == 1), major = as.numeric(cowmaj1 == 1 | cowmaj2 == 1)) %>%
  select(ccode1, ccode2, year, DemLow, NatLow, lndist, lncap, ally, contig, major, lnsys) %>% as.data.frame()
# ---- trade/GDP (lower) at t-1 ----
# Pre-1950 trade/GDP: the local Oneal-Russett JPR file (work/50) covers 1950-92 only and no other pre-1950 nominal trade/GDP series is
#   available, so the rebuild uses IV-years 1950-2000 (DV years 1951-2001) only -- FLAG: sample shorter than the paper's 1900-2001.
gt <- read.table(file.path(WD, "uddtrade_cow.asc"), header = TRUE) %>% filter(year >= 1950, year <= 2000) %>%
  transmute(a = numa, b = numb, year, tr = ifelse(expab < 0 | impab < 0, NA, expab + impab))
gg <- read.table(file.path(WD, "pwt_cow.asc"), header = TRUE) %>% transmute(ccode = statenum, year, gdp = gdppc * pop / 1000)
stopifnot(!anyDuplicated(gg[, 1:2]))
gt <- gt %>% left_join(gg %>% rename(ga = gdp), by = c("a" = "ccode", "year")) %>% left_join(gg %>% rename(gb = gdp), by = c("b" = "ccode", "year")) %>%
  transmute(ccode1 = pmin(a, b), ccode2 = pmax(a, b), year, dep = pmin(tr / ga, tr / gb)) %>% filter(!is.na(dep)) %>%
  distinct(ccode1, ccode2, year, .keep_all = TRUE)
cat("Gleditsch 1950-2000 dyad-years with dependence:", nrow(gt), "\n")
dep <- gt %>% rename(TradeLow = dep)
stopifnot(!anyDuplicated(dep[, 1:3]))
ctl <- ctl %>% left_join(dep, by = c("ccode1", "ccode2", "year")) %>%
  mutate(trade_x_nationalist_leader = TradeLow * NatLow)
d <- nd %>% mutate(yl = year - 1) %>% left_join(ctl %>% rename(yl = year), by = c("ccode1", "ccode2", "yl")) %>%
  mutate(py2 = py^2, py3 = py^3, dyad = ccode1 * 1000 + ccode2) %>% as.data.frame()
rm(nd, ctl, dep, gt, gg); invisible(gc())
rhs0 <- c("TradeLow", "NatLow", "trade_x_nationalist_leader", "DemLow", "lndist", "lncap", "ally", "contig", "major", "lnsys", "py", "py2", "py3")
vars <- c("gmlmidonset", rhs0, "dyad")
for (v in rhs0) cat(sprintf("%-28s NA %7d  mean %10.4g  sd %10.4g\n", v, sum(is.na(d[[v]])), mean(d[[v]], na.rm = TRUE), sd(d[[v]], na.rm = TRUE)))
d <- d[is.finite(d$lncap) | is.na(d$lncap), ]
fitter <- function(rhs, data) fit_bin("gmlmidonset", rhs, data, "logit", "~dyad", terms = c("trade_x_nationalist_leader", "TradeLow", "NatLow"))
full <- d[complete.cases(d[, vars]), ]
cat("rebuilt complete-case N =", nrow(full), "(published 455,963)\n")
r_full <- fitter(rhs0, full); n_orig <- attr(r_full, "n"); rm(full); invisible(gc())
ds <- common_sample(d, vars, "ccode1", "ccode2", "year"); rm(d); invisible(gc())
M <- run_six(rhs0, ds, fitter)
# ---- sensitivity hook (split only): runs only under run_sensitivity.R (SWEEP_SENS=1), never in the main sweep ----
if (Sys.getenv("SWEEP_SENS") == "1") {
  source(file.path(CODE, "hist/sens_helper.R"))
  sens_split(id = "P8_012", dat = ds,
             fit = function(dat, add) { r <- fitter(c(rhs0, add), dat)$trade_x_nationalist_leader; c(b = unname(r[["b"]]), se = unname(r[["se"]]), p = unname(r[["p"]])) })
}
saveRDS(list(M = M, r_full = r_full, n_orig = n_orig, n_used = nrow(ds)), file.path(W, "M_P8_012.rds"))
meta <- list(pid = "P8_012", short_cite = "Choi 2022", journal = "RIPE", year = 2022, doi = "10.1080/09692290.2022.2130959", family = "TIE",
  headline_table_model = "Table 1 Model 4 (All MIDs, Trade x Nationalism)",
  article_estimator = "logit, robust SE, t/t2/t3; non-directed dyads 1900-2001, ongoing years dropped, IVs t-1; N=455,963",
  est_used_all_steps = "logit (fixest feglm), SE clustered by dyad, cubic peace years, all six steps",
  author_has_history_control = "N (t, t2, t3 only)")
write_p8(meta, "trade_x_nationalist_leader", "abstract: 'nationalist leaders are likely to cause the pacifying effect of trade to fall apart'",
  56.976, 18.055, 455963, TRUE, TRUE, M,
  paste0("REBUILT (no archive). Product term TradeLow x NationalismLow (Model 4). Full-sample rebuilt N=", n_orig, " (pub 455,963): product b=",
         signif(r_full$trade_x_nationalist_leader["b"], 4), " p=", signif(r_full$trade_x_nationalist_leader["p"], 3), "; Trade b=",
         signif(r_full$TradeLow["b"], 4), " p=", signif(r_full$TradeLow["p"], 3), " (pub -63.737 (17.145)); Nationalism b=",
         signif(r_full$NatLow["b"], 4), " p=", signif(r_full$NatLow["p"], 3), " (pub .384 (.170)). Sources: Gleditsch v4.1 trade/GDP 1950-2000 (paper's source), ",
         "NO pre-1950 trade data (O&R pre-war series unavailable) -> DV years 1951-2001 only (paper 1900-2001), V-Dem v16 nationalism (paper used v9), polity2 via ccode_democracy, GML onsets, ln system size, land contiguity; ",
         "SE clustered by dyad (paper: 'robust')."),
  full = r_full)
sink()
