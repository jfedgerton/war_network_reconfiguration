# VERIFY 78: Benson 2005 CMPS -- copy of rebuild_78.R matched to the PDF (Table 5, Sample 14, p.129).
# Paper: duration-dependent logit of MID onset, full sample 1948-1992, Barbieri data; IVs at t-1: Salience, Symmetry, Dyadsize,
#   Salience*Dyadsize (KEY), Contiguity (land or <=400 mi water), Joint democracy ((dem-aut+10)/2 product, Polity III), Alliance
#   (defense/neutrality/entente), ln CINC ratio; peace years + 3 BTSCS splines; robust SE clustered by dyad (p.123).
#   Published Sample 14: Salience*Dyadsize -0.019 (0.011), p = .078; N = 118,964.
# Only change vs rebuild_78.R: polity2 now joined from peacesciencer::ccode_democracy (COW codes). add_democracy() (v1.2.0) returned
#   NA polity2 for USSR 1923-91, Yugoslavia, Italy, etc. (Polity-code keys), silently dropping those dyads, which the paper includes.
# Not fixable: Barbieri trade/GDP (COW v4 trade nominal / SDP GDP 2011$ -> magnitude not comparable); Barbieri's missing-trade
#   coverage (N ~2x); BTSCS knots unreported (cubic kept); MID v2.1 onset -> GML onset.
suppressMessages({library(peacesciencer)})
ROOT <- Sys.getenv("SWEEP_DATA"); OUT <- Sys.getenv("SWEEP_OUT"); CODE <- Sys.getenv("SWEEP_CODE")
stopifnot(nzchar(ROOT), nzchar(OUT), nzchar(CODE))
source(file.path(CODE, "vb4/vhelp.R"))
WO <- file.path(OUT, "work/78"); dir.create(WO, recursive = TRUE, showWarnings = FALSE)  # port: outputs
W <- file.path(ROOT, "work/78"); sink(file.path(WO, "log_verify_78.txt"), split = TRUE)

nd <- create_dyadyears(directed = FALSE, subset_years = 1816:1992) %>% add_gml_mids()
nd <- nd %>% arrange(ccode1, ccode2, year) %>% group_by(ccode1, ccode2) %>%
  mutate(spell = cumsum(lag(gmlmidonset, default = 0) == 1)) %>% group_by(ccode1, ccode2, spell) %>%
  mutate(py = row_number() - 1) %>% ungroup() %>% filter(year >= 1948) %>% select(ccode1, ccode2, year, gmlmidonset, py)
f_tr <- system.file("extdata", "cow_trade_ddy.rds", package = "peacesciencer"); stopifnot(nzchar(f_tr), file.exists(f_tr))
tr <- readRDS(f_tr) %>%
  mutate(flow1 = ifelse(flow1 < 0, NA, flow1), flow2 = ifelse(flow2 < 0, NA, flow2),
         lo = pmin(ccode1, ccode2), hi = pmax(ccode1, ccode2)) %>%
  group_by(lo, hi, year) %>% summarise(trade = mean(flow1 + flow2), .groups = "drop")
pol <- ccode_democracy %>% select(ccode, year, polity2) %>% as.data.frame(); stopifnot(!anyDuplicated(pol[, 1:2]))
iv <- create_dyadyears(directed = FALSE, subset_years = 1947:1991) %>%
  add_contiguity() %>% add_nmc() %>% {suppressWarnings(add_sdp_gdp(.))}
iv <- iv %>% left_join(pol %>% rename(polity21 = polity2), by = c("ccode1" = "ccode", "year")) %>%
  left_join(pol %>% rename(polity22 = polity2), by = c("ccode2" = "ccode", "year"))
ca <- peacesciencer::cow_alliance
iv <- iv %>% left_join(ca, by = c("ccode1","ccode2","year")) %>%
  left_join(tr, by = c("ccode1" = "lo", "ccode2" = "hi", "year" = "year")) %>%
  mutate(gdp1 = exp(wbgdp2011est1) / 1e6, gdp2 = exp(wbgdp2011est2) / 1e6,
         dep1 = trade / gdp1, dep2 = trade / gdp2,
         salience = sqrt(dep1 * dep2), symmetry = 1 - abs(dep1 - dep2),
         dyadsize = (gdp1 + gdp2) / 1e3, sal_x_size = salience * dyadsize,
         contig = as.numeric(conttype %in% 1:5),
         jointdem = ((polity21 + 10) / 2) * ((polity22 + 10) / 2),
         alliance = as.numeric(coalesce(cow_defense, 0) + coalesce(cow_neutral, 0) + coalesce(cow_entente, 0) > 0),
         relcap = log(pmax(cinc1, cinc2) / pmin(cinc1, cinc2)), year = year + 1) %>%
  select(ccode1, ccode2, year, salience, symmetry, dyadsize, sal_x_size, contig, jointdem, alliance, relcap)
d <- nd %>% inner_join(iv, by = c("ccode1","ccode2","year")) %>%
  mutate(py2 = py^2 / 10, py3 = py^3 / 100, dyad = ccode1 * 1000 + ccode2) %>% as.data.frame()
rm(nd, iv, tr); invisible(gc())
rhs0 <- c("sal_x_size","salience","symmetry","dyadsize","contig","jointdem","alliance","relcap","py","py2","py3")
d <- d[is.finite(d$relcap) & d$year >= 1948 & d$year <= 1992, ]
FF <- function(r, s) fit_bin("gmlmidonset", r, s, "logit", "~dyad", NULL, "sal_x_size")
x <- d[complete.cases(d[, c("gmlmidonset", rhs0, "dyad")]), ]
rep <- FF(rhs0, x); cat("rebuilt N =", attr(rep, "n"), "(pub 118,964; previous rebuild 243,377)\n"); rm(x)
s <- common_sample(d, c("gmlmidonset", rhs0, "dyad"), "ccode1", "ccode2", "year"); rm(d); invisible(gc())
M <- run_six(rhs0, s, FF)
# ---- sensitivity hook (split only): runs only under run_sensitivity.R (SWEEP_SENS=1), never in the main sweep ----
if (Sys.getenv("SWEEP_SENS") == "1") {
  source(file.path(CODE, "hist/sens_helper.R"))
  sens_split(id = "78", dat = s,
             fit = function(dat, add) { r <- FF(c(rhs0, add), dat)$sal_x_size; c(b = unname(r[["b"]]), se = unname(r[["se"]]), p = unname(r[["p"]])) })
}
write_verify(78, "Benson 2005", "TIE", "sal_x_size", "Table 5 Sample 14 (Barbieri data, p.129)",
  -0.019, 0.011, 118964, "logit, robust SE clustered by dyad, peace years + 3 BTSCS splines", FALSE,
  "Barbieri trade replaced by COW v4 trade/SDP GDP (nominal/real units); N ~2x (trade coverage); cubic peace years instead of BTSCS splines; GML onset instead of MID onset; add_democracy() NA polity2 for USSR/Yugoslavia/Italy etc. dropped those dyads",
  "PARTIAL", "Barbieri trade/GDP data and coverage unavailable (units differ -> magnitude not comparable); BTSCS knots unreported",
  TRUE, M, sprintf("REBUILD (no author files). Full-sample rebuilt b=%.4g p=%.3g N=%d. pct_vs_pub not meaningful (trade units differ).",
                   rep$sal_x_size["b"], rep$sal_x_size["p"], attr(rep, "n")))
sink()
