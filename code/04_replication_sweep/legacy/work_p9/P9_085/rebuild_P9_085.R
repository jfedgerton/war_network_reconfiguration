# P9_085 Mousseau 2000, JCR 44(4) "Market Prosperity, Democratic Consolidation, and Democratic Peace" (10.1177/0022002700044004004).
#   CAPITALISM. REBUILD (protocol B; PARTIAL; no archive).
# Article (staged scanned PDF, OCR, pp. 483-490, Table 1): logit with cubic splines (BKT), Oneal-Russett (1999) all-dyad data 1950-92,
#   DV = involvement in a MID in a given year (O/R coding); all covariates lagged one year. Democracy low = Polity III (democ - autoc) + 11 of
#   the less democratic state (1-21); development low/high = log PWT 5.6 real GDP per capita of the less/more prosperous state (missing filled
#   from previous years up to 3); developed democracy = development low x democracy low; trade interdependence, relative capability, major power
#   status, contiguity, distance (O/R), hegemonic power, hegemonic status. One-tailed tests.
# HEADLINE (fixed before estimating): Table 1 Model 3, Developed democracy -0.02 (0.01), t -2.7, *** (p < .01 one-tailed), N 114,598.
# PARTIAL: hegemonic power and hegemonic status (custom controls) are DROPPED (stated deviation); PWT 5.6 replaced by the {peacesciencer}
#   simulated GDP per capita (mrgdppc); Polity IV 2014 polity2 + 11 for Polity III; GML MID ongoing for O/R involvement.
# Rebuild inputs: {peacesciencer} COW dyad-years 1946-1992 (GML MID ongoing, contiguity 1-4, capital distance, CINC, majors, COW trade,
#   sim GDP); Polity IV 2014 (work_p8/P8_049). DV in t+1 against covariates in t (the one-year lag); history at outcome_year.
ROOT <- Sys.getenv("SWEEP_DATA"); OUT <- Sys.getenv("SWEEP_OUT"); CODE <- Sys.getenv("SWEEP_CODE")
stopifnot(nzchar(ROOT), nzchar(OUT), nzchar(CODE))
source(file.path(CODE, "work_p8/common_b3.R")); source(file.path(CODE, "work_p9/rebuild_common_p9r.R"))
set.seed(123)
W <- file.path(OUT, "work_p9/P9_085"); dir.create(W, recursive = TRUE, showWarnings = FALSE)
sink(file.path(W, "log_P9_085.txt"), split = TRUE)
d <- ps_dyads(1946:1992)
d <- d %>% mutate(dyadid = 1000 * ccode1 + ccode2, inv = as.integer(gmlmidongoing == 1))
lead <- d %>% transmute(ccode1, ccode2, year = year - 1, inv_t1 = inv)
d <- d %>% left_join(lead, by = c("ccode1", "ccode2", "year"))
d <- py_calc(d, "dyadid", "inv")
d <- d[d$year >= 1949 & d$year <= 1991, ]
# GDPpc: fill missing from previous years up to 3 (article)
s <- suppressMessages(create_stateyears(system = "cow", subset_years = 1940:1992) %>% add_sim_gdp_pop()) %>% as.data.frame() %>% select(ccode, year, gpc = mrgdppc) %>%
  arrange(ccode, year) %>% group_by(ccode) %>% mutate(gpcf = coalesce(gpc, dplyr::lag(gpc, 1), dplyr::lag(gpc, 2), dplyr::lag(gpc, 3))) %>% ungroup() %>% select(ccode, year, gpcf)
p <- polity_cow("polity2"); n0 <- nrow(d)
d <- d %>% left_join(p %>% rename(p2a = polity2), by = c("ccode1" = "ccode", "year")) %>% left_join(p %>% rename(p2b = polity2), by = c("ccode2" = "ccode", "year")) %>%
  left_join(s %>% rename(ga = gpcf), by = c("ccode1" = "ccode", "year")) %>% left_join(s %>% rename(gb = gpcf), by = c("ccode2" = "ccode", "year"))
stopifnot(nrow(d) == n0)
d <- d %>% mutate(demlow = pmin(p2a, p2b) + 11, devlow = log(pmin(ga, gb)), devhigh = log(pmax(ga, gb)), devdem = devlow * demlow,
  tdep = pmin((flow1 + flow2) / pwtrgdp1, (flow1 + flow2) / pwtrgdp2), caprat = log(pmax(cinc1, cinc2) / pmin(cinc1, cinc2)),
  major = as.integer(cowmaj1 == 1 | cowmaj2 == 1), contig = as.integer(conttype %in% 1:4), ldist = log(capdist + 1), outcome_year = year + 1)
for (v in c("tdep", "caprat", "devlow", "devhigh")) d[[v]][!is.finite(d[[v]])] <- NA
d <- btscs_splines(d)
rhs0 <- c("demlow", "devhigh", "devlow", "devdem", "tdep", "caprat", "major", "contig", "ldist", "py", "spl1", "spl2", "spl3")
vars <- c("inv_t1", rhs0, "dyadid")
cat("dyad-years 1949-91 (outcome 1950-92):", nrow(d), " missing:\n"); print(sapply(d[, vars], function(x) sum(is.na(x))))
d0 <- d[complete.cases(d[, vars]), ]
cat("complete-case N:", nrow(d0), "(published 114,598)  involvement years", sum(d0$inv_t1), "(published events 856)\n")
TERM <- "devdem"
fitter <- function(r, data) fit_bin("inv_t1", r, data, "logit", "iid", terms = TERM)
full <- fitter(rhs0, d0)
cat("\nREBUILT vs published: devdem b", full[[TERM]][["b"]], "(pub -0.02) se", full[[TERM]][["se"]], "(pub 0.01) N", attr(full, "n"), "(pub 114,598)\n")
ds <- common_sample(d, vars, "ccode1", "ccode2", "outcome_year"); rm(d, d0); invisible(gc())
M <- run_six(rhs0, ds, fitter)
if (Sys.getenv("SWEEP_SENS") == "1") {
  source(file.path(CODE, "hist/sens_helper.R"))
  sens_split(id = "P9_085", dat = ds, fit = function(dat, add) fitter(c(rhs0, add), dat)[[TERM]])
}
saveRDS(list(M = M, full = full), file.path(W, "M_P9_085.rds"))
meta <- list(pid = "P9_085", cite = "Mousseau 2000", journal = "JCR", year = 2000, doi = "10.1177/0022002700044004004", family = "CAPITALISM")
nt <- sprintf(paste0("REBUILT, PARTIAL (no archive): peacesciencer all dyads, covariates 1949-91, DV = GML MID involvement in t+1 (1950-92); Polity IV 2014 polity2+11 low; ",
  "log sim GDPpc low/high (3-yr carry-forward); product; lower COW trade/PWT GDP; log CINC ratio; major; contiguity 1-4; log capital distance; ",
  "btscs splines; conventional SE. Hegemonic power/status DROPPED. Full-sample rebuilt b=%.6f se=%.6f N=%d vs published -0.02 (0.01) N 114,598. ",
  "Common N=%d. Interaction: constituents kept in all steps."), full[[TERM]][["b"]], full[[TERM]][["se"]], attr(full, "n"), nrow(ds))
write_p8(meta, TERM, "abstract: the democratic peace is limited to / stronger among market-developed democracies; Table 1 Model 3",
  "Table 1 Model 3 (logit, MID involvement, all dyads 1950-92; developed democracy)", "logit with cubic splines, one-tailed tests",
  "feglm logit, iid vcov, all six steps", -0.02, 0.01, 114598, TRUE, TRUE, M, "Y (cubic splines)", nt)
dir.create(file.path(OUT, "results_p9"), recursive = TRUE, showWarnings = FALSE)
src <- file.path(OUT, "results_p8", "P9_085__devdem.csv"); stopifnot(file.exists(src), file.rename(src, file.path(OUT, "results_p9", basename(src))))
sink()
