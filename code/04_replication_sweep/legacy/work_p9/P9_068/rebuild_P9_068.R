# P9_068 Henderson 2009, BJPS 39(1) "Disturbing the Peace: African Warfare, Political Inversion and the Universality of the Democratic
#   Peace Thesis" (10.1017/S0007123408000495). DEMOCRACY. REBUILD (protocol B; no archive).
# Article (staged PDF pp. 40-44, Table 1): all non-directed dyad-years 1950-2001 (Gleditsch-Ward system), MID onset (subsequent dispute
#   years not coded as onsets), rare-events logit (relogit) with robust SEs and BKT peace years + 3 cubic splines. DemocracyLO = lower
#   Polity score (Gleditsch modified Polity 4); Economic DevelopmentLO = lower log GDP per capita; Allies; Contiguity; Capability ratio;
#   Trade RatioLO = lower bilateral trade/GDP (Gleditsch); Major (either a major power); Distance = log capital distance;
#   Africa = 1 for dyads of two sub-Saharan African states; Africa x DemocracyLO.
# HEADLINE (fixed before estimating): Table 1 column 2, Africa x DemocracyLO = 0.14 (0.02) ***, N 491,569.
# Rebuild inputs: {peacesciencer} COW dyad-years (stand-in for G-W), GML MID onset, COW contiguity (conttype 1-4 = contiguous),
#   capital distance, NMC CINC, COW majors, COW trade (flow1+flow2, current $m) over PWT real GDP ($m) from add_sim_gdp_pop (stand-in for
#   Gleditsch trade/GDP), Maddison-based GDPpc estimate (mrgdppc); Polity IV 2014 polity2 (work_p8/P8_049; stand-in for the modified
#   Polity); COW alliance v4.1. Sub-Saharan Africa = COW codes 402-599 (North Africa 600-626 excluded).
ROOT <- Sys.getenv("SWEEP_DATA"); OUT <- Sys.getenv("SWEEP_OUT"); CODE <- Sys.getenv("SWEEP_CODE")
stopifnot(nzchar(ROOT), nzchar(OUT), nzchar(CODE))
source(file.path(CODE, "work_p8/common_b3.R")); source(file.path(CODE, "work/relogit_kz_ea1.R")); source(file.path(CODE, "work_p9/rebuild_common_p9r.R"))
set.seed(123)
W <- file.path(OUT, "work_p9/P9_068"); dir.create(W, recursive = TRUE, showWarnings = FALSE)
sink(file.path(W, "log_P9_068.txt"), split = TRUE)
d <- ps_dyads(1946:2001) %>% add_cow_ally()
d <- d %>% mutate(dyadid = 1000 * ccode1 + ccode2, mid = as.integer(gmlmidonset == 1))
d <- py_calc(d, "dyadid", "mid"); d <- d[d$year >= 1950, ]
p <- polity_cow("polity2"); n0 <- nrow(d)
d <- d %>% left_join(p %>% rename(p2a = polity2), by = c("ccode1" = "ccode", "year")) %>% left_join(p %>% rename(p2b = polity2), by = c("ccode2" = "ccode", "year"))
stopifnot(nrow(d) == n0)
ssa <- function(cc) cc >= 402 & cc <= 599
d <- d %>% mutate(demlo = pmin(p2a, p2b), devlo = log(pmin(mrgdppc1, mrgdppc2)), allies = cow_any, contig = as.integer(conttype %in% 1:4),
  caprat = log(pmax(cinc1, cinc2) / pmin(cinc1, cinc2)), trlo = pmin((flow1 + flow2) / pwtrgdp1, (flow1 + flow2) / pwtrgdp2),
  major = as.integer(cowmaj1 == 1 | cowmaj2 == 1), ldist = log(capdist + 1), africa = as.integer(ssa(ccode1) & ssa(ccode2)), afr_x_dem = africa * demlo) %>%
  mutate(caprat = ifelse(is.finite(caprat), caprat, NA_real_), trlo = ifelse(is.finite(trlo), trlo, NA_real_))
d <- btscs_splines(d)
rhs0 <- c("demlo", "devlo", "allies", "contig", "caprat", "trlo", "major", "ldist", "africa", "afr_x_dem", "py", "spl1", "spl2", "spl3")
vars <- c("mid", rhs0, "dyadid")
cat("dyad-years 1950-2001:", nrow(d), " missing:\n"); print(sapply(d[, vars], function(x) sum(is.na(x))))
d0 <- d[complete.cases(d[, vars]), ]
cat("complete-case N:", nrow(d0), "(published 491,569)  onsets", sum(d0$mid), "\n")
TERM <- "afr_x_dem"
fitter <- function(r, data) fit_relogit("mid", r, data, terms = TERM, vcov_type = "hetero")
full <- fitter(rhs0, d0)
cat("\nREBUILT vs published: afr_x_dem b", full[[TERM]][["b"]], "(pub 0.14) se", full[[TERM]][["se"]], "(pub 0.02) N", attr(full, "n"), "(pub 491,569)\n")
ds <- common_sample(d, vars, "ccode1", "ccode2", "year"); rm(d, d0); invisible(gc())
M <- run_six(rhs0, ds, fitter)
saveRDS(list(M = M, full = full), file.path(W, "M_P9_068.rds"))
meta <- list(pid = "P9_068", cite = "Henderson 2009", journal = "BJPS", year = 2009, doi = "10.1017/S0007123408000495", family = "DEMOCRACY")
nt <- sprintf(paste0("REBUILT (no archive): peacesciencer COW all dyads 1950-2001 (paper G-W), GML MID onset; Polity IV 2014 (paper Gleditsch modified Polity); ",
  "lower log GDPpc (mrgdppc), lower COW trade/PWT GDP, COW alliance, conttype 1-4, log CINC ratio, major, log capital distance; ",
  "Africa = both COW 402-599. relogit (King-Zeng), robust (HC) SE, peace years + btscs splines. Full-sample rebuilt b=%.6f se=%.6f N=%d vs ",
  "published 0.14 (0.02) N 491,569. Common 1822-2011 sample N=%d. Interaction: constituents kept in all steps."),
  full[[TERM]][["b"]], full[[TERM]][["se"]], attr(full, "n"), nrow(ds))
write_p8(meta, TERM, "abstract: politically open African states more likely to fight each other; Table 1 col 2 product term",
  "Table 1 column 2 (relogit, MID onset, all dyads 1950-2001; Africa x DemocracyLO)", "rare-events logit (relogit), robust SE, peace years + 3 splines",
  "King-Zeng relogit (work/relogit_kz_ea1.R), HC SE, all six steps", 0.14, 0.02, 491569, TRUE, TRUE, M, "Y (peace years + splines)", nt)
dir.create(file.path(OUT, "results_p9"), recursive = TRUE, showWarnings = FALSE)
src <- file.path(OUT, "results_p8", "P9_068__afr_x_dem.csv"); stopifnot(file.exists(src), file.rename(src, file.path(OUT, "results_p9", basename(src))))
sink()
