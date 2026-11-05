# P9_067 Mousseau 2005, CMPS 22(1) "Comparing New Theory with Prior Beliefs: Market Civilization and the Democratic Peace"
#   (10.1080/07388940590915327). CAPITALISM (democracy conditional on development). REBUILD (protocol B; no archive).
# Article (staged PDF pp. 68-71, Table 1): land-contiguous non-directed dyad-years 1920-2000; DV = onset of a fatal MID among
#   direct dyadic participants (Maoz DYMID 1.1); DemocracyL = Polity2 of the less democratic state, centred at the median;
#   DevelopmentL = log energy consumption per capita (COW CINC data) of the less developed state, centred at the median;
#   dyad-years with Polity DURABLE = 0 for either state omitted; Capability ratio = log(stronger CINC / weaker CINC);
#   Conflict history = 2^(-y/5), y = years without a militarized dispute. Logit (ML; SEs as printed).
# HEADLINE (fixed before estimating): Table 1 Model 2 (market civilization), DemocracyL x DevelopmentL = -0.527 (0.109),
#   p < .0001; Cases 10,503, events 353.
# Rebuild inputs: {peacesciencer} COW dyad-years 1920-2000 with GML dyadic MIDs (stand-in for DYMID; fatality > 0) and COW direct
#   contiguity (conttype 1 = land); Polity IV 2014 polity2/durable and COW NMC v4 pec/tpop (work_p8/P8_049).
ROOT <- Sys.getenv("SWEEP_DATA"); OUT <- Sys.getenv("SWEEP_OUT"); CODE <- Sys.getenv("SWEEP_CODE")
stopifnot(nzchar(ROOT), nzchar(OUT), nzchar(CODE))
source(file.path(CODE, "work_p8/common_b3.R")); source(file.path(CODE, "work_p9/rebuild_common_p9r.R"))
set.seed(123)
W <- file.path(OUT, "work_p9/P9_067"); dir.create(W, recursive = TRUE, showWarnings = FALSE)
sink(file.path(W, "log_P9_067.txt"), split = TRUE)
d <- ps_dyads(1915:2000, trade = FALSE)
d <- d %>% mutate(dyadid = 1000 * ccode1 + ccode2, anymid = as.integer(gmlmidonset == 1),
                  fmid = as.integer(gmlmidonset == 1 & !is.na(fatality) & fatality > 0))
d <- py_calc(d, "dyadid", "anymid")                         # years without a militarized dispute (any MID onset)
d <- d[d$conttype == 1 & d$year >= 1920, ]
p <- polity_cow(c("polity2", "durable"))
n <- read.csv(file.path(ROOT, "work_p8/P8_049/nmc_v4.csv"))[, c("ccode", "year", "pec", "tpop")]
n <- n %>% mutate(epc = ifelse(pec > 0 & tpop > 0, log(pec / tpop), NA_real_)) %>% select(ccode, year, epc)
n0 <- nrow(d)
d <- d %>% left_join(p %>% rename(p2a = polity2, dura = durable), by = c("ccode1" = "ccode", "year")) %>%
  left_join(p %>% rename(p2b = polity2, durb = durable), by = c("ccode2" = "ccode", "year")) %>%
  left_join(n %>% rename(epca = epc), by = c("ccode1" = "ccode", "year")) %>% left_join(n %>% rename(epcb = epc), by = c("ccode2" = "ccode", "year"))
stopifnot(nrow(d) == n0)
d <- d %>% filter(!(dura %in% 0) & !(durb %in% 0)) %>%
  mutate(demL_raw = pmin(p2a, p2b), devL_raw = pmin(epca, epcb),
         caprat = log(pmax(cinc1, cinc2) / pmin(cinc1, cinc2)), chist = 2^(-py / 5))
ok <- complete.cases(d[, c("fmid", "demL_raw", "devL_raw", "caprat", "chist")]) & is.finite(d$caprat)
d$demL <- d$demL_raw - median(d$demL_raw[ok]); d$devL <- d$devL_raw - median(d$devL_raw[ok]); d$dem_x_dev <- d$demL * d$devL
rhs0 <- c("demL", "devL", "dem_x_dev", "caprat", "chist")
vars <- c("fmid", rhs0, "dyadid")
d <- d[is.finite(d$caprat) | is.na(d$caprat), ]
d0 <- d[complete.cases(d[, vars]), ]
cat("complete-case N:", nrow(d0), "(published 10,503)  fatal onsets", sum(d0$fmid), "(published 353)\n")
TERM <- "dem_x_dev"
fitter <- function(r, data) fit_bin("fmid", r, data, "logit", "iid", terms = TERM)
full <- fitter(rhs0, d0)
cat("\nREBUILT vs published: dem_x_dev b", full[[TERM]][["b"]], "(pub -0.527) se", full[[TERM]][["se"]], "(pub 0.109) N", attr(full, "n"), "\n")
ds <- common_sample(d, vars, "ccode1", "ccode2", "year")
M <- run_six(rhs0, ds, fitter)
if (Sys.getenv("SWEEP_SENS") == "1") {
  source(file.path(CODE, "hist/sens_helper.R"))
  sens_split(id = "P9_067", dat = ds, fit = function(dat, add) fitter(c(rhs0, add), dat)[[TERM]])
}
saveRDS(list(M = M, full = full), file.path(W, "M_P9_067.rds"))
meta <- list(pid = "P9_067", cite = "Mousseau 2005", journal = "CMPS", year = 2005, doi = "10.1080/07388940590915327", family = "CAPITALISM")
nt <- sprintf(paste0("REBUILT (no archive): peacesciencer land-contiguous (conttype 1) dyads 1920-2000; fatal GML MID onset (paper Maoz DYMID); ",
  "Polity IV 2014 polity2 min and log NMC energy per capita min, each centred at the sample median; DURABLE=0 rows dropped; ",
  "caprat = log(max/min CINC); chist = 2^(-py/5), py = years since any MID onset. Full-sample rebuilt b=%.6f se=%.6f N=%d (events %d) vs ",
  "published -0.527 (0.109) N 10,503 (353). Common 1822-2011 sample N=%d. Interaction: constituents kept in all steps."),
  full[[TERM]][["b"]], full[[TERM]][["se"]], attr(full, "n"), sum(d0$fmid), nrow(ds))
write_p8(meta, TERM, "abstract: democracy a force for peace only in dyads above the median income; Table 1 Model 2 product term",
  "Table 1 Model 2 (market civilization; fatal MID onset, land-contiguous dyads 1920-2000)", "logit (ML)",
  "feglm logit, iid vcov, all six steps", -0.527, 0.109, 10503, TRUE, TRUE, M, "Y (conflict history 2^(-y/5))", nt)
dir.create(file.path(OUT, "results_p9"), recursive = TRUE, showWarnings = FALSE)
src <- file.path(OUT, "results_p8", "P9_067__dem_x_dev.csv"); stopifnot(file.exists(src), file.rename(src, file.path(OUT, "results_p9", basename(src))))
sink()
