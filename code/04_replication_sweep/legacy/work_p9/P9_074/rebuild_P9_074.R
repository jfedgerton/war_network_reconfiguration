# P9_074 Souva & Prins 2006, International Interactions 32(2) "The Liberal Peace Revisited: The Role of Democracy, Dependence, and
#   Development in Militarized Interstate Dispute Initiation, 1950-1999" (10.1080/03050620600719361). TIE. REBUILD (protocol B).
# Article (staged PDF pp. 190-195, Table 1): monadic state-years 1950-99, GEE logit with AR(1) working correlation (robust SEs);
#   DV = initiation of a fatal MID (revisionist state; MID 3.0); trade dependence = PWT openness (imports + exports over GDP, %);
#   democracy = Polity >= 6 (dummy); development = GDP per capita; democracy x development; CINC; major power; number of neighbours;
#   one-tailed p-values.
# HEADLINE (fixed before estimating): Table 1 Model A (1950-1999), Trade Dependence -0.0093 (0.0036) *** (one-tailed), N 4,586.
# Rebuild inputs: COW MIDB 4.01 participant file (work_p8/P8_049; revisionist participant entering a dispute whose participants record a
#   fatality); openness = COW national trade 3.0 total ($m current) over PWT real GDP ($m) from {peacesciencer} add_sim_gdp_pop (stand-in
#   for PWT openness), x 100; GDPpc = PWT real GDP / population; Polity IV 2014 polity2 >= 6; NMC CINC and COW majors via {peacesciencer};
#   neighbours = number of land-contiguous states (COW direct contiguity via {peacesciencer} directed dyad-years, conttype 1). GEE AR(1) logit with geepack, clustered on state.
ROOT <- Sys.getenv("SWEEP_DATA"); OUT <- Sys.getenv("SWEEP_OUT"); CODE <- Sys.getenv("SWEEP_CODE")
stopifnot(nzchar(ROOT), nzchar(OUT), nzchar(CODE))
source(file.path(CODE, "work_p8/common_b3.R")); source(file.path(CODE, "work_p9/rebuild_common_p9r.R"))
suppressPackageStartupMessages(library(geepack))
set.seed(123)
W <- file.path(OUT, "work_p9/P9_074"); dir.create(W, recursive = TRUE, showWarnings = FALSE)
sink(file.path(W, "log_P9_074.txt"), split = TRUE)
s <- suppressMessages(create_stateyears(system = "cow", subset_years = 1950:1999) %>% add_nmc() %>% add_cow_majors() %>% add_sim_gdp_pop()) %>% as.data.frame()
m <- read.csv(file.path(ROOT, "work_p8/P8_049/MIDB_401.csv"))
fd <- unique(m$DispNum3[m$Fatality > 0])
ini <- m %>% filter(RevState == 1, DispNum3 %in% fd) %>% distinct(ccode, year = StYear) %>% mutate(finit = 1L)
nb <- suppressMessages(create_dyadyears(system = "cow", directed = TRUE, subset_years = 1950:1999) %>% add_contiguity()) %>% as.data.frame() %>%
  group_by(ccode = ccode1, year) %>% summarise(nbors = sum(conttype == 1), .groups = "drop")
tr <- cow_nat_trade(); p <- polity_cow("polity2"); n0 <- nrow(s)
s <- s %>% left_join(ini, by = c("ccode", "year")) %>% left_join(nb, by = c("ccode", "year")) %>% left_join(tr, by = c("ccode", "year")) %>%
  left_join(p, by = c("ccode", "year")) %>%
  mutate(finit = coalesce(finit, 0L), nbors = coalesce(nbors, 0L), tdep = 100 * tot / pwtrgdp, dem = as.integer(polity2 >= 6),
         dev = pwtrgdp / pwtpop, dem_x_dev = dem * dev, major = as.integer(cowmaj == 1))
stopifnot(nrow(s) == n0)
s$tdep[!is.finite(s$tdep)] <- NA
rhs0 <- c("tdep", "dem", "dev", "dem_x_dev", "cinc", "major", "nbors")
vars <- c("finit", rhs0, "ccode")
cat("state-years 1950-99:", nrow(s), " fatal initiations:", sum(s$finit), " missing:\n"); print(sapply(s[, vars], function(x) sum(is.na(x))))
d0 <- s[complete.cases(s[, vars]), ]
cat("complete-case N:", nrow(d0), "(published 4,586)\n")
TERM <- "tdep"
fitter <- function(r, data) {
  data <- data[order(data$ccode, data$year), ]
  g <- geeglm(as.formula(paste("finit ~", paste(r, collapse = " + "))), id = ccode, waves = year, data = data, family = binomial("logit"), corstr = "ar1")
  ct <- summary(g)$coefficients
  cat("\n---- GEE AR1 logit finit ~", paste(r, collapse = " + "), "| N =", nrow(data), "| alpha =", round(summary(g)$corr[1, 1], 4), "\n"); print(round(ct, 6))
  out <- list(c(b = ct[TERM, "Estimate"], se = ct[TERM, "Std.err"], p = ct[TERM, "Pr(>|W|)"])); names(out) <- TERM; attr(out, "n") <- nrow(data); out
}
full <- fitter(rhs0, d0)
cat("\nREBUILT vs published: tdep b", full[[TERM]][["b"]], "(pub -0.0093) se", full[[TERM]][["se"]], "(pub 0.0036) N", attr(full, "n"), "(pub 4,586)\n")
ds <- common_sample(s, vars, "ccode", NULL, "year", monadic = TRUE)
M <- run_six(rhs0, ds, fitter, monadic = TRUE)
saveRDS(list(M = M, full = full), file.path(W, "M_P9_074.rds"))
meta <- list(pid = "P9_074", cite = "Souva & Prins 2006", journal = "II", year = 2006, doi = "10.1080/03050620600719361", family = "TIE")
nt <- sprintf(paste0("REBUILT (no archive): COW state-years 1950-99; fatal MID initiation = revisionist participant (MIDB 4.01) entering a dispute with a ",
  "recorded fatality; trade dependence = 100 x COW national trade / PWT real GDP (paper PWT openness); dem = polity2 >= 6; GDPpc = PWT GDP/pop; ",
  "CINC; COW major; land neighbours. geepack GEE logit AR(1), robust SE. Full-sample rebuilt b=%.6f se=%.6f N=%d vs published -0.0093 (0.0036) N 4,586. ",
  "Common sample N=%d (monadic history)."), full[[TERM]][["b"]], full[[TERM]][["se"]], attr(full, "n"), nrow(ds))
write_p8(meta, TERM, "abstract and Table 1: trade dependence reduces fatal MID initiation (monadic)",
  "Table 1 Model A (GEE AR1 logit, fatal MID initiation, state-years 1950-99)", "GEE logit, AR(1) working correlation, robust SE",
  "geepack geeglm logit AR(1), clustered on state, all six steps", -0.0093, 0.0036, 4586, TRUE, TRUE, M, "N (AR(1) correlation only)", nt)
dir.create(file.path(OUT, "results_p9"), recursive = TRUE, showWarnings = FALSE)
src <- file.path(OUT, "results_p8", "P9_074__tdep.csv"); stopifnot(file.exists(src), file.rename(src, file.path(OUT, "results_p9", basename(src))))
sink()
