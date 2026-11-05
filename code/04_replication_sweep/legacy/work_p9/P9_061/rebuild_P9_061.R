# P9_061 Cederman 2001 APSR 95(1) "Back to Kant: Reinterpreting the Democratic Peace as a Macrohistorical Learning Process"
#   (10.1017/S0003055401000028). REBUILD (protocol B; no author code or archive).
# Article (text extract of the staged PDF): ordinary logit (no clustering, no peace years; fn. 13 explicitly declines the BKT fix),
#   "risky" (politically relevant) dyad-years, DV = MID (COW) per dyad-year. D = both states >= 6 on Polity III DEMOC (Gurr democracy).
#   Learning curves: X'b = b0 + b1 D + b2 D*Y + b3 (1-D)*Y [+ b4 INTER + b5 COLD + b6 ALL + b7 CAP + b8 DEV], Y = calendar year.
#   Tables 2-3 drop the world wars and their two-year "shadows" (1914-1920, 1939-1947; fn. 17, Table 4 header); INTER = 1921-38,
#   COLD = 1948-92. CAP = capability balance, 1 = perfect symmetry, 0 = total asymmetry (weaker/stronger CINC). DEV = lowest dyadic
#   energy consumption per capita (COW NMC). ALL = common COW alliance.
# HEADLINE (fixed before estimating): Table 3, dynamic model (full specification: epoch dummies + controls), key term b2 = D x Y
#   ("Year (Dem.)"): -0.0269 (SE 0.0052, p .0001), N = 43,783 (risky dyads 1837-1992 without world wars). The estimator is a
#   single-equation logit with a product term, NOT a time-varying-parameter model (TVP work is Cederman & Penubarti 1999, fn. 22),
#   so the design is portable.
# Rebuild inputs: Russett-Oneal TRIANGLE.tab (PR dyads 1885-1991: dispute1, allies, lcaprat2) + Polity IV 2014 DEMOC (stand-in for
#   Polity III) + COW NMC v4 energy and population (stand-in for the Bennett-Stam 1999 EUGene NMC). Files already on disk (P8_049).
ROOT <- Sys.getenv("SWEEP_DATA"); OUT <- Sys.getenv("SWEEP_OUT"); CODE <- Sys.getenv("SWEEP_CODE")
stopifnot(nzchar(ROOT), nzchar(OUT), nzchar(CODE))
source(file.path(CODE, "work_p8/common_b3.R"))   # dplyr, fixest, add_history, common_sample, run_six, write_p8, fit_bin; set.seed(123)
set.seed(123)
W <- file.path(OUT, "work_p9/P9_061"); dir.create(W, recursive = TRUE, showWarnings = FALSE)
sink(file.path(W, "log_P9_061.txt"), split = TRUE)

# ---- 1. TRIANGLE (PR dyads 1885-1991) ----
t <- read.delim(file.path(ROOT, "recA/or/TRIANGLE.tab"))
stopifnot(nrow(t) == 39996, all(t$statea < t$stateb), anyDuplicated(t[, c("statea", "stateb", "year")]) == 0,
          !anyNA(t[, c("dispute1", "allies", "lcaprat2")]), all(t$lcaprat2 >= 0), all(t$dispute1 %in% 0:1))
t <- t %>% mutate(dyadid = 1000 * statea + stateb)
n0 <- nrow(t)

# ---- 2. Polity DEMOC (Gurr democracy scale), mapped to COW codes ----
p <- read.csv(file.path(ROOT, "work_p8/P8_049/p4v2014.csv"))[, c("ccode", "year", "democ")]
p <- p %>% filter(year >= 1880, year <= 1992) %>%
  mutate(democ = ifelse(democ %in% c(-66, -77, -88), NA_real_, democ),
         ccode = case_when(ccode == 364 ~ 365L, ccode == 769 ~ 770L, ccode == 818 ~ 816L, ccode == 342 & year < 1921 ~ 345L, TRUE ~ as.integer(ccode))) %>%
  filter(!(ccode == 816 & year == 1976 & is.na(democ))) %>%
  group_by(ccode, year) %>% summarise(democ = if (all(is.na(democ))) NA_real_ else max(democ, na.rm = TRUE), n = n(), .groups = "drop")
cat("Polity ccode-years with >1 source row after mapping:", sum(p$n > 1), "\n"); print(as.data.frame(p[p$n > 1, ]))
stopifnot(sum(p$n > 1) <= 5); p$n <- NULL
t <- t %>% left_join(p %>% rename(democa = democ), by = c("statea" = "ccode", "year")) %>%
  left_join(p %>% rename(democb = democ), by = c("stateb" = "ccode", "year"))
stopifnot(nrow(t) == n0)
t <- t %>% mutate(D = as.integer(democa >= 6 & democb >= 6))
cat("D missing:", sum(is.na(t$D)), " joint democracies:", sum(t$D, na.rm = TRUE), "\n")
cat("cross-check vs TRIANGLE demaut (both >= 6):", round(mean(t$D == as.integer(t$demauta >= 6 & t$demautb >= 6), na.rm = TRUE), 3), "agreement\n")

# ---- 3. NMC v4: lower energy consumption per capita (DEV), raw scale (paper gives no transform) ----
n <- read.csv(file.path(ROOT, "work_p8/P8_049/nmc_v4.csv"))[, c("ccode", "year", "pec", "tpop")]
stopifnot(anyDuplicated(n[, c("ccode", "year")]) == 0)
n <- n %>% mutate(epc = ifelse(pec >= 0 & tpop > 0, pec / tpop, NA_real_)) %>% select(ccode, year, epc)   # thousand coal-ton / thousand people
t <- t %>% left_join(n %>% rename(epca = epc), by = c("statea" = "ccode", "year")) %>% left_join(n %>% rename(epcb = epc), by = c("stateb" = "ccode", "year"))
stopifnot(nrow(t) == n0)
t <- t %>% mutate(dev = pmin(epca, epcb), cap = exp(-lcaprat2), all = as.integer(allies == 1),
  dem_year = D * year, oth_year = (1 - D) * year,
  inter = as.integer(year >= 1921 & year <= 1938), cold = as.integer(year >= 1948 & year <= 1992)) %>% as.data.frame()
stopifnot(all(t$cap > 0 & t$cap <= 1))

# ---- 4. sample: drop world wars + two-year shadows ----
d <- t[!(t$year %in% c(1914:1920, 1939:1947)), ]
cat("rows without world wars:", nrow(d), "of", n0, "\n")
rhs0 <- c("D", "dem_year", "oth_year", "inter", "cold", "all", "cap", "dev")
vars <- c("dispute1", rhs0, "dyadid")
cat("missing by var:\n"); print(sapply(d[, vars], function(x) sum(is.na(x))))
d0 <- d[complete.cases(d[, vars]), ]
cat("complete-case N:", nrow(d0), "(published 43,783, 1837-1992)  MID rate", round(mean(d0$dispute1), 4), "\n")

# ---- 5. headline fit: logit, conventional SEs ----
fitter <- function(r, data) fit_bin("dispute1", r, data, "logit", "iid", terms = "dem_year")
full <- fitter(rhs0, d0)
cat("\nREBUILT vs published: dem_year b", full$dem_year[["b"]], "(pub -0.0269) se", full$dem_year[["se"]], "(pub 0.0052) p",
    full$dem_year[["p"]], "(pub .0001) N", attr(full, "n"), "(pub 43,783)\n")
ds <- common_sample(d, vars, "statea", "stateb", "year")
M <- run_six(rhs0, ds, fitter)
saveRDS(list(M = M, full = full), file.path(W, "M_P9_061.rds"))

# sensitivity hook: NOT wired. The key term is a product (D x Year) entered alongside D and (1-D) x Year (separate-slopes
# parametrization), i.e. an interaction headline -> "split" class under PORT_SPEC_P9, not a single additive column.

meta <- list(pid = "P9_061", cite = "Cederman 2001", journal = "APSR", year = 2001, doi = "10.1017/S0003055401000028", family = "DEMOCRACY")
art <- "logit, conventional SEs (no clustering, no peace years; fn. 13), Table 3 dynamic model"
used <- "feglm logit, iid vcov, all six steps"
ahc <- "N (no temporal/history controls; BKT correction explicitly declined)"
nt <- sprintf(paste0("REBUILT (no archive): TRIANGLE PR dyads 1885-1991 (paper: risky dyads 1837-1992) minus 1914-20 and 1939-47; dispute1 = MID involvement; ",
  "D = both Polity IV 2014 DEMOC >= 6 (paper Polity III; -66/-77/-88 missing); D x Year and (1-D) x Year raw calendar year; INTER 1921-38, COLD 1948-92; ",
  "ALL = TRIANGLE allies; CAP = exp(-lcaprat2) (weaker/stronger CINC); DEV = lower NMC v4 energy per capita, raw (paper: EUGene NMC, transform not stated). ",
  "Full-sample rebuilt b=%.6f se=%.6f p=%.4g N=%d vs published -0.0269 (0.0052) N 43,783. Common 1822-2011 sample N=%d."),
  full$dem_year[["b"]], full$dem_year[["se"]], full$dem_year[["p"]], attr(full, "n"), nrow(ds))
write_p8(meta, "dem_year", "\"the conflict propensities among democracies exhibit a steadily falling trend\"",
  "Table 3 dynamic model (logit, MID, risky dyads without world wars; Year x joint democracy)", art, used,
  -0.0269, 0.0052, 43783, TRUE, TRUE, M, ahc, nt)
dir.create(file.path(OUT, "results_p9"), recursive = TRUE, showWarnings = FALSE)
src <- file.path(OUT, "results_p8", "P9_061__dem_year.csv"); dst <- file.path(OUT, "results_p9", "P9_061__dem_year.csv")
stopifnot(file.exists(src), file.rename(src, dst))
sink()
