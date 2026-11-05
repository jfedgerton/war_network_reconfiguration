# P8_015 Isgandarov 2026 II "Preferential trade agreements and interstate conflict" (10.1080/03050629.2026.2715473)
# Data + code: Harvard Dataverse doi:10.7910/DVN/Q9SX4V (prd_data.csv, data_analysis.R). Data prep below is copied from data_analysis.R.
# HEADLINE: Table 3, Model 6 (full model: depweaklag + intglag + nsclag + jblag + controls + gmlmidspell cubic), logit,
#   vcovCL HC1 clustered by dyadno; politically relevant dyads 1960-2010 (undirected).
# KEY TERMS (a priori, abstract): "both institutional structures of security cooperation and joint bodies significantly reduce
#   the likelihood of conflict onset" -> nsclag (security-cooperation depth, t-1) and jblag (joint-body strength, t-1).
CODE <- Sys.getenv("SWEEP_CODE"); stopifnot(nzchar(CODE))
source(file.path(CODE, "work_p8/common_b3.R"))
suppressPackageStartupMessages({library(tidyr); library(countrycode)})
WD <- file.path(ROOT, "work_p8/P8_015"); W <- file.path(OUT, "work_p8/P8_015"); dir.create(W, recursive = TRUE, showWarnings = FALSE); sink(file.path(W, "log_P8_015.txt"), split = TRUE)
df <- read.csv(file.path(WD, "prd_data.csv"))
df$gdp1[which(df$gdp1 == 0)] <- 1; df$gdp2[which(df$gdp2 == 0)] <- 1
df <- df %>% mutate(depend1 = dyadtrade / gdp1, depend2 = dyadtrade / gdp2) %>%
  mutate(depweak = log(ifelse(depend1 > depend2, depend2, depend1) + 1), distance = log(distance + 1),
         jointdemo = ifelse(polity21 >= 6 & polity22 >= 6, 1, 0),
         capratio = log(ifelse(cinc1 > cinc2, cinc1 / cinc2, cinc2 / cinc1)),
         alliance = ifelse(cow_defense == 1 | cow_entente == 1 | cow_neutral == 1 | cow_nonagg == 1, 1, 0),
         coldwar = ifelse(year < 1992, 1, 0),
         nsc_depth = case_when(nsc01 == 0.0 ~ 0, nsc01 == 0.33 ~ 1, nsc01 == 0.66 ~ 2, nsc01 == 1.0 ~ 3),
         dyadno = as.factor(ifelse(ccode1 < ccode2, ccode1 * 1000 + ccode2, ccode2 * 1000 + ccode1)))
eu <- tribble(~country, ~joining_year, "Belgium", 1958, "France", 1958, "Germany", 1958, "Italy", 1958, "Luxembourg", 1958,
  "Netherlands", 1958, "Denmark", 1973, "Ireland", 1973, "United Kingdom", 1973, "Greece", 1981, "Portugal", 1986, "Spain", 1986,
  "Finland", 1995, "Austria", 1995, "Sweden", 1995, "Czech Republic", 2004, "Estonia", 2004, "Hungary", 2004, "Lithuania", 2004,
  "Latvia", 2004, "Malta", 2004, "Poland", 2004, "Slovakia", 2004, "Slovenia", 2004, "Cyprus", 2004, "Bulgaria", 2007, "Romania", 2007)
eu <- eu %>% mutate(cowcode = countrycode(country, "country.name", "cown"), year = purrr::map2(joining_year, 2010, ~seq(.x, .y))) %>%
  unnest(year) %>% filter(year >= 1960) %>% select(cowcode, year)
df <- df %>% left_join(eu %>% mutate(eu_member1 = 1), by = c("ccode1" = "cowcode", "year")) %>%
  left_join(eu %>% mutate(eu_member2 = 1), by = c("ccode2" = "cowcode", "year")) %>%
  mutate(eu_member1 = replace_na(eu_member1, 0), eu_member2 = replace_na(eu_member2, 0), eu_dyad = ifelse(eu_member1 == 1 & eu_member2 == 1, 1, 0))
df_lag <- df %>% group_by(dyadno) %>%
  mutate(jblag = dplyr::lag(jb_strength), nsclag = dplyr::lag(nsc_depth), intglag = dplyr::lag(intg_depth), depweaklag = dplyr::lag(depweak)) %>%
  as.data.frame() %>% drop_na(jblag, nsclag, intglag)
df_lag$sp2 <- df_lag$gmlmidspell^2; df_lag$sp3 <- df_lag$gmlmidspell^3
rhs0 <- c("depweaklag", "intglag", "nsclag", "jblag", "distance", "capratio", "jointdemo", "alliance", "dyadigos", "eu_dyad", "coldwar",
          "gmlmidspell", "sp2", "sp3")
KEY <- c("nsclag", "jblag")
fitter <- function(rhs, data) fit_bin("gmlmidonset", rhs, data, "logit", "~dyadno", terms = KEY)
vars <- c("gmlmidonset", rhs0, "dyadno")
full <- df_lag[complete.cases(df_lag[, vars]), ]
r_full <- fitter(rhs0, full); n_orig <- attr(r_full, "n")
ds <- common_sample(df_lag, vars, "ccode1", "ccode2", "year")
M <- run_six(rhs0, ds, fitter)
# ---- sensitivity hook: runs only under run_sensitivity.R (SWEEP_SENS=1), never in the main sweep ----
if (Sys.getenv("SWEEP_SENS") == "1") {
  source(file.path(CODE, "hist/sens_helper.R"))
  sens_run(id = "P8_015", dat = ds, y = "gmlmidonset", tie = "jblag", ctrl = setdiff(rhs0, "jblag"), fe = NULL,
           fit = function(dat, add) fitter(c(rhs0, add), dat)[["jblag"]])
}
saveRDS(list(M = M, r_full = r_full, n_orig = n_orig, n_used = nrow(ds)), file.path(W, "M_P8_015.rds"))
# ---- port: result-CSV write step for this study, copied verbatim from work_p8/write_b3.R (paths: saved model list read from OUT) ----
NOPDF <- "Published PDF not accessible (paywall; the device that stages sweep PDFs was disconnected), and the archive has no log with the table numbers -> UNVERIFIED."
R <- readRDS(file.path(W, "M_P8_015.rds"))
meta <- list(pid = "P8_015", cite = "Isgandarov 2026", journal = "II", year = 2026, doi = "10.1080/03050629.2026.2715473", family = "TIE")
for (k in c("nsclag", "jblag")) write_p8(meta, k,
  "abstract: 'both institutional structures of security cooperation and joint bodies significantly reduce the likelihood of conflict onset'",
  "Table 3 Model 6 (full model; data_analysis.R reg6)", "logit, vcovCL HC1 clustered by dyad, gmlmidspell cubic; PRD 1960-2010",
  "logit (fixest), SE clustered by dyad (same small-sample adj. as vcovCL HC1), all six steps", NA, NA, NA, FALSE, FALSE, R$M, "N (peace-spell cubic only)",
  sprintf(paste("%s Authors' own data + code (Dataverse Q9SX4V) rerun: full-sample N=%d, %s b=%.4f p=%.3g (author script output = Table 3 by construction).",
    "Common 1822-2011 sample N=%d."), NOPDF, R$n_orig, k, R$r_full[[k]][["b"]], R$r_full[[k]][["p"]], R$n_used), verdict = "UNVERIFIED")
sink()
