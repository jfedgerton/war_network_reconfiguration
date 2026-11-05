# P9_075 Helfstein 2012, International Finance 15(3) "Liabilities of Globalization: Sovereign Debt, International Investors, and Interstate
#   Conflict" (10.1111/j.1468-2362.2013.12004.x). CAPITALISM (foreign-held public debt). REBUILD (protocol B; no archive).
# Article (staged PDF pp. 283-285, Table 1): time-series cross-section logit, state-years 1970-2001 (World Bank debt-reporting countries);
#   DV = state initiates a MID (MID SideA = 1); key IV = external debt stocks, public and publicly guaranteed (World Bank);
#   controls GDP, exports, CINC, Polity, inflation (Model 2). Units: the text describes the debt SD as $7.3bn-$23.9bn -> US$ billions.
#   Stars are lost in the PDF text layer; 0.015 (0.008) implies p ~ .06; legend "P < 0.10, P < 0.05, P < 0.01".
# HEADLINE (fixed before estimating): Table 1 Model 2, External Public Debt Stocks 0.015 (0.008); N not reported.
# Rebuild inputs: WDI API series DT.DOD.DPPG.CD, NY.GDP.MKTP.CD, NE.EXP.GNFS.CD, FP.CPI.TOTL.ZG (work_p9/P9_075/wdi_P9_075.csv; US$ -> bn);
#   MID initiation = SideA participant entering a dispute that year (COW MIDB 4.01, work_p8/P8_049); NMC CINC via {peacesciencer};
#   Polity IV 2014 polity2. Logit with SEs clustered on state (TSCS logit; clustering assumed). ISO3 -> COW codes via {countrycode}.
ROOT <- Sys.getenv("SWEEP_DATA"); OUT <- Sys.getenv("SWEEP_OUT"); CODE <- Sys.getenv("SWEEP_CODE")
stopifnot(nzchar(ROOT), nzchar(OUT), nzchar(CODE))
source(file.path(CODE, "work_p8/common_b3.R")); source(file.path(CODE, "work_p9/rebuild_common_p9r.R"))
set.seed(123)
W <- file.path(OUT, "work_p9/P9_075"); dir.create(W, recursive = TRUE, showWarnings = FALSE)
sink(file.path(W, "log_P9_075.txt"), split = TRUE)
w <- read.csv(file.path(ROOT, "work_p9/P9_075/wdi_P9_075.csv"))
w$ccode <- suppressWarnings(countrycode::countrycode(w$iso3c, "iso3c", "cown"))
w$ccode[w$iso3c == "SRB"] <- 345; w$ccode[w$iso3c == "VNM"] <- 816; w$ccode[w$iso3c == "YEM"] <- 679
w <- w[!is.na(w$ccode), ] %>% group_by(ccode, year) %>% summarise(across(c(debt_ppg, gdp, exports, inflation), ~ if (all(is.na(.x))) NA_real_ else max(.x, na.rm = TRUE)), .groups = "drop")
s <- suppressMessages(create_stateyears(system = "cow", subset_years = 1970:2001) %>% add_nmc()) %>% as.data.frame()
m <- read.csv(file.path(ROOT, "work_p8/P8_049/MIDB_401.csv"))
ini <- m %>% filter(SideA == 1) %>% distinct(ccode, year = StYear) %>% mutate(init = 1L)
p <- polity_cow("polity2"); n0 <- nrow(s)
s <- s %>% left_join(w, by = c("ccode", "year")) %>% left_join(ini, by = c("ccode", "year")) %>% left_join(p, by = c("ccode", "year")) %>%
  mutate(init = coalesce(init, 0L), debt = debt_ppg / 1e9, gdpb = gdp / 1e9, expb = exports / 1e9)
stopifnot(nrow(s) == n0)
rhs0 <- c("debt", "gdpb", "expb", "cinc", "polity2", "inflation")
vars <- c("init", rhs0, "ccode")
cat("state-years 1970-2001:", nrow(s), " initiations:", sum(s$init), " missing:\n"); print(sapply(s[, vars], function(x) sum(is.na(x))))
d0 <- s[complete.cases(s[, vars]), ]
cat("complete-case N:", nrow(d0), "(not reported)  countries", length(unique(d0$ccode)), "(paper: 125 World Bank countries)\n")
TERM <- "debt"
fitter <- function(r, data) fit_bin("init", r, data, "logit", "~ccode", terms = TERM)
full <- fitter(rhs0, d0)
cat("\nREBUILT vs published: debt b", full[[TERM]][["b"]], "(pub 0.015) se", full[[TERM]][["se"]], "(pub 0.008) N", attr(full, "n"), "\n")
ds <- common_sample(s, vars, "ccode", NULL, "year", monadic = TRUE)
M <- run_six(rhs0, ds, fitter, monadic = TRUE)
if (Sys.getenv("SWEEP_SENS") == "1") {
  source(file.path(CODE, "hist/sens_helper.R"))
  sens_run(id = "P9_075", dat = ds, y = "init", tie = TERM, ctrl = setdiff(rhs0, TERM), fe = NULL,
           fit = function(dat, add) fitter(c(rhs0, add), dat)[[TERM]])
}
saveRDS(list(M = M, full = full), file.path(W, "M_P9_075.rds"))
meta <- list(pid = "P9_075", cite = "Helfstein 2012", journal = "IF", year = 2012, doi = "10.1111/j.1468-2362.2013.12004.x", family = "CAPITALISM")
nt <- sprintf(paste0("REBUILT (no archive): COW state-years 1970-2001 with WDI PPG external debt (US$bn); MID initiation = MIDB SideA entry year; ",
  "GDP and exports (US$bn), CINC, polity2, CPI inflation; logit, SE clustered on state. Published p ~ .06 (not significant at .05). ",
  "Full-sample rebuilt b=%.6f se=%.6f N=%d vs published 0.015 (0.008), N not reported. Common sample N=%d (monadic history)."),
  full[[TERM]][["b"]], full[[TERM]][["se"]], attr(full, "n"), nrow(ds))
write_p8(meta, TERM, "abstract: strong relationship between foreign debt ownership and conflict; Table 1 Model 2",
  "Table 1 Model 2 (logit, MID initiation, state-years 1970-2001)", "TSCS logit", "feglm logit, SE clustered on state, all six steps",
  0.015, 0.008, NA, FALSE, TRUE, M, "N (no temporal terms)", nt)
dir.create(file.path(OUT, "results_p9"), recursive = TRUE, showWarnings = FALSE)
src <- file.path(OUT, "results_p8", "P9_075__debt.csv"); stopifnot(file.exists(src), file.rename(src, file.path(OUT, "results_p9", basename(src))))
sink()
