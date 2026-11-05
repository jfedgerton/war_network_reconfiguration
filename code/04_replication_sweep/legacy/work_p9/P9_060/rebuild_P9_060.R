# P9_060 Debs & Goemans 2010 APSR 104(3) "Regime Type, the Fate of Leaders, and War" (10.1017/S0003055410000195).
#   REBUILD (protocol B; no author code or archive used).
# Article (text extract of the staged PDF, pp. 440-441): Table 5 "Nondemocratic Regime Types and Onset of War": logit, SE clustered
#   on country code; unit = leader-year (leader-year-spell) of nondemocracies post-WWII, Archigos leaders merged with Cheibub & Gandhi
#   (2004) regime types (autocracies 1946-96: monarchic / military / civilian; excluded category = military). DV = 1 if the leader
#   "launched a challenge or became a target" of a war that year. RHS: Civilian, Monarch, Civil War, GDP per capita (logged), GDP Growth,
#   Trade Openness, change in Trade Openness, Population (logged), Age, Times in Office (number of previous times in office), CINC,
#   Military Mobilization, Number of Borders, cubic time since previous onset. Published Civilian -0.351 (SE 0.196, p .074 two-tailed),
#   N 4,040, LL -144.88.
# KEY (a priori): civilian -- "civilian autocrats are indeed significantly less likely to become involved in war than military
#   autocrats". Family DEMOCRACY (regime type).
# Rebuild inputs: DD/CGV "Democracy and Dictatorship Revisited" v1 (Harvard Dataverse file 4258459; regime 3 civilian, 4 military,
#   5 royal), Archigos 4.1 via peacesciencer (paper: Archigos 2.9), COW Inter-State War v4.0 (participant entry dates, correlatesofwar.org),
#   peacesciencer: SDP GDP, NMC, COW trade, COW intra-state wars, COW direct contiguity.
ROOT <- Sys.getenv("SWEEP_DATA"); OUT <- Sys.getenv("SWEEP_OUT"); CODE <- Sys.getenv("SWEEP_CODE")
stopifnot(nzchar(ROOT), nzchar(OUT), nzchar(CODE))
source(file.path(CODE, "work_p8/common_b3.R"))   # dplyr, fixest, add_history, common_sample, run_six, write_p8, fit_bin; set.seed(123)
suppressPackageStartupMessages({library(peacesciencer); library(haven)})
set.seed(123)
W <- file.path(OUT, "work_p9/P9_060"); dir.create(W, recursive = TRUE, showWarnings = FALSE)
D9 <- file.path(ROOT, "work_p9/P9_060")
sink(file.path(W, "log_P9_060.txt"), split = TRUE)

# ---- 1. leader-year-spells (Archigos via peacesciencer, COW codes) ----
ly <- as.data.frame(create_leaderyears(standardize = "cow"))
data(archigos, package = "peacesciencer"); a <- as.data.frame(archigos)
stopifnot(anyDuplicated(a$obsid) == 0, anyDuplicated(ly[, c("obsid", "ccode", "year")]) == 0)   # GFR-1982 / YEM-1978 appear under two COW codes in 1990
a <- a %>% group_by(leadid) %>% arrange(startdate, .by_group = TRUE) %>% mutate(prevtimes = row_number() - 1L) %>% ungroup() %>%
  select(obsid, leadid, startdate, enddate, prevtimes) %>% as.data.frame()
n_ly <- nrow(ly)
ly <- ly %>% left_join(a, by = "obsid"); stopifnot(nrow(ly) == n_ly, !anyNA(ly$startdate))
ly <- ly %>% filter(year >= 1946, year <= 1996)
cat("leader-year-spells 1946-96:", nrow(ly), "\n")

# ---- 2. CGV regime type (country-year, end-of-year rule) ----
dd <- zap_labels(read_dta(file.path(D9, "ddrevisited_data_v1.dta"), col_select = c("cowcode", "year", "regime"))) %>% as.data.frame()
dd <- dd[!is.na(dd$cowcode) & !is.na(dd$regime), ]; stopifnot(anyDuplicated(dd[, c("cowcode", "year")]) == 0)
n1 <- nrow(ly)
ly <- ly %>% left_join(dd, by = c("ccode" = "cowcode", "year")); stopifnot(nrow(ly) == n1)
cat("leader-years without CGV regime:", sum(is.na(ly$regime)), "\n")
ly <- ly %>% filter(regime %in% 3:5) %>% mutate(civilian = as.integer(regime == 3), monarch = as.integer(regime == 5))
cat("nondemocratic leader-years:", nrow(ly), " civilian", sum(ly$civilian), " military", sum(ly$regime == 4), " monarch", sum(ly$monarch), "\n")

# ---- 3. war onset: COW inter-state war entry dates assigned to the leader in office on that date ----
w <- data.table::fread(file.path(D9, "Inter-StateWarData_v4.0.csv"), data.table = FALSE)
stopifnot(all(c("WarNum", "ccode", "StartYear1", "StartMonth1", "StartDay1") %in% names(w)), nrow(w) > 300)
w <- w %>% mutate(m = ifelse(StartMonth1 < 1, 1L, StartMonth1), dday = ifelse(StartDay1 < 1, 1L, StartDay1),
  wdate = as.Date(sprintf("%d-%02d-%02d", StartYear1, m, dday))) %>% select(WarNum, ccode, wyear = StartYear1, wdate) %>% distinct()
stopifnot(!anyNA(w$wdate))
hit <- ly %>% select(obsid, ccode, year, startdate, enddate) %>% inner_join(w, by = c("ccode", "year" = "wyear"), relationship = "many-to-many") %>%
  filter(wdate >= startdate, wdate <= enddate) %>% distinct(obsid, ccode, year) %>% mutate(war_onset = 1L)
ly <- ly %>% left_join(hit, by = c("obsid", "ccode", "year")) %>% mutate(war_onset = coalesce(war_onset, 0L)); stopifnot(nrow(ly) == sum(ly$regime %in% 3:5))
we <- w %>% filter(wyear >= 1946, wyear <= 1996) %>% distinct(ccode, wyear)
cat("war onsets on nondemocratic leader-years:", sum(ly$war_onset), " (state war entries 1946-96, all regimes:", nrow(we), ")\n")

# ---- 4. state-year controls ----
sy <- as.data.frame(create_stateyears(system = "cow", subset_years = 1816:1997))
sy0 <- nrow(sy)
sy <- sy %>% add_sdp_gdp() %>% add_nmc() %>% as.data.frame(); stopifnot(nrow(sy) == sy0)
data(cow_trade_sy, package = "peacesciencer")
tr <- as.data.frame(cow_trade_sy) %>% mutate(trade = ifelse(imports < 0 | exports < 0, NA_real_, imports + exports)) %>% select(ccode, year, trade)
stopifnot(anyDuplicated(tr[, c("ccode", "year")]) == 0)
sy <- sy %>% left_join(tr, by = c("ccode", "year")); stopifnot(nrow(sy) == sy0)
# civil war: state is a party (side a or b) to an ongoing COW intra-state war
data(cow_war_intra, package = "peacesciencer")
cw <- as.data.frame(cow_war_intra) %>% filter(cowintraongoing == 1) %>% { rbind(data.frame(ccode = .$ccodea, year = .$year), data.frame(ccode = .$ccodeb, year = .$year)) } %>%
  filter(ccode > 0) %>% distinct() %>% mutate(civilwar = 1L)
sy <- sy %>% left_join(cw, by = c("ccode", "year")) %>% mutate(civilwar = coalesce(civilwar, 0L)); stopifnot(nrow(sy) == sy0)
# number of borders: distinct directly land-contiguous (type 1) states in force during the year
data(cow_contdir, package = "peacesciencer")
cd <- as.data.frame(cow_contdir) %>% filter(conttype == 1) %>% mutate(y0 = as.integer(format(stdate, "%Y")), y1 = as.integer(format(enddate, "%Y")))
nb <- sy %>% select(ccode, year) %>% inner_join(cd, by = c("ccode" = "ccode1"), relationship = "many-to-many") %>% filter(year >= y0, year <= y1) %>%
  distinct(ccode, year, ccode2) %>% count(ccode, year, name = "borders")
sy <- sy %>% left_join(nb, by = c("ccode", "year")) %>% mutate(borders = coalesce(borders, 0L)); stopifnot(nrow(sy) == sy0)
# time since previous war onset (state level; clock starts at 1816 or system entry), scaled by 10
wy <- w %>% distinct(ccode, wyear)
sy <- sy %>% arrange(ccode, year) %>% group_by(ccode) %>%
  mutate(ons = as.integer(paste(ccode, year) %in% paste(wy$ccode, wy$wyear)),
         last = cummax(ifelse(ons == 1, year, first(year))),
         last_prev = lag(last, default = first(year)), tsince = (year - last_prev) / 10) %>% ungroup() %>% as.data.frame()
sy <- sy %>% arrange(ccode, year) %>% group_by(ccode) %>%
  mutate(gdp_l = lag(wbgdp2011est), open = trade / (exp(wbgdp2011est) / 1e6), open_l = lag(open), yr_l = lag(year)) %>% ungroup() %>%
  mutate(lgdppc = wbgdppc2011est, growth = ifelse(yr_l == year - 1, exp(wbgdp2011est - gdp_l) - 1, NA_real_),
         dopen = ifelse(yr_l == year - 1, open - open_l, NA_real_), lpop = ifelse(tpop > 0, log(tpop), NA_real_),
         milmob = ifelse(tpop > 0 & milper >= 0, milper / tpop, NA_real_), cinc = ifelse(cinc >= 0, cinc, NA_real_),
         tsince2 = tsince^2, tsince3 = tsince^3) %>%
  select(ccode, year, civilwar, lgdppc, growth, open, dopen, lpop, cinc, milmob, borders, tsince, tsince2, tsince3) %>% as.data.frame()
stopifnot(anyDuplicated(sy[, c("ccode", "year")]) == 0)
n2 <- nrow(ly)
d <- ly %>% left_join(sy, by = c("ccode", "year")) %>% as.data.frame(); stopifnot(nrow(d) == n2)
d$age <- d$leaderage; d$prevtimes <- as.numeric(d$prevtimes)

# ---- 5. headline fit: logit, SE clustered on country ----
rhs0 <- c("civilian", "monarch", "civilwar", "lgdppc", "growth", "open", "dopen", "lpop", "age", "prevtimes", "cinc", "milmob", "borders",
          "tsince", "tsince2", "tsince3")
vars <- c("war_onset", rhs0, "ccode")
cat("missing by var:\n"); print(sapply(d[, vars], function(x) sum(is.na(x))))
d0 <- d[complete.cases(d[, vars]), ]
cat("complete-case N:", nrow(d0), "(published 4,040)  onsets", sum(d0$war_onset), "\n")
fitter <- function(r, data) fit_bin("war_onset", r, data, "logit", "~ccode", terms = "civilian")
full <- fitter(rhs0, d0)
cat("\nREBUILT vs published: civilian b", full$civilian[["b"]], "(pub -0.351) se", full$civilian[["se"]], "(pub 0.196) p",
    full$civilian[["p"]], "(pub .074) N", attr(full, "n"), "(pub 4,040)\n")
ds <- common_sample(d, vars, "ccode", NULL, "year", monadic = TRUE)
M <- run_six(rhs0, ds, fitter, monadic = TRUE)
saveRDS(list(M = M, full = full), file.path(W, "M_P9_060.rds"))

# sensitivity hook: NOT wired. Monadic (leader-year) design: sens_run() requires a dyadic history column (hist_dyad_w15), which the
# monadic convention does not produce (same as P9_024).

meta <- list(pid = "P9_060", cite = "Debs & Goemans 2010", journal = "APSR", year = 2010, doi = "10.1017/S0003055410000195", family = "DEMOCRACY")
art <- "logit, SE clustered on country code; leader-years of nondemocracies 1946-96; cubic time since previous onset (Table 5)"
used <- "feglm logit, vcov clustered ccode, all six steps (monadic history convention)"
ahc <- "Y (cubic time since previous war onset only)"
nt <- sprintf(paste0("REBUILT (no archive): Archigos 4.1 leader-year-spells (peacesciencer; paper v2.9) x DD/CGV revisited regime (3 civ/4 mil/5 royal; paper Cheibub-Gandhi 2004) ",
  "1946-96; war onset = COW inter-state war v4.0 participant entry date (StartDate1) within the leader's spell (paper: Chiozza-Goemans challenger/target coding); ",
  "civil war = party to ongoing COW intra-state war; GDPpc (log), GDP growth from SDP; openness = COW trade / SDP GDP (price bases differ) and its 1-yr change; ",
  "log NMC population, CINC, milper/tpop; land borders from COW contdir; Age = current leader age; Times in office = prior Archigos spells; ",
  "time since previous onset /10 and its square and cube (scaling not stated). Full-sample rebuilt b=%.6f se=%.6f p=%.4g N=%d vs published -0.351 (0.196, p .074) N 4,040. ",
  "Common 1822-2011 sample N=%d."), full$civilian[["b"]], full$civilian[["se"]], full$civilian[["p"]], attr(full, "n"), nrow(ds))
write_p8(meta, "civilian", "\"civilian autocrats are indeed significantly less likely to become involved in war than military autocrats\"",
  "Table 5 (logit, war onset, nondemocratic leader-years 1946-96; civilian vs military)", art, used,
  -0.351, 0.196, 4040, FALSE, TRUE, M, ahc, nt)
dir.create(file.path(OUT, "results_p9"), recursive = TRUE, showWarnings = FALSE)
src <- file.path(OUT, "results_p8", "P9_060__civilian.csv"); dst <- file.path(OUT, "results_p9", "P9_060__civilian.csv")
stopifnot(file.exists(src), file.rename(src, dst))
sink()
