# P8_022 Bueno de Mesquita, Morrow, Siverson & Smith 2004 WP "Testing Novel Implications from the Selectorate Theory of War" (10.1353/wp.2004.0017)
# Data: authors' dempeace.zip (Wayback 2016-06-04 capture of A. Smith's NYU page): ratio_data.dta (state-MID-year, built by
#   selection_effect_in_who_fights_may2004.do) + disputes_temp.dta (EUGene directed dispute dyads, for opponents);
#   Logic of Political Survival Dataverse (doi:10.7910/DVN/1LVHEU): Bdm2s2_nation_year_data_may2002 (W, irst) and NationYearCOWcapabilities (tpop).
# Headline (a priori): Table 3 Model 1, logit, DV = war involvement (hostlev==5) given dispute, unit = state-MID-year, plain logit SE, N = 6,640.
#   Key term: Proportion of Capabilities x Coalition Size (Wratioany) b = 1.486 (0.417). do-file line: logit WAR W ratioany Wratioany lnirstpc Wlnirstpc lpop Wlpop
# History (monadic state-in-dispute adaptation): mD = log1p(sum over this MID's opponents of dyadic MID onsets t-1..t-5);
#   mDM = + log1p(own MID onsets t-1..t-5); m26 = t-2..t-6 versions of both. bC = own mean lagged CINC (t-1..t-5); bD = share of t-1..t-5 with own polity2>=6.
CODE <- Sys.getenv("SWEEP_CODE"); stopifnot(nzchar(CODE))
source(file.path(CODE, "work_p8/p8_common.R"))
suppressMessages(library(haven))
WD <- file.path(ROOT, "work_p8/P8_022"); W <- file.path(OUT, "work_p8/P8_022"); dir.create(W, recursive = TRUE, showWarnings = FALSE); sink(file.path(W, "log_P8_022.txt"), split = TRUE)

# ---- rebuild the authors' estimation file ----
r  <- zap_labels(read_dta(file.path(WD, "dp/ratio_data.dta")))
cp <- zap_labels(read_dta(file.path(WD, "f_654924")))[, c("ccode", "year", "tpop")]
ny <- zap_labels(read_dta(file.path(WD, "f_654913"), col_select = c(ccode, year, W, irst)))
stopifnot(!anyDuplicated(cp[, c("ccode", "year")]))
cat("nation-year dup keys:", sum(duplicated(ny[, c("ccode", "year")])), "\n")
ny <- ny[!duplicated(ny[, c("ccode", "year")]), ]
d <- r %>% left_join(cp, by = c("ccode", "year")) %>% left_join(ny, by = c("ccode", "year"))
stopifnot(nrow(d) == nrow(r))
d <- d %>% mutate(WAR = as.numeric(!is.na(hostlev) & hostlev == 5),        # Stata: gen WAR=(hostlev==5), missing -> 0
                  lnirstpc = log(1 + irst / tpop), lpop = log(tpop + 1),
                  Wratioany = W * ratioany, Wlnirstpc = W * lnirstpc, Wlpop = W * lpop) %>% as.data.frame()
rhs <- c("W", "ratioany", "Wratioany", "lnirstpc", "Wlnirstpc", "lpop", "Wlpop")
fitter <- function(rr, data) fit_bin("WAR", rr, data, "logit", "iid", "Wratioany")
d0 <- d[complete.cases(d[, c("WAR", rhs)]), ]
cat("\n==== authors' full sample N:", nrow(d0), "(published 6,640)\n")
full <- fitter(rhs, d0)

# ---- history for a state in a dispute ----
H <- readRDS(Sys.getenv("SWEEP_HIST"))
dt <- zap_labels(read_dta(file.path(WD, "dp/disputes_temp.dta")))
opp <- bind_rows(dt %>% transmute(ccode = ccode1, opp = ccode2, midnum, year),
                 dt %>% transmute(ccode = ccode2, opp = ccode1, midnum, year)) %>% distinct()
dy <- opp %>% mutate(c1 = pmin(ccode, opp), c2 = pmax(ccode, opp)) %>%
  left_join(H$dyh, by = c("c1", "c2", "year")) %>% mutate(d15 = coalesce(d15, 0), d26 = coalesce(d26, 0)) %>%
  group_by(ccode, midnum, year) %>% summarise(od15 = sum(d15), od26 = sum(d26), n_opp = n(), .groups = "drop")
d <- d %>% left_join(dy, by = c("ccode", "midnum", "year")) %>%
  left_join(H$sth, by = c("ccode", "year")) %>% mutate(s15 = coalesce(s15, 0), s26 = coalesce(s26, 0))
cat("state-MID-years without an opponent match:", sum(is.na(d$n_opp)), "\n")
d <- d %>% mutate(hist_dyad_w15 = log1p(od15), hist_dyad_w26 = log1p(od26),
                  hist_mono_w15 = log1p(s15), hist_mono_w26 = log1p(s26))
bm <- lapply(1:5, function(L) data.frame(ccode = d$ccode, yl = d$year - L) %>%
  left_join(H$nmc, by = c("ccode", "yl" = "year")) %>% left_join(H$pol, by = c("ccode", "yl" = "year")))
d$bm_cinc_w15  <- rowMeans(sapply(bm, function(x) x$cinc), na.rm = TRUE)
d$bm_democ_w15 <- rowMeans(sapply(bm, function(x) as.numeric(x$polity2 >= 6)), na.rm = TRUE)
ds <- d[d$year >= 1822 & d$year <= 2011, ]
ds <- ds[complete.cases(ds[, c("WAR", rhs, HV)]) & is.finite(ds$bm_cinc_w15) & is.finite(ds$bm_democ_w15), ]
cat("complete-case N:", nrow(d0), " common N:", nrow(ds), "\n"); stopifnot(nrow(ds) > 0.8 * nrow(d0), sd(ds$hist_dyad_w15) > 0)
cat("share with any dyadic history w15:", round(mean(ds$hist_dyad_w15 > 0), 3), " own history w15:", round(mean(ds$hist_mono_w15 > 0), 3), "\n")
M <- run_six(rhs, ds, fitter)

meta <- list(pid = "P8_022", short_cite = "Bueno de Mesquita et al. 2004", journal = "WP", year = 2004, doi = "10.1353/wp.2004.0017",
  family = "DEMOCRACY", headline_table_model = "Table 3, Model 1 (logit, war involvement given dispute)",
  article_estimator = "logit, conventional SE, state-MID-year (war involvement | dispute), 1816-1992",
  est_used_all_steps = "logit iid SE, all six steps", author_has_history_control = "N")
write_p8(meta, "Wratioany", "Table 3 M1: Proportion of Capabilities x Coalition Size = 1.486 (0.417); H1 selectivity of large-W regimes",
  1.486, 0.417, 6640, TRUE, FALSE, M,
  paste0("Authors' dempeace.zip (Wayback) + LOPS Dataverse W/irst and COW tpop. MONADIC state-in-dispute unit: dyadic history = log1p(sum of ",
         "onsets with this MID's opponents), mono = own onsets; bC/bD monadic analogues (own lagged CINC; own democracy share)."), full = full)
sink()
