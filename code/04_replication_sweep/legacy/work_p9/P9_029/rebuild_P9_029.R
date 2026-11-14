# P9_029 Gartzke 1998 AJPS 42(1) "Kant We All Just Get Along? Opportunity, Willingness, and the Origins of the Democratic Peace"
#   (10.2307/2991745). REBUILD (protocol B): no author code or archive; the affinity index was "available from the author".
# Article (scanned PDF, printed pp. 20-23): logit of MID involvement (dichotomous), politically relevant dyads 1950-85, re-estimating
#   Oneal, Oneal, Maoz & Russett (1996) Table 1 Eqns 1 / 1a and adding UN affinity (Spearman rho of 3-category UNGA votes, lag 1).
#   Table 1 (all PR dyads) is the main table; Table 2 (contiguous dyads only) is the subsample check.
# HEADLINE (fixed from the PDF before estimating): Table 1 Model IIa, key term Affinity:
#   Democrat (dichotomous joint democracy), Growth, Allies, Contig, Capratio (raw), Interdep, Affinity (lag 1); no temporal controls.
#   Published Affinity b = -1.12, SE 0.143 (p<.001), N = 19,253. (Model Ia, Joinreg version: -1.22 (0.142), N 18,286.)
#   Why IIa: Table 1 is the full sample; within it the paper treats Ia/IIa symmetrically for the affinity claim (and the conclusion
#   stresses Democrat staying significant); fn 6 calls Joinreg nonlinear, used only to replicate prior work; Joinreg (Maoz-Russett
#   power-concentration-weighted index) cannot be rebuilt from TRIANGLE, Democrat can.
# SE: table note says "standard errors" (Stata 4/5 logit, no robust/cluster wording) -> conventional ML SEs (vcov iid), chosen a priori.
# Rebuild inputs (reused from X076): Russett-Oneal TRIANGLE PR dyads (dispute1, demaut, allies, lcaprat2, noncontg, depend);
#   Gartzke Affinity of Nations 2010 release s3un4608i (S-score, 3-category, 1964 interpolated) lagged 1 yr;
#   lower GDPpc growth (mean of t-1..t-3) from peacesciencer SDP estimates.
CODE <- Sys.getenv("SWEEP_CODE"); stopifnot(nzchar(CODE))
source(file.path(CODE, "work_p8/common_b3.R"))   # ROOT, OUT, dplyr, fixest, add_history, set.seed(123)
suppressPackageStartupMessages({library(haven); library(peacesciencer)})
W <- file.path(OUT, "work_p9/P9_029"); dir.create(W, recursive = TRUE, showWarnings = FALSE)
sink(file.path(W, "log_P9_029.txt"), split = TRUE)

# ---- 1. Russett-Oneal TRIANGLE, PR dyads 1950-85 ----
t <- read.delim(file.path(ROOT, "recA/or/TRIANGLE.tab"))
stopifnot(nrow(t) == 39996)
t <- t[t$year >= 1950 & t$year <= 1985, ]
stopifnot(nrow(t) == 23589, all(t$statea < t$stateb))
t <- t %>% mutate(c1 = pmin(statea, stateb), c2 = pmax(statea, stateb),
  democrat = as.numeric(demauta >= 6 & demautb >= 6),           # both Polity dem-aut >= 6
  caprat = exp(lcaprat2),                                        # raw ratio stronger/weaker (lcaprat2 = log ratio >= 0)
  contig = -noncontg,                                            # noncontg stored as (dummy - 1): 0 = noncontig, -1 = contig
  smldep = pmin(dependa, dependb))                               # lower bilateral trade/GDP (weak link)
stopifnot(all(t$contig %in% c(0, 1)), all(t$caprat >= 1 - 1e-8), anyDuplicated(t[, c("c1", "c2", "year")]) == 0,
          all(t$dispute1 %in% c(0, 1)), all(t$allies %in% c(0, 1)))
n_t <- nrow(t); cat("TRIANGLE 1950-85 rows:", n_t, " dispute1 mean:", round(mean(t$dispute1), 4), " democrat mean:", round(mean(t$democrat, na.rm = TRUE), 4), "\n")

# ---- 2. UN affinity, lag 1 (affinity in year y -> dyad-year y+1) ----
a <- zap_labels(read_dta(file.path(ROOT, "work_p4/X076/affinity_01242010.dta"))) %>% as.data.frame()
stopifnot("s3un4608i" %in% names(a))
# [2026-11-14] dispute1 is t+1, so affinity in the row year is the one-period lag (was year + 1 = double lag)
a <- a %>% transmute(c1 = pmin(ccodea, ccodeb), c2 = pmax(ccodea, ccodeb), year = year, aff_l1 = s3un4608i) %>% filter(c1 != c2) %>%
  group_by(c1, c2, year) %>% summarise(aff_l1 = mean(aff_l1, na.rm = TRUE), .groups = "drop") %>% mutate(aff_l1 = ifelse(is.nan(aff_l1), NA, aff_l1))
t <- t %>% left_join(a, by = c("c1", "c2", "year"))
stopifnot(nrow(t) == n_t, all(is.na(t$aff_l1) | abs(t$aff_l1) <= 1))
cat("affinity (lag 1) non-missing:", sum(!is.na(t$aff_l1)), "of", n_t, "\n")

# ---- 3. lower economic growth: mean annual % change in GDPpc over t-1..t-3 (SDP), weak link ----
sy <- create_stateyears(system = "cow", subset_years = 1940:1990) %>% add_sdp_gdp()
stopifnot("wbgdppc2011est" %in% names(sy), anyDuplicated(sy[, c("ccode", "year")]) == 0)
sy <- sy %>% arrange(ccode, year) %>% group_by(ccode) %>%
  mutate(gpc = exp(wbgdppc2011est), gr = 100 * (gpc / lag(gpc) - 1), gr3 = (lag(gr, 1) + lag(gr, 2) + lag(gr, 3)) / 3) %>%
  ungroup() %>% select(ccode, year, gr3) %>% as.data.frame()
t <- t %>% left_join(sy %>% rename(ga = gr3), by = c("statea" = "ccode", "year")) %>%
  left_join(sy %>% rename(gb = gr3), by = c("stateb" = "ccode", "year")) %>%
  mutate(lowgrowth = pmin(ga, gb), dyadid = 1000 * c1 + c2) %>% as.data.frame()
stopifnot(nrow(t) == n_t)
cat("lowgrowth non-missing:", sum(!is.na(t$lowgrowth)), "of", n_t, "\n")

# ---- 4. headline model: Table 1 Model IIa ----
rhs0 <- c("aff_l1", "democrat", "lowgrowth", "allies", "contig", "caprat", "smldep")
vars <- c("dispute1", rhs0, "dyadid")
fitter <- function(r, data) fit_bin("dispute1", r, data, "logit", "iid", terms = c("aff_l1", "democrat"))
d0 <- t[complete.cases(t[, vars]), ]
stopifnot(nrow(d0) > 15000)
full <- fitter(rhs0, d0)
cat("\nREBUILD vs published (Table 1 IIa): aff_l1 b", full$aff_l1[["b"]], "se", full$aff_l1[["se"]], "p", full$aff_l1[["p"]],
    "N", attr(full, "n"), "| published -1.12 (0.143) N 19253\n")
cat("democrat b", full$democrat[["b"]], "se", full$democrat[["se"]], "| published -1.16 (0.273)\n")
ref_cl <- fit_bin("dispute1", rhs0, d0, "logit", "~dyadid", terms = "aff_l1")   # reference only: dyad-clustered SE
cat("reference dyad-clustered SE for aff_l1:", ref_cl$aff_l1[["se"]], "\n")

# [2026-11-13 audit] TRIANGLE dispute1 in row year t is MID at t+1 (checked: P(dispute1=1 | MID onset at t+1)=.994) -> anchor on outcome year
t$yr_out <- t$year + 1
ds <- common_sample(t, vars, "statea", "stateb", "yr_out")
M <- run_six(rhs0, ds, fitter)
saveRDS(list(M = M, full = full, n_used = nrow(ds)), file.path(W, "M_P9_029.rds"))

# ---- sensitivity hook: runs only under run_sensitivity.R (SWEEP_SENS=1), never in the main sweep ----
if (Sys.getenv("SWEEP_SENS") == "1") {
  source(file.path(CODE, "hist/sens_helper.R"))
  sens_run(id = "P9_029", dat = ds, y = "dispute1", tie = "aff_l1", ctrl = setdiff(rhs0, "aff_l1"), fe = NULL,
           fit = function(dat, add) fitter(c(rhs0, add), dat)[["aff_l1"]])
}

meta <- list(pid = "P9_029", cite = "Gartzke 1998", journal = "AJPS", year = 1998, doi = "10.2307/2991745", family = "TIE")
art <- "logit, conventional SE, no temporal controls (Oneal et al. 1996 replication + affinity)"
used <- "feglm logit, iid SE, all six steps"
ahc <- "N (no temporal/history controls)"
nt <- sprintf(paste0("REBUILT (no archive): Russett-Oneal TRIANGLE PR dyads 1950-85 (dispute1, Democrat = both demaut>=6, allies, raw capability ratio, ",
  "contiguity, lower trade/GDP) + Gartzke Affinity 2010 s3un4608i lag 1 (paper: Spearman rho on 3-category votes) + lower GDPpc growth ",
  "(SDP; paper: Oneal et al. 1996 growth). Full rebuilt N %d (published 19,253): b=%.5f se=%.5f. Dyad-clustered SE %.5f for reference. ",
  "Common 1822-2011 sample N %d. Table 1 Model IIa (Model Ia Joinreg version not rebuildable)."),
  attr(full, "n"), full$aff_l1[["b"]], full$aff_l1[["se"]], ref_cl$aff_l1[["se"]], nrow(ds))
write_p8(meta, "aff_l1", "abstract: 'national preferences account for the lack of conflict between democracies'",
  "Table 1 Model IIa (all PR dyads 1950-85; Democrat + UN affinity lag 1)", art, used,
  -1.12, 0.143, 19253, TRUE, TRUE, M, ahc, nt)
# phase-9 output location (shared write_p8 writes to results_p8)
dir.create(file.path(OUT, "results_p9"), recursive = TRUE, showWarnings = FALSE)
src <- file.path(OUT, "results_p8", "P9_029__aff_l1.csv"); dst <- file.path(OUT, "results_p9", "P9_029__aff_l1.csv")
stopifnot(file.exists(src), file.rename(src, dst))
sink()
