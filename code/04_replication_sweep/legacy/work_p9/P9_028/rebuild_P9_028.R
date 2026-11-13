# P9_028 Russett, Oneal & Davis 1998 IO 52(3) "The Third Leg of the Kantian Tripod for Peace: International Organizations and
#   Militarized Disputes, 1950-85" (10.1162/002081898550626). REBUILD (protocol B; no author code or archive).
# Article (PDF p.455, Table 1, the only DV=dispute model): probit, Huber robust SE clustered on dyad (Stata cluster(dyad)), two-tailed;
#   politically relevant dyads (contiguous incl. via colonies, or containing a COW major power) 1950-85; DV = involvement in a MID in year t.
#   RHS: IGO_ij,t-1 (shared IGO memberships, lag 1), DEM_L, DEM_H (Polity III democ-autoc), DEPEND_L,t-1 (lower trade/GDP, lag 1),
#   dDEPEND_H (DEPEND_t-1 - DEPEND_t-4 per state, 1-2 yr fallback; the larger absolute change), GROWTH_L (lower 3-yr mean real GDPpc
#   growth, 1-2 yr fallback), CAPRATIO (raw stronger/weaker CINC), ALLIES (direct alliance OR both allied with the US), CONTIG.
#   No peace-year controls in Table 1. Published: IGO -0.008 (SE 0.003, p .01), N = 19,752, LL -3210.2.
# KEY (a priori): igo_l1 -- "We find that shared memberships in IGOs do make a difference".
# Rebuild inputs: Russett-Oneal TRIANGLE.tab (PR dyads 1885-1991; dispute1, demaut, depend, allies, lcaprat2, smigoabi, noncontg,
#   all stored at year t, verified: depend tracks COW trade at t, demaut tracks Polity at t) + lower growth from peacesciencer SDP GDP
#   estimates (X076 pattern). Lags of IGO and dependence built by matching the same dyad's year t-1 row in TRIANGLE.
ROOT <- Sys.getenv("SWEEP_DATA"); OUT <- Sys.getenv("SWEEP_OUT"); CODE <- Sys.getenv("SWEEP_CODE")
stopifnot(nzchar(ROOT), nzchar(OUT), nzchar(CODE))
source(file.path(CODE, "work_p8/common_b3.R"))   # dplyr, fixest, add_history, common_sample, run_six, write_p8, fit_bin; set.seed(123)
suppressPackageStartupMessages(library(peacesciencer))
set.seed(123)
W <- file.path(OUT, "work_p9/P9_028"); dir.create(W, recursive = TRUE, showWarnings = FALSE)
sink(file.path(W, "log_P9_028.txt"), split = TRUE)

# ---- 1. TRIANGLE, all years (lags need pre-1950 rows) ----
t <- read.delim(file.path(ROOT, "recA/or/TRIANGLE.tab"))
stopifnot(nrow(t) == 39996, all(t$statea < t$stateb), anyDuplicated(t[, c("statea", "stateb", "year")]) == 0,
          !anyNA(t[, c("dispute1", "smigoabi", "dependa", "dependb", "demauta", "demautb", "allies", "lcaprat2", "noncontg")]))
t <- t %>% mutate(dyadid = 1000 * statea + stateb)
# [2026-11-13 audit] dispute1 in row year t is MID at t+1 (checked vs COW dyadic MID). The article's DV is MID in year t with lagged IGO/DEPEND,
# so the DV is the same dyad's dispute1 from row t-1 (previously the t+1 outcome was paired with t-1 covariates, a two-year lag).
t <- t %>% left_join(t %>% transmute(statea, stateb, year = year + 1, dispute_t = dispute1), by = c("statea", "stateb", "year"))
stopifnot(nrow(t) == 39996)

# ---- 2. lags of IGO and dependence (same dyad, years t-1 .. t-4) ----
lagk <- function(k) t %>% transmute(statea, stateb, year = year + k, igo = smigoabi, da = dependa, db = dependb) %>%
  rename_with(~paste0(.x, "_l", k), c("igo", "da", "db"))
n0 <- nrow(t)
for (k in 1:4) { t <- t %>% left_join(lagk(k), by = c("statea", "stateb", "year")); stopifnot(nrow(t) == n0) }
dtr <- function(l1, l2, l3, l4) ifelse(!is.na(l4), l1 - l4, ifelse(!is.na(l3), l1 - l3, l1 - l2))   # t-1 minus t-4, fallback t-3, t-2
t <- t %>% mutate(igo_l1 = igo_l1, depl_l1 = pmin(da_l1, db_l1),
  dda = dtr(da_l1, da_l2, da_l3, da_l4), ddb = dtr(db_l1, db_l2, db_l3, db_l4),
  ddeph = ifelse(is.na(dda) | is.na(ddb), NA_real_, ifelse(abs(dda) >= abs(ddb), dda, ddb)),
  deml = pmin(demauta, demautb), demh = pmax(demauta, demautb), caprat = exp(lcaprat2), contig = -noncontg)
stopifnot(all(t$contig %in% c(0, 1)))

# ---- 3. ALLIES: direct (TRIANGLE) or both allied with the US (US dyads are all in the PR set) ----
us <- t %>% filter(statea == 2 | stateb == 2, allies == 1) %>% transmute(o = ifelse(statea == 2, stateb, statea), year, usally = 1L) %>% distinct()
t <- t %>% left_join(us %>% rename(usa = usally), by = c("statea" = "o", "year")) %>% left_join(us %>% rename(usb = usally), by = c("stateb" = "o", "year"))
stopifnot(nrow(t) == n0)
t <- t %>% mutate(allies98 = as.integer(allies == 1 | (coalesce(usa, 0L) == 1 & coalesce(usb, 0L) == 1 & statea != 2)))
cat("allies: direct", sum(t$allies), " direct-or-both-US", sum(t$allies98), "\n")

# ---- 4. GROWTH_L: mean annual % change in SDP GDPpc over t-3..t-1 (fallback 1-2 yrs), weak link ----
sy <- tryCatch(create_stateyears(system = "cow", subset_years = 1940:1992) %>% add_sdp_gdp(), error = function(e) { cat("add_sdp_gdp failed:", conditionMessage(e), "\n"); NULL })
stopifnot(!is.null(sy), "wbgdppc2011est" %in% names(sy))
sy <- sy %>% transmute(ccode, year, gpc = exp(wbgdppc2011est))
sy <- sy %>% left_join(sy %>% transmute(ccode, year = year + 1, gpc_p = gpc), by = c("ccode", "year")) %>% mutate(gr = 100 * (gpc / gpc_p - 1))
g3 <- sy %>% select(ccode, year)
for (k in 1:3) g3 <- g3 %>% left_join(sy %>% transmute(ccode, year = year + k, !!paste0("g", k) := gr), by = c("ccode", "year"))
g3 <- g3 %>% mutate(gr3 = ifelse(is.na(g1) & is.na(g2) & is.na(g3), NA_real_, rowMeans(cbind(g1, g2, g3), na.rm = TRUE))) %>% select(ccode, year, gr3)
stopifnot(anyDuplicated(g3[, c("ccode", "year")]) == 0)
t <- t %>% left_join(g3 %>% rename(ga = gr3), by = c("statea" = "ccode", "year")) %>% left_join(g3 %>% rename(gb = gr3), by = c("stateb" = "ccode", "year"))
stopifnot(nrow(t) == n0)
t <- t %>% mutate(growl = pmin(ga, gb)) %>% as.data.frame()

# ---- 5. sample 1950-85 ----
d <- t[t$year >= 1950 & t$year <= 1985, ]
stopifnot(nrow(d) == 23589)
rhs0 <- c("igo_l1", "deml", "demh", "depl_l1", "ddeph", "growl", "caprat", "allies98", "contig")
vars <- c("dispute_t", rhs0, "dyadid")
cat("1950-85 rows:", nrow(d), "\nmissing by var:\n"); print(sapply(d[, vars], function(x) sum(is.na(x))))
d0 <- d[complete.cases(d[, vars]), ]
cat("complete-case N:", nrow(d0), "(published 19,752)  igo_l1 mean", round(mean(d0$igo_l1), 2), "max", max(d0$igo_l1), "(paper: mean 33, max 129)\n")
cat("cor igo-deml", round(cor(d0$igo_l1, d0$deml), 3), " igo-depl", round(cor(d0$igo_l1, d0$depl_l1), 3), " deml-depl", round(cor(d0$deml, d0$depl_l1), 3), "(paper .52 .47 .32)\n")

# ---- 6. headline fit: probit, SE clustered on dyad ----
fitter <- function(r, data) fit_bin("dispute_t", r, data, "probit", "~dyadid", terms = "igo_l1")
full <- fitter(rhs0, d0)
cat("\nREBUILT vs published: igo_l1 b", full$igo_l1[["b"]], "(pub -0.008) se", full$igo_l1[["se"]], "(pub 0.003) p", full$igo_l1[["p"]],
    "(pub .01) N", attr(full, "n"), "(pub 19,752)\n")
ds <- common_sample(d, vars, "statea", "stateb", "year")
M <- run_six(rhs0, ds, fitter)
saveRDS(list(M = M, full = full), file.path(W, "M_P9_028.rds"))

# ---- sensitivity hook: runs only under run_sensitivity.R (SWEEP_SENS=1), never in the main sweep ----
if (Sys.getenv("SWEEP_SENS") == "1") {
  source(file.path(CODE, "hist/sens_helper.R"))
  sens_run(id = "P9_028", dat = ds, y = "dispute_t", tie = "igo_l1", ctrl = setdiff(rhs0, "igo_l1"), fe = NULL,
           fit = function(dat, add) fitter(c(rhs0, add), dat)[["igo_l1"]])
}

meta <- list(pid = "P9_028", cite = "Russett, Oneal & Davis 1998", journal = "IO", year = 1998, doi = "10.1162/002081898550626", family = "TIE")
art <- "probit, Huber robust SE clustered on dyad, two-tailed; no peace-year controls (Table 1)"
used <- "feglm probit, vcov clustered dyadid, all six steps"
ahc <- "N (no temporal/history controls)"
nt <- sprintf(paste0("REBUILT (no archive): Russett-Oneal TRIANGLE PR dyads 1950-85 (their dispute1/Polity/dependence/allies/caprat/contiguity/IGO coding, ",
  "book-era revision, not the 1998 file). IGO = smigoabi lagged 1 yr (TRIANGLE count differs from 1998 count, e.g. US-PRC 20+ shared IGOs pre-1971 vs paper zero). ",
  "Lower dependence lagged 1 yr; trend in dependence built from TRIANGLE dependence (t-1 minus t-4, fallback t-3/t-2, larger abs change); ",
  "lower growth from peacesciencer SDP GDPpc (paper: PWT), 3-yr mean with 1-2 yr fallback; allies = direct or both allied with US; caprat = exp(lcaprat2). ",
  "Full-sample rebuilt b=%.6f se=%.6f p=%.4f N=%d vs published -0.008 (0.003) N 19,752. Common 1822-2011 sample N=%d."),
  full$igo_l1[["b"]], full$igo_l1[["se"]], full$igo_l1[["p"]], attr(full, "n"), nrow(ds))
write_p8(meta, "igo_l1", "\"We find that shared memberships in IGOs do make a difference\"",
  "Table 1 (probit, involvement in MIDs, PR dyads 1950-85; shared IGOs t-1)", art, used,
  -0.008, 0.003, 19752, TRUE, TRUE, M, ahc, nt)
# phase-9 output location: move the write_p8 row file into results_p9 (helper is shared and writes to results_p8)
dir.create(file.path(OUT, "results_p9"), recursive = TRUE, showWarnings = FALSE)
src <- file.path(OUT, "results_p8", "P9_028__igo_l1.csv"); dst <- file.path(OUT, "results_p9", "P9_028__igo_l1.csv")
stopifnot(file.exists(src), file.rename(src, dst))
sink()
