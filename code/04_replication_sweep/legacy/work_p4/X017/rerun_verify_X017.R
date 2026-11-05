# VERIFY X017 Gartzke & Weisiger 2013 ISQ, Table 2 Model 2.4 (logit, dyad-clustered SE, 1816-2000, N 634,684):
#   Dyadic x Systemic = -0.513 (SE .094). Copy of rebuild_X017.R with three changes, each justified by the paper text:
#   (1) DV = Maoz dyadic MID onset (p.178: "dummy for MID 'onset' as defined by Maoz (1999)") -> Maoz Dyadic MID v4.02
#       (COW site; Maoz 1999 DYMID1.x not posted), undirected, onset = first year of a dyadic MID (durindx==1). Was GML.
#   (2) Temporal dependence = BKT peace years + 3 splines (p.179: "four spline variables using ... Tucker (1999)") -> btscs rule
#       (quartile knots, as in earlier reruns). Was py + py^2 + py^3.
#   (3) Contiguity = 6-category ordinal (p.179: "ordinal variable for six categories of decreasing physical proximity") ->
#       peacesciencer COW direct contiguity conttype 1..5, 6 = none. Was the 2014 file's dummy. (Colonial contiguity not available.)
# RHS otherwise = authors' variables from the GW 2014 ISQ archive (demloi, dydiff, propdem, logdist, allydumy, capratio, onemajor).
# IVs lagged one year (p.178): DV year = year + 1. KEY: dxs = dydiff x propdem. CONTROL: demloi.
ROOT <- Sys.getenv("SWEEP_DATA"); OUT <- Sys.getenv("SWEEP_OUT"); CODE <- Sys.getenv("SWEEP_CODE")
stopifnot(nzchar(ROOT), nzchar(OUT), nzchar(CODE))
source(file.path(CODE, "vb2/vcommon.R"))
suppressPackageStartupMessages({library(readstata13); library(peacesciencer)})
WO <- file.path(OUT, "work_p4/X017"); dir.create(WO, recursive = TRUE, showWarnings = FALSE)  # port: outputs
W <- file.path(ROOT, "work_p4/X017"); sink(file.path(WO, "log_verify_X017.txt"), split = TRUE)

g <- read.dta13(file.path(W, "gw14t3.dta"), convert.factors = FALSE,
                select.cols = c("dyadid","year","demloi","dydiff","propdem","logdist","allydumy","capratio","onemajor"))
g$c1 <- g$dyadid %/% 1000; g$c2 <- g$dyadid %% 1000; g$dvyear <- g$year + 1
stopifnot(mean(g$c1 < g$c2) > .999, anyDuplicated(g[, c("dyadid","year")]) == 0)

# Maoz dyadic MID 4.02 onsets (undirected: keep statea < stateb records)
y <- read.csv(file.path(ROOT, "work_p4/X090/dymid/dyadic_mid_4.02.csv"))[, c("statea","stateb","year","durindx")]
on <- y %>% filter(statea < stateb, durindx == 1) %>% distinct(ccode1 = statea, ccode2 = stateb, year) %>% mutate(onset = 1)
cat("DyMID4.02 undirected onset dyad-years 1817-2001:", sum(on$year >= 1817 & on$year <= 2001), "\n")

# system panel (DV years) + contiguity + peace years from 1816
nd <- create_dyadyears(directed = FALSE, subset_years = 1816:2001) %>% add_contiguity() %>%
  select(ccode1, ccode2, year, conttype) %>% as.data.frame()
stopifnot(all(nd$conttype %in% 0:5))
nd <- nd %>% left_join(on, by = c("ccode1","ccode2","year")) %>% mutate(onset = coalesce(onset, 0)) %>% arrange(ccode1, ccode2, year)
nd$py <- btscs_py(nd$ccode1 * 1000 + nd$ccode2, nd$year, nd$onset)
nd$cont6 <- ifelse(nd$conttype == 0, 6, nd$conttype)
d <- g %>% inner_join(nd %>% select(ccode1, ccode2, dvyear = year, onset, py, cont6), by = c("c1" = "ccode1", "c2" = "ccode2", "dvyear"))
cat("matched:", nrow(d), "of", nrow(g), "\n"); stopifnot(nrow(d) > .9 * nrow(g)); rm(nd, g, y); gc()
d <- d %>% filter(year >= 1816, year <= 2000) %>% mutate(dxs = dydiff * propdem)
d <- btscs_splines(d, "py")
print(table(d$onset)); print(table(d$cont6))

f0 <- onset ~ demloi + dydiff + propdem + dxs + logdist + cont6 + allydumy + capratio + onemajor + py + spl1 + spl2 + spl3
terms <- c("dxs", "demloi"); FF <- ff_logit_cl("dyadid")
dir.create(file.path(WO, "verify_cache"), recursive = TRUE, showWarnings = FALSE)
rep <- vfit(f0, d[complete.cases(d[, c(all.vars(f0), "dyadid")]), ], terms, FF, file.path(WO, "verify_cache/v_repro.rds"))
cat("repro N:", rep$n, "(published 634,684)\n")
ds <- vcommon_sample(d, f0, "dyadid", "c1", "c2", "dvyear"); rm(d); gc()
M <- vrun_six(f0, ds, terms, FF, WO)
# ---- sensitivity hook (split only): runs only under run_sensitivity.R (SWEEP_SENS=1), never in the main sweep ----
if (Sys.getenv("SWEEP_SENS") == "1") {
  source(file.path(CODE, "hist/sens_helper.R"))
  # calls FF directly (what vfit runs): vfit would read the m0/mD/mDM cache files
  sens_split(id = "X017", dat = ds,
             fit = function(dat, add) { ct <- FF(if (length(add)) update(f0, as.formula(paste(". ~ . +", paste(add, collapse = " + ")))) else f0, dat)
               c(b = unname(ct["dxs", 1]), se = unname(ct["dxs", 2]), p = unname(ct["dxs", 3])) })
}
vsave(M, terms, WO, rep)
# port: one row per term from verify_six.csv (long) for ID-keyed extraction; no estimation (same arithmetic as vb2/assemble.py from_six)
V6 <- read.csv(file.path(WO, "verify_six.csv")); g6 <- function(m, t, k) { x <- V6[V6$model == m & V6$term == t, k]; stopifnot(length(x) == 1); x }
w6 <- do.call(rbind, lapply(terms, function(t) { b0 <- g6("m0", t, "b")
  data.frame(term = t, b0 = b0, p0 = g6("m0", t, "p"), b_dm = g6("mDM", t, "b"), p_dm = g6("mDM", t, "p"), b_26 = g6("m26", t, "b"), p_26 = g6("m26", t, "p"),
             pct_bC = 100 * (g6("bC", t, "b") - b0) / abs(b0), pct_bD = 100 * (g6("bD", t, "b") - b0) / abs(b0), n = g6("m0", t, "n")) }))
write.csv(w6, file.path(WO, "verify_six_wide.csv"), row.names = FALSE)
sink()
