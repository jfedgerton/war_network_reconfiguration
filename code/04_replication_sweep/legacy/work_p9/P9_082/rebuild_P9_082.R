# P9_082 Pevehouse & Russett 2006, IO 60(4) "Democratic International Governmental Organizations Promote Peace" (10.1017/S0020818306060322).
#   TIE (joint memberships in densely democratic IGOs). REBUILD (protocol B; no archive).
# Article (staged PDF pp. 980-984, Tables 1-2): all dyads 1885-2000 with data, N 454,380; DV = onset of a MID with at least one fatality in
#   t+1 (Maoz dyadic MIDs); democratic igos = number of IGOs shared by i and j whose members' average Polity IV score is >= 7;
#   democracyL (lower Polity), dependenceS (lower trade/GDP), contiguity (land or < 150 miles water), distance (ln capital distance),
#   major power, cumulative mids (running sum of the dyad's MIDs), joint igos (all shared IGOs), cubic spline (two knots) on years since the
#   last onset. Estimator: logit (SEs in parentheses; clustering not stated - SEs clustered on dyad assumed).
# HEADLINE (fixed before estimating): Table 2 base model, democratic igos -0.079 (0.037) **, N 454,380.
#   Table 2 note: "All significance tests are one-tailed: *** p <= 0.01; ** p <= 0.05; * p <= 0.1" -> alpha_pub 0.10.
#   Verdict rule (PI, 2026-10-11 round 2): match = published sign, significant at alpha_pub when the published estimate is, and
#   rebuilt baseline within 0.5x-2x of the published coefficient.
# Rebuild inputs: DOR/O-R all-dyad file + derived trade file (work_demcap/D1/dor.dta: mzfatal1_orig = Maoz fatal MID onset in t+1, mzmid1_orig;
#   work_p4/X074/or_trade.rds: lowdem, lowdep, contig150, majpow); COW IGO v3 state-year file (work_p9/P9_082/state_year_formatv3.zip; full
#   membership = code 1); Polity IV 2014 polity2 for IGO member averages; capital distance ({peacesciencer}). Cumulative mids = running count of
#   MID onsets (mzmid1_orig, i.e. through t+1 lagged one year = onsets up to t) since 1816; splines: knots at the terciles of years since onset.
ROOT <- Sys.getenv("SWEEP_DATA"); OUT <- Sys.getenv("SWEEP_OUT"); CODE <- Sys.getenv("SWEEP_CODE")
stopifnot(nzchar(ROOT), nzchar(OUT), nzchar(CODE))
source(file.path(CODE, "work_p8/common_b3.R")); source(file.path(CODE, "work_p9/rebuild_common_p9r.R"))
set.seed(123)
W <- file.path(OUT, "work_p9/P9_082"); dir.create(W, recursive = TRUE, showWarnings = FALSE)
sink(file.path(W, "log_P9_082.txt"), split = TRUE)
# ---- 1. IGO memberships: shared democratic IGOs and shared IGOs per dyad-year ----
ig <- read.csv(unz(file.path(ROOT, "work_p9/P9_082/state_year_formatv3.zip"), "state_year_formatv3.csv"), check.names = FALSE)
stopifnot(nrow(ig) == 15557, all(c("ccode", "year", "state") %in% names(ig)))
ig <- ig[ig$year >= 1885 & ig$year <= 2000, ]
io <- setdiff(names(ig), c("ccode", "year", "state"))
p <- polity_cow("polity2")
ig <- ig %>% left_join(p, by = c("ccode", "year"))
mem <- as.matrix(ig[, io]) == 1; mem[is.na(mem)] <- FALSE
shared <- list()
for (y in sort(unique(ig$year))) {
  r <- which(ig$year == y); M <- mem[r, , drop = FALSE]; pol <- ig$polity2[r]
  avg <- colSums(M * ifelse(is.na(pol), 0, pol)) / colSums(M * !is.na(pol))   # mean Polity of members with Polity data
  dem <- !is.na(avg) & avg >= 7
  S_all <- tcrossprod(M * 1); S_dem <- tcrossprod(M[, dem, drop = FALSE] * 1)
  cc <- ig$ccode[r]; ij <- which(upper.tri(S_all), arr.ind = TRUE)
  a <- cc[ij[, 1]]; b <- cc[ij[, 2]]
  shared[[as.character(y)]] <- data.frame(statea = pmin(a, b), stateb = pmax(a, b), year = y, jigo = S_all[ij], demigo = S_dem[ij])
}
sh <- do.call(rbind, shared); rm(shared, mem); invisible(gc())
stopifnot(!anyDuplicated(sh[, c("statea", "stateb", "year")]))
cat("IGO state-years 1885-2000:", nrow(ig), " dyad-years with IGO counts:", nrow(sh), "\n")
# ---- 2. DOR/O-R dyads, cumulative MIDs, splines ----
d <- or_dyads(1816:2000) %>% arrange(dyadid, year) %>% group_by(dyadid) %>%
  mutate(cummid = cumsum(dplyr::lag(coalesce(mzmid1_orig, 0), default = 0))) %>% ungroup() %>% as.data.frame()
d <- d[d$year >= 1885, ]
cd <- suppressMessages(create_dyadyears(system = "cow", directed = FALSE, subset_years = 1885:2000) %>% add_capital_distance()) %>% as.data.frame()
n0 <- nrow(d)
d <- d %>% left_join(sh, by = c("statea", "stateb", "year")) %>% left_join(cd[, c("ccode1", "ccode2", "year", "capdist")], by = c("statea" = "ccode1", "stateb" = "ccode2", "year"))
stopifnot(nrow(d) == n0)
d <- py_calc(d, "dyadid", "mzfatal1_orig")
kn <- quantile(d$py, c(1/3, 2/3), type = 2); Tm <- max(d$py); cat("spline knots (terciles):", kn, " T:", Tm, "\n")
for (i in 1:2) d[[paste0("spl", i)]] <- (pmax(d$py - kn[i], 0)^3 - (Tm - kn[i]) / Tm * d$py^3) / 1000
d <- d %>% mutate(ldcap = log(capdist + 1), outcome_year = year + 1)
rhs0 <- c("demigo", "lowdem", "lowdep", "contig150", "ldcap", "majpow", "cummid", "jigo", "py", "spl1", "spl2")
vars <- c("mzfatal1_orig", rhs0, "dyadid")
cat("dyad-years 1885-2000:", nrow(d), " missing:\n"); print(sapply(d[, vars], function(x) sum(is.na(x))))
d0 <- d[complete.cases(d[, vars]), ]
cat("complete-case N:", nrow(d0), "(published 454,380)  fatal onsets", sum(d0$mzfatal1_orig), "\n")
TERM <- "demigo"
fitter <- function(r, data) fit_bin("mzfatal1_orig", r, data, "logit", "~dyadid", terms = TERM)
full <- fitter(rhs0, d0)
cat("\nREBUILT vs published: demigo b", full[[TERM]][["b"]], "(pub -0.079) se", full[[TERM]][["se"]], "(pub 0.037) N", attr(full, "n"), "(pub 454,380)\n")
ds <- common_sample(d, vars, "statea", "stateb", "outcome_year"); rm(d, d0); invisible(gc())
M <- run_six(rhs0, ds, fitter)
if (Sys.getenv("SWEEP_SENS") == "1") {
  source(file.path(CODE, "hist/sens_helper.R"))
  sens_run(id = "P9_082", dat = ds, y = "mzfatal1_orig", tie = TERM, ctrl = setdiff(rhs0, TERM), fe = NULL,
           fit = function(dat, add) fitter(c(rhs0, add), dat)[[TERM]])
}
saveRDS(list(M = M, full = full), file.path(W, "M_P9_082.rds"))
meta <- list(pid = "P9_082", cite = "Pevehouse & Russett 2006", journal = "IO", year = 2006, doi = "10.1017/S0020818306060322", family = "TIE")
nt <- sprintf(paste0("REBUILT (no archive): DOR/O-R all dyads 1885-2000, DV mzfatal1_orig (Maoz fatal MID onset t+1); democratic IGOs = shared full memberships ",
  "(COW IGO v3) in IGOs whose members' mean Polity IV 2014 polity2 >= 7; joint IGOs = all shared full memberships; lowdem, lowdep, contig150, ",
  "log capital distance, majpow, cumulative MID onsets, years since fatal onset + 2-knot cubic spline; logit clustered on dyad. Full-sample rebuilt ",
  "b=%.6f se=%.6f N=%d vs published -0.079 (0.037) N 454,380. Common sample N=%d (history at outcome year)."),
  full[[TERM]][["b"]], full[[TERM]][["se"]], attr(full, "n"), nrow(ds))
b0 <- M$m0[[TERM]][["b"]]; p0 <- M$m0[[TERM]][["p"]]; pub_p <- 2 * pnorm(-abs(-0.079 / 0.037))
v <- rebuild_verdict("P9_082", TERM, b0, p0, -0.079, NA, 0.037)   # round-2 rule, alpha_pub from registry.csv
write_p8(meta, TERM, "abstract: densely democratic IGOs far more likely to engender peaceful relations; Table 2 base model",
  "Table 2 base model (logit, fatal MID onset t+1, all dyads 1885-2000)", "logit, cubic spline (two knots), cumulative MIDs",
  "feglm logit, SE clustered on dyad, all six steps", -0.079, 0.037, 454380, TRUE, TRUE, M, "Y (cumulative MIDs + spline)", nt, verdict = v)
dir.create(file.path(OUT, "results_p9"), recursive = TRUE, showWarnings = FALSE)
src <- file.path(OUT, "results_p8", "P9_082__demigo.csv"); stopifnot(file.exists(src), file.rename(src, file.path(OUT, "results_p9", basename(src))))
sink()
