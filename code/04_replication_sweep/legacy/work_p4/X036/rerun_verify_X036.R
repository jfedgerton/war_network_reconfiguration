# VERIFY X036 Joshi, Maloy & Peterson 2015 JPR, Table I Model 1 (p.472): RE logit, random intercepts by dyad, N 116,527.
#   Published: initiator IDI(t-1) -0.536 (0.089)***; IDI x target Polity 0.008 (0.008) n.s.
# Copy of rerun_X036.R with ONE change: integration nAGQ = 12 (Stata xtlogit default intpoints(12), adaptive Gauss-Hermite;
#   do-file uses plain `xtlogit`, so the default applies) instead of Laplace (nAGQ = 1). Everything else unchanged.
# X036 Joshi, Maloy & Peterson 2015 JPR "Popular vs. elite democratic structures and international peace" (10.1177/0022343314567722)
# Data + do-file: PRIO JPR archive (2015/52/4 "Joshi, Maloy & Peterson 2015 replication.zip": dem_peace_final.dta, dempeace_do_file_final.do).
# Headline = Table 1 Model 1: xtlogit init_MID c.l.IDI_weighted##c.l.polity2a l.caprat lndist l.s_wt_glo c.init_MID_yr##c.init_MID_yr##c.init_MID_yr
#   coldwar year2001   (Stata xtlogit default = random-effects logit on the xtset panel; directed dyads, established democracy = side A, 1961-2001)
# KEY (a priori, abstract): L.IDI_weighted ("more popular democracies are less warlike with respect to all other regimes, not just other
#   democracies") -> constituent of IDI x target polity; also report the product term (the 'not just other democracies' part).
# Estimation: lme4::glmer random intercept by panel, nAGQ = 12 (matches Stata xtlogit intpoints(12); see lines 3-4).
ROOT <- Sys.getenv("SWEEP_DATA"); OUT <- Sys.getenv("SWEEP_OUT"); CODE <- Sys.getenv("SWEEP_CODE")
stopifnot(nzchar(ROOT), nzchar(OUT), nzchar(CODE))
source(file.path(CODE, "vb2/vcommon.R"))
suppressPackageStartupMessages({library(readstata13); library(lme4)})
WO <- file.path(OUT, "work_p4/X036"); dir.create(WO, recursive = TRUE, showWarnings = FALSE)  # port: outputs
W <- file.path(ROOT, "work_p4/X036"); sink(file.path(WO, "log_verify_X036.txt"), split = TRUE)
fn <- file.path(WO, "Appendix and replication data/dem_peace_final.dta")  # port: zip member unzipped into OUT
if (!file.exists(fn)) unzip(file.path(W, "r.zip"), files = "Appendix and replication data/dem_peace_final.dta", exdir = WO)
stopifnot(file.exists(fn))
h <- read.dta13(fn, select.rows = c(1, 2), convert.factors = FALSE); nm <- names(h)
ef <- attr(h, "expansion.fields"); tsl <- ef[sapply(ef, function(x) any(grepl("^_TS|^iis$|^tis$", x[2])))]; print(tsl); print(nm)
pv <- NULL; for (x in tsl) if (x[2] %in% c("_TSpanel", "iis")) pv <- x[3]
tv <- "year"; for (x in tsl) if (x[2] %in% c("_TStvar", "tis")) tv <- x[3]
c1 <- intersect(c("ccode1", "ccodea", "statea", "ccode_a", "state1"), nm)[1]; c2 <- intersect(c("ccode2", "ccodeb", "stateb", "ccode_b", "state2"), nm)[1]
if (is.null(pv) || !(pv %in% nm)) pv <- intersect(c("ddyad", "dyadid", "dyad", "directed_dyad"), nm)[1]
cat("panel:", pv, " time:", tv, " c1:", c1, " c2:", c2, "\n"); stopifnot(!is.na(c1), !is.na(c2), !is.na(pv))
V <- c("init_MID", "IDI_weighted", "polity2a", "caprat", "lndist", "s_wt_glo", "init_MID_yr", "coldwar", "year2001")
stopifnot(all(V %in% nm))
d <- read.dta13(fn, select.cols = unique(c(V, pv, tv, c1, c2)), convert.factors = FALSE)
names(d)[names(d) == pv] <- "panel"; names(d)[names(d) == tv] <- "year"
if (pv == "year") stop("panel var is year?")
stopifnot(anyDuplicated(d[, c("panel", "year")]) == 0)
cat("rows:", nrow(d), " years:", range(d$year), " panels:", length(unique(d$panel)), "\n")
# Stata L. operator: value at (panel, year-1), NA if that row is absent
lagv <- function(v) { k <- d[, c("panel", "year", v)]; k$year <- k$year + 1; names(k)[3] <- "x"; left_join(d[, c("panel", "year")], k, by = c("panel", "year"))$x }
for (v in c("IDI_weighted", "polity2a", "caprat", "s_wt_glo")) d[[paste0("L_", v)]] <- lagv(v)
d$L_IDIxPol <- d$L_IDI_weighted * d$L_polity2a
d$t1 <- d$init_MID_yr; d$t2 <- d$t1^2; d$t3 <- d$t1^3
print(summary(d[, c("init_MID", "L_IDI_weighted", "L_polity2a", "L_caprat", "lndist", "L_s_wt_glo", "t1")]))
f0 <- init_MID ~ L_IDI_weighted + L_polity2a + L_IDIxPol + L_caprat + lndist + L_s_wt_glo + t1 + t2 + t3 + coldwar + year2001
ff_re <- function(f, data) {
  data <- data[complete.cases(data[, c(all.vars(f), "panel")]), ]
  m <- glmer(update(f, . ~ . + (1 | panel)), data = data, family = binomial, nAGQ = 12,
             control = glmerControl(optimizer = "bobyqa", calc.derivs = FALSE))
  cat("   [RE sd =", round(sqrt(unlist(VarCorr(m))), 4), " conv msgs:", length(m@optinfo$conv$lme4$messages), "]\n")
  ct <- summary(m)$coefficients[, c(1, 2, 4), drop = FALSE]; attr(ct, "n") <- nrow(data); rm(m); gc(); ct
}
terms <- c("L_IDI_weighted", "L_IDIxPol")
dir.create(file.path(WO, "verify_cache"), recursive = TRUE, showWarnings = FALSE)
rep <- vfit(f0, d[complete.cases(d[, c(all.vars(f0), "panel")]), ], terms, ff_re, file.path(WO, "verify_cache/v_repro.rds"))
cat("repro N:", rep$n, "(published 116,527)\n")
ds <- vcommon_sample(d, f0, "panel", c1, c2, "year"); rm(d); gc()
M <- vrun_six(f0, ds, terms, ff_re, WO)
vsave(M, terms, WO, rep)
# port: one row per term from verify_six.csv (long) for ID-keyed extraction; no estimation (same arithmetic as vb2/assemble.py from_six)
V6 <- read.csv(file.path(WO, "verify_six.csv")); g6 <- function(m, t, k) { x <- V6[V6$model == m & V6$term == t, k]; stopifnot(length(x) == 1); x }
w6 <- do.call(rbind, lapply(terms, function(t) { b0 <- g6("m0", t, "b")
  data.frame(term = t, b0 = b0, p0 = g6("m0", t, "p"), b_dm = g6("mDM", t, "b"), p_dm = g6("mDM", t, "p"), b_26 = g6("m26", t, "b"), p_26 = g6("m26", t, "p"),
             pct_bC = 100 * (g6("bC", t, "b") - b0) / abs(b0), pct_bD = 100 * (g6("bD", t, "b") - b0) / abs(b0), n = g6("m0", t, "n")) }))
write.csv(w6, file.path(WO, "verify_six_wide.csv"), row.names = FALSE)
sink()
