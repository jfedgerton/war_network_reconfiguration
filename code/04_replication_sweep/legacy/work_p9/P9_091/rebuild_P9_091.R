# P9_091 Mattes & Vonnahme 2010, JOP 72(4) "Contracting for Peace: Do Nonaggression Pacts Reduce Conflict?" (10.1017/s0022381610000447).
#   TIE (nonaggression pacts). REBUILD (protocol B; replication files "available from the authors", none posted).
# Article (staged PDF pp. 931-935, Tables 1-2, fns 8-13): treatment = politically relevant (PR) dyad-years 1920-2001 in which a bilateral
#   nonaggression or nonaggression/neutrality pact is signed (ATOP 3.0), excluding pacts with defense, offense or consultation obligations
#   and renegotiations/renewals; controls = PR dyad-years with no bilateral or multilateral alliance of any kind; PR = separated by <= 400
#   miles of water or at least one major power. Genetic matching (Diamond & Sekhon) on IIS (Crescenzi-Enterline interstate interaction
#   score), strategic rivalry (Thompson), ln CINC ratio, major power, land contiguity, weakest-link trade dependence (O-R + Gleditsch),
#   lower and higher Polity IV, plus the logit propensity-score linear predictor -> 56 treated matched to 56 controls. Outcome = number of
#   unique MID onsets (initiators or joiners; Maoz dyadic MIDs) in the 5 (10) years after signing, truncated to the pact's life (fn 12),
#   with the matched control observed over the same span. Table 2: OLS on the matched sample with the imbalanced covariates.
# HEADLINE (fixed before estimating): Table 2, All MIDs 5yrs, Nonaggression Pact -0.358 (0.195), p = 0.069, N 112. Two-tailed tests;
#   the text calls this "statistically significant" (p .069) -> alpha_pub 0.10.
#   Verdict rule (PI, 2026-10-11 round 2): match = published sign, significant at alpha_pub when the published estimate is, and
#   rebuilt baseline within 0.5x-2x of the published coefficient.
# Rebuild inputs and departures (disclosed): ATOP 5.1 alliance-level file (work_p9/P9_091/atop5_1_csv.zip; ATOP 3.0 not posted); a renewal
#   is coded as any ATOP alliance between the pair in force in the year before signing; GML dyadic MIDs ({peacesciencer} gml_dirdisp);
#   strategic rivalry = Thompson-Dreyer ({peacesciencer}); trade dependence = O-R/Gleditsch lowdep (work_p4/X074/or_trade.rds, to 2000);
#   IIS is APPROXIMATED as -sum over past MID dyad-years of (max hostility / 5) * 0.5^(years elapsed / 10) (the Crescenzi-Enterline
#   weights are not reported in the article). GenMatch (Matching) with seed 123, default population/generations, 1:1 without replacement.
ROOT <- Sys.getenv("SWEEP_DATA"); OUT <- Sys.getenv("SWEEP_OUT"); CODE <- Sys.getenv("SWEEP_CODE")
stopifnot(nzchar(ROOT), nzchar(OUT), nzchar(CODE))
source(file.path(CODE, "work_p8/common_b3.R")); source(file.path(CODE, "work_p9/rebuild_common_p9r.R"))
suppressPackageStartupMessages(library(Matching)); cat("Matching", as.character(packageVersion("Matching")), "\n")
set.seed(123)
W <- file.path(OUT, "work_p9/P9_091"); dir.create(W, recursive = TRUE, showWarnings = FALSE)
sink(file.path(W, "log_P9_091.txt"), split = TRUE)
# ---- 1. ATOP 5.1: bilateral nonaggression pacts without defense/offense/consultation obligations, signed 1920-2001 ----
a <- read.csv(unz(file.path(ROOT, "work_p9/P9_091/atop5_1_csv.zip"), "ATOP 5.1 (.csv)/atop5_1a.csv"))
stopifnot(nrow(a) > 700, all(c("atopid", "begyr", "endyr", "bilat", "nonagg", "defense", "offense", "consul", "mem1", "mem2") %in% names(a)))
nap <- a %>% filter(bilat == 1, nonagg == 1, defense == 0, offense == 0, consul == 0, begyr >= 1920, begyr <= 2001) %>%
  transmute(atopid, statea = pmin(mem1, mem2), stateb = pmax(mem1, mem2), year = begyr, endyr = ifelse(is.na(endyr) | endyr == 0, 9999, endyr))
cat("ATOP 5.1 bilateral NA pacts (no def/off/consul) signed 1920-2001:", nrow(nap), "\n")
# ---- 2. PR dyad-years 1816-2001 (covariates) with ATOP flags ----
d <- create_dyadyears(system = "cow", directed = FALSE, subset_years = 1816:2001, mry = FALSE)
d <- suppressMessages(d %>% add_contiguity() %>% add_nmc() %>% add_cow_majors() %>% add_atop_alliance() %>% add_strategic_rivalries()) %>% as.data.frame()
stopifnot(!anyDuplicated(d[, c("ccode1", "ccode2", "year")]))
rv <- grep("^ongoingrivalry$|^ongoing", names(d), value = TRUE); cat("rivalry vars:", rv, "\n"); stopifnot(length(rv) >= 1)
d$rival <- as.integer(coalesce(d[[rv[1]]], 0) == 1)
d$ally_any <- as.integer((coalesce(d$atop_defense, 0) + coalesce(d$atop_offense, 0) + coalesce(d$atop_neutral, 0) + coalesce(d$atop_nonagg, 0) +
                          coalesce(d$atop_consul, 0)) > 0)
p <- polity_cow("polity2")
d <- d %>% rename(statea = ccode1, stateb = ccode2) %>%
  left_join(rename(p, p2a = polity2), by = c("statea" = "ccode", "year")) %>% left_join(rename(p, p2b = polity2), by = c("stateb" = "ccode", "year"))
tr <- readRDS(file.path(ROOT, "work_p4/X074/or_trade.rds"))[, c("statea", "stateb", "year", "lowdep")]
stopifnot(!anyDuplicated(tr[, 1:3])); n0 <- nrow(d)
d <- d %>% left_join(tr, by = c("statea", "stateb", "year")); stopifnot(nrow(d) == n0)
# MIDs: unique dyadic dispute onsets (initiators or joiners) and max hostility per dyad-year, GML
g <- as.data.frame(gml_dirdisp) %>% filter(ccode1 < ccode2) %>% group_by(dispnum, ccode1, ccode2) %>%
  summarise(onset = min(year), host = max(hostlev, na.rm = TRUE), .groups = "drop")
gy <- as.data.frame(gml_dirdisp) %>% filter(ccode1 < ccode2) %>% group_by(statea = ccode1, stateb = ccode2, year) %>%
  summarise(hy = max(hostlev, na.rm = TRUE), .groups = "drop")
d <- d %>% left_join(gy, by = c("statea", "stateb", "year")) %>% mutate(hy = coalesce(hy, 0)) %>% arrange(statea, stateb, year)
# approximated IIS: decayed sum of past hostility (years before t), half-life 10 years
d <- d %>% group_by(statea, stateb) %>% mutate(iis = {
  h <- hy / 5; y <- year; s <- numeric(length(y)); acc <- 0; prev <- NA
  for (i in seq_along(y)) { if (i > 1) acc <- (acc + h[i - 1]) * 0.5^((y[i] - y[i - 1]) / 10); s[i] <- -acc }; s }) %>% ungroup() %>% as.data.frame()
d <- d %>% mutate(pr = as.integer(conttype %in% 1:5 | cowmaj1 == 1 | cowmaj2 == 1), majpow = as.integer(cowmaj1 == 1 | cowmaj2 == 1),
                  contig = as.integer(conttype == 1), caprat = log(pmax(cinc1, cinc2) / pmin(cinc1, cinc2)),
                  pollow = pmin(p2a, p2b), polhigh = pmax(p2a, p2b))
# ---- 3. treated and control pools ----
prev <- d %>% transmute(statea, stateb, year = year + 1, ally_prev = ally_any)
nap <- nap %>% left_join(prev, by = c("statea", "stateb", "year")) %>% left_join(d[, c("statea", "stateb", "year", "pr")], by = c("statea", "stateb", "year"))
cat("of which renewals (any ATOP alliance in force year before):", sum(nap$ally_prev %in% 1), " not PR:", sum(!nap$pr %in% 1), "\n")
nap <- nap %>% filter(!ally_prev %in% 1, pr %in% 1) %>% distinct(statea, stateb, year, .keep_all = TRUE)
covs <- c("iis", "rival", "caprat", "majpow", "contig", "lowdep", "pollow", "polhigh")
tre <- nap %>% dplyr::select(statea, stateb, year, endyr) %>% left_join(d[, c("statea", "stateb", "year", covs)], by = c("statea", "stateb", "year")) %>% mutate(Tr = 1L)
ctl <- d %>% filter(pr == 1, ally_any == 0, year >= 1920, year <= 2001) %>% dplyr::select(statea, stateb, year, all_of(covs)) %>% mutate(endyr = NA, Tr = 0L) %>%
  anti_join(nap[, c("statea", "stateb", "year")], by = c("statea", "stateb", "year"))
pool <- rbind(tre[, c("statea", "stateb", "year", "endyr", covs, "Tr")], ctl[, c("statea", "stateb", "year", "endyr", covs, "Tr")])
cat("treated PR NA pacts (before covariate completeness):", nrow(tre), " complete:", sum(complete.cases(tre[, covs])), "\n")
pool <- pool[complete.cases(pool[, covs]) & is.finite(pool$caprat), ]
cat("pool: treated", sum(pool$Tr), " controls", sum(pool$Tr == 0), "(published: 56 treated)\n")
# ---- 4. propensity score + GenMatch (seed 123), cached ----
ps <- glm(reformulate(covs, "Tr"), data = pool, family = binomial()); pool$lp <- predict(ps, type = "link")
X <- as.matrix(pool[, c(covs, "lp")])
gf <- file.path(W, "genmatch_P9_091.rds")
if (file.exists(gf) && Sys.getenv("P9_091_REFIT") != "1") { gm <- readRDS(gf); cat("GenMatch weights read from cache", gf, "\n") } else {
  t0 <- Sys.time()
  gm <- GenMatch(Tr = pool$Tr, X = X, BalanceMatrix = X, estimand = "ATT", M = 1, replace = FALSE, ties = FALSE,
                 unif.seed = 123, int.seed = 123, print.level = 0)
  cat("GenMatch minutes:", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), "\n"); saveRDS(gm, gf)
}
mt <- Match(Tr = pool$Tr, X = X, estimand = "ATT", M = 1, replace = FALSE, ties = FALSE, Weight.matrix = gm)
stopifnot(length(mt$index.treated) == sum(pool$Tr), !anyDuplicated(mt$index.control))
print(MatchBalance(reformulate(c(covs, "lp"), "Tr"), data = pool, match.out = mt, nboots = 0, print.level = 0)$AMsmallest.p.value)
tm <- pool[mt$index.treated, ]; cm <- pool[mt$index.control, ]
L5 <- pmin(5, pmax(1, tm$endyr - tm$year)); L10 <- pmin(10, pmax(1, tm$endyr - tm$year))
cm$L5 <- L5; tm$L5 <- L5; cm$L10 <- L10; tm$L10 <- L10
ds0 <- rbind(tm, cm); ds0$nap <- ds0$Tr
cnt <- function(sa, sb, y, L, minhost = 0) {
  k <- g$ccode1 == sa & g$ccode2 == sb & g$onset >= y + 1 & g$onset <= y + L & g$host >= minhost; sum(k) }
ds0$mid5 <- mapply(cnt, ds0$statea, ds0$stateb, ds0$year, ds0$L5)
ds0$vmid5 <- mapply(cnt, ds0$statea, ds0$stateb, ds0$year, ds0$L5, MoreArgs = list(minhost = 4))
ds0$mid10 <- mapply(cnt, ds0$statea, ds0$stateb, ds0$year, ds0$L10)
cat("matched N:", nrow(ds0), " mean MIDs 5y treated", mean(ds0$mid5[ds0$nap == 1]), "control", mean(ds0$mid5[ds0$nap == 0]),
    " (published t-test means 0.55 vs 0.91)\n")
ds0$outcome_year <- ds0$year + 1
# ---- 5. post-matching OLS (Table 2 controls), six steps ----
TERM <- "nap"; rhs0 <- c("nap", "caprat", "pollow", "polhigh", "iis", "lowdep")
fitter <- function(r, data) {
  m <- lm(reformulate(r, "mid5"), data = data); ct <- summary(m)$coefficients
  cat("\n---- mid5 ~", paste(r, collapse = " + "), "| N =", nobs(m), "\n"); print(round(ct, 6))
  out <- list(nap = c(b = ct["nap", 1], se = ct["nap", 2], p = ct["nap", 4])); attr(out, "n") <- nobs(m); out
}
full <- fitter(rhs0, ds0)
cat("\nREBUILT vs published: nap b", full[[TERM]][["b"]], "(pub -0.358) se", full[[TERM]][["se"]], "(pub 0.195) N", attr(full, "n"), "(pub 112)\n")
vars <- c("mid5", rhs0)
ds <- common_sample(ds0, vars, "statea", "stateb", "outcome_year")
M <- run_six(rhs0, ds, fitter)
if (Sys.getenv("SWEEP_SENS") == "1") {
  source(file.path(CODE, "hist/sens_helper.R"))
  sens_run(id = "P9_091", dat = ds, y = "mid5", tie = TERM, ctrl = setdiff(rhs0, TERM), fe = NULL,
           fit = function(dat, add) fitter(c(rhs0, add), dat)[[TERM]])
}
saveRDS(list(M = M, full = full, ds = ds0), file.path(W, "M_P9_091.rds"))
meta <- list(pid = "P9_091", cite = "Mattes & Vonnahme 2010", journal = "JOP", year = 2010, doi = "10.1017/s0022381610000447", family = "TIE")
nt <- sprintf(paste0("REBUILT (no archive): ATOP 5.1 bilateral nonaggression pacts (no defense/offense/consultation; renewal = any ATOP alliance in force the ",
  "year before) signed by PR dyads (contiguity 1-5 or COW major) 1920-2001; controls = PR dyad-years with no ATOP alliance; GenMatch (seed 123, 1:1, ",
  "no replacement) on approximated IIS (decayed GML hostility, half-life 10y), Thompson-Dreyer rivalry, ln CINC ratio, major power, land contiguity, ",
  "O-R lowdep, lower/higher polity2 and the logit linear predictor; DV = unique GML MID onsets in years t+1..t+min(5, pact life); OLS iid SE. ",
  "Full matched sample b=%.6f se=%.6f N=%d vs published -0.358 (0.195) N 112. Common sample N=%d (history at year after signing)."),
  full[[TERM]][["b"]], full[[TERM]][["se"]], attr(full, "n"), nrow(ds))
b0 <- M$m0[[TERM]][["b"]]; p0 <- M$m0[[TERM]][["p"]]; pub_p <- 0.069
v <- if (sign(b0) == -1 && (pub_p >= .10 || p0 < .10) && b0 / -0.358 >= 0.5 && b0 / -0.358 <= 2) "REBUILT_MATCH" else "REBUILT_MISMATCH"
cat("round-2 verdict:", v, " ratio", round(b0 / -0.358, 3), " p0", signif(p0, 3), "\n")
write_p8(meta, TERM, "abstract: nonaggression pacts reduce conflict between signatories; Table 2 All MIDs 5yrs",
  "Table 2, All MIDs 5yrs (post-matching OLS, N 112)", "genetic matching + OLS on matched sample (IIS control)",
  "lm OLS iid SE on the GenMatch sample, all six steps", -0.358, 0.195, 112, TRUE, TRUE, M, "Y (IIS + strategic rivalry, matched; IIS in OLS)", nt, verdict = v)
dir.create(file.path(OUT, "results_p9"), recursive = TRUE, showWarnings = FALSE)
src <- file.path(OUT, "results_p8", "P9_091__nap.csv"); stopifnot(file.exists(src), file.rename(src, file.path(OUT, "results_p9", basename(src))))
sink()
