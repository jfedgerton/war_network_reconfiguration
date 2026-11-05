# P9_033 Hegre 2004 JCR 48(3) "Size Asymmetry, Trade, and Militarized Conflict" (10.1177/0022002704263954). REBUILD (protocol B).
# No archive on disk (paper points to a defunct personal page). Spec read from the PDF (Table 1, p. 423; pp. 417-421):
#   Cox regression in CALENDAR time (Raknerud & Hegre 1997): at each day on which a fatal militarized action occurs, the risk set is all
#   directed dyads in the system; 5% of the non-failure observations sampled with weight 20; robust SEs. 1951-1992, 438 failures.
#   DV: directed fatal militarized action (actor's militarized action in a MID with >= 1 fatality; Maoz dyadic MID, ROLEA/ROLEB);
#   initiator side at the dyadic start date t, target side at t + epsilon.
#   Model II RHS: s, s^2, e, e*s, e*s^2, democracy actor, democracy target (Polity 0-10, centered), their product, ln GDPpc actor,
#   ln GDPpc target (lag 1), distance (ln), contiguity (land or < 25 nm water), dyad size ln(GDPa + GDPb) (lag 1),
#   Prx(independence) (half-life 2,192 days), Prx(actor action) (2,192), Prx(target action) (14 days).
#   s = GDPa/(GDPa + GDPb) (lag 1, centred on 0.5 for the main terms); e = bilateral trade x (GDPa + GDPb)/(GDPa x GDPb) (lag 1).
# HEADLINE / KEY (fixed from the PDF before estimating): e*s^2 in Model II, published 310.5 (SE 89.7, p < .01). Abstract: trade "reduces
#   the incentives for conflict but ... most clearly seen in relatively symmetric dyads"; Proposition 3 is carried by e*s^2 (and s^2).
#   The e*s^2 coefficient is invariant to centring e or s (centring only moves lower-order terms).
# Rebuild inputs: GML MID 2.2.1 directed dispute-years + participant start dates (peacesciencer), COW directed dyad-years,
#   COW trade (peacesciencer extdata cow_trade_ddy.rds, smoothed flows), SDP GDP / GDPpc (peacesciencer), Polity IV democ
#   (p4v2014.csv already on disk under work_p8/P8_049), COW contiguity (conttype <= 3) and capital distance (peacesciencer).
# Estimator: survival::coxph, strata = event time (= calendar-time risk sets), Breslow ties, case weights (1 failures, 20 sampled
#   non-failures), robust SE clustered on directed dyad. Same in all six steps. Product key term -> no sensitivity hook.
CODE <- Sys.getenv("SWEEP_CODE"); stopifnot(nzchar(CODE))
source(file.path(CODE, "work_p8/common_b3.R"))   # ROOT, OUT, dplyr, fixest, add_history, common_sample, run_six, write_p8; set.seed(123)
suppressPackageStartupMessages({library(peacesciencer); library(survival)})
set.seed(123)
W <- file.path(OUT, "work_p9/P9_033"); dir.create(W, recursive = TRUE, showWarnings = FALSE)
sink(file.path(W, "log_P9_033.txt"), split = TRUE)
Y0 <- 1951; Y1 <- 1992

# ---- 1. directed fatal militarized actions with dates (all years, for the proximity terms) ----
gp <- as.data.frame(gml_part)
stopifnot(all(gp$stmon %in% 1:12), all(gp$stday == -9 | gp$stday %in% 1:31))
gp <- gp %>% mutate(sd = as.Date(sprintf("%d-%02d-%02d", styear, stmon, ifelse(stday == -9, dummy_stday, stday)))) %>%
  group_by(dispnum, ccode) %>% summarise(sd = min(sd), .groups = "drop")
stopifnot(!anyNA(gp$sd))
g <- as.data.frame(gml_dirdisp)
stopifnot(all(g$sidea1 != g$sidea2))
g <- g %>% group_by(dispnum, ccode1, ccode2) %>% arrange(year, .by_group = TRUE) %>% slice(1) %>% ungroup() %>%   # dyad's first dispute-year
  left_join(gp %>% rename(sd1 = sd), by = c("dispnum", "ccode1" = "ccode")) %>% left_join(gp %>% rename(sd2 = sd), by = c("dispnum", "ccode2" = "ccode"))
stopifnot(!anyNA(g$sd1), !anyNA(g$sd2))
# event ccode1 -> ccode2: ccode1 used force (hostlev >= 4) in a dyadic dispute with >= 1 fatality on either side
ev <- g %>% filter(hostlev1 >= 4, pmax(fatality1, fatality2) >= 1) %>%
  mutate(dd = pmax(sd1, sd2), tt = as.numeric(dd) + ifelse(sidea1 == 1, 0, 0.5), year = as.numeric(format(dd, "%Y"))) %>%
  transmute(ccode1, ccode2, tt, year) %>% distinct(ccode1, ccode2, tt, .keep_all = TRUE) %>% as.data.frame()
evs <- ev[ev$year >= Y0 & ev$year <= Y1, ]
cat("fatal directed actions, all years:", nrow(ev), "  ", Y0, "-", Y1, ":", nrow(evs), "(paper: 438 failures)  distinct event times:", length(unique(evs$tt)), "\n")

# ---- 2. directed dyad-year covariates 1951-1992 ----
dy <- create_dyadyears(directed = TRUE, subset_years = Y0:Y1) %>% add_contiguity() %>% add_capital_distance() %>% as.data.frame()
n_dy <- nrow(dy); stopifnot(n_dy == 769966, anyDuplicated(dy[, c("ccode1", "ccode2", "year")]) == 0)
sdp <- as.data.frame(cow_sdp_gdp) %>% transmute(ccode, year = year + 1, gdp = exp(wbgdp2011est), lpc = wbgdppc2011est)   # lag 1
stopifnot(anyDuplicated(sdp[, c("ccode", "year")]) == 0)
dy <- dy %>% left_join(sdp %>% rename(gdp1 = gdp, lpc1 = lpc), by = c("ccode1" = "ccode", "year")) %>%
  left_join(sdp %>% rename(gdp2 = gdp, lpc2 = lpc), by = c("ccode2" = "ccode", "year")); stopifnot(nrow(dy) == n_dy)
trp <- system.file("extdata", "cow_trade_ddy.rds", package = "peacesciencer"); stopifnot(nzchar(trp), file.exists(trp))
tr <- readRDS(trp) %>% as.data.frame() %>% filter(ccode1 < ccode2) %>% transmute(   # file holds both orders of each pair
  a = ccode1, b = ccode2, year = year + 1, trade = smoothflow1 + smoothflow2)   # lag 1, USD mn
stopifnot(all(tr$a < tr$b), anyDuplicated(tr[, c("a", "b", "year")]) == 0)
dy <- dy %>% mutate(a = pmin(ccode1, ccode2), b = pmax(ccode1, ccode2)) %>% left_join(tr, by = c("a", "b", "year")); stopifnot(nrow(dy) == n_dy)
p4 <- read.csv(file.path(ROOT, "work_p8/P8_049/p4v2014.csv"))
p4 <- p4 %>% transmute(ccode = ccode, year = year, democ = ifelse(democ %in% 0:10, democ, NA_real_)) %>%
  group_by(ccode, year) %>% summarise(democ = if (all(is.na(democ))) NA_real_ else mean(democ, na.rm = TRUE), .groups = "drop")   # transition years: mean
dy <- dy %>% left_join(p4 %>% rename(dem1 = democ), by = c("ccode1" = "ccode", "year")) %>%
  left_join(p4 %>% rename(dem2 = democ), by = c("ccode2" = "ccode", "year")); stopifnot(nrow(dy) == n_dy)
dy <- dy %>% mutate(s = gdp1 / (gdp1 + gdp2), sc = s - 0.5, sc2 = sc^2,
  e = 1e6 * trade * (gdp1 + gdp2) / (gdp1 * gdp2), e_s = e * sc, e_s2 = e * sc2,
  contig = as.numeric(conttype %in% 1:3), ldist = log(capdist), dsize = log(gdp1 + gdp2), lpc_a = lpc1, lpc_t = lpc2,
  ddyad = ccode1 * 1000 + ccode2)
dm <- mean(c(dy$dem1, dy$dem2), na.rm = TRUE)
dy <- dy %>% mutate(dem_a = dem1 - dm, dem_t = dem2 - dm, dem_at = dem_a * dem_t)
cat("dyad-years:", n_dy, " missing: s", sum(is.na(dy$s)), " trade", sum(is.na(dy$trade)), " democ", sum(is.na(dy$dem_at)), " capdist", sum(is.na(dy$capdist)),
    "  democ centring mean", round(dm, 3), "\n")
cat("e quantiles (non-missing):\n"); print(signif(quantile(dy$e, c(.5, .9, .99, 1), na.rm = TRUE), 3))

# ---- 3. calendar-time risk sets: every event time; failures + 5% sample of non-failures (weight 20) ----
tts <- sort(unique(evs$tt)); ty <- as.numeric(format(as.Date(floor(tts), origin = "1970-01-01"), "%Y"))
idx_by_year <- split(seq_len(n_dy), dy$year)
rows <- vector("list", length(tts))
for (i in seq_along(tts)) {
  yi <- as.character(ty[i]); ii <- idx_by_year[[yi]]; stopifnot(length(ii) > 0)
  fk <- paste(evs$ccode1[evs$tt == tts[i]], evs$ccode2[evs$tt == tts[i]])
  isf <- paste(dy$ccode1[ii], dy$ccode2[ii]) %in% fk
  if (sum(isf) != length(fk)) cat("  event time", tts[i], ": failing dyads not in system list:", length(fk) - sum(isf), "\n")
  keep <- ii[isf | runif(length(ii)) < 0.05]
  rows[[i]] <- data.frame(r = keep, tt = tts[i], ev = as.numeric(isf[match(keep, ii)]))
}
rs <- do.call(rbind, rows)
d <- cbind(dy[rs$r, ], rs[, c("tt", "ev")]); d$w <- ifelse(d$ev == 1, 1, 20)
cat("risk-set rows:", nrow(d), " failures:", sum(d$ev), " event times:", length(tts), "\n")
stopifnot(nrow(d) > 2e5, sum(d$ev) > 300)

# ---- 4. proximity terms at the event time ----
# Prx(actor action): last fatal action actor -> target strictly before tt; Prx(target action): target -> actor. 0 if none.
lastbefore <- function(a, b, t, E) {
  out <- rep(NA_real_, length(a)); key <- paste(a, b); ek <- split(E$tt, paste(E$ccode1, E$ccode2))
  hit <- which(key %in% names(ek))
  for (j in hit) { x <- ek[[key[j]]]; x <- x[x < t[j]]; if (length(x)) out[j] <- max(x) }
  out
}
la <- lastbefore(d$ccode1, d$ccode2, d$tt, ev); lt <- lastbefore(d$ccode2, d$ccode1, d$tt, ev)
d$prx_act <- ifelse(is.na(la), 0, 2^(-(d$tt - la) / 2192)); d$prx_tgt <- ifelse(is.na(lt), 0, 2^(-(d$tt - lt) / 14))
# Prx(independence): days since the younger state's (current spell) COW system entry
cs <- as.data.frame(cow_states) %>% transmute(ccode, ent = as.numeric(as.Date(sprintf("%d-%02d-%02d", styear, stmonth, stday))),
  ext = as.numeric(as.Date(sprintf("%d-%02d-%02d", endyear, endmonth, endday))))
entry_at <- function(cc, t) { o <- rep(NA_real_, length(cc)); for (k in unique(cc)) { j <- which(cc == k); s <- cs[cs$ccode == k, ]
  for (q in seq_len(nrow(s))) o[j][t[j] >= s$ent[q] - 366] <- s$ent[q] }; o }   # -366: system members entering later in the event year
e1 <- entry_at(d$ccode1, d$tt); e2 <- entry_at(d$ccode2, d$tt); stopifnot(!anyNA(e1), !anyNA(e2))
d$prx_ind <- 2^(-pmax(d$tt - pmax(e1, e2), 0) / 2192)
cat("prx means: actor", round(mean(d$prx_act), 4), " target", round(mean(d$prx_tgt), 4), " indep", round(mean(d$prx_ind), 4),
    "  among failures: actor", round(mean(d$prx_act[d$ev == 1]), 3), " target", round(mean(d$prx_tgt[d$ev == 1]), 3), "\n")
d <- as.data.frame(d); d$t1 <- 1

# ---- 5. headline fit (Table 1 Model II) ----
rhs0 <- c("sc", "sc2", "e", "e_s", "e_s2", "dem_a", "dem_t", "dem_at", "lpc_a", "lpc_t", "ldist", "contig", "dsize", "prx_ind", "prx_act", "prx_tgt")
K <- "e_s2"
fitter <- function(r, data) {
  f <- as.formula(paste("Surv(t1, ev) ~", paste(r, collapse = " + "), "+ strata(tt)"))
  m <- withCallingHandlers(coxph(f, data = data, weights = w, ties = "breslow", cluster = ddyad, control = coxph.control(iter.max = 200)),
    warning = function(w) stop("coxph did not converge: ", conditionMessage(w)))   # fail loudly
  ct <- summary(m)$coefficients; stopifnot(all(c(K, r) %in% rownames(ct)))
  cat("\n---- coxph | N =", m$n, " failures =", m$nevent, " clusters =", length(unique(data$ddyad)), "\n")
  print(signif(ct[, c("coef", "robust se", "Pr(>|z|)")], 5))
  out <- list(c(b = ct[K, "coef"], se = ct[K, "robust se"], p = ct[K, "Pr(>|z|)"])); names(out) <- K
  attr(out, "n") <- m$n; out
}
vars <- c("t1", "ev", "w", "tt", "ddyad", rhs0)
d0 <- d[complete.cases(d[, vars]), ]
full <- fitter(rhs0, d0)
cat("\nREBUILT vs published (Model II): e_s2 b", full[[K]][["b"]], "se", full[[K]][["se"]], "p", full[[K]][["p"]], "N", attr(full, "n"),
    "failures", sum(d0$ev), "| published 310.5 (89.7), p < .01, 438 failures, N not reported\n")
ds <- common_sample(d, vars, "ccode1", "ccode2", "year")
M <- run_six(rhs0, ds, fitter)
saveRDS(list(M = M, full = full, n_used = nrow(ds)), file.path(W, "M_P9_033.rds"))
# no sensitivity hook: key term is a product (e x s^2) in a Cox model

meta <- list(pid = "P9_033", cite = "Hegre 2004", journal = "JCR", year = 2004, doi = "10.1177/0022002704263954", family = "TIE")
art <- "Cox regression in calendar time (event-day risk sets of all directed dyads), 5% non-failure sampling weighted 20, robust SE"
used <- "coxph, strata = event time, Breslow, weights 1/20, robust SE clustered on directed dyad, all six steps"
ahc <- "Y (proximity of past fatal actions actor->target and target->actor, proximity of independence)"
nt <- sprintf(paste0("REBUILT (no archive): calendar-time Cox as conditional risk sets at each fatal-action time 1951-92 (GML 2.2.1 directed dispute data; ",
  "fatal action = actor hostlev >= 4 in a dyadic dispute with fatality >= 1 for either state; target side at t + 0.5 day); 5%% sample of non-failures (seed 123), weight 20. ",
  "s, e, D-size from SDP GDP (2011 USD; paper PWT 5.6 current USD) and COW smoothed bilateral trade (paper Gleditsch 2002), lag 1; democ Polity IV (paper Polity IIId daily); ",
  "ln GDPpc SDP; contiguity COW conttype <= 3; ln capital distance. Failures %d (paper 438). Full rebuilt b=%.2f se=%.2f p=%.4g N=%d rows. Common 1822-2011 N=%d. ",
  "Key term e*s^2 is a product; constituent terms kept in all steps."),
  sum(d0$ev), full[[K]][["b"]], full[[K]][["se"]], full[[K]][["p"]], attr(full, "n"), nrow(ds))
write_p8(meta, K, "abstract: trade reduces conflict incentives 'most clearly seen in relatively symmetric dyads'",
  "Table 1 Model II (Cox, fatal directed militarized actions 1951-92; trade efficiency x size asymmetry squared)", art, used,
  310.5, 89.7, NA, TRUE, TRUE, M, ahc, nt)
dir.create(file.path(OUT, "results_p9"), recursive = TRUE, showWarnings = FALSE)
src <- file.path(OUT, "results_p8", "P9_033__e_s2.csv"); dst <- file.path(OUT, "results_p9", "P9_033__e_s2.csv")
stopifnot(file.exists(src), file.rename(src, dst))
sink()
