# P9_093 Gartzke & Gleditsch 2022, ISQ 66(1) sqab082 "Ties That Bias in International Conflict: A Spatial Approach to Dyadic Dependence
#   from Alliance Ties and Inbetweenness" (10.1093/isq/sqab082). TIE (third-party alliance ties). REBUILD (protocol B; no Dataverse deposit found).
# Article (staged PDF pp. 6-9, Tables 1-2): Oneal-Russett (1999) all dyads 1950-1992, N 271,262; DV = Maoz dyadic MID onset, subsequent
#   dispute years dropped; logit with O-R covariates (joint democracy = lower Polity score, lower/higher dependence, contiguity, logged
#   capital distance, major power, alliance, logged capability ratio, BKT peace years) plus extra-dyadic terms lagged one year relative to
#   the DV: third-order alliance tie (dyad pairs a disputant with an ally of its opponent), fourth-order alliance tie (allies of A paired
#   with allies of B, neither a disputant), MID inbetweenness ratio (third party C to dispute AB: D(AB)/(D(AC)+D(BC)), max over ongoing
#   disputes, 0 for the disputants), Max(A,B IBR) of the two members to disputes not involving them, and its product with logged
#   intercapital distance. SEs are window-subsampling (WSEV.1) estimates.
# HEADLINE (fixed before estimating): Table 2, Third-order alliance tie 0.789 (WSEV.1 SE 0.189), N 271,262.
#   No star legend; the text uses "standard significance tests" -> alpha_pub 0.05 (low confidence, coded from the default).
#   Verdict rule (PI, round 2): published sign, significance at alpha_pub, rebuilt b0 within 0.5x-2x.
# Rebuild inputs and DISCLOSED DEPARTURES:
#   - O-R all-dyad covariates and Maoz onsets from work_p4/X074/or_trade.rds + work_demcap/D1/dor.dta (or_dyads): mzmid1_all = Maoz
#     dyadic MID onset in t+1 including joiners (Maoz dyads are not restricted to originators; the article stresses this).
#   - Ongoing dispute years: GML dyadic disputes ({peacesciencer} gml_dirdisp), used (a) to drop dyad-years whose t+1 is a continuing
#     dispute year, and (b) as the set of disputes ongoing in year t that activates the extra-dyadic terms (Maoz DYDMID 2.0 not on disk).
#   - Alliances for third/fourth-party ties: COW formal alliance v4.1 dyad-year file (any type), as in O-R.
#   - Capital distances: {peacesciencer} capital distance (Gleditsch-Ward style, km) for the inbetweenness ratios.
#   - SE: logit clustered on dyad (the pipeline standard); WSEV.1 is not implemented (novel SE estimator) -> SE comparison not exact.
ROOT <- Sys.getenv("SWEEP_DATA"); OUT <- Sys.getenv("SWEEP_OUT"); CODE <- Sys.getenv("SWEEP_CODE")
stopifnot(nzchar(ROOT), nzchar(OUT), nzchar(CODE))
source(file.path(CODE, "work_p8/common_b3.R")); source(file.path(CODE, "work_p9/rebuild_common_p9r.R"))
set.seed(123)
W <- file.path(OUT, "work_p9/P9_093"); dir.create(W, recursive = TRUE, showWarnings = FALSE)
sink(file.path(W, "log_P9_093.txt"), split = TRUE)
YRS <- 1950:1992
# ---- 1. O-R dyads with Maoz onsets ----
d <- or_dyads(YRS) %>% filter(year %in% YRS) %>% as.data.frame()
stopifnot(all(d$statea < d$stateb), !anyDuplicated(d[, c("statea", "stateb", "year")]))
d <- d %>% mutate(highdep = pmax(dep_a, dep_b))
cat("O-R dyad-years 1950-92:", nrow(d), "(published N 271,262)\n")
# ---- 2. GML dyadic disputes: ongoing sets by year, continuing years ----
g <- as.data.frame(gml_dirdisp) %>% filter(ccode1 < ccode2) %>% transmute(dispnum, a = ccode1, b = ccode2, year, onset = midonset)
stopifnot(nrow(g) > 0)
on_y <- g %>% distinct(a, b, year)                                             # dyads in an ongoing dispute in year
cont <- g %>% group_by(a, b, year) %>% summarise(newd = max(onset), .groups = "drop") %>% filter(newd == 0) %>% transmute(statea = a, stateb = b, year = year - 1, cont1 = 1L)
d <- d %>% left_join(cont, by = c("statea", "stateb", "year")) %>% mutate(cont1 = coalesce(cont1, 0L))
cat("dyad-years dropped as continuing dispute years in t+1 (no new onset):", sum(d$cont1 == 1 & d$mzmid1_all %in% 0), "\n")
d <- d[!(d$cont1 == 1 & d$mzmid1_all %in% 0), ]
# ---- 3. alliances (COW v4.1, any type) ----
al <- read.csv(file.path(ROOT, "work_p8/P8_049/alliance_v4_1_by_dyad_yearly.csv")) %>% transmute(x = pmin(ccode1, ccode2), y = pmax(ccode1, ccode2), year) %>%
  filter(year %in% YRS, x != y) %>% distinct()
# ---- 4. capital distances (km) ----
cd <- suppressMessages(create_dyadyears(system = "cow", directed = FALSE, subset_years = YRS) %>% add_capital_distance()) %>% as.data.frame() %>%
  transmute(x = ccode1, y = ccode2, year, km = capdist)
stopifnot(!anyDuplicated(cd[, c("x", "y", "year")]))
# ---- 5. extra-dyadic terms, year by year ----
ext <- list()
for (yr in YRS) {
  A <- on_y[on_y$year == yr, ]; al_y <- al[al$year == yr, ]; c_y <- cd[cd$year == yr & !is.na(cd$km), ]
  st <- sort(unique(c(c_y$x, c_y$y)))
  Dm <- matrix(NA_real_, length(st), length(st), dimnames = list(st, st)); diag(Dm) <- 0
  Dm[cbind(as.character(c_y$x), as.character(c_y$y))] <- c_y$km; Dm[cbind(as.character(c_y$y), as.character(c_y$x))] <- c_y$km
  allies <- function(s) unique(c(al_y$y[al_y$x == s], al_y$x[al_y$y == s]))
  t3 <- list(); t4 <- list(); ib3 <- list(); ib4 <- setNames(numeric(length(st)), st)
  for (k in seq_len(nrow(A))) {
    a <- A$a[k]; b <- A$b[k]; aa <- setdiff(allies(a), c(a, b)); ab <- setdiff(allies(b), c(a, b))
    if (length(ab)) t3[[length(t3) + 1]] <- data.frame(p = a, q = ab)
    if (length(aa)) t3[[length(t3) + 1]] <- data.frame(p = b, q = aa)
    if (length(aa) && length(ab)) { e <- expand.grid(p = aa, q = ab); e <- e[e$p != e$q, ]; if (nrow(e)) t4[[length(t4) + 1]] <- e }
    ca <- as.character(a); cb <- as.character(b)
    if (ca %in% st && cb %in% st && is.finite(Dm[ca, cb])) {
      r <- Dm[ca, cb] / (Dm[ca, ] + Dm[cb, ]); r[c(ca, cb)] <- NA          # inbetweenness of every third state C to AB
      ok <- is.finite(r)
      ib3[[length(ib3) + 1]] <- rbind(data.frame(p = a, q = st[ok], v = r[ok]), data.frame(p = b, q = st[ok], v = r[ok]))
      ib4[ok] <- pmax(ib4[ok], r[ok])                                         # each state's max IBR to disputes it is not part of
    }
  }
  und <- function(z) if (length(z)) { z <- do.call(rbind, z); data.frame(statea = pmin(z$p, z$q), stateb = pmax(z$p, z$q), v = if (is.null(z$v)) 1 else z$v) } else NULL
  T3 <- und(t3); T4 <- und(t4); I3 <- und(ib3)
  base <- d[d$year == yr, c("statea", "stateb", "year")]
  base$tie3 <- 0L; base$tie4 <- 0L; base$ibr3 <- 0
  key <- paste(base$statea, base$stateb)
  if (!is.null(T3)) base$tie3[key %in% paste(T3$statea, T3$stateb)] <- 1L
  if (!is.null(T4)) base$tie4[key %in% paste(T4$statea, T4$stateb)] <- 1L
  if (!is.null(I3)) { m <- tapply(I3$v, paste(I3$statea, I3$stateb), max); hit <- key %in% names(m); base$ibr3[hit] <- m[key[hit]] }
  # the disputants themselves get 0 on the inbetweenness ratio
  disp <- paste(A$a, A$b); base$ibr3[key %in% disp] <- 0
  base$ibr4 <- pmax(ib4[as.character(base$statea)], ib4[as.character(base$stateb)]); base$ibr4[is.na(base$ibr4)] <- 0
  ext[[as.character(yr)]] <- base
}
ext <- do.call(rbind, ext); n0 <- nrow(d)
d <- d %>% left_join(ext, by = c("statea", "stateb", "year")); stopifnot(nrow(d) == n0)
d <- d %>% mutate(ibr4_x_ldist = ibr4 * ldist, outcome_year = year + 1)
cat("extra-dyadic terms: share tie3", mean(d$tie3), " tie4", mean(d$tie4), " mean ibr3", mean(d$ibr3), " mean ibr4", mean(d$ibr4), "\n")
# ---- 6. peace years and splines on the Maoz onset ----
d <- py_calc(d, "dyadid", "mzmid1_all"); d <- btscs_splines(d)
rhs0 <- c("tie3", "tie4", "ibr3", "ibr4", "ibr4_x_ldist", "lowdem", "lowdep", "highdep", "contig150", "ldist", "majpow", "anyally", "lcaprat",
          "py", "spl1", "spl2", "spl3")
vars <- c("mzmid1_all", rhs0, "dyadid")
cat("missing:\n"); print(sapply(d[, vars], function(x) sum(is.na(x))))
d0 <- d[complete.cases(d[, vars]), ]
cat("complete-case N:", nrow(d0), "(published 271,262)  onsets", sum(d0$mzmid1_all), "\n")
TERM <- "tie3"
fitter <- function(r, data) fit_bin("mzmid1_all", r, data, "logit", "~dyadid", terms = TERM)
full <- fitter(rhs0, d0)
cat("\nREBUILT vs published: tie3 b", full[[TERM]][["b"]], "(pub 0.789) se", full[[TERM]][["se"]], "(pub WSEV.1 0.189) N", attr(full, "n"), "(pub 271,262)\n")
ds <- common_sample(d, vars, "statea", "stateb", "outcome_year"); rm(d, d0); invisible(gc())
M <- run_six(rhs0, ds, fitter)
if (Sys.getenv("SWEEP_SENS") == "1") {
  source(file.path(CODE, "hist/sens_helper.R"))
  sens_run(id = "P9_093", dat = ds, y = "mzmid1_all", tie = TERM, ctrl = setdiff(rhs0, TERM), fe = NULL,
           fit = function(dat, add) fitter(c(rhs0, add), dat)[[TERM]])
}
saveRDS(list(M = M, full = full), file.path(W, "M_P9_093.rds"))
meta <- list(pid = "P9_093", cite = "Gartzke & Gleditsch 2022", journal = "ISQ", year = 2022, doi = "10.1093/isq/sqab082", family = "TIE")
nt <- sprintf(paste0("REBUILT (no archive): O-R all dyads 1950-92 (or_trade/DOR), DV mzmid1_all (Maoz dyadic onset t+1 incl. joiners), continuing dispute ",
  "years dropped (GML dyadic); third/fourth-party alliance ties from GML disputes ongoing in t and COW v4.1 alliances; inbetweenness ratios from ",
  "peacesciencer capital distances; O-R covariates; BTSCS peace years + 3 splines; logit clustered on dyad (WSEV.1 not implemented). Full-sample ",
  "rebuilt b=%.6f se=%.6f N=%d vs published 0.789 (0.189) N 271,262. Common sample N=%d (history at outcome year)."),
  full[[TERM]][["b"]], full[[TERM]][["se"]], attr(full, "n"), nrow(ds))
b0 <- M$m0[[TERM]][["b"]]; p0 <- M$m0[[TERM]][["p"]]; r_ <- b0 / 0.789
v <- rebuild_verdict("P9_093", TERM, b0, p0, 0.789, NA, 0.189)   # round-2 rule, alpha_pub from registry.csv
cat("round-2 verdict:", v, " ratio", round(r_, 3), " p0", signif(p0, 3), "\n")
write_p8(meta, TERM, "abstract and Table 2: third- and fourth-party alliance ties raise the risk of dispute onset",
  "Table 2 (logit, Maoz dyadic MID onset, all dyads 1950-92; third-order alliance tie)", "logit with WSEV.1 window-subsampling SEs, BKT peace years",
  "feglm logit, SE clustered on dyad, all six steps", 0.789, 0.189, 271262, TRUE, TRUE, M, "Y (peace years + splines)", nt, verdict = v)
dir.create(file.path(OUT, "results_p9"), recursive = TRUE, showWarnings = FALSE)
src <- file.path(OUT, "results_p8", "P9_093__tie3.csv"); stopifnot(file.exists(src), file.rename(src, file.path(OUT, "results_p9", basename(src))))
sink()
