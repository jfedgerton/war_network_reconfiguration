# P9_036 Norrevik & Sarwari 2025 CMPS 42(3) "If we cooperate together, we intervene together: Defense cooperation agreements
#   and support to conflict parties" (10.1177/07388942241247953). REBUILD (protocol B): no author archive (Dataverse search
#   2026-10-10: none; SAGE supplement 403 / appendix PDF only; no data-availability statement in the article).
# HEADLINE: Table 1 Model 2. Logit, SE clustered on dyad. DV = both states of an undirected dyad give military support (troops
#   and/or arms) to the same civil-conflict warring party in year t (UCDP ESD). KEY TERM dca (DCAD dcaAnyV1, fn 11).
#   Controls: joint democracy (both Polity >= 6), joint NATO, economic interdependence (any COW trade flow), P5 member (either),
#   divisiveness (UNGA ideal points, Bailey et al.), same region, joint rivalry (MID with the same third state), previous DCA (t-1).
#   Published: dca 0.5196 (0.155) p < .001, N = 320,015. Model 1 (bivariate) N = 497,249 = every undirected DCAD dyad-year 1980-2010.
# Rebuild inputs (DATA/work_p9/P9_036): DCAD v1.0 dyadic (Dataverse doi:10.7910/DVN/DT3VTE), UCDP ESD 18.1 triad-year,
#   UNGA ideal points Jun2024 (Dataverse doi:10.7910/DVN/LEJUQZ); peacesciencer ccode_democracy, cow_mid_dirdisps, cow_trade_ddy.rds.
# Estimator: fixest feglm logit, vcov clustered on dyadid, all six steps. No randomness (seed 123 set in common_b3.R).
CODE <- Sys.getenv("SWEEP_CODE"); stopifnot(nzchar(CODE))
source(file.path(CODE, "work_p8/common_b3.R"))
suppressPackageStartupMessages({library(peacesciencer); library(countrycode)})
WD <- file.path(ROOT, "work_p9/P9_036"); W <- file.path(OUT, "work_p9/P9_036"); dir.create(W, recursive = TRUE, showWarnings = FALSE)
sink(file.path(W, "log_P9_036.txt"), split = TRUE)

# ---- 1. frame: undirected DCAD dyad-years 1980-2010 (reproduces Model 1 N = 497,249 and the 7,463 DCA dyad-years) ----
# Window choice (2026-10-10): the article's text says 1980-2009, but the 1980-2010 frame reproduces the paper's Model 1 N and DCA
# count exactly and its Model 2 N far better (rebuilt 321,603 vs published 320,015; a 1980-2009 frame gives 308,651), so the
# authors evidently kept 2010. Frame and ESD support events therefore use 1980-2010.
x <- as.data.frame(data.table::fread(file.path(WD, "DCAD-v1.0-dyadic.csv")))
stopifnot(nrow(x) == 994498, !anyDuplicated(x[, c("ccode1", "ccode2", "year")]))
d <- x[x$ccode1 < x$ccode2, c("ccode1", "ccode2", "year", "dcaAnyV1")]
names(d)[4] <- "dca"
stopifnot(nrow(d) == 497249, sum(d$dca) == 7463, !anyNA(d$dca), all(d$year >= 1980 & d$year <= 2010))
d$dyadid <- d$ccode1 * 1000 + d$ccode2
lg <- data.frame(dyadid = d$dyadid, year = d$year + 1, dca_l1 = d$dca)
n0 <- nrow(d); d <- left_join(d, lg, by = c("dyadid", "year")); stopifnot(nrow(d) == n0)
cat("frame rows:", nrow(d), " DCA=1:", sum(d$dca), " dca_l1 NA:", sum(is.na(d$dca_l1)), "\n")

# ---- 2. DV from UCDP ESD 18.1 triad-year: state supporters, civil conflicts, active years, troops (ext_x) or weapons (ext_w) ----
e <- as.data.frame(haven::zap_labels(haven::read_dta(file.path(WD, "ucdp-esd-ty-181.dta"))))
stopifnot(nrow(e) == 10852)
e <- e[e$ext_nonstate %in% 0 & e$civil == 1 & e$active == 1 & (e$ext_x %in% 1 | e$ext_w %in% 1) & e$year >= 1980 & e$year <= 2010, ]
cat("ESD supporter-recipient-years kept:", nrow(e), " alleged:", sum(e$ext_alleged == 1), "\n")
nm <- sub("^Government of ", "", e$ext_name)
e$sup <- countrycode(nm, "country.name", "cown", warn = FALSE)
stopifnot(!anyNA(e$sup))
e$sup[e$sup == 255 & e$year <= 1990] <- 260          # West Germany in COW codes through 1990 (DCAD uses 260 in 1990)
e$sup[e$sup == 678 & e$year >= 1990] <- 679          # unified Yemen from 1990
sy <- unique(rbind(data.frame(cc = d$ccode1, year = d$year), data.frame(cc = d$ccode2, year = d$year)))
miss <- anti_join(unique(e[, c("sup", "year")]), sy, by = c("sup" = "cc", "year"))
cat("supporter-years not in DCAD system:\n"); print(miss)
stopifnot(nrow(miss) == 0)
s <- unique(e[, c("sup", "actor_id", "year")])
pp <- inner_join(s, s, by = c("actor_id", "year"), relationship = "many-to-many")
pp <- unique(pp[pp$sup.x < pp$sup.y, c("sup.x", "sup.y", "year")])
names(pp) <- c("ccode1", "ccode2", "year"); pp$milsup <- 1
stopifnot(nrow(anti_join(pp, d, by = c("ccode1", "ccode2", "year"))) == 0)
n0 <- nrow(d); d <- left_join(d, pp, by = c("ccode1", "ccode2", "year")); stopifnot(nrow(d) == n0)
d$milsup[is.na(d$milsup)] <- 0
cat("DV=1 dyad-years:", sum(d$milsup), "(paper reports 7,969)\n")
tb <- table(dca = d$dca, milsup = d$milsup); print(tb)
cat("bivariate log odds ratio (data check vs Model 1 b 3.094):", log(tb[2, 2] * tb[1, 1] / (tb[2, 1] * tb[1, 2])), "\n")

# ---- 3. controls ----
pol <- ccode_democracy[, c("ccode", "year", "polity2")]
stopifnot(!anyDuplicated(pol[, c("ccode", "year")]))
n0 <- nrow(d)
d <- d %>% left_join(rename(pol, p1 = polity2), by = c("ccode1" = "ccode", "year")) %>%
  left_join(rename(pol, p2 = polity2), by = c("ccode2" = "ccode", "year"))
stopifnot(nrow(d) == n0)
d$jointdem <- as.numeric(d$p1 >= 6 & d$p2 >= 6)
# NATO membership by accession year (hard-coded; 260 West Germany 1955-90, 255 Germany from 1990)
nato <- rbind(data.frame(cc = c(2, 20, 200, 220, 210, 211, 212, 385, 390, 395, 235, 325), from = 1949),
              data.frame(cc = c(350, 640), from = 1952), data.frame(cc = c(260, 255), from = 1955),
              data.frame(cc = 230, from = 1982), data.frame(cc = c(316, 310, 290), from = 1999),
              data.frame(cc = c(355, 366, 367, 368, 360, 317, 349), from = 2004), data.frame(cc = c(339, 344), from = 2009))
stopifnot(!anyDuplicated(nato$cc))
in_nato <- function(cc, yr) { f <- nato$from[match(cc, nato$cc)]; !is.na(f) & yr >= f }
d$nato <- as.numeric(in_nato(d$ccode1, d$year) & in_nato(d$ccode2, d$year))
trp <- system.file("extdata", "cow_trade_ddy.rds", package = "peacesciencer"); stopifnot(nzchar(trp), file.exists(trp))
tr <- readRDS(trp)
tr <- as.data.frame(tr[tr$ccode1 < tr$ccode2, c("ccode1", "ccode2", "year", "flow1", "flow2")])
stopifnot(!anyDuplicated(tr[, c("ccode1", "ccode2", "year")]))
d <- left_join(d, tr, by = c("ccode1", "ccode2", "year")); stopifnot(nrow(d) == n0)
d$trade <- ifelse(is.na(d$flow1) & is.na(d$flow2), NA_real_, as.numeric(pmax(d$flow1, d$flow2, na.rm = TRUE) > 0))
p5 <- c(2, 200, 220, 365, 710)
d$p5 <- as.numeric(d$ccode1 %in% p5 | d$ccode2 %in% p5)
ip <- read.csv(file.path(WD, "IdealpointestimatesAll_Jun2024.csv"))
ip <- data.frame(ccode = ip$ccode, year = ip$session + 1945, ideal = ip$IdealPointAll)
stopifnot(!anyDuplicated(ip[, c("ccode", "year")]))
d <- d %>% left_join(rename(ip, i1 = ideal), by = c("ccode1" = "ccode", "year")) %>%
  left_join(rename(ip, i2 = ideal), by = c("ccode2" = "ccode", "year"))
stopifnot(nrow(d) == n0)
d$divis <- abs(d$i1 - d$i2)
reg <- function(cc) cut(cc, c(0, 199, 399, 626, 698, 999), labels = c("Amer", "Eur", "Afr", "ME", "Asia"))
stopifnot(!anyNA(reg(d$ccode1)), !anyNA(reg(d$ccode2)))
d$region <- as.numeric(reg(d$ccode1) == reg(d$ccode2))
# joint rivalry: both states in an ongoing MID (COW MID 5.0, opposite sides) with the same third state in year t
m <- as.data.frame(cow_mid_dirdisps)
m <- unique(m[m$sidea1 != m$sidea2 & m$year >= 1980 & m$year <= 2010, c("ccode1", "ccode2", "year")])
names(m) <- c("st", "opp", "year")
jr <- inner_join(m, m, by = c("opp", "year"), relationship = "many-to-many")
jr <- jr[jr$st.x < jr$st.y & jr$opp != jr$st.x & jr$opp != jr$st.y, ]
jr <- unique(data.frame(ccode1 = jr$st.x, ccode2 = jr$st.y, year = jr$year, jrival = 1))
d <- left_join(d, jr, by = c("ccode1", "ccode2", "year")); stopifnot(nrow(d) == n0)
d$jrival[is.na(d$jrival)] <- 0
cat("MID years in data: up to", max(cow_mid_dirdisps$year), "\n")
cat("means: jointdem", mean(d$jointdem, na.rm = TRUE), "nato", mean(d$nato), "trade", mean(d$trade, na.rm = TRUE), "p5", mean(d$p5),
    "divis", mean(d$divis, na.rm = TRUE), "region", mean(d$region), "jrival", mean(d$jrival), "\n")
cat("NA counts: jointdem", sum(is.na(d$jointdem)), "trade", sum(is.na(d$trade)), "divis", sum(is.na(d$divis)), "dca_l1", sum(is.na(d$dca_l1)), "\n")

# ---- 4. Table 1 Model 2, then the six steps on the common sample ----
term <- "dca"
rhs0 <- c(term, "jointdem", "nato", "trade", "p5", "divis", "region", "jrival", "dca_l1")
vars <- c("milsup", rhs0)
fitter <- function(r, data) fit_bin("milsup", r, data, "logit", "~dyadid", terms = term)
d0 <- d[complete.cases(d[, vars]), ]
full <- fitter(rhs0, d0); n_orig <- attr(full, "n")
cat("\nREBUILT Model 2:", term, "b", full[[term]][["b"]], "se", full[[term]][["se"]], "N", n_orig,
    "| published b 0.5196 se 0.155 N 320,015 | DV=1 in sample:", sum(d0$milsup), "\n")
ds <- common_sample(d, vars, "ccode1", "ccode2", "year")
M <- run_six(rhs0, ds, fitter)
saveRDS(list(M = M, full = full, n_orig = n_orig, n_used = nrow(ds)), file.path(W, "M_P9_036.rds"))

# ---- sensitivity hook: runs only under run_sensitivity.R (SWEEP_SENS=1), never in the main sweep ----
if (Sys.getenv("SWEEP_SENS") == "1") {
  source(file.path(CODE, "hist/sens_helper.R"))
  # dca enters the single-equation logit additively -> full sensitivity
  sens_run(id = "P9_036", dat = ds, y = "milsup", tie = term, ctrl = setdiff(rhs0, term), fe = NULL,
           fit = function(dat, add) fitter(c(rhs0, add), dat)[[term]])
}

meta <- list(pid = "P9_036", cite = "Norrevik & Sarwari 2025", journal = "CMPS", year = 2025, doi = "10.1177/07388942241247953", family = "TIE")
write_p8(meta, term,
  "abstract: 'DCA partners are more likely to provide military assistance to the same warring parties'",
  "Table 1 Model 2 (logit, joint military support to the same civil-conflict party, dyad-years 1980-2010)",
  "logit, SE clustered on dyad",
  "feglm logit, vcov clustered dyadid, all six steps", 0.5196, 0.155, 320015, TRUE, TRUE, M,
  "N (previous-year DCA; joint rivalry = MID with common third state in t)",
  sprintf(paste0("REBUILT (protocol B, no author data): DCAD v1.0 undirected dyad-years 1980-2010 (frame N and DCA count match the paper exactly), ",
    "DV from UCDP ESD 18.1 (paper: ESD 2011), ideal points Jun2024 release, polity2/MID/trade via peacesciencer, NATO and region hard-coded. ",
    "Rebuilt full-sample b=%.4f se=%.4f N=%d (published 0.5196, 0.155, 320,015). Common 1822-2011 sample N=%d."),
    full[[term]][["b"]], full[[term]][["se"]], n_orig, nrow(ds)))
dir.create(file.path(OUT, "results_p9"), recursive = TRUE, showWarnings = FALSE)
src <- file.path(OUT, "results_p8", paste0("P9_036__", term, ".csv")); dst <- file.path(OUT, "results_p9", paste0("P9_036__", term, ".csv"))
stopifnot(file.exists(src), file.rename(src, dst))
sink()
