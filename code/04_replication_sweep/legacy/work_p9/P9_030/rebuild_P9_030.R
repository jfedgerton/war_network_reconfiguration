# P9_030 Bearce, Flanagan & Floros 2006 IO 60(3) "Alliances, Internal Information, and Military Conflict Among Member-States"
#   (10.1017/S0020818306060188). REBUILD (protocol B): no author code or archive found; staged article PDF + pdftotext only.
# Article: logit of MID involvement ("dispute", O&R 1997), politically relevant dyad-years 1951-85 (O&R 1997 sample), robust SE
#   clustered on dyad, peace years + three splines (BKT 1998). Table 2 adds power preponderance x joint alliance to Model 1.3.
# HEADLINE (fixed from the PDF before estimating): Table 2 Model 2.2, key term = product (ln) power preponderance x joint alliance
#   (H2). Published 0.30 (0.10), p < .002; joint alliance (H1) -0.88 (0.27); N = 18,252 (table note).
#   Why 2.2 over 2.1: the text adopts the ln version for "all the statistical models that follow" (Model 2.1 uses raw caprat - 1,
#   interaction 0.006 (0.002)); 2.3/2.5 use a generated alliance instrument with bootstrap SEs; 2.4/2.5 change the DV to use of force.
# Controls (Model 1.3 set): affinity s (Signorino-Ritter S on UNGA votes), democracy high, democracy low, economic growth (lower),
#   contiguity, trade interdependence (lower trade/GDP), ln power preponderance, joint alliance, peace years + 3 splines.
# Rebuild inputs (reused from X076 / P9_029): Russett-Oneal TRIANGLE PR dyads (dispute1, demaut, allies, lcaprat2, noncontg, depend);
#   Gartzke Affinity 2010 s3un4608i lagged 1 yr; lower GDPpc growth (mean t-1..t-3) from peacesciencer SDP.
#   Peace years + splines rebuilt with the btscs rule used in work_demcap/D1 and work/43 (py within dyad, knots = quartiles of py,
#   T = max py), computed on the 1950-85 frame (the O&R 1997 data start in 1950), then the sample is cut to 1951-85.
CODE <- Sys.getenv("SWEEP_CODE"); stopifnot(nzchar(CODE))
source(file.path(CODE, "work_p8/common_b3.R"))   # ROOT, OUT, dplyr, fixest, add_history, common_sample, run_six, write_p8, fit_bin; seed 123
suppressPackageStartupMessages({library(haven); library(peacesciencer)})
set.seed(123)
W <- file.path(OUT, "work_p9/P9_030"); dir.create(W, recursive = TRUE, showWarnings = FALSE)
sink(file.path(W, "log_P9_030.txt"), split = TRUE)

# ---- 1. Russett-Oneal TRIANGLE, PR dyads 1950-85 (1950 kept only to start the peace-year count) ----
t <- read.delim(file.path(ROOT, "recA/or/TRIANGLE.tab"))
stopifnot(nrow(t) == 39996)
t <- t[t$year >= 1950 & t$year <= 1985, ]
stopifnot(nrow(t) == 23589, all(t$statea < t$stateb))
t <- t %>% mutate(c1 = statea, c2 = stateb, dyadid = 1000 * c1 + c2,
  demhigh = pmax(demauta, demautb), demlow = pmin(demauta, demautb),
  contig = -noncontg,                                  # noncontg stored as (dummy - 1): 0 = noncontig, -1 = contig
  smldep = pmin(dependa, dependb),                     # trade interdependence = lower bilateral trade/GDP
  lnprep = lcaprat2,                                   # ln(stronger/weaker capability); 0 = parity
  lnprep_x_ally = lcaprat2 * allies)                   # KEY: (ln) power preponderance x joint alliance
stopifnot(all(t$contig %in% c(0, 1)), all(t$lnprep >= -1e-8, na.rm = TRUE), anyDuplicated(t[, c("c1", "c2", "year")]) == 0,
          all(t$dispute1 %in% c(0, 1)), all(t$allies %in% c(0, 1)))
n_t <- nrow(t)

# ---- 2. peace years + 3 cubic splines (btscs dispute1 year dyadid, nspline(3)) on the 1950-85 frame ----
t <- t[order(t$dyadid, t$year), ]
first <- c(TRUE, t$dyadid[-1] != t$dyadid[-nrow(t)]); ev <- t$dispute1; y <- t$year
py <- numeric(nrow(t))
for (i in seq_len(nrow(t))) py[i] <- if (first[i]) 0 else if (isTRUE(ev[i - 1] == 1)) 0 else py[i - 1] + (y[i] - y[i - 1])
t$py <- py
kn <- quantile(t$py, c(.25, .5, .75), type = 2); Tm <- max(t$py)
for (i in 1:3) t[[paste0("spl", i)]] <- pmax(t$py - kn[i], 0)^3 - (Tm - kn[i]) / Tm * t$py^3
cat("btscs knots:", kn, " T:", Tm, "\n")
stopifnot(Tm <= 35, all(diff(kn) >= 0))

# ---- 3. affinity S, lag 1 (affinity in year y -> dyad-year y+1) ----
a <- zap_labels(read_dta(file.path(ROOT, "work_p4/X076/affinity_01242010.dta"))) %>% as.data.frame()
stopifnot("s3un4608i" %in% names(a))
a <- a %>% transmute(c1 = pmin(ccodea, ccodeb), c2 = pmax(ccodea, ccodeb), year = year + 1, aff_l1 = s3un4608i) %>% filter(c1 != c2) %>%
  group_by(c1, c2, year) %>% summarise(aff_l1 = mean(aff_l1, na.rm = TRUE), .groups = "drop") %>% mutate(aff_l1 = ifelse(is.nan(aff_l1), NA, aff_l1))
t <- t %>% left_join(a, by = c("c1", "c2", "year"))
stopifnot(nrow(t) == n_t, all(is.na(t$aff_l1) | abs(t$aff_l1) <= 1))

# ---- 4. lower economic growth: mean annual % change in GDPpc over t-1..t-3 (SDP), weak link ----
sy <- create_stateyears(system = "cow", subset_years = 1940:1990) %>% add_sdp_gdp()
stopifnot("wbgdppc2011est" %in% names(sy), anyDuplicated(sy[, c("ccode", "year")]) == 0)
sy <- sy %>% arrange(ccode, year) %>% group_by(ccode) %>%
  mutate(gpc = exp(wbgdppc2011est), gr = 100 * (gpc / lag(gpc) - 1), gr3 = (lag(gr, 1) + lag(gr, 2) + lag(gr, 3)) / 3) %>%
  ungroup() %>% select(ccode, year, gr3) %>% as.data.frame()
t <- t %>% left_join(sy %>% rename(ga = gr3), by = c("statea" = "ccode", "year")) %>%
  left_join(sy %>% rename(gb = gr3), by = c("stateb" = "ccode", "year")) %>%
  mutate(lowgrowth = pmin(ga, gb)) %>% as.data.frame()
stopifnot(nrow(t) == n_t)

# ---- 5. estimation sample 1951-85 ----
t <- t[t$year >= 1951, ]
cat("TRIANGLE PR dyads 1951-85:", nrow(t), " dispute1 mean:", round(mean(t$dispute1), 4), " affinity non-missing:", sum(!is.na(t$aff_l1)),
    " lowgrowth non-missing:", sum(!is.na(t$lowgrowth)), "\n")

# ---- 6. headline model: Table 2 Model 2.2 ----
rhs0 <- c("aff_l1", "demhigh", "demlow", "lowgrowth", "contig", "smldep", "lnprep", "allies", "lnprep_x_ally", "py", "spl1", "spl2", "spl3")
vars <- c("dispute1", rhs0, "dyadid")
fitter <- function(r, data) fit_bin("dispute1", r, data, "logit", "~dyadid", terms = c("lnprep_x_ally", "allies", "lnprep"))
d0 <- t[complete.cases(t[, vars]), ]
stopifnot(nrow(d0) > 15000)
full <- fitter(rhs0, d0)
cat("\nREBUILD vs published (Table 2 Model 2.2): lnprep_x_ally b", full$lnprep_x_ally[["b"]], "se", full$lnprep_x_ally[["se"]],
    "p", full$lnprep_x_ally[["p"]], "N", attr(full, "n"), "| published 0.30 (0.10) N 18252\n")
cat("allies b", full$allies[["b"]], "se", full$allies[["se"]], "| published -0.88 (0.27);  lnprep b", full$lnprep[["b"]], "| published -0.30 (0.06)\n")

ds <- common_sample(t, vars, "statea", "stateb", "year")
M <- run_six(rhs0, ds, fitter)
saveRDS(list(M = M, full = full, n_used = nrow(ds), knots = kn, Tmax = Tm), file.path(W, "M_P9_030.rds"))
# Sensitivity hook: not wired. The key term is a product (lnprep x allies) whose constituents must stay in the model -> "split" class.

meta <- list(pid = "P9_030", cite = "Bearce, Flanagan & Floros 2006", journal = "IO", year = 2006, doi = "10.1017/S0020818306060188", family = "TIE")
art <- "logit, robust SE clustered on dyad, peace years + 3 splines"
used <- "feglm logit, vcov clustered dyadid, rebuilt btscs peace years + 3 splines, all six steps"
ahc <- "N (peace years + 3 splines only)"
nt <- sprintf(paste0("REBUILT (no archive): Russett-Oneal TRIANGLE PR dyads 1951-85 (dispute1, demaut high/low, allies, ln capability ratio, ",
  "contiguity, lower trade/GDP) + Gartzke Affinity 2010 s3un4608i lag 1 + lower GDPpc growth (SDP; paper: O&R 1997) + rebuilt btscs ",
  "(knots %s, T %g, counted from 1950). Full rebuilt N %d (published 18,252): b=%.5f se=%.5f; joint alliance %.4f (pub -0.88). ",
  "Common 1822-2011 sample N %d. Table 2 Model 2.2 (ln preponderance version adopted by the authors); interaction constituents kept."),
  paste(kn, collapse = "/"), Tm, attr(full, "n"), full$lnprep_x_ally[["b"]], full$lnprep_x_ally[["se"]], full$allies[["b"]], nrow(ds))
write_p8(meta, "lnprep_x_ally",
  "abstract: 'information provided within the alliance matters most for dyads at or near power parity'",
  "Table 2 Model 2.2 (Model 1.3 + ln power preponderance x joint alliance; PR dyads 1951-85)", art, used,
  0.30, 0.10, 18252, TRUE, TRUE, M, ahc, nt)
# phase-9 output location (shared write_p8 writes to results_p8)
dir.create(file.path(OUT, "results_p9"), recursive = TRUE, showWarnings = FALSE)
src <- file.path(OUT, "results_p8", "P9_030__lnprep_x_ally.csv"); dst <- file.path(OUT, "results_p9", "P9_030__lnprep_x_ally.csv")
stopifnot(file.exists(src), file.rename(src, dst))
sink()
