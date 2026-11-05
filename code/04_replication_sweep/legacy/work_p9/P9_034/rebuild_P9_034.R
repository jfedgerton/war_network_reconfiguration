# P9_034 Barbieri 1996 JPR 33(1) "Economic Interdependence: A Path to Peace or a Source of Interstate Conflict?"
#   (10.1177/0022343396033001003). REBUILD (protocol B; no author code or archive). Scanned PDF read page by page (pp. 35-39, notes p. 45).
# Design: non-directed dyad-years of COW states 1870-1938 excluding WWI, only dyads with import, export and total-trade data
#   (N = 14,341 of 70,293). DV (year t) = MID occurrence, states on opposite sides, coded only in the first year (onset).
#   All IVs at t-1. Trade share_i = dyadic trade_ij / total trade_i; Salience = sqrt(share_i * share_j); Symmetry = 1 - |share_i - share_j|;
#   Interdependence = Salience * Symmetry. Controls: contiguity (direct or indirect, land or sea <= 150 miles), joint democracy
#   (Polity II democracy >= 5 for both), alliance (any COW defense / neutrality / entente), relative capabilities (ln larger/smaller CINC).
#   Logit, conventional SEs. Table I, N = 14,341.
# HEADLINE / KEY (fixed from the PDF before estimating): Table I Model 3 "Interdependence", interdependence_{t-1} = 3.50 (SE 1.25), p <= .01.
#   Abstract: "extensive economic interdependence increases the likelihood that dyads will engage in militarized interstate disputes".
#   Model 3 is the single additive interdependence test; Model 4 (salience + symmetry + product) is the "full model".
# Rebuild inputs: COW trade v4 (peacesciencer extdata cow_trade_ddy.rds raw flows; cow_trade_sy state totals -- Barbieri's own trade
#   database is the pre-1938 basis of COW trade), COW MID via peacesciencer add_cow_mids, COW alliances, COW contiguity, NMC CINC,
#   Polity IV democ (p4v2014.csv on disk under work_p8/P8_049).
CODE <- Sys.getenv("SWEEP_CODE"); stopifnot(nzchar(CODE))
source(file.path(CODE, "work_p8/common_b3.R"))   # ROOT, OUT, dplyr, fixest, add_history, common_sample, run_six, write_p8, fit_bin; seed 123
suppressPackageStartupMessages(library(peacesciencer))
set.seed(123)
W <- file.path(OUT, "work_p9/P9_034"); dir.create(W, recursive = TRUE, showWarnings = FALSE)
sink(file.path(W, "log_P9_034.txt"), split = TRUE)

# ---- 1. non-directed COW dyad-years 1870-1938 with MID onset (DV at t) ----
dy <- create_dyadyears(directed = FALSE, subset_years = 1870:1938) %>% add_cow_mids() %>% as.data.frame()
n0 <- nrow(dy); stopifnot(n0 > 70000)
stopifnot(anyDuplicated(dy[, c("ccode1", "ccode2", "year")]) == 0, all(dy$ccode1 < dy$ccode2), all(dy$cowmidonset %in% 0:1))
dy <- dy[!(dy$year %in% 1914:1918), ]; n0 <- nrow(dy)
cat("population dyad-years 1870-1938 excl. 1914-18:", n0, "(paper 70,293)  MID onsets:", sum(dy$cowmidonset), "\n")

# ---- 2. IVs at t-1 ----
trp <- system.file("extdata", "cow_trade_ddy.rds", package = "peacesciencer"); stopifnot(nzchar(trp), file.exists(trp))
tr <- readRDS(trp) %>% as.data.frame() %>% filter(ccode1 < ccode2) %>% transmute(ccode1, ccode2, year = year + 1, dtrade = flow1 + flow2)
stopifnot(anyDuplicated(tr[, c("ccode1", "ccode2", "year")]) == 0)
sy <- as.data.frame(cow_trade_sy) %>% transmute(ccode, year = year + 1, tot = imports + exports)
stopifnot(anyDuplicated(sy[, c("ccode", "year")]) == 0)
dy <- dy %>% left_join(tr, by = c("ccode1", "ccode2", "year")) %>%
  left_join(sy %>% rename(tot1 = tot), by = c("ccode1" = "ccode", "year")) %>% left_join(sy %>% rename(tot2 = tot), by = c("ccode2" = "ccode", "year"))
stopifnot(nrow(dy) == n0)
dy <- dy %>% mutate(sh1 = ifelse(tot1 > 0, dtrade / tot1, NA_real_), sh2 = ifelse(tot2 > 0, dtrade / tot2, NA_real_),
  salience = sqrt(sh1 * sh2), symmetry = 1 - abs(sh1 - sh2), interdep = salience * symmetry)
cat("trade shares > 1 (data error, kept):", sum(dy$sh1 > 1 | dy$sh2 > 1, na.rm = TRUE), "\n")
# controls at t-1: build on t-1 dyad-years and shift
cy <- create_dyadyears(directed = FALSE, subset_years = 1869:1937) %>% add_contiguity() %>% add_nmc() %>% as.data.frame()
nc <- nrow(cy)
al <- as.data.frame(cow_alliance) %>% filter(ccode1 < ccode2) %>% select(ccode1, ccode2, year, cow_defense, cow_neutral, cow_entente)   # add_cow_alliance() is disabled
stopifnot(anyDuplicated(al[, c("ccode1", "ccode2", "year")]) == 0)
cy <- cy %>% left_join(al, by = c("ccode1", "ccode2", "year")) %>% mutate(across(c(cow_defense, cow_neutral, cow_entente), ~coalesce(.x, 0)))
stopifnot(nrow(cy) == nc)
stopifnot(anyDuplicated(cy[, c("ccode1", "ccode2", "year")]) == 0)
p4 <- read.csv(file.path(ROOT, "work_p8/P8_049/p4v2014.csv")) %>% transmute(ccode, year, democ = ifelse(democ %in% 0:10, democ, NA_real_)) %>%
  group_by(ccode, year) %>% summarise(democ = if (all(is.na(democ))) NA_real_ else mean(democ, na.rm = TRUE), .groups = "drop")   # transition years: mean (as P9_033)
cy <- cy %>% left_join(p4 %>% rename(dem1 = democ), by = c("ccode1" = "ccode", "year")) %>% left_join(p4 %>% rename(dem2 = democ), by = c("ccode2" = "ccode", "year"))
cy <- cy %>% transmute(ccode1, ccode2, year = year + 1,
  contig = as.numeric(conttype %in% 1:4),
  jointdem = ifelse(is.na(dem1) | is.na(dem2), NA_real_, as.numeric(dem1 >= 5 & dem2 >= 5)),
  alliance = as.numeric(cow_defense == 1 | cow_neutral == 1 | cow_entente == 1),
  relcap = ifelse(pmin(cinc1, cinc2) > 0, log(pmax(cinc1, cinc2) / pmin(cinc1, cinc2)), NA_real_))
dy <- dy %>% left_join(cy, by = c("ccode1", "ccode2", "year")) %>% mutate(dyadid = 1000 * ccode1 + ccode2) %>% as.data.frame()
stopifnot(nrow(dy) == n0)
cat("missing by var:\n"); print(sapply(dy[, c("dtrade", "tot1", "tot2", "interdep", "contig", "jointdem", "alliance", "relcap")], function(x) sum(is.na(x))))

# ---- 3. headline fit: Table I Model 3, logit, conventional SE ----
rhs0 <- c("interdep", "contig", "jointdem", "alliance", "relcap")
vars <- c("cowmidonset", rhs0, "dyadid")
d0 <- dy[complete.cases(dy[, vars]), ]
cat("complete-case N:", nrow(d0), "(paper 14,341)  MID onsets:", sum(d0$cowmidonset), "(paper 270)  means: salience", round(mean(d0$salience), 3),
    "symmetry", round(mean(d0$symmetry), 3), "interdep", round(mean(d0$interdep), 3), "(paper .03 .9 .03)  MID rate", round(mean(d0$cowmidonset), 4), "(paper ~.019)\n")
fitter <- function(r, data) fit_bin("cowmidonset", r, data, "logit", "iid", terms = "interdep")
full <- fitter(rhs0, d0)
cat("\nREBUILT vs published (Table I Model 3): interdep b", full$interdep[["b"]], "se", full$interdep[["se"]], "p", full$interdep[["p"]],
    "N", attr(full, "n"), "| published 3.50 (1.25) N 14,341\n")
ref_cl <- fit_bin("cowmidonset", rhs0, d0, "logit", "~dyadid", terms = "interdep")   # reference only
cat("reference dyad-clustered SE:", ref_cl$interdep[["se"]], "\n")
ds <- common_sample(dy, vars, "ccode1", "ccode2", "year")
M <- run_six(rhs0, ds, fitter)
saveRDS(list(M = M, full = full, n_used = nrow(ds)), file.path(W, "M_P9_034.rds"))

# ---- sensitivity hook: runs only under run_sensitivity.R (SWEEP_SENS=1), never in the main sweep ----
if (Sys.getenv("SWEEP_SENS") == "1") {
  source(file.path(CODE, "hist/sens_helper.R"))
  sens_run(id = "P9_034", dat = ds, y = "cowmidonset", tie = "interdep", ctrl = setdiff(rhs0, "interdep"), fe = NULL,
           fit = function(dat, add) fitter(c(rhs0, add), dat)[["interdep"]])
}

meta <- list(pid = "P9_034", cite = "Barbieri 1996", journal = "JPR", year = 1996, doi = "10.1177/0022343396033001003", family = "TIE")
art <- "logit, conventional SE, all IVs lagged 1 yr, no temporal controls (Table I)"
used <- "feglm logit, iid SE, all six steps"
ahc <- "N (no temporal/history controls)"
nt <- sprintf(paste0("REBUILT (no archive): COW non-directed dyad-years 1870-1938 excl. 1914-18 with non-missing t-1 COW trade (raw flows, zeros kept) and state totals ",
  "(cow_trade_sy); DV = COW MID onset (peacesciencer add_cow_mids, opposite sides); contiguity conttype 1-4 (direct only, paper also indirect via colonies); ",
  "joint democracy Polity IV democ >= 5 (paper Polity II); any COW alliance; ln CINC ratio. Full rebuilt b=%.4f se=%.4f p=%.4g N=%d (published 3.50 (1.25) N 14,341). ",
  "Dyad-clustered SE %.4f for reference. Common 1822-2011 N=%d."),
  full$interdep[["b"]], full$interdep[["se"]], full$interdep[["p"]], attr(full, "n"), ref_cl$interdep[["se"]], nrow(ds))
write_p8(meta, "interdep", "abstract: 'extensive economic interdependence increases the likelihood that dyads will engage in militarized interstate disputes'",
  "Table I Model 3 (Interdependence = salience x symmetry, t-1; MID occurrence 1870-1938)", art, used,
  3.50, 1.25, 14341, TRUE, TRUE, M, ahc, nt)
dir.create(file.path(OUT, "results_p9"), recursive = TRUE, showWarnings = FALSE)
src <- file.path(OUT, "results_p8", "P9_034__interdep.csv"); dst <- file.path(OUT, "results_p9", "P9_034__interdep.csv")
stopifnot(file.exists(src), file.rename(src, dst))
sink()
