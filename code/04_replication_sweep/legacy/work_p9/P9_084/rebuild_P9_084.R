# P9_084 Park 2013, CMPS 30(2) "Forward to the Future: The Democratic Peace after the Cold War" (10.1177/0738894212473927). DEMOCRACY.
#   REBUILD (protocol B; no archive).
# Article (staged PDF pp. 180-190, Table 2): post-Cold War equation 1990-2000, Gowa (2011) design (logit, clustered SEs, all dyads, MID onset,
#   contiguity, log distance, capability ratio, one/two major powers, duration controls) plus Cold War common-interest controls built from COW
#   formal alliances (v3.0): Cold War Rival (unallied dyad, one state allied with the US but not the USSR and the other the reverse),
#   US-Bloc Ally, Soviet-Bloc Ally, Dyadic Ally. Democracy: Lower and Higher Polity2 (weak link; the article states the continuous index is
#   preferable to the trichotomy, so Model 2 is the headline column).
# HEADLINE (fixed before estimating): Table 2 panel B (with common-interest controls), Model 2: Lower Democracy -0.0280 (0.0147) *, N 144,861.
# Rebuild inputs: {peacesciencer} COW dyad-years 1946-2000 (GML MID onset, contiguity 1-2, capital distance, CINC, majors); Polity IV 2014;
#   COW alliance v4.1 dyad-years (work_p8/P8_049). Bloc membership = COW alliance with the US (2) / USSR-Russia (365) in 1989 (the article
#   dates the end of the Cold War to 1989), applied to 1990-2000; Dyadic Ally = current COW alliance not covered by the bloc categories.
ROOT <- Sys.getenv("SWEEP_DATA"); OUT <- Sys.getenv("SWEEP_OUT"); CODE <- Sys.getenv("SWEEP_CODE")
stopifnot(nzchar(ROOT), nzchar(OUT), nzchar(CODE))
source(file.path(CODE, "work_p8/common_b3.R")); source(file.path(CODE, "work_p9/rebuild_common_p9r.R"))
set.seed(123)
W <- file.path(OUT, "work_p9/P9_084"); dir.create(W, recursive = TRUE, showWarnings = FALSE)
sink(file.path(W, "log_P9_084.txt"), split = TRUE)
d <- ps_dyads(1946:2000, trade = FALSE) %>% add_cow_ally()
d <- d %>% mutate(dyadid = 1000 * ccode1 + ccode2, mid = as.integer(gmlmidonset == 1))
d <- py_calc(d, "dyadid", "mid")
b89 <- d %>% filter(year == 1989, cow_any == 1, ccode1 %in% c(2, 365) | ccode2 %in% c(2, 365)) %>%
  transmute(sp = ifelse(ccode1 %in% c(2, 365), ccode1, ccode2), st = ifelse(ccode1 %in% c(2, 365), ccode2, ccode1))
us <- unique(c(2, b89$st[b89$sp == 2])); su <- unique(c(365, b89$st[b89$sp == 365]))
usonly <- setdiff(us, su); suonly <- setdiff(su, us)
cat("1989 US bloc:", length(us), " Soviet bloc:", length(su), "\n")
d <- d[d$year >= 1990, ]
p <- polity_cow("polity2"); n0 <- nrow(d)
d <- d %>% left_join(p %>% rename(p2a = polity2), by = c("ccode1" = "ccode", "year")) %>% left_join(p %>% rename(p2b = polity2), by = c("ccode2" = "ccode", "year"))
stopifnot(nrow(d) == n0)
d <- d %>% mutate(lowdem = pmin(p2a, p2b), highdem = pmax(p2a, p2b),
  cwrival = as.integer(cow_any == 0 & ((ccode1 %in% usonly & ccode2 %in% suonly) | (ccode1 %in% suonly & ccode2 %in% usonly))),
  usbloc = as.integer(cow_any == 1 & ccode1 %in% us & ccode2 %in% us), subloc = as.integer(cow_any == 1 & ccode1 %in% su & ccode2 %in% su & !(ccode1 %in% us & ccode2 %in% us)),
  dyally = as.integer(cow_any == 1 & usbloc == 0 & subloc == 0),
  contig = as.integer(conttype %in% 1:2), ldist = log(capdist + 1), ratio = pmax(cinc1, cinc2) / (cinc1 + cinc2),
  onemaj = as.integer(cowmaj1 + cowmaj2 == 1), twomaj = as.integer(cowmaj1 + cowmaj2 == 2))
d <- btscs_splines(d)
rhs0 <- c("lowdem", "highdem", "cwrival", "usbloc", "subloc", "dyally", "contig", "ldist", "ratio", "onemaj", "twomaj", "py", "spl1", "spl2", "spl3")
vars <- c("mid", rhs0, "dyadid")
cat("dyad-years 1990-2000:", nrow(d), " missing:\n"); print(sapply(d[, vars], function(x) sum(is.na(x))))
d0 <- d[complete.cases(d[, vars]), ]
cat("complete-case N:", nrow(d0), "(published 144,861)  onsets", sum(d0$mid), "\n")
TERM <- "lowdem"
fitter <- function(r, data) fit_bin("mid", r, data, "logit", "~dyadid", terms = TERM)
full <- fitter(rhs0, d0)
cat("\nREBUILT vs published: lowdem b", full[[TERM]][["b"]], "(pub -0.0280) se", full[[TERM]][["se"]], "(pub 0.0147) N", attr(full, "n"), "(pub 144,861)\n")
ds <- common_sample(d, vars, "ccode1", "ccode2", "year"); rm(d, d0); invisible(gc())
M <- run_six(rhs0, ds, fitter)
if (Sys.getenv("SWEEP_SENS") == "1") {
  source(file.path(CODE, "hist/sens_helper.R"))
  sens_run(id = "P9_084", dat = ds, y = "mid", tie = TERM, ctrl = setdiff(rhs0, TERM), fe = NULL,
           fit = function(dat, add) fitter(c(rhs0, add), dat)[[TERM]])
}
saveRDS(list(M = M, full = full), file.path(W, "M_P9_084.rds"))
meta <- list(pid = "P9_084", cite = "Park 2013", journal = "CMPS", year = 2013, doi = "10.1177/0738894212473927", family = "DEMOCRACY")
nt <- sprintf(paste0("REBUILT (no archive): peacesciencer all dyads 1990-2000, GML MID onset; Polity IV 2014 lower/higher polity2; Cold War bloc controls from ",
  "COW alliance v4.1 ties with the US/USSR in 1989; Gowa controls; btscs peace years + 3 splines; logit clustered on dyad. Full-sample rebuilt b=%.6f se=%.6f ",
  "N=%d vs published -0.0280 (0.0147) N 144,861. Common N=%d."), full[[TERM]][["b"]], full[[TERM]][["se"]], attr(full, "n"), nrow(ds))
write_p8(meta, TERM, "abstract: joint democracy promotes peace independently of Cold War interests; the democratic peace exists after the Cold War",
  "Table 2 panel B Model 2 (weak link; logit, MID onset, all dyads 1990-2000, Cold War interest controls)", "logit, clustered SE, duration controls",
  "feglm logit, SE clustered on dyad, all six steps", -0.0280, 0.0147, 144861, FALSE, TRUE, M, "Y (peace years + spline)", nt)
dir.create(file.path(OUT, "results_p9"), recursive = TRUE, showWarnings = FALSE)
src <- file.path(OUT, "results_p8", "P9_084__lowdem.csv"); stopifnot(file.exists(src), file.rename(src, file.path(OUT, "results_p9", basename(src))))
sink()
