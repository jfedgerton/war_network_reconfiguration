# P9_086 / P9_087 McDonald 2015, IO 69(3) "Great Powers, Hierarchy, and Endogenous Regimes: Rethinking the Domestic Causes of Peace"
#   (10.1017/S0020818315000120). REBUILD (protocol B; no Dataverse deposit or script-fetchable supplement found 2026-10-11).
#   Two headline terms (one article, two families; precedent Gartzke & Hewitt 2010 = C6/99):
#   P9_086 DEMOCRACY: Table 2 Model 1 (1816-2000, DV MID): DEMOCRACYL -0.0188 (0.0074) **, N 525,114.
#   P9_087 TIE: Table 6 col 1 (1919-2000, DV MID): GPDEM x DEMOCRACYL -0.204 (0.078) **, N 468,189.
# Article (staged PDF pp. 571-573, 580-582): non-directed dyad-years; DV = outbreak of a new MID between dyadic participants that are both
#   originating disputants (MID 3.1); DEMOCRACYL = lower Polity2; GREATPOWER = at least one COW great power; CAPRATIO = ln(stronger CINC /
#   weaker CINC); CONTIGUITY = land border or within 150 miles of water; DISTANCE = ln capital distance; DYADIC ALLIANCE = common COW alliance;
#   INTEREST SIM = unweighted global S (alliance portfolios); Carter-Signorino cubic in peace years; logit, robust SEs clustered on dyad;
#   two-tailed, * p<=.10 ** p<=.05 *** p<=.01.
#   GPDEM (Table 6 note) = both members democracies AND (at least one allied to a great power, or democratized 1918-1925 [1919-45 years], or
#   democratized 1988-1993 [1992-2000 years]).
# Rebuild inputs: {peacesciencer} COW dyad-years (GML MIDs: onset with orig1 = orig2 = 1; contiguity 1-4; capital distance; CINC; majors;
#   Haege srsvaa as the unweighted global S); Polity IV 2014 polity2 (work_p8/P8_049); COW alliance v4.1 dyad-years (common alliance, and
#   defense pacts with great powers). Democracy for GPDEM = polity2 >= 6 (the article's Polity three-tier coding); "democratized" = first year
#   of a democratic spell (polity2 >= 6 after < 6 or system entry) falls in the window.
ROOT <- Sys.getenv("SWEEP_DATA"); OUT <- Sys.getenv("SWEEP_OUT"); CODE <- Sys.getenv("SWEEP_CODE")
stopifnot(nzchar(ROOT), nzchar(OUT), nzchar(CODE))
source(file.path(CODE, "work_p8/common_b3.R")); source(file.path(CODE, "work_p9/rebuild_common_p9r.R"))
set.seed(123)
W <- file.path(OUT, "work_p9/P9_086"); dir.create(W, recursive = TRUE, showWarnings = FALSE)
sink(file.path(W, "log_P9_086.txt"), split = TRUE)
d <- ps_dyads(1816:2000, trade = FALSE, fpsim = TRUE) %>% add_cow_ally()
d <- d %>% mutate(dyadid = 1000 * ccode1 + ccode2, mid = as.integer(gmlmidonset == 1 & orig1 %in% 1 & orig2 %in% 1))
d <- py_calc(d, "dyadid", "mid")
p <- polity_cow("polity2"); n0 <- nrow(d)
d <- d %>% left_join(p %>% rename(p2a = polity2), by = c("ccode1" = "ccode", "year")) %>% left_join(p %>% rename(p2b = polity2), by = c("ccode2" = "ccode", "year"))
stopifnot(nrow(d) == n0)
# great-power defense allies (COW v4.1 defense pacts with a COW major)
a <- read.csv(file.path(ROOT, "work_p8/P8_049/alliance_v4_1_by_dyad_yearly.csv"))
mj <- suppressMessages(create_stateyears(system = "cow", subset_years = 1816:2000) %>% add_cow_majors()) %>% as.data.frame() %>% filter(cowmaj == 1) %>% select(ccode, year)
gpa <- bind_rows(a %>% filter(defense == 1) %>% transmute(s = ccode1, g = ccode2, year), a %>% filter(defense == 1) %>% transmute(s = ccode2, g = ccode1, year)) %>%
  semi_join(mj, by = c("g" = "ccode", "year")) %>% distinct(ccode = s, year) %>% mutate(gpally = 1L)
# democratization years: first year of a democratic spell
pz <- p %>% arrange(ccode, year) %>% group_by(ccode) %>% mutate(dem = as.integer(polity2 >= 6), prev = dplyr::lag(dem), first = !is.na(dem) & dem == 1 & (is.na(prev) | prev == 0)) %>% ungroup()
dz <- pz %>% filter(first) %>% transmute(ccode, dzy = year)
w1 <- unique(dz$ccode[dz$dzy >= 1918 & dz$dzy <= 1925]); w2 <- unique(dz$ccode[dz$dzy >= 1988 & dz$dzy <= 1993])
d <- d %>% left_join(gpa %>% rename(gpa1 = gpally), by = c("ccode1" = "ccode", "year")) %>% left_join(gpa %>% rename(gpa2 = gpally), by = c("ccode2" = "ccode", "year"))
stopifnot(nrow(d) == n0)
d <- d %>% mutate(demL = pmin(p2a, p2b), gp = as.integer(cowmaj1 == 1 | cowmaj2 == 1), caprat = log(pmax(cinc1, cinc2) / pmin(cinc1, cinc2)),
  contig = as.integer(conttype %in% 1:4), ldist = log(capdist + 1), ally = cow_any, ssim = srsvaa, py2 = py^2, py3 = py^3,
  bothdem = as.integer(p2a >= 6 & p2b >= 6),
  dz1 = (year >= 1919 & year <= 1945) & (ccode1 %in% w1 | ccode2 %in% w1), dz2 = (year >= 1992) & (ccode1 %in% w2 | ccode2 %in% w2),
  gpdem = as.integer(bothdem == 1 & (coalesce(gpa1, 0L) == 1 | coalesce(gpa2, 0L) == 1 | dz1 | dz2)), gpdem_x_dem = gpdem * demL)
d$caprat[!is.finite(d$caprat)] <- NA
ctrl <- c("gp", "caprat", "contig", "ldist", "ally", "ssim", "py", "py2", "py3")
SPEC <- list(P9_086 = list(term = "demL", rhs = c("demL", ctrl), yrs = 1816:2000, pub = c(-0.0188, 0.0074, 525114), fam = "DEMOCRACY",
                           hm = "Table 2 Model 1 (logit, MID onset between originators, all dyads 1816-2000)"),
             P9_087 = list(term = "gpdem_x_dem", rhs = c("demL", "gpdem", "gpdem_x_dem", ctrl), yrs = 1919:2000, pub = c(-0.204, 0.078, 468189), fam = "TIE",
                           hm = "Table 6 col 1 (logit, MID onset between originators, all dyads 1919-2000; GPDEM x DEMOCRACYL)"))
dir.create(file.path(OUT, "results_p9"), recursive = TRUE, showWarnings = FALSE)
for (pid in names(SPEC)) {
  sp <- SPEC[[pid]]; TERM <- sp$term
  dd <- d[d$year %in% sp$yrs, ]; vars <- c("mid", sp$rhs, "dyadid")
  cat("\n==== ", pid, " dyad-years:", nrow(dd), " missing:\n"); print(sapply(dd[, vars], function(x) sum(is.na(x))))
  d0 <- dd[complete.cases(dd[, vars]), ]
  cat("complete-case N:", nrow(d0), "(published", sp$pub[3], ")  onsets", sum(d0$mid), " gpdem dyad-years", sum(d0$gpdem), "\n")
  fitter <- function(r, data) fit_bin("mid", r, data, "logit", "~dyadid", terms = TERM)
  full <- fitter(sp$rhs, d0)
  cat("\nREBUILT vs published:", TERM, "b", full[[TERM]][["b"]], "(pub", sp$pub[1], ") se", full[[TERM]][["se"]], "(pub", sp$pub[2], ") N", attr(full, "n"), "\n")
  ds <- common_sample(dd, vars, "ccode1", "ccode2", "year"); rm(d0); invisible(gc())
  M <- run_six(sp$rhs, ds, fitter)
  if (Sys.getenv("SWEEP_SENS") == "1") {
    source(file.path(CODE, "hist/sens_helper.R"))
    if (pid == "P9_086") sens_run(id = pid, dat = ds, y = "mid", tie = TERM, ctrl = setdiff(sp$rhs, TERM), fe = NULL, fit = function(dat, add) fitter(c(sp$rhs, add), dat)[[TERM]])
    else sens_split(id = pid, dat = ds, fit = function(dat, add) fitter(c(sp$rhs, add), dat)[[TERM]])
  }
  b0 <- M$m0[[TERM]][["b"]]; p0 <- M$m0[[TERM]][["p"]]; pub_p <- 2 * pnorm(-abs(sp$pub[1] / sp$pub[2])); r_ <- b0 / sp$pub[1]
  v <- rebuild_verdict(pid, TERM, b0, p0, sp$pub[1], NA, sp$pub[2])   # round-2 rule, alpha_pub from registry.csv (P9_086 unregistered -> 0.05)
  meta <- list(pid = pid, cite = "McDonald 2015", journal = "IO", year = 2015, doi = "10.1017/S0020818315000120", family = sp$fam)
  nt <- sprintf(paste0("REBUILT (no archive): peacesciencer all dyads, GML MID onset with both states originators; lower polity2; COW major; log CINC ratio; contiguity 1-4; ",
    "log capital distance; COW v4.1 common alliance; Haege srsvaa S; peace years cubic; logit clustered on dyad. GPDEM per Table 6 note (polity2 >= 6; COW defense ",
    "pact with a COW major; democratization windows). Full-sample rebuilt b=%.6f se=%.6f N=%d vs published %g (%g) N %d. Common N=%d. Round-2 rule (alpha .05, 0.5-2x)."),
    full[[TERM]][["b"]], full[[TERM]][["se"]], attr(full, "n"), sp$pub[1], sp$pub[2], sp$pub[3], nrow(ds))
  write_p8(meta, TERM, "abstract and Tables 2/6: the democratic peace depends on great-power hierarchy", sp$hm,
    "logit, robust SE clustered on dyad, Carter-Signorino cubic", "feglm logit, SE clustered on dyad, all six steps",
    sp$pub[1], sp$pub[2], sp$pub[3], TRUE, TRUE, M, "Y (peace years cubic)", nt, verdict = v)
  src <- file.path(OUT, "results_p8", paste0(pid, "__", TERM, ".csv")); stopifnot(file.exists(src), file.rename(src, file.path(OUT, "results_p9", basename(src))))
  rm(ds, dd); invisible(gc())
}
sink()
