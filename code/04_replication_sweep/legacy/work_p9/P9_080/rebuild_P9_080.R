# P9_080 Johnson, Joiner et al. 2019 (online; CMPS 38(2), 2021) "Power changes, alliance credibility, and extended deterrence"
#   (10.1177/0738894218824735). TIE (property of the target's alliance ties). REBUILD (protocol B; PARTIAL).
# Article (staged PDF Tables 1-2): probit, SEs clustered by dyad, peace years + square + cube; directed dyad-years 1816-2000 in which the
#   potential target has defensive allies; DV = potential challenger initiates a MID against the potential target (Maoz 2005 dyadic MIDs).
#   Key IV (Table 2 Model 1): any of the target's allies is >= 30% stronger or weaker (CINC) than in the alliance's formation year.
#   Controls: relative strength of the target's side; allies' average S-score; peacetime coordination in the target's alliance (ATOP);
#   challenger has a relevant offensive alliance; challenger has a relevant neutrality pact; ln distance; joint democracy (both Polity2 >= 6);
#   challenger-target S-score. Legend *p<.05 **p<.01.
# HEADLINE (fixed before estimating): Table 2 Model 1, Any ally 30% power change 0.238 (0.042) **, N 559,164.
# PARTIAL rebuild: variable construction copied from P9_035 (same authors' design, same ATOP/peacesciencer inputs) with the 30% threshold.
#   Deviations (minor controls): peacetime coordination (ATOP alliance-level military-coordination coding, not in the dyad-year ATOP file) is
#   DROPPED; "relevant" offensive/neutrality pacts use any ATOP offense/neutrality partner other than the target (Johnson-Leeds applicability
#   conditions not reproducible), as in P9_035. Allies' average S = mean Haege srsvaa between the target and its defense allies (excl. C).
CODE <- Sys.getenv("SWEEP_CODE"); stopifnot(nzchar(CODE))
source(file.path(CODE, "work_p8/common_b3.R"))
suppressPackageStartupMessages(library(peacesciencer))
set.seed(123)
W <- file.path(OUT, "work_p9/P9_080"); dir.create(W, recursive = TRUE, showWarnings = FALSE)
sink(file.path(W, "log_P9_080.txt"), split = TRUE)

# ---- 1. directed dyad-years 1816-2000: MIDs, peace years, distance, CINC, S ----
fp <- system.file("extdata", "dyadic_fp_similarity.rds", package = "peacesciencer"); stopifnot(nzchar(fp), file.exists(fp))
d <- create_dyadyears(directed = TRUE, subset_years = 1816:2000) %>% add_cow_mids() %>% add_peace_years() %>%
  add_capital_distance() %>% add_nmc() %>% add_democracy() %>% as.data.frame()
n0 <- nrow(d); stopifnot(anyDuplicated(d[, c("ccode1", "ccode2", "year")]) == 0)
s <- readRDS(fp) %>% as.data.frame() %>% select(ccode1, ccode2, year, srsvaa); stopifnot(anyDuplicated(s[, c("ccode1", "ccode2", "year")]) == 0)
d <- d %>% left_join(s, by = c("ccode1", "ccode2", "year")); stopifnot(nrow(d) == n0)
sall <- s
d <- d %>% select(ccode1, ccode2, year, cowmidonset, sidea1, orig1, hostlev2, cowmidspell, capdist, cinc1, cinc2, srsvaa, polity21, polity22)
cat("directed dyad-years 1816-2000:", n0, "\n")

# ---- 2. ATOP directed alliance edges (atop_alliance holds both orders) ----
a <- as.data.frame(atop_alliance) %>% filter(year >= 1800, year <= 2000) %>% select(x = ccode1, y = ccode2, year, def = atop_defense, off = atop_offense, neu = atop_neutral)
stopifnot(anyDuplicated(a[, c("x", "y", "year")]) == 0, sum(a$x < a$y) == sum(a$x > a$y))
cinc <- as.data.frame(cow_nmc) %>% select(ccode, year, cinc)
# defense edges T(x) -> ally A(y): formation year = first year of the current consecutive defense spell
de <- a %>% filter(def == 1) %>% arrange(x, y, year) %>% group_by(x, y) %>%
  mutate(sp = cumsum(c(1, diff(year) != 1))) %>% group_by(x, y, sp) %>% mutate(fy = min(year)) %>% ungroup() %>% select(x, y, year, fy)
de <- de %>% left_join(cinc %>% rename(ca = cinc), by = c("y" = "ccode", "year")) %>% left_join(cinc %>% rename(cf = cinc), by = c("y" = "ccode", "fy" = "year")) %>%
  mutate(pc = as.numeric(abs(ca / cf - 1) >= 0.30), age = year - fy) %>%
  left_join(sall %>% rename(sxy = srsvaa), by = c("x" = "ccode1", "y" = "ccode2", "year")) %>% as.data.frame()
cat("defense edges:", nrow(de), " power-change share:", round(mean(de$pc, na.rm = TRUE), 3), " missing pc:", sum(is.na(de$pc)), "\n")
# per target-year aggregates over ALL its defense allies (corrected below when the challenger is one of them)
ta <- de %>% group_by(x, year) %>% arrange(desc(age), .by_group = TRUE) %>%
  summarise(n_al = n(), npc = sum(pc, na.rm = TRUE), nna = sum(is.na(pc)), scinc = sum(ca, na.rm = TRUE), ssum = sum(sxy, na.rm = TRUE), sn = sum(!is.na(sxy)),   # missing CINC: ally contributes 0
             age1 = age[1], arg1 = y[1], age2 = if (n() > 1) age[2] else NA_real_, .groups = "drop")
# challenger offense / neutrality partners and offense partners' CINC
ca_ <- a %>% left_join(cinc, by = c("y" = "ccode", "year")) %>% group_by(x, year) %>%
  summarise(noff = sum(off), nneu = sum(neu), soff = sum(off * coalesce(cinc, 0)), .groups = "drop")
# ---- 3. merge onto directed dyads: T = ccode2, C = ccode1 ----
d <- d %>% left_join(ta, by = c("ccode2" = "x", "year")) %>% left_join(ca_, by = c("ccode1" = "x", "year")) %>%
  left_join(a %>% select(x, y, year, defCT = def, offCT = off, neuCT = neu), by = c("ccode2" = "x", "ccode1" = "y", "year")) %>%   # T's pact with C
  left_join(de %>% select(x, y, year, pcC = pc, ageC = age, caC = ca, sC = sxy), by = c("ccode2" = "x", "ccode1" = "y", "year"))
stopifnot(nrow(d) == n0)
d <- d %>% mutate(across(c(n_al, npc, nna, scinc, noff, nneu, soff, defCT, offCT, neuCT), ~coalesce(.x, 0)), isC = coalesce(defCT, 0) == 1) %>%
  mutate(n_rel = n_al - isC,                                          # relevant allies exclude the challenger
         allypc = ifelse(npc - ifelse(isC, coalesce(pcC, 0), 0) > 0, 1,   # [2026-11-14 P3] unknown ally power change -> NA, not 0
                         ifelse(nna - as.numeric(isC & is.na(pcC)) > 0, NA_real_, 0)),
         allycinc = scinc - ifelse(isC, coalesce(caC, 0), 0),
         allyage = ifelse(isC & arg1 == ccode1, age2, age1),
         chal_off = as.numeric(noff - offCT > 0), chal_neu = as.numeric(nneu - neuCT > 0),
         offcinc = soff - ifelse(offCT == 1, cinc2, 0),
         pwin_t = (cinc2 + allycinc) / (cinc2 + allycinc + cinc1 + offcinc),
         init = as.numeric(coalesce(cowmidonset, 0) == 1 & coalesce(sidea1, 0) == 1 & coalesce(orig1, 0) == 1),
         esc = ifelse(init == 1, as.numeric(hostlev2 >= 3), NA_real_),
         ldist = log(capdist), sscore = srsvaa, py = cowmidspell, py2 = py^2, py3 = py^3, dyadid = 1000 * ccode1 + ccode2,
         allyS = (ssum - ifelse(isC & !is.na(sC), sC, 0)) / (sn - as.numeric(isC & !is.na(sC))),
         jdem = as.numeric(polity21 >= 6 & polity22 >= 6))
d <- d[d$n_rel >= 1, ]
cat("dyad-years with >= 1 relevant defensive ally of the target:", nrow(d), "(paper 559,164 after missing data)\n")
cat("initiations:", sum(d$init), "")

d$allyS[!is.finite(d$allyS)] <- NA
cat("means: allypc", round(mean(d$allypc), 3), " pwin_t", round(mean(d$pwin_t, na.rm = TRUE), 3), " chal_off", round(mean(d$chal_off), 3),
    " chal_neu", round(mean(d$chal_neu), 3), " allyage", round(mean(d$allyage, na.rm = TRUE), 2), "\n")
d <- as.data.frame(d)

# ---- 4. probit: initiation ----
rhs0 <- c("allypc", "pwin_t", "allyS", "chal_off", "chal_neu", "ldist", "jdem", "sscore", "py", "py2", "py3")
vars <- c("init", rhs0, "dyadid")
cat("missing:\n"); print(sapply(d[, vars], function(x) sum(is.na(x))))
d0 <- d[complete.cases(d[, vars]), ]
cat("complete-case N:", nrow(d0), "(published 559,164)  initiations", sum(d0$init), "\n")
TERM <- "allypc"
fitter <- function(r, data) fit_bin("init", r, data, "probit", "~dyadid", terms = TERM)
full <- fitter(rhs0, d0); rm(d0); invisible(gc())
cat("\nREBUILT vs published: allypc b", full[[TERM]][["b"]], "(pub 0.238) se", full[[TERM]][["se"]], "(pub 0.042) N", attr(full, "n"), "(pub 559,164)\n")
ds <- common_sample(d, vars, "ccode1", "ccode2", "year"); rm(d); invisible(gc())
M <- run_six(rhs0, ds, fitter)
if (Sys.getenv("SWEEP_SENS") == "1") {
  source(file.path(CODE, "hist/sens_helper.R"))
  sens_run(id = "P9_080", dat = ds, y = "init", tie = TERM, ctrl = setdiff(rhs0, TERM), fe = NULL,
           fit = function(dat, add) fitter(c(rhs0, add), dat)[[TERM]])
}
saveRDS(list(M = M, full = full), file.path(W, "M_P9_080.rds"))
meta <- list(pid = "P9_080", cite = "Johnson et al. 2019", journal = "CMPS", year = 2019, doi = "10.1177/0738894218824735", family = "TIE")
nt <- sprintf(paste0("REBUILT, PARTIAL (no archive): P9_035 construction (peacesciencer COW directed dyads 1816-2000, ATOP 5 defense allies of the target other than ",
  "the challenger, COW MID 5.0 initiation = onset with C originator on side A) with a 30%% CINC-change threshold; allies' mean Haege S; ",
  "peacetime coordination DROPPED; relevant offense/neutrality approximated; joint democracy polity2 >= 6; ln capital distance; C-T S (srsvaa); ",
  "COW peace years cubic; probit clustered on directed dyad. Full-sample rebuilt b=%.6f se=%.6f N=%d vs published 0.238 (0.042) N 559,164. Common N=%d."),
  full[[TERM]][["b"]], full[[TERM]][["se"]], attr(full, "n"), nrow(ds))
write_p8(meta, TERM, "abstract: states whose allies experienced significant power changes are more likely to be targeted; Table 2 Model 1",
  "Table 2 Model 1 (probit, MID initiation vs targets with defensive allies, directed dyads 1816-2000)", "probit, SE clustered on dyad, peace years cubic",
  "feglm probit, SE clustered on directed dyad, all six steps", 0.238, 0.042, 559164, TRUE, TRUE, M, "Y (peace years cubic)", nt)
dir.create(file.path(OUT, "results_p9"), recursive = TRUE, showWarnings = FALSE)
src <- file.path(OUT, "results_p8", "P9_080__allypc.csv"); stopifnot(file.exists(src), file.rename(src, file.path(OUT, "results_p9", basename(src))))
sink()
