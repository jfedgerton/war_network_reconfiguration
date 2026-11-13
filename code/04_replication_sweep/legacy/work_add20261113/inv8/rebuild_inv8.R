# Rebuild inv8: Long, Nordstrom & Baek 2007 JOP (10.1111/j.1468-2508.2007.00611.x), "Allying for Peace: Treaty Obligations
#   and Conflict between Allies" -- alliance provisions -> hazard of a MID in which both allies use force
# Headline: Table 1 (Cox PH with time-varying covariates, Breslow ties, SE clustered on episode), key term = Non-aggression
#   published b = -0.463 (SE 0.145), p < .05; 3,695 episodes, 135 failures, 57,819 allied dyad-years (1816-2001).
# Protocol B (no author files): built from ATOP v3.0 (the version the authors used) and peacesciencer, round-2 match rule.
# Construction (article pp. 1108-1111; ATOP v3.0 codebook):
#   unit        allied dyad-years from atop3_0dy.csv (all ATOP alliances in force for the pair, ids in atopid1-5)
#   episode     consecutive allied years of a dyad; a gap starts a new episode; years after the first failure are dropped
#   failure     onset of a dyadic MID in which both sides reach hostility level 4+ (use of force); GML dyadic MIDs (peacesciencer)
#               stand in for Maoz's dyadic MID 2.0
#   Non-aggression          0-3: NONAGG (dyad-year) + NOTAIDEN > 0 (either member) + INTERV == 1 (mutual nonintervention)
#   Mediate Disputes        MEDARB == 1
#   Military Institutionalization 0-3: MILCON >= 2 (peacetime contact) + INTCOM == 1 + an organization whose purpose is
#                           military cooperation (ORGPURP1/2 == 1)
#   Organization            ORGAN1/2 == 3 (stand-alone organization with a permanent bureaucracy)
#   Arbitration Commission  ORGPURP1/2 == 5 (organization to resolve disputes among members)
#   provisions are the maximum over all alliances in force for the dyad-year
#   controls    relative capabilities (larger CINC share), shared land border (conttype 1), joint democracy (both polity2 >= 6),
#               similarity (global weighted S, valued alliances: srswvas), wartime alliance (ATOP WARTIME of any alliance in
#               force), major power (either state)
ROOT <- Sys.getenv("SWEEP_DATA"); OUT <- Sys.getenv("SWEEP_OUT"); CODE <- Sys.getenv("SWEEP_CODE")
stopifnot(nzchar(ROOT), nzchar(OUT), nzchar(CODE))
suppressPackageStartupMessages({library(survival); library(dplyr); library(peacesciencer)})
set.seed(123)
WO <- file.path(OUT, "work_add20261113/inv8"); dir.create(WO, recursive = TRUE, showWarnings = FALSE)
sink(file.path(WO, "log_inv8.txt"), split = TRUE)
AD <- file.path(ROOT, "_incoming/Long2007/ATOP v3.0")

# 1. ATOP v3.0: dyad-years, alliance- and member-level provisions
dy <- read.csv(file.path(AD, "atop3_0dy.csv"))
al <- read.csv(file.path(AD, "atop3_0a.csv"))
me <- read.csv(file.path(AD, "atop3_0m.csv"))
stopifnot(all(dy$atopally == 1))
dy <- dy %>% filter(year >= 1816, year <= 2001) %>% mutate(c1 = pmin(mem1, mem2), c2 = pmax(mem1, mem2))
stopifnot(!anyDuplicated(dy[, c("c1", "c2", "year")]))
cat("allied dyad-years 1816-2001:", nrow(dy), "\n")

alv <- al %>% transmute(atopid, wartime = as.integer(wartime == 1), medarb = as.integer(medarb == 1),
  milcon2 = as.integer(milcon >= 2), intcom = as.integer(intcom == 1),
  milorg = as.integer(orgpurp1 %in% 1 | orgpurp2 %in% 1),
  organ = as.integer(organ1 %in% 3 | organ2 %in% 3), arbcom = as.integer(orgpurp1 %in% 5 | orgpurp2 %in% 5),
  interv1 = as.integer(interv %in% 1))
mev <- me %>% group_by(atopid, member) %>% summarise(notaiden = as.integer(any(notaiden > 0, na.rm = TRUE)), .groups = "drop")

long <- bind_rows(lapply(1:5, function(k) dy %>% select(c1, c2, year, atopid = !!paste0("atopid", k)))) %>% filter(!is.na(atopid))
long <- long %>% left_join(alv, by = "atopid") %>%
  left_join(mev %>% rename(c1 = member, na1 = notaiden), by = c("atopid", "c1")) %>%
  left_join(mev %>% rename(c2 = member, na2 = notaiden), by = c("atopid", "c2")) %>%
  mutate(notaid = as.integer(coalesce(na1, 0L) == 1 | coalesce(na2, 0L) == 1))
cat("alliance ids not found in the alliance file:", sum(!long$atopid %in% al$atopid), "\n")
stopifnot(all(long$atopid %in% al$atopid))
mx <- function(x) if (all(is.na(x))) NA_integer_ else max(x, na.rm = TRUE)
prov <- long %>% group_by(c1, c2, year) %>%
  summarise(across(c(wartime, medarb, milcon2, intcom, milorg, organ, arbcom, interv1, notaid), mx), .groups = "drop")
d <- dy %>% select(c1, c2, year, nonagg) %>% inner_join(prov, by = c("c1", "c2", "year")) %>%
  mutate(nonaggression = nonagg + notaid + interv1, milinst = milcon2 + intcom + milorg)
stopifnot(nrow(d) == nrow(dy))

# 2. Failure: onset year of a GML dyadic MID with both sides at hostility level 4+ (an ongoing MID carried into a new
#    alliance is not a failure of that alliance). [2026-11-13] counting ongoing MID-years as failures gives b = -0.170
#    (p = .14, 161 failures); onsets give 149 failures (published 135).
mid <- gml_mid_ddydisps %>% mutate(c1 = pmin(ccode1, ccode2), c2 = pmax(ccode1, ccode2)) %>%
  filter(gmlmidonset == 1, hostlev1 >= 4, hostlev2 >= 4) %>% distinct(c1, c2, year) %>% mutate(fail = 1L)
d <- d %>% left_join(mid, by = c("c1", "c2", "year")) %>% mutate(fail = coalesce(fail, 0L))

# 3. Episodes (consecutive allied years), time at risk, drop years after the first failure in an episode
d <- d %>% arrange(c1, c2, year) %>% group_by(c1, c2) %>%
  mutate(newep = is.na(lag(year)) | year != lag(year) + 1, ep = cumsum(newep)) %>% ungroup() %>%
  mutate(episode = paste(c1, c2, ep, sep = "_")) %>% group_by(episode) %>%
  mutate(t = year - min(year) + 1, prior = cumsum(lag(fail, default = 0L))) %>% ungroup() %>% filter(prior == 0)

# 4. Controls
nmc <- cow_nmc %>% select(ccode, year, cinc)
pol <- ccode_democracy %>% select(ccode, year, polity2)
cont <- cow_contdir %>% mutate(c1 = pmin(ccode1, ccode2), c2 = pmax(ccode1, ccode2), y0 = as.integer(format(as.Date(stdate), "%Y")), y1 = as.integer(format(as.Date(enddate), "%Y"))) %>%
  filter(conttype == 1) %>% select(c1, c2, y0, y1) %>% distinct()
maj <- cow_majors %>% select(ccode, styear, endyear)
ismaj <- function(cc, yr) vapply(seq_along(cc), function(i) any(maj$ccode == cc[i] & maj$styear <= yr[i] & maj$endyear >= yr[i]), logical(1))
fps <- readRDS(system.file("extdata", "dyadic_fp_similarity.rds", package = "peacesciencer"))
stopifnot(nrow(fps) > 0)
fps <- fps %>% mutate(a = pmin(ccode1, ccode2), b = pmax(ccode1, ccode2)) %>% group_by(a, b, year) %>%
  summarise(sim = mean(srswvas), .groups = "drop")
d <- d %>% left_join(nmc %>% rename(c1 = ccode, cinc1 = cinc), by = c("c1", "year")) %>%
  left_join(nmc %>% rename(c2 = ccode, cinc2 = cinc), by = c("c2", "year")) %>%
  left_join(pol %>% rename(c1 = ccode, pol1 = polity2), by = c("c1", "year")) %>%
  left_join(pol %>% rename(c2 = ccode, pol2 = polity2), by = c("c2", "year")) %>%
  left_join(fps %>% rename(c1 = a, c2 = b), by = c("c1", "c2", "year")) %>%
  mutate(relcap = pmax(cinc1, cinc2) / (cinc1 + cinc2), jointdem = as.integer(pol1 >= 6 & pol2 >= 6))
cy <- d %>% select(c1, c2, year) %>% left_join(cont, by = c("c1", "c2"), relationship = "many-to-many") %>%
  mutate(on = !is.na(y0) & year >= y0 & year <= y1) %>% group_by(c1, c2, year) %>% summarise(border = as.integer(any(on)), .groups = "drop")
d <- d %>% left_join(cy, by = c("c1", "c2", "year")) %>% mutate(major = as.integer(ismaj(c1, year) | ismaj(c2, year)))

TIE  <- "nonaggression"
KEYS <- c("arbcom", "organ", "nonaggression", "medarb", "milinst")
CTRL <- c("arbcom", "organ", "medarb", "milinst", "relcap", "border", "jointdem", "sim", "wartime", "major")
d <- d %>% mutate(tstart = t - 1)
cc <- d[complete.cases(d[, c("fail", "t", TIE, CTRL)]), ]
cat("rows", nrow(d), "| complete cases", nrow(cc), "| episodes", n_distinct(cc$episode), "| failures", sum(cc$fail),
    "| published: 57,819 obs, 3,695 episodes, 135 failures\n")

f_rhs <- paste(c(TIE, CTRL), collapse = " + ")
fit <- function(rhs_add, dat) {
  m <- coxph(as.formula(paste("Surv(tstart, t, fail) ~", f_rhs, rhs_add)), data = dat, ties = "breslow", cluster = episode)
  ct <- summary(m)$coefficients
  print(round(ct[TIE, , drop = FALSE], 5))
  list(b = ct[TIE, "coef"], se = ct[TIE, "robust se"], p = ct[TIE, "Pr(>|z|)"], n = m$n, ev = m$nevent, ct = ct)
}

# 5. Rebuilt Table 1 on the full rebuilt sample and the round-2 verdict
cat("\n==== Rebuilt Table 1 ====\n")
r_pub <- fit("", cc)
print(round(r_pub$ct[, c("coef", "robust se", "Pr(>|z|)")], 4))
b_pub <- -0.463; se_pub <- 0.145
source(file.path(CODE, "work_p8/round2_rule.R"))
verdict <- rebuild_verdict("inv8", TIE, r_pub$b, r_pub$p, b_pub, TRUE, se_pub)
cat(sprintf("REBUILD %s: %.4f (se %.4f, p %.4f) N %d events %d | published %.3f (%.3f) | verdict %s\n",
            TIE, r_pub$b, r_pub$se, r_pub$p, r_pub$n, r_pub$ev, b_pub, se_pub, verdict))
n_orig <- r_pub$n

# 6. History (undirected dyad, outcome year = year); restrict to 1822-2011 (sample ends 2001)
source(file.path(CODE, "hist/add_history.R"))
cc <- add_history(cc, "c1", "c2", "year")
cc <- filter(cc, year >= 1822, year <= 2011)
cc <- cc[complete.cases(cc[, c("hist_dyad_w15", "hist_mono_w15", "hist_dyad_w26", "hist_mono_w26", "bm_cinc_w15", "bm_democ_w15")]), ]
cat("\nRows after year restriction + complete cases:", nrow(cc), "(rebuilt N", n_orig, ") | failures", sum(cc$fail), "\n")
stopifnot(nrow(cc) > 0.80 * n_orig)

cat("\n==== m0 ====\n");  r0  <- fit("", cc)
cat("\n==== mD ====\n");  rD  <- fit("+ hist_dyad_w15", cc)
cat("\n==== mDM ====\n"); rDM <- fit("+ hist_dyad_w15 + hist_mono_w15", cc)
cat("\n==== m26 ====\n"); r26 <- fit("+ hist_dyad_w26 + hist_mono_w26", cc)
cat("\n==== bC ====\n");  rC  <- fit("+ bm_cinc_w15", cc)
cat("\n==== bD ====\n");  rBD <- fit("+ bm_democ_w15", cc)
stopifnot(length(unique(c(r0$n, rD$n, rDM$n, r26$n, rC$n, rBD$n))) == 1, r0$n == nrow(cc))
if (Sys.getenv("SWEEP_SENS") == "1") {
  source(file.path(CODE, "hist/sens_helper.R"))
  sens_run(id = "inv8", dat = cc, y = "fail", tie = TIE, ctrl = CTRL, fe = NULL,
           fit = function(dat, add) { r <- fit(if (length(add)) paste("+", paste(add, collapse = " + ")) else "", dat)
             c(b = unname(r$b), se = unname(r$se), p = unname(r$p)) })
}
pct <- function(x) 100 * (x$b - r0$b) / abs(r0$b)
out <- data.frame(cid = "inv8", short_cite = "Long et al. 2007", key_term = TIE, reproduced = verdict, b_pub = b_pub, b_repro = r_pub$b,
  estimator = "Cox PH (Breslow), robust SE clustered by episode (survival::coxph)", tie_var = TIE,
  n_orig = n_orig, n_used = r0$n, b0 = r0$b, se0 = r0$se, p0 = r0$p, b_dm = rDM$b, se_dm = rDM$se, p_dm = rDM$p,
  pct_dm = pct(rDM), b_d = rD$b, pct_d = pct(rD), b_26 = r26$b, se_26 = r26$se, p_26 = r26$p, pct_26 = pct(r26),
  pct_bC = pct(rC), pct_bD = pct(rBD), author_has_history_control = "N",
  notes = sprintf("Table 1 Non-aggression; protocol-B rebuild from ATOP v3.0 + GML MIDs; dropped %d obs vs rebuilt sample.", n_orig - r0$n))
stopifnot(all(is.finite(unlist(out[, c("b0", "se0", "p0", "b_dm", "p_dm", "b_26", "p_26", "pct_bC", "pct_bD")]))))
print(Filter(is.numeric, out))
dir.create(file.path(OUT, "results"), recursive = TRUE, showWarnings = FALSE)
write.csv(out, file.path(OUT, "results/inv8.csv"), row.names = FALSE)
sink()
