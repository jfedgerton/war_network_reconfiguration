# VERIFY 108: Boehmer, Gartzke & Nordstrom 2004 World Politics 57(1) -- matched to the author copy
#   (pages.ucsd.edu/~egartzke/publications/57.1boehmer.pdf), Table 2 (p.24): logit, MID onset, all dyads 1950-91, n = 109,965.
#   Model 2-3 "Joint Membership, Structured" -0.0854 (0.0429)**; Model 2-4 "Joint Membership, Interventionist" -0.0094 (0.0162).
#   Controls: floor/major-power/difference contention, diplomatic missions low, democracy low, interdependence low, capability ratio,
#   allies, contiguity (land or <=150 mi water), ln capital distance, major power dyad, peace years + 3 splines. SE type not reported.
# KEY-TERM FLAG: rebuild_108.R keyed IGOmaxst with b_pub -0.0854, but -0.0854 is the STRUCTURED coefficient (Model 2-3).
#   The authors' file labels IGOminst/IGOmedst/IGOmaxst = min./med./max. institutionalized IGOs, i.e. minimal/structured/
#   interventionist; IGOmaxst (interventionist) is published as -0.0094 (ns) in Model 2-4. Per protocol the existing key
#   (IGOmaxst, Model 2-4) is kept and the paper's headline term (IGOmedst, Model 2-3: "institutionalized IGOs reduce the risk of
#   militarized disputes", p.2) is ALSO run.
# Changes vs rebuild_108.R (paper text): (1) sample = Oneal-Russett all-dyads coverage ("We base our analysis on Oneal and Russett",
#   p.17): rows where the authors' O&R IGO variable smigoabi is non-missing (Table 1 all-dyads n = 149,404 = smigoabi coverage);
#   (2) polity2 from peacesciencer::ccode_democracy (add_democracy() v1.2.0 returns NA for USSR, Yugoslavia, etc.).
# Not fixable: Maoz DYMID onset -> GML onset; O&R trade/GDP -> COW trade (nominal) / SDP GDP (2011$); Polity III -> polity2;
#   BTSCS knots unreported (cubic kept); SE type unreported (dyad clusters kept).
suppressMessages({library(haven); library(peacesciencer)})
ROOT <- Sys.getenv("SWEEP_DATA"); OUT <- Sys.getenv("SWEEP_OUT"); CODE <- Sys.getenv("SWEEP_CODE")
stopifnot(nzchar(ROOT), nzchar(OUT), nzchar(CODE))
source(file.path(CODE, "vb4/vhelp.R"))
W <- file.path(ROOT, "work/108"); WO <- file.path(OUT, "work/108"); dir.create(WO, recursive = TRUE, showWarnings = FALSE); sink(file.path(WO, "log_verify_108.txt"), split = TRUE)
g <- read_dta(file.path(W, "BGN WP05 limited.dta")) %>% zap_labels() %>% as.data.frame() %>%
  select(statea, stateb, year, smigoabi, IGOmedst, IGOmaxst, cavgqsun, cavmjsun, cadifsun, dipsl)
g$lo <- pmin(g$statea, g$stateb); g$hi <- pmax(g$statea, g$stateb)
stopifnot(anyDuplicated(g[, c("lo","hi","year")]) == 0)
g <- g %>% filter(year >= 1950, year <= 1991)
cat("BGN rows 1950-91:", nrow(g), " with O&R coverage (smigoabi):", sum(!is.na(g$smigoabi)), "\n")
g <- g[!is.na(g$smigoabi), ]

nd <- create_dyadyears(directed = FALSE, subset_years = 1816:1991) %>% add_gml_mids()
nd <- nd %>% arrange(ccode1, ccode2, year) %>% group_by(ccode1, ccode2) %>%
  mutate(spell = cumsum(lag(gmlmidonset, default = 0) == 1)) %>% group_by(ccode1, ccode2, spell) %>%
  mutate(py = row_number() - 1) %>% ungroup() %>% filter(year >= 1950) %>% select(ccode1, ccode2, year, gmlmidonset, py)
trf <- system.file("extdata", "cow_trade_ddy.rds", package = "peacesciencer"); stopifnot(nzchar(trf), file.exists(trf))
tr <- readRDS(trf) %>%
  mutate(flow1 = ifelse(flow1 < 0, NA, flow1), flow2 = ifelse(flow2 < 0, NA, flow2),
         lo = pmin(ccode1, ccode2), hi = pmax(ccode1, ccode2)) %>%
  group_by(lo, hi, year) %>% summarise(trade = mean(flow1 + flow2), .groups = "drop")
pol <- ccode_democracy %>% select(ccode, year, polity2) %>% as.data.frame(); stopifnot(!anyDuplicated(pol[, 1:2]))
ctl <- create_dyadyears(directed = FALSE, subset_years = 1950:1991) %>%
  add_contiguity() %>% add_capital_distance() %>% add_nmc() %>% add_cow_majors() %>% {suppressWarnings(add_sdp_gdp(.))}
ctl <- ctl %>% left_join(pol %>% rename(polity21 = polity2), by = c("ccode1" = "ccode", "year")) %>%
  left_join(pol %>% rename(polity22 = polity2), by = c("ccode2" = "ccode", "year"))
ca <- peacesciencer::cow_alliance
ctl <- ctl %>% left_join(ca, by = c("ccode1","ccode2","year")) %>%
  left_join(tr, by = c("ccode1" = "lo", "ccode2" = "hi", "year" = "year")) %>%
  mutate(allies = as.numeric(coalesce(cow_defense, 0) + coalesce(cow_neutral, 0) + coalesce(cow_entente, 0) > 0),
         demlow = pmin(polity21, polity22),
         deplow = pmin(trade * 1e6 / exp(wbgdp2011est1), trade * 1e6 / exp(wbgdp2011est2)),
         lncapratio = log(pmax(cinc1, cinc2) / pmin(cinc1, cinc2)),
         contig = as.numeric(conttype %in% 1:4), lndist = log(capdist + 1),
         majdyad = as.numeric(cowmaj1 == 1 | cowmaj2 == 1)) %>%
  select(ccode1, ccode2, year, allies, demlow, deplow, lncapratio, contig, lndist, majdyad) %>%
  left_join(nd, by = c("ccode1","ccode2","year"))
d <- g %>% inner_join(ctl, by = c("lo" = "ccode1", "hi" = "ccode2", "year" = "year")) %>%
  mutate(py2 = py^2 / 10, py3 = py^3 / 100, dyad = lo * 1000 + hi) %>% as.data.frame()
rm(nd, ctl, tr, g); invisible(gc())
d <- d[is.finite(d$lncapratio) | is.na(d$lncapratio), ]
ctrl <- c("cavgqsun","cavmjsun","cadifsun","dipsl","demlow","deplow","lncapratio","allies","contig","lndist","majdyad","py","py2","py3")
allv <- c("gmlmidonset", "IGOmedst", "IGOmaxst", ctrl, "dyad")
x <- d[complete.cases(d[, allv]) & is.finite(d$lncapratio), ]
cat("rebuilt complete-case n =", nrow(x), "(published 109,965; previous rebuild 213,759)\n")
FF <- function(k) function(r, s) fit_bin("gmlmidonset", r, s, "logit", "~dyad", NULL, k)
r3 <- FF("IGOmedst")(c("IGOmedst", ctrl), x); r4 <- FF("IGOmaxst")(c("IGOmaxst", ctrl), x); rm(x)
s <- common_sample(d[is.finite(d$lncapratio), ], allv, "lo", "hi", "year"); rm(d); invisible(gc())
M4 <- run_six(c("IGOmaxst", ctrl), s, FF("IGOmaxst"))
M3 <- run_six(c("IGOmedst", ctrl), s, FF("IGOmedst"))

# ---- sensitivity hook: runs only under run_sensitivity.R (SWEEP_SENS=1), never in the main sweep ----
if (Sys.getenv("SWEEP_SENS") == "1") {
  source(file.path(CODE, "hist/sens_helper.R"))
  sens_run(id = "108", dat = s, y = "gmlmidonset", tie = "IGOmedst", ctrl = ctrl, fe = NULL,
           fit = function(dat, add) { r <- FF("IGOmedst")(c("IGOmedst", ctrl, add), dat)$IGOmedst; c(b = unname(r[["b"]]), se = unname(r[["se"]]), p = unname(r[["p"]])) })
}

mism <- "b_pub (-0.0854) belonged to Structured IGOs (Model 2-3), not IGOmaxst/interventionist (Model 2-4: -0.0094); sample not restricted to O&R all-dyads coverage (n 213,759 vs 109,965); add_democracy() NA polity2 for USSR/Yugoslavia etc.; GML onset instead of Maoz DYMID onset; COW/SDP trade-to-GDP; cubic peace years instead of BTSCS splines"
unf <- "Maoz DYMID onset, O&R trade/GDP and Polity III not available (GML, COW/SDP, polity2 used); BTSCS knots and SE type unreported"
write_verify(108, "Boehmer, Gartzke & Nordstrom 2004", "TIE", "IGOmaxst", "Table 2 Model 2-4 (interventionist IGOs, p.24)",
  -0.0094, 0.0162, 109965, "logit, MID onset, all dyads 1950-91, peace years + 3 splines; SE type not reported", FALSE, mism, "PARTIAL", unf, TRUE, M4,
  sprintf("REBUILD on authors' IGO file + peacesciencer controls. Existing key kept (interventionist = IGOmaxst). Full-sample rebuilt b=%.4g p=%.3g (pub -0.0094, ns). See IGOmedst row for the paper's headline term.",
          r4$IGOmaxst["b"], r4$IGOmaxst["p"]))
write_verify(108, "Boehmer, Gartzke & Nordstrom 2004", "TIE", "IGOmedst", "Table 2 Model 2-3 (structured IGOs, p.24)",
  -0.0854, 0.0429, 109965, "logit, MID onset, all dyads 1950-91, peace years + 3 splines; SE type not reported", TRUE, mism, "PARTIAL", unf, TRUE, M3,
  sprintf("PAPER'S KEY TERM (flag: differs from the key previously used). Full-sample rebuilt b=%.4g p=%.3g (pub -0.0854, p<.05). Same common sample as IGOmaxst row.",
          r3$IGOmedst["b"], r3$IGOmedst["p"]))
sink()
