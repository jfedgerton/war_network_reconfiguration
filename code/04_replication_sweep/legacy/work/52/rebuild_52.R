# REBUILD 52: Benson 2004 JPR "Dyadic Hostility and the Ties That Bind" (10.1177/0022343304047431)
# Recovery protocol B. Authors' posted file (PRIO) has only the 7 MDS tie IVs + ccode1/ccode2/year/dyadid.
# DV (HOSTLEVD) and controls (EUGene) are rebuilt with peacesciencer (GML MIDs, COW NMC/contiguity/alliance/Polity).
# Published design (paper p.662-668, do-file file45999_benson.do): mlogit HOSTLEVD (0 none,1 threat,2 show,3 use,4 war)
#   on DyadEconDissim IntOrderEconDist DyadicSecDissim IntOrderSecDist Democ Capratio Geodistance Contiguity Alliance Peaceyears,
#   robust cluster(dyadid); 1951-85; n = 189,214 (rows in both orderings of each pair, as in the posted file).
# KEY TERM (fixed before estimation): Model 1, DyadEconDissim (intra-dyadic economic dissimilarity, trade-based tie),
#   contrast 0 -> use of force (modal hostile outcome, 452/734 events). Published factor change 1.295 (z=3.132) -> b_pub = ln(1.295).
suppressMessages({library(dplyr); library(haven); library(peacesciencer)})
ROOT <- Sys.getenv("SWEEP_DATA"); OUT <- Sys.getenv("SWEEP_OUT"); CODE <- Sys.getenv("SWEEP_CODE")
stopifnot(nzchar(ROOT), nzchar(OUT), nzchar(CODE))
source(file.path(CODE, "work/rerun_helpers_b118.R"))
source(file.path(CODE, "recB/mlogit_helpers.R"))
source(file.path(CODE, "work_p8/round2_rule.R"))   # rebuild_verdict(): round-2 rule with registry alpha_pub (2026-10-11)
WO <- file.path(OUT, "work/52"); dir.create(WO, recursive = TRUE, showWarnings = FALSE)  # port: outputs
W <- file.path(ROOT, "work/52"); logf <- file.path(WO, "log_52.txt"); cat("", file = logf)
lg <- function(...) { cat(..., "\n"); cat(..., "\n", file = logf, append = TRUE) }
set.seed(123)

# ---- 1. authors' tie variables ----
b <- read_dta(file.path(W, "BensonUnravelingTies.dta")) %>% zap_labels() %>% as.data.frame()
b <- b[!is.na(b$year) & !is.na(b$ccode1) & !is.na(b$ccode2), ]
stopifnot(anyDuplicated(b[, c("ccode1","ccode2","year")]) == 0)
lg("study 52 tie file rows:", nrow(b))

# ---- 2. rebuild DV: GML MID onsets with dyadic hostility level, nondirected, 1816-1991 (for peace years) ----
nd <- create_dyadyears(directed = FALSE, subset_years = 1816:1991) %>% add_gml_mids(keep = c("hostlev"))
nd <- nd %>% mutate(hostd = ifelse(gmlmidonset == 1 & !is.na(hostlev), pmax(hostlev - 1, 0), 0))  # MID hostlev 2..5 -> 1..4
stopifnot(all(nd$hostd %in% 0:4))
# peace years: years since last onset of any hostile action (counter starts at 0 in dyad's first system year)
nd <- nd %>% arrange(ccode1, ccode2, year) %>% group_by(ccode1, ccode2) %>%
  mutate(spell = cumsum(lag(gmlmidonset, default = 0) == 1)) %>% group_by(ccode1, ccode2, spell) %>%
  mutate(peaceyrs = row_number() - 1) %>% ungroup() %>% select(ccode1, ccode2, year, hostd, peaceyrs)
nd <- nd %>% filter(year >= 1950)

# ---- 3. controls ----
ctl <- create_dyadyears(directed = FALSE, subset_years = 1950:1991) %>%
  add_contiguity() %>% add_capital_distance() %>% add_democracy() %>% add_nmc()
ca <- peacesciencer::cow_alliance  # deprecated in package but present; COW alliance v4.1 dyad-years
ctl <- ctl %>% left_join(ca, by = c("ccode1","ccode2","year")) %>%
  mutate(across(c(cow_defense, cow_neutral, cow_nonagg, cow_entente), ~coalesce(.x, 0)),
         Democ = pmin(polity21, polity22),
         Capratio = pmax(cinc1, cinc2) / pmin(cinc1, cinc2),
         Geodistance = capdist,
         Contiguity = ifelse(conttype == 0, 6, conttype),                       # 1 land ... 5 400mi water, 6 none
         Alliance = case_when(cow_defense == 1 ~ 1, cow_neutral == 1 ~ 2, cow_entente == 1 ~ 3, TRUE ~ 4)) %>%
  select(ccode1, ccode2, year, Democ, Capratio, Geodistance, Contiguity, Alliance)
ctl <- ctl %>% left_join(nd, by = c("ccode1","ccode2","year"))
stopifnot(anyDuplicated(ctl[, c("ccode1","ccode2","year")]) == 0)

# merge onto both orderings of Benson rows
b$lo <- pmin(b$ccode1, b$ccode2); b$hi <- pmax(b$ccode1, b$ccode2)
d <- b %>% left_join(ctl, by = c("lo" = "ccode1", "hi" = "ccode2", "year" = "year"))
stopifnot(nrow(d) == nrow(b))
d <- d %>% filter(year >= 1951, year <= 1985)
vars <- c("hostd","DyadEconDissim","IntOrderEconDist","DyadicSecDissim","IntOrderSecDist",
          "Democ","Capratio","Geodistance","Contiguity","Alliance","peaceyrs")
d0 <- d[complete.cases(d[, vars]) & is.finite(d$Capratio), ]
lg("rebuilt 1951-85 complete-case n =", nrow(d0), "(published 189,214)")
lg("DV distribution:", paste(names(table(d0$hostd)), table(d0$hostd), collapse = "; "), "(pub: 38 threat, 202 show, 452 use, 42 war; counts include both orderings)")

rhs <- c("DyadEconDissim","IntOrderEconDist","DyadicSecDissim","IntOrderSecDist","Democ","Capratio","Geodistance","Contiguity","Alliance","peaceyrs")
m_pub <- mlogit_fit(d0$hostd, d0[, rhs], cluster = d0$dyadid)
lg("\n==== Reproduction attempt (Model 1, cluster = authors' dyadid) ====")
capture.output(print(round(m_pub$coeftable, 4)), file = logf, append = TRUE); print(round(m_pub$coeftable, 4))
st <- mlogit_stat(m_pub, "DyadEconDissim", 3)
b_pub <- log(1.295)
lg(sprintf("DyadEconDissim 0->use: b_repro=%.4f (se %.4f, p %.4g) vs b_pub=%.4f (z=3.132)", st["b"], st["se"], st["p"], b_pub))
rel <- abs(st["b"] - b_pub) / abs(b_pub)
status <- rebuild_verdict("52", "DyadEconDissim", unname(st["b"]), unname(st["p"]), b_pub, TRUE)   # round-2 rule (was: sign and p < .05)
lg("rel diff:", round(rel, 3), "->", status, "(sign & significance vs published)")

# ---- 4. history battery ----
d0 <- add_history(d0, "ccode1", "ccode2", "year"); d0$year_h <- d0$year
db <- battery_sample(d0, c(vars, "dyadid"))
fitter <- function(r, data) mlogit_fit(data$hostd, data[, r], cluster = data$dyadid)
stat <- function(m) mlogit_stat(m, "DyadEconDissim", 3)
R <- run_battery(db, rhs, fitter, stat, logf)
write_result(52, "Benson 2004", "REBUILT", b_pub, unname(st["b"]), "multinomial logit (HOSTLEVD 0-4), cluster dyadid; key contrast 0->use of force",
  "DyadEconDissim (0->use of force)", 189214, nrow(db), R, "Y",
  paste0(status, ": tie IVs are the authors' (PRIO file); DV/controls rebuilt: GML MID onset hostlev (2-5 -> 1-4), ongoing years = 0; peace years from GML onsets since 1816; ",
         "Democ=min polity2; Capratio=CINC high/low; Geodistance=capdist; Contiguity 1-5 (6=none); Alliance COW 1 def/2 neut/3 ent/4 none; ",
         "rows in both orderings as in posted file (n_pub 189,214 implies this); cluster = authors' dyadid (=1000+ccode1+ccode2, non-unique across pairs). ",
         "Key fixed a priori: Model 1 DyadEconDissim, 0->use contrast; published OR 1.295 (z 3.13). Rel diff ", round(rel, 2), "."))
