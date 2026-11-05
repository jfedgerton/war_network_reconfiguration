# VERIFY X091: Mousseau, Hegre & Oneal 2003 EJIR -- rebuild matched to the PDF (Table 2 Model 2, the preferred model:
# "model 2 provides a significantly better account", p.297). Changes vs rebuild_X091.R, each justified by the PDF:
#  (a) sample = politically relevant dyads: contiguous or >=1 major power (p.289, p.294 "contiguous states and pairs that
#      include at least one major power"); pub N = 39,584.                                         [was: all dyads]
#  (b) years: 1885-1992, "All but the first year of the First and Second World Wars are omitted ... immediate postwar
#      years, 1919-20 and 1946-9, are excluded" (p.288-289) -> drop 1915-1920 and 1940-1949.   [was: all years]
#  (c) temporal dependence: "Brevity of Peace ... equals 2^(-y/5)" + "Past Dispute" (spell initiated by a fatal dispute
#      rather than by independence) (p.294).                                                    [was: py + py^2 + py^3]
#  (d) DemocracyL, DependencyL, DevelopmentL "centered each constituent term on its median ... divided by the standard
#      deviation" (p.291-293); interaction = product of the centered/scaled terms.              [was: raw terms]
#  (e) Major Power = 1 if dyad includes at least one great power (p.294).                       [was: both-minor dummy]
# Kept (paper silent or data unavailable): IVs at t-1; GML dyadic fatal onsets (paper: Maoz 1999 DyMID);
#  polity2 (paper: Polity III); COW trade current $ / SDP-WB 2011$ GDP (paper: Gleditsch 2002 trade & GDP);
#  direct contiguity only (paper also counts colonial contiguity); dyad-clustered SE (paper: Huber/White, dyadic).
suppressMessages(library(peacesciencer))
ROOT <- Sys.getenv("SWEEP_DATA"); OUT <- Sys.getenv("SWEEP_OUT"); CODE <- Sys.getenv("SWEEP_CODE")
stopifnot(nzchar(ROOT), nzchar(OUT), nzchar(CODE))
source(file.path(CODE, "work_p4/b3/common_b3.R"))
source(file.path(CODE, "vb3/vhelp.R"))
WO <- file.path(OUT, "work_p4/X091/verify"); dir.create(WO, recursive = TRUE, showWarnings = FALSE)  # port: outputs
W <- file.path(ROOT, "work_p4/X091/verify")  # port: outputs go to WO; nothing is read from W
sink(file.path(OUT, "work_p4/X091/log_verify_X091.txt"), split = TRUE)

fd <- peacesciencer::gml_dirdisp %>% filter(midonset == 1, ccode1 < ccode2) %>%
  mutate(fatal = as.numeric(pmax(fatality1, fatality2) > 0)) %>%
  group_by(ccode1, ccode2, year) %>% summarise(fatalonset = max(fatal), .groups = "drop") %>% filter(fatalonset == 1)
nd <- create_dyadyears(directed = FALSE, subset_years = 1816:1992) %>% select(ccode1, ccode2, year) %>%
  left_join(fd, by = c("ccode1","ccode2","year")) %>% mutate(fatalonset = coalesce(fatalonset, 0)) %>%
  arrange(ccode1, ccode2, year) %>% group_by(ccode1, ccode2) %>%
  mutate(spell = cumsum(dplyr::lag(fatalonset, default = 0) == 1)) %>% group_by(ccode1, ccode2, spell) %>%
  mutate(py = row_number() - 1) %>% ungroup() %>%
  mutate(brevity = 2^(-py / 5), pastdisp = as.numeric(spell > 0)) %>%               # (c)
  filter(year >= 1885, !(year %in% c(1915:1920, 1940:1949))) %>%                      # (b)
  select(ccode1, ccode2, year, fatalonset, brevity, pastdisp)
stopifnot(all(nd$brevity > 0 & nd$brevity <= 1))

tr_f <- system.file("extdata", "cow_trade_ddy.rds", package = "peacesciencer"); stopifnot(nzchar(tr_f), file.exists(tr_f))  # port
tr <- readRDS(tr_f) %>%
  mutate(flow1 = ifelse(flow1 < 0, NA, flow1), flow2 = ifelse(flow2 < 0, NA, flow2),
         lo = pmin(ccode1, ccode2), hi = pmax(ccode1, ccode2)) %>%
  group_by(lo, hi, year) %>% summarise(trade = mean(flow1 + flow2), .groups = "drop")
iv <- create_dyadyears(directed = FALSE, subset_years = 1884:1991) %>%
  add_contiguity() %>% add_capital_distance() %>% add_democracy() %>% add_nmc() %>% add_cow_majors() %>%
  {suppressWarnings(add_sdp_gdp(.))}
iv <- iv %>% left_join(peacesciencer::cow_alliance, by = c("ccode1","ccode2","year")) %>%
  left_join(tr, by = c("ccode1" = "lo", "ccode2" = "hi", "year" = "year")) %>%
  mutate(allies = as.numeric(coalesce(cow_defense, 0) + coalesce(cow_neutral, 0) + coalesce(cow_entente, 0) > 0),
         demlo_raw = pmin(polity21, polity22),
         deplo_raw = pmin(trade * 1e6 / exp(wbgdp2011est1), trade * 1e6 / exp(wbgdp2011est2)),
         gdplo_raw = pmin(wbgdppc2011est1, wbgdppc2011est2),
         lncaprat = log(pmax(cinc1, cinc2) / pmin(cinc1, cinc2)),
         contig = as.numeric(conttype %in% 1:4), lndist = log(capdist + 1),
         majpow = as.numeric(cowmaj1 == 1 | cowmaj2 == 1),                               # (e)
         year = year + 1) %>%
  select(ccode1, ccode2, year, allies, demlo_raw, deplo_raw, gdplo_raw, lncaprat, contig, lndist, majpow)
d <- nd %>% inner_join(iv, by = c("ccode1","ccode2","year")) %>%
  filter(contig == 1 | majpow == 1) %>%                                                  # (a)
  mutate(dyad = ccode1 * 1000 + ccode2) %>% as.data.frame()
d <- d[is.finite(d$lncaprat) | is.na(d$lncaprat), ]
rm(nd, iv, tr, fd); invisible(gc())
rawv <- c("fatalonset","demlo_raw","deplo_raw","gdplo_raw","lncaprat","allies","contig","lndist","majpow","brevity","pastdisp")
d <- d[complete.cases(d[, rawv]), ]
# (d) center on median, divide by SD of the uncentered variable, on the estimation sample
zs <- function(x) (x - median(x)) / sd(x)
d$DemocracyL <- zs(d$demlo_raw); d$DependencyL <- zs(d$deplo_raw); d$DevelopmentL <- zs(d$gdplo_raw)
d$DemXDev <- d$DemocracyL * d$DevelopmentL
cat("PR sample N =", nrow(d), "(pub 39,584); fatal onsets:", sum(d$fatalonset), "\n")
f0 <- fatalonset ~ DemocracyL + DependencyL + DevelopmentL + DemXDev + lncaprat + majpow + allies + contig + lndist + brevity + pastdisp
terms <- c("DependencyL", "DemXDev", "DemocracyL")
rep <- fit_keys(f0, d, terms, "~dyad", "logit", file.path(WO, "cache_repro.rds"))
# Table 1 (no interaction) as a sanity check on the controls
t1 <- fit_keys(update(f0, . ~ . - DemXDev), d, c("DependencyL","DemocracyL"), "~dyad", "logit", file.path(WO, "cache_t1.rds"))
ds <- common_sample(d, f0, "dyad", "ccode1", "ccode2", "year"); n_pr <- nrow(d); rm(d); invisible(gc())
M <- run_six(f0, ds, terms, "~dyad", "logit", WO)
R <- lapply(setNames(terms, terms), function(t) lapply(M, function(z) z$k[[t]]))

# ---- sensitivity hook: runs only under run_sensitivity.R (SWEEP_SENS=1), never in the main sweep ----
if (Sys.getenv("SWEEP_SENS") == "1") {
  source(file.path(CODE, "hist/sens_helper.R"))
  sens_run(id = "X091", dat = ds, y = "fatalonset", tie = "DependencyL",
           ctrl = c("DemocracyL", "DevelopmentL", "DemXDev", "lncaprat", "majpow", "allies", "contig", "lndist", "brevity", "pastdisp"), fe = NULL,
           fit = function(dat, add) { r <- fit_keys(if (length(add)) update(f0, as.formula(paste(". ~ . +", paste(add, collapse = " + ")))) else f0, dat, "DependencyL", "~dyad", "logit")$k$DependencyL
             c(b = unname(r[["b"]]), se = unname(r[["se"]]), p = unname(r[["p"]])) })
}


mism <- "sample all dyads vs politically relevant (contiguous or major power);war years 1915-20 and 1940-49 not dropped;peace years cubic vs brevity-of-peace 2^(-y/5) + past dispute;raw terms vs median-centered/SD-scaled DemocracyL/DependencyL/DevelopmentL;both-minor dummy vs major-power dummy"
unf <- "trade/GDP data vintage (COW current $ / SDP 2011$ vs Gleditsch 2002);Polity III vs polity2;Maoz DyMID vs GML fatal onsets;colonial contiguity not in PR definition;GDPpc source (SDP vs Maddison/Summers-Heston);lag of IVs not stated (kept t-1)"
nt <- sprintf("REBUILD matched to PDF (peacesciencer). PR N=%d vs pub 39,584; full rebuilt b DependencyL=%.3f (p %.3g), DemXDev=%.3f (p %.3g), DemocracyL=%.3f (p %.3g). Table 1 analogue DependencyL=%.3f (pub -.38), DemocracyL=%.3f (pub -.35).",
  n_pr, rep$k$DependencyL["b"], rep$k$DependencyL["p"], rep$k$DemXDev["b"], rep$k$DemXDev["p"], rep$k$DemocracyL["b"], rep$k$DemocracyL["p"],
  t1$k$DependencyL["b"], t1$k$DemocracyL["b"])
hdr <- "Table 2 Model 2 (PR dyads 1885-1992, Democracy x Development only)"
est <- "logit, Huber/White SE (dyadic); brevity of peace + past dispute"
write_verify("X091", "Mousseau, Hegre & Oneal 2003", "TIE", "deplo", hdr, -0.34, 0.14, 39584, est, mism, "PARTIAL", unf, TRUE,
  R$DependencyL, M$m0$n, paste("deplo = DependencyL (centered/scaled lower trade/GDP).", nt), rebuilt = TRUE)
write_verify("X091", "Mousseau, Hegre & Oneal 2003", "DEMOCRACY", "demlo_x_gdplo", hdr, -0.33, 0.09, 39584, est, mism, "PARTIAL", unf, TRUE,
  R$DemXDev, M$m0$n, paste("demlo_x_gdplo = DemocracyL x DevelopmentL (centered/scaled). bD not a placebo.", nt), rebuilt = TRUE)
sink()
