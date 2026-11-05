# X076 Gartzke 2000 ISQ "Preferences and the Democratic Peace" (10.1111/0020-8833.00155). REBUILD (recovery protocol B).
# No article archive (author site lists no replication file; Dataverse none). Paper PDF (author site) read: Table 1 Model II =
#   logit, robust SE, O&R (1998) politically relevant dyads 1950-85, N=17,910; DV MID dummy; IVs: lower & higher democracy (Polity III
#   dem-aut), lower economic growth (avg GDPpc change over previous 3 yrs), allies, capability ratio (raw), contiguity, lower trade/GDP (lag 1),
#   affinity (UNGA voting similarity, lag 1). Published: affinity -1.12755 (z -4.117, p<.001); lower democracy -0.02105 (p=.209).
# Rebuild inputs: Russett-Oneal's own PR-dyad data (TRIANGLE.tab, Dataverse doi:10.7910/DVN/OMJYE5; O&R coding of dispute1, democracy,
#   dependence, allies, capability ratio, contiguity) + Gartzke's own Affinity of Nations 2010 release (s3un4608i: S-based 3-category, 1964 interpolated;
#   the 2000 paper used an earlier Spearman-based version) + lower GDPpc growth from peacesciencer SDP GDP estimates.
# KEY (a priori): affinity (lag 1) -- "residuals from ... regression of affinity remain highly significant as a predictor of the absence of MIDs";
#   the paper's tie term (Model II uses raw affinity; Model IV residualized). CONTROL: lower democracy (weak-link joint democracy).
ROOT <- Sys.getenv("SWEEP_DATA"); OUT <- Sys.getenv("SWEEP_OUT"); CODE <- Sys.getenv("SWEEP_CODE")
stopifnot(nzchar(ROOT), nzchar(OUT), nzchar(CODE))
source(file.path(CODE, "work_p4/common_p4_b04.R"))
source(file.path(CODE, "work_p8/round2_rule.R"))   # rebuild_verdict(): round-2 rule with registry alpha_pub (2026-10-11)
suppressPackageStartupMessages({library(haven); library(peacesciencer)})
WO <- file.path(OUT, "work_p4/X076"); dir.create(WO, recursive = TRUE, showWarnings = FALSE)  # port: outputs
W <- file.path(ROOT, "work_p4/X076"); sink(file.path(WO, "log_X076.txt"), split = TRUE)
t <- read.delim(file.path(ROOT, "recA/or/TRIANGLE.tab"))
t <- t[t$year >= 1950 & t$year <= 1985, ]
t <- t %>% mutate(c1 = pmin(statea, stateb), c2 = pmax(statea, stateb), lowdem = pmin(demauta, demautb), highdem = pmax(demauta, demautb),
  smldep = pmin(dependa, dependb), caprat = exp(lcaprat2), contig = -noncontg)   # noncontg stored as (dummy-1): 0 = noncontig, -1 = contig
stopifnot(all(t$contig %in% c(0, 1)), anyDuplicated(t[, c("c1", "c2", "year")]) == 0)
cat("TRIANGLE 1950-85 rows:", nrow(t), " dispute1 mean:", mean(t$dispute1, na.rm = TRUE), "\n")
# affinity, lag 1
a <- zap_labels(read_dta(file.path(W, "affinity_01242010.dta"))) %>% as.data.frame(); print(names(a))
a <- a %>% transmute(c1 = pmin(ccodea, ccodeb), c2 = pmax(ccodea, ccodeb), year = year + 1, aff_l1 = s3un4608i) %>% filter(c1 != c2) %>%
  group_by(c1, c2, year) %>% summarise(aff_l1 = mean(aff_l1, na.rm = TRUE), .groups = "drop")
t <- t %>% left_join(a, by = c("c1", "c2", "year"))
# lower growth: mean annual % change in GDPpc over t-3..t-1 (SDP estimates), weak link
sy <- tryCatch(create_stateyears(system = "cow", subset_years = 1940:1990) %>% add_sdp_gdp(), error = function(e) { cat("add_sdp_gdp failed:", conditionMessage(e), "\n"); NULL })
stopifnot(!is.null(sy), "wbgdppc2011est" %in% names(sy))
sy <- sy %>% arrange(ccode, year) %>% group_by(ccode) %>%
  mutate(gpc = exp(wbgdppc2011est), gr = 100 * (gpc / lag(gpc) - 1), gr3 = (lag(gr, 1) + lag(gr, 2) + lag(gr, 3)) / 3) %>% ungroup() %>% select(ccode, year, gr3)
t <- t %>% left_join(sy %>% rename(ga = gr3), by = c("statea" = "ccode", "year")) %>% left_join(sy %>% rename(gb = gr3), by = c("stateb" = "ccode", "year")) %>%
  mutate(lowgrowth = pmin(ga, gb), dyadid = 1000 * c1 + c2) %>% as.data.frame()
f0 <- dispute1 ~ lowdem + highdem + lowgrowth + allies + caprat + contig + smldep + aff_l1
FIT <- ff_logit_cl("dyadid")   # published robust SEs (affinity .274, lowdem .0168) match dyad-clustered, not HC1 (.109/.007)
r0 <- FIT(f0, t[complete.cases(t[, c(all.vars(f0), "dyadid")]), ]); print(round(r0, 5)); n_orig <- attr(r0, "n"); cat("rebuilt N:", n_orig, "(published 17,910)\n")
rc <- ff_logit_hc1(f0, t); cat("HC1 (unclustered) SEs for reference:\n"); print(round(rc[c("aff_l1", "lowdem"), ], 5))
st <- function(term, bp, sigp) rebuild_verdict("X076", term, r0[term, 1], r0[term, 3], bp, sigp)   # round-2 rule (was: sign and same .05 significance)
s_aff <- st("aff_l1", -1.12755, TRUE); s_dem <- st("lowdem", -0.02105, FALSE); cat(s_aff, s_dem, "\n")
ds <- common_sample(t, f0, "dyadid", "statea", "stateb", "year")
M <- run_six(f0, ds, c("aff_l1", "lowdem"), FIT, WO)
# ---- sensitivity hook: runs only under run_sensitivity.R (SWEEP_SENS=1), never in the main sweep ----
if (Sys.getenv("SWEEP_SENS") == "1") {
  source(file.path(CODE, "hist/sens_helper.R"))
  sens_run(id = "X076", dat = ds, y = "dispute1", tie = "lowdem",
           ctrl = c("highdem", "lowgrowth", "allies", "caprat", "contig", "smldep", "aff_l1"), fe = NULL,
           fit = function(dat, add) { ct <- FIT(if (length(add)) update(f0, as.formula(paste(". ~ . +", paste(add, collapse = " + ")))) else f0, dat)
             c(b = unname(ct["lowdem", 1]), se = unname(ct["lowdem", 2]), p = unname(ct["lowdem", 3])) })
}
meta <- list(xid = "X076", short_cite = "Gartzke 2000", journal = "ISQ", year = 2000, doi = "10.1111/0020-8833.00155")
est <- "logit, SE clustered on dyad (published 'Robust S.E.' magnitudes match dyad clustering)"
nt <- paste0("REBUILT: O&R TRIANGLE PR dyads 1950-85 (their dispute1/democracy/dependence/allies/caprat/contiguity coding) + Gartzke Affinity (2010 release) s3un4608i lagged 1 yr ",
  "(S-based; paper used earlier Spearman-based affinity) + lower GDPpc growth from peacesciencer SDP estimates (paper: PWT). Table 1 Model II spec, no peace years (as published). ",
  "Rebuilt N ", n_orig, " vs published 17,910. HC1 SEs shown in log for reference (far smaller than published).")
write_p4(meta, "aff_l1", "TIE", "KEY", "\"residuals from ... regression of affinity remain highly significant as a predictor of the absence of MIDs\"",
  s_aff, -1.12755, r0["aff_l1", 1], est, n_orig, M, "N (no temporal/history controls)", paste(nt, "Key = UN affinity (lag 1), Model II raw-affinity version."))
write_p4(meta, "lowdem", "DEMOCRACY", "CONTROL", "\"significance for democracy is shown to be fragile\" (lower democracy in Model II)",
  s_dem, -0.02105, r0["lowdem", 1], est, n_orig, M, "N (no temporal/history controls)", paste(nt, "Lower (weak-link) democracy; published n.s. in Model II."))
sink()
