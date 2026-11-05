# P8_048 Ekmekci 2014 Peace and Conflict Studies 21(1) "Democratic vs. Capitalist Peace: A Test in the Developing World"
# (NSUWorks open access, nsuworks.nova.edu/pcs/vol21/iss1/5; PDF behind Cloudflare -> design read via WebFetch summary of the PDF).
# No replication archive found (Dataverse title/author search, WebSearch) -> REBUILD with peacesciencer (standard dyad-year MID logit).
# Paper design (Table I): undirected dyad-years 1951-2000; contiguous dyads (land or <150 miles water) with at least one developing state
#   (dyads of two developed states excluded; developed = US, Canada, UK, NL, BEL, LUX, FRA, SWZ, W.Germany, AUS(305), ITA, FIN, SWE, NOR,
#   DEN, ICE, JPN, AUL, NZ; microstates dropped); DV = Maoz dyadic MID occurrence (onset and continuation treated alike);
#   logit, robust SE clustered on dyad, peace years + 3 cubic splines (BKT). IVs lagged one year, weak link:
#   Development_low = lower real GDP per capita (Gleditsch 2002, PPP 1996$); Democracy_low = lower Polity (1-21 scale);
#   controls: ln capability ratio, alliance, major power in dyad, state age (younger state's years since independence). N = 12,175.
# Published (Table I): Development_low -0.0743 (robust SE .0238, p<.01); Democracy_low -0.0096 (.0113, ns).
# KEY TERMS (a priori, abstract/title "Democratic vs. Capitalist Peace"): development_low (CAPITALISM) and democracy_low (DEMOCRACY).
# Rebuild substitutions: GML dyadic MID ongoing (gmlmidongoing) for Maoz DYMID; GDP pc = exp(SDP wbgdppc2011est)/1000 (thousand 2011$)
#   for Gleditsch PPP 1996$ (scale differs -> compare sign/significance only); polity2+11 from ccode_democracy (known add_democracy() bug);
#   COW alliance any type; state age = year - first COW system-entry year; cubic peace years (py, py^2, py^3) instead of BKT splines.
suppressMessages({library(peacesciencer)})
CODE <- Sys.getenv("SWEEP_CODE"); stopifnot(nzchar(CODE))
source(file.path(CODE, "work_p8/common_b3.R"))
WD <- file.path(ROOT, "work_p8/P8_048"); W <- file.path(OUT, "work_p8/P8_048"); dir.create(W, recursive = TRUE, showWarnings = FALSE); sink(file.path(W, "log_P8_048.txt"), split = TRUE)
dev <- c(2, 20, 200, 210, 211, 212, 220, 225, 260, 305, 325, 375, 380, 385, 390, 395, 740, 900, 920)
micro <- c(221, 223, 232, 331)   # Monaco, Liechtenstein, Andorra, San Marino
# DV + peace years (years since last MID year in the dyad, from 1816)
nd <- create_dyadyears(directed = FALSE, subset_years = 1816:2000) %>% add_gml_mids() %>%
  mutate(mid = as.numeric(gmlmidongoing == 1 | gmlmidonset == 1)) %>% arrange(ccode1, ccode2, year) %>% group_by(ccode1, ccode2) %>%
  mutate(spell = cumsum(lag(mid, default = 0) == 1)) %>% group_by(ccode1, ccode2, spell) %>% mutate(py = row_number() - 1) %>% ungroup() %>%
  filter(year >= 1951) %>% select(ccode1, ccode2, year, mid, py)
pol <- ccode_democracy %>% select(ccode, year, polity2) %>% as.data.frame(); stopifnot(!anyDuplicated(pol[, 1:2]))
entry <- cow_states %>% group_by(ccode) %>% summarise(entry = min(styear), .groups = "drop")
ctl <- create_dyadyears(directed = FALSE, subset_years = 1950:2000) %>% add_contiguity() %>% add_nmc() %>% add_cow_majors() %>%
  {suppressWarnings(add_sdp_gdp(.))} %>%
  left_join(pol %>% rename(polity21 = polity2), by = c("ccode1" = "ccode", "year")) %>%
  left_join(pol %>% rename(polity22 = polity2), by = c("ccode2" = "ccode", "year")) %>%
  left_join(peacesciencer::cow_alliance, by = c("ccode1", "ccode2", "year")) %>%
  left_join(entry %>% rename(entry1 = entry), by = c("ccode1" = "ccode")) %>% left_join(entry %>% rename(entry2 = entry), by = c("ccode2" = "ccode")) %>%
  mutate(development_low_t = pmin(exp(wbgdppc2011est1), exp(wbgdppc2011est2)) / 1000,
         democracy_low_t = pmin(polity21, polity22) + 11,
         lncapratio = log(pmax(cinc1, cinc2) / pmin(cinc1, cinc2)),
         alliance = as.numeric(coalesce(cow_defense, 0) + coalesce(cow_neutral, 0) + coalesce(cow_entente, 0) > 0),
         majpow = as.numeric(cowmaj1 == 1 | cowmaj2 == 1), stateage = year - pmax(entry1, entry2),
         contig = as.numeric(conttype %in% 1:4)) %>%
  arrange(ccode1, ccode2, year) %>% group_by(ccode1, ccode2) %>%
  mutate(development_low = ifelse(lag(year) == year - 1, lag(development_low_t), NA),
         democracy_low = ifelse(lag(year) == year - 1, lag(democracy_low_t), NA)) %>% ungroup() %>%
  select(ccode1, ccode2, year, contig, development_low, democracy_low, lncapratio, alliance, majpow, stateage)
d <- ctl %>% inner_join(nd, by = c("ccode1", "ccode2", "year")) %>% filter(year >= 1951, contig == 1) %>%
  filter(!(ccode1 %in% dev & ccode2 %in% dev), !(ccode1 %in% micro), !(ccode2 %in% micro)) %>%
  mutate(py2 = py^2, py3 = py^3, dyad = ccode1 * 1000 + ccode2) %>% as.data.frame()
rm(nd, ctl); invisible(gc())
d <- d[is.finite(d$lncapratio) | is.na(d$lncapratio), ]
rhs0 <- c("development_low", "democracy_low", "lncapratio", "alliance", "majpow", "stateage", "py", "py2", "py3")
KEY <- c("development_low", "democracy_low")
fitter <- function(rhs, data) fit_bin("mid", rhs, data, "logit", "~dyad", terms = KEY)
vars <- c("mid", rhs0, "dyad")
full <- d[complete.cases(d[, vars]), ]
cat("rebuilt complete-case N =", nrow(full), "(published 12,175); MID dyad-years:", sum(full$mid), "\n")
r_full <- fitter(rhs0, full); n_orig <- attr(r_full, "n")
ds <- common_sample(d, vars, "ccode1", "ccode2", "year")
M <- run_six(rhs0, ds, fitter)
# ---- sensitivity hook: runs only under run_sensitivity.R (SWEEP_SENS=1), never in the main sweep ----
if (Sys.getenv("SWEEP_SENS") == "1") {
  source(file.path(CODE, "hist/sens_helper.R"))
  sens_run(id = "P8_048", dat = ds, y = "mid", tie = "development_low", ctrl = setdiff(rhs0, "development_low"), fe = NULL,
           fit = function(dat, add) fitter(c(rhs0, add), dat)[["development_low"]])
}
saveRDS(list(M = M, r_full = r_full, n_orig = n_orig, n_used = nrow(ds)), file.path(W, "M_P8_048.rds"))
meta <- list(pid = "P8_048", cite = "Ekmekci 2014", journal = "PCS", year = 2014, doi = "", family = "CAPITALISM")
nt <- sprintf(paste0("REBUILT (no archive). GML MID occurrence (onset+ongoing) for Maoz DYMID; dev_low = min exp(SDP wbgdppc2011est)/1000 (2011$) for Gleditsch PPP 1996$ ",
  "(scale differs); dem_low = min polity2+11 (ccode_democracy); COW alliance; state age from COW entry; cubic peace years for BKT splines; contiguity conttype 1-4; ",
  "dyads of two developed states and microstates excluded. Rebuilt full N=%d (pub 12,175); full-sample b(dev_low)=%.4f p=%.3g, b(dem_low)=%.4f p=%.3g. Common 1822-2011 sample N=%d."),
  n_orig, r_full$development_low[["b"]], r_full$development_low[["p"]], r_full$democracy_low[["b"]], r_full$democracy_low[["p"]], nrow(ds))
write_p8(modifyList(meta, list(family = "CAPITALISM")), "development_low", "title/abstract: development -> peace ('Democratic vs. Capitalist Peace')",
  "Table I (logit, contiguous developing-world dyads 1951-2000)", "logit, robust SE clustered on dyad, peace years + 3 splines",
  "logit (fixest), SE clustered on dyad, all six steps", -0.0743, 0.0238, 12175, TRUE, TRUE, M, "N (peace years + splines only)", nt)
write_p8(modifyList(meta, list(family = "DEMOCRACY")), "democracy_low", "title/abstract: democratic peace term tested against development",
  "Table I (logit, contiguous developing-world dyads 1951-2000)", "logit, robust SE clustered on dyad, peace years + 3 splines",
  "logit (fixest), SE clustered on dyad, all six steps", -0.0096, 0.0113, 12175, FALSE, TRUE, M, "N (peace years + splines only)", nt)
sink()
