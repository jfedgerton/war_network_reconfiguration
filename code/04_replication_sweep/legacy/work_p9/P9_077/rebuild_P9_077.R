# P9_077 Tschantret 2020, EJIR 26(4) "The Democratic Peace's Hidden Costs: Democratic Breakdown and International Aggression"
#   (10.1177/1354066120916525). DEMOCRACY. REBUILD (protocol B; no archive).
# Article (staged PDF pp. 9-12, Table 1): politically relevant directed dyad-years (1946-2010), logit, clustered SEs, cubic time since last
#   MID; DV = MID initiation by side A; Post-Democracy = A is an autocracy whose regime emerged from democratic breakdown (GWF);
#   Other Autocracy = other autocratic A (democracies the base); controls CINC A, CINC B, A's share of dyadic capabilities, global S-score,
#   civil war A and B, major power A and B, regime duration, t, t^2, t^3.
# HEADLINE (fixed before estimating): Table 1 Model 1, Post-Democracy 0.613 (0.149) ***, N 118,141.
# Rebuild inputs: {peacesciencer} COW directed dyad-years 1946-2010 (GML MID onset with A initiator), CINC, majors, contiguity for the PR rule,
#   Haege S (srsvaa: unweighted, valued alliances, absolute distances) via add_fpsim; UCDP/PRIO intrastate conflict (ucdpongoing, G-W codes
#   mapped to COW); GWF 1.2 All Political Regimes (work_p9/gwf): autocracy = non-missing gwf_regimetype, post-democracy = gwf_prior == democracy,
#   regime duration = gwf_duration of A. Clustered on directed dyad (cluster unit not stated).
ROOT <- Sys.getenv("SWEEP_DATA"); OUT <- Sys.getenv("SWEEP_OUT"); CODE <- Sys.getenv("SWEEP_CODE")
stopifnot(nzchar(ROOT), nzchar(OUT), nzchar(CODE))
source(file.path(CODE, "work_p8/common_b3.R")); source(file.path(CODE, "work_p9/rebuild_common_p9r.R"))
set.seed(123)
W <- file.path(OUT, "work_p9/P9_077"); dir.create(W, recursive = TRUE, showWarnings = FALSE)
sink(file.path(W, "log_P9_077.txt"), split = TRUE)
d <- ps_dyads(1946:2010, directed = TRUE, trade = FALSE, fpsim = TRUE)
d <- d %>% mutate(ddyad = 1000 * ccode1 + ccode2, init = as.integer(gmlmidonset == 1 & init1 %in% 1))
d <- py_calc(d, "ddyad", "init")
d <- d[pr_flag(d) == 1, ]
g <- as.data.frame(haven::zap_labels(haven::read_dta(file.path(ROOT, "work_p9/gwf/GWF_AllPoliticalRegimes.dta")))) %>%
  transmute(ccode = as.integer(cowcode), year = as.integer(year), autoc = as.integer(!is.na(gwf_regimetype) & gwf_regimetype != ""),
            postdem = as.integer(autoc == 1 & gwf_prior == "democracy"), regdur = gwf_duration)
stopifnot(!anyDuplicated(g[, c("ccode", "year")]))
cw <- suppressMessages(create_stateyears(system = "gw", subset_years = 1946:2010) %>% add_ccode_to_gw() %>% add_ucdp_acd(type = "intrastate")) %>%
  as.data.frame() %>% filter(!is.na(ccode)) %>% group_by(ccode, year) %>% summarise(cw = max(ucdpongoing), .groups = "drop")
n0 <- nrow(d)
d <- d %>% left_join(g %>% rename_with(~paste0(.x, "A"), c(autoc, postdem, regdur)), by = c("ccode1" = "ccode", "year")) %>%
  left_join(cw %>% rename(cwA = cw), by = c("ccode1" = "ccode", "year")) %>% left_join(cw %>% rename(cwB = cw), by = c("ccode2" = "ccode", "year"))
stopifnot(nrow(d) == n0)
d <- d %>% mutate(otherautA = as.integer(autocA == 1 & postdemA == 0), shareA = cinc1 / (cinc1 + cinc2),
  cwA = coalesce(cwA, 0L), cwB = coalesce(cwB, 0L), py2 = py^2, py3 = py^3)
rhs0 <- c("postdemA", "otherautA", "cinc1", "cinc2", "shareA", "srsvaa", "cwA", "cwB", "cowmaj1", "cowmaj2", "regdurA", "py", "py2", "py3")
vars <- c("init", rhs0, "ddyad")
cat("PR directed dyad-years 1946-2010:", nrow(d), " missing:\n"); print(sapply(d[, vars], function(x) sum(is.na(x))))
d0 <- d[complete.cases(d[, vars]), ]
cat("complete-case N:", nrow(d0), "(published 118,141)  initiations", sum(d0$init), "\n")
TERM <- "postdemA"
fitter <- function(r, data) fit_bin("init", r, data, "logit", "~ddyad", terms = TERM)
full <- fitter(rhs0, d0)
cat("\nREBUILT vs published: postdemA b", full[[TERM]][["b"]], "(pub 0.613) se", full[[TERM]][["se"]], "(pub 0.149) N", attr(full, "n"), "(pub 118,141)\n")
ds <- common_sample(d, vars, "ccode1", "ccode2", "year"); rm(d, d0); invisible(gc())
M <- run_six(rhs0, ds, fitter)
if (Sys.getenv("SWEEP_SENS") == "1") {
  source(file.path(CODE, "hist/sens_helper.R"))
  sens_run(id = "P9_077", dat = ds, y = "init", tie = TERM, ctrl = setdiff(rhs0, TERM), fe = NULL,
           fit = function(dat, add) fitter(c(rhs0, add), dat)[[TERM]])
}
saveRDS(list(M = M, full = full), file.path(W, "M_P9_077.rds"))
meta <- list(pid = "P9_077", cite = "Tschantret 2020", journal = "EJIR", year = 2020, doi = "10.1177/1354066120916525", family = "DEMOCRACY")
nt <- sprintf(paste0("REBUILT (no archive): peacesciencer PR directed dyads 1946-2010, GML MID onset initiated by A; GWF 1.2 post-democracy (gwf_prior democracy) ",
  "and other autocracy for A, regime duration; CINC A/B and share; Haege S (srsvaa); UCDP intrastate conflict A/B; COW majors; cubic time; clustered on ",
  "directed dyad. Full-sample rebuilt b=%.6f se=%.6f N=%d vs published 0.613 (0.149) N 118,141. Common sample N=%d."),
  full[[TERM]][["b"]], full[[TERM]][["se"]], attr(full, "n"), nrow(ds))
write_p8(meta, TERM, "abstract: autocracies established through democratic breakdown are more belligerent; Table 1 Model 1",
  "Table 1 Model 1 (logit, MID initiation, PR directed dyads 1946-2010)", "logit, clustered SE, cubic time", "feglm logit, SE clustered on directed dyad, all six steps",
  0.613, 0.149, 118141, TRUE, TRUE, M, "Y (cubic time since MID)", nt)
dir.create(file.path(OUT, "results_p9"), recursive = TRUE, showWarnings = FALSE)
src <- file.path(OUT, "results_p8", "P9_077__postdemA.csv"); stopifnot(file.exists(src), file.rename(src, file.path(OUT, "results_p9", basename(src))))
sink()
