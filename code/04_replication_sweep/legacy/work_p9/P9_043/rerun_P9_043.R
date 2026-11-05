# P9_043 Peterson & Thies 2012 BJPS 42(4) "Beyond Ricardo: The Link between Intra-Industry Trade and Peace" (10.1017/S0007123412000129)
# Data + do-file: author's replication zip (timothypeterson.org/research, "BJPS-beyond-ricardo-replication.zip"): intra_cut_1_5.dta
#   (staged as peterson_data.dta), intra_models_1_5.do. No log.
# Undirected dyads 1962-2001 (dyadid = ccode1*1000 + ccode2, ccode1 < ccode2). DV ffatalmid = f.fataltest (fatal MID in t+1: a LEAD).
# Do-file "Final models", model 3 ("our noint"): logit ffatalmid lowerdev lowerdep intra_index bothaut bothdem missingpol onemp jointmp lndist
#   atopally lncapratio fpyr fpy2 fpy3 if fcwongo==0, cluster(dyadid). Stub intra_index resolves uniquely to intra_index_final (asserted).
# HEADLINE: model 3; KEY TERM intra_index_final (horizontal intra-industry trade index). Published b/SE/N NOT AVAILABLE here (BJPS closed
#   access; archive has no log) -> pub values NA, verdict NO_PUBLISHED pending a table lookup.
# Estimator: fixest feglm logit, vcov ~dyadid, all six steps. History/benchmarks attached at outcome_year = year + 1 (lead DV, as P9_008).
CODE <- Sys.getenv("SWEEP_CODE"); stopifnot(nzchar(CODE))
source(file.path(CODE, "work_p8/common_b3.R"))
WD <- file.path(ROOT, "work_p9/P9_043"); W <- file.path(OUT, "work_p9/P9_043"); dir.create(W, recursive = TRUE, showWarnings = FALSE)
sink(file.path(W, "log_P9_043.txt"), split = TRUE)
TERM <- "intra_index_final"
rhs0 <- c("lowerdev", "lowerdep", TERM, "bothaut", "bothdem", "missingpol", "onemp", "jointmp", "lndist", "atopally", "lncapratio",
          "fpyr", "fpy2", "fpy3")
d <- as.data.frame(haven::zap_labels(haven::read_dta(file.path(WD, "peterson_data.dta"))))
stopifnot(nrow(d) == 511445, !anyDuplicated(d[, c("dyadid", "year")]))
m <- grep("^intra_index", names(d), value = TRUE); cat("intra_index ->", m, "\n"); stopifnot(identical(m, TERM))
d$ccode1 <- d$dyadid %/% 1000; d$ccode2 <- d$dyadid %% 1000
stopifnot(all(d$ccode1 < d$ccode2), mean(d$ccode1 %in% peacesciencer::cow_states$ccode) > .99)
d <- d[d$fcwongo %in% 0, ]                                         # if fcwongo==0
cat("year range:", range(d$year), "\n")
# lead diagnostic: after ffatalmid[t] == 1 the forward peace-years counter at t+1 (fpyr[t+1]) should be 0
o <- d[order(d$dyadid, d$year), ]; no <- nrow(o)
cs <- o$dyadid[-1] == o$dyadid[-no] & o$year[-1] == o$year[-no] + 1
ev <- which(cs & o$ffatalmid[-no] %in% 1)
cat("P(fpyr[t+1] == 0 | ffatalmid[t] == 1):", mean(o$fpyr[ev + 1] == 0, na.rm = TRUE), " on", length(ev), "events\n"); rm(o, cs, ev, no)
vars <- c("ffatalmid", rhs0, "dyadid")
fitter <- function(r, data) fit_bin("ffatalmid", r, data, "logit", "~dyadid", terms = TERM)
d0 <- d[complete.cases(d[, vars]), ]
full <- fitter(rhs0, d0)
cat("\nFULL-SAMPLE (no published value available):", TERM, "b", full[[TERM]][["b"]], "se", full[[TERM]][["se"]], "N", attr(full, "n"), "\n")
d$outcome_year <- d$year + 1                                       # DV is fatal MID in t+1
ds <- common_sample(d, vars, "ccode1", "ccode2", "outcome_year")
M <- run_six(rhs0, ds, fitter)
saveRDS(list(M = M, full = full, n_used = nrow(ds)), file.path(W, "M_P9_043.rds"))

# ---- sensitivity hook: runs only under run_sensitivity.R (SWEEP_SENS=1), never in the main sweep ----
if (Sys.getenv("SWEEP_SENS") == "1") {
  source(file.path(CODE, "hist/sens_helper.R"))
  sens_run(id = "P9_043", dat = ds, y = "ffatalmid", tie = TERM, ctrl = setdiff(rhs0, TERM), fe = NULL,
           fit = function(dat, add) fitter(c(rhs0, add), dat)[[TERM]])
}

meta <- list(pid = "P9_043", cite = "Peterson & Thies 2012", journal = "BJPS", year = 2012, doi = "10.1017/S0007123412000129", family = "TIE")
write_p8(meta, TERM, "abstract: 'horizontal intra-industry trade is associated with reduced conflict propensity within dyads'",
  "Model 3 (no interactions; do-file 'our noint')", "logit, SE clustered on dyad, peace years + square + cube (forward)",
  "feglm logit, vcov clustered dyadid, all six steps", -4.778, 1.830, 191175, TRUE, FALSE, M, "N (peace years + square + cube only)",
  sprintf(paste0("Undirected dyads, DV fatal MID in t+1, fcwongo==0 (authors' data + do-file). Published = article Table 2 MODEL 2 (no interactions; do-file calls it model 3): -4.778** (1.830), N 191,175. PI DECISION 2026-10-10: the author's archived do-file is accepted as the port. SPEC CAVEAT: the archived do-file differs from the published table. Published controls: both democracies, both autocracies, contiguity, ln distance, alliance, higher capability score, relative capabilities, peace years cubic. Archived do-file controls: bothaut, bothdem, missingpol, onemp, jointmp, lndist, atopally, lncapratio, fpyr cubic (no contiguity, no higher-capability score; adds missing-polity, one-major-power and joint-major-power dummies; log capability ratio); other coefficients differ (e.g. ln distance -0.239 published vs -0.336 here; lower dependence -2.18 vs -4.30), so the near-match on the key term is partly coincidental. ",
    "Full-sample b=%.6f se=%.6f N=%d. Common sample N=%d (history attached at outcome year t+1)."),
    full[[TERM]][["b"]], full[[TERM]][["se"]], attr(full, "n"), nrow(ds)))
dir.create(file.path(OUT, "results_p9"), recursive = TRUE, showWarnings = FALSE)
src <- file.path(OUT, "results_p8", paste0("P9_043__", TERM, ".csv")); stopifnot(file.exists(src), file.rename(src, file.path(OUT, "results_p9", basename(src))))
sink()
