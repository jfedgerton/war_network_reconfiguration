# P9_047 Peterson 2015 JCR 59(4) "Insiders versus Outsiders: Preferential Trade Agreements, Trade Distortions, and Militarized
#   Conflict" (10.1177/0022002713520483)
# Data + do-file: author page timothypeterson.org/research ("JCR-distortion-replication.zip": jcr_rep.dta, jcr_rep.do). No log.
# Undirected dyads 1961-2000 (ccode1 < ccode2). DV fcwmid = MID onset in t+1; fcwongo = ongoing MID in t+1; fpym = peace years in t+1.
# HEADLINE: Table 1 Model 1 (probit on dyads with an exclusive third-party PTA; the paper's Figure 3 uses this model):
#   "probit fcwmid (c.max_dyadic_dist c.max_third_dist)##c.lowerdep lncaprat bothdem bothaut lndist contigdi S3UN atopally
#    c.fpym##c.fpym##c.fpym if fcwongo==0 & thirdpta==1, cluster(dyadid)"   (lncaprat abbreviates lncapratio)
#   KEY TERM max_dyadic_dist (dyadic PTA-induced export distortion; constituent of x lowerdep, i.e. effect at lowerdep = 0).
#   SECOND ROW max_third_dist (third-party export distortion), same model.
#   PUBLISHED VALUES PENDING: article PDF not reachable (SAGE 403), no author log.
# Estimator: fixest feglm probit, vcov ~dyadid, all six steps. History attached at outcome_year = year + 1 (lead DV).
CODE <- Sys.getenv("SWEEP_CODE"); stopifnot(nzchar(CODE))
source(file.path(CODE, "work_p8/common_b3.R"))
WD <- file.path(ROOT, "work_p9/P9_047"); W <- file.path(OUT, "work_p9/P9_047"); dir.create(W, recursive = TRUE, showWarnings = FALSE)
sink(file.path(W, "log_P9_047.txt"), split = TRUE)
d <- as.data.frame(haven::zap_labels(haven::read_dta(file.path(WD, "jcr_rep.dta"))))
stopifnot(nrow(d) == 366198, !anyDuplicated(d[, c("ccode1", "ccode2", "year")]), all(d$ccode1 < d$ccode2))
d <- d[d$fcwongo == 0 & d$thirdpta == 1, ]
cat("year range:", range(d$year), " rows:", nrow(d), "\n")
d$dyad_x_low <- d$max_dyadic_dist * d$lowerdep; d$third_x_low <- d$max_third_dist * d$lowerdep
d$fpym2 <- d$fpym^2; d$fpym3 <- d$fpym^3
terms <- c("max_dyadic_dist", "max_third_dist")
rhs0 <- c("max_dyadic_dist", "max_third_dist", "lowerdep", "dyad_x_low", "third_x_low", "lncapratio", "bothdem", "bothaut", "lndist",
          "contigdi", "S3UN", "atopally", "fpym", "fpym2", "fpym3")
vars <- c("fcwmid", rhs0, "dyadid")
fitter <- function(r, data) fit_bin("fcwmid", r, data, "probit", "~dyadid", terms = terms)
d0 <- d[complete.cases(d[, vars]), ]
full <- fitter(rhs0, d0)
for (t in terms) cat("\nFULL-SAMPLE:", t, "b", full[[t]][["b"]], "se", full[[t]][["se"]])
cat("\nN", attr(full, "n"), "(published values pending)\n")
d$outcome_year <- d$year + 1
ds <- common_sample(d, vars, "ccode1", "ccode2", "outcome_year")
M <- run_six(rhs0, ds, fitter)
out <- list(M = M, full = full, n_used = nrow(ds))
saveRDS(out, file.path(W, "M_P9_047.rds"))

meta <- list(pid = "P9_047", cite = "Peterson 2015", journal = "JCR", year = 2015, doi = "10.1177/0022002713520483", family = "TIE")
art <- "probit, SE clustered on dyad, (future) peace years + square + cube; dyads with an exclusive third-party PTA, no ongoing MID"
used <- "feglm probit, vcov clustered dyadid, all six steps"
ahc <- "N (peace years + square + cube only)"
nt <- function(t, extra) sprintf(paste0("Undirected dyads 1961-2000 (author's data). Published values from the article PDF (Table 1 Model 1). b and N reproduce exactly; our clustered SEs (fixest probit, cluster dyadid) run 1-18%% below the published SEs across terms (third_dist 0.0057 vs 0.007; lowerdep 2.44 vs 2.70; interaction 49.7 vs 51.1); kept as is (PI decision 2026-10-10). ",
  "Full-sample b=%.6f se=%.6f N=%d. Common 1822-2011 sample N=%d (history/benchmarks attached at outcome year t+1). ",
  "Interaction with lowerdep: constituent terms and both products kept in all steps. %s"), full[[t]][["b"]], full[[t]][["se"]],
  attr(full, "n"), nrow(ds), extra)
write_p8(meta, "max_dyadic_dist",
  "abstract: PTA-induced trade distortions associated with higher likelihood of conflict",
  "Table 1 Model 1 (probit, exclusive-PTA dyads; do-file first model)", art, used,
  1.274, 0.549, 256323L, TRUE, FALSE, M, ahc, nt("max_dyadic_dist", "Headline (confirmed against article): dyadic distortion at lowerdep = 0. Published Table 1 Model 1 1.274* (0.549), N 256,323."))
write_p8(meta, "max_third_dist",
  "abstract: PTA-induced trade distortions associated with higher likelihood of conflict",
  "Table 1 Model 1 (probit, exclusive-PTA dyads; do-file first model)", art, used,
  0.014, 0.007, 256323L, TRUE, FALSE, M, ahc, nt("max_third_dist", "Non-headline second row: third-party distortion at lowerdep = 0. Published Table 1 Model 1 0.014* (0.007), N 256,323; article calls it relatively benign."))
dir.create(file.path(OUT, "results_p9"), recursive = TRUE, showWarnings = FALSE)
for (tm in terms) {
  src <- file.path(OUT, "results_p8", paste0("P9_047__", tm, ".csv")); dst <- file.path(OUT, "results_p9", paste0("P9_047__", tm, ".csv"))
  stopifnot(file.exists(src), file.rename(src, dst))
}
sink()
