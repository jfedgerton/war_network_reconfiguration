# P9_041 Fausett & Volgy 2010 ISQ 54(1) "Intergovernmental Organizations (IGOs) and Interstate Conflict" (10.1111/j.1468-2478.2009.00578.x)
# Data + do-file: Harvard Dataverse doi:10.7910/DVN/NSEOHP ("FausettVolgy.dta", "FausettVolgy.do"); staged as fausett.dta / fausett.do. No log.
# Undirected dyads of post-communist states, 1990-2004. DV Count_sum = count of conflict events in the dyad-year (contemporaneous).
# Do-file Table 2 model [1] (all conflict): nbreg Count_sum lagCount russia cap atop_ally ab_ungadiff mingdp100 distance10 dyadmin_p4demo0
#   loglagT, robust cluster(dyadid_max). Stata abbreviations resolve uniquely: cap -> capratio, loglagT -> loglagTIGO (asserted below).
# PI DECISION 2026-10-10: HEADLINE is Table 2 model [2] (low-severity conflict), do-file line 11:
#   nbreg low68 laglow high22 russia cap atop_ally ab_ungadiff mingdp100 distance10 dyadmin_p4demo0 loglagT, robust cluster(dyadid_max)
#   Model [1] (all conflict) kept as a non-headline second row in a separate output file. Published values from the article PDF.
# (superseded) HEADLINE: Table 2 model [1]; KEY TERM loglagTIGO (log of lagged shared IGO count). Published b/SE/N NOT AVAILABLE here (OUP PDF behind a
#   Cloudflare challenge; archive has no log) -> pub values NA, verdict NO_PUBLISHED pending a table lookup.
# Estimator: NB2 ML with SE clustered on dyadid_max -> fixest::fenegbin, vcov ~dyadid_max, all six steps (as P8_027 / P9_015).
CODE <- Sys.getenv("SWEEP_CODE"); stopifnot(nzchar(CODE))
source(file.path(CODE, "work_p8/common_b3.R"))
WD <- file.path(ROOT, "work_p9/P9_041"); W <- file.path(OUT, "work_p9/P9_041"); dir.create(W, recursive = TRUE, showWarnings = FALSE)
sink(file.path(W, "log_P9_041.txt"), split = TRUE)
d <- as.data.frame(haven::zap_labels(haven::read_dta(file.path(WD, "fausett.dta"))))
stopifnot(nrow(d) == 4583)
# Stata abbreviation check: each do-file stub matches exactly one variable
for (st in c("cap", "loglagT", "high22")) { m <- grep(paste0("^", st), names(d), value = TRUE); cat(st, "->", m, "\n"); stopifnot(length(m) == 1) }
d <- d[!is.na(d$state1), ]                                        # 1,469 all-missing padding rows
stopifnot(!anyDuplicated(d[, c("state1", "state2", "year")]), all(d$state1 < d$state2), all(d$dyadid_max == d$dyadid))
# lagged DV check: lagCount[t] == Count_sum[t-1] on consecutive dyad-years
o <- d[order(d$state1, d$state2, d$year), ]; no <- nrow(o)
cs <- o$state1[-1] == o$state1[-no] & o$state2[-1] == o$state2[-no] & o$year[-1] == o$year[-no] + 1
cat("lagCount == lag(Count_sum) on", sum(cs), "pairs:", mean(o$lagCount[-1][cs] == o$Count_sum[-no][cs], na.rm = TRUE), "\n"); rm(o, cs, no)
cat("laglow == lag(low68) on", sum(cs <- (function(o){no<-nrow(o); o$state1[-1]==o$state1[-no] & o$state2[-1]==o$state2[-no] & o$year[-1]==o$year[-no]+1})(o2 <- d[order(d$state1, d$state2, d$year), ])), "pairs:",
    mean(o2$laglow[-1][cs] == o2$low68[-nrow(o2)][cs], na.rm = TRUE), "\n"); rm(cs, o2)
cat("year range:", range(d$year), "\n")
TERM <- "loglagTIGO"
ctrl <- c("russia", "capratio", "atop_ally", "ab_ungadiff", "mingdp100", "distance10", "dyadmin_p4demo0")
specs <- list(
  M2 = list(y = "low68", rhs = c("laglow", "high22_68", ctrl, TERM), pub = c(0.4085, 0.1307, 2425), ahc = "Y (lagged DV laglow)",
            hm = "Table 2 model 2 (low-severity conflict events; do-file)", out = TERM,
            note = "HEADLINE (PI decision 2026-10-10): abstract/primary hypothesis concern low-severity conflict. Published Table 2 Model 2 0.4085** (0.1307), N 2,425. Author controls include contemporaneous high-severity count high22 (Stata stub, resolves to high22_68)."),
  M1 = list(y = "Count_sum", rhs = c("lagCount", ctrl, TERM), pub = c(0.2619, 0.1284, 2425), ahc = "Y (lagged DV lagCount)",
            hm = "Table 2 model 1 (all conflict events; do-file)", out = "loglagTIGO_allconflict",
            note = "Non-headline second row (all conflict). Published Table 2 Model 1 0.2619** (0.1284), N 2,425."))
mk_fit <- function(y) function(rhs, data) {
  m <- fenegbin(as.formula(paste(y, "~", paste(rhs, collapse = " + "))), data = data, vcov = ~dyadid_max, notes = FALSE)
  stopifnot(m$convStatus)
  ct <- coeftable(m); cat("\n---- nbreg", y, "~", paste(rhs, collapse = " + "), "| N =", nobs(m), "| theta =", round(m$theta, 4), "\n"); print(round(ct, 6))
  out <- list(c(b = ct[TERM, 1], se = ct[TERM, 2], p = ct[TERM, 4])); names(out) <- TERM; attr(out, "n") <- nobs(m); out
}
meta <- list(pid = "P9_041", cite = "Fausett & Volgy 2010", journal = "ISQ", year = 2010, doi = "10.1111/j.1468-2478.2009.00578.x", family = "TIE")
dir.create(file.path(OUT, "results_p9"), recursive = TRUE, showWarnings = FALSE)
res <- list()
for (s in names(specs)) {
  sp <- specs[[s]]; fit_nb <- mk_fit(sp$y); vars <- c(sp$y, sp$rhs, "dyadid_max")
  cat("\n==================== spec", s, "DV", sp$y, "\n")
  d0 <- d[complete.cases(d[, vars]), ]
  full <- fit_nb(sp$rhs, d0)
  cat("\nFULL-SAMPLE:", TERM, "b", full[[TERM]][["b"]], "se", full[[TERM]][["se"]], "N", attr(full, "n"), "| published", sp$pub, "\n")
  ds <- common_sample(d, vars, "state1", "state2", "year")
  M <- run_six(sp$rhs, ds, fit_nb)
  res[[s]] <- list(M = M, full = full, n_used = nrow(ds))
  # ---- sensitivity hook (headline spec only): runs only under run_sensitivity.R (SWEEP_SENS=1), never in the main sweep ----
  if (s == "M2" && Sys.getenv("SWEEP_SENS") == "1") {
    source(file.path(CODE, "hist/sens_helper.R"))
    sens_run(id = "P9_041", dat = ds, y = sp$y, tie = TERM, ctrl = setdiff(sp$rhs, TERM), fe = NULL,
             fit = function(dat, add) fit_nb(c(sp$rhs, add), dat)[[TERM]])
  }
  # write_p8 names the file by term; the M1 row is relabelled so the two rows get separate files
  Mw <- M; fw <- full
  if (sp$out != TERM) Mw <- lapply(M, function(x) { if (is.null(x)) return(x); n <- attr(x, "n"); names(x)[names(x) == TERM] <- sp$out; attr(x, "n") <- n; x })
  write_p8(meta, sp$out, "abstract: 'strong association between shared IGO membership and low severity conflict'",
    sp$hm, "negative binomial (NB2), robust SE clustered on dyad; lagged DV",
    "fenegbin (NB2 ML), vcov clustered dyadid_max, all six steps", sp$pub[1], sp$pub[2], sp$pub[3], TRUE, FALSE, Mw, sp$ahc,
    sprintf("Undirected post-communist dyads 1990-2004 (authors' data + do-file). %s Full-sample b=%.6f se=%.6f N=%d. Common sample N=%d.",
      sp$note, full[[TERM]][["b"]], full[[TERM]][["se"]], attr(full, "n"), nrow(ds)))
  src <- file.path(OUT, "results_p8", paste0("P9_041__", sp$out, ".csv")); stopifnot(file.exists(src), file.rename(src, file.path(OUT, "results_p9", basename(src))))
}
saveRDS(res, file.path(W, "M_P9_041.rds"))
sink()
