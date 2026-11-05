# P9_050 Li & Yang 2026 JPR "Extended trade dependence, alliance network centrality, and conflict initiation" (10.1093/jopres/xjaf030)
# Data + do-file: PRIO JPR replication archive ("LiYang_ReplicationFiles.zip": data_1870_2012.dta, 00Analysis.do, Table1.tex,
#   main_doc.pdf). Only the analysis inputs were extracted (archive 557 MB, zip md5 recorded in tmp/p9/P9_050_sources.csv).
# Directed dyads, DV midsidea = cc1 is on side A of a MID in year t (contemporaneous; IVs lagged in the cleaning script).
# HEADLINE: Table 1 Model 1 (post-1951; abstract: extended dependence deters MID initiation "conditionally between 1951 and 2012"):
#   "logit midsidea c.extdeporig##c.c2_tradecentrality_max_scaled dyaddep0 dyaddefpact ccode1_dem ccode2_dem c2maj c1capshare sscore
#    contig capdist CW post911 peaceyrmidsidea1 peaceyrmidsidea2 peaceyrmidsidea3 if 1951 <= year, cluster(ddyad)"
#   KEY TERM ext_x_tradecent = extdeporig x c2_tradecentrality_max_scaled (product term).
#   Published (manuscript Table 1 = author's esttab Table1.tex): -0.0281 (0.00812), N = 1,454,403.
# Estimator: fixest feglm logit, vcov ~ddyad, all six steps.
CODE <- Sys.getenv("SWEEP_CODE"); stopifnot(nzchar(CODE))
source(file.path(CODE, "work_p8/common_b3.R"))
WD <- file.path(ROOT, "work_p9/P9_050"); W <- file.path(OUT, "work_p9/P9_050"); dir.create(W, recursive = TRUE, showWarnings = FALSE)
sink(file.path(W, "log_P9_050.txt"), split = TRUE)
rhs_raw <- c("extdeporig", "c2_tradecentrality_max_scaled", "dyaddep0", "dyaddefpact", "ccode1_dem", "ccode2_dem", "c2maj", "c1capshare",
             "sscore", "contig", "capdist", "CW", "post911", "peaceyrmidsidea1", "peaceyrmidsidea2", "peaceyrmidsidea3")
d <- as.data.frame(haven::zap_labels(haven::read_dta(file.path(WD, "data_1870_2012.dta"),
       col_select = dplyr::all_of(c("ccode1", "ccode2", "year", "ddyad", "midsidea", rhs_raw)))))
stopifnot(nrow(d) == 1699036, !anyDuplicated(d[, c("ccode1", "ccode2", "year")]))
d <- d[d$year >= 1951, ]
cat("year range:", range(d$year), " rows:", nrow(d), "\n")
d$ext_x_tradecent <- d$extdeporig * d$c2_tradecentrality_max_scaled
term <- "ext_x_tradecent"
rhs0 <- c("extdeporig", "c2_tradecentrality_max_scaled", "ext_x_tradecent", setdiff(rhs_raw, c("extdeporig", "c2_tradecentrality_max_scaled")))
vars <- c("midsidea", rhs0, "ddyad")
pub <- c(b = -0.0281, se = 0.00812); pub_n <- 1454403
fitter <- function(r, data) fit_bin("midsidea", r, data, "logit", "~ddyad", terms = term)
d0 <- d[complete.cases(d[, vars]), ]
full <- fitter(rhs0, d0); rm(d0); invisible(gc())
cat("\nREPRO vs published:", term, "b", full[[term]][["b"]], "(pub", pub[["b"]], ") se", full[[term]][["se"]], "(pub", pub[["se"]],
    ") N", attr(full, "n"), "(pub", pub_n, ")\n")
# published b is rounded to 3 significant digits: require |b - pub| <= 0.00005 and the exact N
stopifnot(attr(full, "n") == pub_n, abs(full[[term]][["b"]] - pub[["b"]]) <= 0.00005)
ds <- common_sample(d, vars, "ccode1", "ccode2", "year"); rm(d); invisible(gc())
M <- run_six(rhs0, ds, fitter)
out <- list(M = M, full = full, n_used = nrow(ds))
saveRDS(out, file.path(W, "M_P9_050.rds"))

meta <- list(pid = "P9_050", cite = "Li & Yang 2026", journal = "JPR", year = 2026, doi = "10.1093/jopres/xjaf030", family = "TIE")
art <- "logit, SE clustered on directed dyad, peace years + square + cube; 1951-2012"
used <- "feglm logit, vcov clustered ddyad, all six steps"
ahc <- "N (peace years + square + cube only)"
write_p8(meta, term,
  "abstract: extended dependence deters MID initiations 'conditionally between 1951 and 2012'",
  "Table 1 Model 1 (extended dependence x trade-dependence-weighted alliance centrality, MID initiation, 1951-2012)", art, used,
  pub[["b"]], pub[["se"]], pub_n, TRUE, FALSE, M, ahc,
  sprintf(paste0("Directed dyads 1951-2012 (author's data). Published values from the manuscript Table 1 (= author's esttab Table1.tex). ",
    "Full-sample repro b=%.6f se=%.6f N=%d vs pub %s (%s) N=%d. Common 1822-2011 sample N=%d. ",
    "Interaction: constituent terms extdeporig and c2_tradecentrality_max_scaled kept in all steps."),
    full[[term]][["b"]], full[[term]][["se"]], attr(full, "n"), format(pub[["b"]]), format(pub[["se"]]), pub_n, nrow(ds)))
dir.create(file.path(OUT, "results_p9"), recursive = TRUE, showWarnings = FALSE)
src <- file.path(OUT, "results_p8", paste0("P9_050__", term, ".csv")); dst <- file.path(OUT, "results_p9", paste0("P9_050__", term, ".csv"))
stopifnot(file.exists(src), file.rename(src, dst))
sink()
