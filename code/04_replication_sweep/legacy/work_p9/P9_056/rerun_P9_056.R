# P9_056 Choi 2025 FPA "Do Alliances Deter War?" (10.1093/fpa/oraf027)
# Data + do-file + log: Harvard Dataverse doi:10.7910/DVN/BJFFCW ("jl_modified2.dta", "alliance and war.do", "alliance and war.smcl");
#   staged as jl_modified2.dta / choi.do / choi.smcl. Directed challenger-target dyad-years (Kenwick & Vasquez 2017 build of the
#   Johnson & Leeds data).
# HEADLINE: Table 2 Model 7 (war DV, full sample; the war models are the ones plotted in Figure 1):
#   logit war fbdef nfbdef pchaloff pchalneu ln_distance capprop jdem s_un_glo yearout1w yearout2w yearout3w, cluster(ddyad)
#   yearout1-3w = war peace-years cubic generated in the do-file (time since last war within ddyad; counter resets after a war year).
# KEY TERM fbdef (target holds a defense alliance from its birth). SECOND ROW nfbdef (target defense alliance formed after birth).
#   Both are properties of the TARGET's alliances (not a challenger-target tie).
#   Author's log: fbdef -2.062863 (.6133846), nfbdef -1.308343 (.2955895), N = 1,077,992.
# Estimator: fixest feglm logit, vcov ~ddyad, all six steps. DV contemporaneous (no lead).
CODE <- Sys.getenv("SWEEP_CODE"); stopifnot(nzchar(CODE))
source(file.path(CODE, "work_p8/common_b3.R"))
WD <- file.path(ROOT, "work_p9/P9_056"); W <- file.path(OUT, "work_p9/P9_056"); dir.create(W, recursive = TRUE, showWarnings = FALSE)
sink(file.path(W, "log_P9_056.txt"), split = TRUE)
keep <- c("challenger", "target", "ddyad", "year", "war", "fbdef", "nfbdef", "pchaloff", "pchalneu", "ln_distance", "capprop", "jdem", "s_un_glo")
d <- as.data.frame(haven::zap_labels(haven::read_dta(file.path(WD, "jl_modified2.dta"), col_select = all_of(keep))))
stopifnot(!anyDuplicated(d[, c("ddyad", "year")]))
cat("rows", nrow(d), " year range:", range(d$year), "\n")
# do-file peace-years for war: sort ddyad year; time1 = war[_n-1] (0 for first row); time2 = running sum(time1) (Stata sum: missing -> 0);
# yearspeacew = year - first year of the (ddyad, time2) spell + 1
d <- d[order(d$ddyad, d$year), ]
w0 <- ifelse(is.na(d$war), 0, d$war)
first <- !duplicated(d$ddyad)
t1 <- c(0, w0[-nrow(d)]); t1[first] <- 0
d$time2 <- ave(t1, d$ddyad, FUN = cumsum)
d$yearspeacew <- d$year - ave(d$year, d$ddyad, d$time2, FUN = min) + 1
d$yearout1w <- d$yearspeacew; d$yearout2w <- d$yearspeacew^2; d$yearout3w <- d$yearspeacew^3
terms <- c("fbdef", "nfbdef")
rhs0 <- c("fbdef", "nfbdef", "pchaloff", "pchalneu", "ln_distance", "capprop", "jdem", "s_un_glo", "yearout1w", "yearout2w", "yearout3w")
vars <- c("war", rhs0, "ddyad")
pub <- list(fbdef = c(b = -2.062863, se = 0.6133846), nfbdef = c(b = -1.308343, se = 0.2955895)); pub_n <- 1077992
fitter <- function(r, data) fit_bin("war", r, data, "logit", "~ddyad", terms = terms)
d0 <- d[complete.cases(d[, vars]), ]
full <- fitter(rhs0, d0); rm(d0); invisible(gc())
for (t in terms) cat("\nREPRO vs log:", t, "b", full[[t]][["b"]], "(log", pub[[t]][["b"]], ") se", full[[t]][["se"]], "(log", pub[[t]][["se"]], ")")
cat("\nN", attr(full, "n"), "(log", pub_n, ")\n")
stopifnot(attr(full, "n") == pub_n, all(sapply(terms, function(t) abs(full[[t]][["b"]] - pub[[t]][["b"]]) / abs(pub[[t]][["b"]]) <= 0.01)))
ds <- common_sample(d, vars, "challenger", "target", "year"); rm(d); invisible(gc())
M <- run_six(rhs0, ds, fitter)
out <- list(M = M, full = full, n_used = nrow(ds))
saveRDS(out, file.path(W, "M_P9_056.rds"))

# ---- sensitivity hook: runs only under run_sensitivity.R (SWEEP_SENS=1), never in the main sweep ----
if (Sys.getenv("SWEEP_SENS") == "1") {
  source(file.path(CODE, "hist/sens_helper.R"))
  sens_run(id = "P9_056", dat = ds, y = "war", tie = "fbdef", ctrl = setdiff(rhs0, "fbdef"), fe = NULL,
           fit = function(dat, add) fitter(c(rhs0, add), dat)[["fbdef"]])
}

meta <- list(pid = "P9_056", cite = "Choi 2025", journal = "FPA", year = 2025, doi = "10.1093/fpa/oraf027", family = "TIE")
art <- "logit, SE clustered on directed dyad, war peace-years cubic"
used <- "feglm logit, vcov clustered ddyad, all six steps"
ahc <- "N (war peace-years cubic only)"
nt <- function(t, extra) sprintf(paste0("Directed challenger-target dyad-years (author's data). Published values from the author's Stata log (Table 2 Model 7). ",
  "Full-sample repro b=%.6f se=%.6f N=%d vs log %s (%s) N=%d. Common 1822-2011 sample N=%d. %s"), full[[t]][["b"]], full[[t]][["se"]],
  attr(full, "n"), format(pub[[t]][["b"]]), format(pub[[t]][["se"]]), pub_n, nrow(ds), extra)
write_p8(meta, "fbdef", "title: 'Do Alliances Deter War?'; Figure 1 plots fbdef and nfbdef from the war models",
  "Table 2 Model 7 (war, full sample; author's log)", art, used, pub$fbdef[["b"]], pub$fbdef[["se"]], pub_n, TRUE, FALSE, M, ahc,
  nt("fbdef", "Headline (target's birth-era defense alliance). nfbdef kept in all steps. Key term is a property of the target's alliances."))
write_p8(meta, "nfbdef", "title: 'Do Alliances Deter War?'; Figure 1 plots fbdef and nfbdef from the war models",
  "Table 2 Model 7 (war, full sample; author's log)", art, used, pub$nfbdef[["b"]], pub$nfbdef[["se"]], pub_n, TRUE, FALSE, M, ahc,
  nt("nfbdef", "Non-headline second row (target's later-formed defense alliance) from the same model."))
dir.create(file.path(OUT, "results_p9"), recursive = TRUE, showWarnings = FALSE)
for (tm in terms) {
  src <- file.path(OUT, "results_p8", paste0("P9_056__", tm, ".csv")); dst <- file.path(OUT, "results_p9", paste0("P9_056__", tm, ".csv"))
  stopifnot(file.exists(src), file.rename(src, dst))
}
sink()
