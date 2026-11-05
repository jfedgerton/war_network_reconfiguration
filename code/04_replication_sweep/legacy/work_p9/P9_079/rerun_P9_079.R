# P9_079 Peterson & Drury 2011, JCR 55(4) "Sanctioning Violence: The Effect of Third-Party Economic Coercion on Militarized Conflict"
#   (10.1177/0022002710393915). TIE (target sanctioned by third parties = single-state relational attribute of the target). PORT.
# Data: authors' replication zip (Dropbox; SV_rep.dta + SV_rep_do_file.do, no log). The script reads SV_rep_subset.rds, a column subset of
#   SV_rep.dta built by build_subset_P9_079.R (estimation filter fcwongo == 0 & active == 1 applied; 331,880 rows = published N).
# Author's model (do-file, Table 3 Model 1): logit fcwinit anysanctB sanctionAtoB allyAsBdpsio caprat jdem lndist s_un_glo coldwar fcwpceyrs
#   fpysq fpycu if fcwongo==0 & active==1, cluster(dyadid). DV and peace years are t+1 leads ("rather than lag most IVs, we take the future
#   (t+1) DV and peace years variables") -> history keyed on outcome_year = year + 1. Directed dyad: dyad = 1000 * ccode_A + ccode_B.
# HEADLINE (article Table 3 Model 1, read 2026-10-11): Third-party sanction (anysanctB) 0.354 (0.0909) **, N 331,880.
#   Legend *p < .05 **p < .01 ***p < .001; robust SEs; explanatory variables lagged one year.
CODE <- Sys.getenv("SWEEP_CODE"); stopifnot(nzchar(CODE))
source(file.path(CODE, "work_p8/common_b3.R"))
W <- file.path(OUT, "work_p9/P9_079"); dir.create(W, recursive = TRUE, showWarnings = FALSE)
sink(file.path(W, "log_P9_079.txt"), split = TRUE)
d <- readRDS(file.path(ROOT, "work_p9/P9_079/SV_rep_subset.rds"))
stopifnot(nrow(d) == 331880, !anyDuplicated(d[, c("dyad", "year")]))
d <- d %>% mutate(ccode1 = dyad %/% 1000, ccode2 = dyad %% 1000, outcome_year = year + 1)
stopifnot(all(pmin(d$ccode1, d$ccode2) * 1000 + pmax(d$ccode1, d$ccode2) == d$dyadid))
rhs0 <- c("anysanctB", "sanctionAtoB", "allyAsBdpsio", "caprat", "jdem", "lndist", "s_un_glo", "coldwar", "fcwpceyrs", "fpysq", "fpycu")
vars <- c("fcwinit", rhs0, "dyadid")
TERM <- "anysanctB"
fitter <- function(r, data) fit_bin("fcwinit", r, data, "logit", "~dyadid", terms = TERM)
full <- fitter(rhs0, d[complete.cases(d[, vars]), ])
cat("FULL-SAMPLE Model 1:", TERM, "b", full[[TERM]][["b"]], "se", full[[TERM]][["se"]], "N", attr(full, "n"), "(published 0.354 (0.0909) N 331,880)\n")
stopifnot(attr(full, "n") == 331880, abs(full[[TERM]][["b"]] - 0.354) < 0.0005)
ds <- common_sample(d, vars, "ccode1", "ccode2", "outcome_year")
M <- run_six(rhs0, ds, fitter)
if (Sys.getenv("SWEEP_SENS") == "1") {
  source(file.path(CODE, "hist/sens_helper.R"))
  sens_run(id = "P9_079", dat = ds, y = "fcwinit", tie = TERM, ctrl = setdiff(rhs0, TERM), fe = NULL,
           fit = function(dat, add) fitter(c(rhs0, add), dat)[[TERM]])
}
saveRDS(list(M = M, full = full), file.path(W, "M_P9_079.rds"))
meta <- list(pid = "P9_079", cite = "Peterson & Drury 2011", journal = "JCR", year = 2011, doi = "10.1177/0022002710393915", family = "TIE")
nt <- sprintf("Authors' data and do-file (Table 3 Model 1). Full-sample b=%.6f se=%.6f N=%d (published 0.354 (0.0909) N 331,880). Common sample N=%d (history at outcome year t+1).",
  full[[TERM]][["b"]], full[[TERM]][["se"]], attr(full, "n"), nrow(ds))
write_p8(meta, TERM, "abstract and Table 3 Model 1: third-party sanctions against the target increase MID initiation",
  "Table 3 Model 1 (logit, MID initiation t+1, directed dyad-years 1914-2000; including IO sanctions)",
  "logit, SE clustered on dyad, future peace years + square + cube", "feglm logit, SE clustered on dyadid, all six steps",
  0.354, 0.0909, 331880, TRUE, FALSE, M, "Y (future peace years cubic)", nt)
dir.create(file.path(OUT, "results_p9"), recursive = TRUE, showWarnings = FALSE)
src <- file.path(OUT, "results_p8", "P9_079__anysanctB.csv"); stopifnot(file.exists(src), file.rename(src, file.path(OUT, "results_p9", basename(src))))
sink()
