# P8_041 Barnhart, Trager, Saunders & Dafoe 2020 IO "The Suffragist Peace" (10.1017/S0020818320000508)
# Data + do-files: Dataverse doi:10.7910/DVN/U6ZBWT (DyadicData.dta, MonadicData.dta, "Dyadic Models.do", "Monadic Models .do").
# Published PDF NOT accessible (helper did not stage it; only a 2018 working paper with different models found at ceris.be/ORA),
#   and no log is posted -> published coefficients unavailable -> verdict UNVERIFIED.
# Headline models (do-files, "Main Models"):
#   Dyadic Table 5 Model 1: logit dispute2018 i.lagjointautocracy i.lagjointMSdem i.lagmsdaud i.lagmsdwsd i.lagwsdaut dependL s_un_glo
#     contiguity updatedally year distance distancex2 lagcap distancex3 eithernuclear jointnuclear minorpower peaceyears(^2,^3), robust
#     (omitted category = joint women's-suffrage democracy).
#   Monadic Table 4 Model 1 (directed dyads, initiator's regime): logit altinit i.democracynosuff i.autocracy minorpower relcap eithernuclear
#     jointnuclear ldistancef distancefx2 distancefx3 noncontig dependL s_wt_glo ally year altpeaceyears(^2,^3), robust cluster(dyadid)
#     (omitted category = democracy with women's suffrage).
# KEY TERMS (a priori, abstract "women's more pacific preferences generate a dyadic democratic peace ... as well as a monadic peace"):
#   dyadic: lagjointMSdem (joint male-suffrage-only democracy vs joint suffrage democracy); monadic: democracynosuff (non-suffrage
#   democracy vs suffrage democracy). Positive = suffrage democracies more peaceful.
# Usage: Rscript rerun_P8_041.R dyadic | monadic
args <- commandArgs(TRUE); part <- args[1]
# port: run_one.R calls the script without arguments -> run both parts, each in its own R process as in the sweep
if (length(args) == 0) {
  self <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)); stopifnot(length(self) == 1)
  for (p in c("dyadic", "monadic")) stopifnot(system2(file.path(R.home("bin"), "Rscript"), c(shQuote(self), p)) == 0)
  quit(save = "no")
}
CODE <- Sys.getenv("SWEEP_CODE"); stopifnot(nzchar(CODE))
source(file.path(CODE, "work_p8/common_b2.R"))
suppressPackageStartupMessages(library(haven))
WD <- file.path(ROOT, "work_p8/P8_041"); W <- file.path(OUT, "work_p8/P8_041"); dir.create(W, recursive = TRUE, showWarnings = FALSE); sink(file.path(W, paste0("log_P8_041_", part, ".txt")), split = TRUE)
meta <- list(pid = "P8_041", cite = "Barnhart et al. 2020", journal = "IO", year = 2020, doi = "10.1017/S0020818320000508", family = "DEMOCRACY")
if (part == "dyadic") {
  f <- file.path(WD, "DyadicData.dta"); dv <- "dispute2018"; term <- "lagjointMSdem"; vc <- "hetero"
  rhs <- c("lagjointautocracy","lagjointMSdem","lagmsdaud","lagmsdwsd","lagwsdaut","dependL","s_un_glo","contiguity","updatedally","year",
           "distance","distancex2","lagcap","distancex3","eithernuclear","jointnuclear","minorpower","peaceyears","peaceyearsx2","peaceyearsx3")
  hm <- "Table 5 Model 1 (dyadic, do-file)"; ae <- "logit, Stata robust (HC) SE, non-directed dyad-years, cubic peace-years"
  ev <- "\"women's more pacific preferences generate a dyadic democratic peace\""
} else {
  f <- file.path(WD, "MonadicData.dta"); dv <- "altinit"; term <- "democracynosuff"; vc <- "~dyadid"
  rhs <- c("democracynosuff","autocracy","minorpower","relcap","eithernuclear","jointnuclear","ldistancef","distancefx2","distancefx3",
           "noncontig","dependL","s_wt_glo","ally","year","altpeaceyears","altpeaceyears2","altpeaceyears3")
  hm <- "Table 4 Model 1 (monadic/directed, do-file)"; ae <- "logit, robust SE clustered by dyad, directed dyad-years, cubic peace-years"
  ev <- "\"as well as a monadic peace\""
}
nm <- names(read_dta(f, n_max = 1)); cat(length(nm), "vars\n")
cc <- intersect(c("ccode1","ccode2","statea","stateb","ccodea","ccodeb","state1","state2","ccode_1","ccode_2"), nm); cat("ccode cols:", cc, "\n")
if (length(cc) != 2) { print(nm); stop("ccode columns not found") }
keep <- unique(c(dv, rhs, cc, "dyadid"))
d <- zap_labels(read_dta(f, col_select = all_of(keep))) %>% as.data.frame(); invisible(gc())
cat("rows", nrow(d), "years", range(d$year), "\n")
for (v in c(term)) print(table(d[[v]], useNA = "ifany"))
fitter <- function(r, data) fit_bin(dv, r, data, "logit", vc, terms = term)
r0 <- fitter(rhs, d); n0 <- attr(r0, "n")
d <- d[complete.cases(d[, c(dv, rhs)]), ]; invisible(gc())
ds <- common_sample(d, c(dv, rhs), cc[1], cc[2], "year"); rm(d); invisible(gc())
M <- run_six(rhs, ds, fitter)
# ---- sensitivity hook: runs only under run_sensitivity.R (SWEEP_SENS=1), never in the main sweep ----
if (Sys.getenv("SWEEP_SENS") == "1" && part == "dyadic") {
  source(file.path(CODE, "hist/sens_helper.R"))
  sens_run(id = "P8_041", dat = ds, y = "dispute2018", tie = "lagjointMSdem", ctrl = setdiff(rhs, "lagjointMSdem"), fe = NULL,
           fit = function(dat, add) fitter(c(rhs, add), dat)[["lagjointMSdem"]])
}
write_p8(meta, term, ev, hm, ae, paste0("logit (fixest feglm), SE ", if (vc == "hetero") "heteroskedasticity-robust" else "clustered by dyadid (G/(G-1))", ", all steps"),
  NA, NA, NA, NA, FALSE, M, "N (peace-years cubic only)",
  sprintf("UNVERIFIED: published IO PDF not accessible and no log posted; b_repro from authors' data + do-file = %.4f (se %.4f, N=%d). Common sample N=%d. Key term is a regime-type category dummy (reference = suffrage democracy / joint suffrage democracy).",
          r0[[term]]["b"], r0[[term]]["se"], n0, nrow(ds)), verdict = "UNVERIFIED")
sink()
