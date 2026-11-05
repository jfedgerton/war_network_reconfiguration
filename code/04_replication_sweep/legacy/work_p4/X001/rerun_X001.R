# X001 Reiter & Stam 2003 APSR "Identifying the Culprit". Directed dyad-year, DV sideaa (A initiates MID vs B).
# Headline: Table 1 Model 4 (logit, robust SE clustered on idyr = dyad-year), authors' log (rep.txt) gives b_pub.
# KEY terms (decided from abstract before estimating): pdemdtar (personalist initiator -> democratic target),
#   pdemdin (democratic initiator -> personalist target); secondary abstract claims: mdemtar, sdemtar (military /
#   single-party initiator -> democratic target). CONTROL: democ (joint democracy).
ROOT <- Sys.getenv("SWEEP_DATA"); OUT <- Sys.getenv("SWEEP_OUT"); CODE <- Sys.getenv("SWEEP_CODE")
stopifnot(nzchar(ROOT), nzchar(OUT), nzchar(CODE))
source(file.path(CODE, "work_p4/xb1/p4c.R"))
suppressPackageStartupMessages(library(readstata13))
WO <- file.path(OUT, "work_p4/X001"); dir.create(WO, recursive = TRUE, showWarnings = FALSE)  # port: outputs
W <- file.path(ROOT, "work_p4/X001"); sink(file.path(WO, "log_X001.txt"), split = TRUE)

v <- c("sideaa","pdemdtar","pdemdin","sdemin","sdemtar","mdemin","mdemtar","personal","military","single","democ","contig",
       "majpow","ally","loglsrat","advanced","dispyrs","dspline1","dspline2","dspline3","idyr","statea","stateb","year")
d <- read.dta13(file.path(W, "doublereduce.dta"), convert.factors = FALSE, select.cols = v)
cat("rows:", nrow(d), " years:", range(d$year), "\n"); print(table(d$sideaa, useNA = "ifany"))

rhs <- c("pdemdtar","pdemdin","sdemin","sdemtar","mdemin","mdemtar","personal","military","single","democ","contig",
         "majpow","ally","loglsrat","advanced","dispyrs","dspline1","dspline2","dspline3")
f0 <- reformulate(rhs, "sideaa")
terms <- c("pdemdtar","pdemdin","mdemtar","sdemtar","democ")

# 1. reproduce on full original sample
d0 <- d[complete.cases(d[, c(all.vars(f0), "idyr")]), ]; n_orig <- nrow(d0); cat("n_orig", n_orig, "(log: 753456)\n")
rep <- fit_keys(f0, d0, terms, "idyr", "logit", file.path(WO, "cache_rep.rds"))
PUB <- list(pdemdtar = c(1.04494, 1), pdemdin = c(.1052964, 0), mdemtar = c(.6421512, 1), sdemtar = c(.178899, 0), democ = c(-1.04653, 1))

# 2-4. history, common sample, six models
ds <- common_sample(d, f0, "idyr", "statea", "stateb", "year"); rm(d, d0); gc()
M <- run_six(f0, ds, terms, "idyr", WO)

# ---- sensitivity hook: runs only under run_sensitivity.R (SWEEP_SENS=1), never in the main sweep ----
if (Sys.getenv("SWEEP_SENS") == "1") {
  source(file.path(CODE, "hist/sens_helper.R"))
  sens_run(id = "X001", dat = ds, y = "sideaa", tie = "pdemdtar", ctrl = setdiff(rhs, "pdemdtar"), fe = NULL,
           fit = function(dat, add) { r <- fit_keys(if (length(add)) update(f0, as.formula(paste(". ~ . +", paste(add, collapse = " + ")))) else f0, dat, "pdemdtar", "idyr")$k[["pdemdtar"]]
             c(b = unname(r["b"]), se = unname(r["se"]), p = unname(r["p"])) })
}

EV <- "abstract: 'personalist dictatorships are more likely to challenge democracies, but not vice versa'"
for (t in terms) {
  st <- repro(PUB[[t]][1], PUB[[t]][2] == 1, rep$k[[t]]["b"], rep$k[[t]]["p"])
  write_p4("X001", "DEMOCRACY", ifelse(t == "democ", "CONTROL", "KEY"), t,
    ifelse(t %in% c("mdemtar","sdemtar"), "abstract: 'military regimes and single-party regimes, are more likely to challenge democracies'", EV),
    st, PUB[[t]][1], rep$k[[t]]["b"], "logit, SE clustered on idyr (dyad-year); directed dyads", n_orig, M, "N (peace years + 3 splines only)",
    paste0("Table 1 Model 4 (logit; Model 5 relogit ~identical); b_pub from authors' Stata log. Directed dyads ", paste(range(ds$year), collapse="-"), " (Peceny et al. data doubled). ",
           "History keyed on undirected pair. n_orig=", n_orig, ", dropped ", n_orig - M$m0$n, " (history/benchmark missing or outside 1822-2011)."))
}
sink()
