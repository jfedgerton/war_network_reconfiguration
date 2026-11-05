# P8_032 Choi 2015/2016 ISQ "A Menace to the Democratic Peace? Dyadic and Systemic Difference" (10.1111/isqu.12214)
# Data: Choi's replication data as re-posted in the ISQ symposium dataverse doi:10.7910/DVN/G1Y4D0
#   (gartzke_weisiger_ISQ_2014_table_3_data_012014_skewness.dta; readme: "based on the data used by Choi in his response").
# PDF: author copy (whanchoi.github.io), Table 1 p.575. Headline = Model T2 ("With skewness"): Gartzke & Weisiger (2014) Table 3
#   Model 14 with Dyadic Difference replaced by Dyadic Democracy (high) and Systemic Difference (SD of polity) replaced by its skewness.
#   G&W Model 14 (gw_t3.do): relogit deadlyl polave pcenerg diff1 demloi engypop dydiff logdist cntgdumy allydumy capratio onemajor
#   deadyrs deadyer*, cluster(dyadid). T2: diff1 -> skewness (sysdems), dydiff -> demhii. N 619,102.
#   Published T2: demloi -0.170 (0.034); Systemic Difference (skewness) 0.348 (0.270); demhii 0.124 (0.025).
# KEY TERMS (a priori, abstract): demloi -- "democracy emerges as a pacifying factor"; systemic difference (skewness) --
#   "no dyadic and systemic difference persists".
# Estimator all steps: King-Zeng rare-events logit (relogit, bias-corrected; work/relogit_kz_ea1.R), SE clustered by dyadid.
CODE <- Sys.getenv("SWEEP_CODE"); stopifnot(nzchar(CODE))
source(file.path(CODE, "work_p8/common_b2.R"))
source(file.path(CODE, "work/relogit_kz_ea1.R"))
suppressPackageStartupMessages(library(haven))
WD <- file.path(ROOT, "work_p8/P8_032"); W <- file.path(OUT, "work_p8/P8_032"); dir.create(W, recursive = TRUE, showWarnings = FALSE); sink(file.path(W, "log_P8_032.txt"), split = TRUE)
nm <- names(read_dta(file.path(WD, "choi_skew.dta"), n_max = 1)); print(nm)
spl <- grep("^deadyer", nm, value = TRUE); cat("splines:", spl, "\n")
cc <- intersect(c("ccode1","ccode2","statea","stateb","state1","state2"), nm); cat("ccode cols:", cc, "\n")
keep <- c("deadlyl","polave","pcenerg","sysdems","demloi","demhii","engypop","logdist","cntgdumy","allydumy","capratio","onemajor",
          "deadyrs", spl, "dyadid", "year", cc)
d <- zap_labels(read_dta(file.path(WD, "choi_skew.dta"), col_select = all_of(keep))) %>% as.data.frame()
cat("rows", nrow(d), "years", range(d$year), "\n")
# dyad codes: dyadid = ccode1*1000 + ccode2 (as in P8_031, same G&W data base); checked against COW state codes
d$c1 <- d$dyadid %/% 1000; d$c2 <- d$dyadid %% 1000
cs <- unique(peacesciencer::cow_states$ccode)
cat("share c1 in COW:", mean(d$c1 %in% cs), " c2 in COW:", mean(d$c2 %in% cs), " c1<c2:", mean(d$c1 < d$c2), "\n")
stopifnot(mean(d$c1 %in% cs) > .99, mean(d$c2 %in% cs) > .99, !anyDuplicated(d[, c("dyadid", "year")]))
# DV timing diagnostic (same as P8_031): after an event in t, does peace-years reset at t+1 (contemporaneous) or t+2 (lead)?
d <- d[order(d$dyadid, d$year), ]
nx <- function(x, k) ave(x, d$dyadid, FUN = function(z) c(z[-(1:k)], rep(NA, k)))
ev <- d$deadlyl == 1
r1 <- mean(nx(d$deadyrs, 1)[ev] == 0, na.rm = TRUE); r2 <- mean(nx(d$deadyrs, 2)[ev] == 0, na.rm = TRUE)
cat("P(deadyrs=0 at t+1 | event t) =", round(r1, 3), "; at t+2 =", round(r2, 3), "\n")
lead_dv <- r2 > r1; d$yr_dv <- d$year + as.integer(lead_dv)
cat("DV treated as", if (lead_dv) "LEAD (t+1): history keyed to year+1" else "contemporaneous: history keyed to year", "\n")
rhs <- c("polave","pcenerg","sysdems","demloi","demhii","engypop","logdist","cntgdumy","allydumy","capratio","onemajor","deadyrs", spl)
terms <- c("demloi", "sysdems")
fitter <- function(r, data) {
  f <- as.formula(paste("deadlyl ~", paste(r, collapse = " + ")))
  ct <- relogit_kz(f, data, cluster = "dyadid")
  cat("\n----", deparse(f, width.cutoff = 500), "| N =", attr(ct, "n"), "\n"); print(round(ct[c(terms, paste0(terms, "__unc")), ], 6))
  out <- lapply(setNames(terms, terms), function(t) c(b = ct[t, "b"], se = ct[t, "se"], p = ct[t, "p"])); attr(out, "n") <- attr(ct, "n"); out
}
r0 <- fitter(rhs, d)
ds <- common_sample(d, c("deadlyl", rhs, "dyadid"), "c1", "c2", "yr_dv")
M <- run_six(rhs, ds, fitter)

# ---- sensitivity hook: runs only under run_sensitivity.R (SWEEP_SENS=1), never in the main sweep ----
if (Sys.getenv("SWEEP_SENS") == "1") {
  source(file.path(CODE, "hist/sens_helper.R"))
  sens_run(id = "P8_032", dat = ds, y = "deadlyl", tie = "sysdems", ctrl = setdiff(rhs, "sysdems"), fe = NULL,
           fit = function(dat, add) fitter(c(rhs, add), dat)[["sysdems"]])
}

meta <- list(pid = "P8_032", cite = "Choi 2015", journal = "ISQ", year = 2015, doi = "10.1111/isqu.12214", family = "DEMOCRACY")
ae <- "rare-events logit (relogit; G&W Model 14 command), SE clustered by dyad, peace-years + splines"
eu <- "King-Zeng bias-corrected logit (relogit_kz), SE cluster(dyadid) x n/(n+k), all steps"
nt <- sprintf("Full-sample repro N=%d (pub 619,102): demloi %.4f (%.4f), sysdems %.4f (%.4f). DV deadlyl = fatal MID onset (lead) per G&W; history keyed on %s (peace-year reset diagnostic).",
  attr(r0, "n"), r0$demloi["b"], r0$demloi["se"], r0$sysdems["b"], r0$sysdems["se"], if (lead_dv) "year+1 (DV lead)" else "data year (DV contemporaneous)")
write_p8(meta, "demloi", "\"democracy emerges as a pacifying factor that reduces the likelihood of interstate dispute\"", "Table 1 Model T2",
  ae, eu, -0.170, 0.034, 619102, TRUE, FALSE, M, "N (peace-years splines only)", nt,
  verdict = verdict_fn(r0$demloi["b"], r0$demloi["p"], -0.170, TRUE, FALSE))
write_p8(meta, "sysdems", "\"no dyadic and systemic difference persists\"", "Table 1 Model T2",
  ae, eu, 0.348, 0.270, 619102, FALSE, FALSE, M, "N (peace-years splines only)", paste(nt, "Systemic difference measured as skewness of polity distribution (system-year variable)."),
  verdict = verdict_fn(r0$sysdems["b"], r0$sysdems["p"], 0.348, FALSE, FALSE))
sink()
