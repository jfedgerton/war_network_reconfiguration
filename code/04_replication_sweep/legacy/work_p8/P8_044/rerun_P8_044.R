# P8_044 Mousseau 2009 International Security "The Social Market Roots of Democratic Peace" (10.1162/isec.2009.33.4.52)
# Data + do: author's replication zip (sciences.ucf.edu/.../Replication-Files-IS-2009.zip): Mousseau_IS_2009.dta (variables already
#   constructed as in the do-file) + "Mousseau_IS2009_Stata Do File.do". Non-directed dyad-years 1961-2001, all X lagged.
# Headline = Table 1, Model 2 (fatal MIDs): logit mzfatmid dml bCIE dmloCIE oCIE lncprt mjpw cntg dist hazard dvh distdvh, cluster(dyadid)
# KEY TERMS (a priori, abstract): dml (lower Polity2 of dyad) -- "Democracy does not cause peace among nations";
#   bCIE (both contract-intensive) -- "not a single fatal militarized incident occurred between nations with contract-intensive economies":
#   bCIE perfectly predicts no fatal MID, so (as in Stata, which drops the regressor and its 1-cases) it is not estimable in the logit.
# PDF: green OA copy (JSTOR download via contentdm, found through Semantic Scholar). Table 1 Model 2 (p.67, rotated image):
#   DemocracyL -0.03 (0.02); DemocracyL x One state CIE -0.20 (0.04)***; One state CIE -0.88 (0.30)**; N 276,133; 11,994 clusters;
#   Both states CIE: 'Variable predicts peace perfectly and 10,866 observations, containing 185 clusters, not used.'
CODE <- Sys.getenv("SWEEP_CODE"); stopifnot(nzchar(CODE))
source(file.path(CODE, "work_p8/common_b2.R"))
suppressPackageStartupMessages(library(haven))
WD <- file.path(ROOT, "work_p8/P8_044"); W <- file.path(OUT, "work_p8/P8_044"); dir.create(W, recursive = TRUE, showWarnings = FALSE); sink(file.path(W, "log_P8_044.txt"), split = TRUE)
rhs0 <- c("dml","bCIE","dmloCIE","oCIE","lncprt","mjpw","cntg","dist","hazard","dvh","distdvh")
d <- zap_labels(read_dta(file.path(WD, "Mousseau_IS_2009.dta"), col_select = all_of(c("c1","c2","dyadid","year","mzfatmid", rhs0)))) %>% as.data.frame()
cat("rows", nrow(d), "years", range(d$year), " c1<c2 share:", mean(d$c1 < d$c2), "\n"); stopifnot(anyDuplicated(d[, c("dyadid","year")]) == 0)
print(with(d, table(bCIE, mzfatmid, useNA = "ifany")))
# Stata logit: bCIE != 0 predicts failure perfectly -> bCIE dropped and bCIE==1 obs dropped. Replicate explicitly.
stopifnot(sum(d$mzfatmid[d$bCIE == 1], na.rm = TRUE) == 0)
d <- d[!is.na(d$bCIE) & d$bCIE == 0, ]
rhs <- setdiff(rhs0, "bCIE"); terms <- c("dml")
fitter <- function(r, data) fit_bin("mzfatmid", r, data, "logit", "~dyadid", terms = terms)
r0 <- fitter(rhs, d)
ds <- common_sample(d, c("mzfatmid", rhs, "dyadid"), "c1", "c2", "year")
M <- run_six(rhs, ds, fitter)
meta <- list(pid = "P8_044", cite = "Mousseau 2009", journal = "IS", year = 2009, doi = "10.1162/isec.2009.33.4.52", family = "CAPITALISM")
write_p8(meta, "dml", "\"Democracy does not cause peace among nations\"", "Table 1 Model 2 (do-file; fatal MIDs)",
  "logit, SE clustered by dyad, hazard = 2^(-peace years/5) decay term", "logit (fixest feglm), SE cluster(dyadid) G/(G-1); bCIE==1 obs dropped as Stata does (perfect prediction), all steps",
  -0.03, 0.02, 276133, FALSE, FALSE, M, "N (decay function of peace years only)",
  sprintf("Full-sample repro from authors' data + do-file = %.4f (se %.4f, p %.3f, N=%d; pub -0.03 (0.02), N 276,133; pub b rounded to 2 dp). bCIE (both CIE) perfectly predicts zero fatal MIDs and is dropped with its cases, as Stata does; dml enters with dml x oCIE, so dml is the democracy slope in dyads without a CIE member. Common sample N=%d.",
          r0$dml["b"], r0$dml["se"], r0$dml["p"], attr(r0, "n"), nrow(ds)), verdict = verdict_fn(r0$dml["b"], r0$dml["p"], -0.03, FALSE, FALSE))
sink()
