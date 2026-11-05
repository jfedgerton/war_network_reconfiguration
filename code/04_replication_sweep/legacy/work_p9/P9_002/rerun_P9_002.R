# P9_002 Conrad 2017 BJPS 47(4) "How Democratic Alliances Solve the Power Parity Problem" (10.1017/S0007123415000538)
# Data + do-file: Harvard Dataverse doi:10.7910/DVN/AMDYRS (ConradBJPS-Master Data.dta, "Tables 1-2 Figure 1.do").
# HEADLINE: Table 2, heteroskedastic probit (do-file second Table 2 model; the model Figure 1 and every appendix robustness table use):
#   hetprob MID powerparity FINALalliesmin2 interaction2 atopally jntd6 majdyad contig150 peaceyrs p2 p3,
#     het(powerparity FINALalliesmin2 interaction2) vce(robust)
#   Directed dyads (cc1 initiates MID vs cc2); FINALalliesmin2 = minimum polity2 of cc2's (target's) allies; powerparity = alliance-inclusive
#   capability distribution; interaction2 = powerparity * FINALalliesmin2 (verified identical in the data).
# KEY TERM (a priori, abstract): "uncertainty generated at power parity is mitigated when a target state's allies are more democratic"
#   -> product term interaction2 (mean equation).
# Published (article Table 2 Model 2, mean equation): interaction2 -0.05 (0.03) ** (p<.05, two-tailed), N = 633,178.
# Estimator: glmx::hetglm binomial probit, scale link log (= Stata hetprob: Pr = Phi(xb / exp(zg))); robust (HC0 sandwich x N/(N-1),
#   Stata's ML robust adjustment). History terms enter the mean equation only; the variance equation is unchanged. All six steps.
CODE <- Sys.getenv("SWEEP_CODE"); stopifnot(nzchar(CODE))
source(file.path(CODE, "work_p8/common_b3.R"))
suppressPackageStartupMessages({library(glmx); library(sandwich)})
WD <- file.path(ROOT, "work_p9/P9_002"); W <- file.path(OUT, "work_p9/P9_002"); dir.create(W, recursive = TRUE, showWarnings = FALSE); sink(file.path(W, "log_P9_002.txt"), split = TRUE)
d <- as.data.frame(haven::zap_labels(haven::read_dta(file.path(WD, "ConradBJPS-Master_Data.dta"))))
stopifnot(!anyDuplicated(d[, c("ccode1", "ccode2", "year")]))
ok <- !is.na(d$interaction2); stopifnot(max(abs(d$interaction2[ok] - d$powerparity[ok] * d$FINALalliesmin2[ok]), na.rm = TRUE) < 1e-5)
rhs0 <- c("powerparity", "FINALalliesmin2", "interaction2", "atopally", "jntd6", "majdyad", "contig150", "peaceyrs", "p2", "p3")
het <- c("powerparity", "FINALalliesmin2", "interaction2")
fit_hp <- function(rhs, data) {
  f <- as.formula(paste("MID ~", paste(rhs, collapse = " + "), "|", paste(het, collapse = " + ")))
  m <- hetglm(f, data = data, family = binomial(link = "probit"), link.scale = "log")
  stopifnot(m$converged, m$optim$convergence == 0)
  n <- nobs(m); V <- sandwich(m) * n / (n - 1); b <- coef(m); se <- sqrt(diag(V)); stopifnot(all(is.finite(se)))
  ct <- cbind(b = b, se = se, z = b / se, p = 2 * pnorm(-abs(b / se)))
  cat("\n---- hetprob", deparse(f, width.cutoff = 500), "| N =", n, "| logLik =", round(as.numeric(logLik(m)), 3), "| robust\n"); print(round(ct, 6))
  out <- list(interaction2 = c(b = ct["interaction2", "b"], se = ct["interaction2", "se"], p = ct["interaction2", "p"]))
  attr(out, "n") <- n; rm(m); invisible(gc()); out
}
vars <- c("MID", rhs0)
full <- d[complete.cases(d[, vars]), ]; cat("author sample years:", range(full$year), "\n")
r_full <- fit_hp(rhs0, full); n_orig <- attr(r_full, "n"); rm(full); invisible(gc())
ds <- common_sample(d, vars, "ccode1", "ccode2", "year"); rm(d); invisible(gc())
M <- run_six(rhs0, ds, fit_hp)
saveRDS(list(M = M, r_full = r_full, n_orig = n_orig, n_used = nrow(ds)), file.path(W, "M_P9_002.rds"))
# ---- result-CSV write step (write_p8 writes to results_p8; the file is moved to results_p9 unchanged) ----
PUBNOTE <- paste("Published values from article Table 2 Model 2 (-0.05, SE 0.03, ** p<.05 two-tailed, N = 633,178); author-sample repro matches b and SE at the",
  "article's two-decimal precision and N exactly. The article's ** is not attainable from -0.05/0.03 two-tailed (z ~ 1.7-1.8, p ~ .07-.09); repro p = .070.",
  "Article prints Power Distribution -0.06 (0.18) in the mean equation; repro gives +0.060 (0.21) at a unique optimum (checked from alternative starts).")
meta <- list(pid = "P9_002", cite = "Conrad 2017", journal = "BJPS", year = 2017, doi = "10.1017/S0007123415000538", family = "TIE")
write_p8(meta, "interaction2", "abstract: 'power parity is mitigated when a target state's allies are more democratic'",
  "Table 2 Model 2 heteroskedastic probit (do-file: hetprob MID ... interaction2 ..., het(powerparity FINALalliesmin2 interaction2) vce(robust))",
  "heteroskedastic probit (hetprob), variance eq. powerparity FINALalliesmin2 interaction2; robust SE",
  "glmx::hetglm probit, log scale link, same variance equation; HC0 sandwich x N/(N-1); all six steps (history in mean equation)",
  -0.05, 0.03, 633178, TRUE, FALSE, M, "N (peace years + squared + cubed only)",
  sprintf("%s Directed dyads; authors' Dataverse master data + do-file. Full N=%d, b=%.4f se=%.4f p=%.3g; common 1822-2011 N=%d.",
    PUBNOTE, n_orig, r_full$interaction2[["b"]], r_full$interaction2[["se"]], r_full$interaction2[["p"]], nrow(ds)))
src <- file.path(OUT, "results_p8", "P9_002__interaction2.csv"); dst <- file.path(OUT, "results_p9", "P9_002__interaction2.csv")
dir.create(dirname(dst), recursive = TRUE, showWarnings = FALSE); stopifnot(file.exists(src), file.rename(src, dst))
sink()
