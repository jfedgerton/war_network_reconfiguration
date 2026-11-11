# Rerun ADD_T05: Choi (Seulah) 2026 JPR (10.1093/jopres/xjag029) -- relational status decline (declining embeddedness in
#   diplomatic-recognition link communities) -> MID initiation, directed dyads 1816-2001
# Headline: Table 1 col 1, author do-file StatusDecline_manuscript.do:
#   logit cwinit ccgap militaryStatusGap cap_1 majpow1 logdistance contig_binary alliance_binary jointdemo sharedIGO
#         spline1 spline2 spline3, cluster(dyadid)
# Key term = ccgap (challenger's relational status decline). Single additive term -> sens_run (scope full).
# Published Table 1 not obtainable here (OUP returns 403 to scripted requests); the reproduction below is the authors' own
# model on the authors' data and is recorded as published_check = unverified. Appendix Table A4 (rare-events logit, same RHS)
# reports ccgap 2.506 (0.503), N = 861,663; the logit is checked against that value (within 10%) and N exactly.
# Data: authors' replication package (Choi_ReplicationFileszip.zip), in _incoming/Choi2026.
ROOT <- Sys.getenv("SWEEP_DATA"); OUT <- Sys.getenv("SWEEP_OUT"); CODE <- Sys.getenv("SWEEP_CODE")
stopifnot(nzchar(ROOT), nzchar(OUT), nzchar(CODE))
suppressPackageStartupMessages({library(fixest); library(dplyr); library(haven)})
set.seed(123)
WO <- file.path(OUT, "work_add20261110/ADD_T05"); dir.create(WO, recursive = TRUE, showWarnings = FALSE)
sink(file.path(WO, "log_ADD_T05.txt"), split = TRUE)

# 1. Load authors' analysis data
TIE <- "ccgap"
CTRL <- c("militaryStatusGap", "cap_1", "majpow1", "logdistance", "contig_binary", "alliance_binary", "jointdemo", "sharedIGO",
          "spline1", "spline2", "spline3")
d <- as.data.frame(zap_labels(read_dta(file.path(ROOT, "_incoming/Choi2026/StatusDecline_jpr.dta"))))
cat("rows in author file:", nrow(d), "\n")
stopifnot(all(c("cwinit", TIE, CTRL, "ccode1", "ccode2", "year", "dyadid") %in% names(d)))

f_rhs <- paste(c(TIE, CTRL), collapse = " + ")
fit <- function(rhs_add, dat) {
  m <- feglm(as.formula(paste("cwinit ~", f_rhs, rhs_add)), data = dat, family = "binomial", cluster = ~dyadid)
  ct <- coeftable(m); print(round(ct[TIE, , drop = FALSE], 5))
  list(b = ct[TIE, 1], se = ct[TIE, 2], p = ct[TIE, 4], n = nobs(m))
}

# 2. Reproduce Table 1 col 1 on the author's full sample
cat("\n==== Reproduction (author sample) ====\n")
r_pub <- fit("", d)
b_chk <- 2.506   # appendix Table A4 (relogit, same RHS and sample); Table 1 logit itself not obtainable
cat(sprintf("REPRO %s: %.4f (se %.4f) N %d | appendix relogit %.3f (0.503) N 861663\n", TIE, r_pub$b, r_pub$se, r_pub$n, b_chk))
stopifnot(r_pub$n == 861663, abs(r_pub$b - b_chk) / abs(b_chk) <= 0.10)   # fail loudly
b_pub <- r_pub$b   # published_check = unverified: our reproduction of the author's own model stands in for Table 1
n_orig <- r_pub$n

# 3. Merge history (directed dyad, outcome year = year); restrict to 1822-2011
source(file.path(CODE, "hist/add_history.R"))
d <- add_history(d, "ccode1", "ccode2", "year")
d <- filter(d, year >= 1822, year <= 2011)
d <- d[complete.cases(d[, c("cwinit", TIE, CTRL, "dyadid", "hist_dyad_w15", "hist_mono_w15", "hist_dyad_w26", "hist_mono_w26",
                            "bm_cinc_w15", "bm_democ_w15")]), ]
cat("\nRows after year restriction + complete cases:", nrow(d), "(author model N", n_orig, ")\n")
stopifnot(nrow(d) > 0.80 * n_orig)

# 4. Models on the same sample
cat("\n==== m0 ====\n");  r0  <- fit("", d)
cat("\n==== mD ====\n");  rD  <- fit("+ hist_dyad_w15", d)
cat("\n==== mDM ====\n"); rDM <- fit("+ hist_dyad_w15 + hist_mono_w15", d)
cat("\n==== m26 ====\n"); r26 <- fit("+ hist_dyad_w26 + hist_mono_w26", d)
cat("\n==== bC ====\n");  rC  <- fit("+ bm_cinc_w15", d)
cat("\n==== bD ====\n");  rBD <- fit("+ bm_democ_w15", d)
stopifnot(length(unique(c(r0$n, rD$n, rDM$n, r26$n, rC$n, rBD$n))) == 1, r0$n == nrow(d))
# ---- sensitivity hook: runs only under run_sensitivity.R (SWEEP_SENS=1), never in the main sweep ----
if (Sys.getenv("SWEEP_SENS") == "1") {
  source(file.path(CODE, "hist/sens_helper.R"))
  sens_run(id = "ADD_T05", dat = d, y = "cwinit", tie = TIE, ctrl = CTRL, fe = NULL,
           fit = function(dat, add) { r <- fit(if (length(add)) paste("+", paste(add, collapse = " + ")) else "", dat)
             c(b = unname(r$b), se = unname(r$se), p = unname(r$p)) })
}
pct <- function(x) 100 * (x$b - r0$b) / abs(r0$b)

out <- data.frame(cid = "ADD_T05", short_cite = "Choi 2026", key_term = TIE, reproduced = "REPRODUCED", b_pub = b_pub, b_repro = r_pub$b,
  estimator = "logit, peace-year splines, SE clustered by dyad (fixest feglm)", tie_var = TIE,
  n_orig = n_orig, n_used = r0$n, b0 = r0$b, se0 = r0$se, p0 = r0$p, b_dm = rDM$b, se_dm = rDM$se, p_dm = rDM$p,
  pct_dm = pct(rDM), b_d = rD$b, pct_d = pct(rD), b_26 = r26$b, se_26 = r26$se, p_26 = r26$p, pct_26 = pct(r26),
  pct_bC = pct(rC), pct_bD = pct(rBD), author_has_history_control = "N",
  notes = sprintf("Table 1 col 1; published Table 1 not obtained (b_pub = our reproduction, unverified); appendix relogit check 2.506; dropped %d obs vs author sample (1822-2011 + history/benchmark completeness).", n_orig - r0$n))
stopifnot(all(is.finite(unlist(out[, c("b0", "se0", "p0", "b_dm", "p_dm", "b_26", "p_26", "pct_bC", "pct_bD")]))))
print(Filter(is.numeric, out))
dir.create(file.path(OUT, "results"), recursive = TRUE, showWarnings = FALSE)
write.csv(out, file.path(OUT, "results/ADD_T05.csv"), row.names = FALSE)
sink()
