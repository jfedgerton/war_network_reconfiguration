# ESTAUDIT 66: Peterson 2011 JPR, Table I Model 1. Copy of rerun_66.R (same data, spec, t+1 history anchoring, common sample)
# with the estimator fixed. Article: rare events logit (King & Zeng 1999, 2001), robust SE clustered on dyad
#   (PDF p.? "I run a series of rare events logit models"; do-file: relogit fcwinit ... , cluster(dyadid)).
# Previous steps all used ordinary logit. Now: relogit_kz() = ML logit + King-Zeng finite-sample bias correction (no tau:
# full dyad-year population, not case-control), SE = cluster-robust (dyadid) * n/(n+k). Implemented in work/relogit_kz_ea1.R
# (validated vs Zelig-style lm.influence formula and Firth in ea1/test_relogit.R). Key term = lnot1xswtglo (results_key/66).
ROOT <- Sys.getenv("SWEEP_DATA"); OUT <- Sys.getenv("SWEEP_OUT"); CODE <- Sys.getenv("SWEEP_CODE")
stopifnot(nzchar(ROOT), nzchar(OUT), nzchar(CODE))
source(file.path(CODE, "hist/add_history.R")); source(file.path(CODE, "work/relogit_kz_ea1.R"))
suppressPackageStartupMessages({library(haven); library(fixest)})
set.seed(123)
WO <- file.path(OUT, "work/66"); dir.create(WO, recursive = TRUE, showWarnings = FALSE)  # port: outputs
W <- file.path(ROOT, "work/66"); sink(file.path(WO, "log_est_66.txt"), split = TRUE)
KEY <- c("lnot1xswtglo", "lnotrade1")
X <- c("lnotrade1","lnotrade2","s_wt_glo","lnot1xswtglo","lnot2xswtglo","dep1_imp","dep2_imp","interdep_imp",
       "polity1_adj","polity2_adj","polity_int","caprat","lndist","lngdppc1","lngdppc2","fpy","fpy2","fpy3")
ADD <- list(m0 = character(0), mD = "hist_dyad_w15", mDM = c("hist_dyad_w15","hist_mono_w15"), m26 = c("hist_dyad_w26","hist_mono_w26"),
            bC = "bm_cinc_w15", bD = "bm_democ_w15")
d <- read_dta(file.path(W, "tpt_rep_data.dta"), col_select = any_of(c("fcwinit", X, "dyadid","ccode1","ccode2","year"))) |>
  zap_labels() |> as.data.frame()
d <- d[complete.cases(d[, c("fcwinit", X, "dyadid")]), ]
n_orig <- nrow(d); cat("n_orig =", n_orig, "\n"); stopifnot(n_orig == 732148)
fit <- function(dat, add, tag) {
  ct <- relogit_kz(reformulate(c(X, add), "fcwinit"), dat, cluster = "dyadid")
  cat("\n====", tag, " N =", attr(ct, "n"), " k =", attr(ct, "k"), "\n"); print(signif(ct[c(KEY, paste0(KEY, "__unc"), add), , drop = FALSE], 6))
  ct
}
rep <- fit(d, character(0), "REPRO full sample (relogit)")
cat("relogit repro lnot1xswtglo =", rep["lnot1xswtglo", "b"], " (pub -0.174; logit repro -0.17353)\n")

d$outyear <- d$year + 1                           # DV and peace years are t+1 (as rerun_66.R)
d <- add_history(d, "ccode1", "ccode2", "outyear")
d <- d[d$outyear >= 1822 & d$outyear <= 2011, ]
d <- d[complete.cases(d[, c("hist_dyad_w15","hist_mono_w15","hist_dyad_w26","hist_mono_w26","bm_cinc_w15","bm_democ_w15")]), ]
n_used <- nrow(d); cat("n_used =", n_used, " (rerun_66: 731634)\n"); stopifnot(n_used == 731634)
R <- list()
for (nm in names(ADD)) { ct <- fit(d, ADD[[nm]], nm); stopifnot(attr(ct, "n") == n_used); R[[nm]] <- ct; invisible(gc()) }
saveRDS(list(R = R, rep = rep, n_orig = n_orig, n_used = n_used), file.path(WO, "est_66.rds"))
print(t(sapply(R, function(ct) c(b = ct["lnot1xswtglo", "b"], p = ct["lnot1xswtglo", "p"], b_unc = ct["lnot1xswtglo__unc", "b"], p_unc = ct["lnot1xswtglo__unc", "p"]))))
source(file.path(CODE, "work/est_write_ea1.R"))
KT <- "lnot1xswtglo"; g <- function(ct, s = "") c(b = unname(ct[paste0(KT, s), "b"]), se = unname(ct[paste0(KT, s), "se"]), p = unname(ct[paste0(KT, s), "p"]))
RK <- lapply(R, g); RU <- lapply(R, g, "__unc")
# ---- sensitivity hook (split only): runs only under run_sensitivity.R (SWEEP_SENS=1), never in the main sweep ----
if (Sys.getenv("SWEEP_SENS") == "1") {
  source(file.path(CODE, "hist/sens_helper.R"))
  sens_split(id = "66", dat = d,
             fit = function(dat, add) g(fit(dat, add, paste("sens", paste(add, collapse = " + ")))))
}
write_est_ea1("66", "Peterson 2011", KT,
  "rare-events logit (Stata relogit, King-Zeng bias correction, no tau), SE clustered on dyad, t+1 DV with cubic t+1 peace years",
  "pdf (Peterson (2011) JPR: 'I run a series of rare events logit models (King & Zeng, 1999, 2001)'; tables 'Robust standard errors') + code (tpt_replication.do: relogit ..., cluster(dyadid))",
  "relogit: ML logit (fixest) + King-Zeng bias correction (work/relogit_kz_ea1.R), SE cluster dyadid x n/(n+k)",
  "N", "all six steps used ordinary logit (fixest feglm, cluster dyadid); King-Zeng rare-events bias correction not applied", "Y", RK, n_used,
  sprintf(paste0("KZ-corrected b0=%.6g (p=%.3g); uncorrected logit b0=%.6g (p=%.3g) = previous results_key/66. Uncorrected mDM %.6g (p=%.3g, pct_dm %.1f). ",
    "Full-sample relogit repro %.6g (pub -0.174). Correction at N=%d is tiny. Constituent lnotrade1 KZ: b0=%.4g -> b_dm=%.4g."),
    RK$m0["b"], RK$m0["p"], RU$m0["b"], RU$m0["p"], RU$mDM["b"], RU$mDM["p"], 100 * (RU$mDM["b"] - RU$m0["b"]) / abs(RU$m0["b"]),
    rep[KT, "b"], n_used, R$m0["lnotrade1", "b"], R$mDM["lnotrade1", "b"]))
sink()
