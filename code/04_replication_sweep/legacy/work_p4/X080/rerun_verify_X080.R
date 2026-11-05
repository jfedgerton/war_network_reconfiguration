# VERIFY X080 Ryou-Ellison & Gold 2020 II "Moral hazard at sea" -- matched to the PDF (Table 2, Model 2, p.123).
# Paper text (pp.119-122): unit = maritime claim-dyad-year, W. Hemisphere + Europe 1900-2001 (ICOW); DV = MID in the claim
#   (midissyr); "We primarily use probit regression models"; key = defensive alliance (ATOP defense);
#   controls: One EEZ, Two EEZs, One UNCLOS, Two UNCLOS, Fish Stocks (mcmigrate), Other Salience (mcrest 0-10),
#   Cap Imbalance (relcaps), Previous MIDs weighted (recmidwt); "In each model, we use dyad fixed effects ... as well as temporal
#   polynomials" (Table 2 note: dyad FE + temporal polynomials not shown). SE type not stated -> fixed a priori to dyad clusters
#   (2026-10-10; previously chosen by closeness to the published SE, which gave the same type).
#   Published M2: defense 1.235 (0.289), N = 1676 (dyads without MIDs drop out under dyad FE).
# Changes vs rerun_X080.R (pre-committed logit spec): probit; dyad FE; EEZ/UNCLOS categories added; jointdem6 removed
#   (only in Table 4 robustness); SE clustered on dyad (fixed a priori) instead of claim clusters. Nothing else changed.
ROOT <- Sys.getenv("SWEEP_DATA"); OUT <- Sys.getenv("SWEEP_OUT"); CODE <- Sys.getenv("SWEEP_CODE")
stopifnot(nzchar(ROOT), nzchar(OUT), nzchar(CODE))
source(file.path(CODE, "vb4/vhelp.R"))
WO <- file.path(OUT, "work_p4/X080"); dir.create(WO, recursive = TRUE, showWarnings = FALSE)  # port: outputs
W <- file.path(ROOT, "work_p4/X080"); sink(file.path(WO, "log_verify_X080.txt"), split = TRUE)
d <- read.delim(file.path(W, "3595566"))
stopifnot(nrow(d) == 3231, all(c("eez", "unclos", "dyad", "time3", "defense", "atopally") %in% names(d)))
d$eez1c <- as.numeric(d$eez == 1); d$eez2c <- as.numeric(d$eez == 2)
d$unc1c <- as.numeric(d$unclos == 1); d$unc2c <- as.numeric(d$unclos == 2)
ctrl <- c("eez1c", "eez2c", "unc1c", "unc2c", "mcmigrate", "mcrest", "relcaps", "recmidwt", "time", "time2", "time3")
rhs0 <- c("defense", ctrl)
x <- d[complete.cases(d[, c("midissyr", rhs0, "dyad")]), ]
# SE type: neither the article nor the author do-file on Dataverse (doi:10.7910/DVN/LFD6QA; file 3595567 is 0 bytes) states it.
# We cluster on the panel unit (dyad), the conventional choice for a dyad-FE panel, fixed a priori (not chosen by fit to the published SE).
VC <- "~dyad"
# diagnostic only (logged, never used to choose VC): key-term SE under alternative SE types
cand <- c("iid", "hetero", "~dyad", "~claim")
se_c <- sapply(cand, function(v) fit_bin("midissyr", rhs0, x, "probit", v, "dyad", "defense")$defense["se"])
names(se_c) <- cand; cat("SE diagnostic (iid/hetero/dyad/claim; pub 0.289; not used for selection):\n"); print(se_c); cat("SE type used:", VC, "\n")
FF <- function(r, s) fit_bin("midissyr", r, s, link = "probit", vc = VC, fe = "dyad", terms = "defense")

# reproduction on the full data (no history restriction); Model 1 (atopally) as a diagnostic
rep2 <- FF(rhs0, x); cat("Model 2 repro: defense b =", rep2$defense["b"], "se =", rep2$defense["se"], " N =", attr(rep2, "n"), "(pub 1.235, N 1676)\n")
rep1 <- fit_bin("midissyr", c("atopally", ctrl), x, "probit", VC, "dyad", "atopally")
cat("Model 1 diag: atopally b =", rep1$atopally["b"], "se =", rep1$atopally["se"], " N =", attr(rep1, "n"), "(pub 0.568 (0.246), N 1676)\n")

s <- common_sample(d, c("midissyr", rhs0, "dyad"), "ccode1", "ccode2", "year")
M <- run_six(rhs0, s, FF)
# ---- sensitivity hook: runs only under run_sensitivity.R (SWEEP_SENS=1), never in the main sweep ----
if (Sys.getenv("SWEEP_SENS") == "1") {
  source(file.path(CODE, "hist/sens_helper.R"))
  sens_run(id = "X080", dat = s, y = "midissyr", tie = "defense", ctrl = ctrl, fe = "dyad", fe_sep = TRUE,
           fit = function(dat, add) { r <- FF(c(rhs0, add), dat)$defense; c(b = unname(r[["b"]]), se = unname(r[["se"]]), p = unname(r[["p"]])) })
}
write_verify("X080", "Ryou-Ellison & Gold 2020", "TIE", "defense", "Table 2 Model 2 (p.123)",
  1.235, 0.289, 1676, "probit, dyad FE, temporal polynomials (cubic claim time), SE type not stated", TRUE,
  "logit instead of probit; no dyad FE; EEZ (one/two) and UNCLOS (one/two) controls omitted; jointdem6 included (only in Table 4); claim-clustered SE (paper: unspecified)",
  "Y", "SE type not reported (clustered on dyad, fixed a priori)", FALSE, M,
  sprintf(paste0("Authors' Dataverse data; spec now documented from PDF (previously pre-committed because paper paywalled). SE type not stated in the article or the (empty) Dataverse do-file; clustered on ", VC, " (panel unit), fixed a priori. Diagnostic only, not used to choose: iid/hetero/dyad/claim SEs ", paste(round(se_c, 3), collapse = "/"), " (pub 0.289). Full-data repro: b=%.3f N=%d ",
                 "(pub 1.235, N 1676); Model 1 diag atopally b=%.3f (pub 0.568). Common sample 1900-2001 with ccode/history non-missing. ",
                 "Authors already control recmidwt (weighted previous claim MIDs)."),
          rep2$defense["b"], attr(rep2, "n"), rep1$atopally["b"]))
sink()
