# Rerun 72: DiLorenzo & Rooney 2018 JPR -- Table 1 Model 2 (rivalry dyad-years)
# logit mid_dummy sols_count_exclusive sols_count_shared ally_count_total stalemate past_mids mindem contiguity peace_years_mid, cluster(rivdyad)
ROOT <- Sys.getenv("SWEEP_DATA"); OUT <- Sys.getenv("SWEEP_OUT"); CODE <- Sys.getenv("SWEEP_CODE")
stopifnot(nzchar(ROOT), nzchar(OUT), nzchar(CODE))
suppressPackageStartupMessages({library(haven); library(dplyr); library(sandwich); library(lmtest)})
WO <- file.path(OUT, "work/72"); dir.create(WO, recursive = TRUE, showWarnings = FALSE)
sink(file.path(WO, "log_72.txt"), split = TRUE)

# 1. Load authors' data
d <- as.data.frame(zap_labels(read_dta(file.path(ROOT, "work/72/dr_dyadic_data.dta"))))
stopifnot(nrow(d) == 3752)
X <- "sols_count_exclusive + sols_count_shared + ally_count_total + stalemate + past_mids + mindem + contiguity + peace_years_mid"
fit <- function(add, dat) {
  m <- glm(as.formula(paste("mid_dummy ~", X, add)), data = dat, family = binomial("logit"))
  used <- as.numeric(rownames(model.frame(m)))
  V <- vcovCL(m, cluster = dat$rivdyad[used], type = "HC1")   # Stata cluster(rivdyad)
  ct <- coeftest(m, vcov. = V); print(ct)
  list(b = ct["sols_count_exclusive", 1], se = ct["sols_count_exclusive", 2], p = ct["sols_count_exclusive", 4], n = nobs(m))
}
# 2. Reproduce. Published: 0.100 (0.0428)* -- Appendix Table IV col 1, identical spec to Table 1 M2 (do-file lines 41 & 111)
cat("==== Reproduction ====\n"); r_pub <- fit("", d); b_pub <- 0.100
n_orig <- r_pub$n; stopifnot(n_orig == 3545)
rep_ok <- abs(r_pub$b - b_pub) / abs(b_pub) <= 0.10 && r_pub$p < .05

# 3. History. Timing check (interactive): mid_dummy matches GML onsets in the SAME year (552/698 vs 218-223 at +-1), so merge on year.
source(file.path(CODE, "hist/add_history.R"))
rownames(d) <- NULL
d <- add_history(d, "ccode1", "ccode2", "year")
d <- filter(d, year >= 1822, year <= 2011)
vars <- c("mid_dummy", strsplit(X, " \\+ ")[[1]], "hist_dyad_w15","hist_mono_w15","hist_dyad_w26","hist_mono_w26","bm_cinc_w15","bm_democ_w15")
d <- d[complete.cases(d[, vars]), ]; rownames(d) <- NULL
cat("Rows used:", nrow(d), "\n")

# 4. Models (author's past_mids + peace_years_mid retained)
cat("==== m0 ====\n");  r0  <- fit("", d)
cat("==== mD ====\n");  rD  <- fit("+ hist_dyad_w15", d)
cat("==== mDM ====\n"); rDM <- fit("+ hist_dyad_w15 + hist_mono_w15", d)
cat("==== m26 ====\n"); r26 <- fit("+ hist_dyad_w26 + hist_mono_w26", d)
cat("==== bC ====\n");  rC  <- fit("+ bm_cinc_w15", d)
cat("==== bD ====\n");  rBD <- fit("+ bm_democ_w15", d)
# ---- sensitivity hook: runs only under run_sensitivity.R (SWEEP_SENS=1), never in the main sweep ----
if (Sys.getenv("SWEEP_SENS") == "1") {
  source(file.path(CODE, "hist/sens_helper.R"))
  sens_run(id = "72", dat = d, y = "mid_dummy", tie = "sols_count_exclusive", ctrl = setdiff(strsplit(X, " \\+ ")[[1]], "sols_count_exclusive"), fe = NULL,
           fit = function(dat, add) { r <- fit(if (length(add)) paste("+", paste(add, collapse = " + ")) else "", dat); c(b = r$b, se = r$se, p = r$p) })
}
pct <- function(x) 100 * (x$b - r0$b) / abs(r0$b)
out <- data.frame(cid = 72, short_cite = "DiLorenzo & Rooney 2018", reproduced = ifelse(rep_ok, "REPRODUCED", "NOT_REPRODUCED"),
  b_pub = b_pub, b_repro = r_pub$b, estimator = "logit, rivalry-dyad clustered SE (Stata HC1)", tie_var = "sols_count_exclusive",
  n_orig = n_orig, n_used = r0$n, b0 = r0$b, se0 = r0$se, p0 = r0$p, b_dm = rDM$b, se_dm = rDM$se, p_dm = rDM$p, pct_dm = pct(rDM),
  b_d = rD$b, pct_d = pct(rD), b_26 = r26$b, se_26 = r26$se, p_26 = r26$p, pct_26 = pct(r26), pct_bC = pct(rC), pct_bD = pct(rBD),
  author_has_history_control = "Y",
  notes = "b_pub from appendix Table IV col1 (same spec as T1 M2; main text paywalled). Rivals-only sample; author already controls past_mids + peace years. DV same-year (history merged on year). 175 obs dropped: bm_cinc_w15 NA (no lagged NMC; hist vars never NA).")
print(Filter(is.numeric, out))   # port: no citation/text columns in logs
dir.create(file.path(OUT, "results"), recursive = TRUE, showWarnings = FALSE)
write.csv(out, file.path(OUT, "results/72.csv"), row.names = FALSE)
sink()
