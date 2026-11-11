# Rerun ADD_T01: Mazumder 2017 JPR (10.1177/0022343316687018) -- PTA network centrality x challenger regime type -> MID initiation
# Headline: main Table, Model 4 "Full Specification" (directed dyads 1965-1999), author do-file IOConflictMain.do:
#   logit mzinit PTACent1_lag democracy_1_lag PTACentDemoc_1_lag democracy_2_lag igos dependlow_lag requivlag5 PTASameGT0_lag
#         oil_1_lag contigdum logdist rgdp1_lag rgdp2_lag s_wt_glo initshare majmaj minmaj majmin pcyrsmzinit pcyrsmzinits*,
#         robust cluster(dirdyadid)
# Focal (tie-containing) term = PTACent1_lag: PTA degree centrality of the challenger, i.e. its effect for AUTOCRATIC challengers
#   (democracy_1_lag = 0), the article's headline claim ("Autocracies and the International Sources of Cooperation"). It is a
#   constituent of the product PTACentDemoc_1_lag (= PTACent1_lag * democracy_1_lag), which is carried as a second (non-headline) row.
# Published point estimates: online appendix Table 2 (bootstrapped re-estimation of the main table; point estimates are the
#   original-sample logit estimates) Model 4: PTACent1_lag -7.533 (boot SE 1.693), PTACentDemoc_1_lag 7.040 (1.763), N = 452,415.
#   The main-text (clustered) SEs are not reproduced here (article paywalled); reproduction is judged on b (within 10%).
# Data: authors' replication package (Harvard Dataverse doi:10.7910/DVN/TE2SQ9, Soumyajit_Mazumder.zip), in _incoming/Mazumder2017.
ROOT <- Sys.getenv("SWEEP_DATA"); OUT <- Sys.getenv("SWEEP_OUT"); CODE <- Sys.getenv("SWEEP_CODE")
stopifnot(nzchar(ROOT), nzchar(OUT), nzchar(CODE))
suppressPackageStartupMessages({library(fixest); library(dplyr); library(haven)})
set.seed(123)
WO <- file.path(OUT, "work_add20261110/ADD_T01"); dir.create(WO, recursive = TRUE, showWarnings = FALSE)
sink(file.path(WO, "log_ADD_T01.txt"), split = TRUE)

# 1. Load authors' analysis data (only the model columns; the .dta has 261 columns, 474 MB)
SPL <- c("pcyrsmzinit", "pcyrsmzinits1", "pcyrsmzinits2", "pcyrsmzinits3")
CTRL <- c("democracy_1_lag", "PTACentDemoc_1_lag", "democracy_2_lag", "igos", "dependlow_lag", "requivlag5", "PTASameGT0_lag",
          "oil_1_lag", "contigdum", "logdist", "rgdp1_lag", "rgdp2_lag", "s_wt_glo", "initshare", "majmaj", "minmaj", "majmin", SPL)
TIE <- "PTACent1_lag"; PROD <- "PTACentDemoc_1_lag"
d <- read_dta(file.path(ROOT, "_incoming/Mazumder2017/Datasets/IOConflictWorking.dta"),
              col_select = all_of(c("mzinit", TIE, CTRL, "ccode1", "ccode2", "year", "dirdyadid")))
d <- as.data.frame(zap_labels(d))
cat("rows in author file:", nrow(d), "\n")
stopifnot(nrow(d) == 679931)                                     # appendix Table 1: MID Initiation N = 679,931
# the product term is the author's own column; check it is PTACent1_lag * democracy_1_lag where both are observed
ok <- complete.cases(d[, c(TIE, "democracy_1_lag", PROD)])
stopifnot(max(abs(d[[PROD]][ok] - d[[TIE]][ok] * d$democracy_1_lag[ok])) < 1e-6)

f_rhs <- paste(c(TIE, CTRL), collapse = " + ")
fit <- function(rhs_add, dat, key = TIE) {
  fml <- as.formula(paste("mzinit ~", f_rhs, rhs_add))
  m <- feglm(fml, data = dat, family = "binomial", cluster = ~dirdyadid)
  ct <- coeftable(m)
  print(round(ct[c(TIE, PROD), ], 5))
  list(b = ct[key, 1], se = ct[key, 2], p = ct[key, 4], n = nobs(m), m = m)
}
kv <- function(r, key) { ct <- coeftable(r$m); list(b = ct[key, 1], se = ct[key, 2], p = ct[key, 4], n = r$n) }

# 2. Reproduce Model 4 on the author's full sample
cat("\n==== Reproduction (author sample) ====\n")
r_pub <- fit("", d)
b_pub <- -7.533; b_pub_prod <- 7.040
cat(sprintf("REPRO %s: %.4f (se %.4f) N %d | published %.3f N 452415\n", TIE, r_pub$b, r_pub$se, r_pub$n, b_pub))
cat(sprintf("REPRO %s: %.4f | published %.3f\n", PROD, kv(r_pub, PROD)$b, b_pub_prod))
stopifnot(abs(r_pub$b - b_pub) / abs(b_pub) <= 0.10)             # fail loudly if not reproduced
stopifnot(abs(kv(r_pub, PROD)$b - b_pub_prod) / abs(b_pub_prod) <= 0.10)
n_orig <- r_pub$n

# 3. Merge history (directed dyad; covariates already lagged by the author, outcome year = year); restrict to 1822-2011
source(file.path(CODE, "hist/add_history.R"))
d <- add_history(d, "ccode1", "ccode2", "year")
d <- filter(d, year >= 1822, year <= 2011)
d <- d[complete.cases(d[, c("mzinit", TIE, CTRL, "dirdyadid", "hist_dyad_w15", "hist_mono_w15", "hist_dyad_w26", "hist_mono_w26",
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
stopifnot(length(unique(c(r0$n, rD$n, rDM$n, r26$n, rC$n, rBD$n))) == 1)
# ---- sensitivity hook: runs only under run_sensitivity.R (SWEEP_SENS=1), never in the main sweep ----
# headline term is a constituent of a product term -> three-way split only (sens_split)
if (Sys.getenv("SWEEP_SENS") == "1") {
  source(file.path(CODE, "hist/sens_helper.R"))
  sens_split(id = "ADD_T01", dat = d,
             fit = function(dat, add) { r <- fit(if (length(add)) paste("+", paste(add, collapse = " + ")) else "", dat)
               c(b = unname(r$b), se = unname(r$se), p = unname(r$p)) })
}

# 5. Results: one row per key term (headline constituent and the product term), same columns as other studies
row_for <- function(key, bp, brep) {
  g <- function(r) kv(r, key); b0 <- g(r0)$b; pct <- function(r) 100 * (g(r)$b - b0) / abs(b0)
  data.frame(cid = "ADD_T01", short_cite = "Mazumder 2017", key_term = key,
    reproduced = ifelse(abs(brep - bp) / abs(bp) <= 0.10, "REPRODUCED", "NOT_REPRODUCED"), b_pub = bp, b_repro = brep,
    estimator = "logit, robust SE clustered by directed dyad (fixest feglm)", tie_var = key,
    n_orig = n_orig, n_used = g(r0)$n, b0 = b0, se0 = g(r0)$se, p0 = g(r0)$p,
    b_dm = g(rDM)$b, se_dm = g(rDM)$se, p_dm = g(rDM)$p, pct_dm = pct(rDM), b_d = g(rD)$b, pct_d = pct(rD),
    b_26 = g(r26)$b, se_26 = g(r26)$se, p_26 = g(r26)$p, pct_26 = pct(r26), pct_bC = pct(rC), pct_bD = pct(rBD),
    author_has_history_control = "N",
    notes = sprintf("Model 4 full spec; published b from appendix Table 2 (bootstrap table, original-sample point estimates); peace years + 3 splines (B_duration); dropped %d obs vs author sample (1822-2011 + history/benchmark completeness).", n_orig - g(r0)$n))
}
out <- rbind(row_for(TIE, b_pub, r_pub$b), row_for(PROD, b_pub_prod, kv(r_pub, PROD)$b))
stopifnot(all(is.finite(as.matrix(out[, c("b0", "se0", "p0", "b_dm", "p_dm", "b_26", "p_26", "pct_bC", "pct_bD")]))))
print(Filter(is.numeric, out))   # no citation/text columns in logs
dir.create(file.path(OUT, "results"), recursive = TRUE, showWarnings = FALSE)
write.csv(out, file.path(OUT, "results/ADD_T01.csv"), row.names = FALSE)
sink()
