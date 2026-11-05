# Rerun 13: Kenwick, Vasquez & Powers 2015 JOP -- Figure 1a / Table A5, J&L data 1816-2000:
#   cem jdem(#0) pchaloff(#0) pchalneu(#0) contig(#0) ccode1_maj(#0) ccode2_maj(#0) rivalry_bef(#0) disputes_before(#0)
#       ptargoff(#0) ptargneu(#0), tr(longitudinal_treatment)
#   logit anydispaft longitudinal_treatment [iweight=cem_weights]
# (#0) = no coarsening -> exact matching on all listed covariates. CEM weights reimplemented (Stata cem formula):
#   matched treated w=1; matched control w = (m_C/m_T) * (m_T^s/m_C^s); unmatched w=0.
# Unit: directed dyad, 10-year before/after windows; 'year' = window reference year.
ROOT <- Sys.getenv("SWEEP_DATA"); OUT <- Sys.getenv("SWEEP_OUT"); CODE <- Sys.getenv("SWEEP_CODE")
stopifnot(nzchar(ROOT), nzchar(OUT), nzchar(CODE))
source(file.path(CODE, "work/rerun_common_k6.R"))
W <- file.path(ROOT, "work/13")
WO <- file.path(OUT, "work/13"); dir.create(WO, recursive = TRUE, showWarnings = FALSE)
sink(file.path(WO, "log_13.txt"), split = TRUE)

# 1. Authors' analysis data (Dataverse 10.7910/DVN/29743, jl_10yrd_1816-2000.dta original format)
d <- as.data.frame(zap_labels(read_dta(file.path(W, "jl_10yrd_1816-2000.dta"))))
stopifnot(nrow(d) == 45721, all(c("ccode1","ccode2","year","anydispaft","longitudinal_treatment","disputes_before") %in% names(d)))
MATCH <- c("jdem","pchaloff","pchalneu","contig","ccode1_maj","ccode2_maj","rivalry_bef","disputes_before","ptargoff","ptargneu")
stopifnot(!anyNA(d[, c(MATCH, "longitudinal_treatment", "anydispaft")]))

cem_weights <- function(dat) {
  tr <- dat$longitudinal_treatment
  s <- interaction(dat[, MATCH], drop = TRUE)
  nT <- tapply(tr == 1, s, sum); nC <- tapply(tr == 0, s, sum)
  matched <- (nT[s] > 0) & (nC[s] > 0)
  mT <- sum(tr == 1 & matched); mC <- sum(tr == 0 & matched)
  w <- ifelse(!matched, 0, ifelse(tr == 1, 1, (mC / mT) * (nT[s] / nC[s])))
  cat(sprintf("CEM: strata=%d matched strata=%d | matched T=%d C=%d | unmatched T=%d C=%d\n",
      nlevels(s), sum(nT > 0 & nC > 0), mT, mC, sum(tr == 1 & !matched), sum(tr == 0 & !matched)))
  as.numeric(w)
}
BASE <- "anydispaft ~ longitudinal_treatment"
MV <- c("anydispaft","longitudinal_treatment", MATCH)
fitter <- function(add, dat) fit_tie(as.formula(paste0(BASE, add)), dat, link = "logit", tie = "longitudinal_treatment", weights = "cemw")

# 2. Reproduce on the full published sample. Author log (KVP_jop_rep.smcl): b = .3532531 (se .125996, p=.005), N = 44929;
#    strata 311, matched 89; matched T=6870 C=38059.
cat("\n==== Reproduction: Figure 1a JL 1816-2000 (full sample) ====\n")
d$cemw <- cem_weights(d)
rp <- fitter("", d)
b_pub <- 0.3532531; p_pub <- 0.005
stopifnot(rp$n == 44929)
status <- repro_status(b_pub, p_pub, rp$b, rp$p)
cat("\nReproduction status:", status, " b_repro =", rp$b, "\n")
n_orig <- rp$n

# 3. Merge history (directed ccode1/ccode2); restrict to 1822-2011 + complete hist/bm; then RE-RUN the authors' CEM on the
#    restricted sample (matching is part of the estimator, so baseline b0 uses weights recomputed on the same sample).
d <- add_history(d, "ccode1", "ccode2", "year")
d <- restrict_sample(d, MV)
cat("\nCEM re-run on restricted sample:\n")
d$cemw <- cem_weights(d)

# 4. Six models on one sample (history/benchmark terms added to the weighted outcome logit; matching set unchanged)
R <- run_six(fitter, d)

# 5. Results row
write_row(13, "Kenwick, Vasquez & Powers 2015", status, b_pub, rp$b, "CEM (exact) + weighted logit, iweight model SEs", "longitudinal_treatment",
  n_orig, R$m0$n, R, ahc = "Y",
  notes = sprintf("Fig 1a/Table A5 JL 1816-2000, DV anydispaft (MID in 10 yrs after). Tie = defense pact formation (treatment). Authors already exact-match on disputes_before (10-yr MID count) and rivalry_bef, but those are directed-initiation counts; GML hist_dyad is undirected (either side), and 567 obs have disputes_before=0 yet hist_dyad_w15>0, so hist adds info. CEM re-run on restricted sample; hist/bm vars enter outcome logit. N counts matched obs (w>0); dropped %d vs orig N=%d.", n_orig - R$m0$n, n_orig))
sink()
