# Rerun 55: Kimball 2006 JPR "Alliance formation and conflict initiation: the missing link"
# Headline: Table I Bivariate Probit Model 2 (recursive; allform2 enters the conflict equation), robust SE:
#   biprobit (allform2 = RELCAP logdist contig jtdem jaut sharerival allyrs _spline1af-3af)
#            (conflict = RELCAP contig jtdem jaut pwrs allform2 peaceyrs _spline1-3), robust
# Tie = allform2 in the conflict equation. Estimator hand-coded (work/55/biprobit.R; passed synthetic tests
# in test_biprobit.R: BVN CDF vs mvtnorm TVPACK < 1e-10, parameter recovery, analytic scores).
# Published table not obtainable (paywalled; Laval repository copy restricted) -> UNVERIFIED.
# History + benchmark variables are added to the CONFLICT equation only.
ROOT <- Sys.getenv("SWEEP_DATA"); OUT <- Sys.getenv("SWEEP_OUT"); CODE <- Sys.getenv("SWEEP_CODE")
stopifnot(nzchar(ROOT), nzchar(OUT), nzchar(CODE))
source(file.path(CODE, "work/rerun_helpers_s68.R"))
WO <- file.path(OUT, "work/55"); dir.create(WO, recursive = TRUE, showWarnings = FALSE)  # port: outputs
W <- file.path(ROOT, "work/55"); sink(file.path(WO, "log_55.txt"), split = TRUE)
source(file.path(CODE, "work/55/biprobit.R")); library(haven)

d <- as.data.frame(zap_labels(read_dta(file.path(W, "file47984_missinglink_jprfinal.dta"))))
names(d) <- sub("^_", "x_", names(d))
stopifnot(nrow(d) == 94163)
eqA <- c("RELCAP","logdist","contig","jtdem","jaut","sharerival","allyrs","x_spline1af","x_spline2af","x_spline3af")
eqC <- c("RELCAP","contig","jtdem","jaut","pwrs","allform2","peaceyrs","x_spline1","x_spline2","x_spline3")
TIE <- "eq2:allform2"

run_bp <- function(dat, extraC = character(0), label = "") {
  vars <- unique(c("allform2","conflict", eqA, eqC, extraC))
  dat <- dat[complete.cases(dat[, vars]), ]
  X1 <- cbind(`_cons` = 1, as.matrix(dat[, eqA])); X2 <- cbind(`_cons` = 1, as.matrix(dat[, c(eqC, extraC)]))
  f <- fit_biprobit(X1, X2, dat$allform2, dat$conflict)
  se <- sqrt(diag(f$V)); z <- f$coef / se
  cat("\n---- [", label, "] biprobit; conflict eq extra:", paste(extraC, collapse = ","), "| N =", f$n,
      " ll =", round(f$ll, 3), " rho =", round(f$rho, 4), "\n")
  print(round(cbind(b = f$coef, se = se, z = z, p = 2 * pnorm(-abs(z))), 5))
  list(b = unname(f$coef[TIE]), se = unname(se[TIE]), p = unname(2 * pnorm(-abs(z[TIE]))), n = f$n)
}

# 2. Reproduce (plus single-probit column from the same do-file for reference)
pr <- glm(reformulate(eqC, "conflict"), data = d, family = binomial("probit"))
cat("Single probit (Table I col 1) allform2:", coef(pr)["allform2"], " N =", nobs(pr), "\n")
rep <- run_bp(d, label = "repro Model 2")
n_orig <- rep$n
status <- "UNVERIFIED"; b_pub <- NA

# 3. History
d <- add_history(d, "ccode1", "ccode2", "year")
s <- d[d$year >= 1822 & d$year <= 2011, ]
s <- s[complete.cases(s[, unique(c("allform2","conflict", eqA, eqC, "hist_dyad_w15","hist_mono_w15","hist_dyad_w26","hist_mono_w26"))]), ]
cat("Estimation sample", nrow(s), "; dropped vs original", n_orig - nrow(s), "(years 1816-1821 removed)\n")

# 4. Models
m0  <- run_bp(s, label = "m0")
mD  <- run_bp(s, "hist_dyad_w15", "mD")
mDM <- run_bp(s, c("hist_dyad_w15","hist_mono_w15"), "mDM")
m26 <- run_bp(s, c("hist_dyad_w26","hist_mono_w26"), "m26")
sC <- s[!is.na(s$bm_cinc_w15), ]; sD <- s[!is.na(s$bm_democ_w15), ]
cat("bC sample", nrow(sC), " bD sample", nrow(sD), "\n")
b0C <- if (nrow(sC) < nrow(s)) run_bp(sC, label = "m0 on bC sample")$b else NULL
bC  <- run_bp(sC, "bm_cinc_w15", "bC")
b0D <- if (nrow(sD) < nrow(s)) run_bp(sD, label = "m0 on bD sample")$b else NULL
bD  <- run_bp(sD, "bm_democ_w15", "bD")

write_row(55, "Kimball 2006", status, b_pub, rep$b, "recursive bivariate probit (hand-coded ML), robust SE", "allform2 (conflict eq)",
  n_orig, m0$n, m0, mD, mDM, m26, bC, bD, b0C, b0D, ahc = "Y",
  notes = sprintf("Table I biprobit M2; pub table unobtainable -> UNVERIFIED (estimator validated on synthetic data). History/benchmarks added to conflict eq only. Single-probit allform2=%.3f. Dropped %d obs (1816-21). bC N=%d, bD N=%d.",
                  coef(pr)["allform2"], n_orig - m0$n, bC$n, bD$n))
sink()
