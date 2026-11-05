# Rerun 86: Hafner-Burton & Montgomery 2012 CMPS -- Table 1 "(3) Network Model" (= do-file T1C5 / T2C1)
# logit dispon_l1 PTASame DependD AccCentD polMin polMax rgdppc_p3Min traBilMax rgdpMin rgdpMax caprat allies gdpHegPerc
#   contigEQ1 majpowers lndist dispon_spl*, cluster(dyadid)
ROOT <- Sys.getenv("SWEEP_DATA"); OUT <- Sys.getenv("SWEEP_OUT"); CODE <- Sys.getenv("SWEEP_CODE")
stopifnot(nzchar(ROOT), nzchar(OUT), nzchar(CODE))
source(file.path(CODE, "work/rerun_common_b8.R"))
library(haven)
WO <- file.path(OUT, "work/86"); dir.create(WO, recursive = TRUE, showWarnings = FALSE)
sink(file.path(WO, "log_86.txt"), split = TRUE)

# 1. Load authors' data (nondirected dyad-years, covariates at t, DV = MID onset at t+1)
d <- as.data.frame(zap_labels(read_dta(file.path(ROOT, "work/86/PTAMIDTables.dta"))))
# dyadid = 1000*ccode_low + ccode_high (verified: all low<high, no dup dyad-years)
d$c1 <- d$dyadid %/% 1000; d$c2 <- d$dyadid %% 1000
stopifnot(all(d$c1 < d$c2), !any(duplicated(d[, c("dyadid","year")])))

TIE <- "PTASame"   # joint PTA membership = the dyadic tie; DependD / AccCentD (network-position terms) reported in notes
RHS <- paste("dispon_l1 ~ PTASame + DependD + AccCentD + polMin + polMax + rgdppc_p3Min + traBilMax + rgdpMin + rgdpMax + caprat +",
             "allies + gdpHegPerc + contigEQ1 + majpowers + lndist + dispon_spl0 + dispon_spl1 + dispon_spl2 + dispon_spl3")
mv <- all.vars(as.formula(RHS))
d <- d[complete.cases(d[, c(mv, "dyadid")]), ]
fit3 <- function(add, dat) {   # returns tie + the other two network terms
  f <- as.formula(paste(RHS, add))
  m <- fixest::feglm(f, data = dat, family = binomial(), cluster = ~dyadid, ssc = fixest::ssc(adj = FALSE, cluster.adj = TRUE), notes = FALSE)
  stopifnot(m$convStatus)
  ct <- fixest::coeftable(m)
  cat("\n---- ", paste(deparse(f, width.cutoff = 500), collapse = ""), " | N =", nobs(m), "\n"); print(signif(ct, 5))
  list(b = ct[TIE, 1], se = ct[TIE, 2], p = ct[TIE, 4], n = nobs(m), dep = ct["DependD", ], acc = ct["AccCentD", ])
}
fitter <- function(add, dat) fit3(add, dat)

# 2. Reproduce (published Table 1 col 3: PTASame -2.99E-02 (1.30E-02)*, DependD 0.22 (0.10)*, AccCentD 2.46 (0.83)**, N = 407,209)
cat("==== Reproduction ====\n")
r_pub <- fitter("", d)
b_pub <- -0.0299
cat("DependD:", r_pub$dep[1:2], " AccCentD:", r_pub$acc[1:2], "\n")
stopifnot(r_pub$n == 407209, abs(r_pub$b - b_pub) / abs(b_pub) <= 0.10, r_pub$p < 0.05)
n_orig <- r_pub$n

# 3. Merge history keyed to the OUTCOME year (DV is onset at t+1), so windows are t..t-4 / t-1..t-5 relative to covariate year
source(file.path(CODE, "hist/add_history.R"))
d$year_dv <- d$year + 1
d <- add_history(d, "c1", "c2", "year_dv")
d <- filter(d, year_dv >= 1822, year_dv <= 2011)
cat("\nNA counts in history/benchmark vars:\n"); print(colSums(is.na(d[, HIST_VARS])))
d <- d[complete.cases(d[, c(mv, HIST_VARS)]), ]
cat("Rows after year restriction + complete cases:", nrow(d), " (orig model N", n_orig, ")\n")

# 4. Models on the same sample
R <- run_six(fitter, d)
pc <- function(x, k) 100 * (x[[k]][1] - R$m0[[k]][1]) / abs(R$m0[[k]][1])
cat(sprintf("\nDependD: b0 %.4f (p %.4f) -> b_dm %.4f (p %.4f), pct_dm %.1f\n", R$m0$dep[1], R$m0$dep[4], R$mDM$dep[1], R$mDM$dep[4], pc(R$mDM, "dep")))
cat(sprintf("AccCentD: b0 %.4f (p %.4f) -> b_dm %.4f (p %.4f), pct_dm %.1f\n", R$m0$acc[1], R$m0$acc[4], R$mDM$acc[1], R$mDM$acc[4], pc(R$mDM, "acc")))

# ---- sensitivity hook: runs only under run_sensitivity.R (SWEEP_SENS=1), never in the main sweep ----
if (Sys.getenv("SWEEP_SENS") == "1") {
  source(file.path(CODE, "hist/sens_helper.R"))
  # key term of this study is DependD (the row extracted below), not the script's TIE (PTASame), which is therefore a control here
  sens_run(id = "86", dat = d, y = "dispon_l1", tie = "DependD",
           ctrl = c("PTASame", "AccCentD", "polMin", "polMax", "rgdppc_p3Min", "traBilMax", "rgdpMin", "rgdpMax", "caprat", "allies",
                    "gdpHegPerc", "contigEQ1", "majpowers", "lndist", "dispon_spl0", "dispon_spl1", "dispon_spl2", "dispon_spl3"), fe = NULL,
           fit = function(dat, add) { r <- fit3(if (length(add)) paste0(" + ", paste(add, collapse = " + ")) else "", dat); c(b = unname(r$dep[1]), se = unname(r$dep[2]), p = unname(r$dep[4])) })
}

# 5. Write row
write_row(86, "Hafner-Burton & Montgomery 2012", "REPRODUCED", b_pub, r_pub$b, "logit (fixest feglm), SE clustered on dyadid (Stata-style)",
  "PTASame", n_orig, R$m0$n, R, ahc = "Y",
  notes = sprintf(paste0("Table 1 col 3; ccodes decoded from dyadid; history keyed to DV year t+1; author controls BKT spline (PEACE_YEARS); ",
    "other network terms: DependD b0 %.3f p %.3f -> b_dm %.3f p %.3f (%.0f%%), AccCentD b0 %.3f p %.3f -> b_dm %.3f p %.3f (%.0f%%); dropped %d obs vs orig"),
    R$m0$dep[1], R$m0$dep[4], R$mDM$dep[1], R$mDM$dep[4], pc(R$mDM, "dep"),
    R$m0$acc[1], R$m0$acc[4], R$mDM$acc[1], R$mDM$acc[4], pc(R$mDM, "acc"), n_orig - R$m0$n))

# port (extraction only, no estimation): the master row for key term DependD was hand-copied by write_key_mine.py from log_86
# (fit3 prints signif(ct, 5)); the pipeline reports full precision instead (differences < 1e-5 relative).
s5 <- function(k, j) unname(R[[k]]$dep[j])   # dep = coeftable row: 1 = b, 2 = se, 4 = p
b0 <- s5("m0", 1)
ex <- data.frame(id = "86", key_term = "DependD", b0 = b0, se0 = s5("m0", 2), p0 = s5("m0", 4), b_dm = s5("mDM", 1), p_dm = s5("mDM", 4),
  b_26 = s5("m26", 1), p_26 = s5("m26", 4), pct_bC = 100 * (s5("bC", 1) - b0) / abs(b0), pct_bD = 100 * (s5("bD", 1) - b0) / abs(b0),
  n_used = R$m0$n, b0_full = unname(R$m0$dep[1]), b_dm_full = unname(R$mDM$dep[1]))
stopifnot(all(is.finite(unlist(ex[, -(1:2)]))))
write.csv(ex, file.path(WO, "key_86_extract.csv"), row.names = FALSE)
sink()
