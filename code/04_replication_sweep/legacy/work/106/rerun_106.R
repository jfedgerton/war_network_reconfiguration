# Rerun 106: Choi 2023 II -- Table 3 col 1
# logit mzinit lnbet1 cap1 majpow1 logdistance binarycontig binaryalliance jointIGO jointdemo peaceyears _spline1-3, cluster(dyadid)
ROOT <- Sys.getenv("SWEEP_DATA"); OUT <- Sys.getenv("SWEEP_OUT"); CODE <- Sys.getenv("SWEEP_CODE")
stopifnot(nzchar(ROOT), nzchar(OUT), nzchar(CODE))
source(file.path(CODE, "work/rerun_common_b8.R"))
library(haven)
WO <- file.path(OUT, "work/106"); dir.create(WO, recursive = TRUE, showWarnings = FALSE)
sink(file.path(WO, "log_106.txt"), split = TRUE)

# 1. Load authors' data (directed dyad-years; ccode1 = potential initiator)
v <- c("dyadid","ccode1","ccode2","year","mzinit","lnbet1","cap1","majpow1","logdistance","binarycontig","binaryalliance",
       "jointIGO","jointdemo","peaceyears","_spline1","_spline2","_spline3")
d <- as.data.frame(zap_labels(read_dta(file.path(ROOT, "work/106/DiploMid_violent.dta"), col_select = all_of(v))))
names(d) <- sub("^_", "", names(d))   # _spline1 -> spline1 (R-safe)
cat("rows:", nrow(d), " years:", range(d$year), "\n")

TIE <- "lnbet1"
RHS <- "mzinit ~ lnbet1 + cap1 + majpow1 + logdistance + binarycontig + binaryalliance + jointIGO + jointdemo + peaceyears + spline1 + spline2 + spline3"
mv <- all.vars(as.formula(RHS))
d <- d[complete.cases(d[, c(mv, "dyadid")]), ]
fitter <- function(add, dat) fit_tie_fe(as.formula(paste(RHS, add)), dat, binomial(), TIE, cluster = ~dyadid)

# 2. Reproduce. Published article table not obtainable (paywalled, no OA copy; the 2021 BU dissertation
#    chapter uses a different spec: relogit + capability ratio + cubic peace years), no author log -> UNVERIFIED.
cat("==== Reproduction (authors' data + do-file spec) ====\n")
r_pub <- fitter("", d)
b_pub <- NA
n_orig <- r_pub$n

# 3. Merge history; restrict to 1822-2011
source(file.path(CODE, "hist/add_history.R"))
d <- add_history(d, "ccode1", "ccode2", "year")
n_pre <- nrow(d)
d <- filter(d, year >= 1822, year <= 2011)
cat("Dropped by 1822-2011 year restriction:", n_pre - nrow(d), "\n")
cat("\nNA counts in history/benchmark vars:\n"); print(colSums(is.na(d[, HIST_VARS])))
d <- d[complete.cases(d[, c(mv, HIST_VARS)]), c(mv, "dyadid", HIST_VARS)]   # keep only model columns (memory)
gc(verbose = FALSE)
cat("Rows after year restriction + complete cases:", nrow(d), " (orig model N", n_orig, ")\n")

# 4. Models on the same sample
R <- run_six(fitter, d)

# ---- sensitivity hook: runs only under run_sensitivity.R (SWEEP_SENS=1), never in the main sweep ----
if (Sys.getenv("SWEEP_SENS") == "1") {
  source(file.path(CODE, "hist/sens_helper.R"))
  sens_run(id = "106", dat = d, y = "mzinit", tie = TIE,
           ctrl = c("cap1", "majpow1", "logdistance", "binarycontig", "binaryalliance", "jointIGO", "jointdemo", "peaceyears",
                    "spline1", "spline2", "spline3"), fe = NULL,
           fit = function(dat, add) { r <- fitter(if (length(add)) paste("+", paste(add, collapse = " + ")) else "", dat); c(b = r$b, se = r$se, p = r$p) })
}

# 5. Write row
write_row(106, "Choi 2023", "UNVERIFIED", b_pub, r_pub$b, "logit (fixest feglm), SE clustered on dyadid (Stata-style)",
  "lnbet1 (initiator betweenness, diplomatic network)", n_orig, R$m0$n, R, ahc = "Y",
  notes = sprintf("Table 3 col 1 spec from do-file; published coef not obtainable (paywalled); author controls peaceyears+3 splines (PEACE_YEARS); tie is monadic node position of ccode1; dropped %d obs vs orig (year restriction + missing benchmarks)", n_orig - R$m0$n))
sink()
