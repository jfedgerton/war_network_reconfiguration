# Rerun 25: Malis 2021 IO -- Table 1 IV col (mod_iv1): 2SLS start_mid_US ~ turnover (instr. lag3enter)
#   + lag4 covariates + prior vacancy cubic + lag4-6 MID onset | FE cow_name + year x prior_pol; cluster cow_name
# Unit: US-country i dyad-years 1960-2010 (US implicit other side, ccode 2).
ROOT <- Sys.getenv("SWEEP_DATA"); OUT <- Sys.getenv("SWEEP_OUT"); CODE <- Sys.getenv("SWEEP_CODE")
stopifnot(nzchar(ROOT), nzchar(OUT), nzchar(CODE))
source(file.path(CODE, "work/rerun_common_b8.R"))
library(fixest)
WO <- file.path(OUT, "work/25"); dir.create(WO, recursive = TRUE, showWarnings = FALSE)
sink(file.path(WO, "log_25.txt"), split = TRUE)

# 1. Load authors' analysis data and apply their sample filter (mid_regressions_main.R)
d <- read.csv(file.path(ROOT, "work/25/Malis_DelegatedDiplomacy_IO_Replication/analysis/data_for_analysis.csv"), stringsAsFactors = FALSE)
d <- subset(d, year %in% 1960:2010 & lag3ineligible == 0 & !is.na(totalvac))
d$us <- 2
TIE <- "fit_turnover"
COV <- paste("lag4logimp_imputed + lag4logexp_imputed + lag4idealdiff_imputed + lag4polity_imputed + lag4deltapolity_imputed",
             "+ lag4cap + lag4ally + totalvac + totalvac2 + totalvac3 + lag4startmid + lag5startmid + lag6startmid")
# factor(cow_name) + factor(year)*factor(prior_pol) == FE cow_name + (year x prior_pol) cells.
# SEs: sandwich::vcovCL on an ivreg object defaults to HC0 meat + G/(G-1) (HC1 only for class "lm");
# coeftest p uses t with residual df (n - K incl. all dummies). Mimic: ssc(adj=FALSE, cluster.adj=TRUE, t.df="conventional", fixef.K="full").
# fixef.rm="none" keeps the one singleton-country obs that lm/ivreg keeps (N = 6,279).
fitter <- function(add, dat) {
  f <- as.formula(paste("start_mid_US ~", COV, add, "| cow_name + year^prior_pol | turnover ~ lag3enter"))
  m <- feols(f, data = dat, cluster = ~cow_name, fixef.rm = "none",
             ssc = ssc(adj = FALSE, fixef.K = "full", cluster.adj = TRUE, t.df = "conventional"))
  print(summary(m)); ct <- coeftable(m)
  stopifnot(TIE %in% rownames(ct))
  list(b = ct[TIE, 1], se = ct[TIE, 2], p = ct[TIE, 4], n = nobs(m))
}

# 2. Reproduce (authors' mid_IV.tex col 3: turnover 0.022 (0.011), p = 0.040, N = 6,279)
cat("==== Reproduction ====\n")
r_pub <- fitter("", d)
b_pub <- 0.022
stopifnot(r_pub$n == 6279, abs(r_pub$b - b_pub) / abs(b_pub) <= 0.10, r_pub$p < 0.05)
n_orig <- r_pub$n

# 3. Merge history (dyad = US(2) x cow_code); years 1960-2010 inside 1822-2011
source(file.path(CODE, "hist/add_history.R"))
d <- add_history(d, "us", "cow_code", "year")
d <- filter(d, year >= 1822, year <= 2011)
cat("\nNA counts in history/benchmark vars:\n"); print(colSums(is.na(d[, HIST_VARS])))
need <- c("start_mid_US", "turnover", "lag3enter", "cow_name", "prior_pol", all.vars(as.formula(paste("~", COV))), HIST_VARS)
d <- d[complete.cases(d[, need]), ]
cat("Rows after year restriction + complete cases:", nrow(d), " (orig model N", n_orig, ")\n")

# 4. Models on the same sample
R <- run_six(fitter, d)

# 5. Write row
write_row(25, "Malis 2021", "REPRODUCED", b_pub, r_pub$b, "2SLS (fixest feols IV), FE country + year x prior_pol, cluster country",
  "turnover (instr. lag3enter)", n_orig, R$m0$n, R, ahc = "Y",
  notes = sprintf("Table 1 IV col 3; author already controls lag4-6 US-i MID onset (LAGGED_CONFLICT) + country FE; hist_dyad = US-i MIDs t-1..t-5 overlaps IV timing (lags 1-3 are post-instrument); tie = ambassador turnover (borderline scope); dropped %d obs vs orig (missing Polity benchmark bm_democ_w15)", n_orig - R$m0$n))
sink()
