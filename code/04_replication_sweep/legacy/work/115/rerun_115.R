# Rerun 115: Edgerton 2024 JOP -- Table 3 col 1, FE logit fatal_mid ~ out_group + controls
ROOT <- Sys.getenv("SWEEP_DATA"); OUT <- Sys.getenv("SWEEP_OUT"); CODE <- Sys.getenv("SWEEP_CODE")
stopifnot(nzchar(ROOT), nzchar(OUT), nzchar(CODE))
suppressPackageStartupMessages({library(fixest); library(dplyr)})
WO <- file.path(OUT, "work/115"); dir.create(WO, recursive = TRUE, showWarnings = FALSE)
sink(file.path(WO, "log_115.txt"), split = TRUE)

# 1. Load authors' analysis data
war.data <- readRDS(file.path(ROOT, "_parent", "coop/war_data_replication.rds"))
stopifnot(nrow(war.data) == 993975)

# 2. Reproduce: authors' seed 1 defines the 75% training split (published sample)
set.seed(1)
training_rows <- sample(nrow(war.data), round(0.75 * nrow(war.data)))
d <- war.data[training_rows, ]
d <- filter(d, !is.na(fatal_mid))
CL <- ~year + dyad + ccode1 + ccode2
f0 <- fatal_mid ~ out_group + democ_joint + autoc_joint + major_power_joint + cinc_ratio

fe <- "| ccode1 + ccode2 + dyad + year"
fit <- function(rhs_add, dat) {
  fml <- as.formula(paste("fatal_mid ~ out_group + democ_joint + autoc_joint + major_power_joint + cinc_ratio", rhs_add, fe))
  m <- feglm(fml, data = dat, family = "binomial", cluster = CL)
  print(summary(m))
  ct <- coeftable(m)
  list(b = ct["out_group", 1], se = ct["out_group", 2], p = ct["out_group", 4], n = nobs(m))
}
cat("\n==== Reproduction (full training sample) ====\n")
r_pub <- fit("", d)
b_pub <- 0.847
stopifnot(abs(r_pub$b - b_pub) / abs(b_pub) <= 0.10)  # fail loudly if not reproduced
n_orig <- r_pub$n

# 3. Merge history; restrict to 1822-2011

source(file.path(CODE, "hist/add_history.R"))
set.seed(1); tr <- sample(nrow(war.data), round(0.75 * nrow(war.data)))
d <- filter(war.data[tr, ], !is.na(fatal_mid))
d <- add_history(d, "ccode1", "ccode2", "year")
d <- filter(d, year >= 1822, year <= 2011)
d <- d[complete.cases(d[, c("fatal_mid","out_group","democ_joint","autoc_joint","major_power_joint","cinc_ratio",
                            "hist_dyad_w15","hist_mono_w15","hist_dyad_w26","hist_mono_w26","bm_cinc_w15","bm_democ_w15")]), ]
cat("\nRows after year restriction + complete cases:", nrow(d), "\n")

# 4. Models on the same sample
cat("\n==== m0 ====\n");  r0  <- fit("", d)
cat("\n==== mD ====\n");  rD  <- fit("+ hist_dyad_w15", d)
cat("\n==== mDM ====\n"); rDM <- fit("+ hist_dyad_w15 + hist_mono_w15", d)
cat("\n==== m26 ====\n"); r26 <- fit("+ hist_dyad_w26 + hist_mono_w26", d)
cat("\n==== bC ====\n");  rC  <- fit("+ bm_cinc_w15", d)
cat("\n==== bD ====\n");  rBD <- fit("+ bm_democ_w15", d)
# ---- sensitivity hook: runs only under run_sensitivity.R (SWEEP_SENS=1), never in the main sweep ----
if (Sys.getenv("SWEEP_SENS") == "1") {
  source(file.path(CODE, "hist/sens_helper.R"))
  sens_run(id = "115", dat = d, y = "fatal_mid", tie = "out_group", ctrl = c("democ_joint", "autoc_joint", "major_power_joint", "cinc_ratio"),
           fe = c("ccode1", "ccode2", "dyad", "year"), fe_sep = TRUE,
           fit = function(dat, add) { r <- fit(if (length(add)) paste("+", paste(add, collapse = " + ")) else "", dat); c(b = r$b, se = r$se, p = r$p) })
}
pct <- function(x) 100 * (x$b - r0$b) / abs(r0$b)

out <- data.frame(cid = 115, short_cite = "Edgerton 2024", reproduced = "REPRODUCED", b_pub = b_pub, b_repro = r_pub$b,
  estimator = "fixest feglm logit, FE ccode1+ccode2+dyad+year, 4-way clustered", tie_var = "out_group",
  n_orig = n_orig, n_used = r0$n, b0 = r0$b, se0 = r0$se, p0 = r0$p, b_dm = rDM$b, se_dm = rDM$se, p_dm = rDM$p,
  pct_dm = pct(rDM), b_d = rD$b, pct_d = pct(rD), b_26 = r26$b, se_26 = r26$se, p_26 = r26$p, pct_26 = pct(r26),
  pct_bC = pct(rC), pct_bD = pct(rBD), author_has_history_control = "N",
  notes = sprintf("75%% training split with author seed 1; nobs counts exclude FE-perfectly-predicted obs; dropped %d obs vs orig (year restriction 1822-2011 + FE separation). fixest reports cluster as ccode1&dyad&ccode1&ccode2: year cluster silently dropped (identical to author code SEs). b_pub from task spec (Table 3 col 1).", n_orig - r0$n))
print(Filter(is.numeric, out))   # port: no citation/text columns in logs
dir.create(file.path(OUT, "results"), recursive = TRUE, showWarnings = FALSE)
write.csv(out, file.path(OUT, "results/115.csv"), row.names = FALSE)
sink()
