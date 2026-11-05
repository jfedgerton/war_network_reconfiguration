# Rerun 5: Fuhrmann & Sechser 2014 AJPS -- Table 2 Model 1
# probit military_conflict ~ defense_pact_nuclear_cow + controls + time cubic, if politically_relevant==1, cluster(dyad_id)
ROOT <- Sys.getenv("SWEEP_DATA"); OUT <- Sys.getenv("SWEEP_OUT"); CODE <- Sys.getenv("SWEEP_CODE")
stopifnot(nzchar(ROOT), nzchar(OUT), nzchar(CODE))
source(file.path(CODE, "work/rerun_common_b8.R"))
library(haven)
WO <- file.path(OUT, "work/5"); dir.create(WO, recursive = TRUE, showWarnings = FALSE)
sink(file.path(WO, "log_5.txt"), split = TRUE)

# 1. Load authors' data (directed dyads: state_a = challenger, state_b = target; 1950-2000)
d <- read_dta(file.path(ROOT, "work/5/Replication Material to Post/Fuhrmann-Sechser-Appendices.dta"))
d <- zap_labels(d)
d <- as.data.frame(d[d$politically_relevant == 1, ])
TIE <- "defense_pact_nuclear_cow"
RHS <- paste("military_conflict ~ defense_pact_nuclear_cow + nuclear_deployment + defense_pact_nonnuclear_cow + us_troops +",
  "challenger_nuclear_weapons + target_nuclear_weapons + contiguity + alliance_with_challenger_cow + foreign_policy_similarity +",
  "power_ratio + challenger_polity + target_polity + polity_interaction + time_conflict + time_conflict2 + time_conflict3")
fitter <- function(add, dat) fit_tie(as.formula(paste(RHS, add)), dat, binomial(link = "probit"), TIE, cluster = "dyad_id")

# 2. Reproduce (published: b = -0.388, se 0.0938, N = 85,306)
cat("==== Reproduction ====\n")
r_pub <- fitter("", d)
b_pub <- -0.388
stopifnot(r_pub$n == 85306)
stopifnot(abs(r_pub$b - b_pub) / abs(b_pub) <= 0.10, r_pub$p < 0.05)
n_orig <- r_pub$n

# 3. Merge history; restrict to 1822-2011 (all obs 1950-2000, so no year loss expected)
source(file.path(CODE, "hist/add_history.R"))
d <- add_history(d, "state_a", "state_b", "year")
d <- filter(d, year >= 1822, year <= 2011)
vars <- c(all.vars(as.formula(RHS)), "dyad_id", HIST_VARS)
cat("\nNA counts in history/benchmark vars:\n"); print(colSums(is.na(d[, HIST_VARS])))
d <- d[complete.cases(d[, vars]), ]
cat("Rows after year restriction + complete cases:", nrow(d), " (orig model N", n_orig, ")\n")

# 4. Models on the same sample
R <- run_six(fitter, d)

# ---- sensitivity hook: runs only under run_sensitivity.R (SWEEP_SENS=1), never in the main sweep ----
if (Sys.getenv("SWEEP_SENS") == "1") {
  source(file.path(CODE, "hist/sens_helper.R"))
  sens_run(id = "5", dat = d, y = "military_conflict", tie = TIE,
           ctrl = c("nuclear_deployment", "defense_pact_nonnuclear_cow", "us_troops", "challenger_nuclear_weapons", "target_nuclear_weapons",
                    "contiguity", "alliance_with_challenger_cow", "foreign_policy_similarity", "power_ratio", "challenger_polity",
                    "target_polity", "polity_interaction", "time_conflict", "time_conflict2", "time_conflict3"), fe = NULL,
           fit = function(dat, add) { r <- fitter(if (length(add)) paste0(" + ", paste(add, collapse = " + ")) else "", dat); c(b = r$b, se = r$se, p = r$p) })
}

# 5. Write row
write_row(5, "Fuhrmann & Sechser 2014", "REPRODUCED", b_pub, r_pub$b, "probit, SE clustered on dyad_id",
  TIE, n_orig, R$m0$n, R, ahc = "Y",
  notes = sprintf("Table 2 M1; author controls cubic peace-time (PEACE_YEARS); DV = violent MID initiation; dropped %d obs vs orig (missing benchmark CINC/Polity lags)", n_orig - R$m0$n))
sink()
