# [2026-09-26] Headline system-level DV switched from KL to Jensen-Shannon (JS)
# divergence. KL results are retained in results/appendix/kl_canonical/.
###############################################################################
# 01_system_level_analysis.R
# System-level analysis: War → Network Reconfiguration
#
# DV: KL divergence (main). Alt DVs → appendix script.
#
# 4 IVs (counts of ongoing conflicts, raw + log):
#   (a) mid_count, (b) log1p(mid_count),
#   (c) fatal_mid_count, (d) log1p(fatal_mid_count)
#
# Onset measures dropped from system level: ongoing wars continue to
# reshape network topology, so cumulative counts are the appropriate
# IV. Onset is used at the state level where it captures the discrete
# shock to a specific state's network.
# Log(1+x) dampens influence of outlier years (WWI, WWII).
#
# Structure:
#   1. OLS with controls for each IV
#   2. VAR / Granger causality for each IV
#   3. Forecast Error Variance Decomposition (FEVD) for each IV
#
# Appendix (09_appendix_irf.R): IRFs and Local Projections
#
# Controls (full models): n_democracies, n_autocracies (log1p)
###############################################################################

library(tidyverse)
library(vars)
library(strucchange)
library(tseries)
library(fixest)
library(sandwich)
library(texreg)

select <- dplyr::select

OUTPUT_DIR <- here::here("data", "processed")
FIG_DIR <- here::here("code", "03_figures")
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

# --------------------------------------------------------------------------- #
# 0. Load data
# --------------------------------------------------------------------------- #
set.seed(123)
system_data <- readRDS(file.path(OUTPUT_DIR, "system_analysis_data.rds"))

# --------------------------------------------------------------------------- #
# 1. Prepare time series
# --------------------------------------------------------------------------- #

ts_data <- system_data |>
  filter(!is.na(js_divergence)) |>
  select(
    year,
    kl_divergence,
    js_divergence,
    edge_turnover_pct,
    density_change,
    centralization_change,
    modularity_change,
    degree_kl,
    degree_wasserstein,
    mid_count,
    fatal_mid_count,
    any_of(c("n_states", "n_democracies", "n_autocracies"))
  ) |>
  arrange(year) |>
  mutate(
    log_js_divergence = log1p(js_divergence),
    log_mid_count = log1p(mid_count),
    log_fatal_mid_count = log1p(fatal_mid_count),
    # Lagged IVs for OLS (avoids simultaneity)
    mid_count_lag1 = lag(mid_count),
    log_mid_count_lag1 = lag(log_mid_count),
    fatal_mid_count_lag1 = lag(fatal_mid_count),
    log_fatal_mid_count_lag1 = lag(log_fatal_mid_count),
    # AR(1) lagged DV — used in the main-table dynamics column
    log_js_divergence_lag1 = lag(log_js_divergence),
    # Period FE: 4 international-order regimes (matches 12_order_period_fe.R)
    order_period = factor(case_when(
      year <= 1913 ~ "Concert of Europe",
      year <= 1945 ~ "WWI / Interwar / WWII",
      year <= 1991 ~ "Cold War",
      TRUE         ~ "Post-Cold War"
    ),
    levels = c("Concert of Europe", "WWI / Interwar / WWII",
               "Cold War", "Post-Cold War")),
    # `unit` column required by feols panel.id when computing NW HAC SEs
    # on a single time series (treats it as a panel of one unit).
    unit = 1L
  )

# IVs used inside VAR (contemporaneous)
war_vars <- c("mid_count", "log_mid_count", "fatal_mid_count", "log_fatal_mid_count")
# IVs used inside OLS (lagged to avoid simultaneity)
ols_ivs <- c("mid_count_lag1", "log_mid_count_lag1",
             "fatal_mid_count_lag1", "log_fatal_mid_count_lag1")

cat("Stationarity tests (ADF):\n")
for (var in c("js_divergence", war_vars)) {
  adf_result <- adf.test(na.omit(ts_data[[var]]))
  cat(sprintf("  %s: ADF stat = %.3f, p = %.4f\n", var, adf_result$statistic, adf_result$p.value))
}

# =========================================================================== #
# MAIN MODELS — DV: KL DIVERGENCE
# =========================================================================== #
#
# OLS: uses lagged IVs (t-1) to avoid simultaneity.
# VAR/Granger/FEVD: uses contemporaneous IVs (VAR handles endogeneity
#   internally via its own lag structure).
# =========================================================================== #

CONTROL_TERMS_LOG <- c("log1p(n_democracies)",
                        "log1p(n_autocracies)",
                        "log1p(n_states)")

# ----------------------------------------------------------------------- #
# OLS: 2 sets x 2 controls x 2 IVs = 8 specs
# ----------------------------------------------------------------------- #
cat("\n========== OLS (lagged IVs, t-1) ==========\n")

build_ols_formula <- function(set, controls, iv) {
  dv <- if (set == "logged") "log_js_divergence" else "js_divergence"
  iv_var <- if (set == "logged") {
    if (iv == "mid") "log_mid_count_lag1" else "log_fatal_mid_count_lag1"
  } else {
    if (iv == "mid") "mid_count_lag1" else "fatal_mid_count_lag1"
  }
  parts <- iv_var
  if (controls == "full") parts <- c(parts, CONTROL_TERMS_LOG)
  parts <- c(parts, "poly(year, 2)")
  as.formula(paste(dv, "~", paste(parts, collapse = " + ")))
}

ols_grid <- expand.grid(
  set      = c("logged", "raw"),
  controls = c("none", "full"),
  iv       = c("mid", "fatal"),
  stringsAsFactors = FALSE
)

ols_models <- list()
for (i in seq_len(nrow(ols_grid))) {
  spec <- ols_grid[i, ]
  key  <- sprintf("ols__%s__%s__%s", spec$set, spec$controls, spec$iv)
  fml  <- build_ols_formula(spec$set, spec$controls, spec$iv)
  fit  <- lm(fml, data = ts_data)
  ols_models[[key]] <- list(
    set = spec$set, controls = spec$controls, iv = spec$iv,
    formula = fml, fit = fit
  )
  cat(sprintf("\nOLS [%s | %s | %s]\n  %s\n",
              spec$set, spec$controls, spec$iv, deparse1(fml)))
  print(summary(fit))
}

# --- VAR / Granger / FEVD with contemporaneous IVs ---
#
# TRANSFORMATION: js_divergence is log1p-transformed before entering the
# VAR (LHS only). KL is non-negative and right-skewed, and right-skewness
# on the LHS drives instability of IRF / FEVD coefficients via outlier
# leverage. Logging fixes that and lets IRF coefficients read as
# approximate elasticities -- the natural unit for a bidirectional
# reciprocal-dynamics comparison alongside script 05_svar.R.
#
# The OLS section above (lagged-IV main-text spec) is intentionally NOT
# transformed. OLS slope estimates are unbiased under a skewed LHS, the
# four-IV ladder already includes raw and log(1+x) variants on the IV
# side, and changing the OLS LHS would force a refactor of every main-
# text coefficient table for no analytical gain.
#
# The 4-IV sweep (mid / log_mid / fatal / log_fatal on the IV side) is
# preserved here. Variants with raw IVs estimate the response of log(KL)
# to a unit change in MID counts (semi-elasticity); variants with log IVs
# estimate the elasticity. Both are interpretable; the writeup flags the
# unit difference.

# ----------------------------------------------------------------------- #
# VAR / Granger / FEVD: 2 sets x 2 controls x 2 IVs = 8 specs
# ----------------------------------------------------------------------- #
cat("\n========== VAR / Granger / FEVD ==========\n")

build_var_inputs <- function(set, controls, iv) {
  endo_iv <- if (set == "logged") {
    if (iv == "mid") "log_mid_count" else "log_fatal_mid_count"
  } else {
    if (iv == "mid") "mid_count" else "fatal_mid_count"
  }
  endo_dv <- if (set == "logged") "log_js_divergence" else "js_divergence"

  cols_needed <- c(endo_iv, endo_dv)
  if (controls == "full") {
    cols_needed <- c(cols_needed,
                      "n_democracies", "n_autocracies", "n_states")
  }
  blk <- ts_data |> select(all_of(cols_needed)) |> na.omit()

  endo <- blk |> select(all_of(c(endo_iv, endo_dv))) |> as.ts()
  exog <- if (controls == "full") {
    blk |>
      select(n_democracies, n_autocracies, n_states) |>
      mutate_all(log1p) |>
      as.ts()
  } else {
    NULL
  }
  list(endo = endo, exog = exog,
       endo_iv = endo_iv, endo_dv = endo_dv)
}

var_models <- list()
granger_results_list <- list()
fevd_results_list <- list()

var_grid <- expand.grid(
  set      = c("logged", "raw"),
  controls = c("none", "full"),
  iv       = c("mid", "fatal"),
  stringsAsFactors = FALSE
)

for (i in seq_len(nrow(var_grid))) {
  spec <- var_grid[i, ]
  key  <- sprintf("var__%s__%s__%s", spec$set, spec$controls, spec$iv)

  cat(sprintf("\n%s\n", strrep("=", 70)))
  cat(sprintf("  VAR [%s | %s | %s]\n", spec$set, spec$controls, spec$iv))
  cat(sprintf("%s\n", strrep("=", 70)))

  inp <- build_var_inputs(spec$set, spec$controls, spec$iv)

  optimal_lag <- tryCatch({
    sel_args <- list(y = inp$endo, lag.max = 10, type = "both")
    if (!is.null(inp$exog)) sel_args$exogen <- inp$exog
    lag_select <- do.call(VARselect, sel_args)
    cat(sprintf("  Optimal lag (AIC): %d, (BIC): %d\n",
                lag_select$selection["AIC(n)"],
                lag_select$selection["SC(n)"]))
    max(1L, as.integer(lag_select$selection["AIC(n)"]))
  }, error = function(e) {
    cat("  VARselect failed:", conditionMessage(e), "; defaulting to 2\n")
    2L
  })

  var_args <- list(y = inp$endo, p = as.integer(optimal_lag), type = "both")
  if (!is.null(inp$exog)) var_args$exogen <- inp$exog
  var_fit <- do.call(VAR, var_args)
  var_models[[key]] <- list(
    set = spec$set, controls = spec$controls, iv = spec$iv,
    endo_iv = inp$endo_iv, endo_dv = inp$endo_dv,
    fit = var_fit, lag = optimal_lag
  )
  cat(sprintf("\nVAR(%d) summary (%s equation):\n",
              optimal_lag, inp$endo_dv))
  print(summary(var_fit$varresult[[inp$endo_dv]]))

  # Granger both directions
  gc_to_kl <- causality(var_fit, cause = inp$endo_iv, boot.runs = 1000)
  gc_kl_to <- causality(var_fit, cause = inp$endo_dv, boot.runs = 1000)
  cat(sprintf("\nGranger: %s -> %s:  F=%.3f, p=%.4f\n",
              inp$endo_iv, inp$endo_dv,
              gc_to_kl$Granger$statistic, gc_to_kl$Granger$p.value))
  cat(sprintf("Granger: %s -> %s:  F=%.3f, p=%.4f\n",
              inp$endo_dv, inp$endo_iv,
              gc_kl_to$Granger$statistic, gc_kl_to$Granger$p.value))

  granger_results_list[[key]] <- list(
    set = spec$set, controls = spec$controls, iv = spec$iv,
    endo_iv = inp$endo_iv, endo_dv = inp$endo_dv,
    fwd = gc_to_kl, rev = gc_kl_to
  )

  # FEVD
  fevd_fit <- fevd(var_fit, n.ahead = 10)
  fevd_dv <- as.data.frame(fevd_fit[[inp$endo_dv]])
  fevd_dv$horizon  <- 1:nrow(fevd_dv)
  fevd_dv$set      <- spec$set
  fevd_dv$controls <- spec$controls
  fevd_dv$iv       <- spec$iv
  fevd_dv$endo_iv  <- inp$endo_iv
  fevd_dv$endo_dv  <- inp$endo_dv
  fevd_results_list[[key]] <- fevd_dv

  cat(sprintf("\nFEVD: %s variance explained by %s shocks\n",
              inp$endo_dv, inp$endo_iv))
  print(fevd_dv)
}

# --------------------------------------------------------------------------- #
# Save results
# --------------------------------------------------------------------------- #

cat("\nSaving results...\n")

results <- list(
  ols = ols_models,
  var = var_models,
  granger = granger_results_list,
  fevd = fevd_results_list,
  system_data = ts_data,
  war_vars = war_vars,
  ols_ivs = ols_ivs
)

saveRDS(results, file.path(OUTPUT_DIR, "system_level_results.rds"))

# --- CSV exports (new 2 sets x 2 controls x 2 IVs structure) ---
write_csv(ts_data, file.path(OUTPUT_DIR, "system_analysis_data.csv"))

# OLS coefficients: one combined CSV with set / controls / iv columns.
# Plain OLS (iid) standard errors. SE robustness comparisons (HC1, HC3, HAC,
# bootstrap) live exclusively in 19_appendix_hac.R.
ols_coefs <- list()
for (key in names(ols_models)) {
  m <- ols_models[[key]]
  sm <- summary(m$fit)$coefficients
  ols_coefs[[key]] <- tibble(
    set      = m$set,
    controls = m$controls,
    iv       = m$iv,
    spec_key = key,
    term     = rownames(sm),
    estimate = sm[, "Estimate"],
    se       = sm[, "Std. Error"],
    t_value  = sm[, "t value"],
    p_value  = sm[, "Pr(>|t|)"]
  )
}
ols_coefs_df <- bind_rows(ols_coefs)
write_csv(ols_coefs_df,
          file.path(OUTPUT_DIR, "system_ols_coefs_combined.csv"))
# Convenience splits (one CSV per set) for the manuscript tables.
write_csv(ols_coefs_df |> filter(set == "logged"),
          file.path(OUTPUT_DIR, "system_ols_coefs_logged.csv"))
write_csv(ols_coefs_df |> filter(set == "raw"),
          file.path(OUTPUT_DIR, "system_ols_coefs_raw.csv"))

# Granger summary: one combined CSV; columns set / controls / iv plus
# F-stats and p-values for both directions.
granger_rows <- list()
for (key in names(granger_results_list)) {
  g <- granger_results_list[[key]]
  granger_rows[[key]] <- tibble(
    set      = g$set,
    controls = g$controls,
    iv       = g$iv,
    spec_key = key,
    endo_iv  = g$endo_iv,
    endo_dv  = g$endo_dv,
    fwd_F    = as.numeric(g$fwd$Granger$statistic),
    fwd_p    = as.numeric(g$fwd$Granger$p.value),
    rev_F    = as.numeric(g$rev$Granger$statistic),
    rev_p    = as.numeric(g$rev$Granger$p.value)
  )
}
granger_summary <- bind_rows(granger_rows)
write_csv(granger_summary,
          file.path(OUTPUT_DIR, "granger_results_combined.csv"))
write_csv(granger_summary |> filter(set == "logged"),
          file.path(OUTPUT_DIR, "granger_results_logged.csv"))
write_csv(granger_summary |> filter(set == "raw"),
          file.path(OUTPUT_DIR, "granger_results_raw.csv"))

# FEVD: one combined CSV.
fevd_combined <- bind_rows(fevd_results_list)
write_csv(fevd_combined, file.path(OUTPUT_DIR, "fevd_results_combined.csv"))
write_csv(fevd_combined |> filter(set == "logged"),
          file.path(OUTPUT_DIR, "fevd_results_logged.csv"))
write_csv(fevd_combined |> filter(set == "raw"),
          file.path(OUTPUT_DIR, "fevd_results_raw.csv"))

# =========================================================================== #
# LATEX TABLES (texreg)
# =========================================================================== #

TABLE_DIR <- here::here("results", "main")
dir.create(TABLE_DIR, recursive = TRUE, showWarnings = FALSE)

# Legacy coefficient map — retained as a reference for the simple-OLS
# CSV exports above. The 4-col main table uses ols_coef_map_4col defined
# below (which keeps the AR(1) lagged DV row).
ols_coef_map <- list(
  "mid_count_lag1"              = "Conflict measure ($t{-}1$)",
  "log_mid_count_lag1"          = "Conflict measure ($t{-}1$)",
  "fatal_mid_count_lag1"        = "Conflict measure ($t{-}1$)",
  "log_fatal_mid_count_lag1"    = "Conflict measure ($t{-}1$)",
  "log1p(n_democracies)"        = "$\\log(1 + \\text{Democracies})$",
  "log1p(n_autocracies)"        = "$\\log(1 + \\text{Autocracies})$",
  "log1p(n_states)"             = "$\\log(1 + \\text{States})$",
  "(Intercept)"                 = "Intercept"
)

# =========================================================================== #
# MAIN-TEXT 4-COLUMN SYSTEM TABLE (HAC inference throughout)
# =========================================================================== #
# Anchored on log(1+KL) ~ log(1+MID count_lag1) — same DV/IV pairing in
# every column. Each column adds one identification challenge:
#   (1) bivariate                          — unconditional association
#   (2) + full controls                    — observable confounders + smooth
#                                            time trend
#   (3) + AR(1) lagged DV                  — conditional on KL momentum
#   (4) + order-period FE                  — within-regime identification
#
# All four columns use Newey-West HAC SEs (lag = 5) for consistent
# inference. Cluster-on-period (G=4) is too few clusters for reliable
# cluster inference; HAC handles within-period serial correlation
# without that pitfall.
#
# Same 4-col structure for the fatal-MID appendix companion.
# --------------------------------------------------------------------------- #

cat("\n========== Main system table: 4-col HAC ladder ==========\n")

NW_LAG <- 5L
SYS_CTRL_FULL <- "log1p(n_democracies) + log1p(n_autocracies) + log1p(n_states) + poly(year, 2)"

# Drop the 1 obs lost to the AR(1) lag so all 4 cols are on the same sample.
# (Without this, col 3 would have n=196 while others have n=197 — readers
# don't appreciate the inconsistency, and it's a costless harmonization.)
ts_data_main <- ts_data |> filter(!is.na(log_js_divergence_lag1),
                                    !is.na(mid_count_lag1))

fit_main_table_4col <- function(iv_family) {
  stopifnot(iv_family %in% c("mid", "fatal"))
  iv_var <- if (iv_family == "mid") "log_mid_count_lag1"
                                     else "log_fatal_mid_count_lag1"

  # Col 1: IV + quadratic year trend (bivariate-with-trend baseline).
  # Including `poly(year, 2)` here keeps col 1 honest — with 1816--2014
  # annual data and clear secular drift in both KL and MID, a strict
  # bivariate would conflate the MID effect with the trend. All four
  # columns therefore share `poly(year, 2)` as a baseline; the ladder
  # varies by what's added on top.
  m1 <- feols(as.formula(paste("log_js_divergence ~", iv_var,
                                 "+ poly(year, 2)")),
              data = ts_data_main,
              panel.id = ~unit + year, vcov = NW(lag = NW_LAG))

  # Col 2: + full controls (adds n_dem, n_aut, n_states logged on top of
  # the year trend already in col 1).
  m2 <- feols(as.formula(paste("log_js_divergence ~", iv_var, "+",
                                 SYS_CTRL_FULL)),
              data = ts_data_main,
              panel.id = ~unit + year, vcov = NW(lag = NW_LAG))

  # Col 3: + AR(1) lagged DV
  m3 <- feols(as.formula(paste("log_js_divergence ~", iv_var, "+",
                                 SYS_CTRL_FULL,
                                 "+ log_js_divergence_lag1")),
              data = ts_data_main,
              panel.id = ~unit + year, vcov = NW(lag = NW_LAG))

  # Col 4: + order-period FE (passed via | order_period in feols formula)
  m4 <- feols(as.formula(paste("log_js_divergence ~", iv_var, "+",
                                 SYS_CTRL_FULL,
                                 "| order_period")),
              data = ts_data_main,
              panel.id = ~unit + year, vcov = NW(lag = NW_LAG))

  list(m1 = m1, m2 = m2, m3 = m3, m4 = m4)
}

main_4col_mid   <- fit_main_table_4col("mid")
main_4col_fatal <- fit_main_table_4col("fatal")

cat("\n  MID 4-col estimates (HAC NW(5)):\n")
for (k in names(main_4col_mid)) {
  m <- main_4col_mid[[k]]
  iv <- "log_mid_count_lag1"
  if (iv %in% rownames(coeftable(m))) {
    ct <- coeftable(m)
    cat(sprintf("    %s: beta=%.4f  HAC SE=%.4f  p=%.4g\n",
                k, ct[iv, "Estimate"], ct[iv, "Std. Error"],
                ct[iv, "Pr(>|t|)"]))
  }
}

cat("\n  FATAL 4-col estimates (HAC NW(5)):\n")
for (k in names(main_4col_fatal)) {
  m <- main_4col_fatal[[k]]
  iv <- "log_fatal_mid_count_lag1"
  if (iv %in% rownames(coeftable(m))) {
    ct <- coeftable(m)
    cat(sprintf("    %s: beta=%.4f  HAC SE=%.4f  p=%.4g\n",
                k, ct[iv, "Estimate"], ct[iv, "Std. Error"],
                ct[iv, "Pr(>|t|)"]))
  }
}

# --- texreg helper: 4-col main table, MID or fatal ---
# Coefficient map. Both IV variants map to one "Conflict measure" row;
# the AR(1) lagged DV gets its own row that only shows in col 3.
ols_coef_map_4col <- list(
  "log_mid_count_lag1"           = "$\\log(1{+}\\text{MID Count})_{t{-}1}$",
  "log_fatal_mid_count_lag1"     = "$\\log(1{+}\\text{Fatal MID Count})_{t{-}1}$",
  "log1p(n_democracies)"         = "$\\log(1{+}\\text{Democracies})$",
  "log1p(n_autocracies)"         = "$\\log(1{+}\\text{Autocracies})$",
  "log1p(n_states)"              = "$\\log(1{+}\\text{States})$",
  "log_js_divergence_lag1"       = "$\\log(1{+}\\text{JS})_{t{-}1}$ (AR(1))",
  "(Intercept)"                  = "Intercept"
)

write_main_table_4col <- function(models_list, iv_family,
                                    file_path, table_dir_label) {
  stopifnot(iv_family %in% c("mid", "fatal"))
  iv_short <- if (iv_family == "mid")
    "$\\log(1{+}\\text{MID})$"
  else
    "$\\log(1{+}\\text{Fatal MID})$"

  caption_lead <- if (iv_family == "mid")
    paste0("System-Level OLS: $\\log(1+\\text{JS})$ on ", iv_short,
           " (HAC inference throughout). All four columns include a ",
           "quadratic year trend. Each subsequent column tightens one ",
           "identification challenge: column 2 adds observable controls; ",
           "column 3 additionally controls for JS momentum (AR(1) lagged ",
           "DV); column 4 absorbs four international-order regimes via ",
           "fixed effects.")
  else
    paste0("Appendix. System-Level OLS: $\\log(1+\\text{JS})$ on ",
           iv_short, " (companion to main-text Table~\\ref{tab:system_main}).")

  label_suffix <- if (iv_family == "mid") "system_main" else "system_appx_fatal"

  texreg(
    list(models_list$m1, models_list$m2, models_list$m3, models_list$m4),
    custom.model.names = c("(1) IV + trend",
                            "(2) + Full ctrls",
                            "(3) + AR(1)",
                            "(4) + Period FE"),
    custom.coef.map = ols_coef_map_4col,
    caption = caption_lead,
    label   = sprintf("tab:%s", label_suffix),
    custom.note = paste0(
      "\\parbox{\\linewidth}{\\vspace{2pt}\\textit{Note:} ",
      "DV is $\\log(1+\\text{JS})$. IV is ", iv_short, " lagged one ",
      "period ($t-1$). Newey-West HAC standard errors (lag $=$ 5) in ",
      "parentheses. \\textbf{All four columns include a quadratic year ",
      "trend (\\texttt{poly(year, 2)}); terms suppressed for legibility.} ",
      "Column 1 is the IV plus the year trend (no other controls). ",
      "Columns 2--4 add the full controls covariate set ",
      "($\\log(1+\\text{Democracies})$, $\\log(1+\\text{Autocracies})$, ",
      "$\\log(1+\\text{States})$). Column 3 additionally controls for ",
      "$\\log(1+\\text{JS})_{t-1}$ (AR(1) term). Column 4 absorbs ",
      "fixed effects for four international-order regimes (Concert of ",
      "Europe / WWI--Interwar--WWII / Cold War / Post-Cold War). ",
      "COW system, 1816--2014. ",
      "$^{***}p<0.001$; $^{**}p<0.01$; $^{*}p<0.05$; ",
      "$^{\\dagger}p<0.1$.}"
    ),
    stars = c(0.001, 0.01, 0.05, 0.1),
    symbol = "\\dagger",
    digits = 3,  # JS divergence is on a [0, log 2] scale; 2 digits hides the estimates
    float.pos = "htbp",
    booktabs = TRUE,
    use.packages = FALSE,
    include.rsquared = TRUE,
    include.adjrs = FALSE,
    include.nobs = TRUE,
    file = file_path
  )
  cat(sprintf("  System main 4-col [%s] -> %s (%s)\n",
              iv_family, basename(file_path), table_dir_label))
}

# Main text: MID only. Appendix: Fatal MID.
write_main_table_4col(main_4col_mid, "mid",
                        file.path(TABLE_DIR, "table1_system_main.tex"),
                        "main")
write_main_table_4col(main_4col_fatal, "fatal",
                        file.path(here::here("results", "appendix"),
                                  "tableA_system_fatal.tex"),
                        "appendix")

# --- Granger + FEVD tables (one per set) -----------------------------

star_fn <- function(p) {
  if (is.na(p)) return("")
  if (p < 0.001) return("$^{***}$")
  if (p < 0.01)  return("$^{**}$")
  if (p < 0.05)  return("$^{*}$")
  if (p < 0.10)  return("$^{\\dagger}$")
  return("")
}

# --- Granger/FEVD table by IV family (MID or fatal MID) -------
# 2 columns, BOTH endogenous variables logged to mirror Table 1's
# anchoring on log(1+KL) ~ log(1+MID count). Cols are no-controls vs.
# full-controls VARs.
#
# Three panels:
#   A: Forward Granger (Conflict -> KL)
#   B: Reverse Granger (KL -> Conflict)  [shows asymmetry]
#   C: FEVD share of KL variance attributable to conflict shocks at
#      multiple horizons (h = 1, 3, 5, 10) — dynamic profile, not just one
#      number.
#
# Raw-scale variants from the legacy table are dropped; they remain
# available in the underlying granger_results_combined.csv /
# fevd_results_combined.csv files for any reader who wants them.
write_granger_fevd_by_iv <- function(iv_family, file_path, table_dir_label) {
  stopifnot(iv_family %in% c("mid", "fatal"))
  col_specs <- list(
    list(controls = "none"),
    list(controls = "full")
  )

  iv_short <- if (iv_family == "mid")
    "$\\log(1{+}\\text{MID})$"
  else
    "$\\log(1{+}\\text{Fatal MID})$"

  col_headers <- vapply(col_specs, function(s) {
    sprintf("%s, %s", iv_short,
            ifelse(s$controls == "none", "no ctrls", "full ctrls"))
  }, character(1))

  fwd_F  <- fwd_p  <- rev_F  <- rev_p  <- character(length(col_specs))
  fevd_h <- list(h1 = character(length(col_specs)),
                  h3 = character(length(col_specs)),
                  h5 = character(length(col_specs)),
                  h10 = character(length(col_specs)))

  for (j in seq_along(col_specs)) {
    s   <- col_specs[[j]]
    key <- sprintf("var__logged__%s__%s", s$controls, iv_family)
    g   <- granger_results_list[[key]]
    fF  <- as.numeric(g$fwd$Granger$statistic)
    fP  <- as.numeric(g$fwd$Granger$p.value)
    rF  <- as.numeric(g$rev$Granger$statistic)
    rP  <- as.numeric(g$rev$Granger$p.value)
    fwd_F[j] <- sprintf("%.3f", fF)
    fwd_p[j] <- sprintf("%.4f%s", fP, star_fn(fP))
    rev_F[j] <- sprintf("%.3f", rF)
    rev_p[j] <- sprintf("%.4f%s", rP, star_fn(rP))

    fd <- fevd_results_list[[key]]
    fevd_h$h1[j]  <- sprintf("%.1f\\%%", fd[[g$endo_iv]][1]  * 100)
    fevd_h$h3[j]  <- sprintf("%.1f\\%%", fd[[g$endo_iv]][3]  * 100)
    fevd_h$h5[j]  <- sprintf("%.1f\\%%", fd[[g$endo_iv]][5]  * 100)
    fevd_h$h10[j] <- sprintf("%.1f\\%%", fd[[g$endo_iv]][10] * 100)
  }

  caption_lead <- if (iv_family == "mid")
    paste0("Granger Causality and FEVD: $\\log(1+\\text{MID Count})$ ",
           "and $\\log(1+\\text{JS Divergence})$. ",
           "VAR with both endogenous variables logged (mirrors the ",
           "main-text Table~\\ref{tab:system_main} anchoring).")
  else
    paste0("Appendix. Granger Causality and FEVD: ",
           "$\\log(1+\\text{Fatal MID Count})$ and ",
           "$\\log(1+\\text{JS Divergence})$. ",
           "Companion to main-text Table~\\ref{tab:granger_fevd_mid}.")

  tex_lines <- c(
    "\\begin{table}[htbp]",
    "\\centering",
    sprintf("\\caption{%s}", caption_lead),
    sprintf("\\label{tab:granger_fevd_%s}", iv_family),
    "\\begin{tabular}{lcc}",
    "\\toprule",
    paste0(" & ", paste(col_headers, collapse = " & "), " \\\\"),
    "\\midrule",
    "\\multicolumn{3}{l}{\\textit{Panel A: Granger Causality (Conflict $\\rightarrow$ JS)}} \\\\[3pt]",
    paste0("$F$-statistic & ", paste(fwd_F, collapse = " & "), " \\\\"),
    paste0("$p$-value & ",     paste(fwd_p, collapse = " & "), " \\\\[6pt]"),
    "\\multicolumn{3}{l}{\\textit{Panel B: Granger Causality (JS $\\rightarrow$ Conflict)}} \\\\[3pt]",
    paste0("$F$-statistic & ", paste(rev_F, collapse = " & "), " \\\\"),
    paste0("$p$-value & ",     paste(rev_p, collapse = " & "), " \\\\[6pt]"),
    "\\multicolumn{3}{l}{\\textit{Panel C: FEVD --- \\% of $\\log(1+\\text{JS})$ variance from conflict shocks}} \\\\[3pt]",
    paste0("$h = 1$ & ",  paste(fevd_h$h1,  collapse = " & "), " \\\\"),
    paste0("$h = 3$ & ",  paste(fevd_h$h3,  collapse = " & "), " \\\\"),
    paste0("$h = 5$ & ",  paste(fevd_h$h5,  collapse = " & "), " \\\\"),
    paste0("$h = 10$ & ", paste(fevd_h$h10, collapse = " & "), " \\\\"),
    "\\bottomrule",
    "\\end{tabular}",
    "\\par\\vspace{4pt}",
    paste0(
      "{\\footnotesize \\textit{Note:} Granger tests from VAR with ",
      "AIC-selected lag order and 1,000 bootstrap replications. Both ",
      "endogenous variables logged ($\\log(1+\\text{JS})$ and ",
      iv_short, ") to match the main-text Table~",
      "\\ref{tab:system_main} anchoring. ``no ctrls'' VARs include only ",
      "constant and trend; ``full ctrls'' VARs additionally include ",
      "exogenous $\\log(1+\\text{Democracies})$, ",
      "$\\log(1+\\text{Autocracies})$, $\\log(1+\\text{States})$. Panel C ",
      "reports the share of forecast variance in $\\log(1+\\text{JS})$ ",
      "attributable to conflict shocks at horizons $h = 1, 3, 5, 10$ years. ",
      "Raw-scale variants reported in the underlying ",
      "\\texttt{granger\\_results\\_combined.csv} and ",
      "\\texttt{fevd\\_results\\_combined.csv}. COW system, 1816--2014. ",
      "$^{***}p<0.001$; $^{**}p<0.01$; $^{*}p<0.05$; $^{\\dagger}p<0.1$.}"
    ),
    "\\end{table}"
  )

  writeLines(tex_lines, file_path)
  cat(sprintf("  Granger + FEVD [%s] -> %s (%s)\n",
              iv_family, basename(file_path), table_dir_label))
}

# Main text: MID only. Fatal-MID Granger/FEVD appendix companion is
# subsumed by 03_robustness.R's multiverse (VAR_Granger forward + reverse
# for fatal IVs and FEVD computed across horizons 1-15) and by 03's
# fig_granger_forward_vs_reverse.png and fig_fevd_horizon.png. Cite those
# in the paper instead of generating a redundant standalone fatal table.
write_granger_fevd_by_iv("mid",
                          file.path(TABLE_DIR, "table2_granger_fevd_mid.tex"),
                          "main")

# =========================================================================== #
# FEVD Decomposition Plot — REMOVED (subsumed by 03_robustness.R)
# =========================================================================== #
# 03's `fig_fevd_horizon.png` shows FEVD horizon profiles for all main
# specs in one figure. The standalone logged/raw FEVD plots that lived
# here were redundant. The underlying `fevd_results_combined.csv` is
# still saved (line ~367) for any reader who wants to plot this themselves.
cat("\nFEVD decomposition plots: NOT generated by 01 (see 03's ",
    "fig_fevd_horizon.png).\n", sep = "")

# =========================================================================== #
# SUMMARY STATISTICS TABLE
# =========================================================================== #
cat("\n========== Summary Statistics ==========\n")

# Aligned with the variables that appear in the main-text Tables 1 and 2
# full-control regressions: DV (KL divergence, raw + logged), IV
# (MID count, raw + logged; fatal-MID variants for the appendix
# companion), and the exogenous controls (n_democracies, n_autocracies,
# n_states — reported on the raw scale; regressions use log1p of these).
# Alt DVs (JS, edge turnover, degree KL/Wasserstein) and lagged versions
# of IVs / DV are dropped — they appear only in 03's multiverse, not in
# the main-text regressions.
sumstat_vars <- c(
  # Network outcomes (DV in Tables 1 and 2)
  "js_divergence", "log_js_divergence",
  # Conflict measures (IV in Tables 1 and 2; fatal in appendix companion)
  "mid_count", "log_mid_count",
  "fatal_mid_count", "log_fatal_mid_count",
  # Exogenous controls (raw form for interpretability; regressions use log1p)
  "n_democracies", "n_autocracies", "n_states"
)

sumstat_labels <- c(
  js_divergence       = "JS Divergence",
  log_js_divergence   = "$\\log(1 + \\text{JS Divergence})$",
  mid_count           = "MID Count",
  log_mid_count       = "$\\log(1 + \\text{MID Count})$",
  fatal_mid_count     = "Fatal MID Count",
  log_fatal_mid_count = "$\\log(1 + \\text{Fatal MID Count})$",
  n_democracies       = "N Democracies",
  n_autocracies       = "N Autocracies",
  n_states            = "N States"
)

sumstat_df <- ts_data |>
  select(any_of(sumstat_vars)) |>
  pivot_longer(everything(), names_to = "variable", values_to = "value") |>
  group_by(variable) |>
  summarise(
    N      = sum(!is.na(value)),
    Mean   = mean(value, na.rm = TRUE),
    SD     = sd(value, na.rm = TRUE),
    Min    = min(value, na.rm = TRUE),
    Q25    = quantile(value, 0.25, na.rm = TRUE),
    Median = median(value, na.rm = TRUE),
    Q75    = quantile(value, 0.75, na.rm = TRUE),
    Max    = max(value, na.rm = TRUE),
    .groups = "drop"
  ) |>
  mutate(
    variable = factor(variable, levels = sumstat_vars),
    Label    = sumstat_labels[as.character(variable)]
  ) |>
  arrange(variable) |>
  select(Label, N, Mean, SD, Min, Q25, Median, Q75, Max)

# CSV export to data/processed for analyst use
write_csv(sumstat_df, file.path(OUTPUT_DIR, "summary_stats_system.csv"))
cat("  CSV -> summary_stats_system.csv\n")

# LaTeX export to results/main alongside the main-text tables
tex_lines <- c(
  "\\begin{table}[htbp]",
  "\\centering",
  "\\caption{Summary Statistics: System-Level Variables}",
  "\\label{tab:sumstats_system}",
  "\\small",
  "\\begin{tabular}{lrrrrrrrr}",
  "\\toprule",
  "Variable & N & Mean & SD & Min & Q25 & Median & Q75 & Max \\\\",
  "\\midrule"
)
for (i in seq_len(nrow(sumstat_df))) {
  row <- sumstat_df[i, ]
  tex_lines <- c(tex_lines, sprintf(
    "%s & %d & %.4f & %.4f & %.4f & %.4f & %.4f & %.4f & %.4f \\\\",
    row$Label, row$N, row$Mean, row$SD, row$Min, row$Q25,
    row$Median, row$Q75, row$Max
  ))
}
tex_lines <- c(tex_lines,
  "\\bottomrule",
  "\\end{tabular}",
  "\\par\\vspace{4pt}",
  paste0(
    "{\\footnotesize \\textit{Note:} Variables appearing in the main-text ",
    "Tables~\\ref{tab:system_main} and \\ref{tab:granger_fevd_mid} ",
    "full-control regressions. Both raw and $\\log(1+x)$ forms reported ",
    "for the DV (JS divergence) and IVs (MID count, fatal MID count) so ",
    "readers can see the distributional impact of the transform. ",
    "Exogenous controls (N Democracies, N Autocracies, N States) shown ",
    "on the raw scale; the regressions use $\\log(1+x)$ versions. ",
    "COW system, 1816--2014.}"
  ),
  "\\end{table}"
)
writeLines(tex_lines, file.path(TABLE_DIR, "summary_stats_system.tex"))
cat("  LaTeX -> summary_stats_system.tex (in results/main)\n")

cat("\nDone. System-level analysis complete.\n")