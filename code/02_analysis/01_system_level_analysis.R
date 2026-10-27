# [2026-10-01] Headline DV is RAW (unlogged) JS divergence; conflict counts stay logged.
# JS is bounded in [0, log 2] with mean ~0.016, so log(1 + JS) ~ JS and the log adds nothing.
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
cat(sprintf("R %s | dynamicmultiplex %s | igraph %s | fixest %s\n",
            getRversion(), packageVersion("dynamicmultiplex"),
            packageVersion("igraph"), packageVersion("fixest")))

# Inputs must come from the same (most recent) Stage 1 run: print their mtimes
# and stop if they were written more than 6 hours apart.
INPUT_FILES <- file.path(OUTPUT_DIR, c("system_analysis_data.rds", "dyad_year_network_layers.rds"))
stopifnot(all(file.exists(INPUT_FILES)))
input_mtime <- file.mtime(INPUT_FILES)
for (i in seq_along(INPUT_FILES)) cat(sprintf("  input %-32s mtime %s\n", basename(INPUT_FILES[i]), format(input_mtime[i], "%Y-%m-%d %H:%M:%S")))
if (diff(range(input_mtime)) > as.difftime(6, units = "hours"))
  stop("input files were written more than 6 hours apart; they are not from one Stage 1 run")

system_data <- readRDS(file.path(OUTPUT_DIR, "system_analysis_data.rds"))

# --------------------------------------------------------------------------- #
# 0a. Peace-dyad JS divergence [2026-10-12]
# --------------------------------------------------------------------------- #
# JS divergence between the year t-1 and year t edge-weight distributions after
# dropping, from BOTH years, every dyad that is in an ongoing MID in t-1 or t.
# Computed here from dyad_year_network_layers.rds so that Stage 1 is untouched.
# kl_divergence() and js_divergence() are copied verbatim from
# code/01_data_cleaning/02_compute_network_measures.R, and the no-drop version
# of the same code must reproduce the saved js_divergence exactly.
kl_divergence <- function(p, q, eps = 1e-10) {
  p <- p + eps
  q <- q + eps
  p <- p / sum(p)
  q <- q / sum(q)
  sum(p * log(p / q))
}
js_divergence <- function(p, q, eps = 1e-10) {
  p <- p + eps
  q <- q + eps
  p <- p / sum(p)
  q <- q / sum(q)
  m <- 0.5 * (p + q)
  0.5 * kl_divergence(p, m, eps = 0) + 0.5 * kl_divergence(q, m, eps = 0)
}

cat("Computing peace-dyad JS divergence from dyad_year_network_layers.rds...\n")
dyad_data <- readRDS(file.path(OUTPUT_DIR, "dyad_year_network_layers.rds"))
MAIN_WEIGHT_VAR <- "main_weight"
stopifnot(MAIN_WEIGHT_VAR %in% names(dyad_data), "cowmidongoing" %in% names(dyad_data))
js_years <- sort(unique(dyad_data$year))
stopifnot(identical(js_years, sort(unique(system_data$year))))

# year -> edges with positive weight (both stored directions, as in Stage 1)
js_edges <- lapply(js_years, function(yr) dyad_data |>
  filter(year == yr, .data[[MAIN_WEIGHT_VAR]] > 0) |>
  select(ccode1, ccode2, w = !!sym(MAIN_WEIGHT_VAR)))
names(js_edges) <- js_years
# year -> belligerent dyads (both directions) with an ongoing MID in that year
js_opp <- lapply(js_years, function(yr) {
  o <- dyad_data |> filter(year == yr, cowmidongoing == 1) |> select(ccode1, ccode2)
  bind_rows(o, o |> rename(ccode1 = ccode2, ccode2 = ccode1)) |> distinct()
})
names(js_opp) <- js_years
# [2026-10-12] year -> dyads (both directions) with an ongoing FATAL MID in that year (headline peace rule)
js_opp_fatal <- lapply(js_years, function(yr) {
  o <- dyad_data |> filter(year == yr, cowmidongoing == 1, cowmidfatal == 1) |> select(ccode1, ccode2)
  bind_rows(o, o |> rename(ccode1 = ccode2, ccode2 = ccode1)) |> distinct()
})
names(js_opp_fatal) <- js_years

js_pair <- function(e_t, e_t1) {
  merged <- full_join(e_t |> rename(w_t = w), e_t1 |> rename(w_t1 = w), by = c("ccode1", "ccode2")) |>
    mutate(w_t = replace_na(w_t, 0), w_t1 = replace_na(w_t1, 0))
  if (min(merged$w_t) < 0) merged$w_t <- merged$w_t - min(merged$w_t)
  if (min(merged$w_t1) < 0) merged$w_t1 <- merged$w_t1 - min(merged$w_t1)
  js_divergence(merged$w_t, merged$w_t1)
}

js_peace_rows <- vector("list", length(js_years))
for (i in seq_along(js_years)) {
  if (i == 1) { js_peace_rows[[i]] <- tibble(year = js_years[i], js_check = NA_real_, js_divergence_peace = NA_real_, bellig_weight_share = NA_real_,
                                             js_divergence_peace_ft1 = NA_real_, bellig_weight_share_ft1 = NA_real_); next }
  e_t  <- js_edges[[i]]; e_t1 <- js_edges[[i - 1]]
  drop <- bind_rows(js_opp[[i]], js_opp[[i - 1]]) |> distinct() |> mutate(drop = TRUE)
  e_t_p  <- e_t  |> anti_join(drop, by = c("ccode1", "ccode2"))
  e_t1_p <- e_t1 |> anti_join(drop, by = c("ccode1", "ccode2"))
  js_peace_rows[[i]] <- tibble(
    year = js_years[i],
    js_check = js_pair(e_t, e_t1),                      # no dyads dropped: must equal the saved js_divergence
    js_divergence_peace = js_pair(e_t_p, e_t1_p),
    bellig_weight_share = 1 - sum(e_t_p$w) / sum(e_t$w), # share of year-t tie weight on dyads in a MID in t-1 or t
    # [2026-10-12] headline peace rule: drop (from both years) dyads in an ongoing FATAL MID in t-1 only
    js_divergence_peace_ft1 = js_pair(e_t  |> anti_join(js_opp_fatal[[i - 1]], by = c("ccode1", "ccode2")),
                                      e_t1 |> anti_join(js_opp_fatal[[i - 1]], by = c("ccode1", "ccode2"))),
    bellig_weight_share_ft1 = 1 - sum((e_t |> anti_join(js_opp_fatal[[i - 1]], by = c("ccode1", "ccode2")))$w) / sum(e_t$w))
}
js_peace_df <- bind_rows(js_peace_rows)
stopifnot(nrow(js_peace_df) == nrow(system_data), !anyDuplicated(js_peace_df$year))
js_val <- system_data |> select(year, js_divergence) |> inner_join(js_peace_df, by = "year") |> filter(!is.na(js_divergence))
js_max_diff <- max(abs(js_val$js_check - js_val$js_divergence))
cat(sprintf("  Validation: no-drop JS vs saved js_divergence, max |diff| = %.3g over %d years\n", js_max_diff, nrow(js_val)))
if (!(js_max_diff < 1e-10)) stop("peace-dyad JS code does not reproduce the saved js_divergence (max |diff| = ", js_max_diff, ")")
stopifnot(all(is.na(js_val$js_divergence_peace) == is.na(js_val$js_divergence)),
          all(is.na(js_val$js_divergence_peace_ft1) == is.na(js_val$js_divergence)))
system_data <- system_data |> left_join(js_peace_df |> select(-js_check), by = "year")

# Footnote statistics (printed and saved; no prose drafted here)
js_fn <- system_data |> filter(!is.na(js_divergence), !is.na(js_divergence_peace))
js_footnote_stats <- tibble(
  statistic = c("cor_js_all_peace", "mean_js_all", "mean_js_peace",
                "bellig_weight_share_median", "bellig_weight_share_max", "bellig_weight_share_max_year", "n_years",
                "cor_js_all_peace_ft1", "mean_js_peace_ft1", "bellig_weight_share_ft1_median", "bellig_weight_share_ft1_max"),
  value = c(cor(js_fn$js_divergence, js_fn$js_divergence_peace), mean(js_fn$js_divergence), mean(js_fn$js_divergence_peace),
            median(js_fn$bellig_weight_share), max(js_fn$bellig_weight_share),
            js_fn$year[which.max(js_fn$bellig_weight_share)], nrow(js_fn),
            cor(js_fn$js_divergence, js_fn$js_divergence_peace_ft1), mean(js_fn$js_divergence_peace_ft1),
            median(js_fn$bellig_weight_share_ft1), max(js_fn$bellig_weight_share_ft1)))
cat("  JS peace-dyad footnote statistics:\n"); print(as.data.frame(js_footnote_stats), digits = 6)
dir.create(here::here("results", "main"), recursive = TRUE, showWarnings = FALSE)
write_csv(js_footnote_stats, here::here("results", "main", "js_peace_footnote_stats.csv"))
rm(dyad_data, js_edges, js_opp, js_opp_fatal); invisible(gc())

# --------------------------------------------------------------------------- #
# 1. Prepare time series
# --------------------------------------------------------------------------- #

# [2026-10-12] Conflict lags are built on the full annual series BEFORE dropping years without JS
# (1816), so 1817 keeps its 1816 conflict count. The AR(1) DV lag is still undefined in 1817.
stopifnot(all(diff(sort(system_data$year)) == 1), !anyDuplicated(system_data$year))
sys_conflict_lags <- system_data |>
  arrange(year) |>
  transmute(year,
            mid_count_lag1 = lag(mid_count),
            log_mid_count_lag1 = lag(log1p(mid_count)),
            fatal_mid_count_lag1 = lag(fatal_mid_count),
            log_fatal_mid_count_lag1 = lag(log1p(fatal_mid_count)))

ts_data <- system_data |>
  filter(!is.na(js_divergence)) |>
  select(
    year,
    kl_divergence,
    js_divergence,
    js_divergence_peace,
    bellig_weight_share,
    js_divergence_peace_ft1,
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
    # AR(1) lagged DV — used in the main-table dynamics column
    log_js_divergence_lag1 = lag(log_js_divergence),
    js_divergence_lag1 = lag(js_divergence),            # AR(1) term for the raw-JS main table
    js_divergence_peace_lag1 = lag(js_divergence_peace), # AR(1) term for the peace-dyad panels
    js_divergence_peace_ft1_lag1 = lag(js_divergence_peace_ft1), # AR(1) term for the headline peace panel
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
  ) |>
  left_join(sys_conflict_lags, by = "year")
stopifnot(nrow(ts_data) == sum(!is.na(system_data$js_divergence)))

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
  dv <- "js_divergence"   # raw JS in both sets; "logged" now means logged IVs only [2026-10-01]
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
# TRANSFORMATION [2026-10-01]: js_divergence now enters the VAR UNLOGGED (the text below
# describes the earlier log1p choice, made when the DV was KL). KL is non-negative and right-skewed, and right-skewness
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

build_var_inputs <- function(set, controls, iv, dv = "js_divergence") {
  stopifnot(dv %in% c("js_divergence", "js_divergence_peace_ft1"))
  endo_iv <- if (set == "logged") {
    if (iv == "mid") "log_mid_count" else "log_fatal_mid_count"
  } else {
    if (iv == "mid") "mid_count" else "fatal_mid_count"
  }
  endo_dv <- dv   # raw JS in both sets; "logged" now means logged conflict counts only [2026-10-01]
                  # [2026-10-14] dv = "js_divergence_peace_ft1" drops dyads in an ongoing fatal MID in t-1 (main-text Table 3, Panel B)

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
var_grid$dv <- "js_divergence"
# [2026-10-14] Main-text Table 3, Panel B: the same logged fatal-MID VARs with JS computed after dropping
# dyads in an ongoing fatal MID in t-1 (the peace rule used in Table 2, Panel B).
var_grid <- rbind(var_grid,
                  data.frame(set = "logged", controls = c("none", "full"), iv = "fatal",
                             dv = "js_divergence_peace_ft1", stringsAsFactors = FALSE))
stopifnot(nrow(var_grid) == 10, !anyDuplicated(var_grid))

for (i in seq_len(nrow(var_grid))) {
  spec <- var_grid[i, ]
  key  <- sprintf("var__%s__%s__%s%s", spec$set, spec$controls, spec$iv,
                  if (spec$dv == "js_divergence") "" else "__peacefatal")

  cat(sprintf("\n%s\n", strrep("=", 70)))
  cat(sprintf("  VAR [%s | %s | %s | %s]\n", spec$set, spec$controls, spec$iv, spec$dv))
  cat(sprintf("%s\n", strrep("=", 70)))

  inp <- build_var_inputs(spec$set, spec$controls, spec$iv, spec$dv)

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
  # [2026-10-12] boot = TRUE: without it vars::causality() ignores boot.runs and returns the asymptotic F-test,
  # while the table note reports 1,000 bootstrap replications. Seeded for reproducibility.
  set.seed(123)
  gc_to_kl <- causality(var_fit, cause = inp$endo_iv, boot = TRUE, boot.runs = 1000)
  set.seed(123)
  gc_kl_to <- causality(var_fit, cause = inp$endo_dv, boot = TRUE, boot.runs = 1000)
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
ts_data_main <- ts_data |> filter(!is.na(js_divergence_lag1), !is.na(js_divergence_peace_lag1), !is.na(js_divergence_peace_ft1_lag1),
                                    !is.na(mid_count_lag1))

# dv: "js_divergence" (all dyads) or "js_divergence_peace" (belligerent dyads dropped);
# the AR(1) column uses the lag of the same DV.
fit_main_table_4col <- function(iv_family, dv = "js_divergence") {
  stopifnot(iv_family %in% c("mid", "fatal"), dv %in% c("js_divergence", "js_divergence_peace", "js_divergence_peace_ft1"))
  iv_var <- if (iv_family == "mid") "log_mid_count_lag1"
                                     else "log_fatal_mid_count_lag1"
  dv_lag <- paste0(dv, "_lag1")

  # Col 1: IV + quadratic year trend (bivariate-with-trend baseline).
  # Including `poly(year, 2)` here keeps col 1 honest — with 1816--2014
  # annual data and clear secular drift in both KL and MID, a strict
  # bivariate would conflate the MID effect with the trend. All four
  # columns therefore share `poly(year, 2)` as a baseline; the ladder
  # varies by what's added on top.
  m1 <- feols(as.formula(paste(dv, "~", iv_var,
                                 "+ poly(year, 2)")),
              data = ts_data_main,
              panel.id = ~unit + year, vcov = NW(lag = NW_LAG))

  # Col 2: + full controls (adds n_dem, n_aut, n_states logged on top of
  # the year trend already in col 1).
  m2 <- feols(as.formula(paste(dv, "~", iv_var, "+",
                                 SYS_CTRL_FULL)),
              data = ts_data_main,
              panel.id = ~unit + year, vcov = NW(lag = NW_LAG))

  # Col 3: + AR(1) lagged DV
  m3 <- feols(as.formula(paste(dv, "~", iv_var, "+",
                                 SYS_CTRL_FULL,
                                 "+", dv_lag)),
              data = ts_data_main,
              panel.id = ~unit + year, vcov = NW(lag = NW_LAG))

  # Col 4: + order-period FE (passed via | order_period in feols formula)
  m4 <- feols(as.formula(paste(dv, "~", iv_var, "+",
                                 SYS_CTRL_FULL,
                                 "| order_period")),
              data = ts_data_main,
              panel.id = ~unit + year, vcov = NW(lag = NW_LAG))

  list(m1 = m1, m2 = m2, m3 = m3, m4 = m4)
}

main_4col_mid   <- fit_main_table_4col("mid")
main_4col_fatal <- fit_main_table_4col("fatal")
main_4col_mid_peace   <- fit_main_table_4col("mid",   dv = "js_divergence_peace")
main_4col_fatal_peace <- fit_main_table_4col("fatal", dv = "js_divergence_peace")
sys_panels <- list(A = main_4col_mid, B = main_4col_mid_peace, C = main_4col_fatal, D = main_4col_fatal_peace)
for (pn in names(sys_panels)) stopifnot(all(sapply(sys_panels[[pn]], nobs) == 197))
for (pn in names(sys_panels)) {
  iv <- if (pn %in% c("A", "B")) "log_mid_count_lag1" else "log_fatal_mid_count_lag1"
  cat(sprintf("  Panel %s [%s]: %s\n", pn, iv, paste(sapply(sys_panels[[pn]], function(m) {
    ct <- coeftable(m); sprintf("%.4f (%.4f) p=%.3g", ct[iv, 1], ct[iv, 2], ct[iv, 4]) }), collapse = " | ")))
}

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
  "js_divergence_lag1"           = "$\\text{JS}_{t{-}1}$ (AR(1))",
  "(Intercept)"                  = "Intercept"
)

# Tables report coefficients and standard errors multiplied by 100 (two decimals, so no cell rounds to 0.00).
# The AR(1) lagged-DV coefficient is unchanged by rescaling both the DV and its lag, so it is not multiplied.
TAB_SCALE <- 100
tab_coef <- function(ms) lapply(ms, function(m) { b <- coef(m); b * ifelse(grepl("divergence_lag1$", names(b)), 1, TAB_SCALE) })
tab_se   <- function(ms) lapply(ms, function(m) { s <- se(m);   s * ifelse(grepl("divergence_lag1$", names(s)), 1, TAB_SCALE) })
TAB_SCALE_NOTE <- "Coefficients and standard errors are multiplied by 100 (the AR(1) coefficient is unaffected by this rescaling). "

write_main_table_4col <- function(models_list, iv_family,
                                    file_path, table_dir_label) {
  stopifnot(iv_family %in% c("mid", "fatal"))
  iv_short <- if (iv_family == "mid")
    "$\\log(1{+}\\text{MID})$"
  else
    "$\\log(1{+}\\text{Fatal MID})$"

  caption_lead <- if (iv_family == "mid")
    paste0("System-Level OLS: JS divergence on ", iv_short,
           " (HAC inference throughout). All four columns include a ",
           "quadratic year trend. Each subsequent column tightens one ",
           "identification challenge: column 2 adds observable controls; ",
           "column 3 additionally controls for JS momentum (AR(1) lagged ",
           "DV); column 4 absorbs four international-order regimes via ",
           "fixed effects.")
  else
    paste0("Appendix. System-Level OLS: JS divergence on ",
           iv_short, " (companion to main-text Table~\\ref{tab:system_main}).")

  label_suffix <- if (iv_family == "mid") "system_main" else "system_appx_fatal"

  ms <- list(models_list$m1, models_list$m2, models_list$m3, models_list$m4)
  texreg(
    ms, override.coef = tab_coef(ms), override.se = tab_se(ms),
    custom.model.names = c("(1) IV + trend",
                            "(2) + Full ctrls",
                            "(3) + AR(1)",
                            "(4) + Period FE"),
    custom.coef.map = ols_coef_map_4col,
    caption = caption_lead,
    label   = sprintf("tab:%s", label_suffix),
    custom.note = paste0(
      "\\parbox{\\linewidth}{\\vspace{2pt}\\textit{Note:} ",
      "DV is JS divergence (not logged). IV is ", iv_short, " lagged one ",
      "period ($t-1$). Newey-West HAC standard errors (lag $=$ 5) in ",
      "parentheses. ", TAB_SCALE_NOTE, "\\textbf{All four columns include a quadratic year ",
      "trend (\\texttt{poly(year, 2)}); terms suppressed for legibility.} ",
      "Column 1 is the IV plus the year trend (no other controls). ",
      "Columns 2--4 add the full controls covariate set ",
      "($\\log(1+\\text{Democracies})$, $\\log(1+\\text{Autocracies})$, ",
      "$\\log(1+\\text{States})$). Column 3 additionally controls for ",
      "$\\text{JS}_{t-1}$ (AR(1) term). Column 4 absorbs ",
      "fixed effects for four international-order regimes (Concert of ",
      "Europe / WWI--Interwar--WWII / Cold War / Post-Cold War). ",
      "COW system, 1816--2014. ",
      "$^{***}p<0.001$; $^{**}p<0.01$; $^{*}p<0.05$; ",
      "$^{\\dagger}p<0.1$.}"
    ),
    stars = c(0.001, 0.01, 0.05, 0.1),
    symbol = "\\dagger",
    digits = 2,  # coefficients are x 100 (TAB_SCALE), so two decimals show every estimate
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

# [2026-10-12] Main text: four stacked panels (MID / fatal MID x all dyads / peace dyads).
# tableA_system_fatal.tex is retired; the fatal results are Panels C-D here.
# write_main_table_4col() is kept for reference and is no longer called.
sys_rows <- function(iv_var, dv_lag) list(
  c(ols_coef_map_4col[[iv_var]], iv_var),
  c("$\\log(1{+}\\text{Democracies})$", "log1p(n_democracies)"),
  c("$\\log(1{+}\\text{Autocracies})$", "log1p(n_autocracies)"),
  c("$\\log(1{+}\\text{States})$",      "log1p(n_states)"),
  c("$\\text{JS}_{t{-}1}$ (AR(1))",        dv_lag),
  c("Intercept",                           "(Intercept)"))
sys_star <- function(p) ifelse(p < 0.001, "^{***}", ifelse(p < 0.01, "^{**}", ifelse(p < 0.05, "^{*}", ifelse(p < 0.1, "^{\\dagger}", ""))))
sys_panel <- function(models, iv_var, dv_lag) {
  out <- character(0)
  for (r in sys_rows(iv_var, dv_lag)) {
    k <- if (r[2] == dv_lag) 1 else TAB_SCALE     # the AR(1) coefficient is not rescaled
    est <- se <- character(length(models))
    for (m in seq_along(models)) {
      ct <- coeftable(models[[m]])
      if (r[2] %in% rownames(ct)) {
        est[m] <- sprintf("$%.2f%s$", k * ct[r[2], 1], sys_star(ct[r[2], 4]))
        se[m]  <- sprintf("$(%.2f)$", k * ct[r[2], 2])
      } else { est[m] <- "--"; se[m] <- "" }
    }
    out <- c(out, paste(r[1], "&", paste(est, collapse = " & "), "\\\\"), paste(" &", paste(se, collapse = " & "), "\\\\"))
  }
  w <- function(m, type) { x <- tryCatch(fixest::r2(m, type), error = function(e) NA); if (is.na(x)) "--" else sprintf("%.2f", x) }
  g <- function(m) { n <- m$fixef_sizes; if (is.null(n) || !"order_period" %in% names(n)) "--" else as.character(n[["order_period"]]) }
  c(out, "\\midrule",
    paste("Num. obs. &", paste(sapply(models, nobs), collapse = " & "), "\\\\"),
    paste("R$^2$ (full model) &", paste(sapply(models, w, "r2"), collapse = " & "), "\\\\"),
    paste("R$^2$ (proj model) &", paste(sapply(models, w, "wr2"), collapse = " & "), "\\\\"),
    paste("Num. groups: Order period &", paste(sapply(models, g), collapse = " & "), "\\\\"))
}
sys_panel_titles <- c(
  A = "Panel A. $\\log(1+\\text{MID count})_{t-1}$; JS divergence over all dyads",
  B = "Panel B. $\\log(1+\\text{MID count})_{t-1}$; JS divergence over peace dyads only",
  C = "Panel C. $\\log(1+\\text{fatal MID count})_{t-1}$; JS divergence over all dyads",
  D = "Panel D. $\\log(1+\\text{fatal MID count})_{t-1}$; JS divergence over peace dyads only")
sys_ivs  <- c(A = "log_mid_count_lag1", B = "log_mid_count_lag1", C = "log_fatal_mid_count_lag1", D = "log_fatal_mid_count_lag1")
sys_lags <- c(A = "js_divergence_lag1", B = "js_divergence_peace_lag1", C = "js_divergence_lag1", D = "js_divergence_peace_lag1")
sys_body <- character(0)
for (pn in names(sys_panels)) {
  sys_body <- c(sys_body, if (pn != "A") "\\midrule",
                sprintf("\\multicolumn{5}{l}{\\textit{%s}} \\\\", sys_panel_titles[[pn]]), "\\midrule",
                sys_panel(unname(sys_panels[[pn]]), sys_ivs[[pn]], sys_lags[[pn]]))
}
sys_note <- paste0(
  "\\multicolumn{5}{l}{\\scriptsize{\\parbox{\\linewidth}{\\vspace{2pt}\\textit{Note:} ",
  "DV is JS divergence between the year $t-1$ and year $t$ tie-weight distributions (not logged). ",
  "Panels B and D drop from both years every dyad that was in an ongoing MID in $t-1$ or $t$ before computing JS; ",
  "their AR(1) term is the lag of the peace-dyad JS. ",
  "The conflict measure is $\\log(1+\\text{MID count})$ in Panels A--B and $\\log(1+\\text{fatal MID count})$ in Panels C--D, ",
  "lagged one period ($t-1$). Newey-West HAC standard errors (lag $=$ 5) in parentheses. ", TAB_SCALE_NOTE,
  "All four columns include a quadratic year trend (\\texttt{poly(year, 2)}); terms suppressed for legibility. ",
  "Column 1 is the conflict measure plus the year trend. Columns 2--4 add $\\log(1+\\text{Democracies})$, ",
  "$\\log(1+\\text{Autocracies})$ and $\\log(1+\\text{States})$. Column 3 additionally controls for the lagged DV (AR(1)). ",
  "Column 4 absorbs fixed effects for four international-order regimes (Concert of Europe / WWI--Interwar--WWII / Cold War / Post-Cold War). ",
  "COW system, 1816--2014. $^{***}p<0.001$; $^{**}p<0.01$; $^{*}p<0.05$; $^{\\dagger}p<0.1$.}}} \\\\")
# [2026-10-12] Split a four-panel table into two floats (Panels A-B; Panels C-D) so neither exceeds a
# page. Part 2 repeats the column header and note; its caption is prefixed "(Continued)" and its label
# gets the suffix "_cd". \ContinuedFloat (caption package) keeps one table number for both parts.
# Fails loudly if the table does not have the expected layout.
split_ab_cd <- function(tex) {
  i_top <- which(tex == "\\toprule"); i_bot <- which(tex == "\\bottomrule")
  i_C   <- grep("textit\\{Panel C\\.", tex)
  i_cap <- grep("^\\\\caption\\{", tex); i_lab <- grep("^\\\\label\\{", tex)
  stopifnot(length(i_top) == 1, length(i_bot) == 1, length(i_C) == 1, length(i_cap) == 1, length(i_lab) == 1,
            tex[i_top + 2] == "\\midrule", tex[i_C - 1] == "\\midrule", i_top + 2 < i_C - 2, i_C < i_bot)
  pre   <- tex[1:(i_top + 2)]
  post  <- tex[i_bot:length(tex)]
  post2 <- post
  post2[post2 == tex[i_cap]] <- paste0("\\ContinuedFloat ", sub("^\\\\caption\\{", "\\\\caption{(Continued) ", tex[i_cap]))
  post2[post2 == tex[i_lab]] <- sub("\\}$", "_cd}", tex[i_lab])
  part1 <- c(pre, tex[(i_top + 3):(i_C - 2)], post)
  part2 <- c(pre, tex[i_C:(i_bot - 1)], post2)
  stopifnot(!any(grepl("Panel C\\.|Panel D\\.", part1)), !any(grepl("Panel A\\.|Panel B\\.", part2)),
            sum(grepl("Panel [AB]\\.", part1)) == 2, sum(grepl("Panel [CD]\\.", part2)) == 2,
            any(grepl("_cd\\}$", part2)), any(grepl("\\(Continued\\)", part2)), any(grepl("ContinuedFloat", part2)))
  c(part1, "", part2)
}
sys_tex <- c(
  "\\begin{table}[p]", "\\scriptsize", "\\begin{center}", "\\begin{tabular}{l c c c c}", "\\toprule",
  " & (1) IV + trend & (2) + Full ctrls & (3) + AR(1) & (4) + Period FE \\\\", "\\midrule",
  sys_body, "\\bottomrule", sys_note,
  "\\end{tabular}", "\\end{center}",
  "\\caption{System-level regressions for conflict and network reconfiguration. MID count (Panels A--B) and fatal MID count (Panels C--D), each with JS divergence over all dyads and over peace dyads only.}",
  "\\label{tab:system_allmid}", "\\end{table}")
sys_tex[grepl("^\\\\caption\\{System-level regressions", sys_tex)] <- paste0(
  "\\caption{Appendix. System-level regressions for all MIDs and fatal MIDs, each with JS divergence over all dyads and over peace dyads ",
  "(dyads in any MID in $t-1$ or $t$ dropped), companion to main-text Table~\\ref{tab:system_main}.}")
writeLines(split_ab_cd(sys_tex), file.path(here::here("results", "appendix"), "tableA_system_allmid.tex"))
cat("  tab:system_allmid -> results/appendix/tableA_system_allmid.tex (Panels A-D, split A-B / C-D)\n")

# =========================================================================== #
# [2026-10-12] HEADLINE, main-text system table: log(1 + fatal MID count) at t-1.
#   Panel A: JS over all dyads. Panel B: JS over peace dyads, dropping from both years every dyad in an
#   ongoing FATAL MID in t-1 (js_divergence_peace_ft1).
# =========================================================================== #
hl_sys_panels <- list(A = main_4col_fatal, B = fit_main_table_4col("fatal", dv = "js_divergence_peace_ft1"))
for (pn in names(hl_sys_panels)) stopifnot(all(sapply(hl_sys_panels[[pn]], nobs) == 197))
hl_sys_lags <- c(A = "js_divergence_lag1", B = "js_divergence_peace_ft1_lag1")
# [2026-10-27] Main-text Table 1 as one 8-column table: (1)-(4) Panel A (all dyads), (5)-(8) Panel B (peace dyads).
#   Variable labels carry a superscript a; the log(1+x) transform is stated once in the note.
hl8_models <- c(unname(hl_sys_panels$A), unname(hl_sys_panels$B))
hl8_lags   <- c(rep(hl_sys_lags[["A"]], 4), rep(hl_sys_lags[["B"]], 4))
stopifnot(length(hl8_models) == 8, all(sapply(hl8_models, nobs) == 197))
hl8_rows <- list(
  c("Fatal MID Count$_{t-1}^{a}$", "log_fatal_mid_count_lag1"),
  c("Democracies$^{a}$",           "log1p(n_democracies)"),
  c("Autocracies$^{a}$",           "log1p(n_autocracies)"),
  c("States$^{a}$",                "log1p(n_states)"),
  c("JS Divergence$_{t-1}$",       "DV_LAG"),
  c("Intercept",                   "(Intercept)"))
hl_sys_body <- character(0)
for (r in hl8_rows) {
  est <- se <- character(8)
  for (m in 1:8) {
    term <- if (r[2] == "DV_LAG") hl8_lags[m] else r[2]
    k <- if (r[2] == "DV_LAG") 1 else TAB_SCALE        # the AR(1) coefficient is not rescaled
    ct <- coeftable(hl8_models[[m]])
    if (term %in% rownames(ct)) {
      est[m] <- sprintf("$%.2f%s$", k * ct[term, 1], sys_star(ct[term, 4]))
      se[m]  <- sprintf("$(%.2f)$", k * ct[term, 2])
    } else { est[m] <- "--"; se[m] <- "" }
  }
  hl_sys_body <- c(hl_sys_body, paste(r[1], "&", paste(est, collapse = " & "), "\\\\"),
                   paste(" &", paste(se, collapse = " & "), "\\\\"))
}
hl8_r2  <- function(m, type) { x <- tryCatch(fixest::r2(m, type), error = function(e) NA); if (is.na(x)) "--" else sprintf("%.2f", x) }
hl8_grp <- function(m) { n <- m$fixef_sizes; if (is.null(n) || !"order_period" %in% names(n)) "--" else as.character(n[["order_period"]]) }
hl_sys_body <- c(hl_sys_body, "\\midrule",
  paste("Num. obs. &", paste(sapply(hl8_models, nobs), collapse = " & "), "\\\\"),
  paste("R$^2$ (full model) &", paste(sapply(hl8_models, hl8_r2, "r2"), collapse = " & "), "\\\\"),
  paste("R$^2$ (proj model) &", paste(sapply(hl8_models, hl8_r2, "wr2"), collapse = " & "), "\\\\"),
  paste("Num. groups: Order period &", paste(sapply(hl8_models, hl8_grp), collapse = " & "), "\\\\"))
hl_sys_note <- paste0(
  "\\multicolumn{9}{l}{\\scriptsize{\\parbox{\\linewidth}{\\vspace{2pt}\\textit{Note:} ",
  "DV is JS divergence between the year $t-1$ and year $t$ tie-weight distributions (not logged). ",
  "Panel B drops from both years every dyad that was in an ongoing fatal MID in $t-1$ before computing JS; ",
  "its lagged DV is the lag of the same peace-dyad JS. ",
  "$^{a}$Entered as $\\log(1+x)$. The conflict measure is lagged one period ($t-1$). ",
  "Newey-West HAC standard errors (lag $=$ 5) in parentheses. ", TAB_SCALE_NOTE,
  "All eight columns include a quadratic year trend (\\texttt{poly(year, 2)}); terms suppressed for legibility. ",
  "Columns 1 and 5 include the conflict measure and the year trend. Columns 2--4 and 6--8 add the counts of democracies, ",
  "autocracies and states. Columns 3 and 7 additionally control for the lagged DV. ",
  "Columns 4 and 8 absorb fixed effects for four international-order regimes (Concert of Europe / WWI--Interwar--WWII / Cold War / Post-Cold War). ",
  "All MIDs and the any-MID peace rule are in Appendix Table~\\ref{tab:system_allmid}. ",
  "COW system, 1816--2014. $^{***}p<0.001$; $^{**}p<0.01$; $^{*}p<0.05$; $^{\\dagger}p<0.1$.}}} \\\\")
hl_sys_tex <- c(
  "\\begin{table}[p]", "\\scriptsize", "\\begin{center}", "\\setlength{\\tabcolsep}{4pt}",
  "\\begin{tabular*}{\\linewidth}{@{\\extracolsep{\\fill}}l c c c c c c c c}", "\\toprule",
  " & \\multicolumn{4}{c}{Panel A. JS Divergence, all dyads} & \\multicolumn{4}{c}{Panel B. JS Divergence, peace dyads} \\\\",
  "\\cmidrule(lr){2-5}\\cmidrule(lr){6-9}",
  " & (1) & (2) & (3) & (4) & (5) & (6) & (7) & (8) \\\\", "\\midrule",
  hl_sys_body, "\\bottomrule", "\\end{tabular*}", "\\end{center}", "\\vspace{-8pt}",
  sub("}}} \\\\", "\\par}", sub("\\multicolumn{9}{l}{\\scriptsize{\\parbox{\\linewidth}{\\vspace{2pt}", "{\\scriptsize ", hl_sys_note, fixed = TRUE), fixed = TRUE),
  "\\caption{System-level regressions for fatal conflict and network reconfiguration. Fatal MID count in $t-1$, with JS divergence over all dyads (Panel A) and over peace dyads (Panel B).}",
  "\\label{tab:system_main}", "\\end{table}")
stopifnot(sum(grepl("Fatal MID Count", hl_sys_tex) & grepl("&", hl_sys_tex)) == 1, !any(grepl("log\\(1\\{\\+\\}", hl_sys_tex)),
          all(lengths(regmatches(hl_sys_body[1:12], gregexpr("&", hl_sys_body[1:12]))) == 8))
writeLines(hl_sys_tex, file.path(TABLE_DIR, "table1_system_main.tex"))
cat("  tab:system_main -> table1_system_main.tex (fatal MID count t-1; Panels A-B)\n")
write_csv(bind_rows(lapply(names(hl_sys_panels), function(pn) bind_rows(lapply(seq_along(hl_sys_panels[[pn]]), function(m)
  broom::tidy(hl_sys_panels[[pn]][[m]]) |> mutate(panel = pn, column = m, iv = "log_fatal_mid_count_lag1",
                                                  dv = if (pn == "A") "js_divergence" else "js_divergence_peace_ft1",
                                                  nobs = nobs(hl_sys_panels[[pn]][[m]])))))),
  file.path(OUTPUT_DIR, "table1_system_main_coefs.csv"))
write_csv(bind_rows(lapply(names(sys_panels), function(pn) bind_rows(lapply(seq_along(sys_panels[[pn]]), function(m)
  broom::tidy(sys_panels[[pn]][[m]]) |> mutate(panel = pn, column = m, iv = sys_ivs[[pn]],
                                               dv = if (pn %in% c("A", "C")) "js_divergence" else "js_divergence_peace",
                                               nobs = nobs(sys_panels[[pn]][[m]])))))),
  file.path(OUTPUT_DIR, "tableA_system_allmid_coefs.csv"))

# Appendix companion [2026-10-02]: the main table with log(1 + JS) as the DV.
# The main text and every other appendix model use raw JS; this is the one logged-DV system table.
stopifnot(all(c("log_js_divergence", "log_js_divergence_lag1") %in% names(ts_data_main)))
logjs_4col <- list(
  m1 = feols(log_js_divergence ~ log_fatal_mid_count_lag1 + poly(year, 2),
             data = ts_data_main, panel.id = ~unit + year, vcov = NW(lag = NW_LAG)),
  m2 = feols(as.formula(paste("log_js_divergence ~ log_fatal_mid_count_lag1 +", SYS_CTRL_FULL)),
             data = ts_data_main, panel.id = ~unit + year, vcov = NW(lag = NW_LAG)),
  m3 = feols(as.formula(paste("log_js_divergence ~ log_fatal_mid_count_lag1 +", SYS_CTRL_FULL, "+ log_js_divergence_lag1")),
             data = ts_data_main, panel.id = ~unit + year, vcov = NW(lag = NW_LAG)),
  m4 = feols(as.formula(paste("log_js_divergence ~ log_fatal_mid_count_lag1 +", SYS_CTRL_FULL, "| order_period")),
             data = ts_data_main, panel.id = ~unit + year, vcov = NW(lag = NW_LAG)))
stopifnot(all(sapply(logjs_4col, nobs) == nobs(main_4col_mid$m1)))
texreg(
  logjs_4col, override.coef = tab_coef(logjs_4col), override.se = tab_se(logjs_4col),
  custom.model.names = c("(1) IV + trend", "(2) + Full ctrls", "(3) + AR(1)", "(4) + Period FE"),
  custom.coef.map = c(ols_coef_map_4col[c("log_fatal_mid_count_lag1", "log1p(n_democracies)", "log1p(n_autocracies)", "log1p(n_states)")],
                      list("log_js_divergence_lag1" = "$\\log(1{+}\\text{JS})_{t{-}1}$ (AR(1))", "(Intercept)" = "Intercept")),
  caption = paste0("Appendix. System-Level OLS with $\\log(1+\\text{JS})$ as the dependent variable ",
                   "(companion to main-text Table~\\ref{tab:system_main})."),
  label = "tab:system_main_logjs",
  custom.note = paste0(
    "\\parbox{\\linewidth}{\\vspace{2pt}\\textit{Note:} DV is $\\log(1+\\text{JS divergence})$. ",
    "IV is $\\log(1{+}\\text{fatal MID count})$ lagged one period ($t-1$). Newey-West HAC standard errors (lag $=$ 5) in parentheses. ", TAB_SCALE_NOTE,
    "All four columns include a quadratic year trend. Columns and sample match Table~\\ref{tab:system_main}. ",
    "COW system, 1816--2014. $^{***}p<0.001$; $^{**}p<0.01$; $^{*}p<0.05$; $^{\\dagger}p<0.1$.}"),
  stars = c(0.001, 0.01, 0.05, 0.1), symbol = "\\dagger", digits = 2, float.pos = "htbp", booktabs = TRUE,
  use.packages = FALSE, include.rsquared = TRUE, include.adjrs = FALSE, include.nobs = TRUE,
  file = file.path(here::here("results", "appendix"), "tableA_system_main_logjs.tex"))
cat("  System main 4-col [log JS companion] -> tableA_system_main_logjs.tex (appendix)\n")

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
# [2026-10-14] Rows for one outcome (all dyads, or belligerent dyads dropped): forward and reverse
# Granger (bootstrap p) and FEVD at h = 1, 3, 5, 10, for the no-controls and full-controls VARs.
gf_rows <- function(iv_family, suffix) {
  out <- list(fwd_F = character(2), fwd_p = character(2), rev_F = character(2), rev_p = character(2),
              h1 = character(2), h3 = character(2), h5 = character(2), h10 = character(2))
  ctrls <- c("none", "full")
  for (j in 1:2) {
    key <- sprintf("var__logged__%s__%s%s", ctrls[j], iv_family, suffix)
    g   <- granger_results_list[[key]]
    fd  <- fevd_results_list[[key]]
    stopifnot(!is.null(g), !is.null(fd), nrow(fd) >= 10)
    fF <- as.numeric(g$fwd$Granger$statistic); fP <- as.numeric(g$fwd$Granger$p.value)
    rF <- as.numeric(g$rev$Granger$statistic); rP <- as.numeric(g$rev$Granger$p.value)
    out$fwd_F[j] <- sprintf("%.3f", fF)
    out$fwd_p[j] <- if (!is.na(fP) && fP < 0.001) paste0("$<$0.001", star_fn(fP)) else sprintf("%.4f%s", fP, star_fn(fP))   # bootstrap p
    out$rev_F[j] <- sprintf("%.3f", rF)
    out$rev_p[j] <- if (!is.na(rP) && rP < 0.001) paste0("$<$0.001", star_fn(rP)) else sprintf("%.4f%s", rP, star_fn(rP))
    out$h1[j]  <- sprintf("%.1f\\%%", fd[[g$endo_iv]][1]  * 100)
    out$h3[j]  <- sprintf("%.1f\\%%", fd[[g$endo_iv]][3]  * 100)
    out$h5[j]  <- sprintf("%.1f\\%%", fd[[g$endo_iv]][5]  * 100)
    out$h10[j] <- sprintf("%.1f\\%%", fd[[g$endo_iv]][10] * 100)
  }
  out
}
gf_block <- function(r) c(
  "\\multicolumn{3}{l}{\\quad Granger causality, conflict $\\rightarrow$ JS} \\\\",
  paste0("\\quad $F$-statistic & ", paste(r$fwd_F, collapse = " & "), " \\\\"),
  paste0("\\quad $p$-value & ",     paste(r$fwd_p, collapse = " & "), " \\\\[3pt]"),
  "\\multicolumn{3}{l}{\\quad Granger causality, JS $\\rightarrow$ conflict} \\\\",
  paste0("\\quad $F$-statistic & ", paste(r$rev_F, collapse = " & "), " \\\\"),
  paste0("\\quad $p$-value & ",     paste(r$rev_p, collapse = " & "), " \\\\[3pt]"),
  "\\multicolumn{3}{l}{\\quad FEVD, \\% of JS variance from conflict shocks} \\\\",
  paste0("\\quad $h = 1$ & ",  paste(r$h1,  collapse = " & "), " \\\\"),
  paste0("\\quad $h = 3$ & ",  paste(r$h3,  collapse = " & "), " \\\\"),
  paste0("\\quad $h = 5$ & ",  paste(r$h5,  collapse = " & "), " \\\\"),
  paste0("\\quad $h = 10$ & ", paste(r$h10, collapse = " & "), " \\\\"))

write_granger_fevd_by_iv <- function(iv_family, file_path, table_dir_label) {
  stopifnot(iv_family %in% c("mid", "fatal"))
  # [2026-10-27] Main-text fatal table: one 4-column table, (1)-(2) Panel A (all dyads), (3)-(4) Panel B (peace dyads),
  #   no-controls and full-controls VARs in each panel. The log(1+x) transform is in the note (superscript a).
  if (iv_family == "fatal") {
    rA <- gf_rows("fatal", ""); rB <- gf_rows("fatal", "__peacefatal")
    rr <- function(lab, key, end = " \\\\") paste0(lab, " & ", paste(c(rA[[key]], rB[[key]]), collapse = " & "), end)
    body4 <- c(
      "\\multicolumn{5}{l}{\\textit{Granger causality, fatal MID$^{a}$ $\\rightarrow$ JS divergence}} \\\\",
      rr("\\quad $F$-statistic", "fwd_F"), rr("\\quad $p$-value", "fwd_p", " \\\\[3pt]"),
      "\\multicolumn{5}{l}{\\textit{Granger causality, JS divergence $\\rightarrow$ fatal MID$^{a}$}} \\\\",
      rr("\\quad $F$-statistic", "rev_F"), rr("\\quad $p$-value", "rev_p", " \\\\[3pt]"),
      "\\multicolumn{5}{l}{\\textit{FEVD, \\% of JS divergence variance from fatal MID$^{a}$ shocks}} \\\\",
      rr("\\quad $h = 1$", "h1"), rr("\\quad $h = 3$", "h3"), rr("\\quad $h = 5$", "h5"), rr("\\quad $h = 10$", "h10"),
      "\\midrule",
      "Controls$^{a}$ & No & Yes & No & Yes \\\\")
    stopifnot(sum(grepl("F\\$-statistic", body4)) == 2,
              all(lengths(regmatches(body4[grepl("^\\\\quad", body4)], gregexpr("&", body4[grepl("^\\\\quad", body4)]))) == 4))
    tex4 <- c(
      "\\begin{table}[htbp]", "\\centering",
      "\\caption{Granger causality and FEVD: fatal MID count and JS divergence, over all dyads (Panel A) and over peace dyads (Panel B).}",
      "\\label{tab:granger_fevd_fatal}",
      "\\begin{tabular}{lcccc}", "\\toprule",
      " & \\multicolumn{2}{c}{Panel A. JS Divergence, all dyads} & \\multicolumn{2}{c}{Panel B. JS Divergence, peace dyads} \\\\",
      "\\cmidrule(lr){2-3}\\cmidrule(lr){4-5}",
      " & (1) & (2) & (3) & (4) \\\\", "\\midrule",
      body4, "\\bottomrule", "\\end{tabular}", "\\par\\vspace{4pt}",
      paste0(
        "{\\footnotesize \\textit{Note:} Granger tests from VARs with AIC-selected lag order (chosen separately for each VAR) ",
        "and 1,000 bootstrap replications. Endogenous variables are JS divergence (not logged) and the fatal MID count. ",
        "$^{a}$Entered as $\\log(1+x)$; this applies to the fatal MID count and to the exogenous controls. ",
        "Panel B computes JS after dropping, from both years, every dyad that was in an ongoing fatal MID in $t-1$, ",
        "as in Table~\\ref{tab:system_main}, Panel B. Columns 1 and 3 include only a constant and trend. ",
        "Columns 2 and 4 additionally include the exogenous counts of democracies, autocracies and states. ",
        "The FEVD rows report the share of forecast variance in JS divergence attributable to fatal MID shocks ",
        "at horizons $h = 1, 3, 5, 10$ years. Raw-scale variants are reported in the underlying ",
        "\\texttt{granger\\_results\\_combined.csv} and \\texttt{fevd\\_results\\_combined.csv}. COW system, 1816--2014. ",
        "$^{***}p<0.001$; $^{**}p<0.01$; $^{*}p<0.05$; $^{\\dagger}p<0.1$.}"),
      "\\end{table}")
    stopifnot(!any(grepl("log(1{+}", tex4, fixed = TRUE)))
    writeLines(tex4, file_path)
    cat(sprintf("  Granger + FEVD [fatal, 4-col] -> %s (%s)\n", basename(file_path), table_dir_label))
    return(invisible(NULL))
  }
  iv_short <- if (iv_family == "mid") "$\\log(1{+}\\text{MID})$" else "$\\log(1{+}\\text{Fatal MID})$"
  col_headers <- c(paste0(iv_short, ", no ctrls"), paste0(iv_short, ", full ctrls"))

  if (iv_family == "fatal") {
    # main text: Panel A all dyads, Panel B dyads in an ongoing fatal MID in t-1 dropped
    body <- c("\\multicolumn{3}{l}{\\textit{Panel A. JS divergence over all dyads}} \\\\[3pt]",
              gf_block(gf_rows("fatal", "")),
              "\\midrule",
              "\\multicolumn{3}{l}{\\textit{Panel B. JS divergence over peace dyads (dyads in an ongoing fatal MID in $t-1$ dropped)}} \\\\[3pt]",
              gf_block(gf_rows("fatal", "__peacefatal")))
    caption_lead <- paste0("Granger causality and FEVD: $\\log(1+\\text{Fatal MID Count})$ and JS divergence, ",
                           "over all dyads (Panel A) and over peace dyads (Panel B).")
    panel_note <- paste0("Panel B computes JS after dropping, from both years, every dyad that was in an ongoing fatal MID in $t-1$, ",
                         "as in Table~\\ref{tab:system_main}, Panel B. ")
  } else {
    body <- c("\\multicolumn{3}{l}{\\textit{JS divergence over all dyads}} \\\\[3pt]",
              gf_block(gf_rows("mid", "")))
    caption_lead <- paste0("Appendix. Granger causality and FEVD: $\\log(1+\\text{MID Count})$ and JS divergence. ",
                           "Companion to main-text Table~\\ref{tab:granger_fevd_fatal}.")
    panel_note <- ""
  }
  stopifnot(sum(grepl("F\\$-statistic", body)) == if (iv_family == "fatal") 4 else 2)

  tex_lines <- c(
    "\\begin{table}[htbp]",
    "\\centering",
    sprintf("\\caption{%s}", caption_lead),
    sprintf("\\label{tab:granger_fevd_%s}", iv_family),
    "\\begin{tabular}{lcc}",
    "\\toprule",
    paste0(" & ", paste(col_headers, collapse = " & "), " \\\\"),
    "\\midrule",
    body,
    "\\bottomrule",
    "\\end{tabular}",
    "\\par\\vspace{4pt}",
    paste0(
      "{\\footnotesize \\textit{Note:} Granger tests from VAR with ",
      "AIC-selected lag order (chosen separately for each VAR) and 1,000 bootstrap replications. ",
      "Endogenous variables are JS divergence (not logged) and ",
      iv_short, ", to match the ",
      if (iv_family == "fatal") "main-text Table~\\ref{tab:system_main}" else "Appendix Table~\\ref{tab:system_allmid}",
      " anchoring. ", panel_note,
      "``no ctrls'' VARs include only ",
      "constant and trend; ``full ctrls'' VARs additionally include ",
      "exogenous $\\log(1+\\text{Democracies})$, ",
      "$\\log(1+\\text{Autocracies})$, $\\log(1+\\text{States})$. The FEVD rows ",
      "report the share of forecast variance in JS divergence ",
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
write_granger_fevd_by_iv("fatal",
                          file.path(TABLE_DIR, "table2_granger_fevd_fatal.tex"),
                          "main")
write_granger_fevd_by_iv("mid",
                          file.path(here::here("results", "appendix"), "tableA_granger_fevd_mid.tex"),
                          "appendix")

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
    "Tables~\\ref{tab:system_main} and \\ref{tab:granger_fevd_fatal} ",
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