###############################################################################
# 20_oos_prediction_period.R
# Robustness: out-of-sample prediction by period.
#
# Motivation. The outline flags a known scope condition: KL/JS divergence
# coefficients flip sign in the post-1945 period and weaken further
# post-1991. If the war -> network reconfiguration relationship genuinely
# breaks down in the postwar institutional order, then a model fit on
# the early period should fail to predict the later period -- and
# specifically, errors should be biased systematically (over-predicting
# reconfiguration when it didn't happen).
#
# Procedure. For each cutoff c in {1945, 1980, 1990}:
#
#   SYSTEM
#     Fit:    js_divergence (raw, as of 2026-10-02) ~ <iv>_lag1 + log1p(n_dem) + log1p(n_aut)
#             on year <= c
#     Predict: year > c
#     Report: RMSE, MAE, bias (mean residual), R^2 vs realized series
#
#   STATE
#     Asymptotic TWFE cannot extrapolate state or year FE to unseen
#     state-years. We report two variants per cutoff to bracket the
#     answer:
#       (a) POOLED OLS (no FE, same controls). Cleanest OOS exercise --
#           does the conditional-mean coefficient on conflict generalize?
#       (b) STATE-FE only (no year FE), training-set state FE re-used at
#           prediction. Tests generalization conditional on time-invariant
#           state heterogeneity.
#     Either way, only state-years whose state appeared in the training
#     window get a prediction in (b); (a) predicts everywhere.
#
# All four IVs are run separately (mid_count, log, fatal, log fatal at
# the system level; mid_onset_count and friends at the state level).
#
# Outputs:
#   data/processed/robustness_oos_prediction_metrics.csv
#   data/processed/robustness_oos_prediction_series.csv
#   code/03_figures/output/manuscript/fig_ms_oos_metrics.png
#   code/03_figures/output/manuscript/fig_ms_oos_kl_series.png
###############################################################################

library(tidyverse)
library(fixest)
library(peacesciencer)
select <- dplyr::select
set.seed(123)

OUTPUT_DIR <- here::here("data", "processed")
FIG_DIR    <- here::here("code", "03_figures", "output", "manuscript")
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

ACCENT  <- "#2C6CB0"
ACCENT2 <- "#B03366"

theme_pres <- function(base_size = 12) {
  theme_bw(base_size = base_size) %+replace%
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "gray90", linewidth = 0.3),
      plot.title   = element_text(size = base_size + 2, face = "bold",
                                  hjust = 0, margin = margin(b = 8)),
      plot.caption = element_text(size = base_size - 2, color = "gray30",
                                  hjust = 1, margin = margin(t = 8)),
      axis.text    = element_text(color = "black", size = base_size - 1),
      axis.title   = element_text(color = "black", size = base_size),
      strip.background = element_rect(fill = "gray92", color = "gray70"),
      strip.text   = element_text(face = "bold", size = base_size,
                                  color = "black"),
      legend.position   = "bottom",
      plot.margin       = margin(10, 12, 8, 10)
    )
}

CUTOFFS <- c(1945L, 1980L, 1990L)

# ============================================================
# Load + prepare data
# ============================================================

cat("Loading inputs...\n")

dyad_data        <- readRDS(file.path(OUTPUT_DIR, "dyad_year_network_layers.rds"))
system_data_raw  <- readRDS(file.path(OUTPUT_DIR, "system_analysis_data.rds"))
state_churn_base <- readRDS(file.path(OUTPUT_DIR, "state_network_churn.rds"))

system_data <- system_data_raw |>
  arrange(year) |>
  # [2026-10-12] conflict lags are built BEFORE dropping years without JS, so 1817 keeps its 1816 conflict
  (\(d) { stopifnot(all(diff(d$year) == 1)); d })() |>   # consecutive years
  mutate(
    log_mid_count            = log1p(mid_count),
    log_fatal_mid_count      = log1p(fatal_mid_count),
    mid_count_lag1           = dplyr::lag(mid_count,           1),
    log_mid_count_lag1       = dplyr::lag(log_mid_count,       1),
    fatal_mid_count_lag1     = dplyr::lag(fatal_mid_count,     1),
    log_fatal_mid_count_lag1 = dplyr::lag(log_fatal_mid_count, 1)
  ) |>
  filter(!is.na(mid_count_lag1), !is.na(js_divergence))

state_controls_needed <- c("partner_war_frac", "democracy", "autocracy",
                            "cinc", "major_power")
if (!all(state_controls_needed %in% names(state_churn_base))) {
  cat("Preparing state-level controls...\n")

  if (!"main_weight" %in% names(dyad_data)) {
    ML <- c("layer_offensive", "layer_defensive", "layer_igo")
    mat <- as.matrix(dyad_data[, ML, drop = FALSE]); mat[is.na(mat)] <- 0
    dyad_data$main_raw <- rowSums(mat)
    dyad_data <- dyad_data |>
      group_by(year) |>
      mutate(main_weight = if (max(main_raw, na.rm = TRUE) > 0) {
        main_raw / max(main_raw, na.rm = TRUE)
      } else 0) |>
      ungroup()
    dyad_data$main_raw <- NULL
  }

  state_attrs <- dyad_data |>
    select(ccode = ccode1, year, cinc = cinc1, major_power = cowmaj1) |>
    distinct(ccode, year, .keep_all = TRUE) |>
    mutate(major_power = as.integer(coalesce(major_power, 0L) == 1))
  state_churn_base <- state_churn_base |>
    left_join(state_attrs, by = c("ccode", "year"))

  state_years <- create_stateyears(
    subset_years = min(state_churn_base$year):max(state_churn_base$year))
  attr(state_years, "ps_data_type") <- "state_year"
  attr(state_years, "ps_system")    <- "cow"
  state_years <- state_years |> add_democracy()
  polity <- state_years |>
    select(ccode, year, polity2) |>
    mutate(democracy = as.integer(polity2 >= 6),
           autocracy = as.integer(polity2 <= -6))
  state_churn_base <- state_churn_base |>
    left_join(polity, by = c("ccode", "year"))

  partner_sets <- dyad_data |>
    filter(main_weight > 0) |>
    select(ccode1, ccode2, year, w = main_weight) |>
    bind_rows(dyad_data |> filter(main_weight > 0) |>
                select(ccode1 = ccode2, ccode2 = ccode1, year, w = main_weight))
  at_war <- dyad_data |>
    filter(cowmidongoing == 1) |>
    select(ccode1, ccode2, year) |>
    pivot_longer(c(ccode1, ccode2), values_to = "ccode") |>
    distinct(ccode, year) |>
    mutate(at_war = 1L)
  at_war_fatal <- dyad_data |>
    filter(cowmidfatal == 1) |>
    select(ccode1, ccode2, year) |>
    pivot_longer(c(ccode1, ccode2), values_to = "ccode") |>
    distinct(ccode, year) |>
    mutate(at_war_fatal = 1L)
  pwe <- partner_sets |>
    mutate(year_next = year + 1L) |>
    left_join(at_war,       by = c("ccode2" = "ccode", "year_next" = "year")) |>
    left_join(at_war_fatal, by = c("ccode2" = "ccode", "year_next" = "year")) |>
    mutate(
      at_war       = replace_na(at_war, 0L),
      at_war_fatal = replace_na(at_war_fatal, 0L)
    ) |>
    group_by(ccode1, year_next) |>
    summarise(
      partner_war_frac                = mean(at_war),
      partner_war_frac_weighted       = sum(at_war * w) / sum(w),
      partner_fatal_war_frac          = mean(at_war_fatal),
      partner_fatal_war_frac_weighted = sum(at_war_fatal * w) / sum(w),
      .groups = "drop"
    ) |>
    rename(ccode = ccode1, year = year_next)
  state_churn_base <- state_churn_base |>
    left_join(pwe, by = c("ccode", "year")) |>
    mutate(partner_war_frac                = replace_na(partner_war_frac, 0),
           partner_war_frac_weighted       = replace_na(partner_war_frac_weighted, 0),
           partner_fatal_war_frac          = replace_na(partner_fatal_war_frac, 0),
           partner_fatal_war_frac_weighted = replace_na(partner_fatal_war_frac_weighted, 0),
           major_power = replace_na(major_power, 0))
}

state_churn_base <- state_churn_base |> mutate(log_churn = log1p(churn))
# [2026-10-12] Headline state IVs are onset counts at t-1, taken from the full COW state-year universe in
# dyad_data (same construction as 01_data_cleaning/02_compute_network_measures.R and 02). The churn panel
# holds transition rows only, so lagging the previous ROW would drop every state's first transition year.
# A t-1 value is NA only if the state was not in the system in t-1 (after a COW membership gap).
onset_univ <- dplyr::bind_rows(dyad_data |> dplyr::distinct(ccode = ccode1, year),
                               dyad_data |> dplyr::distinct(ccode = ccode2, year)) |>
  dplyr::distinct() |>
  dplyr::left_join(
    dyad_data |>
      dplyr::filter(cowmidongoing == 1, ccode1 < ccode2) |>
      dplyr::select(ccode1, ccode2, year, cowmidfatal, cowmidonset) |>
      tidyr::pivot_longer(c(ccode1, ccode2), values_to = "ccode") |>
      dplyr::group_by(ccode, year) |>
      dplyr::summarise(u_onset = sum(cowmidonset == 1, na.rm = TRUE),
                       u_fatal_onset = sum(cowmidfatal == 1 & cowmidonset == 1, na.rm = TRUE),
                       .groups = "drop"),
    by = c("ccode", "year")) |>
  dplyr::mutate(u_onset = as.numeric(tidyr::replace_na(u_onset, 0L)),
                u_fatal_onset = as.numeric(tidyr::replace_na(u_fatal_onset, 0L)))
chk_univ <- dplyr::inner_join(state_churn_base |> dplyr::select(ccode, year, mid_onset_count, fatal_mid_onset_count),
                              onset_univ, by = c("ccode", "year"))
stopifnot(nrow(chk_univ) == nrow(state_churn_base), all(chk_univ$mid_onset_count == chk_univ$u_onset),
          all(chk_univ$fatal_mid_onset_count == chk_univ$u_fatal_onset))
rm(chk_univ)
state_churn_base <- state_churn_base |>
  dplyr::select(-dplyr::any_of(c("onset_lag1", "fatal_onset_lag1", "log_onset_lag1", "log_fatal_onset_lag1"))) |>
  dplyr::left_join(onset_univ |>
                     dplyr::transmute(ccode, year = year + 1L,
                                      onset_lag1 = u_onset, fatal_onset_lag1 = u_fatal_onset,
                                      log_onset_lag1 = log1p(u_onset), log_fatal_onset_lag1 = log1p(u_fatal_onset)),
                   by = c("ccode", "year"))
stopifnot(!anyDuplicated(state_churn_base[c("ccode", "year")]))
system_data <- system_data |> mutate(log_kl_divergence = log1p(kl_divergence))

state_data <- state_churn_base |>
  mutate(log_mid_onset_count       = log1p(mid_onset_count),
         log_fatal_mid_onset_count = log1p(fatal_mid_onset_count)) |>
  filter(!is.na(churn))

# ============================================================
# IVs and control strings
# ============================================================

sys_ivs <- list(
  list(col = "mid_count_lag1",           label = "MID count (t-1)"),
  list(col = "log_mid_count_lag1",       label = "log(1 + MID count) (t-1)"),
  list(col = "fatal_mid_count_lag1",     label = "Fatal MID count (t-1)"),
  list(col = "log_fatal_mid_count_lag1", label = "log(1 + Fatal MID count) (t-1)")
)
# [2026-10-12] state onset IVs at t-1 (headline timing)
state_ivs <- list(
  list(col = "onset_lag1",           label = "MID onset (t-1)"),
  list(col = "log_onset_lag1",       label = "log(1 + MID onset) (t-1)"),
  list(col = "fatal_onset_lag1",     label = "Fatal MID onset (t-1)"),
  list(col = "log_fatal_onset_lag1", label = "log(1 + Fatal MID onset) (t-1)")
)
# --- State-level control set ---
# [2026-10-12] partner war exposure removed from the control set (H3 dropped from the paper).
# CTRL_MID and CTRL_FATAL are now identical (democracy + autocracy + cinc +
# major_power + is_isolated_t1); both names and ctrl_for_iv() are kept so the
# call sites are unchanged.
CTRL_MID   <- "democracy + autocracy + cinc + major_power + is_isolated_t1"
CTRL_FATAL <- "democracy + autocracy + cinc + major_power + is_isolated_t1"
ctrl_for_iv <- function(iv) if (grepl("fatal", iv)) CTRL_FATAL else CTRL_MID

sys_ctrl_str   <- "log1p(n_democracies) + log1p(n_autocracies) + log1p(n_states) + poly(year, 2)"
# state_ctrl_str is now routed per IV via ctrl_for_iv() inside fit_predict_state().

# ============================================================
# Helpers
# ============================================================

pred_metrics <- function(actual, predicted) {
  ok <- !is.na(actual) & !is.na(predicted)
  if (sum(ok) < 2) {
    return(tibble(rmse = NA_real_, mae = NA_real_, bias = NA_real_,
                  r2 = NA_real_, n = sum(ok)))
  }
  a <- actual[ok]; p <- predicted[ok]
  resid <- a - p
  ss_res <- sum(resid^2)
  ss_tot <- sum((a - mean(a))^2)
  r2 <- if (ss_tot > 0) 1 - ss_res / ss_tot else NA_real_
  tibble(rmse = sqrt(mean(resid^2)),
         mae  = mean(abs(resid)),
         bias = mean(resid),
         r2   = r2,
         n    = sum(ok))
}

# ---------- System level ----------
fit_predict_system <- function(iv_col, iv_label, cutoff) {
  train <- system_data |> filter(year <= cutoff)
  test  <- system_data |> filter(year >  cutoff)
  if (nrow(train) < 30 || nrow(test) < 5) return(NULL)

  f <- as.formula(paste0("js_divergence ~ ", iv_col, " + ", sys_ctrl_str))
  m <- tryCatch(lm(f, data = train), error = function(e) NULL)
  if (is.null(m)) return(NULL)

  pred_test  <- suppressWarnings(predict(m, newdata = test))
  pred_train <- suppressWarnings(predict(m, newdata = train))

  cf <- coef(m); vv <- vcov(m)
  est <- cf[iv_col]; se <- sqrt(diag(vv)[iv_col])

  metrics <- bind_rows(
    pred_metrics(train$js_divergence, pred_train) |>
      mutate(sample = "train", .before = 1),
    pred_metrics(test$js_divergence,  pred_test)  |>
      mutate(sample = "test",  .before = 1)
  ) |>
    mutate(level = "system", iv = iv_col, iv_label = iv_label,
           cutoff = cutoff, model = "ols",
           estimate = unname(est), se = unname(se),
           ci_low  = unname(est - 1.96 * se),
           ci_high = unname(est + 1.96 * se))

  series <- bind_rows(
    train |> select(year) |> mutate(actual    = train$js_divergence,
                                     predicted = pred_train,
                                     sample    = "train"),
    test  |> select(year) |> mutate(actual    = test$js_divergence,
                                     predicted = pred_test,
                                     sample    = "test")
  ) |>
    mutate(level = "system", iv = iv_col, iv_label = iv_label,
           cutoff = cutoff, model = "ols")

  list(metrics = metrics, series = series)
}

# ---------- State level ----------
# (a) Pooled OLS, no FE
# (b) State FE only (one-way), via fixest. Predict on test using
#     newdata; fixest::predict handles seen-state extrapolation. For
#     test states absent from training, prediction is NA.
fit_predict_state <- function(iv_col, iv_label, cutoff) {
  train <- state_data |> filter(year <= cutoff)
  test  <- state_data |> filter(year >  cutoff)
  if (nrow(train) < 100 || nrow(test) < 50) return(NULL)

  out_metrics <- list()
  out_series  <- list()

  # (a) Pooled OLS
  f_a <- as.formula(paste0("churn ~ ", iv_col, " + ", ctrl_for_iv(iv_col)))
  m_a <- tryCatch(lm(f_a, data = train), error = function(e) NULL)
  if (!is.null(m_a)) {
    p_test  <- suppressWarnings(predict(m_a, newdata = test))
    p_train <- suppressWarnings(predict(m_a, newdata = train))

    cf <- coef(m_a); vv <- vcov(m_a)
    est <- cf[iv_col]; se <- sqrt(diag(vv)[iv_col])

    out_metrics[["pooled"]] <- bind_rows(
      pred_metrics(train$churn, p_train) |>
        mutate(sample = "train", .before = 1),
      pred_metrics(test$churn,  p_test)  |>
        mutate(sample = "test",  .before = 1)
    ) |>
      mutate(level = "state", iv = iv_col, iv_label = iv_label,
             cutoff = cutoff, model = "pooled_ols",
             estimate = unname(est), se = unname(se),
             ci_low  = unname(est - 1.96 * se),
             ci_high = unname(est + 1.96 * se))
  }

  # (b) State FE only
  f_b <- as.formula(paste0("churn ~ ", iv_col, " + ",
                            ctrl_for_iv(iv_col), " | ccode"))
  m_b <- tryCatch(feols(f_b, data = train,
                          cluster = ~ ccode),
                   error = function(e) NULL)
  if (!is.null(m_b)) {
    # fixest predict returns NA for state FE not in training
    p_test  <- tryCatch(predict(m_b, newdata = test),
                         error = function(e) rep(NA_real_, nrow(test)))
    p_train <- tryCatch(predict(m_b, newdata = train),
                         error = function(e) rep(NA_real_, nrow(train)))

    cf <- coef(m_b)
    if (iv_col %in% names(cf)) {
      vv <- vcov(m_b, cluster = ~ ccode)
      est <- cf[iv_col]; se <- sqrt(diag(vv)[iv_col])

      out_metrics[["state_fe"]] <- bind_rows(
        pred_metrics(train$churn, p_train) |>
          mutate(sample = "train", .before = 1),
        pred_metrics(test$churn,  p_test)  |>
          mutate(sample = "test",  .before = 1)
      ) |>
        mutate(level = "state", iv = iv_col, iv_label = iv_label,
               cutoff = cutoff, model = "state_fe_only",
               estimate = unname(est), se = unname(se),
               ci_low  = unname(est - 1.96 * se),
               ci_high = unname(est + 1.96 * se))
    }
  }

  list(metrics = bind_rows(out_metrics),
       series  = NULL)  # state-level series too large; skip raw saving
}

# ============================================================
# Run all combinations
# ============================================================

cat("\n--- System OOS ---\n")
sys_runs <- expand_grid(iv_idx = seq_along(sys_ivs), cutoff = CUTOFFS) |>
  pmap(function(iv_idx, cutoff) {
    iv <- sys_ivs[[iv_idx]]
    cat(sprintf("  cutoff=%d  iv=%s\n", cutoff, iv$col))
    fit_predict_system(iv$col, iv$label, cutoff)
  }) |> compact()

cat("\n--- State OOS ---\n")
st_runs <- expand_grid(iv_idx = seq_along(state_ivs), cutoff = CUTOFFS) |>
  pmap(function(iv_idx, cutoff) {
    iv <- state_ivs[[iv_idx]]
    cat(sprintf("  cutoff=%d  iv=%s\n", cutoff, iv$col))
    fit_predict_state(iv$col, iv$label, cutoff)
  }) |> compact()

metrics <- bind_rows(
  purrr::map_dfr(sys_runs, "metrics"),
  purrr::map_dfr(st_runs, "metrics")
)
write_csv(metrics, file.path(OUTPUT_DIR,
                              "robustness_oos_prediction_metrics.csv"))
cat(sprintf("\nSaved metrics: %d rows\n", nrow(metrics)))

sys_series <- purrr::map_dfr(sys_runs, "series")
write_csv(sys_series, file.path(OUTPUT_DIR,
                                  "robustness_oos_prediction_series.csv"))
cat(sprintf("Saved system series: %d rows\n", nrow(sys_series)))

# ============================================================
# Plots
# ============================================================

cat("\n--- Plotting ---\n")

iv_levels_sys <- sapply(sys_ivs,   `[[`, "label")
iv_levels_st  <- sapply(state_ivs, `[[`, "label")

panel_theme <- theme_pres() +
  theme(panel.spacing = unit(1.2, "lines"),
        panel.spacing.x = unit(1.5, "lines"),
        strip.placement = "outside",
        strip.text = element_text(color = "black", size = 9, face = "bold",
                                    margin = margin(4, 4, 4, 4)),
        axis.text.x = element_text(angle = 0),
        plot.margin = margin(t = 10, r = 12, b = 14, l = 10))

# (1) Train vs test RMSE bar chart -- separate KL vs churn plots
metrics_plot <- metrics |>
  mutate(cutoff = factor(cutoff, levels = CUTOFFS),
         sample = factor(sample, levels = c("train", "test")))

# --- System (JS) ---
metrics_sys <- metrics_plot |>
  filter(level == "system") |>
  mutate(iv_label = factor(iv_label, levels = iv_levels_sys))

p_metrics_sys <- ggplot(metrics_sys,
                         aes(x = cutoff, y = rmse, fill = sample)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.7) +
  scale_fill_manual(values = c("train" = ACCENT, "test" = ACCENT2)) +
  facet_wrap(~ iv_label, scales = "free_y", ncol = 2,
             labeller = label_wrap_gen(width = 22)) +
  labs(
    title = "System: out-of-sample RMSE",
    x = "Training cutoff (train: <= cutoff, test: > cutoff)",
    y = "RMSE",
    fill = NULL
  ) +
  panel_theme

ggsave(file.path(FIG_DIR, "fig_ms_oos_metrics_kl.png"),
       p_metrics_sys, width = 9.5, height = 6, dpi = 400, bg = "white")
cat("  -> fig_ms_oos_metrics_kl.png\n")

# --- State (churn) ---
metrics_st <- metrics_plot |>
  filter(level == "state") |>
  mutate(iv_label = factor(iv_label, levels = iv_levels_st))

p_metrics_st <- ggplot(metrics_st,
                        aes(x = cutoff, y = rmse,
                            fill = interaction(sample, model, sep = " - "))) +
  geom_col(position = position_dodge(width = 0.85), width = 0.8) +
  scale_fill_manual(values = c(
    "train - pooled_ols"      = ACCENT,
    "test - pooled_ols"       = ACCENT2,
    "train - state_fe_only"   = "#7FA9D9",
    "test - state_fe_only"    = "#D98AA6"
  )) +
  facet_wrap(~ iv_label, scales = "free_y", ncol = 2,
             labeller = label_wrap_gen(width = 22)) +
  labs(
    title = "State: out-of-sample RMSE",
    x = "Training cutoff (train: <= cutoff, test: > cutoff)",
    y = "RMSE",
    fill = NULL
  ) +
  panel_theme +
  theme(legend.position = "bottom")

ggsave(file.path(FIG_DIR, "fig_ms_oos_metrics_churn.png"),
       p_metrics_st, width = 10, height = 6, dpi = 400, bg = "white")
cat("  -> fig_ms_oos_metrics_churn.png\n")

# (2) System JS: actual vs predicted time series, one panel per cutoff
p_kl <- sys_series |>
  filter(iv == "log_fatal_mid_count_lag1") |>   # representative IV [2026-10-12: fatal headline]
  pivot_longer(c(actual, predicted),
               names_to = "series", values_to = "value") |>
  mutate(cutoff = factor(cutoff, levels = CUTOFFS)) |>
  ggplot(aes(x = year, y = value, color = series, linetype = series)) +
  geom_vline(aes(xintercept = as.numeric(as.character(cutoff))),
             linetype = "dashed", color = "gray40", linewidth = 0.4) +
  geom_line(linewidth = 0.6) +
  scale_color_manual(values = c("actual" = "black", "predicted" = ACCENT)) +
  scale_linetype_manual(values = c("actual" = "solid", "predicted" = "dashed")) +
  facet_wrap(~ cutoff, ncol = 1,
              labeller = labeller(cutoff = function(x) sprintf("Train <= %s", x))) +
  labs(
    title = "JS: actual vs predicted",
    x = NULL, y = "JS divergence",
    color = NULL, linetype = NULL
  ) +
  theme_pres() +
  theme(strip.text = element_text(size = 9),
        panel.spacing = unit(1.2, "lines"),
        panel.spacing.y = unit(1.5, "lines"),
        plot.margin = margin(t = 10, r = 12, b = 14, l = 10))

ggsave(file.path(FIG_DIR, "fig_ms_oos_kl_series.png"),
       p_kl, width = 9, height = 7, dpi = 400, bg = "white")
cat("  -> fig_ms_oos_kl_series.png\n")

cat("\nDone.\n")
