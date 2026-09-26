# [2026-09-26] Headline system-level DV switched from KL to Jensen-Shannon (JS)
# divergence. KL results are retained in results/appendix/kl_canonical/.
###############################################################################
# 18_igo_density_break.R
# Robustness: post-1980 IGO proliferation regime as a theoretical period.
#
# Motivation: the global IGO ecosystem grew explosively after ~1980
# (third-wave IGO formation, regional trade blocs, end-of-Cold-War
# institution building). If the conflict -> network-reconfiguration
# finding is driven by that secular IGO densification rather than by
# conflict per se, the coefficient on conflict should differ across
# the pre-1980 and post-1980 periods.
#
# Procedure (THEORETICAL break, not data-driven):
#   1. Hard-code BREAK_YEAR = 1980 (configurable below). Optionally
#      visualize where the densification accelerates as a sanity check;
#      do not use a data-driven break.
#   2. Build a `post` indicator = (year > BREAK_YEAR).
#   3. Fit a SINGLE interacted model per IV at each level:
#        system: js_divergence ~ IV * post + controls
#        state : churn ~ IV * post + controls | ccode + year
#      Year FE absorb the post main effect at state level; the
#      IV main effect is the pre-period slope, IV:post is the
#      differential, post-period slope = IV main + IV:post.
#      This is more powerful than splitting the sample because it
#      uses all observations to identify the controls.
#   4. Report pre-, post-, and interaction (post - pre) coefficients
#      with two-sided p-value on the interaction (the test of regime
#      heterogeneity).
#
# Outputs:
#   data/processed/robustness_igo_break_results.csv
#   code/03_figures/output/manuscript/fig_ms_igo_density_break.png
#   code/03_figures/output/manuscript/fig_ms_igo_break_coefs_kl.png
#   code/03_figures/output/manuscript/fig_ms_igo_break_coefs_churn.png
###############################################################################

library(tidyverse)
library(fixest)
library(peacesciencer)
select <- dplyr::select
set.seed(123)

# ============================================================
# Theoretical break year. Default: 1980 (start of the post-1980
# IGO proliferation period). Adjust if you want a different
# substantive cut (e.g., 1989 for end of Cold War).
# ============================================================
BREAK_YEAR <- 1980L

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

# ============================================================
# Load inputs
# ============================================================

cat("Loading inputs...\n")
dyad_data        <- readRDS(file.path(OUTPUT_DIR, "dyad_year_network_layers.rds"))
system_data_raw  <- readRDS(file.path(OUTPUT_DIR, "system_analysis_data.rds"))
state_churn_base <- readRDS(file.path(OUTPUT_DIR, "state_network_churn.rds"))

# ============================================================
# 1. Per-year IGO density series
# ============================================================
# Use the layer_igo carried in dyad_data (already normalized 0-1 within
# year, so cross-year levels reflect the share of dyads with non-zero
# shared-IGO ties). Mean of layer_igo over active dyad-years gives a
# scalar per year that tracks IGO densification.

igo_density <- dyad_data |>
  filter(!is.na(layer_igo)) |>
  group_by(year) |>
  summarise(
    igo_density = mean(layer_igo, na.rm = TRUE),
    n_pairs_with_igo = sum(layer_igo > 0, na.rm = TRUE),
    n_pairs_total    = n(),
    pct_pairs_igo    = n_pairs_with_igo / n_pairs_total,
    .groups = "drop"
  ) |>
  arrange(year)

cat(sprintf("  IGO density series: %d years (%d-%d)\n",
            nrow(igo_density),
            min(igo_density$year), max(igo_density$year)))

# ============================================================
# 2. Theoretical break year
# ============================================================
# No detection step. The break is a theoretical claim: the post-1980
# IGO proliferation regime is qualitatively different from what came
# before. The density series is shown as a descriptive sanity check
# that 1980 is a reasonable cut, not as a justification of the cut.

break_year <- BREAK_YEAR
cat(sprintf("Using theoretical break year: %d (post-period starts at %d)\n",
            break_year, break_year + 1L))

# ============================================================
# 2a. Plot the density series with the theoretical break line
# ============================================================

p_density <- igo_density |>
  ggplot(aes(x = year, y = igo_density)) +
  geom_line(color = ACCENT, linewidth = 0.7) +
  geom_vline(xintercept = break_year, color = ACCENT2,
             linetype = "dashed", linewidth = 0.6) +
  annotate("text", x = break_year,
           y = max(igo_density$igo_density) * 1.02,
           label = sprintf("Theoretical break: %d", break_year),
           color = ACCENT2, size = 3, vjust = 0, hjust = -0.05) +
  labs(
    title = "IGO network density",
    x = NULL, y = "Mean layer_igo across dyad-years"
  ) +
  theme_pres() +
  theme(plot.margin = margin(t = 10, r = 12, b = 14, l = 10))

ggsave(file.path(FIG_DIR, "fig_ms_igo_density_break.png"),
       p_density, width = 9, height = 4.5, dpi = 400, bg = "white")
cat("  -> fig_ms_igo_density_break.png\n")

# ============================================================
# 3. Re-fit main regressions on pre-break and post-break subsamples
# ============================================================

# --- System data: same setup as 01_system_level_analysis.R ---
system_data <- system_data_raw |>
  filter(!is.na(js_divergence)) |>
  arrange(year) |>
  mutate(
    log_mid_count            = log1p(mid_count),
    log_fatal_mid_count      = log1p(fatal_mid_count),
    mid_count_lag1           = dplyr::lag(mid_count,           1),
    log_mid_count_lag1       = dplyr::lag(log_mid_count,       1),
    fatal_mid_count_lag1     = dplyr::lag(fatal_mid_count,     1),
    log_fatal_mid_count_lag1 = dplyr::lag(log_fatal_mid_count, 1)
  ) |>
  filter(!is.na(mid_count_lag1))

# --- State controls (mirrors 11_igo_membership_sliding.R) ---
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
system_data <- system_data |> mutate(log_js_divergence = log1p(js_divergence))

state_data <- state_churn_base |>
  mutate(log_mid_onset_count       = log1p(mid_onset_count),
         log_fatal_mid_onset_count = log1p(fatal_mid_onset_count)) |>
  filter(!is.na(churn))

# --- Partner-war control routing by IV severity ---
# Non-fatal IVs (any MID exposure) -> partner_war_frac_weighted
# Fatal IVs (battle-death-bearing MIDs) -> partner_fatal_war_frac_weighted
# State-level churn regressions in this script use the tie-WEIGHTED partner-
# war exposure (alliance-strength-weighted share of partners at war), which is
# the harmonized robustness-script default. The unweighted variants are
# reserved for the main 02_state_level_analysis.R headline tables.
CTRL_MID   <- "partner_war_frac_weighted + democracy + autocracy + cinc + major_power + is_isolated_t1"
CTRL_FATAL <- "partner_fatal_war_frac_weighted + democracy + autocracy + cinc + major_power + is_isolated_t1"
ctrl_for_iv <- function(iv) if (grepl("fatal", iv)) CTRL_FATAL else CTRL_MID

sys_ctrl_str   <- "log1p(n_democracies) + log1p(n_autocracies) + log1p(n_states) + poly(year, 2)"
# state_ctrl_str is now routed per IV via ctrl_for_iv() inside fit_interaction().

sys_ivs <- list(
  list(col = "mid_count_lag1",           label = "MID count (t-1)"),
  list(col = "log_mid_count_lag1",       label = "log(1 + MID count) (t-1)"),
  list(col = "fatal_mid_count_lag1",     label = "Fatal MID count (t-1)"),
  list(col = "log_fatal_mid_count_lag1", label = "log(1 + Fatal MID count) (t-1)")
)
state_ivs <- list(
  list(col = "mid_onset_count",           label = "MID onset"),
  list(col = "log_mid_onset_count",       label = "log(1 + MID onset)"),
  list(col = "fatal_mid_onset_count",     label = "Fatal MID onset"),
  list(col = "log_fatal_mid_onset_count", label = "log(1 + Fatal MID onset)")
)

# ============================================================
# 4. Interaction model: IV * post on the full sample
# ============================================================
# Single model per IV per level. Pre-period slope is the IV main
# effect; post-period slope is (IV main + IV:post) via linear
# combination from coef + vcov; the interaction term is the
# differential and its p-value is the regime-heterogeneity test.

system_int <- system_data |>
  mutate(post = as.integer(year > break_year))
state_int <- state_data |>
  mutate(post = as.integer(year > break_year))

# Linear combination helper for post-period slope.
lincom_post <- function(m, iv_term, int_term) {
  cf <- coef(m); vv <- vcov(m)
  if (!iv_term %in% names(cf) || !int_term %in% names(cf)) return(NULL)
  L <- numeric(length(cf)); names(L) <- names(cf)
  L[iv_term] <- 1; L[int_term] <- 1
  est <- sum(L * cf)
  var <- as.numeric(t(L) %*% vv %*% L)
  list(est = est, se = sqrt(var))
}

fit_interaction <- function(level, iv_col, iv_lab, df) {
  if (level == "system") {
    f <- as.formula(sprintf("log_js_divergence ~ %s * post + %s",
                              iv_col, sys_ctrl_str))
    m <- tryCatch(lm(f, data = df), error = function(e) NULL)
  } else {
    f <- as.formula(sprintf("log_churn ~ %s * post + %s | ccode + year",
                              iv_col, ctrl_for_iv(iv_col)))
    m <- tryCatch(feols(f, data = df, cluster = ~ ccode + year),
                   error = function(e) NULL)
  }
  if (is.null(m)) return(NULL)
  cf <- coef(m); vv <- vcov(m)
  iv_term  <- iv_col
  int_term <- paste0(iv_col, ":post")
  if (!iv_term %in% names(cf) || !int_term %in% names(cf)) return(NULL)

  pre_est <- unname(cf[iv_term])
  pre_se  <- sqrt(diag(vv)[iv_term])

  pl <- lincom_post(m, iv_term, int_term)
  post_est <- pl$est; post_se <- pl$se

  int_est <- unname(cf[int_term])
  int_se  <- sqrt(diag(vv)[int_term])
  int_p   <- 2 * (1 - pnorm(abs(int_est / int_se)))

  bind_rows(
    tibble(regime = "pre",  estimate = pre_est,  se = pre_se),
    tibble(regime = "post", estimate = post_est, se = post_se),
    tibble(regime = "diff", estimate = int_est,  se = int_se)
  ) |>
    mutate(level = level, iv = iv_col, iv_label = iv_lab,
           int_p = int_p,
           ci_low  = estimate - 1.96 * se,
           ci_high = estimate + 1.96 * se,
           n_obs   = nobs(m))
}

cat(sprintf("\n--- Fitting interaction models (break = %d) ---\n", break_year))
all_results <- bind_rows(
  purrr::map_dfr(sys_ivs, function(iv)
    fit_interaction("system", iv$col, iv$label, system_int)),
  purrr::map_dfr(state_ivs, function(iv)
    fit_interaction("state",  iv$col, iv$label, state_int))
)

iv_levels_all <- c(sapply(sys_ivs,   `[[`, "label"),
                    sapply(state_ivs, `[[`, "label"))

all_results <- all_results |>
  mutate(iv_label    = factor(iv_label, levels = iv_levels_all),
         level_label = factor(level, levels = c("system", "state"),
                                labels = c("System (JS)", "State (churn)")),
         regime      = factor(regime,
                                levels = c("pre", "post", "diff"),
                                labels = c(sprintf("≤ %d",  break_year),
                                           sprintf("> %d",  break_year),
                                           "Post - Pre (interaction)")))

write_csv(all_results,
          file.path(OUTPUT_DIR, "robustness_igo_break_results.csv"))
cat(sprintf("\nSaved: robustness_igo_break_results.csv (%d rows)\n",
            nrow(all_results)))

# ============================================================
# 5. Coefficient-comparison plot
# ============================================================

panel_theme <- theme_pres() +
  theme(panel.spacing = unit(1.2, "lines"),
        panel.spacing.x = unit(1.5, "lines"),
        strip.placement = "outside",
        strip.text = element_text(color = "black", size = 9, face = "bold",
                                    margin = margin(4, 4, 4, 4)),
        plot.margin = margin(t = 10, r = 12, b = 14, l = 10))

iv_sys_levels   <- sapply(sys_ivs,   `[[`, "label")
iv_state_levels <- sapply(state_ivs, `[[`, "label")

# Pre/post point estimates (drop the 'diff' rows for the side-by-side plot;
# they remain in the CSV for inspection). Annotate each facet with the
# interaction p-value -- the test of regime heterogeneity.
plot_df <- all_results |>
  filter(regime != "Post - Pre (interaction)")

ann_df <- all_results |>
  filter(regime == "Post - Pre (interaction)") |>
  select(level, iv_label, int_p) |>
  distinct() |>
  mutate(p_label = ifelse(int_p < 0.001,
                          "p[interaction] < 0.001",
                          sprintf("p[interaction] = %.3f", int_p)))

color_scale <- scale_color_manual(
  values = c("gray30", ACCENT, ACCENT2)[1:2] |>
              setNames(c(sprintf("≤ %d", break_year),
                          sprintf("> %d", break_year))),
  guide = "none")

# --- System (KL) ---
plot_sys <- plot_df |>
  filter(level == "system") |>
  mutate(iv_label = factor(iv_label, levels = iv_sys_levels))
ann_sys <- ann_df |>
  filter(level == "system") |>
  mutate(iv_label = factor(iv_label, levels = iv_sys_levels))

p_coef_sys <- ggplot(plot_sys,
                       aes(x = regime, y = estimate, color = regime)) +
  geom_hline(yintercept = 0, color = "gray55", linewidth = 0.3) +
  geom_errorbar(aes(ymin = ci_low, ymax = ci_high),
                width = 0.15, linewidth = 0.8) +
  geom_point(size = 2.8) +
  color_scale +
  geom_text(data = ann_sys,
            aes(x = 1.5, y = Inf, label = p_label),
            inherit.aes = FALSE,
            vjust = 1.5, hjust = 0.5, size = 3, color = "gray25") +
  facet_wrap(~ iv_label, scales = "free_y", ncol = 2,
             labeller = label_wrap_gen(width = 22)) +
  labs(
    title = sprintf("System: pre/post-%d", break_year),
    x = NULL, y = "Coefficient (95% CI)"
  ) +
  panel_theme

ggsave(file.path(FIG_DIR, "fig_ms_igo_break_coefs_kl.png"),
       p_coef_sys, width = 9.5, height = 5.5, dpi = 400, bg = "white")
cat("  -> fig_ms_igo_break_coefs_kl.png\n")

# --- State (churn) ---
plot_st <- plot_df |>
  filter(level == "state") |>
  mutate(iv_label = factor(iv_label, levels = iv_state_levels))
ann_st <- ann_df |>
  filter(level == "state") |>
  mutate(iv_label = factor(iv_label, levels = iv_state_levels))

p_coef_st <- ggplot(plot_st,
                      aes(x = regime, y = estimate, color = regime)) +
  geom_hline(yintercept = 0, color = "gray55", linewidth = 0.3) +
  geom_errorbar(aes(ymin = ci_low, ymax = ci_high),
                width = 0.15, linewidth = 0.8) +
  geom_point(size = 2.8) +
  color_scale +
  geom_text(data = ann_st,
            aes(x = 1.5, y = Inf, label = p_label),
            inherit.aes = FALSE,
            vjust = 1.5, hjust = 0.5, size = 3, color = "gray25") +
  facet_wrap(~ iv_label, scales = "free_y", ncol = 2,
             labeller = label_wrap_gen(width = 22)) +
  labs(
    title = sprintf("State: pre/post-%d", break_year),
    x = NULL, y = "Coefficient (95% CI)"
  ) +
  panel_theme

ggsave(file.path(FIG_DIR, "fig_ms_igo_break_coefs_churn.png"),
       p_coef_st, width = 9.5, height = 5.5, dpi = 400, bg = "white")
cat("  -> fig_ms_igo_break_coefs_churn.png\n")

cat("\nDone.\n")
