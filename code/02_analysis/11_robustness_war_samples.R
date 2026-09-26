# [2026-09-26] Headline system-level DV switched from KL to Jensen-Shannon (JS)
# divergence. KL results are retained in results/appendix/kl_canonical/.
###############################################################################
# 12_robustness_war_samples.R
# Robustness: refit the preliminary OLS (system) and TWFE (state) models on
# two alternative samples:
#
#   (A) DROP WWI + WWII years  (year in 1914-1918 or 1939-1945)
#       Tests whether the world-wars period is driving the headline result.
#
#   (B) LARGE CONFLICTS ONLY
#       System level: restrict to years with fatal_mid_count >= 5
#       State level:  restrict to state-years where the state participated
#                     in a fatal MID (fatal_mid_onset_count > 0)
#       Tests whether the effect survives when we *only* look at observations
#       with substantive conflict content.
#
# Both samples are re-run alongside the full-sample baseline so coefficients
# are comparable side-by-side.
#
# Outputs:
#   data/processed/robustness_war_samples.csv
#   code/03_figures/output/presentation/pres_fig16_war_samples.png
###############################################################################

library(tidyverse)
library(fixest)
library(peacesciencer)
select <- dplyr::select

set.seed(123)

OUTPUT_DIR <- here::here("data", "processed")
FIG_DIR    <- here::here("code", "03_figures", "output", "manuscript")
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------
# Manuscript theme — white background, publication-suitable
# ------------------------------------------------------------
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
      axis.ticks   = element_line(color = "gray40"),
      strip.background = element_rect(fill = "gray92", color = "gray70"),
      strip.text   = element_text(face = "bold", size = base_size,
                                    color = "black"),
      legend.position   = "bottom",
      legend.background = element_blank(),
      legend.text       = element_text(size = base_size - 1),
      legend.title      = element_text(size = base_size, face = "bold"),
      plot.margin       = margin(10, 12, 8, 10)
    )
}

save_pres <- function(filename, plot, width = 7, height = 5, dpi = 400) {
  ggsave(file.path(FIG_DIR, filename), plot,
         width = width, height = height, dpi = dpi, bg = "white")
  cat(sprintf("  -> %s\n", filename))
}

# ============================================================
# Config — thresholds for large-conflict sample
# ============================================================

# System-level: sweep the minimum fatal-MID threshold used to define the
# "large conflicts only" subsample. 3 is permissive (most years retained);
# 5 is moderate; 10 is strict (only the most conflict-intensive years).
SYSTEM_FATAL_MINS <- c(3, 5, 10)

# World-war year ranges
WW1_YEARS <- 1914:1918
WW2_YEARS <- 1939:1945
WAR_YEARS <- c(WW1_YEARS, WW2_YEARS)

# ============================================================
# Load data
# ============================================================

cat("Loading data...\n")

# Match main-text baseline: lagged IVs, plain OLS, no post_1945.
system_data <- readRDS(file.path(OUTPUT_DIR, "system_analysis_data.rds")) |>
  filter(!is.na(js_divergence)) |>
  arrange(year) |>
  mutate(
    log_js_divergence    = log1p(js_divergence),
    log_mid_count        = log1p(mid_count),
    log_fatal_mid_count  = log1p(fatal_mid_count),
    mid_count_lag1           = dplyr::lag(mid_count,           1),
    log_mid_count_lag1       = dplyr::lag(log_mid_count,       1),
    fatal_mid_count_lag1     = dplyr::lag(fatal_mid_count,     1),
    log_fatal_mid_count_lag1 = dplyr::lag(log_fatal_mid_count, 1)
  ) |>
  filter(!is.na(mid_count_lag1))

state_churn <- readRDS(file.path(OUTPUT_DIR, "state_network_churn.rds"))
dyad_data   <- readRDS(file.path(OUTPUT_DIR, "dyad_year_network_layers.rds"))

# Attach state-level controls if the RDS doesn't already have them
needs_prep <- !all(c("cinc", "major_power", "democracy", "autocracy",
                      "partner_war_frac") %in% names(state_churn))
if (needs_prep) {
  cat("  Preparing state-level controls...\n")
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
  state_churn <- state_churn |> left_join(state_attrs, by = c("ccode", "year"))

  state_years <- create_stateyears(
    subset_years = min(state_churn$year):max(state_churn$year))
  attr(state_years, "ps_data_type") <- "state_year"
  attr(state_years, "ps_system")    <- "cow"
  state_years <- state_years |> add_democracy()
  polity <- state_years |>
    select(ccode, year, polity2) |>
    mutate(democracy = as.integer(polity2 >= 6),
           autocracy = as.integer(polity2 <= -6))
  state_churn <- state_churn |> left_join(polity, by = c("ccode", "year"))

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
  state_churn <- state_churn |>
    left_join(pwe, by = c("ccode", "year")) |>
    mutate(partner_war_frac                = replace_na(partner_war_frac, 0),
           partner_war_frac_weighted       = replace_na(partner_war_frac_weighted, 0),
           partner_fatal_war_frac          = replace_na(partner_fatal_war_frac, 0),
           partner_fatal_war_frac_weighted = replace_na(partner_fatal_war_frac_weighted, 0),
           major_power = replace_na(major_power, 0))
}

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

state_churn <- state_churn |>
  mutate(
    log_churn                 = log1p(churn),
    log_mid_onset_count       = log1p(mid_onset_count),
    log_fatal_mid_onset_count = log1p(fatal_mid_onset_count)
  )

# ============================================================
# Model specs and fit helpers — four IVs, each fit separately
# ============================================================

# IV sweep (each assessed as its own spec)
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

sys_controls_str <- "log1p(n_democracies) + log1p(n_autocracies) + log1p(n_states) + poly(year, 2)"
# state_controls is now routed per IV via ctrl_for_iv() inside fit_state().

fit_sys <- function(d, sample_label) {
  purrr::map_dfr(sys_ivs, function(iv) {
    iv_col <- iv$col; iv_lab <- iv$label
    f <- as.formula(paste0("log_js_divergence ~ ", iv_col, " + ", sys_controls_str))
    m <- tryCatch(lm(f, data = d), error = function(e) NULL)
    if (is.null(m)) return(NULL)
    cf <- coef(m); vv <- vcov(m)
    tibble(level = "system", sample = sample_label,
           iv = iv_col, iv_label = iv_lab,
           estimate = unname(cf[iv_col]),
           se       = sqrt(diag(vv)[iv_col]),
           n_obs    = nobs(m))
  })
}

fit_state <- function(d, sample_label) {
  purrr::map_dfr(state_ivs, function(iv) {
    iv_col <- iv$col; iv_lab <- iv$label
    f <- as.formula(paste0("log_churn ~ ", iv_col, " + ",
                            ctrl_for_iv(iv_col), " | ccode + year"))
    m <- tryCatch(feols(f, data = d, cluster = ~ ccode + year),
                   error = function(e) NULL)
    if (is.null(m)) return(NULL)
    cf <- coef(m); vv <- vcov(m)
    tibble(level = "state", sample = sample_label,
           iv = iv_col, iv_label = iv_lab,
           estimate = unname(cf[iv_col]),
           se       = sqrt(diag(vv)[iv_col]),
           n_obs    = nobs(m))
  })
}

# ============================================================
# Samples
# ============================================================

cat("\nFitting models across samples...\n")

# --- SYSTEM ---
sys_full  <- system_data
sys_no_ww <- system_data |> filter(!year %in% WAR_YEARS)

cat(sprintf("  System  full:         %d years\n", nrow(sys_full)))
cat(sprintf("  System  no WWI/WWII:  %d years (-%d)\n",
            nrow(sys_no_ww), nrow(sys_full) - nrow(sys_no_ww)))

# Baseline + drop-war rows
sys_rows <- bind_rows(
  fit_sys(sys_full,  "(1) Full sample"),
  fit_sys(sys_no_ww, "(2) Drop WWI + WWII")
)

# Sweep the large-conflict threshold
for (thr in SYSTEM_FATAL_MINS) {
  # Threshold applied to the contemporaneous fatal MID count, but IVs
  # remain lagged (subsampling the system-year panel still retains the
  # main-text temporal structure).
  d <- system_data |> filter(fatal_mid_count >= thr)
  cat(sprintf("  System  large confl:  %d years (>= %d fatal MIDs)\n",
              nrow(d), thr))
  sys_rows <- bind_rows(
    sys_rows,
    fit_sys(d, sprintf("(3) Fatal MIDs >= %d", thr))
  )
}

# --- STATE ---
st_full  <- state_churn
st_no_ww <- state_churn |> filter(!year %in% WAR_YEARS)
st_big   <- state_churn |> filter(mid_onset_count > 0)  # any MID participant

cat(sprintf("  State   full:         %d state-years\n", nrow(st_full)))
cat(sprintf("  State   no WWI/WWII:  %d state-years (-%d)\n",
            nrow(st_no_ww), nrow(st_full) - nrow(st_no_ww)))
cat(sprintf("  State   large confl:  %d state-years (fatal participants only)\n",
            nrow(st_big)))

st_rows <- bind_rows(
  fit_state(st_full,  "(1) Full sample"),
  fit_state(st_no_ww, "(2) Drop WWI + WWII"),
  fit_state(st_big,   "(3) MID participants only")
)

results <- bind_rows(sys_rows, st_rows) |>
  mutate(ci_low  = estimate - 1.96 * se,
         ci_high = estimate + 1.96 * se)

print(results, digits = 4)

write_csv(results,
          file.path(OUTPUT_DIR, "robustness_war_samples.csv"))
cat(sprintf("\nSaved: robustness_war_samples.csv (%d rows)\n", nrow(results)))

# ============================================================
# Plot
# ============================================================

cat("\n--- Plotting ---\n")

plot_df <- results |>
  mutate(level_label = factor(level, levels = c("system", "state"),
                               labels = c("System (JS divergence)",
                                           "State (churn)")),
         sample = factor(sample, levels = unique(sample)),
         iv_label = factor(iv_label,
                           levels = c(sapply(sys_ivs, `[[`, "label"),
                                       sapply(state_ivs, `[[`, "label"))))

p16_system <- ggplot(subset(plot_df, level_label == "System (JS divergence)" &
                              sample != "(2) Drop WWI + WWII"),
              aes(x = estimate, y = forcats::fct_rev(sample))) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray55") +
  geom_errorbarh(aes(xmin = ci_low, xmax = ci_high),
                 height = 0.15, linewidth = 0.7, color = ACCENT) +
  geom_point(size = 2.6, color = ACCENT) +
  facet_wrap( ~ iv_label, scales = "free_x",
              labeller = label_wrap_gen(width = 22)) +
  labs(
    title = "System: sample restrictions",
    x = "Coefficient (95% CI)", y = NULL
  ) +
  theme_pres() +
  theme(panel.spacing = unit(1.2, "lines"),
        panel.spacing.x = unit(1.5, "lines"),
        strip.text = element_text(color = "black", size = 9,
                                    face = "bold",
                                    margin = margin(4, 4, 4, 4)),
        strip.placement = "outside",
        plot.margin = margin(t = 10, r = 12, b = 14, l = 10))

save_pres("fig_ms_war_samples_system.png", p16_system,
          width = 9, height = 5.5)


p16_state <- ggplot(subset(plot_df, level_label != "System (JS divergence)" &
                             sample == "(1) Full sample"),
                     aes(x = estimate, y = forcats::fct_rev(sample))) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray55") +
  geom_errorbarh(aes(xmin = ci_low, xmax = ci_high),
                 height = 0.15, linewidth = 0.7, color = ACCENT) +
  geom_point(size = 2.6, color = ACCENT) +
  facet_wrap( ~ iv_label, scales = "free_x",
              labeller = label_wrap_gen(width = 22)) +
  labs(
    title = "State: sample restrictions",
    x = "Coefficient (95% CI)", y = NULL
  ) +
  theme_pres() +
  theme(panel.spacing = unit(1.2, "lines"),
        panel.spacing.x = unit(1.5, "lines"),
        strip.text = element_text(color = "black", size = 9,
                                  face = "bold",
                                  margin = margin(4, 4, 4, 4)),
        strip.placement = "outside",
        plot.margin = margin(t = 10, r = 12, b = 14, l = 10))

save_pres("fig_ms_war_samples_state.png", p16_state,
          width = 9, height = 5.5)

cat("\nDone.\n")
