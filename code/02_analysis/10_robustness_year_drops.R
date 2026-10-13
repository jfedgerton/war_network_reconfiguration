# [2026-09-26] Headline system-level DV switched from KL to Jensen-Shannon (JS)
# divergence. KL results are retained in results/appendix/kl_canonical/.
###############################################################################
# 10_robustness_year_drops.R
# Robustness: drop years from the OLS and re-estimate, at BOTH levels.
#
# SYSTEM LEVEL
#   DV: js_divergence
#   IVs: mid_count, fatal_mid_count at t-1 (raw and log) [2026-10-12]
#   Controls: log1p(n_democracies), log1p(n_autocracies)
#   SEs: Newey-West HAC (lag = 5)
#
# STATE LEVEL
#   DV: churn
#   IVs: mid_onset_count, fatal_mid_onset_count at t-1 (raw and log) [2026-10-12]
#   Controls: democracy, autocracy, cinc, major_power, is_isolated_t1
#   FE: ccode + year  (TWFE)
#   SEs: two-way clustered (~ ccode + year)
#
# Two variants at each level:
#   (A) RANDOM year drops of size 10, 15, 20, 25 — 500 reps each
#   (B) SLIDING-WINDOW exclusion of 10, 15, 20, 25 consecutive years
#
# Outputs:
#   data/processed/robustness_year_drops_baseline.csv
#   data/processed/robustness_year_drops_random.csv
#   data/processed/robustness_year_drops_sliding.csv
#   code/03_figures/output/presentation/
#       pres_fig12_random_year_drops.png    (system + state, faceted)
#       pres_fig13_sliding_window.png       (system + state, faceted)
###############################################################################

library(tidyverse)
library(fixest)
library(patchwork)
library(scales)
library(peacesciencer)
select <- dplyr::select

set.seed(123)

OUTPUT_DIR <- here::here("data", "processed")
FIG_DIR    <- here::here("code", "03_figures", "output", "manuscript")
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# Manuscript theme — white background, publication-suitable
# ============================================================

ACCENT  <- "#2C6CB0"   # muted blue (darker than slide blue for white bg)
ACCENT2 <- "#B03366"   # muted red/magenta for contrast on white

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
# Load data
# ============================================================

cat("Loading data...\n")

# System-level. Match main-text baseline in 01_system_level_analysis.R:
# KL_t = alpha + beta*IV_{t-1} + log1p(n_democracies) + log1p(n_autocracies),
# plain OLS (no HAC), lagged IVs, no post_1945 indicator.
system_data <- readRDS(file.path(OUTPUT_DIR, "system_analysis_data.rds")) |>
  arrange(year) |>
  # [2026-10-12] conflict lags are built BEFORE dropping years without JS, so 1817 keeps its 1816 conflict
  (\(d) { stopifnot(all(diff(d$year) == 1)); d })() |>   # consecutive years
  mutate(
    log_js_divergence    = log1p(js_divergence),
    log_mid_count        = log1p(mid_count),
    log_fatal_mid_count  = log1p(fatal_mid_count),
    mid_count_lag1           = dplyr::lag(mid_count,           1),
    log_mid_count_lag1       = dplyr::lag(log_mid_count,       1),
    fatal_mid_count_lag1     = dplyr::lag(fatal_mid_count,     1),
    log_fatal_mid_count_lag1 = dplyr::lag(log_fatal_mid_count, 1)
  ) |>
  filter(!is.na(mid_count_lag1), !is.na(js_divergence))

sys_years <- sort(unique(system_data$year))
cat(sprintf("  System: %d years (%d-%d)\n",
            length(sys_years), min(sys_years), max(sys_years)))

# State-level
state_churn <- readRDS(file.path(OUTPUT_DIR, "state_network_churn.rds"))
dyad_data   <- readRDS(file.path(OUTPUT_DIR, "dyad_year_network_layers.rds"))

needs_prep <- !all(c("cinc", "major_power", "democracy", "autocracy",
                      "partner_war_frac") %in% names(state_churn))

if (needs_prep) {
  cat("  Preparing state-level controls (cinc, democracy, partner_war_frac...)...\n")

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

# --- State-level control set ---
# Non-fatal and fatal IVs share the same state-level controls:
# democracy + autocracy + cinc + major_power + is_isolated_t1.
# [2026-10-12] partner war exposure removed from the control set (H3 dropped from the paper).
CTRL_MID   <- "democracy + autocracy + cinc + major_power + is_isolated_t1"
CTRL_FATAL <- "democracy + autocracy + cinc + major_power + is_isolated_t1"
ctrl_for_iv <- function(iv) if (grepl("fatal", iv)) CTRL_FATAL else CTRL_MID

state_churn <- state_churn |>
  mutate(log_churn                 = log1p(churn),
         log_mid_onset_count       = log1p(mid_onset_count),
         log_fatal_mid_onset_count = log1p(fatal_mid_onset_count))
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
chk_univ <- dplyr::inner_join(state_churn |> dplyr::select(ccode, year, mid_onset_count, fatal_mid_onset_count),
                              onset_univ, by = c("ccode", "year"))
stopifnot(nrow(chk_univ) == nrow(state_churn), all(chk_univ$mid_onset_count == chk_univ$u_onset),
          all(chk_univ$fatal_mid_onset_count == chk_univ$u_fatal_onset))
rm(chk_univ)
state_churn <- state_churn |>
  dplyr::select(-dplyr::any_of(c("onset_lag1", "fatal_onset_lag1", "log_onset_lag1", "log_fatal_onset_lag1"))) |>
  dplyr::left_join(onset_univ |>
                     dplyr::transmute(ccode, year = year + 1L,
                                      onset_lag1 = u_onset, fatal_onset_lag1 = u_fatal_onset,
                                      log_onset_lag1 = log1p(u_onset), log_fatal_onset_lag1 = log1p(u_fatal_onset)),
                   by = c("ccode", "year"))
stopifnot(!anyDuplicated(state_churn[c("ccode", "year")]))

st_years <- sort(unique(state_churn$year))
cat(sprintf("  State: %d years, %d state-year obs\n",
            length(st_years), nrow(state_churn)))

# ============================================================
# Model specs + fit helpers — four IVs, each fit separately
# ============================================================

# System-level IVs: lagged (t-1) to match 01_system_level_analysis.R
sys_ivs <- list(
  list(col = "mid_count_lag1",           label = "All MIDs"),
  list(col = "log_mid_count_lag1",       label = "log All MIDs"),
  list(col = "fatal_mid_count_lag1",     label = "Fatal MIDs"),
  list(col = "log_fatal_mid_count_lag1", label = "log Fatal MIDs")
)

# [2026-10-12] state onset IVs at t-1 (headline timing); labels shared with the system IVs.
state_ivs <- list(
  list(col = "onset_lag1",           label = "All MIDs"),
  list(col = "log_onset_lag1",       label = "log All MIDs"),
  list(col = "fatal_onset_lag1",     label = "Fatal MIDs"),
  list(col = "log_fatal_onset_lag1", label = "log Fatal MIDs")
)

sys_controls_str <- "log1p(n_democracies) + log1p(n_autocracies) + log1p(n_states) + poly(year, 2)"
# state_controls is now routed per IV via ctrl_for_iv() inside fit_state().

# fit_sys/fit_state now return 4 rows (one per IV) per call.
# System level uses plain OLS (lm) to match main-text; state uses feols TWFE.
fit_sys <- function(d) {
  purrr::map_dfr(sys_ivs, function(iv) {
    f <- as.formula(paste0("js_divergence ~ ", iv$col, " + ", sys_controls_str))
    m <- tryCatch(lm(f, data = d), error = function(e) NULL)
    if (is.null(m)) return(NULL)
    cf <- coef(m); vv <- vcov(m)
    tibble(level = "system",
           term = iv$col, iv_label = iv$label,
           estimate = unname(cf[iv$col]),
           se       = sqrt(diag(vv)[iv$col]),
           n_obs    = nobs(m))
  })
}

fit_state <- function(d) {
  purrr::map_dfr(state_ivs, function(iv) {
    f <- as.formula(paste0("churn ~ ", iv$col, " + ", ctrl_for_iv(iv$col),
                            " | ccode + year"))
    m <- tryCatch(feols(f, data = d, cluster = ~ ccode + year),
                   error = function(e) NULL)
    if (is.null(m)) return(NULL)
    cf <- coef(m); vv <- vcov(m)
    tibble(level = "state",
           term = iv$col, iv_label = iv$label,
           estimate = unname(cf[iv$col]),
           se       = sqrt(diag(vv)[iv$col]),
           n_obs    = nobs(m))
  })
}

# ============================================================
# Baselines
# ============================================================

cat("\n--- Baselines (full sample) ---\n")
baseline <- bind_rows(fit_sys(system_data),
                      fit_state(state_churn)) |>
  mutate(ci_low  = estimate - 1.96 * se,
         ci_high = estimate + 1.96 * se)
print(baseline, digits = 4)

write_csv(baseline,
          file.path(OUTPUT_DIR, "robustness_year_drops_baseline.csv"))

# ============================================================
# (A) RANDOM DROPS
# ============================================================

drop_sizes <- c(5, 10, 15, 20, 25)
n_reps     <- 500

cat(sprintf("\n--- Random drops (%d reps, drop sizes %s) ---\n",
            n_reps, paste(drop_sizes, collapse = ", ")))

random_rows <- list(); idx <- 1L

# -- system
for (ds in drop_sizes) {
  for (r in seq_len(n_reps)) {
    dropped <- sample(sys_years, ds, replace = FALSE)
    d <- system_data |> filter(!year %in% dropped)
    res <- fit_sys(d); if (is.null(res)) next
    res$drop_size <- ds; res$rep <- r
    for (i in seq_len(nrow(res))) {
      random_rows[[idx]] <- res[i, ]; idx <- idx + 1L
    }
  }
  cat(sprintf("  system   drop %2d: %d reps done\n", ds, n_reps))
}

# -- state
for (ds in drop_sizes) {
  for (r in seq_len(n_reps)) {
    dropped <- sample(st_years, ds, replace = FALSE)
    d <- state_churn |> filter(!year %in% dropped)
    res <- fit_state(d); if (is.null(res)) next
    res$drop_size <- ds; res$rep <- r
    for (i in seq_len(nrow(res))) {
      random_rows[[idx]] <- res[i, ]; idx <- idx + 1L
    }
  }
  cat(sprintf("  state    drop %2d: %d reps done\n", ds, n_reps))
}

random_df <- bind_rows(random_rows) |>
  mutate(ci_low = estimate - 1.96 * se,
         ci_high = estimate + 1.96 * se)

write_csv(random_df,
          file.path(OUTPUT_DIR, "robustness_year_drops_random.csv"))
cat(sprintf("  Saved random drops: %d rows\n", nrow(random_df)))

# ============================================================
# (B) SLIDING WINDOWS
# ============================================================

cat("\n--- Sliding-window exclusion ---\n")

sliding_rows <- list(); idx <- 1L

# -- system
for (ws in drop_sizes) {
  starts <- seq(min(sys_years), max(sys_years) - ws + 1)
  for (start in starts) {
    excluded <- start:(start + ws - 1)
    d <- system_data |> filter(!year %in% excluded)
    res <- fit_sys(d); if (is.null(res)) next
    res$window_size <- ws; res$start_year <- start
    res$end_year <- start + ws - 1
    res$mid_year <- start + (ws - 1) / 2
    for (i in seq_len(nrow(res))) {
      sliding_rows[[idx]] <- res[i, ]; idx <- idx + 1L
    }
  }
  cat(sprintf("  system   window %2d: %d positions\n", ws, length(starts)))
}

# -- state
for (ws in drop_sizes) {
  starts <- seq(min(st_years), max(st_years) - ws + 1)
  for (start in starts) {
    excluded <- start:(start + ws - 1)
    d <- state_churn |> filter(!year %in% excluded)
    res <- fit_state(d); if (is.null(res)) next
    res$window_size <- ws; res$start_year <- start
    res$end_year <- start + ws - 1
    res$mid_year <- start + (ws - 1) / 2
    for (i in seq_len(nrow(res))) {
      sliding_rows[[idx]] <- res[i, ]; idx <- idx + 1L
    }
  }
  cat(sprintf("  state    window %2d: %d positions\n", ws, length(starts)))
}

sliding_df <- bind_rows(sliding_rows) |>
  mutate(ci_low = estimate - qnorm(0.95) * se)

write_csv(sliding_df,
          file.path(OUTPUT_DIR, "robustness_year_drops_sliding.csv"))
cat(sprintf("  Saved sliding window: %d rows\n", nrow(sliding_df)))

# ============================================================
# PLOT helpers
# ============================================================

iv_levels_all <- unique(c(sapply(sys_ivs, `[[`, "label"),
                          sapply(state_ivs, `[[`, "label")))

random_df <- random_df |>
  mutate(iv_label = factor(iv_label, levels = iv_levels_all),
         drop_size_label = factor(paste0("Drop ", drop_size, " years"),
                                   levels = paste0("Drop ", drop_sizes,
                                                    " years")),
         level_label = factor(level, levels = c("system", "state"),
                               labels = c("System (JS)",
                                           "State (churn)")))

sliding_df <- sliding_df |>
  mutate(iv_label = factor(iv_label, levels = iv_levels_all),
         window_label = factor(paste0("Window = ", window_size, " yrs"),
                                levels = paste0("Window = ", drop_sizes,
                                                 " yrs")),
         level_label = factor(level, levels = c("system", "state"),
                               labels = c("System (JS)",
                                           "State (churn)")))

baseline_plot <- baseline |>
  mutate(iv_label = factor(iv_label, levels = iv_levels_all),
         level_label = factor(level, levels = c("system", "state"),
                                labels = c("System (JS)", "State (churn)")))

# ============================================================
# FIGURE 12: Random drops — one row, two panels per level:
#   (a) box plot of the coefficient across the 500 reps, and
#   (b) stacked bar of the share of reps significant at p < 0.05,
# with the number of years dropped on the x-axis, faceted by IV.
# ============================================================

cat("\n--- Plotting ---\n")

# Restrict the random-drop figures to the log(1 + fatal MID) IV only (fatal t-1 headline).
# Category = positive AND two-sided significant at 0.05 (blue) vs null (red),
# with a clean numeric x-axis label for the number of years dropped.
random_df <- random_df |>
  filter(iv_label == "log Fatal MIDs") |>
  mutate(
    drop_n  = factor(drop_size, levels = drop_sizes),
    sig_pos = (abs(estimate / se) > qnorm(0.975)) & (estimate > 0),
    sig_cat = factor(if_else(sig_pos,
                             "Positive & significant (p < 0.05)", "Null"),
                     levels = c("Positive & significant (p < 0.05)", "Null"))
  )

make_drop_panel <- function(lvl_label, file_out) {
  d_lvl <- random_df |> filter(level_label == lvl_label)

  # (a) coefficient distribution across reps
  p_box <- ggplot(d_lvl, aes(x = drop_n, y = estimate)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray60") +
    geom_boxplot(fill = ACCENT, alpha = 0.7, width = 0.6,
                 outlier.size = 0.4, linewidth = 0.35) +
    facet_wrap(~ iv_label, scales = "free_y", ncol = 1,
               labeller = label_wrap_gen(width = 22)) +
    labs(x = "Years dropped (random)", y = "Coefficient",
         title = "(a) Coefficient distribution") +
    theme_pres(base_size = 11) +
    theme(strip.text = element_text(size = 8, face = "bold"))

  # (b) share of reps significant at p < 0.05
  share_df <- d_lvl |>
    count(iv_label, drop_n, sig_cat, .drop = FALSE) |>
    mutate(iv_label = as.character(iv_label)) |>
    filter(iv_label == "log Fatal MIDs") |> 
    #mutate(iv_label = ifelse(iv_label == "log Fatal MIDs", "Log count of MIDs"))
    group_by(iv_label, drop_n) |>
    mutate(share = n / sum(n)) |>
    ungroup()

   p_bar <- ggplot(share_df, aes(x = drop_n, y = share, fill = sig_cat)) +
     geom_col(width = 0.6) +
     geom_text(aes(label = ifelse(sig_cat == "Positive & significant (p < 0.05)" & share > 0.06,
                                   scales::percent(share, accuracy = 1), "")),
               position = position_stack(vjust = 0.5),
               color = "white", size = 2.8, fontface = "bold") +
     facet_wrap(~ iv_label, ncol = 1, labeller = label_wrap_gen(width = 22)) +
     scale_fill_manual(values = c("Positive & significant (p < 0.05)" = "#2E75B6",
                                   "Null"                             = "#C0392B"),
                       name = NULL) +
     scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                        expand = expansion(mult = c(0, 0.02))) +
     labs(x = "Years dropped (random)", y = "Share of 500 reps",
          title = "(b) % positive & significant at p < 0.05") +
     theme_pres(base_size = 11) +
     theme(legend.position = "bottom",
           strip.text = element_text(size = 8, face = "bold"))
   p_combined <- p_box + p_bar +
     plot_layout(widths = c(1, 1)) +
     plot_annotation(
       title = sprintf("%s: random-year drops (500 reps per drop size)", lvl_label)
     )

  save_pres(file_out, p_combined, width = 11, height = 5)
}

make_drop_panel("System (JS)",  "fig_ms_random_year_drops_system.png")
make_drop_panel("State (churn)", "fig_ms_random_year_drops_state.png")

# ============================================================
# FIGURE 13: Sliding window (both levels, 4 IVs)
# ============================================================
# Color by window size, facet by level × IV.
# Baseline (full-sample) estimate as dashed reference line per panel.

p13system <- ggplot(subset(sliding_df, level_label == "System (JS)" &
                             iv_label == "log Fatal MIDs"),
              aes(x = mid_year, y = estimate,
                  color = window_label, fill = window_label)) +
  geom_hline(yintercept = 0, color = "gray55", linewidth = 0.3) +
  geom_ribbon(aes(ymin = ci_low, ymax = estimate),
              alpha = 0.12, color = NA) +
  geom_line(linewidth = 0.6) +
  scale_color_brewer(palette = "Dark2", name = NULL) +
  scale_fill_brewer(palette = "Dark2", name = NULL) +
  facet_wrap( ~ iv_label, scales = "free_y",
              labeller = label_wrap_gen(width = 22)) +
  labs(
    title = "System: sliding-window exclusion",
    x = "Midpoint year of excluded window",
    y = "Coefficient (95% CI One-Sided)"
  ) +
  theme_pres() +
  theme(panel.spacing = unit(1.2, "lines"),
        panel.spacing.x = unit(1.5, "lines"),
        strip.placement = "outside",
        legend.text = element_text(size = 9),
        strip.text = element_text(size = 9, face = "bold",
                                    margin = margin(4, 4, 4, 4)),
        plot.margin = margin(t = 10, r = 12, b = 14, l = 10))

save_pres("fig_ms_sliding_window_system.png", p13system,
          width = 8, height = 5)

p13sstate <- ggplot(subset(sliding_df, level_label != "System (JS)" &
                             iv_label == "log Fatal MIDs"),
                    aes(x = mid_year, y = estimate,
                        color = window_label, fill = window_label)) +
  geom_hline(yintercept = 0, color = "gray55", linewidth = 0.3) +
  geom_ribbon(aes(ymin = ci_low, ymax = estimate),
              alpha = 0.12, color = NA) +
  geom_line(linewidth = 0.6) +
  scale_color_brewer(palette = "Dark2", name = NULL) +
  scale_fill_brewer(palette = "Dark2", name = NULL) +
  facet_wrap( ~ iv_label, scales = "free_y",
              labeller = label_wrap_gen(width = 22)) +
  labs(
    title = "State: sliding-window exclusion",
    x = "Midpoint year of excluded window",
    y = "Coefficient (95% CI One-Sided)"
  ) +
  theme_pres() +
  theme(panel.spacing = unit(1.2, "lines"),
        panel.spacing.x = unit(1.5, "lines"),
        strip.placement = "outside",
        legend.text = element_text(size = 9),
        strip.text = element_text(size = 9, face = "bold",
                                  margin = margin(4, 4, 4, 4)),
        plot.margin = margin(t = 10, r = 12, b = 14, l = 10))

save_pres("fig_ms_sliding_window_state.png", p13sstate,
          width = 8, height = 5)



cat("\nDone.\n")