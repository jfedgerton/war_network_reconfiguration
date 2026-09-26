###############################################################################
# 26_multiverse_summary_plots.R
# Summary visualizations for the layer combo × weighting multiverse from
# scripts 24 (excl-belligerent) and 25 (strength control).
#
# Produces:
#   (A) Box plots of coefficient distributions (canonical vs alt)
#   (B) Stacked bar plots of p-value significance bins
#
# Inputs:
#   data/processed/robustness_excl_belligerent_multiverse.csv
#   data/processed/robustness_strength_control_multiverse.csv (if available)
# Outputs:
#   code/03_figures/output/manuscript/fig_ms_multiverse_coef_boxplot.png
#   code/03_figures/output/manuscript/fig_ms_multiverse_pval_bars.png
###############################################################################

library(tidyverse)

OUTPUT_DIR <- here::here("data", "processed")
FIG_DIR    <- here::here("code", "03_figures", "output", "manuscript")
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

ACCENT  <- "#2C6CB0"
ACCENT2 <- "#B03366"
NEUTRAL <- "#888888"

theme_pres <- function(base_size = 14) {
  theme_bw(base_size = base_size) %+replace%
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.major.y = element_line(color = "gray90", linewidth = 0.3),
      plot.title       = element_text(size = base_size + 2, face = "bold",
                                       hjust = 0, margin = margin(b = 10)),
      plot.subtitle    = element_text(size = base_size - 1, color = "gray30",
                                       margin = margin(b = 10)),
      axis.text        = element_text(color = "black", size = base_size - 1),
      axis.title       = element_text(color = "black", size = base_size),
      axis.ticks       = element_line(color = "gray40"),
      strip.background = element_rect(fill = "gray92", color = "gray70"),
      strip.text       = element_text(face = "bold", size = base_size, color = "black"),
      legend.position  = "bottom",
      legend.text      = element_text(size = base_size - 1),
      legend.title     = element_blank(),
      plot.margin      = margin(12, 14, 10, 12)
    )
}

save_fig <- function(filename, plot, width = 10, height = 7, dpi = 400) {
  ggsave(file.path(FIG_DIR, filename), plot,
         width = width, height = height, dpi = dpi, bg = "white")
  cat(sprintf("  -> %s\n", filename))
}

# Helper: bin p-values into significance categories
bin_pval <- function(p) {
  cut(p, breaks = c(0, 0.01, 0.05, 0.10, 1),
      labels = c("p < 0.01", "p < 0.05", "p < 0.10", "p >= 0.10"),
      include.lowest = TRUE, right = FALSE)
}

# ============================================================================
# 1. Load data
# ============================================================================

panels <- list()

# Script 24: excl-belligerent
f24 <- file.path(OUTPUT_DIR, "robustness_excl_belligerent_multiverse.csv")
if (file.exists(f24)) {
  d24 <- read_csv(f24, show_col_types = FALSE) |>
    mutate(panel = "Excl. belligerent dyads",
           group = outcome,
           pval_bin = bin_pval(pval))
  panels[["excl_bell"]] <- d24
  cat(sprintf("Loaded script 24 multiverse: %d rows\n", nrow(d24)))
}

# Script 25: strength control
f25 <- file.path(OUTPUT_DIR, "robustness_strength_control_multiverse.csv")
if (file.exists(f25)) {
  d25 <- read_csv(f25, show_col_types = FALSE) |>
    mutate(panel = "Strength control",
           group = spec,
           pval_bin = bin_pval(pval))
  # If 25 has both MID onset and fatal MID, keep only MID onset for cleaner plot
  if ("iv" %in% names(d25)) {
    d25 <- d25 |> filter(iv == "log_mid_onset_count")
  }
  panels[["strength"]] <- d25
  cat(sprintf("Loaded script 25 multiverse: %d rows\n", nrow(d25)))
}

if (length(panels) == 0) {
  stop("No multiverse CSV files found. Run scripts 24 and/or 25 first.")
}

# ============================================================================
# 2. Coefficient box plots
# ============================================================================

cat("\n=== Coefficient box plots ===\n")

plot_list_coef <- list()

if ("excl_bell" %in% names(panels)) {
  p24_box <- ggplot(panels$excl_bell,
                     aes(x = group, y = estimate, fill = group)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    geom_boxplot(width = 0.6, outlier.size = 1.2, outlier.alpha = 0.5) +
    scale_fill_manual(values = c("Canonical" = ACCENT,
                                  "Excl. belligerent" = ACCENT2)) +
    labs(x = NULL, y = "Coefficient on log(1+MID Onset)",
         title = "A. Excl. belligerent dyads",
         subtitle = sprintf("120 layer/weighting combos. Canonical median: %.4f | Excl. median: %.4f",
                            median(panels$excl_bell$estimate[panels$excl_bell$group == "Canonical"], na.rm = TRUE),
                            median(panels$excl_bell$estimate[panels$excl_bell$group == "Excl. belligerent"], na.rm = TRUE))) +
    theme_pres(base_size = 13) +
    theme(legend.position = "none")
  plot_list_coef[["excl_bell"]] <- p24_box
}

if ("strength" %in% names(panels)) {
  p25_box <- ggplot(panels$strength,
                     aes(x = group, y = estimate, fill = group)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    geom_boxplot(width = 0.6, outlier.size = 1.2, outlier.alpha = 0.5) +
    scale_fill_manual(values = c("Without strength" = ACCENT,
                                  "With strength" = ACCENT2)) +
    labs(x = NULL, y = "Coefficient on log(1+MID Onset)",
         title = "B. Strength control",
         subtitle = sprintf("120 layer/weighting combos. Without median: %.4f | With median: %.4f",
                            median(panels$strength$estimate[panels$strength$group == "Without strength"], na.rm = TRUE),
                            median(panels$strength$estimate[panels$strength$group == "With strength"], na.rm = TRUE))) +
    theme_pres(base_size = 13) +
    theme(legend.position = "none")
  plot_list_coef[["strength"]] <- p25_box
}

if (length(plot_list_coef) == 2) {
  library(patchwork)
  p_coef_combined <- plot_list_coef[[1]] + plot_list_coef[[2]] +
    plot_layout(ncol = 2)
  save_fig("fig_ms_multiverse_coef_boxplot.png", p_coef_combined,
           width = 12, height = 6)
} else if (length(plot_list_coef) == 1) {
  save_fig("fig_ms_multiverse_coef_boxplot.png", plot_list_coef[[1]],
           width = 7, height = 6)
}

# ============================================================================
# 3. P-value stacked bar plots
# ============================================================================

cat("\n=== P-value stacked bar plots ===\n")

# Colors for significance bins (green = strong, red = not sig)
pval_colors <- c("p < 0.01" = "#1a9641",
                 "p < 0.05" = "#a6d96a",
                 "p < 0.10" = "#fdae61",
                 "p >= 0.10" = "#d7191c")

plot_list_pval <- list()

if ("excl_bell" %in% names(panels)) {
  # Only positive coefficients count as "significant in the right direction"
  d24_pval <- panels$excl_bell |>
    mutate(pval_bin_dir = if_else(estimate > 0, as.character(pval_bin), "p >= 0.10"),
           pval_bin_dir = factor(pval_bin_dir,
                                 levels = c("p < 0.01", "p < 0.05", "p < 0.10", "p >= 0.10")))

  pval_summary_24 <- d24_pval |>
    count(group, pval_bin_dir) |>
    group_by(group) |>
    mutate(pct = 100 * n / sum(n)) |>
    ungroup()

  p24_pval <- ggplot(pval_summary_24,
                      aes(x = group, y = pct, fill = pval_bin_dir)) +
    geom_col(width = 0.6, color = "white", linewidth = 0.3) +
    scale_fill_manual(values = pval_colors) +
    scale_y_continuous(labels = function(x) paste0(x, "%")) +
    labs(x = NULL, y = "% of specifications",
         title = "A. Excl. belligerent dyads",
         subtitle = "Significance of positive MID-onset coefficients across 120 combos") +
    theme_pres(base_size = 13) +
    theme(legend.position = "bottom")
  plot_list_pval[["excl_bell"]] <- p24_pval
}

if ("strength" %in% names(panels)) {
  d25_pval <- panels$strength |>
    mutate(pval_bin_dir = if_else(estimate > 0, as.character(pval_bin), "p >= 0.10"),
           pval_bin_dir = factor(pval_bin_dir,
                                 levels = c("p < 0.01", "p < 0.05", "p < 0.10", "p >= 0.10")))

  pval_summary_25 <- d25_pval |>
    count(group, pval_bin_dir) |>
    group_by(group) |>
    mutate(pct = 100 * n / sum(n)) |>
    ungroup()

  p25_pval <- ggplot(pval_summary_25,
                      aes(x = group, y = pct, fill = pval_bin_dir)) +
    geom_col(width = 0.6, color = "white", linewidth = 0.3) +
    scale_fill_manual(values = pval_colors) +
    scale_y_continuous(labels = function(x) paste0(x, "%")) +
    labs(x = NULL, y = "% of specifications",
         title = "B. Strength control",
         subtitle = "Significance of positive MID-onset coefficients across 120 combos") +
    theme_pres(base_size = 13) +
    theme(legend.position = "bottom")
  plot_list_pval[["strength"]] <- p25_pval
}

if (length(plot_list_pval) == 2) {
  library(patchwork)
  p_pval_combined <- plot_list_pval[[1]] + plot_list_pval[[2]] +
    plot_layout(ncol = 2, guides = "collect") &
    theme(legend.position = "bottom")
  save_fig("fig_ms_multiverse_pval_bars.png", p_pval_combined,
           width = 12, height = 6)
} else if (length(plot_list_pval) == 1) {
  save_fig("fig_ms_multiverse_pval_bars.png", plot_list_pval[[1]],
           width = 7, height = 6)
}

cat("\nDone. 26_multiverse_summary_plots.R complete.\n")
