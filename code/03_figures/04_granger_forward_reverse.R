###############################################################################
# 04_granger_forward_reverse.R
# Rebuilds fig_granger_forward_vs_reverse.png from the saved multiverse
# specifications, without rerunning 07_robustness_multiverse.R (~70 min).
# The plotting block is copied from 07 (keep the two in sync).
# Input:  data/processed/robustness_specifications.csv (written by 07)
# Output: code/03_figures/output/fig_granger_forward_vs_reverse.png
# [2026-11-05]
###############################################################################
suppressPackageStartupMessages({library(dplyr); library(readr); library(ggplot2)})
setwd(here::here())
FIG_DIR <- here::here("code", "03_figures", "output")
specs_df <- read_csv(here::here("data", "processed", "robustness_specifications.csv"), show_col_types = FALSE)
stopifnot(all(c("model", "dv", "war_variable", "pval") %in% names(specs_df)))
stopifnot(all(c("VAR_Granger", "VAR_Granger_reverse") %in% specs_df$model))

granger_compare <- specs_df |>
  filter(model %in% c("VAR_Granger", "VAR_Granger_reverse")) |>
  mutate(
    direction = ifelse(model == "VAR_Granger", "War → Network", "Network → War"),
    sig = pval < 0.05
  )

if (nrow(granger_compare) > 0) {
  granger_wv_labels <- c(
    mid_count = "MID\nCount",
    log_mid_count = "Log(1+MID\nCount)",
    fatal_mid_count = "Fatal MID\nCount",
    log_fatal_mid_count = "Log(1+Fatal\nMID Count)",
    mid_onset_count = "MID Onset\nCount",
    log_mid_onset_count = "Log(1+MID\nOnset Count)",
    fatal_mid_onset_count = "Fatal MID\nOnset Count",
    log_fatal_mid_onset_count = "Log(1+Fatal\nMID Onset Count)"
  )
  granger_dv_labels <- c(
    kl_divergence = "KL\nDivergence",
    js_divergence = "JS\nDivergence",
    edge_turnover_pct = "Edge\nTurnover (%)",
    degree_kl = "Degree\nKL",
    degree_wasserstein = "Degree\nWasserstein",
    log_kl_divergence = "Log KL\nDivergence",
    log_js_divergence = "Log JS\nDivergence",
    log_degree_kl = "Log Degree\nKL",
    log_degree_wasserstein = "Log Degree\nWasserstein"
  )

  granger_summary <- granger_compare |>
    group_by(direction, dv, war_variable) |>
    summarise(
      pct_sig = mean(sig) * 100,
      median_pval = median(pval, na.rm = TRUE),
      n = n(), .groups = "drop"
    ) |>
    mutate(
      wv_label = recode(war_variable, !!!granger_wv_labels),
      dv_label = recode(dv, !!!granger_dv_labels),
      wv_label = factor(wv_label, levels = c(
        "MID\nCount", "Log(1+MID\nCount)",
        "Fatal MID\nCount", "Log(1+Fatal\nMID Count)",
        "MID Onset\nCount", "Log(1+MID\nOnset Count)",
        "Fatal MID\nOnset Count", "Log(1+Fatal\nMID Onset Count)"
      ))
    )

  p_granger <- ggplot(granger_summary, aes(x = dv_label, y = pct_sig, fill = direction)) +
    geom_col(position = position_dodge(width = 0.7), width = 0.6) +
    facet_wrap(~ wv_label, ncol = 2,
               labeller = label_wrap_gen(width = 22)) +
    scale_fill_manual(values = c("War \u2192 Network" = "steelblue",
                                 "Network \u2192 War" = "firebrick")) +
    geom_hline(yintercept = 5, linetype = "dashed", color = "gray40") +
    labs(x = NULL, y = "% Significant (p < 0.05)",
         title = "Forward vs. reverse Granger",
         fill = "Direction") +
    theme_minimal(base_size = 13) +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5, size = 9, colour = "black"),   # horizontal two-line labels [2026-11-05]
          strip.text = element_text(face = "bold", size = 9),
          legend.position = "top",
          legend.title = element_text(size = 13, face = "bold"),
          legend.text = element_text(size = 12),
          legend.key.size = unit(1, "cm"),
          panel.grid.minor = element_blank(),
          panel.spacing = unit(1.2, "lines"),
          plot.margin = margin(t = 10, r = 12, b = 14, l = 10))

  ggsave(file.path(FIG_DIR, "fig_granger_forward_vs_reverse.png"), p_granger,
         width = 14, height = 8, dpi = 300)
}


stopifnot(exists("granger_summary"), all(granger_summary$pct_sig >= 0 & granger_summary$pct_sig <= 100),
          !any(is.na(granger_summary$dv_label)), !any(grepl("_", granger_summary$dv_label)))
cat("fig_granger_forward_vs_reverse.png written\n")
