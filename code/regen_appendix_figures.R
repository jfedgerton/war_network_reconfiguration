###############################################################################
# code/regen_appendix_figures.R
#
# Intermediary script that:
#   (1) Regenerates Figures 1, 2, and 8 with the styling fixes applied.
#       Figures 2 and 8 are rebuilt from cached multiverse data so the
#       full 03_robustness.R does NOT need to be re-run.
#   (2) Re-runs short analysis scripts whose figures are missing:
#         - 14_igo_membership_sliding.R (upper-tail + two-way IGO panels)
#         - 15_igo_membership_status.R  (IGO membership-status moderator)
#   (3) Mirrors all figures into figs/ for upload.
#   (4) Mirrors every PNG/PDF figure from the scattered output dirs
#       into a single ottoman_byzantine/figs/ folder for one-shot upload
#       to Overleaf.
#
# Usage: from project root
#   Rscript code/regen_appendix_figures.R
#
# Set seed inside re-sourced scripts (they each set their own seed = 123).
###############################################################################

suppressPackageStartupMessages({
  library(here)
  library(tidyverse)
  library(patchwork)
})

set.seed(123)

OUTPUT_DIR  <- here("data", "processed")
FIG_DIR_03  <- here("code", "03_figures", "output")
FIG_DIR_MS  <- here("code", "03_figures", "output", "manuscript")
FIG_DIR_RES <- here("results", "main")
FIGS_OUT    <- here("figs")
dir.create(FIGS_OUT,    recursive = TRUE, showWarnings = FALSE)
dir.create(FIG_DIR_03,  recursive = TRUE, showWarnings = FALSE)
dir.create(FIG_DIR_MS,  recursive = TRUE, showWarnings = FALSE)

cat("=========================================================\n")
cat("regen_appendix_figures.R\n")
cat("=========================================================\n\n")

# ──────────────────────────────────────────────────────────────────────────────
# (1a) FIGURE 1: fig_edge_weight_validation.png
#      Two-line y-axis labels. Re-source 08 to regenerate.
# ──────────────────────────────────────────────────────────────────────────────
cat("(1a) Figure 1: edge-weight validation\n")
val_script <- here("code", "02_analysis", "08_validate_edge_weights.R")
if (file.exists(val_script)) {
  tryCatch(source(val_script, echo = FALSE),
           error = function(e)
             cat("    ERROR:", conditionMessage(e), "\n"))
  cat("    -> fig_edge_weight_validation.png\n")
} else {
  cat("    SKIP — 08_validate_edge_weights.R not found\n")
}

# ──────────────────────────────────────────────────────────────────────────────
# (1b) FIGURE 8: fig_meta_coefplot.png
#      Two-line y-axis labels + larger text + taller. Rebuild from
#      cached robustness_meta_coefs.csv without re-running the multiverse.
# ──────────────────────────────────────────────────────────────────────────────
cat("\n(1b) Figure 8: meta-coefplot\n")
meta_csv <- file.path(OUTPUT_DIR, "robustness_meta_coefs.csv")
if (file.exists(meta_csv)) {
  meta_coefs_all <- read_csv(meta_csv, show_col_types = FALSE)

  p_meta <- meta_coefs_all |>
    ggplot(aes(x = estimate, y = reorder(term_clean, estimate),
               color = significant)) +
    geom_point(size = 2.8) +
    geom_errorbarh(aes(xmin = ci_lower, xmax = ci_upper), height = 0.25) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    facet_wrap(~ model, scales = "free_x", ncol = 3,
               labeller = label_wrap_gen(width = 22)) +
    scale_color_manual(values = c("FALSE" = "gray50", "TRUE" = "black")) +
    scale_y_discrete(labels = function(x) str_wrap(x, width = 20)) +
    labs(x = "Coefficient", y = NULL,
         title = "Meta-analysis", color = "p < 0.05") +
    theme_minimal() +
    theme(panel.grid.minor = element_blank(),
          legend.position  = "bottom",
          axis.text.x      = element_text(size = 12),
          axis.text.y      = element_text(size = 12, lineheight = 0.9),
          axis.title.x     = element_text(size = 14),
          strip.text       = element_text(size = 11, face = "bold"),
          legend.text      = element_text(size = 12),
          legend.title     = element_text(size = 12),
          panel.spacing    = unit(1.4, "lines"),
          panel.spacing.x  = unit(1.8, "lines"),
          plot.margin      = margin(t = 10, r = 12, b = 14, l = 10))

  ggsave(file.path(FIG_DIR_03, "fig_meta_coefplot.png"), p_meta,
         width = 15, height = 14, dpi = 300)
  cat("    -> fig_meta_coefplot.png (regenerated from cached meta_coefs)\n")
} else {
  cat("    SKIP — robustness_meta_coefs.csv not found\n")
  cat("    Run 03_robustness.R first to generate meta_coefs_all\n")
}

# ──────────────────────────────────────────────────────────────────────────────
# (1c) FIGURE 2: fig_forest_dvs.png
#      Larger y-axis text + taller. Reconstruct forest_all from
#      cached robustness_specifications.csv without re-running the multiverse.
# ──────────────────────────────────────────────────────────────────────────────
cat("\n(1c) Figure 2: forest plot\n")
specs_csv <- file.path(OUTPUT_DIR, "robustness_specifications.csv")
if (file.exists(specs_csv)) {
  specs_df <- read_csv(specs_csv, show_col_types = FALSE)

  ALPHA_ONESIDE <- 0.10
  is_sig_pos <- function(coef, pval) coef > 0 & (pval / 2) < ALPHA_ONESIDE

  forest_all <- specs_df |>
    filter(model %in% c("OLS_lagged", "FE_full_ccode_year"),
           weighting == "unweighted", transform == "raw",
           layers == "defensive+igo") |>
    mutate(
      ci_lo = coef - 1.96 * se,
      ci_hi = coef + 1.96 * se,
      dv_label = case_when(
        dv == "kl_divergence"     ~ "KL Divergence",
        dv == "js_divergence"     ~ "JS Divergence",
        dv == "edge_turnover_pct" ~ "Edge Turnover %",
        dv == "degree_kl"         ~ "Degree KL",
        dv == "degree_wasserstein"~ "Degree Wasserstein",
        dv == "churn_alt"         ~ "Churn",
        TRUE ~ dv
      ),
      row_label = paste0(dv_label, "\n(",
                         ifelse(model == "OLS_lagged", "OLS", "FE"), ")"),
      level = str_to_title(level),
      sig   = is_sig_pos(coef, pval)
    )

  forest_theme <- theme_minimal(base_size = 13) +
    theme(legend.position    = "bottom",
          panel.grid.minor   = element_blank(),
          strip.text         = element_text(face = "bold", size = 12),
          plot.title         = element_text(size = 13, face = "bold"),
          axis.text.y        = element_text(size = 15),
          axis.text.x        = element_text(size = 11),
          axis.title.x       = element_text(size = 13))

  make_forest_panel <- function(data, wv, title_text) {
    d <- data |> filter(war_variable == wv)
    if (nrow(d) == 0) return(NULL)
    ggplot(d, aes(x = coef, y = reorder(row_label, coef))) +
      geom_point(aes(color = sig), size = 3) +
      geom_errorbarh(aes(xmin = ci_lo, xmax = ci_hi, color = sig),
                     height = 0.2) +
      geom_vline(xintercept = 0, linetype = "dashed") +
      scale_color_manual(values = c("FALSE" = "gray50",
                                     "TRUE" = "forestgreen"),
                         name = sprintf("Sig. positive (one-sided, alpha = %.2f)",
                                         ALPHA_ONESIDE)) +
      labs(x = "Coefficient (95% CI)", y = NULL, title = title_text) +
      forest_theme
  }

  p_sys_mid     <- make_forest_panel(forest_all |> filter(level == "System"),
                                      "mid_count",       "System: MID Count")
  p_state_mid   <- make_forest_panel(forest_all |> filter(level == "State"),
                                      "mid_count",       "State: MID Count")
  p_sys_fatal   <- make_forest_panel(forest_all |> filter(level == "System"),
                                      "fatal_mid_count", "System: Fatal MID Count")
  p_state_fatal <- make_forest_panel(forest_all |> filter(level == "State"),
                                      "fatal_mid_count", "State: Fatal MID Count")

  `%||%` <- function(a, b) if (!is.null(a)) a else b

  if (!is.null(p_sys_mid)) {
    p_forest <- (p_sys_mid | (p_state_mid %||% plot_spacer())) /
                (p_sys_fatal | (p_state_fatal %||% plot_spacer())) +
      plot_layout(guides = "collect") +
      plot_annotation(
        title = "Effect sizes across DVs",
        theme = theme(plot.title = element_text(face = "bold", size = 16))
      ) &
      theme(legend.position = "bottom")

    ggsave(file.path(FIG_DIR_03, "fig_forest_dvs.png"), p_forest,
           width = 14, height = 14, dpi = 300)
    cat("    -> fig_forest_dvs.png (regenerated from cached specs_df)\n")
  } else {
    cat("    WARN — no system-MID rows in specs_df\n")
  }
} else {
  cat("    SKIP — robustness_specifications.csv not found\n")
  cat("    Run 03_robustness.R first to generate specs_df\n")
}

# ──────────────────────────────────────────────────────────────────────────────
# (2a) Re-run 14_igo_membership_sliding.R
#      Produces fig_ms_igo_share_lower_*, fig_ms_igo_share_upper, _both
# ──────────────────────────────────────────────────────────────────────────────
cat("\n(2a) IGO sliding figures (upper-tail + two-way)\n")
igo_sliding <- here("code", "02_analysis", "14_igo_membership_sliding.R")
if (file.exists(igo_sliding)) {
  tryCatch(source(igo_sliding, echo = FALSE),
           error = function(e)
             cat("    ERROR:", conditionMessage(e), "\n"))
  cat("    -> fig_ms_igo_share_upper.png\n")
  cat("    -> fig_ms_igo_share_both.png\n")
} else {
  cat("    SKIP — 14_igo_membership_sliding.R not found\n")
}

# [REMOVED 2026-07-14] (2b) 15_igo_membership_status.R -- archived. The file was
# truncated and the robustness (alternative IGO membership codings) changed
# nothing. fig_ms_igo_status_kl.png is no longer regenerated.

# [REMOVED 2026-07-14] (3) ABM OFAT sensitivity (Python).
#   The ABM was cut from the manuscript and the appendix.
#   Code archived at archive/abm_removed_20260714/code/03_abm/.

# ──────────────────────────────────────────────────────────────────────────────
# (3) Mirror every PNG/PDF into ottoman_byzantine/figs/ for one-shot upload
# ──────────────────────────────────────────────────────────────────────────────
cat("\n(3) Mirroring figures into ottoman_byzantine/figs/\n")

source_dirs <- c(
  FIG_DIR_MS,
  FIG_DIR_03
  FIG_DIR_RES
)

total_copied <- 0L
for (sd in source_dirs) {
  if (!dir.exists(sd)) next
  files <- list.files(sd, pattern = "\\.(png|pdf)$", full.names = TRUE,
                       ignore.case = TRUE, recursive = FALSE)
  for (f in files) {
    dest <- file.path(FIGS_OUT, basename(f))
    if (file.copy(f, dest, overwrite = TRUE)) {
      total_copied <- total_copied + 1L
    }
  }
  cat(sprintf("    %s -> %d files\n", sd, length(files)))
}

cat(sprintf("\nTotal files mirrored: %d\n", total_copied))
cat(sprintf("Mirror location: %s\n", FIGS_OUT))
cat("Upload contents to Overleaf's figs/ folder.\n")

cat("\n=========================================================\n")
cat("Done.\n")
cat("=========================================================\n")
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       