###############################################################################
# 07c_multiverse_figures_bellig.R   [2026-11-10]
#
# Rebuilds the two multiverse appendix figures from SAVED 07 output with the belligerent-dyad rule as a decision.
# No model is re-estimated. Copied from the figure block at the end of 07_robustness_multiverse.R; changes are
# marked [07c].
#
# Changes from 07's figure block
#   * Specs come from robustness_specifications_bellig.csv.gz (all nine dyad rules: "none" = all dyads, plus eight
#     rules that drop belligerent dyads before the outcome is computed). The "none" subset must reproduce
#     robustness_specifications.csv exactly (checked).
#   * Heatmap: new row (k) "Dyads dropped" with the nine rules. Every other row pools over the nine rules.
#     Row (h) (wild bootstrap) exists for all dyads only, so it compares like-for-like on rule "none".
#     Rows (i)-(j) (timing grid) pool over the nine rules in the timing file.
#   * Overall bars: all dyads AND pooled over the nine rules, for each hypothesis.
#   * Coefficient plot: baseline = model 3 of the main tables (H1: Table 2 column 3, OLS with the lagged DV,
#     model "sys_T1_ar1"; H2: Table 4 column 3, state + year TWFE, "FE_full_ccode_year"), dyad rule "none",
#     JS divergence / churn, raw, unweighted, main layers, fatal MIDs. Every coefficient is a deviation from that
#     baseline, which sits at 0 (dashed line; the baseline is also drawn as a labelled point at 0).
#     The dyad rule enters as a factor in all meta-regressions, including the timing grid.
#   * Writes support shares by rule to data/processed/multiverse_support_by_bellig_rule.csv.
#
# Outputs (code/03_figures/output/): fig_multiverse_hypotheses_bellig.png, fig_multiverse_meta_by_hyp_bellig.png
# Run from the repository root:  Rscript code/02_analysis/07c_multiverse_figures_bellig.R
###############################################################################


library(tidyverse)
library(patchwork)
set.seed(123)

OUTPUT_DIR <- here::here("data", "processed")
FIG_DIR    <- here::here("code", "03_figures", "output")
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

# [07c] all nine dyad rules
specs <- readr::read_csv(file.path(OUTPUT_DIR, "robustness_specifications_bellig.csv.gz"),
                         show_col_types = FALSE) |>
  mutate(dv_logged = grepl("^log_", dv), iv_is_fatal = grepl("fatal", war_variable),
         iv_is_log = grepl("^log_", war_variable))   # same definitions as 07 (dv_logged, iv_is_fatal, iv_is_log)
BELLIG_LEVELS <- c("none", "ongoing_any_t1", "ongoing_any_t1t", "ongoing_fatal_t1", "ongoing_fatal_t1t",
                   "onset_any_t1", "onset_any_t1t", "onset_fatal_t1", "onset_fatal_t1t")
stopifnot(setequal(unique(specs$bellig_rule), BELLIG_LEVELS))

# Python wrote "True"/"False"; coerce the design booleans to logical to be safe.
bool_cols <- c("has_offensive", "has_defensive", "has_igo", "has_trade",
               "has_dca", "iv_is_fatal", "iv_is_log", "dv_logged")
specs <- specs |> mutate(across(any_of(bool_cols), ~ as.logical(.x)))

# ===========================================================================
# 1. Classify specs into hypotheses; flag one-sided positive significance
# ===========================================================================
specs <- specs |>
  filter(!str_starts(model, "VAR_")) |>           # drop VAR / Granger / FEVD cells
  filter(!is.na(coef), !is.na(pval)) |>
  # [2026-11-05] all directional specs, including the log-GDP rows (GDP crosses the full state grid)

  mutate(
    hypothesis = case_when(
      level == "system"                            ~ "H1",
      level == "state"                             ~ "H2"
    ),
    sig_pos = coef > 0 & (pval / 2) < 0.10          # one-sided positive, p < 0.10
  ) |>
  filter(!is.na(hypothesis))

HYP <- c("H1", "H2")
specs <- specs |> mutate(extra_control = dplyr::coalesce(extra_control, "none"))
stopifnot(all(specs$extra_control %in% c("none", "gdp")))

# [07c] CHECK: the all-dyads subset reproduces the main 07 file (same rows, coefficients and p-values).
CHK_KEYS <- c("level", "dv", "war_variable", "weighting", "transform", "layers", "model", "extra_control")
main_chk <- readr::read_csv(file.path(OUTPUT_DIR, "robustness_specifications.csv"), show_col_types = FALSE) |>
  mutate(extra_control = dplyr::coalesce(extra_control, "none")) |>
  filter(!str_starts(model, "VAR_"), !is.na(coef), !is.na(pval)) |>
  select(all_of(CHK_KEYS), coef_m = coef, pval_m = pval)
none_chk <- specs |> filter(bellig_rule == "none") |> select(all_of(CHK_KEYS), coef, pval)
stopifnot(!anyDuplicated(main_chk[, CHK_KEYS]), !anyDuplicated(none_chk[, CHK_KEYS]), nrow(main_chk) == nrow(none_chk))
chk <- inner_join(none_chk, main_chk, by = CHK_KEYS)
stopifnot(nrow(chk) == nrow(main_chk), max(abs(chk$coef - chk$coef_m)) < 1e-10, max(abs(chk$pval - chk$pval_m)) < 1e-10)
cat(sprintf("  [07c] check passed: rule 'none' reproduces robustness_specifications.csv (%d specs)\n", nrow(chk)))
# [2026-11-12 churn_alt] State outcome = churn recomputed on each layer combination and weighting (churn_alt,
# log_churn_alt). The fixed churn rows use 02's main network in every cell: they do not vary with the weighting and
# vary across layer sets only through the trade-years sample, so they are dropped here (after the check above).
specs <- specs |> filter(!(level == "state" & dv %in% c("churn", "log_churn")))
stopifnot(all(specs$dv[specs$level == "state"] %in% c("churn_alt", "log_churn_alt")),
          sum(specs$level == "state") == 95040)
# Every rule must carry the same grid as "none".
grid_n <- specs |> count(hypothesis, bellig_rule)
stopifnot(all(table(grid_n$hypothesis) == length(BELLIG_LEVELS)),
          all(tapply(grid_n$n, grid_n$hypothesis, function(x) length(unique(x))) == 1))
specs_base <- specs |> filter(extra_control == "none")
N_SPEC <- table(specs$hypothesis)                                   # [07c] pooled over the nine rules

# ===========================================================================
# 2. Readable labels for each decision value
# ===========================================================================
dv_measure <- c(
  churn = "Churn", log_churn = "Churn",
  churn_alt = "Churn", log_churn_alt = "Churn",   # [2026-11-12 churn_alt] the only state DV
  degree_kl = "Degree KL", log_degree_kl = "Degree KL",
  degree_wasserstein = "Degree Wasserstein", log_degree_wasserstein = "Degree Wasserstein",
  edge_turnover_pct = "Edge turnover", log_edge_turnover_pct = "Edge turnover",
  js_divergence = "JS divergence", log_js_divergence = "JS divergence",
  kl_divergence = "KL divergence", log_kl_divergence = "KL divergence"
)
weight_lab <- c(
  contig_euclidean = "Contiguity (Euclidean)", contig_manhattan = "Contiguity (Manhattan)",
  dist_euclidean = "Distance (Euclidean)", dist_manhattan = "Distance (Manhattan)",
  unweighted = "Unweighted"
)
transform_lab <- c(log1p = "log(1 + x)", raw = "Raw", zscore = "z-score")
# [07c] dyad rules: event (MID onset / ongoing MID), severity (any / fatal), window (t-1 / t-1 and t)
bellig_lab <- c(none = "None (all dyads)",
                ongoing_any_t1 = "Ongoing MID, t-1", ongoing_any_t1t = "Ongoing MID, t-1 and t",
                ongoing_fatal_t1 = "Ongoing fatal MID, t-1", ongoing_fatal_t1t = "Ongoing fatal MID, t-1 and t",
                onset_any_t1 = "MID onset, t-1", onset_any_t1t = "MID onset, t-1 and t",
                onset_fatal_t1 = "Fatal MID onset, t-1", onset_fatal_t1t = "Fatal MID onset, t-1 and t")
stopifnot(setequal(names(bellig_lab), BELLIG_LEVELS))

# Display order (top -> bottom) within each category.
ITEM_ORDER <- c(
  # (a) reconfiguration measure
  "Churn", "Degree KL", "Degree Wasserstein",
  "Edge turnover", "JS divergence", "KL divergence",
  # (b) network layer
  "Defensive alliances", "Offensive alliances", "Shared IGOs", "Trade",
  # (c) transform
  "log(1 + x)", "Raw", "z-score",
  # (d) conflict measure
  "All MIDs", "Fatal MIDs",
  # (e) edge weighting
  "Contiguity (Euclidean)", "Contiguity (Manhattan)",
  "Distance (Euclidean)", "Distance (Manhattan)", "Unweighted",
  # (f) extra control
  "No GDP", "Log GDP",
  # (g)-(j) [2026-11-08]
  "Outcome raw", "Outcome log(1 + x)",
  "Clustered SE (state + year)", "Wild cluster bootstrap",
  "Conflict in year t", "Conflict in year t-1", "Years t and t-1 (sum)",
  "State + year FE", "State + year + community FE",
  # (k) dyads dropped [07c]
  unname(bellig_lab)
)
# [07c] short two-line labels for the vertical strips
CAT_A <- "(a) Reconfig-\nuration"
CAT_B <- "(b)\nLayer"
CAT_C <- "(c)\nTransform"
CAT_D <- "(d)\nCon-\nflict"   # [07c] three lines so it fits two rows at size 11
CAT_E <- "(e)\nWeighting"
CAT_F <- "(f)\nControl"
CAT_G <- "(g)\nScale"                 # [2026-11-08]
CAT_H <- "(h)\nInfer-\nence"   # [07c] three lines so it fits two rows at size 11             # [2026-11-08]
CAT_I <- "(i)\nTiming"                # [2026-11-08] timing grid
CAT_J <- "(j)\nFE"                    # [2026-11-08] timing grid
CAT_K <- "(k)\nDyads dropped"                  # [07c] belligerent-dyad rule

pct <- function(x) 100 * mean(x)   # share TRUE over a non-empty logical vector

# ===========================================================================
# 3. Heatmap data: % significant-positive by (decision value, hypothesis)
# ===========================================================================
d_recon <- specs |>
  mutate(item = dv_measure[dv]) |>
  group_by(hypothesis, item) |>
  summarise(val = pct(sig_pos), .groups = "drop") |>
  mutate(category = CAT_A)

layer_map <- tibble(
  col  = c("has_offensive", "has_defensive", "has_igo", "has_trade"),
  item = c("Offensive alliances", "Defensive alliances", "Shared IGOs", "Trade")
)
d_layer <- purrr::pmap_dfr(layer_map, function(col, item) {
  specs |> filter(.data[[col]]) |>
    group_by(hypothesis) |>
    summarise(val = pct(sig_pos), .groups = "drop") |>
    mutate(item = item)
}) |> mutate(category = CAT_B)

d_trans <- specs |>
  mutate(item = transform_lab[transform]) |>
  group_by(hypothesis, item) |>
  summarise(val = pct(sig_pos), .groups = "drop") |>
  mutate(category = CAT_C)

d_conf <- specs |>
  mutate(item = if_else(iv_is_fatal, "Fatal MIDs", "All MIDs")) |>
  group_by(hypothesis, item) |>
  summarise(val = pct(sig_pos), .groups = "drop") |>
  mutate(category = CAT_D)

d_weight <- specs |>
  mutate(item = weight_lab[weighting]) |>
  group_by(hypothesis, item) |>
  summarise(val = pct(sig_pos), .groups = "drop") |>
  mutate(category = CAT_E)

# (f) like-for-like: only the model classes that carry a GDP variant (state level), so None vs GDP
# compares the same models. H1 has no GDP rows and shows dashes.
CELL_KEYS <- c("level", "dv", "war_variable", "weighting", "transform", "layers", "model")
gdp_cells <- specs |> filter(extra_control == "gdp") |> distinct(across(all_of(CELL_KEYS)))
d_ctrl <- specs |>
  semi_join(gdp_cells, by = CELL_KEYS) |>   # only cells estimated both with and without GDP
  mutate(item = if_else(extra_control == "gdp", "Log GDP", "No GDP")) |>
  group_by(hypothesis, item) |>
  summarise(val = pct(sig_pos), .groups = "drop") |>
  mutate(category = CAT_F)

# [2026-11-08] (g) outcome scale: the same specifications with the dependent variable raw or log(1 + x).
d_scale <- specs |>
  mutate(item = if_else(dv_logged, "Outcome log(1 + x)", "Outcome raw")) |>
  group_by(hypothesis, item) |>
  summarise(val = pct(sig_pos), .groups = "drop") |>
  mutate(category = CAT_G)

# [2026-11-08] (h) inference, like-for-like: state TWFE cells estimated with both the state + year clustered SE
# and the wild cluster bootstrap p (same coefficient, different p-value).
BOOT_KEYS <- c(setdiff(CELL_KEYS, "model"), "extra_control")
wb <- readr::read_csv(file.path(OUTPUT_DIR, "robustness_specifications_wildboot.csv"), show_col_types = FALSE) |>
  mutate(across(any_of(bool_cols), ~ as.logical(.x)),
         extra_control = dplyr::coalesce(extra_control, "none"),
         hypothesis = "H2", sig_pos = coef > 0 & (pval / 2) < 0.10,
         dv_logged = grepl("^log_", dv), iv_is_fatal = grepl("fatal", war_variable)) |>
  filter(!is.na(coef), !is.na(pval), dv %in% c("churn_alt", "log_churn_alt"))   # [2026-11-12 churn_alt]
stopifnot(nrow(wb) == 1920)   # [2026-11-12 churn_alt] appendix Table 1 row (h)
stopifnot(nrow(wb) > 0, all(wb$model == "FE_full_wildboot"), !any(specs$model == "FE_full_wildboot"))
boot_cells <- wb |> distinct(across(all_of(BOOT_KEYS)))
d_inf <- bind_rows(specs |> filter(model == "FE_full_ccode_year", bellig_rule == "none"), wb) |>   # [07c] bootstrap exists for all dyads only
  semi_join(boot_cells, by = BOOT_KEYS) |>
  mutate(item = if_else(model == "FE_full_wildboot", "Wild cluster bootstrap", "Clustered SE (state + year)")) |>
  group_by(hypothesis, item) |>
  summarise(val = pct(sig_pos), .groups = "drop") |>
  mutate(category = CAT_H)
stopifnot(all(d_inf$hypothesis == "H2"), nrow(d_inf) == 2)

# [2026-11-08] (i)-(j) state timing x fixed-effects grid (kept out of specs_df and the overall bars above).
# [07c] All nine dyad rules, no extra control; the IV enters at t, at t-1, or as the t + (t-1) sum from the joint model.
# [2026-11-09] The timing file holds two versions (community = "canonical": 02's community on the main network for
# every layer combination; "combo": community detected on each layer combination). MV_TFE_COMMUNITY picks one
# (default combo). A file written before 2026-11-09 has no community column and only the canonical version.
TFE_COMMUNITY <- Sys.getenv("MV_TFE_COMMUNITY", "combo")
stopifnot(TFE_COMMUNITY %in% c("combo", "canonical"))
tfe <- readr::read_csv(file.path(OUTPUT_DIR, "robustness_specifications_timing_fe.csv.gz"), show_col_types = FALSE)
if (!"community" %in% names(tfe)) {
  if (TFE_COMMUNITY == "combo") stop("timing file has no combo-community rows (written before 2026-11-09); rerun 07 or set MV_TFE_COMMUNITY=canonical")
  tfe$community <- "canonical"
}
tfe <- tfe |>
  filter(community == TFE_COMMUNITY, dplyr::coalesce(extra_control, "none") == "none",   # [07c] all nine dyad rules
         iv_term %in% c("t", "t1", "joint_sum"), !is.na(coef), !is.na(pval),
         dv %in% c("churn_alt", "log_churn_alt")) |>   # [2026-11-12 churn_alt]
  mutate(hypothesis = "H2", sig_pos = coef > 0 & (pval / 2) < 0.10)
TFE_COMM_LABEL <- if (TFE_COMMUNITY == "combo") "community detected on each layer combination" else "canonical community (main network)"
stopifnot(nrow(tfe) > 0, all(tfe$fe_structure %in% c("state_year", "state_year_comm")),
          setequal(unique(tfe$bellig_rule), BELLIG_LEVELS))   # [07c]
N_TFE <- nrow(tfe)
if (TFE_COMMUNITY == "combo") stopifnot(N_TFE == 207360)   # [2026-11-12 churn_alt] appendix Table 1 timing grid
d_timing <- tfe |>
  mutate(item = dplyr::recode(iv_term, t = "Conflict in year t", t1 = "Conflict in year t-1",
                              joint_sum = "Years t and t-1 (sum)")) |>
  group_by(hypothesis, item) |>
  summarise(val = pct(sig_pos), .groups = "drop") |>
  mutate(category = CAT_I)
d_fe <- tfe |>
  mutate(item = if_else(fe_structure == "state_year", "State + year FE", "State + year + community FE")) |>
  group_by(hypothesis, item) |>
  summarise(val = pct(sig_pos), .groups = "drop") |>
  mutate(category = CAT_J)

# [07c] (k) dyads dropped: share significant-positive under each rule (same grid under every rule)
d_bellig <- specs |>
  mutate(item = bellig_lab[bellig_rule]) |>
  group_by(hypothesis, item) |>
  summarise(val = pct(sig_pos), .groups = "drop") |>
  mutate(category = CAT_K)

heat <- bind_rows(d_recon, d_layer, d_trans, d_conf, d_weight, d_ctrl, d_scale, d_inf, d_timing, d_fe, d_bellig) |>
  # Fill the (item x hypothesis) grid so undefined combinations become NA (dash),
  # NOT a mislabeled numeric. This is what removes the stray "NA" cell.
  complete(nesting(category, item), hypothesis = HYP) |>
  mutate(
    category   = factor(category, levels = c(CAT_A, CAT_B, CAT_C, CAT_D, CAT_E, CAT_F, CAT_G, CAT_H, CAT_I, CAT_J, CAT_K)),
    item       = factor(item, levels = rev(ITEM_ORDER)),   # first listed on top
    hypothesis = factor(hypothesis, levels = HYP),
    label      = ifelse(is.na(val), "—", sprintf("%.0f", val))
  )

p_heat <- ggplot(heat, aes(x = hypothesis, y = item, fill = val)) +
  geom_tile(color = "white", linewidth = 0.6) +
  geom_text(aes(label = label, color = !is.na(val) & val > 60), size = 3.8, show.legend = FALSE) +
  scale_color_manual(values = c(`FALSE` = "black", `TRUE` = "white"), guide = "none") +
  facet_grid(category ~ ., scales = "free_y", space = "free_y", switch = "y") +
  # [2026-11-13] cividis (sequential, colorblind-safe) replaces red-white-blue diverging at 50:
  # scale_fill_gradient2(low = "#C0392B", mid = "white", high = "#2C6CB0",
  #                      midpoint = 50, limits = c(0, 100),
  #                      na.value = "grey90", name = "% significant positive") +
  scale_fill_viridis_c(option = "cividis", direction = -1, limits = c(0, 100),
                       na.value = "grey90", name = "% significant positive") +
  scale_x_discrete(position = "top") +
  labs(title = "Which specification choices matter?",
       subtitle = "Percent significant-positive, pooled over nine dyad rules. Dashes = not defined at that level.",   # [07c] shortened; row notes moved to LaTeX caption
       # old_subtitle = sprintf("Share of specifications significant-positive, by decision, pooled over nine dyad rules\n(H1: %s specs; H2: %s specs). Dashes = not defined at that level.\nRow (f) compares the same models with and without log GDP. Row (k) shows each dyad rule.",
       #                   format(N_SPEC[["H1"]], big.mark = ","), format(N_SPEC[["H2"]], big.mark = ",")) |>
       #  paste0(sprintf("\nRow (h) compares the same TWFE models (all dyads) with clustered and wild bootstrap p-values.\nRows (i)-(j) use the separate timing grid (%s specs, nine dyad rules; %s).",
       #                 format(N_TFE, big.mark = ","), TFE_COMM_LABEL)),
       x = NULL, y = NULL) +
  guides(fill = guide_colorbar(barwidth = 16, barheight = 0.8, title.position = "top")) +
  theme_minimal(base_size = 11) +
  theme(panel.grid = element_blank(),
        strip.text.y.left = element_text(angle = 90, hjust = 0.5, face = "bold", size = 11),   # [07c] vertical to save width
        strip.placement = "outside",
        legend.position = "bottom",
        axis.text = element_text(color = "black"), axis.text.y = element_text(color = "black", size = 11),   # [07c] larger labels
        axis.text.x.top = element_text(color = "black", size = 15, face = "bold"),   # [07c] larger H1/H2 column headers (requested)
        plot.title = element_text(face = "bold"),
        legend.title = element_text(size = 13), legend.text = element_text(size = 12),   # [07c] larger scale legend (requested)
        plot.subtitle = element_text(color = "gray30", size = 9))

# ===========================================================================
# 4. Overall bars: share of all specs supporting each hypothesis
# ===========================================================================
# [07c] two bars per hypothesis: all dyads and pooled over the nine dyad rules
BAR_LEVELS <- c("H1\nall dyads", "H1\nnine dyad rules", "H2\nall dyads", "H2\nnine dyad rules")   # [07c] two-line labels
overall <- bind_rows(
  specs |> filter(bellig_rule == "none") |> group_by(hypothesis) |>
    summarise(sig = 100 * mean(sig_pos), .groups = "drop") |> mutate(scope = "all dyads"),
  specs |> group_by(hypothesis) |>
    summarise(sig = 100 * mean(sig_pos), .groups = "drop") |> mutate(scope = "nine dyad rules")) |>
  mutate(null = 100 - sig, bar = paste0(hypothesis, "\n", scope)) |>
  pivot_longer(c(sig, null), names_to = "outcome", values_to = "share") |>
  mutate(
    outcome = factor(outcome, levels = c("sig", "null"),
                     labels = c("Significant and positive",
                                "Null or significant negative")),
    hypothesis = factor(bar, levels = rev(BAR_LEVELS))   # H1 at top; column name kept for the plot code below
  )
stopifnot(nrow(overall) == 8)

p_bar <- ggplot(overall, aes(x = share, y = hypothesis, fill = outcome)) +
  geom_col(width = 0.62) +
  geom_text(aes(label = sprintf("%.1f%%", share)),
            position = position_stack(vjust = 0.5), color = "white",
            fontface = "bold", size = 3.4) +
  scale_fill_manual(values = c("Significant and positive" = "#1F2C74",
                               "Null or significant negative" = "#8B1E1E"),
                    name = NULL) +
  scale_x_continuous(labels = function(x) paste0(x, "%"),
                     expand = expansion(mult = c(0, 0))) +
  labs(title = "Overall multiverse result",
       subtitle = " ",   # [07c] dropped; blank line keeps panel tops aligned with the heatmap
       # old_subtitle = "Share of all specifications supporting each hypothesis,\nall dyads and pooled over the nine dyad rules.",
       x = "Share of specifications", y = NULL) +
  theme_minimal(base_size = 11) +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "bottom",
        legend.direction = "vertical",
        axis.text = element_text(color = "black"), axis.text.y = element_text(angle = 90, hjust = 0.5, vjust = 0.5, size = 13),   # [07c] vertical H1/H2 labels (requested)
        axis.text.x = element_text(size = 12), axis.title.x = element_text(size = 14),   # [07c] larger (requested)
        plot.title = element_text(face = "bold"),
        plot.subtitle = element_text(color = "gray30", size = 7),
        legend.text = element_text(size = 13),   # [07c] larger bar legend (requested)
        legend.key.size = unit(0.7, "cm"),
        plot.margin = margin(5.5, 18, 5.5, 5.5))   # [07c] room for the 100% tick label

p_fig1 <- p_heat + p_bar + plot_layout(widths = c(1.35, 1))
ggsave(file.path(FIG_DIR, "fig_multiverse_hypotheses_bellig.png"), p_fig1,
       width = 11.5, height = 13.5, dpi = 400, bg = "white")   # [07c] taller for row (k); narrower with vertical strips
cat("  -> fig_multiverse_hypotheses_bellig.png\n")

# ===========================================================================
# 5. Coefficient plot: how each decision shifts the downstream estimate,
#    faceted across H1 / H2 (free scales).
#    Outcome = the spec's coefficient, z-scored WITHIN hypothesis, so the
#    meta-regression coefficients read as SD-unit shifts in the downstream
#    estimate and are comparable within each facet.
# ===========================================================================
# [07c] The timing x FE grid enters the H2 meta-regression as estimated terms. Its rows are added with indicators for
# conflict timing (t, t-1, sum; t-1 = baseline), community FE, state + year + community clustering, and a timing-grid
# sample indicator (the grid is estimated on its own sample, so its t-1 TWFE cells differ slightly from the main grid's;
# the indicator absorbs that, and the timing / FE / clustering terms are identified within the grid).
tfe_rows <- tfe |>
  mutate(across(any_of(bool_cols), ~ as.logical(.x)),
         dv_logged = grepl("^log_", dv), iv_is_fatal = grepl("fatal", war_variable),
         tfe_sample = TRUE, timing = iv_term,
         fe_comm = fe_structure == "state_year_comm",
         clust_comm = clustering == "ccode_year_comm",
         model = dplyr::recode(clustering, ccode_year = "FE_full_ccode_year", year_only = "FE_full_year_only",
                               ccode_only = "FE_full_ccode_only", ccode_year_comm = "FE_full_ccode_year"),
         extra_control = "none", hypothesis = "H2")
stopifnot(!anyNA(tfe_rows$model), all(tfe_rows$timing %in% c("t", "t1", "joint_sum")),
          all(tfe_rows$clustering %in% c("ccode_year", "year_only", "ccode_only", "ccode_year_comm")))
main_flags <- function(d) d |> mutate(tfe_sample = FALSE, timing = "t1", fe_comm = FALSE, clust_comm = FALSE)
meta_df <- bind_rows(main_flags(specs), main_flags(wb |> mutate(bellig_rule = "none")), tfe_rows) |>   # [07c] bootstrap rows are all-dyads   # [2026-11-05] all specs; GDP enters as extra_ctl. [2026-11-08] + bootstrap rows (H2 inference term)
  group_by(hypothesis) |>
  mutate(coef_z = as.numeric(scale(coef))) |>
  ungroup() |>
  mutate(conflict = if_else(iv_is_fatal, "Fatal MIDs", "All MIDs"),
         recon    = dv_measure[dv],
         # Layer choices recoded as DEVIATIONS from the main spec, which includes
         # offensive + defensive + IGO and excludes trade. Each indicator is
         # FALSE in the main spec, so the main layer set is the baseline.
         drop_offensive = !has_offensive,
         drop_defensive = !has_defensive,
         drop_igo       = !has_igo,
         add_trade      = has_trade,
         # [2026-11-13] wild bootstrap p has a floor of 1/999 (999 draws); 20 rows store p = 0, which the 1e-300 floor
         # turned into -log10 p = 300 and inflated the inference term. Floor those at 1/999.
         neglog10p      = -log10(pmax(pval, if_else(dplyr::coalesce(model == "FE_full_wildboot", FALSE), 1 / 999, 1e-300))),   # significance, continuous
         sig_num        = as.integer(sig_pos),          # 0/1 for the LPM
         dv_log         = dv_logged,                    # [2026-11-08] outcome scale, raw = baseline
         bellig         = bellig_rule)                  # [07c] dyad rule, "none" = baseline
stopifnot(!anyNA(meta_df$tfe_sample), all(meta_df$timing[meta_df$hypothesis == "H1"] == "t1"))

# Reference = the MAIN specification for each hypothesis, so the intercept IS the
# main model and every coefficient reads as a deviation from it:
#   H1 (system) : [07c] sys_T1_ar1 (model 3 of Table 2: OLS with lagged DV), JS divergence, raw, unweighted, main layers, Fatal MIDs
#   H2 (state): FE_full_ccode_year (state+year TWFE, main-text Table 4 column 3; IV at t-1), Churn, raw,
#                  unweighted, main layers, Fatal MIDs   [2026-10-12: reference switched to fatal MIDs]
relevel_to <- function(x, ref) {
  x <- droplevels(factor(x))
  if (ref %in% levels(x)) relevel(x, ref = ref) else x
}

# Wald tidy: estimate +/- 1.96*SE (fast; avoids profiling confint() on 100k rows).
tidy_wald <- function(m) {
  ct <- summary(m)$coefficients
  tibble(term = rownames(ct), estimate = ct[, 1], se = ct[, 2]) |>
    filter(term != "(Intercept)") |>
    mutate(conf.low = estimate - 1.96 * se, conf.high = estimate + 1.96 * se) |>
    select(term, estimate, conf.low, conf.high)
}

run_meta <- function(df, hyp) {
  ref_model   <- if (hyp == "H1") "sys_T1_ar1" else "FE_full_ccode_year"   # [07c] model 3 of Tables 2 and 4
  ref_measure <- if (hyp == "H1") "JS divergence" else "Churn"  # JS = main-text DV (2026-09-26)
  
  df <- df |> mutate(
    model_f   = relevel_to(model, ref_model),
    transform = relevel_to(transform, "raw"),
    weighting = relevel_to(weighting, "unweighted"),
    conflict  = relevel_to(conflict, "Fatal MIDs"),   # [2026-10-12] main text reports fatal MIDs only
    recon     = relevel_to(recon, ref_measure),
    extra_ctl = relevel_to(dplyr::coalesce(extra_control, "none"), "none"),  # [2026-10-13]
    bellig    = relevel_to(bellig, "none"),                                    # [07c]
    timing_f  = relevel_to(timing, "t1")                                       # [07c] timing grid
  )
  stopifnot(levels(df$model_f)[1] == ref_model, levels(df$bellig)[1] == "none")   # [07c] baseline really is model 3
  # Keep only predictors that actually vary within this hypothesis.
  keep <- function(nm, x) if (length(unique(x[!is.na(x)])) > 1) nm else NULL
  preds <- c(keep("model_f", df$model_f),
             keep("transform", df$transform), keep("weighting", df$weighting),
             keep("conflict", df$conflict), keep("recon", df$recon),
             keep("drop_offensive", df$drop_offensive),
             keep("drop_defensive", df$drop_defensive),
             keep("drop_igo", df$drop_igo), keep("add_trade", df$add_trade),
             keep("extra_ctl", df$extra_ctl),   # [2026-10-13]
             keep("dv_log", df$dv_log),         # [2026-11-08]
             keep("bellig", df$bellig),         # [07c]
             keep("timing_f", df$timing_f), keep("fe_comm", df$fe_comm),         # [07c] timing grid terms
             keep("clust_comm", df$clust_comm), keep("tfe_sample", df$tfe_sample))
  rhs <- paste(preds, collapse = " + ")
  
  # THREE downstream outcomes, same predictors, same main-spec baseline:
  #   (1) the estimate (z-scored within hypothesis),
  #   (2) statistical significance as -log10(p),
  #   (3) P(significant-positive) via an LPM (avoids logit separation; reads as a
  #       change in the share of significant-positive specs).
  m_coef <- lm(as.formula(paste("coef_z ~",    rhs)), data = df)
  m_p    <- lm(as.formula(paste("neglog10p ~", rhs)), data = df)
  
  bind_rows(
    tidy_wald(m_coef) |> mutate(outcome = "Coefficient (SD)"),
    tidy_wald(m_p)    |> mutate(outcome = "-log10(p-value)")
  )   # [07c] P(significant-positive) column dropped (redundant with the heatmap)
}

model_labels <- c(
  OLS_lagged          = "Model: OLS with controls (no lagged DV)",   # [07c] no longer the H1 baseline
  sys_T1_ar1          = "Model: AR(1)",
  sys_T1_iv_only      = "Model: IV only",
  sys_T1_order_period = "Model: order/period",
  LP_h5               = "Model: Local projection (h5)",
  FE_full_year_only   = "Cluster: year only (TWFE)",    # [07c] same TWFE as the baseline, different clustering
  FE_full_ccode_only  = "Cluster: state only (TWFE)",
  full_noFE           = "Model: No FE (full controls)",
  full_noFE_year_poly = "Model: No FE + year poly",
  baseline_noFE       = "Model: Baseline (no FE)",
  FE_full_wildboot    = "Inference: Wild cluster bootstrap"   # [2026-11-08]
)

clean_term <- function(t) dplyr::case_when(
  str_starts(t, "model_f")         ~ unname(model_labels[str_remove(t, "^model_f")]),
  t == "transformlog1p"            ~ "Transform: log(1+x)",
  t == "transformzscore"           ~ "Transform: z-score",
  t == "weightingcontig_manhattan" ~ "Weight: Contiguity (Manhattan)",
  t == "weightingcontig_euclidean" ~ "Weight: Contiguity (Euclidean)",
  t == "weightingdist_manhattan"   ~ "Weight: Distance (Manhattan)",
  t == "weightingdist_euclidean"   ~ "Weight: Distance (Euclidean)",
  t == "conflictFatal MIDs"        ~ "Conflict: Fatal MIDs",
  t == "conflictAll MIDs"          ~ "Conflict: All MIDs",   # [2026-10-12] reference is now Fatal MIDs
  t == "drop_offensiveTRUE"        ~ "Drop offensive alliances",
  t == "drop_defensiveTRUE"        ~ "Drop defensive alliances",
  t == "drop_igoTRUE"              ~ "Drop shared IGOs",
  t == "add_tradeTRUE"             ~ "Add trade",
  str_starts(t, "recon")           ~ paste0("Measure: ", str_remove(t, "^recon")),
  str_starts(t, "extra_ctl")       ~ paste0("Extra control: ", str_remove(t, "^extra_ctl")),
  t == "dv_logTRUE"                ~ "Outcome: log(1+x)",                 # [2026-11-08]
  t == "timing_ft"                 ~ "Timing: conflict in year t",        # [2026-11-08] timing grid terms
  t == "timing_fjoint_sum"         ~ "Timing: years t and t-1 (sum)",
  t == "fe_fstate_year_comm"       ~ "FE: add community FE",
  t == "fe_commTRUE"               ~ "FE: add community FE",                 # [07c]
  t == "clust_commTRUE"            ~ "Cluster: state + year + community",    # [07c]
  t == "tfe_sampleTRUE"            ~ "Timing-grid sample",                   # [07c]
  t == "clust_fyear_only"          ~ "Cluster: year only",
  t == "clust_fccode_only"         ~ "Cluster: state only",
  t == "clust_fccode_year_comm"    ~ "Cluster: state + year + community",
  str_starts(t, "bellig")          ~ paste0("Dyads dropped: ", unname(bellig_lab[str_remove(t, "^bellig")])),   # [07c]
  TRUE ~ t
)

meta_coef <- purrr::imap_dfr(split(meta_df, meta_df$hypothesis),
                             ~ run_meta(.x, .y) |> mutate(hypothesis = .y)) |>
  mutate(term_clean = clean_term(term))
# [07c] Bootstrap rows exist for all dyads only, so the bootstrap term is estimated like-for-like on rule "none"
# (in the pooled fit it would be compared with a TWFE reference averaged over nine rules).
boot_none <- run_meta(meta_df |> filter(hypothesis == "H2", bellig_rule == "none", !tfe_sample), "H2") |>
  filter(term == "model_fFE_full_wildboot") |>
  mutate(hypothesis = "H2", term_clean = clean_term(term))
stopifnot(nrow(boot_none) == 2)
meta_coef <- meta_coef |> filter(term != "model_fFE_full_wildboot") |> bind_rows(boot_none)

# [07c] The separate timing-grid meta-regression is replaced by the timing terms inside the H2 regression above.
stopifnot(all(c("timing_ft", "timing_fjoint_sum", "fe_commTRUE", "clust_commTRUE", "tfe_sampleTRUE") %in% meta_coef$term[meta_coef$hypothesis == "H2"]))

# Order y-axis terms by decision category (model -> transform -> weighting ->
# conflict -> measure -> layers), consistent across all facets.
cat_rank <- function(tc) dplyr::case_when(
  str_starts(tc, "Model:")     ~ 1,
  str_starts(tc, "Transform:") ~ 2,
  str_starts(tc, "Weight:")    ~ 3,
  str_starts(tc, "Conflict:")  ~ 4,
  str_starts(tc, "Measure:")   ~ 5,
  str_starts(tc, "Timing:") | str_starts(tc, "FE:") | str_starts(tc, "Cluster:") | str_starts(tc, "Inference:") ~ 5.5,   # [07c]
  str_starts(tc, "Dyads dropped:") ~ 7,       # [07c]
  str_starts(tc, "Baseline")   ~ 0,           # [07c]
  TRUE                         ~ 6            # layer deviations (Drop.../Add...)
)
# [07c] draw the baseline (model 3 of Tables 2 and 4) as a point at 0 in every panel
BASE_LAB <- "Baseline: model 3 (main specification)"
meta_coef <- bind_rows(meta_coef,
                       tidyr::expand_grid(hypothesis = unique(meta_coef$hypothesis), outcome = unique(meta_coef$outcome)) |>
                         mutate(term = "baseline", term_clean = BASE_LAB, estimate = 0, conf.low = 0, conf.high = 0))
stopifnot(!anyNA(meta_coef$term_clean))
term_levels <- meta_coef |> distinct(term_clean) |>
  mutate(r = cat_rank(term_clean)) |> arrange(r, term_clean) |> pull(term_clean)

meta_coef <- meta_coef |>
  mutate(hypothesis = factor(hypothesis, levels = HYP),   # [07c] timing grid folded into H2
         term_clean = factor(term_clean, levels = rev(term_levels)),
         outcome    = factor(outcome, levels = c("Coefficient (SD)",
                                                 "-log10(p-value)")))

# facet_grid free scales: x free per column (the 3 outcomes have very different
# scales), y free per row (each hypothesis shows only its own terms).
p_fig2 <- ggplot(meta_coef, aes(x = estimate, y = term_clean)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray55") +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), height = 0.2,
                 color = "#2C6CB0") +
  geom_point(size = 1.9, color = "#2C6CB0") +
  facet_grid(hypothesis ~ outcome, scales = "free", space = "free_y") +
  labs(x = "Effect vs. main specification", y = NULL) +   # [2026-11-13] in-plot title/subtitle removed (caption carries them)
  theme_bw(base_size = 10) +
  theme(panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "gray92", color = "gray70"),
        strip.text = element_text(face = "bold", size = 12),
        axis.text.y = element_text(size = 7.5),
        axis.text.x = element_text(size = 11),
        axis.title.x = element_text(size = 13),   # [07c] larger (requested)
        plot.title = element_text(face = "bold", size = 15),
        plot.subtitle = element_text(color = "gray30", size = 11))

ggsave(file.path(FIG_DIR, "fig_multiverse_meta_by_hyp_bellig.png"), p_fig2,
       width = 10, height = 13, dpi = 400, bg = "white")   # [07c] two outcome columns, two hypothesis rows
cat("  -> fig_multiverse_meta_by_hyp_bellig.png\n")

# [07c] support shares by rule, all dyads and pooled, for the text
support <- bind_rows(
  specs |> group_by(hypothesis, bellig_rule) |> summarise(n = n(), pct_sig_pos = 100 * mean(sig_pos),
                                                         pct_pos = 100 * mean(coef > 0), .groups = "drop"),
  specs |> group_by(hypothesis) |> summarise(n = n(), pct_sig_pos = 100 * mean(sig_pos),
                                             pct_pos = 100 * mean(coef > 0), .groups = "drop") |>
    mutate(bellig_rule = "pooled_nine_rules"))
readr::write_csv(support, file.path(OUTPUT_DIR, "multiverse_support_by_bellig_rule.csv"))
cat("\nDone. Support by hypothesis and dyad rule (one-sided positive, p<0.10):\n")
print(support |> mutate(across(c(pct_sig_pos, pct_pos), ~ round(.x, 1))), n = Inf)
