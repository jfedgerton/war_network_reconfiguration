###############################################################################
# 27_sensitivity_plots.R
# Standalone plotting script for the Sensitivity Analyses section.
# Loads the CSVs written by the analysis scripts and produces four figures:
#
#   (1) Reverse-timing placebo   <- robustness_placebo_reverse_timing.csv
#   (2) Permutation null          <- robustness_placebo_permutation.csv
#   (3) Alternative-mechanism med <- robustness_altmech_mediation.csv  (run 16 first)
#   (4) Multiverse coefficient box<- robustness_specifications.csv
#
# Conventions:
#   - Variable-of-interest labels DROP the "log(1 + .)" wrapper (noted in text).
#   - Plots 1-3 are facet-wrapped by level (System / State).
#   - x-axis text is horizontal; long box labels are broken onto two lines.
###############################################################################

library(tidyverse)
select <- dplyr::select

set.seed(123)

DATA_DIR <- here::here("data", "processed")
FIG_DIR  <- here::here("code", "03_figures", "output", "manuscript")
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

ACCENT  <- "#2C6CB0"   # blue  (observed / significant-positive)
ACCENT2 <- "#B03366"   # red   (null)
GREY    <- "gray55"

theme_sens <- function(base_size = 12) {
  theme_bw(base_size = base_size) %+replace%
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "gray90", linewidth = 0.3),
      plot.title    = element_text(size = base_size + 2, face = "bold",
                                     hjust = 0, margin = margin(b = 4)),
      plot.subtitle = element_text(size = base_size - 2, color = "gray30",
                                     hjust = 0, margin = margin(b = 8)),
      axis.text     = element_text(color = "black", size = base_size - 1),
      axis.text.x   = element_text(angle = 0, hjust = 0.5),   # horizontal
      axis.title    = element_text(color = "black", size = base_size),
      strip.background = element_rect(fill = "gray92", color = "gray70"),
      strip.text    = element_text(face = "bold", size = base_size, color = "black"),
      legend.position = "bottom",
      panel.spacing = unit(1.1, "lines"),
      plot.margin   = margin(10, 12, 8, 10)
    )
}

save_fig <- function(filename, plot, width = 8, height = 5, dpi = 400) {
  ggsave(file.path(FIG_DIR, filename), plot,
         width = width, height = height, dpi = dpi, bg = "white")
  cat(sprintf("  -> %s\n", filename))
}

lvl_factor <- function(x) factor(if_else(x == "system", "System level", "State level"),
                                 levels = c("System level", "State level"))

# ===========================================================================
# (1) Reverse-timing placebo
# ===========================================================================
cat("(1) Placebo plot...\n")

# Panels are labeled by HYPOTHESIS and ordered H1 -> H2. Labels drop the
# "log(1 + .)" wrapper (noted in the text instead). [2026-10-12] H3 (partner war exposure) removed.
hyp_levels <- c("H1: Systemic fatal MID count (t-1)",
                "H2: State fatal MID onset (t-1)")

placebo <- read_csv(file.path(DATA_DIR, "robustness_placebo_reverse_timing.csv"),
                    show_col_types = FALSE) |>
  mutate(
    panel = case_when(
      term == "log_fatal_mid_count_lag1"  ~ "H1: Systemic fatal MID count (t-1)",
      term == "log_fatal_onset_lag1"      ~ "H2: State fatal MID onset (t-1)",
      TRUE                                ~ NA_character_
    ),
    panel = factor(panel, levels = hyp_levels),
    # [2026-10-12] the observed (horizon 0) estimate is the fatal t-1 focal IV at both levels
    horizon_lab = factor(if_else(horizon == 0, "t-1", paste0("t+", horizon)),
                         levels = c("t-1", "t+1", "t+3", "t+5", "t+10")),
    kind = if_else(horizon == 0, "Observed (t-1)", "Placebo (future)")
  ) |>
  filter(!is.na(panel))

p1 <- ggplot(placebo, aes(x = horizon_lab, y = estimate, color = kind)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = GREY) +
  geom_errorbar(aes(ymin = ci_low, ymax = ci_high), width = 0.15, linewidth = 0.7) +
  geom_point(size = 2.6) +
  scale_color_manual(values = c("Observed (t-1)" = ACCENT, "Placebo (future)" = GREY),
                     name = NULL) +
  facet_wrap(~ panel, scales = "free_y", ncol = 2) +
  labs(
    title = "Reverse-timing placebo tests",
    subtitle = "Observed conflict effect vs. future placebo leads; future leads should be null.",
    x = "Conflict lead", y = "Coefficient (95% CI)"
  ) +
  theme_sens()

save_fig("fig_sens_placebo.png", p1, width = 7.5, height = 4.5)

# ===========================================================================
# (2) Permutation null distribution + observed value
# ===========================================================================
cat("(2) Permutation plot...\n")

perm <- read_csv(file.path(DATA_DIR, "robustness_placebo_permutation.csv"),
                 show_col_types = FALSE) |>
  mutate(level_label = lvl_factor(level))

# Observed coefficient = reverse-timing estimate at horizon 0 (the focal IV, at t-1).
# Join on level AND term so the join stays one-to-one if more IVs are added.
# The permutation test only covers MID onset (state) and MID count (system).
observed <- placebo |>
  filter(horizon == 0) |>
  distinct(level, term, estimate) |>
  rename(obs_est = estimate)

perm <- perm |> left_join(observed, by = c("level", "term"))
perm$level_label <- lvl_factor(perm$level)

# One-sided permutation p: share of permuted coefficients >= observed.
perm_p <- perm |>
  group_by(level_label) |>
  summarise(obs_est = first(obs_est),
            p_perm  = mean(estimate >= first(obs_est), na.rm = TRUE),
            x_text  = max(estimate, na.rm = TRUE),
            .groups = "drop")

p2 <- ggplot(perm, aes(x = estimate)) +
  geom_histogram(bins = 40, fill = "gray75", color = "white") +
  geom_vline(data = perm_p, aes(xintercept = obs_est),
             color = ACCENT, linewidth = 1) +
  geom_text(data = perm_p,
            aes(x = x_text, y = Inf,
                label = sprintf("observed\n(perm. p = %.3f)", p_perm)),
            hjust = 1, vjust = 1.4, size = 3, color = ACCENT) +
  facet_wrap(~ level_label, scales = "free", ncol = 2) +
  labs(
    title = "Permutation tests",
    subtitle = "Null distribution of the coefficient under randomly reassigned conflict timing; line = observed.",
    x = "Coefficient under permuted conflict timing", y = "Count"
  ) +
  theme_sens()

save_fig("fig_sens_permutation.png", p2, width = 8.5, height = 5)

# ===========================================================================
# (3) Alternative-mechanism mediation (excl. battle deaths)
# ===========================================================================
cat("(3) Mediation plot...\n")

med <- read_csv(file.path(DATA_DIR, "robustness_altmech_mediation.csv"),
                show_col_types = FALSE)

# Exclude battle deaths; restrict to the log focal IVs if the column is present.
med <- med |> filter(!grepl("batdeath", mediator, ignore.case = TRUE))
if ("iv" %in% names(med)) {
  log_ivs <- c("log_fatal_onset_lag1", "log_fatal_mid_count_lag1")
  if (any(med$iv %in% log_ivs)) med <- med |> filter(iv %in% log_ivs)
}
stopifnot("indirect_ci_low"  %in% names(med),
          "indirect_ci_high" %in% names(med))

# Recode on the `mediator` column (the variable name) rather than the label
# string, so these hold regardless of how script 16 happens to label them.
med_lbl <- c(
  # state level
  abs_d_polity2  = "|Δ Polity|",
  mp_flip        = "Major power\nflip",
  tc_event       = "Territorial\nchange",
  delta_cinc     = "Δ CINC",
  # system level
  abs_d_n_states = "|Δ States|",
  n_tc_events    = "Territorial\nchange events"
)

med <- med |>
  mutate(level_label = lvl_factor(level),
         mediator_label = if_else(mediator %in% names(med_lbl),
                                  unname(med_lbl[mediator]),
                                  as.character(mediator_label)),
         mediator_label = fct_inorder(mediator_label))

p3 <- ggplot(med, aes(x = indirect_ab, y = mediator_label)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = GREY) +
  geom_errorbarh(aes(xmin = indirect_ci_low, xmax = indirect_ci_high),
                 height = 0.16, linewidth = 0.7, color = ACCENT) +
  geom_point(size = 2.7, color = ACCENT) +
  facet_wrap(~ level_label, scales = "free", ncol = 2) +
  # Fixed notation instead of the default "0e+00 / 2e-04" scientific labels.
  # accuracy = 1e-5 keeps the small state-level effects legible; drop0trailing
  # keeps the larger system-level breaks clean (0.02, not 0.02000).
  scale_x_continuous(labels = scales::label_number(accuracy = 1e-5,
                                                   drop0trailing = TRUE)) +
  labs(
    title = "Alternative-mechanism mediation",
    subtitle = "Indirect effect (a b path product) with 95% bootstrap CI; a CI overlapping 0 rules out that channel.",
    x = "Indirect effect", y = NULL
  ) +
  theme_sens()

save_fig("fig_sens_mediation.png", p3, width = 8, height = 3)

# ===========================================================================
# (4) Multiverse coefficient boxplot (system + state), sig-colored jitter
# ===========================================================================
cat("(4) Multiverse boxplot (reading large specifications file)...\n")

library(patchwork)

# ---- p-value column ---------------------------------------------------------
# "pval"    = iid  (what the spec-curve percentages in the text are based on)
# "nw_pval" = Newey-West;  "hc3_pval" = heteroskedasticity-robust
# NOTE: this choice moves the numbers a lot. Under iid the three coefficients
# are 52.9% / 81.5% / 76.2%; under NW they are 59.7% / 99.7% / 98.4%. Keep this
# consistent with whatever the manuscript reports.
P_COL <- "pval"

specs <- read_csv(file.path(DATA_DIR, "robustness_specifications.csv"),
                  col_select = c(level, war_variable, coef, dv, weighting,
                                 transform, iv_is_fatal, model,
                                 has_offensive, has_defensive,
                                 has_igo, has_trade, has_dca,
                                 any_of("extra_control"),
                                 all_of(P_COL)),
                  show_col_types = FALSE) |>
  rename(pval = all_of(P_COL)) |>
  # [2026-10-12] VAR Granger rows carry F-statistics (and the reverse test is Network -> Conflict); they are
  # not conflict coefficients, so they are excluded here exactly as in 07's hypothesis figures.
  filter(!grepl("^VAR_", model)) |>
  filter(!is.na(pval), !is.na(coef))

# [2026-10-13] Extra controls (GDP, strength, spillover, community spillover) are run only on FE_full x two IVs x
# four churn DVs, so pooling them would overweight those cells (about a third of state rows). The hypothesis
# shares and figures use the baseline rows (extra_control == "none"); each control's share is printed
# separately on the cells where it is estimated (07 also writes robustness_extra_control_shares.csv).
if (!"extra_control" %in% names(specs)) specs$extra_control <- "none"
specs <- specs |> mutate(extra_control = coalesce(extra_control, "none"))
extra_ivs <- unique(specs$war_variable[specs$level == "state" & specs$extra_control != "none"])   # IVs with a control arm
extra_cells <- specs |>
  filter(level == "state", grepl("^FE_full_", model), war_variable %in% extra_ivs)
if (n_distinct(extra_cells$extra_control) > 1) {
  cat("  State H2 share significant-positive (one-sided p < 0.1) by extra control, on the cells where controls run:\n")
  print(extra_cells |>
          group_by(extra_control) |>
          summarise(n = n(), pct_sig_pos = 100 * mean(coef > 0 & pval / 2 < 0.1), .groups = "drop"),
        n = Inf)
}
specs <- specs |> filter(extra_control == "none")

# The three hypotheses. Each POOLS the conflict-measure variants that the
# multiverse varies (MID vs fatal MID; raw vs logged; unweighted vs weighted),
# which is measurement decision (d) in the manuscript. This reproduces the
# reported figures. [2026-10-12] H3 (partner war exposure) removed; the multiverse no longer contains it.
measure_map <- tibble::tribble(
  ~level,   ~war_variable,                     ~hyp,
  # H1 — systemic conflict
  "system", "mid_count",                       "H1: Systemic MID Count",
  "system", "log_mid_count",                   "H1: Systemic MID Count",
  "system", "fatal_mid_count",                 "H1: Systemic MID Count",
  "system", "log_fatal_mid_count",             "H1: Systemic MID Count",
  # H2 — belligerent's own conflict
  "state",  "mid_onset_count",                 "H2: State MID Onset",
  "state",  "log_mid_onset_count",             "H2: State MID Onset",
  "state",  "fatal_mid_onset_count",           "H2: State MID Onset",
  "state",  "log_fatal_mid_onset_count",       "H2: State MID Onset"
)
# unique(): measure_map now has 4 IV variants per hypothesis, so the hyp column
# repeats each label four times. factor(levels = ) rejects duplicates.
hyp_levels  <- unique(measure_map$hyp)
hyp_short   <- c("H1: Systemic MID Count" = "H1",
                 "H2: State MID Onset" = "H2")

# Significance rule: ONE-SIDED test at the 0.1 level.
#   BLUE = supports the theory (positive AND significant)
#   RED  = does not (null OR significant negative)
specs2 <- specs |>
  inner_join(measure_map, by = c("level", "war_variable")) |>
  mutate(
    p_one_sided = if_else(coef > 0, pval / 2, 1 - pval / 2),
    sig_pos     = (coef > 0) & (p_one_sided < 0.1),
    hyp         = factor(hyp, levels = hyp_levels),
    hyp_s       = factor(hyp_short[as.character(hyp)], levels = unname(hyp_short))
  )

# ===========================================================================
# (4a) SPECIFICATION-INFLUENCE HEATMAP
#      Five varied features become five ROW-GROUPS; x = hypothesis;
#      fill = share of specs with that choice that are significant-positive.
# ===========================================================================
mk <- function(d, feature, lvl) {
  d |> mutate(feature = feature, lvl = lvl) |>
    select(hyp_s, feature, lvl, sig_pos)
}

# The five row-groups map ONE-TO-ONE onto the five measurement decisions
# described in the manuscript:
#   (a) how to measure network reconfiguration   -> "(a) Reconfiguration measure"
#   (b) which network layers are included        -> "(b) Network layer included"
#   (c) transformation of the variables          -> "(c) Variable transform"
#   (d) MIDs as the conflict measure             -> "(d) Conflict measure"
#   (e) how tie strength is aggregated/weighted  -> "(e) Edge weighting"
# Hyphenated so the longest LINE is short: with vertical strip text, the label's
# longest line is what has to fit the facet height, not the total character count.
F_A <- "(a) Reconfig-\nuration measure"
F_B <- "(b) Network layer\nincluded"
F_C <- "(c) Variable\ntransform"
F_D <- "(d) Conflict\nmeasure"
F_E <- "(e) Edge\nweighting"

hm <- bind_rows(
  # (a) reconfiguration measure — collapse the log_ prefix (transform is its own group)
  mk(specs2, F_A, str_remove(specs2$dv, "^log_")),
  # (b) network layers included
  mk(filter(specs2, has_offensive), F_B, "Offensive alliances"),
  mk(filter(specs2, has_defensive), F_B, "Defensive alliances"),
  mk(filter(specs2, has_igo),       F_B, "Shared IGOs"),
  mk(filter(specs2, has_trade),     F_B, "Trade"),
  # DCA facet removed [2026-07-15]: 07 no longer emits DCA specs (has_dca is
  # all FALSE), so this row would be empty. DCA lives only in the contamination
  # multiverse.
  # (c) variable transformation
  mk(specs2, F_C, specs2$transform),
  # (d) conflict measure: MID vs fatal MID
  mk(specs2, F_D, if_else(specs2$iv_is_fatal, "Fatal MIDs", "All MIDs")),
  # (e) edge weighting / aggregation
  mk(specs2, F_E, specs2$weighting)
)

# Human-readable labels: no underscores anywhere in the rendered figure.
lbl_map <- c(
  # (a) reconfiguration measures
  "churn"              = "Churn",
  "churn_alt"          = "Churn (alternate)",
  "kl_divergence"      = "KL divergence",
  "js_divergence"      = "JS divergence",
  "edge_turnover_pct"  = "Edge turnover",
  "degree_kl"          = "Degree KL",
  "degree_wasserstein" = "Degree Wasserstein",
  # (c) transforms
  "raw"                = "Raw",
  "log1p"              = "log(1 + x)",
  "zscore"             = "z-score",
  # (e) edge weighting
  "unweighted"         = "Unweighted",
  "dist_manhattan"     = "Distance (Manhattan)",
  "dist_euclidean"     = "Distance (Euclidean)",
  "contig_manhattan"   = "Contiguity (Manhattan)",
  "contig_euclidean"   = "Contiguity (Euclidean)"
)
clean_lbl <- function(x) {
  out <- unname(lbl_map[x])
  ifelse(is.na(out), x, out)   # layer + conflict-measure labels are already clean
}

hm <- hm |>
  mutate(lvl = clean_lbl(lvl)) |>
  group_by(hyp_s, feature, lvl) |>
  summarise(pct = 100 * mean(sig_pos), n = n(), .groups = "drop") |>
  # Complete the grid so structurally inapplicable cells become explicit NAs
  # rather than blank gaps. The reconfiguration measures are LEVEL-SPECIFIC:
  # churn exists only at the state level (H2); the divergence measures only
  # at the system level (H1). Those cells are "not applicable", not "missing",
  # and are drawn in grey with an em dash.
  tidyr::complete(hyp_s, tidyr::nesting(feature, lvl)) |>
  mutate(feature = factor(feature, levels = c(F_A, F_B, F_C, F_D, F_E)))

p4a <- ggplot(hm, aes(x = hyp_s, y = fct_rev(lvl), fill = pct)) +
  geom_tile(color = "white", linewidth = 0.6) +
  geom_text(data = filter(hm, !is.na(pct)),
            aes(label = sprintf("%.0f", pct)), size = 2.9, color = "black") +
  geom_text(data = filter(hm, is.na(pct)),
            aes(label = "—"), size = 3.2, color = "grey45") +
  scale_fill_gradient2(low = ACCENT2, mid = "grey92", high = ACCENT,
                       midpoint = 50, limits = c(0, 100),
                       na.value = "grey85",          # inapplicable cells = grey
                       name = "% significant positive") +
  # Strips on the RIGHT (facet_grid default), text rotated vertical.
  # NOTE: no `space = "free_y"`. With free space, each facet's height is
  # proportional to its row count, so the 2-row (d) block gets ~9% of the panel
  # and the rotated strip label gets CLIPPED ("Conf. n sures"). Equal-height
  # facets give every rotated label room to render. Trade-off: tiles are taller
  # in the blocks with fewer rows.
  facet_grid(feature ~ ., scales = "free_y") +
  guides(fill = guide_colorbar(direction = "horizontal",
                               title.position = "top", title.hjust = 0.5,
                               barwidth = 11, barheight = 0.6)) +
  labs(title = "Which specification choices matter?",
       subtitle = "Share of specifications significant-positive, by decision. Dashes = measure not defined at that level.",
       x = NULL, y = NULL) +
  theme_sens() +
  theme(legend.position = "bottom",
        strip.text.y = element_text(angle = -90, size = 8, face = "bold"),
        # extra vertical breathing room between the decision blocks
        panel.spacing.y = unit(0.65, "lines"),
        panel.grid = element_blank())

# ===========================================================================
# (4b) OVERALL STACKED PROPORTION BAR
# ===========================================================================
DARKBLUE <- "#00008B"
DARKRED  <- "#8B0000"

SIG_LAB  <- "Significant positive (one-sided p < 0.1)"
NULL_LAB <- "Null or significant negative"

prop_df <- specs2 |>
  mutate(status = factor(if_else(sig_pos, SIG_LAB, NULL_LAB),
                         levels = c(SIG_LAB, NULL_LAB))) |>
  count(hyp, hyp_s, status, name = "n") |>
  group_by(hyp_s) |>
  mutate(prop = n / sum(n), n_total = sum(n)) |>
  ungroup()

cat(sprintf("  Share significant-positive (one-sided p < 0.1), using '%s':\n", P_COL))
print(prop_df |>
        filter(status == SIG_LAB) |>
        select(hyp, n_sig_pos = n, n_total, prop),
      n = Inf, digits = 3)

# reverse = TRUE puts the FIRST factor level (blue) on the LEFT of the stack.
p4b <- ggplot(prop_df, aes(x = prop, y = fct_rev(hyp_s), fill = status)) +
  geom_col(width = 0.62, position = position_stack(reverse = TRUE)) +
  geom_text(aes(label = if_else(prop >= 0.05,
                                scales::percent(prop, accuracy = 0.1), "")),
            position = position_stack(vjust = 0.5, reverse = TRUE),
            color = "white", fontface = "bold", size = 3.4) +
  scale_x_continuous(labels = scales::percent_format(accuracy = 1),
                     limits = c(0, 1), breaks = seq(0, 1, 0.25),
                     expand = expansion(mult = c(0, 0))) +
  scale_fill_manual(values = c(DARKBLUE, DARKRED),
                    breaks = c(SIG_LAB, NULL_LAB),   # blue first in the legend
                    name = NULL) +
  guides(fill = guide_legend(nrow = 2)) +
  labs(title = "Overall multiverse result",
       subtitle = "Share of all specifications supporting each hypothesis.",
       x = "Share of specifications", y = NULL) +
  theme_sens()

p4 <- p4a + p4b + plot_layout(widths = c(1.25, 1))

save_fig("fig_sens_multiverse.png", p4, width = 13, height = 8.3)

cat("\nDone. Four figures written to", FIG_DIR, "\n")
