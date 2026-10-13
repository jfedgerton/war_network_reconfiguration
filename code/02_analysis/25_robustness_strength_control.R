###############################################################################
# 25_robustness_strength_control.R
# Robustness: re-run the MID-onset state-level spec ladder with
# log1p(strength_t1) as an additional control, across the full multiverse
# of layer combos × weighting schemes.
#
# RATIONALE: Weighted Jaccard churn is already normalized by network size
# (sum_max in the denominator), but states with different network positions
# may still have systematically different churn dynamics. Including lagged
# network strength (sum of edge weights at t-1) as a control absorbs any
# residual confounding from network size / density and asks: does MID onset
# increase churn BEYOND what the state's prior connectivity predicts?
#
# DESIGN:
#   1. Loop over 24 layer combos × 5 weighting methods (same grid as 03).
#   2. For each cell, compute edge weights, then:
#      a. Canonical churn (weighted Jaccard)
#      b. Strength at t-1 (sum of edge weights to all partners at t-1)
#   3. Run the MID-onset spec ladder WITH and WITHOUT log1p(strength_t1):
#      - (1) Bivariate TWFE
#      - (2) +Controls TWFE
#      - (3) +Controls TWFE + strength
#   4. Record coefficients to compare: does the MID coefficient change
#      when conditioning on prior network strength?
#
# All helper functions defined locally. No shared code is modified.
#
# Inputs:  data/processed/dyad_year_network_layers.rds
#          data/processed/state_network_churn.rds
#          data/processed/network_metadata.rds
# Outputs:
#   data/processed/robustness_strength_control_results.csv
#   data/processed/robustness_strength_control_multiverse.csv
#   code/04_figures/output/manuscript/fig_ms_strength_control.png
#   code/04_figures/output/manuscript/fig_ms_strength_multiverse.png
###############################################################################

library(tidyverse)
library(fixest)
library(peacesciencer)
library(patchwork)
select <- dplyr::select

set.seed(123)

OUTPUT_DIR <- here::here("data", "processed")
FIG_DIR    <- here::here("code", "04_figures", "output", "manuscript")
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

ACCENT  <- "#2C6CB0"
ACCENT2 <- "#B03366"

theme_pres <- function(base_size = 12) {
  theme_bw(base_size = base_size) %+replace%
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "gray90", linewidth = 0.3),
      plot.title       = element_text(size = base_size + 2, face = "bold",
                                       hjust = 0, margin = margin(b = 8)),
      axis.text        = element_text(color = "black", size = base_size - 1),
      axis.title       = element_text(color = "black", size = base_size),
      axis.ticks       = element_line(color = "gray40"),
      strip.background = element_rect(fill = "gray92", color = "gray70"),
      strip.text       = element_text(face = "bold", size = base_size, color = "black"),
      legend.position  = "bottom",
      legend.text      = element_text(size = base_size - 1),
      legend.title     = element_text(size = base_size, face = "bold"),
      plot.margin      = margin(10, 12, 8, 10)
    )
}

save_pres <- function(filename, plot, width = 9, height = 6, dpi = 400) {
  ggsave(file.path(FIG_DIR, filename), plot,
         width = width, height = height, dpi = dpi, bg = "white")
  cat(sprintf("  -> %s\n", filename))
}

# ============================================================================
# 1. Load data
# ============================================================================

cat("Loading data...\n")
dyad_data   <- readRDS(file.path(OUTPUT_DIR, "dyad_year_network_layers.rds"))
state_churn <- readRDS(file.path(OUTPUT_DIR, "state_network_churn.rds"))
metadata    <- readRDS(file.path(OUTPUT_DIR, "network_metadata.rds"))

# ============================================================================
# 2. Build analysis frame (same controls as main spec)
# ============================================================================

state_attrs <- dyad_data |>
  select(ccode = ccode1, year, cinc = cinc1, major_power = cowmaj1) |>
  distinct(ccode, year, .keep_all = TRUE) |>
  mutate(major_power = as.integer(coalesce(major_power, 0L) == 1))

state_churn <- state_churn |>
  left_join(state_attrs, by = c("ccode", "year"))

# Polity
state_years <- create_stateyears(
  subset_years = min(state_churn$year):max(state_churn$year))
attr(state_years, "ps_data_type") <- "state_year"
attr(state_years, "ps_system")    <- "cow"
state_years <- state_years |> add_democracy()

polity <- state_years |>
  select(ccode, year, polity2) |>
  mutate(democracy = as.integer(polity2 >= 6),
         autocracy = as.integer(polity2 <= -6))

state_churn <- state_churn |>
  left_join(polity, by = c("ccode", "year"))

# Partner war exposure
partner_sets <- dyad_data |>
  filter(main_weight > 0) |>
  select(ccode1, ccode2, year, w = main_weight) |>
  bind_rows(
    dyad_data |> filter(main_weight > 0) |>
      select(ccode1 = ccode2, ccode2 = ccode1, year, w = main_weight)
  )

at_war_dyad <- dyad_data |>
  filter(cowmidongoing == 1) |>
  select(ccode1, ccode2, year) |>
  pivot_longer(c(ccode1, ccode2), values_to = "ccode") |>
  distinct(ccode, year) |>
  mutate(at_war = 1L)

at_war_dyad_fatal <- dyad_data |>
  filter(cowmidongoing == 1, cowmidfatal == 1) |>
  select(ccode1, ccode2, year) |>
  pivot_longer(c(ccode1, ccode2), values_to = "ccode") |>
  distinct(ccode, year) |>
  mutate(at_war_fatal = 1L)

pwe <- partner_sets |>
  mutate(year_next = year + 1L) |>
  left_join(at_war_dyad,       by = c("ccode2" = "ccode", "year_next" = "year")) |>
  left_join(at_war_dyad_fatal, by = c("ccode2" = "ccode", "year_next" = "year")) |>
  mutate(at_war       = replace_na(at_war, 0L),
         at_war_fatal = replace_na(at_war_fatal, 0L)) |>
  group_by(ccode1, year_next) |>
  summarise(
    partner_war_frac_weighted       = sum(at_war * w)       / sum(w),
    partner_fatal_war_frac_weighted = sum(at_war_fatal * w) / sum(w),
    .groups = "drop"
  ) |>
  rename(ccode = ccode1, year = year_next)

state_churn <- state_churn |>
  left_join(pwe, by = c("ccode", "year")) |>
  mutate(partner_war_frac_weighted       = replace_na(partner_war_frac_weighted, 0),
         partner_fatal_war_frac_weighted = replace_na(partner_fatal_war_frac_weighted, 0))

# IV transforms
state_churn <- state_churn |>
  mutate(
    log_churn                 = log1p(churn),
    log_mid_onset_count       = log1p(mid_onset_count),
    log_fatal_mid_onset_count = log1p(fatal_mid_onset_count)
  )
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

# [2026-10-12] partner war exposure removed from the control set (H3 dropped from the paper).
CTRL_MID   <- "democracy + autocracy + cinc + major_power + is_isolated_t1"
CTRL_FATAL <- "democracy + autocracy + cinc + major_power + is_isolated_t1"

# ============================================================================
# 3. Headline results: main spec with and without strength control
# ============================================================================
# Use canonical main_weight for the headline, then multiverse below.

cat("\n=== Headline: canonical main_weight with strength control ===\n")

# Compute strength at t-1 from main_weight
# [2026-10-12] dyad_data stores both directions of each pair, so collapse to one edge per unordered
# pair (max weight, as in 07 and the community graph) before summing endpoint strengths. Otherwise a
# single tie of weight w contributes 2w, which log1p() does not absorb as a constant.
pairs_main <- dyad_data |>
  filter(main_weight > 0, ccode1 != ccode2) |>
  transmute(a = pmin(ccode1, ccode2), b = pmax(ccode1, ccode2), year, w = main_weight) |>
  group_by(a, b, year) |>
  summarise(w = max(w, na.rm = TRUE), .groups = "drop")
strength_main <- bind_rows(
  pairs_main |> transmute(ccode = a, year, w),
  pairs_main |> transmute(ccode = b, year, w)
) |>
  group_by(ccode, year) |>
  summarise(strength = sum(w, na.rm = TRUE), .groups = "drop") |>
  mutate(year_next = year + 1L) |>
  select(ccode, year = year_next, strength_t1 = strength)

state_churn <- state_churn |>
  left_join(strength_main, by = c("ccode", "year")) |>
  mutate(strength_t1 = replace_na(strength_t1, 0),
         log_strength_t1 = log1p(strength_t1))

# [2026-10-12] IVs at t-1 (headline timing), from the full state-year universe (defined for every transition row).
analysis_sample <- state_churn |>
  filter(!is.na(churn),
         !is.na(log_onset_lag1),
         !is.na(log_fatal_onset_lag1))

cat(sprintf("Analysis sample: %s state-years\n",
            format(nrow(analysis_sample), big.mark = ",")))

# Fit one spec — IV name is parameterised so the same helper handles
# both MID-onset and fatal-MID-onset.
fit_spec <- function(formula, data, vcov_arg, spec_label, iv_var) {
  fit <- tryCatch(feols(formula, data = data, vcov = vcov_arg),
                  error = function(e) NULL)
  if (is.null(fit)) return(NULL)
  ct <- coeftable(fit)
  if (!iv_var %in% rownames(ct)) return(NULL)
  tibble(
    iv       = iv_var,
    spec     = spec_label,
    estimate = ct[iv_var, "Estimate"],
    se       = ct[iv_var, "Std. Error"],
    pval     = ct[iv_var, "Pr(>|t|)"],
    n_obs    = nobs(fit),
    ci_low   = ct[iv_var, "Estimate"] - 1.96 * ct[iv_var, "Std. Error"],
    ci_high  = ct[iv_var, "Estimate"] + 1.96 * ct[iv_var, "Std. Error"]
  )
}

# Run the 5-spec ladder for one IV (severity-routed control set).
run_ladder <- function(iv_var, ctrl_set) {
  # [2026-10-12] Specs 3 and 5 ADD log_strength_t1 and keep is_isolated_t1. The two are not collinear
  # (strength varies among non-isolates; the dummy is 0 for all of them), so dropping the dummy changed two
  # controls at once. Any genuine cell-level collinearity is left to the estimator.
  bind_rows(
    fit_spec(as.formula(paste0("churn ~ ", iv_var, " | ccode + year")),
             analysis_sample, ~ ccode + year, "(1) Bivariate TWFE", iv_var),
    fit_spec(as.formula(paste0("churn ~ ", iv_var, " + ", ctrl_set,
                                " | ccode + year")),
             analysis_sample, ~ ccode + year, "(2) +Controls TWFE", iv_var),
    # Spec (3): add log_strength_t1 to the spec-2 controls.
    fit_spec(as.formula(paste0("churn ~ ", iv_var, " + log_strength_t1 + ",
                                ctrl_set, " | ccode + year")),
             analysis_sample, ~ ccode + year,
             "(3) +Controls +Strength TWFE", iv_var),
    fit_spec(as.formula(paste0("churn ~ ", iv_var, " + ", ctrl_set,
                                " + poly(year, 2)")),
             analysis_sample, ~ ccode, "(4) +Controls, no FE", iv_var),
    # Spec (5): add log_strength_t1 to the spec-4 controls.
    fit_spec(as.formula(paste0("churn ~ ", iv_var, " + log_strength_t1 + ",
                                ctrl_set, " + poly(year, 2)")),
             analysis_sample, ~ ccode,
             "(5) +Controls +Strength, no FE", iv_var)
  )
}

# Both IVs use the same state-level control set (CTRL_MID == CTRL_FATAL).
headline_results <- bind_rows(
  run_ladder("log_onset_lag1",       CTRL_MID),
  run_ladder("log_fatal_onset_lag1", CTRL_FATAL)
)

# Filter analysis_sample to require fatal_mid_onset_count too — drop
# rows where fatal IV would be NA so both ladders are on the same panel
# IF we want apples-to-apples. (Currently we leave them on per-IV
# samples; both will print N for transparency.)

cat("\n========== HEADLINE RESULTS ==========\n")
print(headline_results, digits = 4, n = Inf)

# --- Diagnostic: is the log(1+strength) control itself significant? ---
# Reports the log_strength_t1 coefficient + p-value in the strength specs
# (isolation dummy kept [2026-10-12]), for both IVs and both FE / no-FE forms.
cat("\n--- log_strength_t1 coefficient (own significance) ---\n")
strength_diag <- bind_rows(
  fit_spec(churn ~ log_onset_lag1 + log_strength_t1 +
             democracy + autocracy + cinc +
             major_power + is_isolated_t1 | ccode + year,
           analysis_sample, ~ ccode + year, "(3) TWFE, MID", "log_strength_t1"),
  fit_spec(churn ~ log_fatal_onset_lag1 + log_strength_t1 +
             democracy + autocracy + cinc +
             major_power + is_isolated_t1 | ccode + year,
           analysis_sample, ~ ccode + year, "(3) TWFE, Fatal", "log_strength_t1"),
  fit_spec(churn ~ log_onset_lag1 + log_strength_t1 +
             democracy + autocracy + cinc +
             major_power + is_isolated_t1 + poly(year, 2),
           analysis_sample, ~ ccode, "(5) no FE, MID", "log_strength_t1"),
  fit_spec(churn ~ log_fatal_onset_lag1 + log_strength_t1 +
             democracy + autocracy + cinc +
             major_power + is_isolated_t1 + poly(year, 2),
           analysis_sample, ~ ccode, "(5) no FE, Fatal", "log_strength_t1")
) |>
  mutate(sig_05 = pval < 0.05)
print(strength_diag |> select(spec, estimate, se, pval, sig_05), digits = 4)

write_csv(headline_results, file.path(OUTPUT_DIR,
          "robustness_strength_control_results.csv"))

# Headline plot — focal IV only: fatal MID onset at t-1 (main-text Table 3).
plot_headline <- bind_rows(
  headline_results |> filter(iv == "log_fatal_onset_lag1")
) |>
  mutate(iv_label = factor(
    case_when(
      iv == "log_fatal_onset_lag1"      ~ "log(1+Fatal MID Onset), t-1",
      TRUE                              ~ iv
    ),
    levels = c("log(1+Fatal MID Onset), t-1")
  ))

p_headline <- ggplot(plot_headline,
                      aes(x = estimate, y = fct_rev(spec))) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray60") +
  geom_errorbarh(aes(xmin = ci_low, xmax = ci_high),
                 height = 0.2, color = ACCENT, linewidth = 0.6) +
  geom_point(size = 2.8, color = ACCENT) +
  facet_wrap(~ iv_label, ncol = 1, scales = "free_x") +
  labs(x = "Coefficient (95% CI)", y = NULL,
       title = "Conflict-IV effect with and without strength control",
       subtitle = "DV: churn (not logged). Strength = sum of edge weights at t-1.") +
  theme_pres(base_size = 11)
save_pres("fig_ms_strength_control.png", p_headline, width = 9, height = 6)

# ============================================================================
# 4. Multiverse: layer combos × weighting with strength control
# ============================================================================

cat("\n=== MULTIVERSE: Layer combos x Weighting with strength control ===\n")

all_layer_vars <- metadata$layers
all_layer_vars <- all_layer_vars[sapply(all_layer_vars, function(lv) {
  sd(dyad_data[[lv]], na.rm = TRUE) > 0
})]
MAIN_LAYERS <- c("layer_offensive", "layer_defensive", "layer_igo")

# Build layer combos (same as 03)
all_combos_raw <- list()
for (k in seq_along(all_layer_vars)) {
  all_combos_raw <- c(all_combos_raw, combn(all_layer_vars, k, simplify = FALSE))
}
extra_layers <- setdiff(all_layer_vars, MAIN_LAYERS)
if (length(extra_layers) > 0) {
  for (k in seq_along(extra_layers)) {
    for (ex in combn(extra_layers, k, simplify = FALSE)) {
      all_combos_raw <- c(all_combos_raw, list(sort(c(MAIN_LAYERS, ex))))
    }
  }
}
all_combos_raw <- unique(lapply(all_combos_raw, sort))
all_combos <- Filter(function(combo) {
  "layer_igo" %in% combo | "layer_defensive" %in% combo
}, all_combos_raw)

# [2026-10-12] Coverage windows copied from 07_robustness_multiverse.R. Each combo is estimated only
# on the years where all its layers are measured (trade begins 1870; DCA 1980-2010), and churn_alt
# in the first covered year is dropped because it compares with a pre-coverage year.
PANEL_END <- 2014L
LAYER_COVERAGE <- list(
  layer_offensive = c(1816L, PANEL_END),
  layer_defensive = c(1816L, PANEL_END),
  layer_igo       = c(1816L, PANEL_END),
  layer_trade     = c(1870L, PANEL_END),
  layer_dca       = c(1980L, 2010L)
)
get_valid_year_range <- function(layers) {
  rng <- lapply(layers, function(lv) {
    if (!is.null(LAYER_COVERAGE[[lv]])) return(LAYER_COVERAGE[[lv]])
    ay <- dyad_data$year[!is.na(dyad_data[[lv]]) & dyad_data[[lv]] > 0]
    c(min(ay), max(ay))
  })
  c(max(sapply(rng, `[`, 1)), min(sapply(rng, `[`, 2)))
}

weighting_methods <- c("contig_manhattan", "contig_euclidean",
                       "dist_manhattan", "dist_euclidean", "unweighted")

cat(sprintf("  Layer combos: %d | Weighting: %d | Total: %d\n",
            length(all_combos), length(weighting_methods),
            length(all_combos) * length(weighting_methods)))

# --- Local helpers ---
range01_local <- function(x) {
  r <- range(x, na.rm = TRUE)
  if (r[2] == r[1]) return(rep(0, length(x)))
  (x - r[1]) / (r[2] - r[1])
}

compute_weight_local <- function(data, layers, method) {
  mat <- as.matrix(data[, layers, drop = FALSE])
  mat[is.na(mat)] <- 0
  if (method == "contig_manhattan") {
    wcols <- intersect(paste0(layers, "_weighted"), names(data))
    if (length(wcols) == 0) return(rep(0, nrow(data)))
    mat_w <- as.matrix(data[, wcols, drop = FALSE]); mat_w[is.na(mat_w)] <- 0
    result <- rowSums(mat_w) / length(wcols)
  } else if (method == "contig_euclidean") {
    wcols <- intersect(paste0(layers, "_weighted"), names(data))
    if (length(wcols) == 0) return(rep(0, nrow(data)))
    mat_w <- as.matrix(data[, wcols, drop = FALSE]); mat_w[is.na(mat_w)] <- 0
    mat_w[mat_w < 0] <- 0
    result <- sqrt(rowSums(mat_w^2)) / length(wcols)
  } else if (method == "dist_manhattan") {
    dcols <- intersect(paste0(layers, "_dist_weighted"), names(data))
    if (length(dcols) == 0) return(rep(0, nrow(data)))
    mat_d <- as.matrix(data[, dcols, drop = FALSE]); mat_d[is.na(mat_d)] <- 0
    result <- rowSums(mat_d) / length(dcols)
  } else if (method == "dist_euclidean") {
    dcols <- intersect(paste0(layers, "_dist_weighted"), names(data))
    if (length(dcols) == 0) return(rep(0, nrow(data)))
    mat_d <- as.matrix(data[, dcols, drop = FALSE]); mat_d[is.na(mat_d)] <- 0
    mat_d[mat_d < 0] <- 0
    result <- sqrt(rowSums(mat_d^2)) / length(dcols)
  } else if (method == "unweighted") {
    result <- rowSums(mat)
  } else {
    return(rep(0, nrow(data)))
  }
  years <- data$year
  for (yr in unique(years)) {
    idx <- which(years == yr)
    result[idx] <- range01_local(result[idx])
  }
  result
}

compute_churn_strength_pwf <- function(data, w_vec,
                                         at_war_lookup,
                                         at_war_fatal_lookup) {
  # Computes churn, strength_t1, and BOTH severity-routed partner_war_frac
  # variants from the same edge weights, so all four are consistent
  # within a cell.
  #
  # at_war_lookup:       (ccode, year, at_war       = 1L) — any-MID
  # at_war_fatal_lookup: (ccode, year, at_war_fatal = 1L) — fatal-MID

  data$w_alt <- w_vec

  ties <- data |>
    select(ccode1, ccode2, year, w_alt) |>
    filter(w_alt > 0) |>
    bind_rows(
      data |> filter(w_alt > 0) |>
        select(ccode1 = ccode2, ccode2 = ccode1, year, w_alt)
    ) |>
    group_by(ccode1, ccode2, year) |>
    summarise(w = max(w_alt, na.rm = TRUE), .groups = "drop")

  # --- Strength at t-1: sum of edge weights ---
  strength <- ties |>
    group_by(ccode1, year) |>
    summarise(strength = sum(w, na.rm = TRUE), .groups = "drop") |>
    rename(ccode = ccode1) |>
    mutate(year_next = year + 1L) |>
    select(ccode, year = year_next, strength_t1 = strength)

  # --- Combo-specific PWF: any-MID and fatal-MID ---
  pwf <- ties |>
    mutate(year_next = year + 1L) |>
    left_join(at_war_lookup,
              by = c("ccode2" = "ccode", "year_next" = "year")) |>
    mutate(at_war = replace_na(at_war, 0L)) |>
    group_by(ccode1, year_next) |>
    summarise(pwf_combo = sum(at_war * w) / sum(w), .groups = "drop") |>
    rename(ccode = ccode1, year = year_next)

  pwf_fatal <- ties |>
    mutate(year_next = year + 1L) |>
    left_join(at_war_fatal_lookup,
              by = c("ccode2" = "ccode", "year_next" = "year")) |>
    mutate(at_war_fatal = replace_na(at_war_fatal, 0L)) |>
    group_by(ccode1, year_next) |>
    summarise(pwf_fatal_combo = sum(at_war_fatal * w) / sum(w),
              .groups = "drop") |>
    rename(ccode = ccode1, year = year_next)

  # --- Churn per (state, year) ---
  # [2026-10-12] Iterate over the full state-year universe (every state-year in the dyad data), so the
  # isolation branches below (both isolated -> 0; one side isolated -> 1) actually fire, matching
  # 07_robustness_multiverse.R and the canonical rule in 01_data_cleaning/02_compute_network_measures.R.
  universe <- dplyr::bind_rows(data |> select(ccode = ccode1, year), data |> select(ccode = ccode2, year)) |>
    dplyr::distinct(ccode, year) |>
    dplyr::arrange(ccode, year)
  states <- unique(universe$ccode)
  result_list <- vector("list", nrow(universe))
  idx <- 0L
  for (st in states) {
    st_ties <- ties |> filter(ccode1 == st)
    st_years <- universe$year[universe$ccode == st]
    if (length(st_years) < 2) next
    for (j in 2:length(st_years)) {
      yr_t  <- st_years[j]
      yr_t1 <- st_years[j - 1]
      if (yr_t - yr_t1 != 1) next
      w_t  <- st_ties |> filter(year == yr_t)  |> select(ccode2, w)
      w_t1 <- st_ties |> filter(year == yr_t1) |> select(ccode2, w)
      if (nrow(w_t) == 0 && nrow(w_t1) == 0) { ch <- 0
      } else if (nrow(w_t) == 0 || nrow(w_t1) == 0) { ch <- 1
      } else {
        w_both <- full_join(w_t, w_t1, by = "ccode2", suffix = c("_t", "_t1")) |>
          mutate(w_t = replace_na(w_t, 0), w_t1 = replace_na(w_t1, 0))
        sum_min <- sum(pmin(w_both$w_t, w_both$w_t1))
        sum_max <- sum(pmax(w_both$w_t, w_both$w_t1))
        ch <- if (sum_max > 0) 1 - sum_min / sum_max else NA_real_
      }
      idx <- idx + 1L
      result_list[[idx]] <- tibble(ccode = st, year = yr_t, churn_alt = ch)
    }
  }
  churn_df <- bind_rows(result_list[seq_len(idx)]) |>
    mutate(log_churn_alt = log1p(churn_alt))

  list(churn = churn_df, strength = strength,
       pwf = pwf, pwf_fatal = pwf_fatal)
}

# Pre-compute at-war lookups (used by all cells; severity-routed PWF)
at_war_lookup <- dyad_data |>
  filter(cowmidongoing == 1) |>
  select(ccode1, ccode2, year) |>
  pivot_longer(c(ccode1, ccode2), values_to = "ccode") |>
  distinct(ccode, year) |>
  mutate(at_war = 1L)

at_war_fatal_lookup <- dyad_data |>
  filter(cowmidongoing == 1, cowmidfatal == 1) |>
  select(ccode1, ccode2, year) |>
  pivot_longer(c(ccode1, ccode2), values_to = "ccode") |>
  distinct(ccode, year) |>
  mutate(at_war_fatal = 1L)

# --- Cell worker: returns up to 4 rows per cell (2 IVs × 2 specs) ---
run_mv_cell <- function(combo_i, method) {
  layers <- all_combos[[combo_i]]
  combo_label <- paste(gsub("layer_", "", layers), collapse = "+")

  w_vec <- tryCatch(compute_weight_local(dyad_data, layers, method),
                    error = function(e) NULL)
  if (is.null(w_vec) || all(w_vec == 0)) return(NULL)

  cs <- tryCatch(
    compute_churn_strength_pwf(dyad_data, w_vec,
                                at_war_lookup, at_war_fatal_lookup),
    error = function(e) { warning("25 cell ", combo_label, "/", method, ": ", conditionMessage(e)); NULL }
  )
  if (is.null(cs)) return(NULL)

  yr_rng <- get_valid_year_range(layers)
  # [2026-10-12] strength is LEFT-joined and zero-filled: cs$strength only has states with ties in t-1,
  # so an inner join dropped the first-partner-gain rows (churn_alt = 1) and the both-isolated rows.
  merged <- analysis_sample |>
    filter(year >= yr_rng[1], year <= yr_rng[2]) |>
    select(ccode, year,
           log_onset_lag1, log_fatal_onset_lag1,
           democracy, autocracy, cinc, major_power, is_isolated_t1) |>
    inner_join(cs$churn |> select(ccode, year, churn_alt),
               by = c("ccode", "year")) |>
    mutate(churn_alt = if_else(year == yr_rng[1], NA_real_, churn_alt)) |>
    filter(!is.na(churn_alt)) |>
    left_join(cs$strength, by = c("ccode", "year")) |>
    left_join(cs$pwf,       by = c("ccode", "year")) |>
    left_join(cs$pwf_fatal, by = c("ccode", "year")) |>
    mutate(strength_t1     = replace_na(strength_t1, 0),
           log_strength_t1 = log1p(strength_t1),
           pwf_combo       = replace_na(pwf_combo, 0),
           pwf_fatal_combo = replace_na(pwf_fatal_combo, 0))

  if (nrow(merged) < 100) return(NULL)

  # [2026-10-12] partner war exposure removed from the control set (H3 dropped from the paper).
  CTRL_COMBO_MID   <- "democracy + autocracy + cinc + major_power + is_isolated_t1"
  CTRL_COMBO_FATAL <- "democracy + autocracy + cinc + major_power + is_isolated_t1"

  fit_iv <- function(iv_var, ctrl) {
    fit_no <- tryCatch(
      feols(as.formula(paste0("churn_alt ~ ", iv_var, " + ", ctrl,
                               " | ccode + year")),
            data = merged, vcov = ~ ccode + year),
      error = function(e) NULL
    )
    fit_str <- tryCatch(
      feols(as.formula(paste0("churn_alt ~ ", iv_var,
                               " + log_strength_t1 + ", ctrl,
                               " | ccode + year")),
            data = merged, vcov = ~ ccode + year),
      error = function(e) NULL
    )
    rows <- list()
    for (fp in list(list(fit = fit_no,  label = "Without strength"),
                     list(fit = fit_str, label = "With strength"))) {
      f <- fp$fit
      if (is.null(f)) next
      ct <- coeftable(f)
      if (!iv_var %in% rownames(ct)) next
      rows[[length(rows) + 1]] <- tibble(
        layers   = combo_label,
        method   = method,
        iv       = iv_var,
        spec     = fp$label,
        estimate = ct[iv_var, "Estimate"],
        se       = ct[iv_var, "Std. Error"],
        pval     = ct[iv_var, "Pr(>|t|)"],
        n_obs    = nobs(f)
      )
    }
    bind_rows(rows)
  }

  bind_rows(
    fit_iv("log_onset_lag1",       CTRL_COMBO_MID),
    fit_iv("log_fatal_onset_lag1", CTRL_COMBO_FATAL)
  )
}

# --- Parallel dispatch (mirrors 03's pattern) ---
spec_grid_mv <- expand.grid(
  combo_i = seq_along(all_combos),
  method  = weighting_methods,
  stringsAsFactors = FALSE
)
total_cells <- nrow(spec_grid_mv)
cat(sprintf("  Total cells: %d (%d combos x %d methods)\n",
            total_cells, length(all_combos), length(weighting_methods)))

.use_parallel <- requireNamespace("furrr", quietly = TRUE) &&
                 requireNamespace("future", quietly = TRUE)

# Defensive parallel setup:
#   - Hard-cap workers at 6 to avoid memory exhaustion (each multisession
#     worker gets a copy of dyad_data + analysis_sample + lookups, so 19
#     workers blew through RAM and crashed with FutureError on the
#     previous run).
#   - Bump future.globals.maxSize to 8GB so large globals can be
#     transferred to workers without hitting the 500MB default cap.
#   - Wrap the parallel dispatch in tryCatch so a worker crash falls
#     back to sequential rather than killing the whole script.
run_multiverse <- function(use_parallel) {
  if (use_parallel) {
    .n_workers_cap <- 6L
    .n_workers_env <- as.integer(Sys.getenv("ROBUSTNESS_WORKERS", as.character(.n_workers_cap)))
    .n_workers <- min(.n_workers_env, .n_workers_cap,
                      max(1L, future::availableCores() - 1L))
    options(future.globals.maxSize = 8 * 1024^3)  # 8 GB
    future::plan(future::multisession, workers = .n_workers)
    cat(sprintf("  Parallel: %d workers (capped to avoid memory pressure)\n",
                future::nbrOfWorkers()))
    on.exit(future::plan(future::sequential), add = TRUE)
    furrr::future_map2_dfr(
      spec_grid_mv$combo_i, spec_grid_mv$method,
      run_mv_cell,
      .options = furrr::furrr_options(
        seed = TRUE,
        # [2026-10-12] with an explicit globals list the workers did not attach dplyr/tidyr/fixest, every cell
        # failed inside a tryCatch and the multiverse silently returned 0 rows. Attach them explicitly.
        packages = c("dplyr", "tidyr", "fixest"),
        globals = c("all_combos", "dyad_data",
                    "at_war_lookup", "at_war_fatal_lookup",
                    "analysis_sample",
                    "compute_weight_local",
                    "compute_churn_strength_pwf",
                    "range01_local",
                    "LAYER_COVERAGE", "get_valid_year_range")
      )
    )
  } else {
    cat("  Sequential (furrr/future not installed)\n")
    purrr::map2_dfr(
      spec_grid_mv$combo_i, spec_grid_mv$method,
      run_mv_cell
    )
  }
}

mv_df <- tryCatch(
  run_multiverse(.use_parallel),
  error = function(e) {
    cat(sprintf("  Parallel run failed: %s\n", conditionMessage(e)))
    cat("  Falling back to sequential execution...\n")
    future::plan(future::sequential)
    purrr::map2_dfr(
      spec_grid_mv$combo_i, spec_grid_mv$method,
      run_mv_cell
    )
  }
)

cat(sprintf("\n  Multiverse complete: %d rows\n", nrow(mv_df)))
# [2026-10-12] fail loudly: an empty multiverse means every cell errored (see the warnings above)
if (nrow(mv_df) == 0) stop("25: strength multiverse returned 0 rows; every cell failed")

if (nrow(mv_df) > 0) {
  # Attenuation from adding strength — computed PER IV so the spread
  # for any-MID and fatal-MID can be reported separately.
  mv_wide <- mv_df |>
    select(iv, layers, method, spec, estimate) |>
    pivot_wider(names_from = spec, values_from = estimate) |>
    mutate(change_pct = 100 * (`With strength` - `Without strength`) /
             abs(`Without strength`))

  for (iv_lbl in unique(mv_wide$iv)) {
    sub <- mv_wide |> filter(iv == iv_lbl)
    cat(sprintf("\n--- Strength control: change in %s coefficient ---\n", iv_lbl))
    cat(sprintf("  Mean change:   %.1f%%\n", mean(sub$change_pct, na.rm = TRUE)))
    cat(sprintf("  Median change: %.1f%%\n", median(sub$change_pct, na.rm = TRUE)))
    cat(sprintf("  Range: %.1f%% to %.1f%%\n",
                min(sub$change_pct, na.rm = TRUE),
                max(sub$change_pct, na.rm = TRUE)))

    str_rows <- mv_df |> filter(iv == iv_lbl, spec == "With strength")
    cat(sprintf("  With strength: %d / %d sig at p < 0.05 (%.1f%%)\n",
                sum(str_rows$pval < 0.05, na.rm = TRUE), nrow(str_rows),
                100 * mean(str_rows$pval < 0.05, na.rm = TRUE)))
  }

  write_csv(mv_df, file.path(OUTPUT_DIR,
            "robustness_strength_control_multiverse.csv"))
  write_csv(mv_wide, file.path(OUTPUT_DIR,
            "robustness_strength_control_change.csv"))

  # Robustness check visualization: boxplot of coefficients +
  #   stacked barplot of share positive & sig at one-sided 0.05.
  #   Replaces the % change tile, which was unstable when the
  #   without-strength coefficient was near zero. One-sided test is
  #   appropriate: the multiverse substantiates the main MID-onset > 0
  #   result, so we ask whether each cell agrees with that direction.
  iv_lbl <- c(log_onset_lag1       = "log(1+MID Onset), t-1",
              log_fatal_onset_lag1 = "log(1+Fatal MID Onset), t-1")

  mv_box_df <- mv_df |>
    mutate(spec     = factor(spec,
                             levels = c("Without strength", "With strength")),
           iv_label = factor(iv_lbl[iv],
                             levels = c("log(1+MID Onset), t-1",
                                        "log(1+Fatal MID Onset), t-1")))

  mv_class_df <- mv_box_df |>
    mutate(
      pos_sig_1side = (estimate > 0) & ((pval / 2) < 0.05),
      cat = factor(if_else(pos_sig_1side,
                           "Positive & sig (1-sided p<0.05)",
                           "Not positive & sig"),
                   levels = c("Positive & sig (1-sided p<0.05)",
                              "Not positive & sig"))
    ) |>
    count(iv_label, spec, cat, .drop = FALSE) |>
    group_by(iv_label, spec) |>
    mutate(share = n / sum(n)) |>
    ungroup()

  n_per_variant <- mv_box_df |>
    count(iv_label, spec) |>
    summarise(n = max(n)) |>
    pull(n)

  p_box <- ggplot(mv_box_df,
                  aes(x = spec, y = estimate, fill = spec)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray60") +
    geom_boxplot(width = 0.5, outlier.size = 0.6,
                 alpha = 0.75, linewidth = 0.4) +
    facet_wrap(~ iv_label, ncol = 2, scales = "free_y") +
    scale_fill_manual(values = c("Without strength" = ACCENT,
                                  "With strength"    = ACCENT2),
                      guide = "none") +
    scale_x_discrete(labels = c("Without strength" = "Without\nstrength",
                                 "With strength"    = "With\nstrength")) +
    labs(x = NULL, y = "Coefficient on conflict IV",
         title = "(a) Coefficient distribution across multiverse cells") +
    theme_pres(base_size = 11)

  p_bar <- ggplot(mv_class_df,
                  aes(x = spec, y = share, fill = cat)) +
    geom_col(width = 0.5) +
    geom_text(aes(label = ifelse(share > 0.04,
                                  scales::percent(share, accuracy = 1),
                                  "")),
              position = position_stack(vjust = 0.5),
              color = "white", size = 3.4, fontface = "bold") +
    facet_wrap(~ iv_label, ncol = 2) +
    scale_fill_manual(values = c("Positive & sig (1-sided p<0.05)" = "#C0392B",
                                  "Not positive & sig"             = "#2E75B6"),
                      name = NULL) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                       expand = expansion(mult = c(0, 0.02))) +
    scale_x_discrete(labels = c("Without strength" = "Without\nstrength",
                                 "With strength"    = "With\nstrength")) +
    labs(x = NULL, y = "Share of multiverse cells",
         title = "(b) Positive & significant at one-sided p<0.05") +
    theme_pres(base_size = 11) +
    theme(legend.position = "bottom",
          legend.key.size = unit(0.4, "cm"))

  p_combined <- p_box / p_bar +
    plot_annotation(
      title    = "Strength-control multiverse: with vs without log1p(strength_t1)",
      subtitle = sprintf(
        "%d cells per variant per IV. Boxplot shows coefficient spread; bars show 1-sided positive-significance share.",
        n_per_variant)
    )

  save_pres("fig_ms_strength_multiverse.png", p_combined,
            width = 11, height = 10)
}

cat("\nDone. 25_robustness_strength_control.R complete.\n")
