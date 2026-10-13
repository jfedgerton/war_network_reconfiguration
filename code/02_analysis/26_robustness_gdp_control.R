###############################################################################
# 22_robustness_gdp_control.R
# Robustness: replicate the main-text Table 3 state-level ladder with
# log GDP added as an extra control. Tests whether the H2 (own MID)
# effects survive when we condition on
# economic capacity in addition to CINC.
#
# Data: GDP (already in natural log) from Anders/Fariss/Markowitz extended
# series, accessed through peacesciencer::add_sdp_gdp(). Coverage is
# 1816-2014, matching the panel.
#
# Specifications mirror Table 3 columns 3-4 (the controlled regressions):
#   (3) +Ctrls, no FE  +log GDP
#   (4) TWFE           +log GDP
# (Column 5 from the main-text ladder — TWFE + community-clustered SE — is
# omitted here. This script is a control-set robustness check, not an SE
# robustness check, so requiring the multilayer Leiden community label
# would just add overhead without changing the substantive comparison.)
# Each is also re-run WITHOUT log GDP on the same restricted sample
# (state-years where GDP is non-missing) so the coefficient comparison
# is apples-to-apples.
#
# Inputs:  data/processed/state_network_churn.rds
#          data/processed/dyad_year_network_layers.rds
# Outputs:
#   data/processed/robustness_gdp_control.csv
#   results/appendix/tableA_gdp_robustness.tex
#   code/03_figures/output/manuscript/fig_ms_gdp_robustness.png
###############################################################################

library(tidyverse)
library(fixest)
library(peacesciencer)
library(texreg)
select <- dplyr::select

set.seed(123)

OUTPUT_DIR <- here::here("data", "processed")
APPX_DIR   <- here::here("results", "appendix")
FIG_DIR    <- here::here("code", "03_figures", "output", "manuscript")
dir.create(APPX_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG_DIR,  recursive = TRUE, showWarnings = FALSE)

# ----------------------------------------------------------------------------
# Manuscript theme (matches 13_robustness_hostlev pattern)
# ----------------------------------------------------------------------------
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
      strip.text       = element_text(face = "bold", size = base_size,
                                       color = "black"),
      legend.position  = "bottom",
      legend.text      = element_text(size = base_size - 1),
      legend.title     = element_text(size = base_size, face = "bold"),
      plot.margin      = margin(10, 12, 8, 10)
    )
}

save_pres <- function(filename, plot, width = 8, height = 6, dpi = 400) {
  ggsave(file.path(FIG_DIR, filename), plot,
         width = width, height = height, dpi = dpi, bg = "white")
  cat(sprintf("  -> %s\n", filename))
}

# ============================================================================
# 1. Load data
# ============================================================================

cat("Loading state-level churn data...\n")
state_churn <- readRDS(file.path(OUTPUT_DIR, "state_network_churn.rds"))
dyad_data   <- readRDS(file.path(OUTPUT_DIR, "dyad_year_network_layers.rds"))

# State attributes already merged from 02_state_level_analysis.R aren't
# saved alongside churn — re-merge the controls we need (CINC, major_power,
# polity, partner-war exposure, is_isolated_t1) using the same logic.
# This mirrors the prep block in 13_robustness_hostlev.R.

state_attrs <- dyad_data |>
  select(ccode = ccode1, year, cinc = cinc1, major_power = cowmaj1) |>
  distinct(ccode, year, .keep_all = TRUE) |>
  mutate(major_power = as.integer(coalesce(major_power, 0L) == 1))

state_churn <- state_churn |>
  left_join(state_attrs, by = c("ccode", "year"))

# Polity (democracy / autocracy)
cat("Adding polity2 democracy / autocracy...\n")
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

# ============================================================================
# 2. Add log GDP via peacesciencer::add_sim_gdp_pop()
# ============================================================================
# add_sim_gdp_pop() (replaces deprecated add_sdp_gdp) returns:
#   mrgdppc   — Markowitz/Fariss simulated real GDP per capita (raw units)
#   pwtrgdp   — PWT real GDP (raw units)
#   pwtpop    — PWT population (raw units)
#   sd_*      — posterior SDs of the simulated quantities
# We use pwtrgdp (total real GDP, not per-capita) and log it so it
# parallels the headline "log GDP" interpretation. Use log1p() to be
# defensive against any near-zero values.

cat("Adding log GDP via peacesciencer::add_sim_gdp_pop()...\n")

state_years_gdp <- create_stateyears(
  subset_years = min(state_churn$year):max(state_churn$year))
attr(state_years_gdp, "ps_data_type") <- "state_year"
attr(state_years_gdp, "ps_system")    <- "cow"
state_years_gdp <- state_years_gdp |> add_sim_gdp_pop()

# Keep just (ccode, year, log_gdp).
gdp_panel <- state_years_gdp |>
  select(ccode, year, pwtrgdp) |>
  mutate(log_gdp = log1p(pwtrgdp)) |>
  select(ccode, year, log_gdp)

state_churn <- state_churn |>
  left_join(gdp_panel, by = c("ccode", "year"))

cat(sprintf("  log GDP coverage: %d / %d state-years (%.1f%%)\n",
            sum(!is.na(state_churn$log_gdp)),
            nrow(state_churn),
            100 * mean(!is.na(state_churn$log_gdp))))
cat(sprintf("  Year range with GDP: %d - %d\n",
            min(state_churn$year[!is.na(state_churn$log_gdp)]),
            max(state_churn$year[!is.na(state_churn$log_gdp)])))

# ============================================================================
# 3. Compute partner-war exposure (severity-routed, weighted)
# ============================================================================
# Same logic as 02_state_level_analysis.R, repeated here so the script is
# stand-alone. main_weight = alliance + IGO composite.

if (!"main_weight" %in% names(dyad_data)) {
  stop("dyad_data is missing main_weight. Run 02_compute_network_measures.R first.")
}

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
  filter(cowmidfatal == 1) |>
  select(ccode1, ccode2, year) |>
  pivot_longer(c(ccode1, ccode2), values_to = "ccode") |>
  distinct(ccode, year) |>
  mutate(at_war_fatal = 1L)

pwe <- partner_sets |>
  mutate(year_next = year + 1L) |>
  left_join(at_war_dyad,       by = c("ccode2" = "ccode", "year_next" = "year")) |>
  left_join(at_war_dyad_fatal, by = c("ccode2" = "ccode", "year_next" = "year")) |>
  mutate(at_war = replace_na(at_war, 0L),
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

# Final harmonization: log_churn, log_mid_onset, etc.
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

# ============================================================================
# 4. Restrict to GDP-non-missing sample and run apples-to-apples comparisons
# ============================================================================
# Run each specification TWICE on the same sample:
#   (a) baseline controls (no GDP)
#   (b) baseline + log GDP
# Sample is restricted to state-years where log_gdp is non-missing so the
# coefficient comparison reflects the marginal effect of adding GDP, not
# a sample composition shift.

state_gdp_sample <- state_churn |>
  filter(!is.na(log_gdp))

cat(sprintf("\n  GDP-non-missing sample: %s state-years (%.1f%% of full)\n",
            format(nrow(state_gdp_sample), big.mark = ","),
            100 * nrow(state_gdp_sample) / nrow(state_churn)))

# Control sets (matches 02 + downstream conventions); MID and fatal-MID
# sets are now identical.
# [2026-10-12] partner war exposure removed from the control set (H3 dropped from the paper).
CTRL_MID_NO_GDP   <- "democracy + autocracy + cinc + major_power + is_isolated_t1"
CTRL_MID_GDP      <- paste(CTRL_MID_NO_GDP, "+ log_gdp")
CTRL_FATAL_NO_GDP <- "democracy + autocracy + cinc + major_power + is_isolated_t1"
CTRL_FATAL_GDP    <- paste(CTRL_FATAL_NO_GDP, "+ log_gdp")

# Pretty-label lookup for every coefficient we extract
COEF_LABELS <- c(
  "log_onset_lag1"                  = "log(1+MID Onset), t-1",
  "log_fatal_onset_lag1"            = "log(1+Fatal MID Onset), t-1"
)

# Helper: fit one spec, return one tidy row per coefficient in `coefs_to_extract`
fit_one <- function(formula, data, vcov_arg, label,
                    coefs_to_extract, control_label, model_label) {
  fit <- tryCatch(
    feols(formula, data = data, vcov = vcov_arg),
    error = function(e) {
      cat(sprintf("    ! %s failed: %s\n", label, conditionMessage(e)))
      return(NULL)
    }
  )
  if (is.null(fit)) return(NULL)
  ct <- coeftable(fit)
  rows <- lapply(coefs_to_extract, function(v) {
    if (!v %in% rownames(ct)) return(NULL)
    tibble(
      iv       = v,
      iv_label = unname(COEF_LABELS[v]),
      spec     = model_label,
      control  = control_label,
      estimate = ct[v, "Estimate"],
      se       = ct[v, "Std. Error"],
      pval     = ct[v, "Pr(>|t|)"],
      n_obs    = nobs(fit),
      ci_low   = ct[v, "Estimate"] - 1.96 * ct[v, "Std. Error"],
      ci_high  = ct[v, "Estimate"] + 1.96 * ct[v, "Std. Error"]
    )
  })
  bind_rows(rows)
}

# Build the regressions: 2 IVs × 2 specs × 2 control sets = 8 fits.
# Each fit produces 1 tidy row (focal IV),
# so results_df is 8 rows total.
# state-clustered SE for col 3 (no FE), 2-way (state+year) for col 4 —
# matches Table 3 conventions for cols 3 and 4.

results <- list()

for (iv_pair in list(
  list(iv  = "log_onset_lag1",
       ctrl_no  = CTRL_MID_NO_GDP,   ctrl_yes = CTRL_MID_GDP),
  list(iv  = "log_fatal_onset_lag1",
       ctrl_no  = CTRL_FATAL_NO_GDP, ctrl_yes = CTRL_FATAL_GDP)
)) {
  iv  <- iv_pair$iv
  # Each fit returns the row for the focal IV (H2).
  coefs_to_extract <- c(iv)

  for (gdp_variant in list(
    list(label = "No GDP",  ctrl = iv_pair$ctrl_no),
    list(label = "+log GDP", ctrl = iv_pair$ctrl_yes)
  )) {
    cl <- gdp_variant$label
    ctrl <- gdp_variant$ctrl

    # Spec 3: full controls + poly(year, 2), no FE
    fml_3 <- as.formula(paste0("churn ~ ", iv, " + ", ctrl, " + poly(year, 2)"))
    results[[length(results) + 1]] <- fit_one(
      fml_3, state_gdp_sample, ~ ccode,
      sprintf("(%s, %s) +Ctrls,noFE", iv, cl),
      coefs_to_extract, cl, "(2) +Ctrls, no FE")

    # Spec 4: TWFE
    fml_4 <- as.formula(paste0("churn ~ ", iv, " + ", ctrl, " | ccode + year"))
    results[[length(results) + 1]] <- fit_one(
      fml_4, state_gdp_sample, ~ ccode + year,
      sprintf("(%s, %s) TWFE", iv, cl),
      coefs_to_extract, cl, "(3) TWFE")
  }
}

results_df <- bind_rows(results)
print(results_df, n = Inf, digits = 4)

write_csv(results_df, file.path(OUTPUT_DIR, "robustness_gdp_control.csv"))
cat(sprintf("\nSaved: data/processed/robustness_gdp_control.csv (%d rows)\n",
            nrow(results_df)))

# ============================================================================
# 5. Coefficient-comparison plot
# ============================================================================

# [2026-10-12] focal IV only: fatal MID onset at t-1 (main-text Table 3).
plot_df <- results_df |>
  filter(iv_label %in% c("log(1+Fatal MID Onset), t-1")) |>
  mutate(
    spec     = factor(spec, levels = c("(2) +Ctrls, no FE",
                                        "(3) TWFE")),
    control  = factor(control, levels = c("No GDP", "+log GDP")),
    iv_label = factor(iv_label,
                      levels = c("log(1+Fatal MID Onset), t-1"))
  )

p_gdp <- ggplot(plot_df,
                aes(x = estimate, y = spec, color = control, shape = control)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray60") +
  geom_errorbarh(aes(xmin = ci_low, xmax = ci_high),
                 height = 0.18, position = position_dodge(width = 0.5),
                 linewidth = 0.6) +
  geom_point(size = 2.6, position = position_dodge(width = 0.5)) +
  facet_wrap(~ iv_label, scales = "free_x", ncol = 2,
             labeller = label_wrap_gen(width = 24)) +
  scale_color_manual(values = c("No GDP" = ACCENT, "+log GDP" = ACCENT2),
                     name = NULL) +
  scale_shape_manual(values = c("No GDP" = 16, "+log GDP" = 17),
                     name = NULL) +
  labs(x = "Coefficient (95% CI)", y = NULL,
       title = "GDP-control robustness (H2 IV)") +
  theme_pres(base_size = 11) +
  theme(panel.spacing = unit(1.2, "lines"),
        strip.text = element_text(size = 10))

save_pres("fig_ms_gdp_robustness.png", p_gdp, width = 9, height = 4)

# ============================================================================
# 6. LaTeX summary table → results/appendix
# ============================================================================

tex <- c(
  "\\begin{table}[htbp]", "\\centering",
  "\\caption{State-level Table~\\ref{tab:state_main} ladder (columns 2--3) with log GDP added as a control. ",
  "Reports the conflict IV ($\\log(1+\\text{MID onset})_{t-1}$ and its fatal-MID variant, the headline IV). ",
  "Same restricted sample (state-years with non-missing GDP) ",
  "is used in both the no-GDP and +log GDP columns to keep the coefficient ",
  "comparison apples-to-apples.}",
  "\\label{tab:gdp_robustness}",
  "\\small",
  "\\begin{tabular}{lllrrrr}", "\\toprule",
  "Coefficient & Spec & Control set & $\\hat{\\beta}$ & SE & $p$ & N \\\\",
  "\\midrule"
)
for (i in seq_len(nrow(results_df))) {
  r <- results_df[i, ]
  tex <- c(tex, sprintf(
    "%s & %s & %s & %.4f & %.4f & %.3f & %s \\\\",
    r$iv_label, r$spec, r$control,
    r$estimate, r$se, r$pval,
    format(r$n_obs, big.mark = ",")
  ))
}
tex <- c(tex, "\\bottomrule", "\\end{tabular}",
  "\\begin{tablenotes}\\small",
  "\\item DV: churn (not logged). Each row is the conflict coefficient from a ",
  "single fit, for $\\log(1+\\text{MID onset})_{t-1}$ or its fatal-MID variant ",
  "(the headline IV). Specs (2)--(3) replicate columns 2--3 of ",
  "Table~\\ref{tab:state_main} with log GDP added as an extra control. ",
  "Log GDP is $\\log(1+\\text{PWT real GDP})$ from the ",
  "Anders/Fariss/Markowitz extended series, accessed via ",
  "\\texttt{peacesciencer::add\\_sim\\_gdp\\_pop()}. The same restricted ",
  "sample (state-years with non-missing GDP) is used in both the no-GDP ",
  "and +log GDP columns. Spec (2) clusters SEs by state; (3) by state $+$ ",
  "year. The community-clustered SE spec (column 4 in the main-text table) ",
  "is not reproduced here since this is a control-set robustness check, ",
  "not a SE-clustering robustness check.",
  "\\end{tablenotes}",
  "\\end{table}")
writeLines(tex, file.path(APPX_DIR, "tableA_gdp_robustness.tex"))
cat(sprintf("  LaTeX -> %s\n",
            file.path("results", "appendix", "tableA_gdp_robustness.tex")))

cat("\nDone. 22_robustness_gdp_control.R complete.\n")
