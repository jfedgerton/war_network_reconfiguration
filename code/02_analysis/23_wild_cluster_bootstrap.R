###############################################################################
# 19_wild_cluster_bootstrap.R
# Robustness: wild-cluster bootstrap p-values for the state-level TWFE
# main models.
#
# The main-text TWFE uses Cameron-Gelbach-Miller two-way clustered
# standard errors (~ ccode + year). Asymptotic two-way clustered
# inference can over-reject when cluster sizes are unequal -- and
# they are in this panel: state cluster sizes range from a handful of
# years for short-lived polities (Republic of Texas, post-WWII new
# states observed late) to ~200 years for the major powers.
#
# Wild-cluster bootstrap (Cameron, Gelbach & Miller 2008;
# MacKinnon & Webb 2018, "The wild bootstrap for few (treated) clusters")
# resamples residuals with cluster-level Rademacher weights. It is
# robust to few / unbalanced clusters and gives valid inference when
# the asymptotic clustered distribution is poorly approximated.
#
# Implemented via fwildclusterboot::boottest, which works directly on
# fixest::feols objects. Two clustering schemes are reported per IV:
#   one-way ~ ccode    (reasonable when year-FE absorbs time shocks)
#   two-way ~ ccode + year   (matches the main-text spec)
#
# Outputs:
#   data/processed/robustness_wild_bootstrap.csv
#   code/03_figures/output/manuscript/fig_ms_wild_bootstrap.png
#
# Reading:
#   Cameron, Gelbach & Miller (2008) JBES; MacKinnon & Webb (2018) JAE;
#   Roodman, Nielsen, MacKinnon & Webb (2019) Stata Journal.
###############################################################################

library(tidyverse)
library(fixest)
library(peacesciencer)
if (!requireNamespace("fwildclusterboot", quietly = TRUE)) {
  stop("Install fwildclusterboot: install.packages('fwildclusterboot')")
}
library(fwildclusterboot)
select <- dplyr::select
set.seed(123)

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
# Load data + state-level controls (mirrors scripts 11/13/14/...)
# ============================================================

cat("Loading inputs...\n")

state_churn <- readRDS(file.path(OUTPUT_DIR, "state_network_churn.rds"))
dyad_data   <- readRDS(file.path(OUTPUT_DIR, "dyad_year_network_layers.rds"))

state_controls_needed <- c("partner_war_frac", "democracy", "autocracy",
                            "cinc", "major_power")
if (!all(state_controls_needed %in% names(state_churn))) {
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
  state_churn <- state_churn |>
    left_join(state_attrs, by = c("ccode", "year"))

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

state_churn <- state_churn |> mutate(log_churn = log1p(churn))

state_data <- state_churn |>
  mutate(log_mid_onset_count       = log1p(mid_onset_count),
         log_fatal_mid_onset_count = log1p(fatal_mid_onset_count)) |>
  filter(!is.na(churn))

cat(sprintf("  State-year obs: %d, states: %d, years: %d-%d\n",
            nrow(state_data),
            n_distinct(state_data$ccode),
            min(state_data$year), max(state_data$year)))

# ============================================================
# State-level IVs (match scripts 02 / 13)
# ============================================================

# Focal IVs only -- avoid running the wild bootstrap (B = 9999) over the
# full conflict-measure set. Figure shows just these two.
state_ivs <- list(
  list(col = "log_mid_onset_count", label = "log(1 + MID onset)"),
  # Weighted partner-war exposure as focal IV: keep the conflict IV as a
  # control but drop partner_war_frac_weighted from the control set (it is
  # the IV here) to match the main H3 specification.
  list(col = "partner_war_frac_weighted", label = "Partner War Exposure",
       ctrl = "log_mid_onset_count + democracy + autocracy + cinc + major_power + is_isolated_t1")
)

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
# state_ctrl_str is now routed per IV via ctrl_for_iv() inside run_one_iv().

# ============================================================
# Helper: asymptotic clustered SE/p (matches main spec) +
# wild-cluster bootstrap p
# ============================================================

B_REPS <- 9999L

run_one_iv <- function(iv) {
  iv_col <- iv$col; iv_lab <- iv$label
  # Allow a per-IV control override (e.g. partner-war exposure as the focal
  # IV needs a control set that excludes itself); else severity-route.
  ctrl <- if (!is.null(iv$ctrl)) iv$ctrl else ctrl_for_iv(iv_col)
  cat(sprintf("\n[IV] %s\n", iv_col))

  f <- as.formula(paste0("log_churn ~ ", iv_col, " + ", ctrl,
                          " | ccode + year"))
  m <- feols(f, data = state_data, cluster = ~ ccode + year)
  cf <- coef(m)[iv_col]
  vv <- vcov(m, cluster = ~ ccode + year)
  asy_se <- sqrt(diag(vv)[iv_col])
  asy_t  <- cf / asy_se
  asy_p  <- 2 * pnorm(-abs(asy_t))

  vv2 <- vcov(m, cluster = ~ ccode)
  asy_se_1way <- sqrt(diag(vv2)[iv_col])
  asy_p_1way  <- 2 * pnorm(-abs(cf / asy_se_1way))

  # Wild-cluster bootstrap, one-way (~ ccode)
  cat("  wild bootstrap (one-way: ccode)...\n")
  bt1 <- tryCatch(
    boottest(m, param = iv_col, B = B_REPS,
             clustid = "ccode",
             type = "rademacher", impose_null = TRUE),
    error = function(e) { cat("    failed: ", conditionMessage(e), "\n"); NULL }
  )

  # Wild-cluster bootstrap, two-way (~ ccode + year)
  cat("  wild bootstrap (two-way: ccode + year)...\n")
  bt2 <- tryCatch(
    boottest(m, param = iv_col, B = B_REPS,
             clustid = c("ccode", "year"),
             type = "rademacher", impose_null = TRUE),
    error = function(e) { cat("    failed: ", conditionMessage(e), "\n"); NULL }
  )

  tibble(
    iv             = iv_col,
    iv_label       = iv_lab,
    estimate       = unname(cf),
    se_asy_1way    = asy_se_1way,
    p_asy_1way     = asy_p_1way,
    se_asy_2way    = unname(asy_se),
    p_asy_2way     = unname(asy_p),
    p_wcb_1way     = if (!is.null(bt1)) unname(bt1$p_val) else NA_real_,
    ci_low_wcb_1way  = if (!is.null(bt1)) unname(bt1$conf_int[1]) else NA_real_,
    ci_high_wcb_1way = if (!is.null(bt1)) unname(bt1$conf_int[2]) else NA_real_,
    p_wcb_2way     = if (!is.null(bt2)) unname(bt2$p_val) else NA_real_,
    ci_low_wcb_2way  = if (!is.null(bt2)) unname(bt2$conf_int[1]) else NA_real_,
    ci_high_wcb_2way = if (!is.null(bt2)) unname(bt2$conf_int[2]) else NA_real_,
    n_obs          = nobs(m),
    n_clusters_state = n_distinct(state_data$ccode),
    n_clusters_year  = n_distinct(state_data$year),
    B              = B_REPS
  )
}

results <- purrr::map_dfr(state_ivs, run_one_iv)

write_csv(results,
          file.path(OUTPUT_DIR, "robustness_wild_bootstrap.csv"))
cat(sprintf("\nSaved: robustness_wild_bootstrap.csv (%d rows)\n",
            nrow(results)))

# ============================================================
# Plot: side-by-side comparison of asymptotic vs WCB p-values
# ============================================================

cat("\n--- Plotting ---\n")

# Levels and labels supplied inline so the script is robust against
# a polluted workspace where iv_levels_all may have been redefined
# as the 8-element sys+state set by another script in the session.
iv_levels_state <- c("log(1 + MID onset)", "Partner War Exposure")
iv_facet_labels <- c("log(1 +\nMID onset)", "Partner War\nExposure")

p_dat <- results |>
  filter(iv_label %in% iv_levels_state) |>
  select(iv_label, p_asy_1way, p_asy_2way, p_wcb_1way, p_wcb_2way) |>
  pivot_longer(-iv_label, names_to = "method", values_to = "pval") |>
  mutate(
    cluster_scheme = case_when(
      str_detect(method, "1way") ~ "Cluster: state",
      str_detect(method, "2way") ~ "Cluster: state + year"),
    inference = case_when(
      str_detect(method, "asy") ~ "Asymptotic",
      str_detect(method, "wcb") ~ "Wild-cluster bootstrap"),
    iv_label = factor(iv_label,
                       levels = iv_levels_state,
                       labels = iv_facet_labels)
  )

panel_theme <- theme_pres() +
  theme(panel.spacing = unit(1.2, "lines"),
        panel.spacing.x = unit(1.5, "lines"),
        strip.placement = "outside",
        strip.text = element_text(color = "black", size = 9, face = "bold",
                                    margin = margin(4, 4, 4, 4)),
        legend.position = "bottom",
        plot.margin = margin(t = 10, r = 12, b = 14, l = 10))

# State-level (churn) only -- this script tests TWFE clustering at the state
# panel; no system-level analogue. Naming kept consistent with the
# split-by-outcome convention used in scripts 16-23.
p_pvals <- p_dat |>
  ggplot(aes(x = iv_label, y = pval/2,
             color = inference, shape = inference, group = inference)) +
  geom_hline(yintercept = 0.05, color = "black",
             linetype = "dashed", linewidth = 0.8) +
  geom_hline(yintercept = 0.01, color = "black",
             linetype = "dotted", linewidth = 0.8) +
  geom_point(size = 3) +
  scale_color_manual(values = c("Asymptotic" = ACCENT2,
                                "Wild-cluster bootstrap" = ACCENT)) +
  scale_y_continuous(limits = c(0, NA)) +
  facet_wrap(~ cluster_scheme, ncol = 2,
             labeller = label_wrap_gen(width = 22)) +
  labs(
    title = "Wild-cluster bootstrap",
    x = NULL, y = "p-value (two-sided)",
    color = NULL, shape = NULL
  ) +
  panel_theme +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 1))

ggsave(file.path(FIG_DIR, "fig_ms_wild_bootstrap_churn.png"),
       p_pvals, width = 10, height = 6, dpi = 400, bg = "white")
