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

# Focal IV only -- avoid running the wild bootstrap (B = 9999) over the
# full conflict-measure set. [2026-10-12] The partner-war-exposure (H3) IV was removed.
# [2026-10-12] focal IV follows the fatal t-1 headline (main-text Table 3, column 3).
state_ivs <- list(
  list(col = "log_fatal_onset_lag1", label = "log(1 + Fatal MID onset) (t-1)")
)

# [2026-10-12] partner war exposure removed from the control set (H3 dropped from the paper).
# Same control set as main-text Table 3 column 3 (state + year FE, SEs clustered by state and year).
CTRL_MID   <- "democracy + autocracy + cinc + major_power + is_isolated_t1"
CTRL_FATAL <- "democracy + autocracy + cinc + major_power + is_isolated_t1"
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

  f <- as.formula(paste0("churn ~ ", iv_col, " + ", ctrl,
                          " | ccode + year"))
  # [2026-10-12] defensive: drop any row where the t-1 IV is undefined (state not in the system in t-1)
  # up front so boottest() sees exactly the estimation sample.
  d_iv <- state_data |> filter(!is.na(.data[[iv_col]]))
  m <- feols(f, data = d_iv, cluster = ~ ccode + year)
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
    n_clusters_state = n_distinct(d_iv$ccode),
    n_clusters_year  = n_distinct(d_iv$year),
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
iv_levels_state <- c("log(1 + Fatal MID onset) (t-1)")
iv_facet_labels <- c("log(1 + Fatal\nMID onset) (t-1)")

p_dat <- results |>
  filter(iv_label %in% iv_levels_state) |>
  select(iv_label, estimate, p_asy_1way, p_asy_2way, p_wcb_1way, p_wcb_2way) |>
  pivot_longer(-c(iv_label, estimate), names_to = "method", values_to = "pval") |>
  mutate(
    # [2026-10-12] sign-aware one-sided p-value for H: effect > 0 (was pval/2 under a two-sided label)
    pval_one = if_else(estimate > 0, pval / 2, 1 - pval / 2),
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
  ggplot(aes(x = iv_label, y = pval_one,
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
    x = NULL, y = "One-sided p-value (H: effect > 0)",
    color = NULL, shape = NULL
  ) +
  panel_theme +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 1))

ggsave(file.path(FIG_DIR, "fig_ms_wild_bootstrap_churn.png"),
       p_pvals, width = 10, height = 6, dpi = 400, bg = "white")
