###############################################################################
# 23_robustness_nonparticipant_kl.R
# H1 robustness check: does the SYSTEM-LEVEL network reconfigure even
# among states that did not participate in the war?
#
# RATIONALE: The main-text H1 finding is that system-level KL divergence
# (year-over-year network distribution change) elevates during war.
# A predictable critique is that this is mechanical — warring states
# are reshuffling their ties, and aggregating up gives elevated KL by
# construction. This script defuses that critique by restricting the
# system to NON-PARTICIPANT states only and checking whether the
# subnetwork's year-over-year KL still elevates during war.
#
# If non-participant subnetworks ALSO show elevated KL during war
# years, system-level reconfiguration is genuinely systemic — the
# cascade reaches even the states standing on the sidelines. That is
# the strongest evidence for H1 and is consistent with the paper's
# information-revelation mechanism: war reveals information about
# alignment beyond the combatants themselves.
#
# DESIGN: For each of three war episodes we
#   (1) Define non-participants: states present in COW at the pre-war
#       year with NO MID involvement during the war years.
#   (2) Hold the non-participant set fixed across the analysis window
#       (no churn-by-construction from changing membership).
#   (3) Build the per-year subnetwork on `main_weight` (alliance + IGO,
#       same as the manuscript main spec) restricted to dyads where
#       BOTH endpoints are non-participants.
#   (4) Compute year-over-year KL divergence on the subnetwork using
#       the same `kl_divergence()` helper as 02_compute_network_measures.R.
#   (5) Plot full-system vs subnetwork KL time series with war shading.
#   (6) Test for elevation in two ways:
#       (a) Descriptive: pre-war vs during-war mean comparison (Welch t)
#           plus Andrews-Quandt structural-break test on the subnetwork
#           series.
#       (b) Inferential: 3-col HAC NW(5) ladder MIRRORING SCRIPT 01's
#           main system-level OLS specification, with the DV swapped
#           from full-system KL to non-participant subnetwork KL:
#             (1) log(1+KL_subnet) ~ log(1+MID)_{t-1} + poly(year, 2)
#             (2) + log1p(n_democracies, n_autocracies, n_states)
#             (3) + log(1+KL_subnet)_{t-1}  (AR(1))
#           Same panel.id structure (~unit + year, unit = 1L), same
#           Newey-West lag (5). Col 4 (order-period FE) is dropped
#           because each per-war window sits inside a single regime.
#
# WAR EPISODES:
#   - WWII             1939-1945
#   - Korean War       1950-1953
#   - Vietnam War      1965-1975 (US ground combat era)
#
# WHY THESE THREE: each is a major systemic shock, each has a clean
# pre-war baseline period, each has a defensible non-participant set
# from COW data. Three replications turn "an interesting WWII
# observation" into "a generalized cascade pattern."
#
# CAVEATS (state in manuscript):
#   - Non-participants are not "untreated" — they are at network
#     distance >= 2 from the combatants, but the cascade still reaches
#     them. That is the WHOLE POINT of this robustness check: under H1's
#     cascade story, even non-participants reconfigure.
#   - The pre-war baseline is contaminated by other shocks (Spanish
#     Civil War, Munich pre-WWII; etc.). We address this by also showing
#     the full time series rather than relying on a single pre/post test.
#   - The subnetwork is smaller than the full system; KL estimates are
#     noisier. We hold the non-participant set fixed at pre-war to
#     remove dimensional changes from the comparison.
#
# Inputs:  data/processed/dyad_year_network_layers.rds
#          data/processed/system_network_measures.rds  (full-system KL)
# Outputs:
#   data/processed/nonparticipant_kl_timeseries.csv
#   data/processed/nonparticipant_kl_tests.csv
#   data/processed/nonparticipant_kl.rds
#   results/appendix/tableA_nonparticipant_kl.tex
#   code/03_figures/output/manuscript/fig_ms_nonparticipant_kl_<war>.png
#   code/03_figures/output/manuscript/fig_ms_nonparticipant_kl_combined.png
###############################################################################

library(tidyverse)
library(igraph)
library(fixest)
library(strucchange)
select <- dplyr::select

set.seed(123)

OUTPUT_DIR <- here::here("data", "processed")
APPX_DIR   <- here::here("results", "appendix")
FIG_DIR    <- here::here("code", "03_figures", "output", "manuscript")
dir.create(APPX_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG_DIR,  recursive = TRUE, showWarnings = FALSE)

# ----------------------------------------------------------------------------
# Theme + colors
# ----------------------------------------------------------------------------
ACCENT_FULL <- "#2C6CB0"   # full-system KL
ACCENT_SUB  <- "#B03366"   # non-participant subnetwork KL
SHADE_WAR   <- "#888888"

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

save_pres <- function(filename, plot, width = 8.5, height = 5, dpi = 400) {
  ggsave(file.path(FIG_DIR, filename), plot,
         width = width, height = height, dpi = dpi, bg = "white")
  cat(sprintf("  -> %s\n", filename))
}

# ----------------------------------------------------------------------------
# KL divergence helper — IDENTICAL to 02_compute_network_measures.R
# ----------------------------------------------------------------------------
kl_divergence <- function(p, q, eps = 1e-10) {
  p <- p + eps
  q <- q + eps
  p <- p / sum(p)
  q <- q / sum(q)
  sum(p * log(p / q))
}

# ============================================================================
# 0. War definitions
# ============================================================================
# Window covers a long pre-war and post-war baseline so the time-series
# plot shows the war episode as a discontinuity rather than the whole
# panel. Each war's pre-war baseline starts ~10-15 years before onset.

WAR_DEFS <- list(
  list(
    name        = "WWII",
    file_prefix = "wwii",
    start_year  = 1939L,
    end_year    = 1945L,
    window_start = 1925L,
    window_end   = 1955L,
    pre_baseline_start = 1925L,
    pre_baseline_end   = 1938L
  ),
  list(
    name        = "Korean_War",
    file_prefix = "korean",
    start_year  = 1950L,
    end_year    = 1953L,
    window_start = 1946L,
    window_end   = 1965L,
    pre_baseline_start = 1946L,
    pre_baseline_end   = 1949L
  ),
  list(
    name        = "Vietnam_War",
    file_prefix = "vietnam",
    start_year  = 1965L,
    end_year    = 1975L,
    window_start = 1955L,
    window_end   = 1985L,
    pre_baseline_start = 1955L,
    pre_baseline_end   = 1964L
  )
)

# ============================================================================
# 1. Load data
# ============================================================================

cat("Loading dyad-year network layers, full-system KL, and system controls...\n")
dyad_data       <- readRDS(file.path(OUTPUT_DIR, "dyad_year_network_layers.rds"))
system_main_ts  <- readRDS(file.path(OUTPUT_DIR, "system_network_measures.rds"))
# system_analysis_data.rds is the same panel script 01 uses for the
# manuscript main HAC NW(5) ladder. Pulling it here so this robustness
# check uses the IDENTICAL IV transforms and control set.
system_data     <- readRDS(file.path(OUTPUT_DIR, "system_analysis_data.rds"))

stopifnot("main_weight"   %in% names(dyad_data))
stopifnot("kl_divergence" %in% names(system_main_ts))
stopifnot(all(c("mid_count", "fatal_mid_count",
                "n_democracies", "n_autocracies", "n_states") %in% names(system_data)))

# ============================================================================
# 2. Helper: identify war participants (FATAL MID involvement during war years)
# ============================================================================
# RATIONALE: An earlier version of this function used `cowmidongoing == 1`
# (any MID, any hostility level). For long windows like WWII (1939-1945)
# this swept up almost every COW state — 64 of 66 1938 members showed
# some MID, leaving only 3 non-participants. Many of those "MIDs" were
# trivial: border disputes, fishing rights, low-level threats with no
# combat. They are not what we mean by "participated in the war."
#
# Restricting to FATAL MIDs (cowmidfatal == 1) gives the canonical
# definition of war participation: a state had at least one militarized
# dispute involving fatalities during the war years. This filters out
# the trivial cases and matches how COW war participation is typically
# operationalized in the literature.

get_war_participants <- function(start_year, end_year) {
  dyad_data |>
    filter(cowmidongoing == 1L, cowmidfatal == 1L,
           year >= start_year, year <= end_year) |>
    select(ccode1, ccode2) |>
    pivot_longer(c(ccode1, ccode2), values_to = "ccode") |>
    distinct(ccode) |>
    pull(ccode) |>
    sort()
}

# COW members in a given year (any state with at least one observed dyad)
cow_members_in_year <- function(yr) {
  dyad_data |>
    filter(year == yr) |>
    select(ccode1, ccode2) |>
    pivot_longer(c(ccode1, ccode2), values_to = "ccode") |>
    distinct(ccode) |>
    pull(ccode) |>
    sort()
}

# ============================================================================
# 3. Helper: compute subnetwork KL time series
# ============================================================================
# For a given non-participant ccode set, build per-year edge frames using
# main_weight, restricted to dyads where BOTH endpoints are in the set.
# Then compute year-over-year KL divergence on the edge-weight distribution.

compute_subnet_kl <- function(np_ccodes, year_range) {
  # Pre-filter dyads to the non-participant set
  np_dyads <- dyad_data |>
    filter(ccode1 %in% np_ccodes, ccode2 %in% np_ccodes,
           year >= min(year_range) - 1L, year <= max(year_range)) |>
    select(ccode1, ccode2, year, w = main_weight) |>
    mutate(w = replace_na(w, 0))

  # For each consecutive (t-1, t) pair in year_range, compute KL on the
  # edge-weight distribution. Mirrors compute_divergence() in
  # 02_compute_network_measures.R lines 195-205.
  yrs <- sort(unique(year_range))
  out <- vector("list", length(yrs))
  for (i in seq_along(yrs)) {
    yr_t  <- yrs[i]
    yr_t1 <- yr_t - 1L
    e_t  <- np_dyads |> filter(year == yr_t)  |> select(ccode1, ccode2, w_t  = w)
    e_t1 <- np_dyads |> filter(year == yr_t1) |> select(ccode1, ccode2, w_t1 = w)
    merged <- full_join(e_t, e_t1, by = c("ccode1", "ccode2")) |>
      mutate(w_t = replace_na(w_t, 0), w_t1 = replace_na(w_t1, 0))
    if (nrow(merged) == 0 || sum(merged$w_t) == 0 || sum(merged$w_t1) == 0) {
      kl <- NA_real_
    } else {
      kl <- kl_divergence(merged$w_t, merged$w_t1)
    }
    out[[i]] <- tibble(year = yr_t, kl_subnet = kl,
                        n_dyads_t = sum(merged$w_t > 0),
                        n_dyads_t1 = sum(merged$w_t1 > 0))
  }
  bind_rows(out)
}

# ============================================================================
# 4. Run the analysis for each war
# ============================================================================

run_one_war <- function(war) {
  # Extract list elements into locals up front. This avoids a dplyr
  # scoping collision: the `mutate(war = war$name)` below creates a
  # column called `war`, and downstream `filter()` calls would then
  # resolve `war` to that atomic column rather than the function-scope
  # list, throwing "$ operator is invalid for atomic vectors."
  w_name        <- war$name
  w_file_prefix <- war$file_prefix
  w_start       <- war$start_year
  w_end         <- war$end_year
  w_win_start   <- war$window_start
  w_win_end     <- war$window_end
  w_pre_start   <- war$pre_baseline_start
  w_pre_end     <- war$pre_baseline_end

  cat(sprintf("\n%s\n  %s (%d-%d)\n%s\n",
              strrep("=", 70), w_name, w_start, w_end, strrep("=", 70)))

  # 4a. Identify participants and non-participants
  participants <- get_war_participants(w_start, w_end)
  cow_pre <- cow_members_in_year(w_start - 1L)
  np_set  <- setdiff(cow_pre, participants)

  cat(sprintf("  Participants (any MID %d-%d): %d states\n",
              w_start, w_end, length(participants)))
  cat(sprintf("  COW members in %d: %d states\n",
              w_start - 1L, length(cow_pre)))
  cat(sprintf("  Non-participant subnetwork: %d states (held fixed at %d composition)\n",
              length(np_set), w_start - 1L))

  if (length(np_set) < 5) {
    cat("  WARNING: non-participant set < 5 states; subnetwork too small.\n")
    return(NULL)
  }

  # 4b. Compute subnetwork KL
  yrs <- w_win_start:w_win_end
  np_kl <- compute_subnet_kl(np_set, yrs)

  # 4c. Pull full-system KL for the same window for side-by-side comparison
  full_kl <- system_main_ts |>
    filter(year >= w_win_start, year <= w_win_end) |>
    select(year, kl_full = kl_divergence)

  combined <- np_kl |>
    left_join(full_kl, by = "year") |>
    mutate(war_period = year >= w_start & year <= w_end,
           war = w_name)

  # 4d. Statistical tests
  pre_kl <- combined |>
    filter(year >= w_pre_start, year <= w_pre_end) |>
    pull(kl_subnet)
  war_kl <- combined |>
    filter(year >= w_start, year <= w_end) |>
    pull(kl_subnet)

  pre_mean <- mean(pre_kl, na.rm = TRUE)
  war_mean <- mean(war_kl, na.rm = TRUE)
  diff_est <- war_mean - pre_mean

  # Welch t-test on the means
  t_test <- tryCatch(
    t.test(war_kl, pre_kl, var.equal = FALSE),
    error = function(e) NULL
  )
  t_stat <- if (!is.null(t_test)) unname(t_test$statistic) else NA_real_
  t_pval <- if (!is.null(t_test)) t_test$p.value else NA_real_

  # Andrews-Quandt structural-break test using strucchange::sctest.
  # Looks for an unknown break in the mean of kl_subnet across the window.
  ts_kl <- combined |>
    filter(!is.na(kl_subnet)) |>
    arrange(year)
  break_test <- tryCatch(
    sctest(ts_kl$kl_subnet ~ 1, type = "supF",
           data = data.frame(kl_subnet = ts_kl$kl_subnet)),
    error = function(e) NULL
  )
  break_stat <- if (!is.null(break_test)) unname(break_test$statistic) else NA_real_
  break_pval <- if (!is.null(break_test)) break_test$p.value else NA_real_

  # 4d-bis. Mirror script 01's main HAC NW(5) ladder with the SUBNETWORK
  # KL as the outcome. The IV is system-wide log(1+MID count)_{t-1};
  # this is the manuscript's H1 IV and is identical to what script 01
  # uses. The DV is kl_subnet (not logged, as of 2026-10-02), the non-participant subnetwork's
  # year-over-year KL divergence. If the manuscript's headline DV-IV
  # relationship survives when the DV is restricted to non-participants,
  # H1 is genuinely systemic (cascade reaches non-fighters), not just
  # mechanical aggregation of belligerent reconfiguration.
  #
  # Specs mirror script 01 cols 1-3:
  #   (1) IV + poly(year, 2)          — bivariate-with-trend baseline
  #   (2) + log1p({n_dem, n_aut, n_states})  — full controls
  #   (3) + kl_subnet_lag1            — AR(1) momentum
  # Col 4 (order_period FE) from script 01 is dropped because each
  # per-war window sits inside a single regime, so the FE has nothing
  # to absorb here. SE is Newey-West HAC with lag = 5, identical to
  # script 01.

  NW_LAG <- 5L
  CTRL_LOG <- "log1p(n_democracies) + log1p(n_autocracies) + log1p(n_states)"

  reg_data <- combined |>
    filter(!is.na(kl_subnet)) |>
    arrange(year) |>
    left_join(
      system_data |>
        select(year, mid_count, fatal_mid_count,
               n_democracies, n_autocracies, n_states),
      by = "year"
    ) |>
    mutate(
      log_kl_subnet            = log1p(kl_subnet),
      log_kl_subnet_lag1       = lag(log_kl_subnet),
      kl_subnet_lag1           = lag(kl_subnet),
      log_mid_count            = log1p(mid_count),
      log_mid_count_lag1       = lag(log_mid_count),
      log_fatal_mid_count      = log1p(fatal_mid_count),
      log_fatal_mid_count_lag1 = lag(log_fatal_mid_count),
      unit                      = 1L
    ) |>
    filter(!is.na(log_mid_count_lag1),
           !is.na(log_fatal_mid_count_lag1),
           !is.na(n_democracies), !is.na(n_autocracies), !is.na(n_states))

  fit_one_ladder <- function(iv_var, data) {
    if (nrow(data) < 10) return(NULL)
    m1 <- tryCatch(
      feols(as.formula(paste("kl_subnet ~", iv_var, "+ poly(year, 2)")),
            data = data,
            panel.id = ~ unit + year, vcov = NW(lag = NW_LAG)),
      error = function(e) { cat(sprintf("    m1 failed (%s): %s\n", iv_var, conditionMessage(e))); NULL })
    m2 <- tryCatch(
      feols(as.formula(paste("kl_subnet ~", iv_var, "+", CTRL_LOG,
                              "+ poly(year, 2)")),
            data = data,
            panel.id = ~ unit + year, vcov = NW(lag = NW_LAG)),
      error = function(e) { cat(sprintf("    m2 failed (%s): %s\n", iv_var, conditionMessage(e))); NULL })
    m3 <- tryCatch(
      feols(as.formula(paste("kl_subnet ~", iv_var, "+", CTRL_LOG,
                              "+ poly(year, 2) + kl_subnet_lag1")),
            data = data |> filter(!is.na(kl_subnet_lag1)),
            panel.id = ~ unit + year, vcov = NW(lag = NW_LAG)),
      error = function(e) { cat(sprintf("    m3 failed (%s): %s\n", iv_var, conditionMessage(e))); NULL })
    list(m1 = m1, m2 = m2, m3 = m3)
  }

  ladder_mid   <- fit_one_ladder("log_mid_count_lag1",       reg_data)
  ladder_fatal <- fit_one_ladder("log_fatal_mid_count_lag1", reg_data)

  extract_iv_row <- function(fit, iv_var) {
    if (is.null(fit) || !iv_var %in% rownames(coeftable(fit))) return(NULL)
    ct <- coeftable(fit)
    tibble(estimate = ct[iv_var, "Estimate"],
           se       = ct[iv_var, "Std. Error"],
           pval     = ct[iv_var, "Pr(>|t|)"],
           n_obs    = nobs(fit))
  }

  ladder_rows <- bind_rows(
    list(
      tibble(iv = "MID",   col = "(1) IV + trend",        extract_iv_row(ladder_mid$m1,   "log_mid_count_lag1")),
      tibble(iv = "MID",   col = "(2) + Full ctrls",      extract_iv_row(ladder_mid$m2,   "log_mid_count_lag1")),
      tibble(iv = "MID",   col = "(3) + AR(1)",           extract_iv_row(ladder_mid$m3,   "log_mid_count_lag1")),
      tibble(iv = "Fatal", col = "(1) IV + trend",        extract_iv_row(ladder_fatal$m1, "log_fatal_mid_count_lag1")),
      tibble(iv = "Fatal", col = "(2) + Full ctrls",      extract_iv_row(ladder_fatal$m2, "log_fatal_mid_count_lag1")),
      tibble(iv = "Fatal", col = "(3) + AR(1)",           extract_iv_row(ladder_fatal$m3, "log_fatal_mid_count_lag1"))
    )
  ) |>
    mutate(war = w_name) |>
    select(war, iv, col, estimate, se, pval, n_obs)

  # 4e. Print summary
  cat(sprintf("\n  Pre-war (%d-%d) mean KL (subnet): %.5f\n",
              w_pre_start, w_pre_end, pre_mean))
  cat(sprintf("  War-period (%d-%d) mean KL (subnet): %.5f\n",
              w_start, w_end, war_mean))
  cat(sprintf("  Difference: %.5f  Welch t = %.3f  p = %.4f\n",
              diff_est, t_stat, t_pval))
  cat(sprintf("  Andrews-Quandt supF: %.3f  p = %.4f\n",
              break_stat, break_pval))
  cat("\n  --- 3-col HAC NW(5) ladder mirroring script 01 ---\n")
  print(ladder_rows, n = Inf, digits = 4)

  # 4f. Time-series plot
  plot_df <- combined |>
    pivot_longer(c(kl_subnet, kl_full),
                 names_to = "series", values_to = "kl") |>
    mutate(series = factor(series,
                           levels = c("kl_full", "kl_subnet"),
                           labels = c("Full system",
                                      sprintf("Non-participants (n = %d)", length(np_set)))))

  series_labs <- levels(plot_df$series)
  color_vals  <- setNames(c(ACCENT_FULL, ACCENT_SUB), series_labs)
  ltype_vals  <- setNames(c("solid", "dashed"),       series_labs)

  p <- ggplot(plot_df, aes(x = year, y = kl, color = series, linetype = series)) +
    annotate("rect", xmin = w_start - 0.5, xmax = w_end + 0.5,
             ymin = -Inf, ymax = Inf, alpha = 0.18, fill = SHADE_WAR) +
    geom_line(linewidth = 0.9, na.rm = TRUE) +
    geom_point(size = 1.6, na.rm = TRUE) +
    scale_color_manual(values = color_vals, name = NULL) +
    scale_linetype_manual(values = ltype_vals, name = NULL) +
    labs(x = "Year", y = "KL divergence (year-over-year)",
         title = sprintf("Non-participant subnetwork KL: %s",
                          gsub("_", " ", w_name)),
         subtitle = sprintf("War period shaded. %d non-participant states held fixed at %d composition.",
                            length(np_set), w_start - 1L)) +
    theme_pres(base_size = 11)

  save_pres(sprintf("fig_ms_nonparticipant_kl_%s.png", w_file_prefix), p)

  # 4g. Return tidy results
  list(
    war          = w_name,
    participants = participants,
    np_set       = np_set,
    timeseries   = combined,
    ladder_rows  = ladder_rows,
    ladder_mid   = ladder_mid,
    ladder_fatal = ladder_fatal,
    tests = tibble(
      war        = w_name,
      pre_mean   = pre_mean,
      war_mean   = war_mean,
      diff       = diff_est,
      t_stat     = t_stat,
      t_pval     = t_pval,
      break_stat = break_stat,
      break_pval = break_pval,
      n_np       = length(np_set),
      n_part     = length(participants)
    )
  )
}

# Iterate sequentially (per coding-style preference). Wrap each call in
# tryCatch so a failure on one war doesn't kill the others.
results_list <- list()
for (war_def in WAR_DEFS) {
  res <- tryCatch(
    run_one_war(war_def),
    error = function(e) {
      cat(sprintf("\n  *** ERROR in %s: %s ***\n", war_def$name, conditionMessage(e)))
      NULL
    }
  )
  if (!is.null(res)) results_list[[war_def$name]] <- res
}

# ============================================================================
# 5. Combined panel figure: all wars side-by-side
# ============================================================================

ts_combined <- bind_rows(lapply(results_list, function(r) r$timeseries))
ts_combined <- ts_combined |>
  pivot_longer(c(kl_subnet, kl_full),
               names_to = "series", values_to = "kl") |>
  mutate(series = factor(series,
                         levels = c("kl_full", "kl_subnet"),
                         labels = c("Full system", "Non-participants")))

# War shading per facet
war_shade_df <- bind_rows(lapply(WAR_DEFS, function(w) {
  tibble(war = w$name, xmin = w$start_year - 0.5, xmax = w$end_year + 0.5)
})) |> filter(war %in% names(results_list))

p_combined <- ggplot(ts_combined, aes(x = year, y = kl, color = series, linetype = series)) +
  geom_rect(data = war_shade_df,
            aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
            inherit.aes = FALSE, alpha = 0.18, fill = SHADE_WAR) +
  geom_line(linewidth = 0.85, na.rm = TRUE) +
  geom_point(size = 1.4, na.rm = TRUE) +
  facet_wrap(~ war, scales = "free", ncol = 1,
             labeller = as_labeller(function(x) gsub("_", " ", x))) +
  scale_color_manual(values = c("Full system" = ACCENT_FULL,
                                  "Non-participants" = ACCENT_SUB),
                      name = NULL) +
  scale_linetype_manual(values = c("Full system" = "solid",
                                     "Non-participants" = "dashed"),
                         name = NULL) +
  labs(x = "Year", y = "KL divergence (year-over-year)",
       title = "Non-participant subnetwork KL across three war episodes",
       subtitle = "War periods shaded. Full system (blue, solid) vs non-participant subnetwork (red, dashed).") +
  theme_pres(base_size = 11)

save_pres("fig_ms_nonparticipant_kl_combined.png", p_combined,
          width = 9, height = 9)

# ============================================================================
# 6. Save artifacts (CSV, RDS, LaTeX)
# ============================================================================

# Time-series CSV (long form)
ts_out <- bind_rows(lapply(results_list, function(r) r$timeseries))
write_csv(ts_out, file.path(OUTPUT_DIR, "nonparticipant_kl_timeseries.csv"))

# Descriptive tests CSV (pre/during means + Welch + Andrews-Quandt)
tests_out <- bind_rows(lapply(results_list, function(r) r$tests))
write_csv(tests_out, file.path(OUTPUT_DIR, "nonparticipant_kl_tests.csv"))

cat("\n========== DESCRIPTIVE TESTS ==========\n")
print(tests_out, digits = 4, n = Inf)

# Ladder-coefficient CSV: per-war script-01-style 3-col HAC NW(5) ladder
ladder_out <- bind_rows(lapply(results_list, function(r) r$ladder_rows))
write_csv(ladder_out, file.path(OUTPUT_DIR, "nonparticipant_kl_ladder.csv"))

cat("\n========== HAC NW(5) LADDER (mirroring script 01) ==========\n")
print(ladder_out, digits = 4, n = Inf)

# RDS
saveRDS(list(
  war_defs    = WAR_DEFS,
  results     = results_list,
  timeseries  = ts_out,
  tests       = tests_out,
  ladder      = ladder_out
), file.path(OUTPUT_DIR, "nonparticipant_kl.rds"))

# ----- LaTeX summary table (descriptive + ladder coefficient on log MID) -----
# Two panels:
#   Panel A: descriptive pre-war vs during-war mean KL on the subnetwork
#   Panel B: log(1+MID count)_{t-1} coefficient from script-01-style 3-col
#            HAC ladder, applied to kl_subnet (not logged)

tex <- c(
  "\\begin{table}[htbp]", "\\centering",
  "\\caption{H1 robustness check: KL divergence on non-participant ",
  "subnetworks during three major wars. \\textbf{Panel A} reports ",
  "descriptive pre-war vs.\\ during-war mean KL on the subnetwork plus a ",
  "Welch $t$-test and an Andrews-Quandt structural-break test. ",
  "\\textbf{Panel B} mirrors the manuscript main system-level OLS ",
  "(Table~\\ref{tab:system_main}) with the DV restricted to the ",
  "non-participant subnetwork: $\\text{KL}_{\\text{subnet}} ",
  "\\sim \\log(1{+}\\text{MID})_{t{-}1}$ with Newey-West HAC SE ",
  "(lag $=$ 5), reporting the IV coefficient from cols 1--3 of script ",
  "01's ladder (col 4 / order-period FE dropped because each per-war ",
  "window sits in a single regime).}",
  "\\label{tab:nonparticipant_kl}",
  "\\small",
  "\\begin{tabular}{lrrrrr}", "\\toprule",
  "\\multicolumn{6}{l}{\\textit{Panel A: Descriptive pre-war vs.\\ during-war}} \\\\",
  "\\midrule",
  "War & $n_{\\text{np}}$ & Pre-war mean KL & War mean KL & ",
  "Diff. & Welch $p$ \\\\", "\\midrule"
)
for (i in seq_len(nrow(tests_out))) {
  r <- tests_out[i, ]
  tex <- c(tex, sprintf(
    "%s & %d & %.4f & %.4f & %.4f & %s \\\\",
    gsub("_", " ", r$war), r$n_np,
    r$pre_mean, r$war_mean, r$diff,
    if (is.na(r$t_pval)) "---" else sprintf("%.3f", r$t_pval)
  ))
}
tex <- c(tex,
  "\\midrule",
  "\\multicolumn{6}{l}{\\textit{Panel B: HAC NW(5) ladder, $\\log(1+\\text{MID})_{t-1}$ coefficient}} \\\\",
  "\\midrule",
  "War & IV & Spec & $\\hat{\\beta}$ & SE & $p$ \\\\", "\\midrule"
)
for (i in seq_len(nrow(ladder_out))) {
  r <- ladder_out[i, ]
  est_str <- if (is.na(r$estimate)) "---" else sprintf("%.4f", r$estimate)
  se_str  <- if (is.na(r$se))       "---" else sprintf("%.4f", r$se)
  p_str   <- if (is.na(r$pval))     "---" else sprintf("%.3f", r$pval)
  tex <- c(tex, sprintf(
    "%s & %s & %s & %s & %s & %s \\\\",
    gsub("_", " ", r$war), r$iv, r$col, est_str, se_str, p_str
  ))
}
tex <- c(tex, "\\bottomrule", "\\end{tabular}",
  "\\begin{tablenotes}\\small",
  "\\item Non-participants are COW member states present in the year ",
  "before the war with no \\emph{fatal} MID involvement during the war ",
  "years. Subnetwork built on \\texttt{main\\_weight} (alliance + IGO ",
  "composite, the same weight used in the manuscript main spec) ",
  "restricted to dyads where both endpoints are non-participants and ",
  "held fixed at pre-war composition for the entire window. Year-over-",
  "year KL divergence on the edge-weight distribution uses the same ",
  "helper function as \\texttt{02\\_compute\\_network\\_measures.R}. ",
  "\\textbf{Panel B} mirrors script 01's main HAC NW(5) ladder ",
  "(\\texttt{table1\\_system\\_main.tex}) with the DV swapped from ",
  "the full-network divergence to $\\text{KL}_{\\text{subnet}}$ (not logged). ",
  "Spec (1) is IV + quadratic year trend; (2) adds full controls ",
  "($\\log1p$ of n\\_democracies, n\\_autocracies, n\\_states); (3) ",
  "additionally controls for the AR(1) lagged subnetwork KL. Both MID ",
  "and Fatal-MID IVs reported. ",
  "\\textbf{Interpretation:} A positive, significant $\\hat{\\beta}$ on ",
  "$\\log(1{+}\\text{MID})_{t-1}$ in Panel B means system-wide war ",
  "activity drives KL divergence even when measured only on states ",
  "that did not fight — the cascade reaches non-participants, ",
  "supporting H1. Null = the system-level signal in the manuscript ",
  "main results is driven primarily by combatant reconfiguration ",
  "rather than systemic cascade.",
  "\\end{tablenotes}",
  "\\end{table}")
writeLines(tex, file.path(APPX_DIR, "tableA_nonparticipant_kl.tex"))
cat(sprintf("\n  LaTeX -> %s\n",
            file.path("results", "appendix", "tableA_nonparticipant_kl.tex")))

cat("\nDone. 23_robustness_nonparticipant_kl.R complete.\n")
