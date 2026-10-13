# [2026-09-26] Headline system-level DV switched from KL to Jensen-Shannon (JS)
# divergence. KL results are retained in results/appendix/kl_canonical/.
###############################################################################
# 14_robustness_placebos.R
# Two placebo tests.
#
# (A) WITHIN-STATE PERMUTATION PLACEBO (state level)
#     For each replication, shuffle (fatal, nonfatal) MID onset counts
#     WITHIN each state across its years. This holds the total amount of
#     conflict each state faced over the period constant, but breaks the
#     year-to-year timing relative to churn. Refit the main TWFE model
#     each rep and store coefficients. The observed coefficient should be
#     extreme relative to the permutation null if the timing relationship
#     is genuine.
#
# (B) REVERSE-TIMING PLACEBO (both levels)
#     Shift the IV forward by one year — does conflict at t+1 predict
#     churn at t? If yes, the main result is confounded by something that
#     predicts both network change and upcoming conflict.
#
# Outputs:
#   data/processed/robustness_placebo_permutation.csv
#   data/processed/robustness_placebo_reverse_timing.csv
#   code/03_figures/output/presentation/
#       pres_fig18_permutation.png
#       pres_fig19_reverse_timing.png
###############################################################################

library(tidyverse)
library(fixest)
library(peacesciencer)
select <- dplyr::select

set.seed(123)

OUTPUT_DIR <- here::here("data", "processed")
FIG_DIR    <- here::here("code", "03_figures", "output", "manuscript")
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------
# Manuscript theme — white background, publication-suitable
# ------------------------------------------------------------
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
      axis.ticks   = element_line(color = "gray40"),
      strip.background = element_rect(fill = "gray92", color = "gray70"),
      strip.text   = element_text(face = "bold", size = base_size,
                                    color = "black"),
      legend.position   = "bottom",
      legend.background = element_blank(),
      legend.text       = element_text(size = base_size - 1),
      legend.title      = element_text(size = base_size, face = "bold"),
      plot.margin       = margin(10, 12, 8, 10)
    )
}

save_pres <- function(filename, plot, width = 7, height = 5, dpi = 400) {
  ggsave(file.path(FIG_DIR, filename), plot,
         width = width, height = height, dpi = dpi, bg = "white")
  cat(sprintf("  -> %s\n", filename))
}

# ============================================================
# Load data
# ============================================================

cat("Loading data...\n")

# Main-text baseline: lagged IVs, plain OLS (no HAC), no post_1945.
system_data <- readRDS(file.path(OUTPUT_DIR, "system_analysis_data.rds")) |>
  arrange(year) |>
  # [2026-10-12] conflict lags are built BEFORE dropping years without JS, so 1817 keeps its 1816 conflict
  (\(d) { stopifnot(all(diff(d$year) == 1)); d })() |>   # consecutive years
  mutate(
    log_mid_count       = log1p(mid_count),
    log_fatal_mid_count = log1p(fatal_mid_count),
    mid_count_lag1           = dplyr::lag(mid_count,           1),
    log_mid_count_lag1       = dplyr::lag(log_mid_count,       1),
    fatal_mid_count_lag1     = dplyr::lag(fatal_mid_count,     1),
    log_fatal_mid_count_lag1 = dplyr::lag(log_fatal_mid_count, 1)
  ) |>
  filter(!is.na(mid_count_lag1), !is.na(js_divergence))

state_churn <- readRDS(file.path(OUTPUT_DIR, "state_network_churn.rds"))

dyad_data <- readRDS(file.path(OUTPUT_DIR, "dyad_year_network_layers.rds"))

# Build state-level controls inline if not already attached.
needs_prep <- !all(c("cinc", "major_power", "democracy", "autocracy",
                      "partner_war_frac") %in% names(state_churn))
if (needs_prep) {
  cat("Preparing state-level controls inline (cinc / democracy / partner_war_frac / ...)...\n")

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
  state_churn <- state_churn |> left_join(state_attrs, by = c("ccode", "year"))

  state_years <- create_stateyears(
    subset_years = min(state_churn$year):max(state_churn$year))
  attr(state_years, "ps_data_type") <- "state_year"
  attr(state_years, "ps_system")    <- "cow"
  state_years <- state_years |> add_democracy()
  polity <- state_years |>
    select(ccode, year, polity2) |>
    mutate(democracy = as.integer(polity2 >= 6),
           autocracy = as.integer(polity2 <= -6))
  state_churn <- state_churn |> left_join(polity, by = c("ccode", "year"))

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

# --- State-level control set ---
# Non-fatal and fatal IVs share the same state-level controls:
# democracy + autocracy + cinc + major_power + is_isolated_t1.
# [2026-10-12] partner war exposure removed from the control set (H3 dropped from the paper).
CTRL_MID   <- "democracy + autocracy + cinc + major_power + is_isolated_t1"
CTRL_FATAL <- "democracy + autocracy + cinc + major_power + is_isolated_t1"
ctrl_for_iv <- function(iv) {
  if (grepl("fatal", iv)) CTRL_FATAL
  else                     CTRL_MID
}

state_churn <- state_churn |> mutate(log_churn = log1p(churn))
system_data <- system_data |>
  mutate(log_js_divergence = log1p(js_divergence),
         # `unit` column required by feols panel.id when computing NW HAC SEs
         # on a single time series (treats it as a panel of one unit).
         # Mirrors 01_system_level_analysis.R.
         unit = 1L)

state_churn <- state_churn |>
  mutate(log_mid_onset_count       = log1p(mid_onset_count),
         log_fatal_mid_onset_count = log1p(fatal_mid_onset_count))
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

# ============================================================
# Baseline models (for reference lines on plots) — 4 IVs each level
# ============================================================

# state_controls is now routed per IV via ctrl_for_iv().
sys_ctrl_str   <- "log1p(n_democracies) + log1p(n_autocracies) + log1p(n_states) + poly(year, 2)"

# Focal IVs only -- the placebo figures show just these, so we avoid
# looping the permutation / reverse-timing over the full 4-IV set.
# LEAD ALIGNMENT (important — this was a bug).
# `col`      = the IV used in the BASELINE model.
# `lead_col` = the UNLAGGED conflict series from which placebo leads are built,
#              so that a lead of h means conflict in year t+h (genuinely future).
#
# The system baseline IV is log_mid_count_LAG1 (conflict at t-1). Taking
# lead(log_mid_count_lag1, 1) therefore returns conflict at year t — which is
# CONTEMPORANEOUS with the KL window (KL_t compares the network at t vs t-1),
# not a placebo at all. That mis-alignment is what produced the spurious
# "significant t+1 placebo" at the system level. Building leads from the
# unlagged series (log_mid_count) fixes it.
#
# [2026-10-12] Both baselines are now t-1 IVs (state: log fatal MID onset; system: log fatal MID count).
# Leads are built from the UNLAGGED series, so lead h is conflict in t+h (genuinely future). Note the
# horse race below now pairs conflict at t-1 with conflict at t+h; conflict at t is not in that model.
sys_ivs <- list(
  list(col = "log_fatal_mid_count_lag1", label = "log(1 + Fatal MID count) (t-1)",
       lead_col = "log_fatal_mid_count")
)

# PLACEBO IV at the state level — the headline IV in the main table:
#   MID onset (H2)             : not lagged, so a plain lead() is genuinely future.
#
# `lead_type`:
#   "shift"    -> build the placebo with lead(lead_col, h)
state_ivs <- list(
  list(col = "log_fatal_onset_lag1", label = "log(1 + Fatal MID onset) (t-1)",
       lead_col = "log_fatal_mid_onset_count", lead_type = "shift")
)

# The permutation test shuffles MID ONSET within state, so it only applies to
# the onset IV. Keep the permutation loop on onset alone.
state_perm_ivs <- state_ivs[1]

fit_state_one <- function(d, iv_col) {
  f <- as.formula(paste0("churn ~ ", iv_col, " + ",
                          ctrl_for_iv(iv_col), " | ccode + year"))
  feols(f, data = d, cluster = ~ ccode + year)
}

# System inference MATCHES THE MAIN TEXT: Newey-West HAC (lag = 5), computed
# via feols by treating the single time series as a panel of one unit (same
# construction as 01_system_level_analysis.R). This matters here: iid SEs
# understate uncertainty on a series whose residualized ACF is ~0.46 at lag 1,
# which is precisely the autocorrelation that drives the t+1 placebo leak.
NW_LAG <- 5L

fit_sys_one <- function(d, iv_col) {
  f <- as.formula(paste0("js_divergence ~ ", iv_col, " + ", sys_ctrl_str))
  feols(f, data = d, panel.id = ~ unit + year, vcov = NW(lag = NW_LAG))
}

baseline_rows <- bind_rows(
  purrr::map_dfr(state_ivs, function(iv) {
    m <- fit_state_one(state_churn, iv$col)
    tibble(level = "state", term = iv$col, iv_label = iv$label,
           estimate = unname(coef(m)[iv$col]),
           se = sqrt(diag(vcov(m)))[iv$col])
  }),
  purrr::map_dfr(sys_ivs, function(iv) {
    m <- fit_sys_one(system_data, iv$col)
    tibble(level = "system", term = iv$col, iv_label = iv$label,
           estimate = unname(coef(m)[iv$col]),
           se = sqrt(diag(vcov(m)))[iv$col])
  })
)
baseline_df <- baseline_rows

# ============================================================
# (A) WITHIN-STATE PERMUTATION PLACEBO
# ============================================================

cat("\n--- (A) Within-state permutation placebo ---\n")

N_PERMS <- 1000

perm_rows <- vector("list", N_PERMS)

# Pre-sort so we can shuffle within ccode by year reliably
state_churn_srt <- state_churn |> arrange(ccode, year)

for (r in seq_len(N_PERMS)) {
  if (r %% 100 == 0) cat(sprintf("  perm %d / %d\n", r, N_PERMS))
  # [2026-10-12] Permute the focal IV (log fatal MID onset at t-1) across years within each state; rows
  # where the t-1 IV is undefined (state not in the system in t-1) stay NA.
  perm_d <- state_churn_srt |>
    group_by(ccode) |>
    mutate(
      log_fatal_onset_lag1 = {
        v_ <- log_fatal_onset_lag1; ok_ <- !is.na(v_)
        v_[ok_] <- v_[ok_][sample.int(sum(ok_))]
        v_
      }
    ) |>
    ungroup()

  # Only the onset IV is permuted (see state_perm_ivs note above).
  for (iv in state_perm_ivs) {
    m <- tryCatch(fit_state_one(perm_d, iv$col), error = function(e) NULL)
    if (is.null(m)) next
    cf <- coef(m); vv <- vcov(m)
    perm_rows[[length(perm_rows) + 1]] <- tibble(
      rep = r, level = "state", term = iv$col, iv_label = iv$label,
      estimate = unname(cf[iv$col]),
      se       = sqrt(diag(vv)[iv$col])
    )
  }
}

perm_df <- bind_rows(perm_rows)

# --- System-level permutation (focal IV: log(1+MID count), lagged) ---
# Only one unit at the system level, so shuffle the lagged conflict IV
# across years to break the IV-DV timing alignment, then refit.
cat("  System-level permutation (log fatal MID count)...\n")
sys_focal <- "log_fatal_mid_count_lag1"
sys_perm_rows <- vector("list", N_PERMS)
for (r in seq_len(N_PERMS)) {
  perm_sys <- system_data
  perm_sys[[sys_focal]] <- sample(perm_sys[[sys_focal]])
  m <- tryCatch(fit_sys_one(perm_sys, sys_focal), error = function(e) NULL)
  if (is.null(m)) next
  sys_perm_rows[[r]] <- tibble(
    rep = r, level = "system", term = sys_focal,
    iv_label = "log(1 + Fatal MID count) (t-1)",
    estimate = unname(coef(m)[sys_focal]),
    se       = sqrt(diag(vcov(m))[sys_focal])
  )
}
perm_df <- bind_rows(perm_df, bind_rows(sys_perm_rows))

write_csv(perm_df,
          file.path(OUTPUT_DIR, "robustness_placebo_permutation.csv"))
cat(sprintf("  Saved %d permutation rows\n", nrow(perm_df)))

# Per-IV permutation p-value (two-sided) at state level
state_perm_pvals <- purrr::map_dfr(state_perm_ivs, function(iv) {
  obs <- baseline_df |>
    filter(level == "state", term == iv$col) |> pull(estimate)
  null_draws <- perm_df |> filter(term == iv$col) |> pull(estimate)
  tibble(term = iv$col, iv_label = iv$label,
         observed = obs,
         p_two_sided = mean(abs(null_draws) >= abs(obs)))
})

cat("\n  Per-IV permutation p-values (state level):\n")
print(state_perm_pvals, digits = 4)

# ============================================================
# (B) REVERSE-TIMING PLACEBO
# ============================================================

cat("\n--- (B) Reverse-timing placebo: leads at t+1, t+3, t+5, t+10 ---\n")

# Multi-horizon placebo. The t+1 lead bleeds signal through state-level
# MID autocorrelation; longer leads (t+3, t+5, t+10) break that
# autocorrelation. If the placebo coefficient shrinks toward zero as
# the horizon grows, the t+1 leak is autocorrelation rather than
# anticipation, and the baseline is identified.
PLACEBO_LEADS <- c(1L, 3L, 5L, 10L)

# For each IV at each level and each lead horizon, build the lead and refit.
rev_rows <- list()

# State-level
for (iv in state_ivs) {
  # Baseline (IV at t-1) — pulled from baseline_df once per IV
  rev_rows[[length(rev_rows) + 1]] <- tibble(
    level = "state", iv_label = iv$label, term = iv$col,
    horizon = 0L,
    spec = "Baseline (IV at t-1)",  # [2026-10-12] both baselines are the t-1 focal IVs; 06 keys on horizon, not this label
    estimate = baseline_df |>
      filter(level == "state", term == iv$col) |> pull(estimate),
    se = baseline_df |>
      filter(level == "state", term == iv$col) |> pull(se)
  )

  for (h in PLACEBO_LEADS) {
    d_lead <- state_churn |>
        group_by(ccode) |>
        arrange(year, .by_group = TRUE) |>
        mutate(iv_lead = if_else(lead(year, h) == year + h, lead(.data[[iv$lead_col]], n = h), NA_real_)) |>
        ungroup() |>
        filter(!is.na(iv_lead))
    f_lead <- as.formula(paste0("churn ~ iv_lead + ", ctrl_for_iv(iv$col),
                                 " | ccode + year"))
    m_lead <- tryCatch(
      feols(f_lead, data = d_lead, cluster = ~ ccode + year),
      error = function(e) NULL)
    if (is.null(m_lead)) next

    rev_rows[[length(rev_rows) + 1]] <- tibble(
      level = "state", iv_label = iv$label, term = iv$col,
      horizon = h,
      spec = sprintf("Placebo (IV at t+%d)", h),
      estimate = unname(coef(m_lead)["iv_lead"]),
      se = sqrt(diag(vcov(m_lead))["iv_lead"])
    )
  }
}

# System-level
for (iv in sys_ivs) {
  rev_rows[[length(rev_rows) + 1]] <- tibble(
    level = "system", iv_label = iv$label, term = iv$col,
    horizon = 0L,
    spec = "Baseline (IV at t-1)",  # [2026-10-12] both baselines are the t-1 focal IVs; 06 keys on horizon, not this label
    estimate = baseline_df |>
      filter(level == "system", term == iv$col) |> pull(estimate),
    se = baseline_df |>
      filter(level == "system", term == iv$col) |> pull(se)
  )

  for (h in PLACEBO_LEADS) {
    d_lead <- system_data |>
      arrange(year) |>
      mutate(iv_lead = if_else(lead(year, h) == year + h, lead(.data[[iv$lead_col]], n = h), NA_real_)) |>
      filter(!is.na(iv_lead))
    f_lead <- as.formula(paste0("js_divergence ~ iv_lead + ", sys_ctrl_str))
    m_lead <- tryCatch(
      feols(f_lead, data = d_lead, panel.id = ~ unit + year,
            vcov = NW(lag = NW_LAG)),
      error = function(e) NULL)
    if (is.null(m_lead)) next

    rev_rows[[length(rev_rows) + 1]] <- tibble(
      level = "system", iv_label = iv$label, term = iv$col,
      horizon = h,
      spec = sprintf("Placebo (IV at t+%d)", h),
      estimate = unname(coef(m_lead)["iv_lead"]),
      se = sqrt(diag(vcov(m_lead))["iv_lead"])
    )
  }
}

reverse_df <- bind_rows(rev_rows) |>
  mutate(ci_low  = estimate - 1.96 * se,
         ci_high = estimate + 1.96 * se)

write_csv(reverse_df,
          file.path(OUTPUT_DIR, "robustness_placebo_reverse_timing.csv"))

# ============================================================
# (C) CONDITIONAL ("HORSE-RACE") PLACEBO
# ============================================================
# WHY THIS EXISTS. The placebo in (B) puts the LEAD in the model ALONE. Conflict
# is persistent — at the system level the IV counts ONGOING MIDs, so a war
# running at t is often still running at t+1 — and the residualized ACF of the
# system IV is ~0.46 at lag 1. So the t+1 lead proxies for contemporaneous
# conflict and inherits its signal. That is exactly the pattern observed: the
# system t+1 lead is significant (~82% of baseline) while t+3/5/10 are null,
# tracking the decay of the residualized ACF (0.46 -> 0.17 -> 0.10 -> -0.07).
# The state level, whose IV is MID ONSET (a discrete, non-persisting event),
# shows no t+1 leak — consistent with persistence, not anticipation.
#
# THE DECISIVE TEST. Put conflict at t AND conflict at t+h in the SAME model.
#   - If war -> ties: conflict at t stays significant and the LEAD collapses to
#     null once t is controlled. The leak was contemporaneous signal bleeding
#     through persistence.
#   - If the LEAD SURVIVES conditional on t: that is genuine anticipation or
#     reverse causation, and the baseline is not identified.
# This converts an assertion about autocorrelation into a demonstration.

cat("\n--- (C) Conditional 'horse-race' placebo: conflict at t AND t+h ---\n")

hr_rows <- list()

# --- State level ---
for (iv in state_ivs) {
  for (h in PLACEBO_LEADS) {
    d_hr <- state_churn |>
        group_by(ccode) |>
        arrange(year, .by_group = TRUE) |>
        mutate(iv_t    = .data[[iv$col]],
               iv_lead = if_else(lead(year, h) == year + h, lead(.data[[iv$lead_col]], n = h), NA_real_)) |>
        ungroup() |>
        filter(!is.na(iv_lead), !is.na(iv_t))

    f_hr <- as.formula(paste0("churn ~ iv_t + iv_lead + ",
                               ctrl_for_iv(iv$col), " | ccode + year"))
    m_hr <- tryCatch(feols(f_hr, data = d_hr, cluster = ~ ccode + year),
                     error = function(e) NULL)
    if (is.null(m_hr)) next

    cf <- coef(m_hr); ses <- sqrt(diag(vcov(m_hr)))
    for (trm in c("iv_t", "iv_lead")) {
      hr_rows[[length(hr_rows) + 1]] <- tibble(
        level = "state", iv_label = iv$label, term = iv$col, horizon = h,
        coef_name = if (trm == "iv_t") "Conflict at t-1"
                    else sprintf("Conflict at t+%d", h),
        is_lead   = (trm == "iv_lead"),
        estimate  = unname(cf[trm]),
        se        = unname(ses[trm])
      )
    }
  }
}

# --- System level ---
for (iv in sys_ivs) {
  for (h in PLACEBO_LEADS) {
    d_hr <- system_data |>
      arrange(year) |>
      mutate(iv_t    = .data[[iv$col]],
             iv_lead = if_else(lead(year, h) == year + h, lead(.data[[iv$lead_col]], n = h), NA_real_)) |>
      filter(!is.na(iv_lead), !is.na(iv_t))

    f_hr <- as.formula(paste0("js_divergence ~ iv_t + iv_lead + ",
                               sys_ctrl_str))
    m_hr <- tryCatch(
      feols(f_hr, data = d_hr, panel.id = ~ unit + year,
            vcov = NW(lag = NW_LAG)),
      error = function(e) NULL)
    if (is.null(m_hr)) next

    cf <- coef(m_hr); ses <- sqrt(diag(vcov(m_hr)))
    for (trm in c("iv_t", "iv_lead")) {
      hr_rows[[length(hr_rows) + 1]] <- tibble(
        level = "system", iv_label = iv$label, term = iv$col, horizon = h,
        coef_name = if (trm == "iv_t") "Conflict at t-1"
                    else sprintf("Conflict at t+%d", h),
        is_lead   = (trm == "iv_lead"),
        estimate  = unname(cf[trm]),
        se        = unname(ses[trm])
      )
    }
  }
}

horserace_df <- bind_rows(hr_rows) |>
  mutate(ci_low  = estimate - 1.96 * se,
         ci_high = estimate + 1.96 * se,
         sig     = (ci_low > 0 | ci_high < 0))

write_csv(horserace_df,
          file.path(OUTPUT_DIR, "robustness_placebo_horserace.csv"))
cat(sprintf("  Saved %d horse-race rows\n", nrow(horserace_df)))

cat("\n  Horse-race results (want: 'Conflict at t-1' SIG, lead NULL):\n")
print(horserace_df |>
        select(level, horizon, coef_name, estimate, ci_low, ci_high, sig),
      n = Inf, digits = 3)

# Explicit verdict per level x IV x horizon.
# NOTE: iterates over `term` as well so that (level, term, horizon) identifies
# exactly one t / lead pair per verdict.
cat("\n  Verdict:\n")
for (lv in unique(horserace_df$level)) {
  terms_lv <- unique(horserace_df$term[horserace_df$level == lv])
  for (tm in terms_lv) {
    for (h in PLACEBO_LEADS) {
      sub <- horserace_df |> filter(level == lv, term == tm, horizon == h)
      t_sig    <- sub |> filter(!is_lead) |> pull(sig)
      lead_sig <- sub |> filter(is_lead)  |> pull(sig)
      if (length(t_sig) != 1L || length(lead_sig) != 1L) next

      verdict <- if (t_sig && !lead_sig) {
        "PASS - lead collapses once conflict at t-1 is controlled"
      } else if (t_sig && lead_sig) {
        "FLAG - lead SURVIVES controlling for t-1 (anticipation / reverse causation)"
      } else if (!t_sig && lead_sig) {
        "FAIL - only the lead is significant"
      } else {
        "both null - uninformative"
      }
      cat(sprintf("    %-7s %-26s t+%-2d : %s\n", lv, tm, h, verdict))
    }
  }
}
cat("  Saved reverse-timing placebo results.\n")
print(reverse_df, digits = 4)

# ============================================================
# PLOT A: permutation null distribution
# ============================================================

cat("\n--- Plotting ---\n")

# Focal IVs only: system log(1+MID count) + state log(1+MID onset).
panel_lvls <- c("System: log(1 + Fatal MID count, t-1)", "State: log(1 + Fatal MID onset, t-1)")

perm_focal <- bind_rows(
  perm_df |> filter(level == "system", term == "log_fatal_mid_count_lag1") |>
    mutate(panel = panel_lvls[1]),
  perm_df |> filter(level == "state",  term == "log_fatal_onset_lag1") |>
    mutate(panel = panel_lvls[2])
) |>
  mutate(panel = factor(panel, levels = panel_lvls))

# Observed coefficients + two-sided permutation p-values for the two focal IVs.
obs_sys <- baseline_df |> filter(level == "system", term == "log_fatal_mid_count_lag1") |> pull(estimate)
obs_st  <- baseline_df |> filter(level == "state",  term == "log_fatal_onset_lag1") |> pull(estimate)
null_sys <- perm_df |> filter(level == "system", term == "log_fatal_mid_count_lag1") |> pull(estimate)
null_st  <- perm_df |> filter(level == "state",  term == "log_fatal_onset_lag1") |> pull(estimate)
obs_lines <- tibble(
  panel    = factor(panel_lvls, levels = panel_lvls),
  observed = c(obs_sys, obs_st),
  p_label  = sprintf("p = %.3f", c(mean(abs(null_sys) >= abs(obs_sys)),
                                   mean(abs(null_st)  >= abs(obs_st))))
)

p18 <- ggplot(perm_focal, aes(x = estimate)) +
  geom_histogram(bins = 40, fill = "gray60", color = "gray30", alpha = 0.7) +
  geom_vline(data = obs_lines,
             aes(xintercept = observed), color = ACCENT, linewidth = 1.0) +
  geom_text(data = obs_lines,
            aes(x = observed, y = Inf, label = p_label),
            hjust = 1.25, vjust = 2, color = "black", size = 4) +
  facet_wrap(~ panel, scales = "free", nrow = 1) +
  labs(
    title = "Permutation placebo",
    x = "Coefficient under permuted conflict timing",
    y = "Count"
  ) +
  theme_pres() +
  theme(strip.text = element_text(color = "black", size = 9,
                                    face = "bold",
                                    margin = margin(b = 4)),
        panel.spacing = unit(1.2, "lines"),
        panel.spacing.x = unit(1.5, "lines"),
        plot.margin = margin(t = 10, r = 12, b = 14, l = 10))

save_pres("fig_ms_permutation.png", p18,
          width = 9, height = 4)

# ============================================================
# PLOT B: reverse-timing coefficient comparison
# ============================================================

panel_theme <- theme_pres() +
  theme(panel.spacing = unit(1.2, "lines"),
        panel.spacing.x = unit(1.5, "lines"),
        strip.placement = "outside",
        strip.text = element_text(color = "black", size = 9,
                                    face = "bold",
                                    margin = margin(4, 4, 4, 4)),
        plot.margin = margin(t = 10, r = 12, b = 14, l = 10))

iv_sys_levels   <- sapply(sys_ivs,   `[[`, "label")
iv_state_levels <- sapply(state_ivs, `[[`, "label")

spec_levels <- c("Baseline (IV at t-1)",
                  sprintf("Placebo (IV at t+%d)", PLACEBO_LEADS))
spec_colors <- c("Baseline (IV at t-1)" = "gray20",
                  setNames(c(ACCENT, ACCENT2, "#5C8C2A", "#A65A1F"),
                            sprintf("Placebo (IV at t+%d)", PLACEBO_LEADS)))

rev_plot <- reverse_df |>
  mutate(spec = factor(spec, levels = spec_levels))

# Focal IVs only, single row: system log(1+MID count) + state log(1+MID onset).
panel_lvls_rev <- c("System: log(1 + Fatal MID count, t-1)", "State: log(1 + Fatal MID onset, t-1)")
rev_focal <- rev_plot |>
  filter((level == "system" & term == "log_fatal_mid_count_lag1") |
         (level == "state"  & term == "log_fatal_onset_lag1")) |>
  mutate(panel = factor(
    if_else(level == "system", panel_lvls_rev[1], panel_lvls_rev[2]),
    levels = panel_lvls_rev))

p19 <- ggplot(rev_focal,
              aes(x = estimate, y = forcats::fct_rev(spec), color = spec)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray55") +
  geom_errorbarh(aes(xmin = ci_low, xmax = ci_high),
                 height = 0.15, linewidth = 0.7) +
  geom_point(size = 2.6) +
  scale_color_manual(values = spec_colors, guide = "none") +
  facet_wrap(~ panel, scales = "free_x", nrow = 1) +
  labs(
    title = "Reverse-timing placebo",
    x = "Coefficient (95% CI)", y = NULL
  ) +
  panel_theme

save_pres("fig_ms_reverse_timing.png", p19, width = 10, height = 4)

cat("\nDone.\n")
