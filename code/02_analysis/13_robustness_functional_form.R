# [2026-09-26] Headline system-level DV switched from KL to Jensen-Shannon (JS)
# divergence. KL results are retained in results/appendix/kl_canonical/.
###############################################################################
# 15_robustness_functional_form.R
# Two functional-form robustness variations.
#
# (A) BINNED TREATMENT (Idea 6)
#     Replace the linear MID count with categorical bins.
#     Fatal:     {0 (ref), 1, 2-3, 4+}
#     Non-fatal: {0 (ref), 1, 2-3, 4+}
#     Shows the dose-response shape rather than imposing linearity.
#     Fit at state (TWFE) and system (NW HAC-5) levels.
#
# (B) LONG-DIFFERENCE DV (Idea 7)
#     Rebuild the DV over a 5-year horizon:
#       State: weighted Jaccard distance between a state's neighbor weights
#              at year t and year t-5
#       System: KL divergence between edge-weight distribution at year t
#               and year t-5
#     IV: 5-year sum of fatal / non-fatal MID counts.
#     Tests whether the effect operates at annual or multi-year cadence.
#
# Outputs:
#   data/processed/robustness_functional_form_binned.csv
#   data/processed/robustness_functional_form_longdiff.csv
#   code/03_figures/output/manuscript/
#       fig_ms_binned_treatment_kl.png      (system KL, 1 row x 2 cols)
#       fig_ms_binned_treatment_churn.png   (state churn, 1 row x 2 cols)
#       fig_ms_longdiff_kl.png              (system KL, 2x2 IV facet)
#       fig_ms_longdiff_churn.png           (state churn, 2x2 IV facet)
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
# Manuscript theme
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

# Main-text baseline: lagged IVs, plain OLS, no post_1945.
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
dyad_data   <- readRDS(file.path(OUTPUT_DIR, "dyad_year_network_layers.rds"))

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
# [2026-10-12] partner war exposure removed from the control set (H3 dropped from the paper).
# CTRL_MID and CTRL_FATAL are now identical (democracy + autocracy + cinc +
# major_power + is_isolated_t1); both names and ctrl_for_iv() are kept so the
# call sites are unchanged.
CTRL_MID   <- "democracy + autocracy + cinc + major_power + is_isolated_t1"
CTRL_FATAL <- "democracy + autocracy + cinc + major_power + is_isolated_t1"
ctrl_for_iv <- function(iv) if (grepl("fatal", iv)) CTRL_FATAL else CTRL_MID

state_churn <- state_churn |> mutate(log_churn = log1p(churn))
system_data <- system_data |> mutate(log_js_divergence = log1p(js_divergence))

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
# (A) BINNED TREATMENT
# ============================================================

cat("\n--- (A) Binned treatment ---\n")

bin_counts <- function(x) {
  cut(x,
      breaks = c(-Inf, 0.5, 1.5, 3.5, Inf),
      labels = c("0", "1", "2-3", "4+"),
      right = TRUE)
}

# STATE LEVEL — bins for MID (total) and fatal MID, each in its own model
# [2026-10-12] bins built from the t-1 onset counts (headline timing)
bin_state_mid   <- state_churn |> mutate(iv_bin = bin_counts(onset_lag1))
bin_state_fatal <- state_churn |> mutate(iv_bin = bin_counts(fatal_onset_lag1))

# Two formulas: bins built from mid_onset_count use CTRL_MID; bins built
# from fatal_mid_onset_count use CTRL_FATAL.
f_state_bin_mid   <- as.formula(paste0(
  "churn ~ iv_bin + ", CTRL_MID,   " | ccode + year"))
f_state_bin_fatal <- as.formula(paste0(
  "churn ~ iv_bin + ", CTRL_FATAL, " | ccode + year"))

m_state_bin_mid   <- feols(f_state_bin_mid,   data = bin_state_mid,
                             cluster = ~ ccode + year)
m_state_bin_fatal <- feols(f_state_bin_fatal, data = bin_state_fatal,
                             cluster = ~ ccode + year)

# SYSTEM LEVEL — wider bins, applied to LAGGED MID count to match main text.
bin_counts_sys <- function(x) {
  cut(x,
      breaks = c(-Inf, 0.5, 2.5, 5.5, Inf),
      labels = c("0", "1-2", "3-5", "6+"),
      right = TRUE)
}

bin_sys_mid   <- system_data |> mutate(iv_bin = bin_counts_sys(mid_count_lag1))
bin_sys_fatal <- system_data |> mutate(iv_bin = bin_counts_sys(fatal_mid_count_lag1))

f_sys_bin <- js_divergence ~ iv_bin +
  log1p(n_democracies) + log1p(n_autocracies) + log1p(n_states) +
  poly(year, 2)

m_sys_bin_mid   <- lm(f_sys_bin, data = bin_sys_mid)
m_sys_bin_fatal <- lm(f_sys_bin, data = bin_sys_fatal)

extract_bin_coefs <- function(m, level, iv_family_label) {
  cf <- coef(m); vv <- vcov(m)
  bin_terms <- grep("^iv_bin", names(cf), value = TRUE)
  if (length(bin_terms) == 0) return(NULL)
  tibble(
    level = level, iv_family = iv_family_label,
    term  = bin_terms,
    estimate = unname(cf[bin_terms]),
    se       = sqrt(diag(vv)[bin_terms])
  ) |>
    mutate(
      bin_label = sub("^iv_bin", "", term),
      ci_low  = estimate - 1.96 * se,
      ci_high = estimate + 1.96 * se
    )
}

bin_df <- bind_rows(
  extract_bin_coefs(m_state_bin_mid,   "state",  "MID onset (t-1)"),
  extract_bin_coefs(m_state_bin_fatal, "state",  "Fatal MID onset (t-1)"),
  extract_bin_coefs(m_sys_bin_mid,     "system", "MID count"),
  extract_bin_coefs(m_sys_bin_fatal,   "system", "Fatal MID count")
)
write_csv(bin_df,
          file.path(OUTPUT_DIR, "robustness_functional_form_binned.csv"))
cat(sprintf("  Saved binned-coef rows: %d\n", nrow(bin_df)))

# ============================================================
# (B) LONG-DIFFERENCE DV (5-year)
# ============================================================

cat("\n--- (B) Long-difference DV (5-year) ---\n")

HORIZON <- 5

# Build dyad-year main_weight if missing
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

years_sorted <- sort(unique(dyad_data$year))

# --- System-level KL divergence over 5-year lag ---
cat("  Computing 5-year system JS divergence...\n")
kl_long <- tibble(year = integer(), js_divergence_long = double())
eps <- 1e-10

for (yr in years_sorted) {
  if (!(yr - HORIZON) %in% years_sorted) next
  d_t  <- dyad_data |> filter(year == yr,         main_weight > 0) |>
    transmute(key = paste0(pmin(ccode1, ccode2), "_",
                            pmax(ccode1, ccode2)),
              w = main_weight)
  d_s  <- dyad_data |> filter(year == yr - HORIZON, main_weight > 0) |>
    transmute(key = paste0(pmin(ccode1, ccode2), "_",
                            pmax(ccode1, ccode2)),
              w = main_weight)
  if (nrow(d_t) == 0 || nrow(d_s) == 0) next
  p_curr <- setNames(d_t$w / sum(d_t$w), d_t$key)
  p_prev <- setNames(d_s$w / sum(d_s$w), d_s$key)
  all_k <- union(names(p_curr), names(p_prev))
  p <- rep(0, length(all_k)); names(p) <- all_k
  q <- rep(0, length(all_k)); names(q) <- all_k
  p[names(p_curr)] <- p_curr
  q[names(p_prev)] <- p_prev
  p <- (p + eps); p <- p / sum(p)
  q <- (q + eps); q <- q / sum(q)
  m <- 0.5 * (p + q)
  kl <- 0.5 * sum(p * log(p / m)) + 0.5 * sum(q * log(q / m))  # JS divergence
  kl_long <- bind_rows(kl_long,
                        tibble(year = yr, js_divergence_long = kl))
}

# --- State-level churn over 5-year lag ---
cat("  Computing 5-year state churn...\n")
churn_long <- tibble(ccode = integer(), year = integer(),
                     churn_long = double())

# Pre-build ccode -> (year -> neighbor weight vector) for speed
# Long neighbor frame per year
nbr_by_year <- dyad_data |>
  filter(main_weight > 0) |>
  transmute(
    year,
    a = ccode1, b = ccode2, w = main_weight
  ) |>
  {\(d) bind_rows(
    d |> transmute(ccode = a, nbr = b, year, w),
    d |> transmute(ccode = b, nbr = a, year, w)
  )}()

for (yr in years_sorted) {
  if (!(yr - HORIZON) %in% years_sorted) next
  t_df <- nbr_by_year |> filter(year == yr)         |> select(ccode, nbr, w_t   = w)
  s_df <- nbr_by_year |> filter(year == yr - HORIZON) |> select(ccode, nbr, w_tm = w)
  if (nrow(t_df) == 0 && nrow(s_df) == 0) next
  combined <- full_join(t_df, s_df, by = c("ccode", "nbr")) |>
    mutate(w_t = replace_na(w_t, 0), w_tm = replace_na(w_tm, 0))
  churn_yr <- combined |>
    group_by(ccode) |>
    summarise(
      num = sum(pmin(w_t, w_tm)),
      den = sum(pmax(w_t, w_tm)),
      churn_long = if (den > 0) 1 - num / den else NA_real_,
      .groups = "drop"
    ) |> mutate(year = yr) |> select(ccode, year, churn_long)
  churn_long <- bind_rows(churn_long, churn_yr)
}

# Merge back: add 5-year sums (years t-4..t) of MID and fatal MID counts.
# [2026-10-12] Sums are built by CALENDAR year on the full annual series (system) and the full COW
# state-year universe (state), and only where all five years exist. The old row-based rollapplyr ran after
# filtering to defined long differences, so after a gap the five rows could span more than five years.
sys_win <- bind_rows(lapply(0:(HORIZON - 1), function(k)
  system_data |> transmute(year = year + k, mid_count, fatal_mid_count))) |>
  group_by(year) |>
  summarise(n_yrs = n(), mid_sum5 = sum(mid_count), fatal_mid_sum5 = sum(fatal_mid_count), .groups = "drop") |>
  filter(n_yrs == HORIZON) |>
  select(-n_yrs)
sys_long <- system_data |>
  left_join(kl_long, by = "year") |>
  filter(!is.na(js_divergence_long)) |>
  left_join(sys_win, by = "year") |>
  arrange(year) |>
  mutate(
    log_js_divergence_long = log1p(js_divergence_long),
    log_mid_sum5       = log1p(mid_sum5),
    log_fatal_mid_sum5 = log1p(fatal_mid_sum5)
  ) |>
  filter(!is.na(mid_sum5))

st_win <- bind_rows(lapply(0:(HORIZON - 1), function(k)
  onset_univ |> transmute(ccode, year = year + k, u_onset, u_fatal_onset))) |>
  group_by(ccode, year) |>
  summarise(n_yrs = n(), mid_onset_sum5 = sum(u_onset), fatal_onset_sum5 = sum(u_fatal_onset),
            .groups = "drop") |>
  filter(n_yrs == HORIZON) |>
  select(-n_yrs)
st_long <- state_churn |>
  left_join(churn_long, by = c("ccode", "year")) |>
  filter(!is.na(churn_long)) |>
  mutate(log_churn_long = log1p(churn_long)) |>
  left_join(st_win, by = c("ccode", "year")) |>
  mutate(
    log_mid_onset_sum5   = log1p(mid_onset_sum5),
    log_fatal_onset_sum5 = log1p(fatal_onset_sum5)
  ) |>
  filter(!is.na(mid_onset_sum5))

# 4-IV sweep at each level × annual vs 5-year long difference DV
sys_ivs_ann <- list(
  list(col = "mid_count_lag1",           label = "MID count (t-1)"),
  list(col = "log_mid_count_lag1",       label = "log(1 + MID count) (t-1)"),
  list(col = "fatal_mid_count_lag1",     label = "Fatal MID count (t-1)"),
  list(col = "log_fatal_mid_count_lag1", label = "log(1 + Fatal MID count) (t-1)")
)
sys_ivs_long <- list(
  list(col = "mid_sum5",           label = "MID count (5-yr sum)"),
  list(col = "log_mid_sum5",       label = "log(1 + MID count) (5-yr sum)"),
  list(col = "fatal_mid_sum5",     label = "Fatal MID count (5-yr sum)"),
  list(col = "log_fatal_mid_sum5", label = "log(1 + Fatal MID count) (5-yr sum)")
)
# [2026-10-12] annual state IVs at t-1 (headline timing); the 5-year sums below are unchanged
state_ivs_ann <- list(
  list(col = "onset_lag1",           label = "MID onset (t-1)"),
  list(col = "log_onset_lag1",       label = "log(1 + MID onset) (t-1)"),
  list(col = "fatal_onset_lag1",     label = "Fatal MID onset (t-1)"),
  list(col = "log_fatal_onset_lag1", label = "log(1 + Fatal MID onset) (t-1)")
)
state_ivs_long <- list(
  list(col = "mid_onset_sum5",       label = "MID onset"),
  list(col = "log_mid_onset_sum5",   label = "log(1 + MID onset)"),
  list(col = "fatal_onset_sum5",     label = "Fatal MID onset"),
  list(col = "log_fatal_onset_sum5", label = "log(1 + Fatal MID onset)")
)

fit_and_row <- function(f, d, cluster_spec = NULL, engine = "state",
                        iv_col, iv_label, level, spec) {
  if (engine == "lm") {
    m <- lm(f, data = d)
  } else {
    m <- feols(f, data = d, cluster = cluster_spec)
  }
  cf <- coef(m); vv <- vcov(m)
  tibble(spec = spec, level = level,
         iv = iv_col, iv_label = iv_label,
         estimate = unname(cf[iv_col]),
         se       = sqrt(diag(vv)[iv_col]),
         ci_low   = unname(cf[iv_col] - 1.96 * sqrt(diag(vv)[iv_col])),
         ci_high  = unname(cf[iv_col] + 1.96 * sqrt(diag(vv)[iv_col])))
}

ctrl_sys  <- "log1p(n_democracies) + log1p(n_autocracies) + log1p(n_states) + poly(year, 2)"

# System × annual (plain lm to match main text; IV already lagged)
long_sys_ann <- purrr::map_dfr(sys_ivs_ann, function(iv) {
  f <- as.formula(paste0("js_divergence ~ ", iv$col, " + ", ctrl_sys))
  fit_and_row(f, system_data, engine = "lm",
              iv_col = iv$col, iv_label = iv$label,
              level = "system", spec = "Annual DV")
})
# System × long-diff (logged 5-year-difference KL on logged 5-yr sums)
long_sys_long <- purrr::map_dfr(sys_ivs_long, function(iv) {
  f <- as.formula(paste0("js_divergence_long ~ ", iv$col, " + ", ctrl_sys))
  fit_and_row(f, sys_long, engine = "lm",
              iv_col = iv$col, iv_label = iv$label,
              level = "system", spec = "5-year long difference")
})
# State × annual
long_st_ann <- purrr::map_dfr(state_ivs_ann, function(iv) {
  f <- as.formula(paste0("churn ~ ", iv$col, " + ",
                          ctrl_for_iv(iv$col), " | ccode + year"))
  fit_and_row(f, state_churn, cluster_spec = ~ ccode + year,
              iv_col = iv$col, iv_label = iv$label,
              level = "state", spec = "Annual DV")
})
# State × long-diff
long_st_long <- purrr::map_dfr(state_ivs_long, function(iv) {
  f <- as.formula(paste0("churn_long ~ ", iv$col, " + ",
                          ctrl_for_iv(iv$col), " | ccode + year"))
  fit_and_row(f, st_long, cluster_spec = ~ ccode + year,
              iv_col = iv$col, iv_label = iv$label,
              level = "state", spec = "5-year long difference")
})

longdiff_df <- bind_rows(long_sys_ann, long_sys_long,
                          long_st_ann,  long_st_long)

write_csv(longdiff_df,
          file.path(OUTPUT_DIR, "robustness_functional_form_longdiff.csv"))
cat(sprintf("  Saved long-difference rows: %d\n", nrow(longdiff_df)))

# ============================================================
# Plots
# ============================================================

cat("\n--- Plotting ---\n")

# --- Helper: standard panel theme ---
panel_theme <- theme_pres() +
  theme(panel.spacing = unit(1.2, "lines"),
        panel.spacing.x = unit(1.5, "lines"),
        strip.placement = "outside",
        strip.text = element_text(size = 9, face = "bold",
                                    color = "black",
                                    margin = margin(4, 4, 4, 4)),
        plot.margin = margin(t = 10, r = 12, b = 14, l = 10))

# ============================================================
# Binned-treatment dose-response: separate KL vs churn plots
# ============================================================

bin_plot <- bin_df |>
  mutate(bin_order = case_when(
           bin_label %in% c("1", "1-2")  ~ 1L,
           bin_label %in% c("2-3", "3-5")~ 2L,
           bin_label %in% c("4+", "6+")  ~ 3L,
           TRUE ~ 0L
         )) |>
  arrange(level, iv_family, bin_order) |>
  mutate(bin_label = factor(bin_label, levels = unique(bin_label)))

# --- System (KL) ---
bin_sys_plot <- bin_plot |>
  filter(level == "system") |>
  mutate(iv_label = factor(iv_family,
                            levels = c("MID count", "Fatal MID count")))

p20_sys <- ggplot(bin_sys_plot,
                   aes(x = bin_label, y = estimate, group = iv_label)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray55") +
  geom_errorbar(aes(ymin = ci_low, ymax = ci_high),
                width = 0.15, linewidth = 0.7, color = ACCENT) +
  geom_point(size = 2.5, color = ACCENT) +
  geom_line(linewidth = 0.5, color = ACCENT, alpha = 0.6) +
  facet_wrap(~ iv_label, scales = "free", ncol = 2,
             labeller = label_wrap_gen(width = 22)) +
  labs(
    title = "System: binned dose-response",
    x = "MID count bin (reference = 0)",
    y = "Coefficient (95% CI)"
  ) +
  panel_theme

save_pres("fig_ms_binned_treatment_kl.png", p20_sys, width = 9, height = 4.5)

# --- State (churn) ---
bin_st_plot <- bin_plot |>
  filter(level == "state") |>
  mutate(iv_label = factor(iv_family,
                            levels = c("MID onset (t-1)", "Fatal MID onset (t-1)")))

p20_st <- ggplot(bin_st_plot,
                  aes(x = bin_label, y = estimate, group = iv_label)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray55") +
  geom_errorbar(aes(ymin = ci_low, ymax = ci_high),
                width = 0.15, linewidth = 0.7, color = ACCENT) +
  geom_point(size = 2.5, color = ACCENT) +
  geom_line(linewidth = 0.5, color = ACCENT, alpha = 0.6) +
  facet_wrap(~ iv_label, scales = "free", ncol = 2,
             labeller = label_wrap_gen(width = 22)) +
  labs(
    title = "State: binned dose-response",
    x = "MID onset bin (reference = 0)",
    y = "Coefficient (95% CI)"
  ) +
  panel_theme

save_pres("fig_ms_binned_treatment_churn.png", p20_st, width = 9, height = 4.5)

# ============================================================
# Annual vs 5-year long-difference DV: single figure, two IVs:
#   log(1+fatal MID count) [system JS] and log(1+fatal MID onset) [state churn].
# One row of two panels; y-axis = annual vs 5-year difference.
# ============================================================

# Normalize IV family labels (strip "(t-1)" / "(5-yr sum)" suffixes).
clean_iv_family <- function(s) {
  s <- sub("\\s*\\(t-1\\)$",       "", s)
  s <- sub("\\s*\\(5-yr sum\\)$",  "", s)
  s
}

# [2026-10-12] focal IVs are the fatal measures (main-text tables are fatal only)
keep_ivs <- c("log(1 + Fatal MID count)", "log(1 + Fatal MID onset)")

ld_plot <- longdiff_df |>
  mutate(iv_family = clean_iv_family(iv_label)) |>
  filter(iv_family %in% keep_ivs) |>
  mutate(
    iv_label = factor(iv_family, levels = keep_ivs),
    spec     = factor(spec, levels = c("Annual DV", "5-year long difference"),
                      labels = c("Annual", "5-year\ndifference"))
  )

p21 <- ggplot(ld_plot,
              aes(x = estimate, y = forcats::fct_rev(spec))) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray55") +
  geom_errorbarh(aes(xmin = ci_low, xmax = ci_high),
                 height = 0.15, linewidth = 0.7, color = ACCENT) +
  geom_point(size = 2.6, color = ACCENT) +
  facet_wrap(~ iv_label, scales = "free_x", nrow = 1,
             labeller = label_wrap_gen(width = 22)) +
  labs(
    title = "Annual vs 5-year long-difference estimates",
    x = "Coefficient (95% CI)", y = NULL
  ) +
  panel_theme

save_pres("fig_ms_longdiff.png", p21, width = 11, height = 3.6)

cat("\nDone.\n")
