# [2026-09-26] Headline system-level DV switched from KL to Jensen-Shannon (JS)
# divergence. KL results are retained in results/appendix/kl_canonical/.
###############################################################################
# 16_robustness_alt_mechanisms.R  (clean rewrite)
#
# Mediation battery: does the conflict -> network-reconfiguration association
# run through alternative DOWNSTREAM channels? If the indirect effect (a*b) is
# ~0 (CI overlaps 0), that channel does not carry the association and the rival
# mechanism is ruled out.
#
# Focal IV (one per level, matching the main-text baseline):
#   State level  : log(1 + MID onset)            -> log(churn)
#   System level : log(1 + MID count), lagged     -> log(KL divergence)
#
# Rival mediators:
#   State  : |Delta polity2|, major-power flip, territorial change, Delta CINC
#   System : |Delta states in system|, territorial-change events
#
# Indirect-effect 95% CIs:
#   State  : cluster bootstrap (resample states, fresh cluster ids)
#   System : case bootstrap over years (approximate; single unit)
#
# Output:
#   data/processed/robustness_altmech_mediation.csv
#   code/03_figures/output/manuscript/fig_ms_altmech_indirect_ci.png
#
# NOTE: battle-deaths (confirmatory) mediation lives in 09; the clean-subsample
# refit and the 4-IV battery were moved to the archived v1 of this script
# (_archive/) and are reported in the appendix.
###############################################################################

library(tidyverse)
library(fixest)
library(peacesciencer)
select <- dplyr::select

set.seed(123)

# Bootstrap replications for the indirect-effect (a*b) CIs. Drop to ~200 for a
# quick check; 1000 for the reported figure.
N_BOOT <- 1000

OUTPUT_DIR <- here::here("data", "processed")
FIG_DIR    <- here::here("code", "03_figures", "output", "manuscript")
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------
# Manuscript theme
# ------------------------------------------------------------
ACCENT  <- "#2C6CB0"

theme_pres <- function(base_size = 12) {
  theme_bw(base_size = base_size) %+replace%
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "gray90", linewidth = 0.3),
      plot.title    = element_text(size = base_size + 2, face = "bold",
                                     hjust = 0, margin = margin(b = 4)),
      plot.subtitle = element_text(size = base_size - 2, color = "gray30",
                                     hjust = 0, margin = margin(b = 8)),
      axis.text     = element_text(color = "black", size = base_size - 1),
      axis.title    = element_text(color = "black", size = base_size),
      axis.ticks    = element_line(color = "gray40"),
      strip.background = element_rect(fill = "gray92", color = "gray70"),
      strip.text    = element_text(face = "bold", size = base_size,
                                     color = "black"),
      panel.spacing = unit(1.1, "lines"),
      plot.margin   = margin(10, 12, 8, 10)
    )
}

save_pres <- function(filename, plot, width = 8, height = 6, dpi = 400) {
  ggsave(file.path(FIG_DIR, filename), plot,
         width = width, height = height, dpi = dpi, bg = "white")
  cat(sprintf("  -> %s\n", filename))
}

# ============================================================
# Load base data
# ============================================================

cat("Loading data...\n")

# Main-text baseline: lagged system IV, plain OLS.
system_data <- readRDS(file.path(OUTPUT_DIR, "system_analysis_data.rds")) |>
  filter(!is.na(js_divergence)) |>
  arrange(year) |>
  mutate(
    log_mid_count      = log1p(mid_count),
    log_mid_count_lag1 = dplyr::lag(log_mid_count, 1)
  ) |>
  filter(!is.na(log_mid_count_lag1))

state_churn <- readRDS(file.path(OUTPUT_DIR, "state_network_churn.rds"))
dyad_data   <- readRDS(file.path(OUTPUT_DIR, "dyad_year_network_layers.rds"))

# Build state-level controls inline if not already attached.
needs_prep <- !all(c("cinc", "major_power", "democracy", "autocracy",
                      "partner_war_frac", "polity2") %in% names(state_churn))
if (needs_prep) {
  cat("Preparing state-level controls inline (cinc / democracy / partner_war_frac / polity2)...\n")

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
  pwe <- partner_sets |>
    mutate(year_next = year + 1L) |>
    left_join(at_war, by = c("ccode2" = "ccode", "year_next" = "year")) |>
    mutate(at_war = replace_na(at_war, 0L)) |>
    group_by(ccode1, year_next) |>
    summarise(
      partner_war_frac          = mean(at_war),
      partner_war_frac_weighted = sum(at_war * w) / sum(w),
      .groups = "drop"
    ) |>
    rename(ccode = ccode1, year = year_next)
  state_churn <- state_churn |>
    left_join(pwe, by = c("ccode", "year")) |>
    mutate(partner_war_frac          = replace_na(partner_war_frac, 0),
           partner_war_frac_weighted = replace_na(partner_war_frac_weighted, 0),
           major_power               = replace_na(major_power, 0))
}

# State-level churn control set (tie-weighted partner-war exposure + attrs).
CTRL <- "partner_war_frac_weighted + democracy + autocracy + cinc + major_power + is_isolated_t1"

state_churn <- state_churn |>
  mutate(log_churn           = log1p(churn),
         log_mid_onset_count = log1p(mid_onset_count))
system_data <- system_data |> mutate(log_js_divergence = log1p(js_divergence))

# --- Territorial change (COW TC v5) ---
TC_PATH <- here::here("data", "raw", "tc", "tc2018.csv")
if (file.exists(TC_PATH)) {
  cat("Loading COW Territorial Change (tc2018.csv)...\n")
  tc_raw <- readr::read_csv(TC_PATH, show_col_types = FALSE)
  tc_state_year <- bind_rows(
    tc_raw |> transmute(ccode = gainer, year, tc_event = 1L),
    tc_raw |> transmute(ccode = loser,  year, tc_event = 1L)
  ) |>
    filter(!is.na(ccode), ccode > 0, !is.na(year)) |>
    distinct(ccode, year, .keep_all = TRUE)
  tc_system_year <- tc_raw |>
    filter(!is.na(year)) |>
    count(year, name = "n_tc_events")
} else {
  cat("TC data not found at ", TC_PATH,
      "; territorial-change mediator will be all-zero (skipped).\n")
  tc_state_year  <- tibble(ccode = integer(), year = integer(),
                            tc_event = integer())
  tc_system_year <- tibble(year = integer(), n_tc_events = integer())
}

state_churn <- state_churn |>
  left_join(tc_state_year, by = c("ccode", "year")) |>
  mutate(tc_event = replace_na(tc_event, 0L))
system_data <- system_data |>
  left_join(tc_system_year, by = "year") |>
  mutate(n_tc_events = replace_na(n_tc_events, 0L))

# ============================================================
# Build the alternative-mechanism (rival) mediators
# ============================================================

cat("Building alt-mechanism covariates...\n")

# --- STATE LEVEL: year-over-year change channels ---
state_alt <- state_churn |>
  arrange(ccode, year) |>
  group_by(ccode) |>
  mutate(
    abs_d_polity2 = abs(polity2 - lag(polity2)),
    mp_flip       = as.integer(major_power != lag(major_power) &
                                 !is.na(lag(major_power))),
    delta_cinc    = cinc - lag(cinc)
  ) |>
  ungroup() |>
  mutate(
    abs_d_polity2 = replace_na(abs_d_polity2, 0),
    mp_flip       = replace_na(mp_flip, 0L),
    delta_cinc    = replace_na(delta_cinc, 0)
  )

# --- SYSTEM LEVEL: year-over-year state-count change ---
# Name the derived count `n_states_obs` and keep ONLY the change measure.
# system_analysis_data.rds already carries `n_states` as a control; joining a
# second `n_states` would collide into n_states.x / n_states.y and break the
# sys_ctrl formula ("object 'n_states' not found").
state_count_by_year <- state_churn |>
  group_by(year) |>
  summarise(n_states_obs = n_distinct(ccode), .groups = "drop") |>
  arrange(year) |>
  mutate(abs_d_n_states = abs(n_states_obs - lag(n_states_obs))) |>
  select(year, abs_d_n_states)

system_alt <- system_data |>
  left_join(state_count_by_year, by = "year") |>
  mutate(abs_d_n_states = replace_na(abs_d_n_states, 0))

# ============================================================
# Baron-Kenny mediation with bootstrap indirect-effect CIs
# ============================================================

coef_se <- function(m) sqrt(diag(vcov(m)))
coef_ps <- function(m) 2 * (1 - pnorm(abs(coef(m) / coef_se(m))))

# --- STATE: TWFE (ccode + year), cluster bootstrap over states ---
bk_state <- function(data, iv, mediator, n_boot = N_BOOT) {
  f_tot <- as.formula(paste0("log_churn ~ ", iv, " + ", CTRL, " | ccode + year"))
  f_a   <- as.formula(paste0(mediator, " ~ ", iv, " + ", CTRL, " | ccode + year"))
  f_bc  <- as.formula(paste0("log_churn ~ ", iv, " + ", mediator, " + ", CTRL,
                             " | ccode + year"))
  d <- data |> filter(!is.na(.data[[mediator]]))

  m_tot <- tryCatch(feols(f_tot, data = d, cluster = ~ ccode + year), error = function(e) NULL)
  m_a   <- tryCatch(feols(f_a,   data = d, cluster = ~ ccode + year), error = function(e) NULL)
  m_bc  <- tryCatch(feols(f_bc,  data = d, cluster = ~ ccode + year), error = function(e) NULL)
  if (is.null(m_tot) || is.null(m_a) || is.null(m_bc)) return(NULL)

  c_tot <- coef(m_tot)[iv]; a <- coef(m_a)[iv]
  b <- coef(m_bc)[mediator]; c_prm <- coef(m_bc)[iv]
  ind <- a * b
  pct <- if (abs(c_tot) > 1e-12) ind / c_tot else NA_real_

  # Cluster bootstrap: resample states, give each drawn cluster a fresh id
  # (so a state drawn twice acts as two independent clusters).
  ci_low <- NA_real_; ci_high <- NA_real_
  if (n_boot > 0) {
    cl_idx <- split(seq_len(nrow(d)), d$ccode)
    n_cl   <- length(cl_idx)
    ind_bt <- rep(NA_real_, n_boot)
    for (brep in seq_len(n_boot)) {
      draw     <- sample.int(n_cl, n_cl, replace = TRUE)
      rows     <- unlist(cl_idx[draw], use.names = FALSE)
      db       <- d[rows, ]
      db$ccode <- rep.int(seq_along(draw), lengths(cl_idx[draw]))
      ma <- tryCatch(feols(f_a,  data = db), error = function(e) NULL)
      mb <- tryCatch(feols(f_bc, data = db), error = function(e) NULL)
      if (is.null(ma) || is.null(mb)) next
      ind_bt[brep] <- tryCatch(coef(ma)[iv] * coef(mb)[mediator],
                               error = function(e) NA_real_)
    }
    qs <- quantile(ind_bt, c(0.025, 0.975), na.rm = TRUE)
    ci_low <- unname(qs[1]); ci_high <- unname(qs[2])
  }

  tibble(
    level = "state", iv = iv, mediator = mediator,
    c_total = unname(c_tot), a_path = unname(a), b_path = unname(b),
    c_prime = unname(c_prm),
    indirect_ab = ind, indirect_ci_low = ci_low, indirect_ci_high = ci_high,
    boot_B = n_boot, pct_mediated = pct, n_obs = nobs(m_tot)
  )
}

# --- SYSTEM: plain OLS, case bootstrap over years (approximate) ---
bk_system <- function(data, iv, mediator, n_boot = N_BOOT) {
  sys_ctrl <- "log1p(n_democracies) + log1p(n_autocracies) + log1p(n_states) + poly(year, 2)"
  f_tot <- as.formula(paste0("log_js_divergence ~ ", iv, " + ", sys_ctrl))
  f_a   <- as.formula(paste0(mediator, " ~ ", iv, " + ", sys_ctrl))
  f_bc  <- as.formula(paste0("log_js_divergence ~ ", iv, " + ", mediator, " + ", sys_ctrl))
  d <- data |> filter(!is.na(.data[[mediator]]))

  m_tot <- tryCatch(lm(f_tot, data = d), error = function(e) NULL)
  m_a   <- tryCatch(lm(f_a,   data = d), error = function(e) NULL)
  m_bc  <- tryCatch(lm(f_bc,  data = d), error = function(e) NULL)
  if (is.null(m_tot) || is.null(m_a) || is.null(m_bc)) return(NULL)

  c_tot <- coef(m_tot)[iv]; a <- coef(m_a)[iv]
  b <- coef(m_bc)[mediator]; c_prm <- coef(m_bc)[iv]
  ind <- a * b
  pct <- if (abs(c_tot) > 1e-12) ind / c_tot else NA_real_

  ci_low <- NA_real_; ci_high <- NA_real_
  if (n_boot > 0) {
    n_row <- nrow(d); ind_bt <- rep(NA_real_, n_boot)
    for (brep in seq_len(n_boot)) {
      db <- d[sample.int(n_row, n_row, replace = TRUE), ]
      ma <- tryCatch(lm(f_a,  data = db), error = function(e) NULL)
      mb <- tryCatch(lm(f_bc, data = db), error = function(e) NULL)
      if (is.null(ma) || is.null(mb)) next
      ind_bt[brep] <- tryCatch(coef(ma)[iv] * coef(mb)[mediator],
                               error = function(e) NA_real_)
    }
    qs <- quantile(ind_bt, c(0.025, 0.975), na.rm = TRUE)
    ci_low <- unname(qs[1]); ci_high <- unname(qs[2])
  }

  tibble(
    level = "system", iv = iv, mediator = mediator,
    c_total = unname(c_tot), a_path = unname(a), b_path = unname(b),
    c_prime = unname(c_prm),
    indirect_ab = ind, indirect_ci_low = ci_low, indirect_ci_high = ci_high,
    boot_B = n_boot, pct_mediated = pct, n_obs = nobs(m_tot)
  )
}

# ============================================================
# Run: one focal IV per level, four/two rival mediators
# ============================================================

cat("\n--- Mediation (log MID onset -> churn; log MID count[t-1] -> JS) ---\n")

STATE_IV <- "log_mid_onset_count"
SYS_IV   <- "log_mid_count_lag1"

state_mediators <- list(
  list(col = "abs_d_polity2", label = "|Δ Polity|"),
  list(col = "mp_flip",       label = "Major-power flip"),
  list(col = "tc_event",      label = "Territorial change"),
  list(col = "delta_cinc",    label = "Δ CINC")
)
system_mediators <- list(
  list(col = "abs_d_n_states", label = "|Δ States in system|"),
  list(col = "n_tc_events",    label = "Territorial-change events")
)

med_rows <- list()
for (m in state_mediators) {
  res <- bk_state(state_alt, STATE_IV, m$col)
  if (!is.null(res)) { res$mediator_label <- m$label; med_rows[[length(med_rows) + 1]] <- res }
}
for (m in system_mediators) {
  res <- bk_system(system_alt, SYS_IV, m$col)
  if (!is.null(res)) { res$mediator_label <- m$label; med_rows[[length(med_rows) + 1]] <- res }
}
med_df <- bind_rows(med_rows)

write_csv(med_df, file.path(OUTPUT_DIR, "robustness_altmech_mediation.csv"))
cat(sprintf("  Saved %d mediation rows.\n", nrow(med_df)))
print(med_df |> select(level, mediator_label, c_total, c_prime, indirect_ab,
                       indirect_ci_low, indirect_ci_high, pct_mediated),
      digits = 4, n = Inf)

# ============================================================
# Single figure: indirect effects (system + state), horizontal
# ============================================================

cat("\n--- Plotting ---\n")

med_plot <- med_df |>
  mutate(
    level_label = factor(if_else(level == "system", "System level", "State level"),
                         levels = c("System level", "State level")),
    mediator_label = factor(mediator_label,
                            levels = c("|Δ States in system|", "Territorial-change events",
                                        "|Δ Polity|", "Major-power flip",
                                        "Territorial change", "Δ CINC"))
  )

p_indirect <- ggplot(med_plot, aes(x = indirect_ab, y = mediator_label)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray55") +
  geom_errorbarh(aes(xmin = indirect_ci_low, xmax = indirect_ci_high),
                 height = 0.16, linewidth = 0.7, color = ACCENT) +
  geom_point(size = 2.7, color = ACCENT) +
  facet_wrap(~ level_label, ncol = 1, scales = "free") +
  labs(
    title    = "Alternative-mechanism mediation: indirect effects",
    subtitle = "Point = indirect effect (a b path product); bar = 95% bootstrap CI. A CI overlapping 0 rules out that channel.",
    x = expression("Indirect effect ("*a %*% b*")"), y = NULL
  ) +
  theme_pres() +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5))

save_pres("fig_ms_altmech_indirect_ci.png", p_indirect, width = 8, height = 6)

cat("\nDone.\n")
