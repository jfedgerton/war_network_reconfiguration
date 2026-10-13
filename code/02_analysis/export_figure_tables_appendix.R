# [2026-09-26] Headline system-level DV switched from KL to Jensen-Shannon (JS)
# divergence. KL results are retained in results/appendix/kl_canonical/.
# =========================================================================== #
# export_figure_tables_appendix.R  (STANDALONE VERSION)
#
# Full regression tables backing the four main-text figures, per APSR
# technical-check request (tables must show DV name, N, and estimates for
# every substantive control; FE may be noted as absorbed).
#
# RUN: source() this file top-to-bottom in a fresh R session from the project
# root. Sections 1-3 load everything they need from data/processed/ (the
# panels your pipeline already writes). Section 4 (Figure 8) needs the
# dyadic estimation sample from 03_contamination_demo.R and will source()
# that script automatically if `est` is not in memory (takes a few minutes).
#
# Output: 12 booktabs .tex tables in results/appendix/.
# Significance codes match the manuscript: *** .001, ** .01, * .05, dagger .1.
# =========================================================================== #

library(tidyverse)
library(fixest)
library(here)

DATA_DIR <- here::here("data", "processed")
APPX_TAB_DIR <- here::here("results", "appendix")
dir.create(APPX_TAB_DIR, recursive = TRUE, showWarnings = FALSE)

SIGNIF <- c("***" = 0.001, "**" = 0.01, "*" = 0.05, "\\dagger" = 0.10)

# Control sets, hard-coded to match the analysis scripts exactly
# (02 line 338; 04 lines 183-192; 05 line 161; 03 lines 304-305).
# [2026-10-12] partner war exposure removed from the control set (H3 dropped from the paper).
CTRL_MID     <- "democracy + autocracy + cinc + major_power + is_isolated_t1"
SYS_CTRL_STR <- "log1p(n_democracies) + log1p(n_autocracies) + log1p(n_states) + poly(year, 2)"
NW_LAG       <- 5L

# --------------------------------------------------------------------------- #
# 0. Load the saved analysis panels (written by 01/02_*_analysis.R)
# --------------------------------------------------------------------------- #

state_panel <- readr::read_csv(file.path(DATA_DIR, "state_churn_data.csv"),
                               show_col_types = FALSE)
system_panel <- readr::read_csv(file.path(DATA_DIR, "system_analysis_data.csv"),
                                show_col_types = FALSE)

stopifnot(all(c("churn", "log_mid_onset_count", "log_fatal_mid_onset_count", "log_fatal_onset_lag1",
                "democracy", "autocracy", "cinc", "major_power",
                "is_isolated_t1") %in% names(state_panel)))
# [2026-10-08] Churn and JS divergence are multiplied by 100 for the appendix tables so that two decimals show
# every estimate (matches main-text Tables 2 and 4). Tables with other dependent variables keep their scale.
state_panel$churn <- 100 * state_panel$churn
system_panel$js_divergence <- 100 * system_panel$js_divergence
stopifnot(all(c("js_divergence", "log_mid_count", "log_mid_count_lag1",
                "log_fatal_mid_count", "log_fatal_mid_count_lag1",
                "n_democracies", "n_autocracies", "n_states", "unit")
              %in% names(system_panel)))

# Shared variable-name dictionary for etable
FIG_DICT <- c(
  log_mid_onset_count       = "log(1 + MID onset)$_{t}$",
  log_mid_count_lag1        = "log(1 + MID count)$_{t-1}$",
  log_fatal_mid_onset_count = "log(1 + fatal MID onset)$_{t}$",
  log_fatal_onset_lag1      = "log(1 + fatal MID onset)$_{t-1}$",
  log_fatal_mid_count_lag1  = "log(1 + fatal MID count)$_{t-1}$",
  iv_lead                   = "IV at lead horizon",
  democracy = "Democracy", autocracy = "Autocracy", cinc = "CINC",
  major_power = "Major power", is_isolated_t1 = "Isolated at $t-1$",
  churn = "Churn", js_divergence = "JS divergence",
  ccode = "State", year = "Year",
  "log1p(n_democracies)" = "log(1 + Democracies)",
  "log1p(n_autocracies)" = "log(1 + Autocracies)",
  "log1p(n_states)"      = "log(1 + States)",
  h0 = "Churn$_{t}$",   h1 = "Churn$_{t+1}$",
  h2 = "Churn$_{t+2}$", h3 = "Churn$_{t+3}$",
  h4 = "Churn$_{t+4}$", h5 = "Churn$_{t+5}$",
  abs_d_polity2  = "$|\\Delta$ Polity$|$",
  mp_flip        = "Major-power flip",
  tc_event       = "Territorial change",
  delta_cinc     = "$\\Delta$ CINC",
  abs_d_n_states = "$|\\Delta$ States in system$|$",
  n_tc_events    = "Territorial-change events",
  tie_alliance_lag_sd  = "Alliance tie$_{t-1}$ (SD)",
  tie_igo_lag_sd       = "IGO tie$_{t-1}$ (SD)",
  tie_multiplex_lag_sd = "Multiplex tie$_{t-1}$ (SD)",
  tie_alliance  = "Alliance tie$_{t}$",
  tie_igo       = "IGO tie$_{t}$",
  tie_multiplex = "Multiplex tie$_{t}$",
  dyad_mid_5yr_lag_sd = "Dyadic conflict history$_{t-2\\ldots t-6}$ (SD)",
  mid5_sum_tp_sd      = "Monadic conflict history$_{t-2\\ldots t-6}$ (SD)",
  cap_between_sd      = "Capability change, between (SD)",
  cap_within_sd       = "Capability change, within (SD)",
  regime_between_sd   = "Regime change, between (SD)",
  regime_within_sd    = "Regime change, within (SD)",
  cap_ratio       = "Capability ratio",
  joint_democracy = "Joint democracy",
  peace_years     = "Peace years",
  "I(peace_years^2)" = "Peace years$^2$",
  "I(peace_years^3)" = "Peace years$^3$",
  cowmidonset = "MID onset$_{t}$",
  dyad_id = "Dyad"
)

# =========================================================================== #
# 1. FIGURE 4 (fig:panel_local_projection) — panel local projections, h = 0..5
# =========================================================================== #

cat("\n== Figure 4 tables: panel local projections ==\n")

lp_data <- state_panel |>
  group_by(ccode) |>
  arrange(year, .by_group = TRUE) |>
  mutate(across(churn, list(
    h0 = ~.x,
    h1 = ~if_else(lead(year, 1) == year + 1, lead(.x, 1), NA_real_),   # [2026-10-12] gap-safe, as in 02
    h2 = ~if_else(lead(year, 2) == year + 2, lead(.x, 2), NA_real_),
    h3 = ~if_else(lead(year, 3) == year + 3, lead(.x, 3), NA_real_),
    h4 = ~if_else(lead(year, 4) == year + 4, lead(.x, 4), NA_real_),
    h5 = ~if_else(lead(year, 5) == year + 5, lead(.x, 5), NA_real_)
  ), .names = "{.fn}")) |>
  ungroup()

lp_mods_full <- list()
lp_mods_none <- list()
for (h in 0:5) {
  # [2026-10-12] fatal MID onset, matching the main LP figure in 02 (figure1_panel_lp.png)
  fml_full <- as.formula(paste0("h", h, " ~ log_fatal_mid_onset_count + ", CTRL_MID,
                                " | ccode + year"))
  fml_none <- as.formula(paste0("h", h, " ~ log_fatal_mid_onset_count | ccode + year"))
  lp_mods_full[[h + 1]] <- feols(fml_full, data = lp_data, vcov = ~ccode + year)
  lp_mods_none[[h + 1]] <- feols(fml_none, data = lp_data, vcov = ~ccode + year)
}

etable(lp_mods_full,
       dict = FIG_DICT, headers = paste0("$h=", 0:5, "$"),
       signif.code = SIGNIF, digits = "r2", digits.stats = 2,
       fitstat = ~ n + r2, tex = TRUE, replace = TRUE,
       title = "Full results for Figure 4: panel local projections of churn$_{t+h}$ on $\\log(1+\\text{fatal MID onset})_t$, full-controls specification. State and year fixed effects absorbed; standard errors two-way clustered on state and year. Churn coefficients and standard errors are multiplied by 100.",
       label = "tab:fig4_lp_full",
       file  = file.path(APPX_TAB_DIR, "tableA_fig4_lp_full.tex"))

etable(lp_mods_none,
       dict = FIG_DICT, headers = paste0("$h=", 0:5, "$"),
       signif.code = SIGNIF, digits = "r2", digits.stats = 2,
       fitstat = ~ n + r2, tex = TRUE, replace = TRUE,
       title = "Full results for Figure 4: panel local projections on $\\log(1+\\text{fatal MID onset})_t$, no-controls specification (IV and state/year fixed effects only). Standard errors two-way clustered on state and year. Churn coefficients and standard errors are multiplied by 100.",
       label = "tab:fig4_lp_none",
       file  = file.path(APPX_TAB_DIR, "tableA_fig4_lp_none.tex"))

cat("  -> tableA_fig4_lp_full.tex, tableA_fig4_lp_none.tex\n")

# =========================================================================== #
# 2. FIGURE 6 (fig:placebo_time) — reverse-timing placebos at t+{1,3,5,10}
# =========================================================================== #

cat("\n== Figure 6 tables: reverse-timing placebos ==\n")

## ---- 2a. State level, fatal MID onset (baseline at t-1, as in 04) ------------
plc_state_mid <- list()
plc_state_mid[["Baseline"]] <- feols(
  as.formula(paste0("churn ~ log_fatal_onset_lag1 + ", CTRL_MID,
                    " | ccode + year")),
  data = state_panel, cluster = ~ ccode + year)

for (h in c(1L, 3L, 5L, 10L)) {
  d_lead <- state_panel |>
    group_by(ccode) |>
    arrange(year, .by_group = TRUE) |>
    mutate(iv_lead = if_else(lead(year, h) == year + h, lead(log_fatal_mid_onset_count, n = h), NA_real_)) |>
    ungroup() |>
    filter(!is.na(iv_lead))
  plc_state_mid[[paste0("t+", h)]] <- feols(
    as.formula(paste0("churn ~ iv_lead + ", CTRL_MID, " | ccode + year")),
    data = d_lead, cluster = ~ ccode + year)
}

etable(plc_state_mid,
       dict = FIG_DICT, headers = names(plc_state_mid),
       signif.code = SIGNIF, digits = "r2", digits.stats = 2,
       fitstat = ~ n + r2, tex = TRUE, replace = TRUE,
       title = "Full results for Figure 5, state level (fatal MID onset): baseline with the IV at $t-1$ and reverse-timing placebos with the IV led to $t+h$. State and year fixed effects absorbed; SEs two-way clustered on state and year. Churn coefficients and standard errors are multiplied by 100.",
       label = "tab:fig6_state_mid",
       file  = file.path(APPX_TAB_DIR, "tableA_fig6_placebo_state_mid.tex"))

## ---- 2c. System level, Newey-West HAC (lag = 5) ---------------------------
# NOTE: restrict to the baseline's estimation sample (years with JS and a t-1 conflict count: 1817-2014,
# N = 198 [2026-10-12: 1817 now keeps its 1816 conflict]). 04_robustness_placebos.R estimates the lead
# models on the same sample, and Figure 6 plots those values. Table 1 is N = 197 because its AR(1) column
# needs JS in 1816; the placebo baseline is Table 1 column 2 on one more year.
# --- Fig 6 system placebo table, on the baseline's estimation sample ---
system_placebo_panel <- system_panel |>
  filter(!is.na(log_mid_count_lag1)) |>   # same first-year drop as 04 (log_fatal_mid_count_lag1 is NA on the same row)
  arrange(year)

plc_sys <- list()
plc_sys[["Baseline"]] <- feols(
  as.formula(paste0("js_divergence ~ log_fatal_mid_count_lag1 + ", SYS_CTRL_STR)),
  data = system_placebo_panel, panel.id = ~ unit + year, vcov = NW(lag = NW_LAG))

for (h in c(1L, 3L, 5L, 10L)) {
  d_lead <- system_placebo_panel |>
    mutate(iv_lead = lead(log_fatal_mid_count, n = h)) |>
    filter(!is.na(iv_lead))
  plc_sys[[paste0("t+", h)]] <- feols(
    as.formula(paste0("js_divergence ~ iv_lead + ", SYS_CTRL_STR)),
    data = d_lead, panel.id = ~ unit + year, vcov = NW(lag = NW_LAG))
}

etable(plc_sys,
       dict = FIG_DICT, headers = names(plc_sys),
       signif.code = SIGNIF, digits = "r2", digits.stats = 2,
       fitstat = ~ n + r2, tex = TRUE, replace = TRUE,
       title = "Full results for Figure 5, system level: baseline (fatal MID count at $t-1$) and reverse-timing placebos with the fatal MID count led to $t+h$. Newey--West HAC standard errors (lag $=5$); quadratic year trend included. JS divergence coefficients and standard errors are multiplied by 100.",
       label = "tab:fig6_system",
       file  = file.path(APPX_TAB_DIR, "tableA_fig6_placebo_system.tex"))

cat("  -> tableA_fig6_placebo_{state_mid,system}.tex\n")

# =========================================================================== #
# 3. FIGURE 7 (fig:sen_mediation) — mediation, 4 state + 2 system channels
# Mediators rebuilt exactly as in 05_robustness_alt_mechanisms.R.
# =========================================================================== #

cat("\n== Figure 7 tables: mediation ==\n")

CTRL_MED <- CTRL_MID       # 05 uses the same state control set
# [2026-10-12] focal IVs follow the fatal t-1 headline (as in 05)
STATE_IV <- "log_fatal_onset_lag1"
SYS_IV   <- "log_fatal_mid_count_lag1"

# --- Territorial change (COW TC v5), as in 05 ---
TC_PATH <- here::here("data", "raw", "tc", "tc2018.csv")
if (file.exists(TC_PATH)) {
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
  stop("COW territorial change data not found at data/raw/tc/tc2018.csv")
}

state_alt <- state_panel |>
  left_join(tc_state_year, by = c("ccode", "year")) |>
  mutate(tc_event = replace_na(tc_event, 0L)) |>
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

state_count_by_year <- state_panel |>
  group_by(year) |>
  summarise(n_states_obs = n_distinct(ccode), .groups = "drop") |>
  arrange(year) |>
  mutate(abs_d_n_states = abs(n_states_obs - lag(n_states_obs))) |>
  dplyr::select(year, abs_d_n_states)

system_alt <- system_panel |>
  left_join(tc_system_year, by = "year") |>
  mutate(n_tc_events = replace_na(n_tc_events, 0L)) |>
  left_join(state_count_by_year, by = "year") |>
  mutate(abs_d_n_states = replace_na(abs_d_n_states, 0)) |>
  filter(!is.na(log_fatal_mid_count_lag1))

STATE_MEDS <- c("abs_d_polity2", "mp_flip", "tc_event", "delta_cinc")
SYS_MEDS   <- c("abs_d_n_states", "n_tc_events")

## ---- 3a. State: mediator equations (a-path) -------------------------------
med_a_state <- list()
for (md in STATE_MEDS) {
  d_md <- state_alt |> filter(!is.na(.data[[md]]))
  med_a_state[[md]] <- feols(
    as.formula(paste0(md, " ~ ", STATE_IV, " + ", CTRL_MED, " | ccode + year")),
    data = d_md, cluster = ~ ccode + year)
}

etable(med_a_state,
       dict = FIG_DICT, headers = unname(FIG_DICT[STATE_MEDS]),
       signif.code = SIGNIF, digits = 4, digits.stats = 3,
       fitstat = ~ n + r2, tex = TRUE, replace = TRUE,
       title = "Full results for Figure 6 (state level), mediator equations (a-path): each candidate mediator regressed on fatal MID onset at $t-1$ and the full control set. State and year fixed effects absorbed; SEs two-way clustered.",
       label = "tab:fig7_state_apath",
       file  = file.path(APPX_TAB_DIR, "tableA_fig7_mediation_state_apath.tex"))

## ---- 3b. State: outcome equations with mediator (b-path + direct) ---------
med_bc_state <- list()
for (md in STATE_MEDS) {
  d_md <- state_alt |> filter(!is.na(.data[[md]]))
  med_bc_state[[md]] <- feols(
    as.formula(paste0("churn ~ ", STATE_IV, " + ", md, " + ", CTRL_MED,
                      " | ccode + year")),
    data = d_md, cluster = ~ ccode + year)
}

etable(med_bc_state,
       dict = FIG_DICT, headers = unname(FIG_DICT[STATE_MEDS]),
       signif.code = SIGNIF, digits = "r2", digits.stats = 2,
       fitstat = ~ n + r2, tex = TRUE, replace = TRUE,
       title = "Full results for Figure 6 (state level), outcome equations: churn on fatal MID onset at $t-1$, the candidate mediator (column header), and the full control set. Indirect effects in Figure 6 are $a \\times b$ with cluster-bootstrap CIs. State and year fixed effects absorbed; SEs two-way clustered. Churn coefficients and standard errors are multiplied by 100.",
       label = "tab:fig7_state_outcome",
       file  = file.path(APPX_TAB_DIR, "tableA_fig7_mediation_state_outcome.tex"))

## ---- 3c. System: mediator and outcome equations ---------------------------
# 05 fits these with lm(); feols with iid vcov reproduces the same point
# estimates and classical SEs in a form etable can print.
med_a_sys  <- list()
med_bc_sys <- list()
for (md in SYS_MEDS) {
  d_md <- system_alt |> filter(!is.na(.data[[md]]))
  med_a_sys[[md]] <- feols(
    as.formula(paste0(md, " ~ ", SYS_IV, " + ", SYS_CTRL_STR)),
    data = d_md, vcov = "iid")
  med_bc_sys[[md]] <- feols(
    as.formula(paste0("js_divergence ~ ", SYS_IV, " + ", md, " + ",
                      SYS_CTRL_STR)),
    data = d_md, vcov = "iid")
}

etable(c(med_a_sys, med_bc_sys),
       dict = FIG_DICT,
       headers = c(paste0("a: ", unname(FIG_DICT[SYS_MEDS])),
                   paste0("outcome: ", unname(FIG_DICT[SYS_MEDS]))),
       signif.code = SIGNIF, digits = "r2", digits.stats = 2,
       fitstat = ~ n + r2, tex = TRUE, replace = TRUE,
       title = "Full results for Figure 6 (system level): mediator equations (first two columns) and outcome equations (last two columns), with full system controls and quadratic year trend. Classical SEs, matching the estimation in the mediation script; indirect-effect CIs in Figure 6 are case-bootstrap over years. Coefficients and standard errors in the JS divergence equations (columns 3--4) are multiplied by 100.",
       label = "tab:fig7_system",
       file  = file.path(APPX_TAB_DIR, "tableA_fig7_mediation_system.tex"))

cat("  -> tableA_fig7_mediation_{state_apath,state_outcome,system}.tex\n")

# =========================================================================== #
# 4. FIGURE 8 (fig:cont_ladder) — contamination step 1 + step 2 ladder
# The 880k-dyad estimation sample `est` is built in 03_contamination_demo.R;
# source it if not already in memory (takes a few minutes).
# =========================================================================== #

cat("\n== Figure 8 tables: contamination ladder ==\n")

if (!exists("est")) {
  cat("  `est` not found -- sourcing 03_contamination_demo.R (this takes a while)...\n")
  source(here::here("code", "02_analysis", "03_contamination_demo.R"))
}
stopifnot(exists("est"))

CTRL_DYAD <- "cap_ratio + joint_democracy + major_power +
              peace_years + I(peace_years^2) + I(peace_years^3)"

ladder_adds <- c("",
                 "+ dyad_mid_5yr_lag_sd",
                 "+ mid5_sum_tp_sd",
                 "+ cap_between_sd", "+ cap_within_sd",
                 "+ regime_between_sd", "+ regime_within_sd")
ladder_labs_short <- c("Baseline", "+ Dyadic hist.", "+ Monadic hist.",
                       "+ Cap. between", "+ Cap. within",
                       "+ Regime between", "+ Regime within")

## ---- 4a. Step 1: tie equations (one column per layer) ---------------------
step1_mods <- list()
for (lyr in c("tie_alliance", "tie_igo", "tie_multiplex")) {
  step1_mods[[lyr]] <- feols(
    as.formula(paste(lyr, "~ dyad_mid_5yr_lag_sd + mid5_sum_tp_sd +",
                     CTRL_DYAD, "| dyad_id + year")),
    data = est, vcov = ~dyad_id)
}

etable(step1_mods,
       dict = FIG_DICT, headers = c("Alliance", "IGO", "Multiplex"),
       signif.code = SIGNIF, digits = 4, digits.stats = 3,
       fitstat = ~ n + r2, tex = TRUE, replace = TRUE,
       title = "Full results for Figure 8, step 1 (tie equations): tie strength at $t$ on dyadic and monadic conflict history over $t-2\\ldots t-6$, with full dyadic controls. Dyad and year fixed effects absorbed; SEs clustered by dyad.",
       label = "tab:fig8_step1",
       file  = file.path(APPX_TAB_DIR, "tableA_fig8_step1_ties.tex"))

## ---- 4b. Step 2: the full 7-rung ladder, one table per layer --------------
LADDER_LAYERS <- c("Alliance"  = "tie_alliance_lag_sd",
                   "IGO"       = "tie_igo_lag_sd",
                   "Multiplex" = "tie_multiplex_lag_sd")

for (k in seq_along(LADDER_LAYERS)) {
  lyr_name <- names(LADDER_LAYERS)[k]
  iv       <- LADDER_LAYERS[[k]]
  ladder_mods <- list()
  for (j in seq_along(ladder_adds)) {
    ladder_mods[[j]] <- feols(
      as.formula(paste("cowmidonset ~", iv, "+", CTRL_DYAD, ladder_adds[j],
                       "| dyad_id + year")),
      data = est, vcov = ~dyad_id)
  }
  etable(ladder_mods,
         dict = FIG_DICT, headers = ladder_labs_short,
         signif.code = SIGNIF, digits = 4, digits.stats = 3,
         fitstat = ~ n + r2, tex = TRUE, replace = TRUE,
         title = paste0("Full results for Figure 8, step 2 (", lyr_name,
                        " layer): MID onset on the standardized ", lyr_name,
                        " tie at $t-1$ across the conditioning ladder. Columns add dyadic conflict history, monadic conflict history, and the four benchmark controls, one at a time, matching the rungs plotted in Figure 8. Dyad and year fixed effects absorbed; SEs clustered by dyad."),
         label = paste0("tab:fig8_ladder_", tolower(lyr_name)),
         file  = file.path(APPX_TAB_DIR,
                           paste0("tableA_fig8_ladder_",
                                  tolower(lyr_name), ".tex")))
  cat(sprintf("  -> tableA_fig8_ladder_%s.tex\n", tolower(lyr_name)))
}

# etable prints the dagger in its significance note outside math mode, which does not compile in LaTeX
# [2026-10-12] explicit list of the tables written above, so retired files (e.g. the PWE placebo table) are not touched
EXPORTED_TABLES <- c("tableA_fig4_lp_full.tex", "tableA_fig4_lp_none.tex",
                     "tableA_fig6_placebo_state_mid.tex", "tableA_fig6_placebo_system.tex",
                     "tableA_fig7_mediation_state_apath.tex", "tableA_fig7_mediation_state_outcome.tex",
                     "tableA_fig7_mediation_system.tex", "tableA_fig8_step1_ties.tex",
                     paste0("tableA_fig8_ladder_", tolower(names(LADDER_LAYERS)), ".tex"))
stopifnot(all(file.exists(file.path(APPX_TAB_DIR, EXPORTED_TABLES))))
for (f in file.path(APPX_TAB_DIR, EXPORTED_TABLES)) {
  x <- readLines(f); writeLines(gsub("\\dagger: 0.1", "$\\dagger$: 0.1", x, fixed = TRUE), f)
}
cat(sprintf("\nAll %d tables written to results/appendix/.\n", length(EXPORTED_TABLES)))
cat("Spot-check: the focal coefficients in each table should reproduce the\n")
cat("point estimates plotted in Figures 4, 6, 7, and 8 (compare with the CSVs written by\n")
cat("02 (LP), 04 (placebos) and 05 (mediation); [2026-10-12] all now use the fatal-MID IVs).\n")

# --------------------------------------------------------------------------- #
# SI: add a section "Full Regression Tables for Main-Text Figures" that
# \input{}s the twelve files. Figure notes to add in main.tex (one sentence):
#   Fig 4: "Full regression results appear in appendix section X."
#   Fig 6: "Full regression results for every horizon appear in appendix section X."
#   Fig 7: "Full mediator and outcome equations appear in appendix section X."
#   Fig 8: "Full results for every ladder rung appear in appendix section X."
# --------------------------------------------------------------------------- #
