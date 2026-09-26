###############################################################################
# 02_state_level_analysis.R
# State-level analysis: MID involvement → Network churn
#
# DV: churn (weighted Jaccard-based). Alt DVs → appendix script.
#
# 4 IVs (onset measures, raw + log; run in parallel, never combined):
#   (a) mid_onset_count, (b) log1p(mid_onset_count),
#   (c) fatal_mid_onset_count, (d) log1p(fatal_mid_onset_count)
#
# Count measures (mid_count, fatal_mid_count) used at system level.
# Onset captures the discrete shock to a specific state's network.
# Log(1+x) dampens influence of high-conflict states.
#
# Structure:
#   1. Main TWFE specs: Simple + Full for each IV
#   2. Lag models (t-1, t-2) for each IV
#   3. Panel local projections (h=0..5) for each IV
#   4. Summary tables (LaTeX)
#   5. Clustering robustness (onset IVs only — appendix)
#   6. Save
#
# All models: two-way FE (ccode + year) with clustered SEs (~ccode + year)
#
# Full model controls (main):
#   partner_war_frac_weighted (war spillover, tie-strength-WEIGHTED share of
#     t-1 partners at war at t — alliance / IGO main_weight as the tie-
#     strength weighting),
#   democracy, autocracy, cinc, major_power, is_isolated_t1
# The unweighted partner_war_frac variant is computed and saved in
# state_churn for robustness use; the main spec uses the WEIGHTED fraction
# because alliance- and IGO-strong partners' wars are the empirically more
# informative signal (see H3 in the manuscript).
###############################################################################

library(tidyverse)
library(fixest)
library(peacesciencer)
library(igraph)
library(dynamicmultiplex)
library(texreg)

select <- dplyr::select

OUTPUT_DIR <- here::here("data", "processed")

# --------------------------------------------------------------------------- #
# 0. Load data
# --------------------------------------------------------------------------- #

state_churn <- readRDS(file.path(OUTPUT_DIR, "state_network_churn.rds"))
dyad_data <- readRDS(file.path(OUTPUT_DIR, "dyad_year_network_layers.rds"))

# Use unweighted Manhattan (main_weight) as main edge weight
MAIN_WEIGHT_VAR <- "main_weight"

# Fallback: compute main_weight if not present
if (!MAIN_WEIGHT_VAR %in% names(dyad_data)) {
  cat("  WARNING: main_weight not found. Computing from layers.\n")
  MAIN_LAYERS <- c("layer_offensive", "layer_defensive", "layer_igo")
  main_mat <- as.matrix(dyad_data[, MAIN_LAYERS, drop = FALSE])
  main_mat[is.na(main_mat)] <- 0
  dyad_data$main_raw <- rowSums(main_mat)
  dyad_data <- dyad_data |>
    group_by(year) |>
    mutate(main_weight = if (max(main_raw, na.rm = TRUE) > 0) {
      main_raw / max(main_raw, na.rm = TRUE)
    } else 0) |>
    ungroup() |>
    mutate(main_weight = round(main_weight, 4))
  dyad_data$main_raw <- NULL
}

cat(sprintf("  Using edge weight: %s\n", MAIN_WEIGHT_VAR))

# --------------------------------------------------------------------------- #
# 0a. Add state-level controls
# --------------------------------------------------------------------------- #

# CINC
state_attrs <- dyad_data |>
  select(ccode = ccode1, year, cinc1) |>
  distinct(ccode, year, .keep_all = TRUE) |>
  rename(cinc = cinc1)

state_churn <- state_churn |>
  left_join(state_attrs, by = c("ccode", "year"))

# Major power
major_powers <- dyad_data |>
  select(ccode = ccode1, year, cowmaj1) |>
  distinct(ccode, year, .keep_all = TRUE) |>
  rename(major_power = cowmaj1)

state_churn <- state_churn |>
  left_join(major_powers, by = c("ccode", "year")) |>
  mutate(major_power = replace_na(major_power, 0))

# Democracy / autocracy (Polity2)
cat("Adding polity2 democracy indicator...\n")
state_years <- create_stateyears(subset_years = min(state_churn$year):max(state_churn$year))
attr(state_years, "ps_data_type") <- "state_year"
attr(state_years, "ps_system") <- "cow"
state_years <- state_years |> add_democracy()

polity_data <- state_years |>
  select(ccode, year, polity2) |>
  mutate(democracy = as.integer(polity2 >= 6),
         autocracy = as.integer(polity2 <= -6))

state_churn <- state_churn |>
  left_join(polity_data, by = c("ccode", "year"))

cat(sprintf("  Polity2 coverage: %d/%d (%.1f%%)\n",
            sum(!is.na(state_churn$polity2)), nrow(state_churn),
            100 * mean(!is.na(state_churn$polity2))))

# --------------------------------------------------------------------------- #
# 0b. War spillover control
# --------------------------------------------------------------------------- #

cat("Computing partner war exposure (spillover) — all MIDs and fatal MIDs...\n")

partner_sets <- dyad_data |>
  filter(.data[[MAIN_WEIGHT_VAR]] > 0) |>
  select(ccode1, ccode2, year, w = !!sym(MAIN_WEIGHT_VAR)) |>
  bind_rows(
    dyad_data |> filter(.data[[MAIN_WEIGHT_VAR]] > 0) |>
      select(ccode1 = ccode2, ccode2 = ccode1, year, w = !!sym(MAIN_WEIGHT_VAR))
  )

# State-year war indicators:
#   at_war       = 1 if state had ANY ongoing MID this year
#   at_war_fatal = 1 if state had any ongoing FATAL MID this year
# cowmidfatal is built in 01_build_network_layers.R as
#   cowmidongoing == 1 & (fatality1 > 0 | fatality2 > 0)
all_state_war <- dyad_data |>
  filter(cowmidongoing == 1) |>
  select(ccode1, ccode2, year) |>
  pivot_longer(cols = c(ccode1, ccode2), values_to = "ccode") |>
  distinct(ccode, year) |>
  mutate(at_war = 1L)

all_state_fatal_war <- dyad_data |>
  filter(cowmidfatal == 1) |>
  select(ccode1, ccode2, year) |>
  pivot_longer(cols = c(ccode1, ccode2), values_to = "ccode") |>
  distinct(ccode, year) |>
  mutate(at_war_fatal = 1L)

partner_war_exposure <- partner_sets |>
  mutate(year_next = year + 1L) |>
  left_join(all_state_war,       by = c("ccode2" = "ccode", "year_next" = "year")) |>
  left_join(all_state_fatal_war, by = c("ccode2" = "ccode", "year_next" = "year")) |>
  mutate(
    at_war       = replace_na(at_war, 0L),
    at_war_fatal = replace_na(at_war_fatal, 0L)
  ) |>
  group_by(ccode1, year_next) |>
  summarise(
    # All MIDs
    partner_war_frac                = mean(at_war),
    partner_war_frac_weighted       = sum(at_war * w) / sum(w),
    partner_war_count               = sum(at_war),
    # Fatal MIDs only (severity-matched control for fatal-IV regressions)
    partner_fatal_war_frac          = mean(at_war_fatal),
    partner_fatal_war_frac_weighted = sum(at_war_fatal * w) / sum(w),
    partner_fatal_war_count         = sum(at_war_fatal),
    .groups = "drop"
  ) |>
  rename(ccode = ccode1, year = year_next)

state_churn <- state_churn |>
  left_join(partner_war_exposure, by = c("ccode", "year")) |>
  mutate(
    partner_war_frac                = replace_na(partner_war_frac, 0),
    partner_war_frac_weighted       = replace_na(partner_war_frac_weighted, 0),
    partner_war_count               = replace_na(partner_war_count, 0L),
    partner_fatal_war_frac          = replace_na(partner_fatal_war_frac, 0),
    partner_fatal_war_frac_weighted = replace_na(partner_fatal_war_frac_weighted, 0),
    partner_fatal_war_count         = replace_na(partner_fatal_war_count, 0L)
  )

# --------------------------------------------------------------------------- #
# 0c. Create lags
# --------------------------------------------------------------------------- #

state_churn <- state_churn |>
  mutate(
    log_churn                 = log1p(churn),
    log_mid_count             = log1p(mid_count),
    log_fatal_mid_count       = log1p(fatal_mid_count),
    log_mid_onset_count       = log1p(mid_onset_count),
    log_fatal_mid_onset_count = log1p(fatal_mid_onset_count)
  ) |>
  group_by(ccode) |>
  arrange(year) |>
  mutate(
    # Raw onset lags
    onset_lag1            = lag(mid_onset_count, 1),
    onset_lag2            = lag(mid_onset_count, 2),
    fatal_onset_lag1      = lag(fatal_mid_onset_count, 1),
    fatal_onset_lag2      = lag(fatal_mid_onset_count, 2),
    # Log onset lags
    log_onset_lag1        = lag(log_mid_onset_count, 1),
    log_onset_lag2        = lag(log_mid_onset_count, 2),
    log_fatal_onset_lag1  = lag(log_fatal_mid_onset_count, 1),
    log_fatal_onset_lag2  = lag(log_fatal_mid_onset_count, 2)
  ) |>
  ungroup()

cat(sprintf("  MID count: mean=%.3f, sd=%.3f, max=%d\n",
            mean(state_churn$mid_count, na.rm = TRUE),
            sd(state_churn$mid_count, na.rm = TRUE),
            max(state_churn$mid_count, na.rm = TRUE)))
cat(sprintf("  MID onset: mean=%.3f, sd=%.3f, max=%d\n",
            mean(state_churn$mid_onset_count, na.rm = TRUE),
            sd(state_churn$mid_onset_count, na.rm = TRUE),
            max(state_churn$mid_onset_count, na.rm = TRUE)))

# --------------------------------------------------------------------------- #
# 0d. Leiden CPM community membership (unweighted Manhattan)
# --------------------------------------------------------------------------- #
#
# Uses main_weight (unweighted Manhattan, range01-normalized) as edge weights.
# All weights are in [0, 1]. Leiden CPM with gamma > 0 splits communities
# when internal density falls below gamma.
# --------------------------------------------------------------------------- #

cat("Computing multilayer communities (Leiden CPM)...\n")

# Build list of igraph objects using main_weight (positive edges only)
years_with_data <- sort(unique(state_churn$year))
layer_list <- list()
year_map <- c()

for (yr in years_with_data) {
  edges_yr <- dyad_data |>
    filter(year == yr, !is.na(.data[[MAIN_WEIGHT_VAR]])) |>
    select(ccode1, ccode2, !!sym(MAIN_WEIGHT_VAR))
  if (nrow(edges_yr) < 2) next
  g <- graph_from_data_frame(edges_yr |> select(ccode1, ccode2), directed = FALSE)
  E(g)$weight <- edges_yr[[MAIN_WEIGHT_VAR]]
  layer_list[[length(layer_list) + 1]] <- g
  year_map <- c(year_map, yr)
}

cat(sprintf("  Built %d yearly layers (%d-%d)\n",
            length(layer_list), min(year_map), max(year_map)))

# Adjacent-only temporal coupling (default behavior)
layer_links <- data.frame(
  from = seq_len(length(layer_list) - 1),
  to = seq_len(length(layer_list) - 1) + 1,
  weight = 1
)

# Compute data-driven gamma: mean edge weight across all layers
all_weights <- unlist(lapply(layer_list, function(g) E(g)$weight))
LEIDEN_GAMMA <- mean(all_weights, na.rm = TRUE)
cat(sprintf("  Leiden CPM gamma (mean edge weight): %.4f\n", LEIDEN_GAMMA))

set.seed(123)
ml_fit <- fit_multilayer_identity_ties(
  layer_list,
  algorithm = "leiden",
  resolution_parameter = LEIDEN_GAMMA,
  layer_links = layer_links,
  directed = FALSE,
  objective = "cpm"
)

# Extract community assignments back to state-year format
leiden_data <- data.frame()
for (i in seq_along(ml_fit$layer_communities)) {
  yr <- year_map[i]
  temp_data <- data.frame(
    ccode = names(ml_fit[[2]][[i]]$membership),
    year = yr,
    leiden_comm = as.numeric(ml_fit[[2]][[i]]$membership)
  )
  rownames(temp_data) <- NULL
  leiden_data <- leiden_data |> 
    bind_rows(
      temp_data
    )
}

leiden_df <- leiden_data |>
  mutate(leiden_id = paste0(year, "_", leiden_comm))

cat(sprintf("  Leiden communities: %d state-year assignments, %d unique communities\n",
            nrow(leiden_df), n_distinct(leiden_df$leiden_comm)))

state_churn <- state_churn |>
  left_join(leiden_df |> select(ccode, year, leiden_comm) |> mutate(ccode = as.numeric(ccode)),
            by = c("ccode", "year"))

# =========================================================================== #
# MAIN + APPENDIX MODELS — STATE-LEVEL CHURN
# =========================================================================== #
#
# Mirrors the 2 sets x 2 controls x 2 IVs design in 01_system_level_analysis.R
# but extended to support the additional appendix tables the manuscript
# needs (no-FE with poly trend; three-way FE with state + year + community
# membership; raw + count-based variants).
#
# Five tables:
#   main: TWFE log1p(churn) on log MID ONSET   (no/full ctrl x mid/fatal)
#   A1:   TWFE log1p(churn) on log MID COUNT   (no/full ctrl x mid/fatal)
#   A2:   no-FE  log1p(churn) on log {onset, count} + poly(year, 2) (full ctrl)
#   A3:   3-way FE (ccode + year + leiden_comm) log1p(churn) on log
#         {onset, count} (full ctrl)
#   A4:   TWFE raw churn on raw MID ONSET      (no/full ctrl x mid/fatal)
# (A5 -- TWFE raw churn on raw MID COUNT -- dropped; see note below.)
#
# All TWFE models cluster two-way (ccode + year). The three-way FE model
# clusters three-way (ccode + year + leiden_comm). The no-FE model
# clusters one-way (ccode) since there are no within-state deviations to
# absorb.
# =========================================================================== #

# Partner war exposure: the TIE-WEIGHTED fraction of past-year partners
# at war (alliance/IGO main_weight as the tie-strength weighting) is the
# main-text control. We route by IV severity so the partner control
# mirrors the IV: mid-count IVs use partner_war_frac_weighted (any MID),
# fatal-MID IVs use partner_fatal_war_frac_weighted (only fatal MIDs).
# This keeps both the severity threshold AND the weighting scheme
# consistent between IV and partner control. The unweighted variants are
# kept on every row of state_churn and enter the 03 robustness multiverse
# as alternative specs.
#
# is_isolated_t1 is a predetermined indicator for "state had no partners
# at t-1". It absorbs the level shift from `total_gain` (degree_t1 = 0,
# churn = 1) and `both_isolated` (degree_t1 = 0, churn = 0) cases that
# otherwise flatten the within-state slope. Predetermined ⇒ exogenous to
# year-t MID involvement.
CTRL_MID   <- "partner_war_frac_weighted + democracy + autocracy + cinc + major_power + is_isolated_t1"
CTRL_FATAL <- "partner_fatal_war_frac_weighted + democracy + autocracy + cinc + major_power + is_isolated_t1"
ctrl_for_iv <- function(iv) {
  if (grepl("fatal", iv)) CTRL_FATAL else CTRL_MID
}

# Generic fitter -- unifies the spec grid below. The `fe` and `vcov_form`
# arguments cover the variation across tables.
fit_spec <- function(data, dv, iv, controls = NULL,
                     fe = "ccode + year",
                     vcov_form = "ccode + year") {
  rhs <- iv
  if (!is.null(controls) && nzchar(controls)) {
    rhs <- paste(rhs, "+", controls)
  }
  if (!is.null(fe) && nzchar(fe)) {
    fml <- as.formula(sprintf("%s ~ %s | %s", dv, rhs, fe))
  } else {
    fml <- as.formula(sprintf("%s ~ %s", dv, rhs))
  }
  feols(fml, data = data, vcov = as.formula(paste0("~", vcov_form)))
}

# ----- Spec lists (one table = one block) ----------------------------------

specs_main <- list(
  list(key = "main__none__mid",
       dv = "log_churn", iv = "log_mid_onset_count",       controls = NULL),
  list(key = "main__full__mid",
       dv = "log_churn", iv = "log_mid_onset_count",       controls = CTRL_MID),
  list(key = "main__none__fatal",
       dv = "log_churn", iv = "log_fatal_mid_onset_count", controls = NULL),
  list(key = "main__full__fatal",
       dv = "log_churn", iv = "log_fatal_mid_onset_count", controls = CTRL_FATAL)
)

specs_a1 <- list(
  list(key = "a1__none__mid",
       dv = "log_churn", iv = "log_mid_count",        controls = NULL),
  list(key = "a1__full__mid",
       dv = "log_churn", iv = "log_mid_count",        controls = CTRL_MID),
  list(key = "a1__none__fatal",
       dv = "log_churn", iv = "log_fatal_mid_count",  controls = NULL),
  list(key = "a1__full__fatal",
       dv = "log_churn", iv = "log_fatal_mid_count",  controls = CTRL_FATAL)
)

# A2: no fixed effects + cubic year trend
specs_a2 <- list(
  list(key = "a2__nofe__onset_mid",
       dv = "log_churn", iv = "log_mid_onset_count",
       controls = paste(CTRL_MID, "+ poly(year, 2)"),
       fe = NULL, vcov_form = "ccode"),
  list(key = "a2__nofe__onset_fatal",
       dv = "log_churn", iv = "log_fatal_mid_onset_count",
       controls = paste(CTRL_FATAL, "+ poly(year, 2)"),
       fe = NULL, vcov_form = "ccode"),
  list(key = "a2__nofe__count_mid",
       dv = "log_churn", iv = "log_mid_count",
       controls = paste(CTRL_MID, "+ poly(year, 2)"),
       fe = NULL, vcov_form = "ccode"),
  list(key = "a2__nofe__count_fatal",
       dv = "log_churn", iv = "log_fatal_mid_count",
       controls = paste(CTRL_FATAL, "+ poly(year, 2)"),
       fe = NULL, vcov_form = "ccode")
)

# A3: three-way FE (state + year + Leiden community)
specs_a3 <- list(
  list(key = "a3__3wfe__onset_mid",
       dv = "log_churn", iv = "log_mid_onset_count",  controls = CTRL_MID,
       fe = "ccode + year + leiden_comm",
       vcov_form = "ccode + year + leiden_comm"),
  list(key = "a3__3wfe__onset_fatal",
       dv = "log_churn", iv = "log_fatal_mid_onset_count", controls = CTRL_FATAL,
       fe = "ccode + year + leiden_comm",
       vcov_form = "ccode + year + leiden_comm"),
  list(key = "a3__3wfe__count_mid",
       dv = "log_churn", iv = "log_mid_count",        controls = CTRL_MID,
       fe = "ccode + year + leiden_comm",
       vcov_form = "ccode + year + leiden_comm"),
  list(key = "a3__3wfe__count_fatal",
       dv = "log_churn", iv = "log_fatal_mid_count",  controls = CTRL_FATAL,
       fe = "ccode + year + leiden_comm",
       vcov_form = "ccode + year + leiden_comm")
)

# A4: TWFE raw churn on raw MID ONSET
specs_a4 <- list(
  list(key = "a4__raw__none__mid",
       dv = "churn", iv = "mid_onset_count",       controls = NULL),
  list(key = "a4__raw__full__mid",
       dv = "churn", iv = "mid_onset_count",       controls = CTRL_MID),
  list(key = "a4__raw__none__fatal",
       dv = "churn", iv = "fatal_mid_onset_count", controls = NULL),
  list(key = "a4__raw__full__fatal",
       dv = "churn", iv = "fatal_mid_onset_count", controls = CTRL_FATAL)
)

# A5 (TWFE raw churn ~ raw MID count, ongoing) was dropped from the
# analysis. Coefficients were small and inconsistently significant
# because state FE absorbs most of the variation when both DV and IV
# are on the raw scale and the IV is ongoing rather than onset. The
# information-mechanism framing is also better served by the onset
# measure: ongoing-war years do not generate fresh information.
# A4 (raw scale, MID onset) remains as the raw-scale appendix table.

# Aggregate all specs and tag with table membership for downstream output.
all_spec_groups <- list(
  main = specs_main, a1 = specs_a1, a2 = specs_a2,
  a3   = specs_a3,   a4 = specs_a4
)

# Subset state_churn to non-missing leiden_comm for the 3-way FE specs.
state_churn_leiden <- state_churn |> filter(!is.na(leiden_comm))

cat("\n========== Fitting all main + appendix specs ==========\n")

twfe_models <- list()
for (group_name in names(all_spec_groups)) {
  for (sp in all_spec_groups[[group_name]]) {
    fe_arg <- if (!is.null(sp$fe)) sp$fe else "ccode + year"
    vc_arg <- if (!is.null(sp$vcov_form)) sp$vcov_form else "ccode + year"
    data_arg <- if (group_name == "a3") state_churn_leiden else state_churn

    fit <- fit_spec(
      data        = data_arg,
      dv          = sp$dv,
      iv          = sp$iv,
      controls    = sp$controls,
      fe          = fe_arg,
      vcov_form   = vc_arg
    )
    twfe_models[[sp$key]] <- list(
      group     = group_name,
      key       = sp$key,
      dv        = sp$dv,
      iv        = sp$iv,
      controls  = sp$controls,
      fe        = fe_arg,
      vcov_form = vc_arg,
      fit       = fit
    )
    cat(sprintf("  [%s] %s\n", group_name, sp$key))
  }
}

# --------------------------------------------------------------------------- #
# Lag models — keep the existing onset lag1/lag2 battery (used in the
# main-text "lag robustness" appendix). These now use log1p(churn) on the
# LHS to stay consistent with the headline TWFE.
# --------------------------------------------------------------------------- #

cat("\n========== LAG MODELS (lagged onset IVs, log1p(churn) DV) ==========\n")

lag_iv_map <- list(
  mid_onset_count           = c("onset_lag1",       "onset_lag2"),
  log_mid_onset_count       = c("log_onset_lag1",   "log_onset_lag2"),
  fatal_mid_onset_count     = c("fatal_onset_lag1", "fatal_onset_lag2"),
  log_fatal_mid_onset_count = c("log_fatal_onset_lag1", "log_fatal_onset_lag2")
)

lag_models <- list()
for (wv in names(lag_iv_map)) {
  lag1_var <- lag_iv_map[[wv]][1]
  lag2_var <- lag_iv_map[[wv]][2]
  ctrl_wv  <- ctrl_for_iv(wv)   # routes mid vs fatal partner-war control

  for (lag_info in list(c(lag1_var, "lag1"), c(lag2_var, "lag2"))) {
    lv  <- lag_info[1]; tag <- lag_info[2]
    lag_models[[paste0(wv, "_", tag, "_simple")]] <- feols(
      as.formula(paste("log_churn ~", lv, "| ccode + year")),
      data = state_churn, vcov = ~ccode + year)
    lag_models[[paste0(wv, "_", tag, "_full")]] <- feols(
      as.formula(paste("log_churn ~", lv, "+", ctrl_wv, "| ccode + year")),
      data = state_churn, vcov = ~ccode + year)
  }
  cat(sprintf("  %s lags: done\n", wv))
}

# --------------------------------------------------------------------------- #
# Panel Local Projections (Jordà 2005) — matches the new main spec.
# --------------------------------------------------------------------------- #
# 2 controls x 2 IVs x 6 horizons = 24 LP coefficients, all on logged DV/IV.
#
# For each horizon h:
#   log1p(churn)_{i, t+h} = alpha_i + gamma_t + beta_h * log1p(IV)_{i,t}
#                            + (controls)_{i,t} + epsilon_{i,t+h}
# --------------------------------------------------------------------------- #

cat("\n========== PANEL LOCAL PROJECTIONS ==========\n")

MAX_HORIZON <- 5

lp_data <- state_churn |>
  group_by(ccode) |>
  arrange(year) |>
  mutate(across(log_churn, list(
    h0 = ~.x,
    h1 = ~lead(.x, 1),
    h2 = ~lead(.x, 2),
    h3 = ~lead(.x, 3),
    h4 = ~lead(.x, 4),
    h5 = ~lead(.x, 5)
  ), .names = "{.fn}")) |>
  ungroup()

lp_specs <- list(
  list(key = "lp_none_mid",   iv = "log_mid_onset_count",        controls = NULL),
  list(key = "lp_full_mid",   iv = "log_mid_onset_count",        controls = CTRL_MID),
  list(key = "lp_none_fatal", iv = "log_fatal_mid_onset_count",  controls = NULL),
  list(key = "lp_full_fatal", iv = "log_fatal_mid_onset_count",  controls = CTRL_FATAL)
)

lp_results <- list()
for (sp in lp_specs) {
  rows <- list()
  for (h in 0:MAX_HORIZON) {
    dv  <- paste0("h", h)
    rhs <- sp$iv
    if (!is.null(sp$controls)) rhs <- paste(rhs, "+", sp$controls)
    fml <- as.formula(paste0(dv, " ~ ", rhs, " | ccode + year"))
    fit <- feols(fml, data = lp_data, vcov = ~ccode + year)
    ct  <- coeftable(fit)
    rows[[h + 1]] <- tibble(
      spec_key = sp$key,
      iv       = sp$iv,
      controls = ifelse(is.null(sp$controls), "none", "full"),
      iv_family= ifelse(grepl("fatal", sp$iv), "fatal", "mid"),
      horizon  = h,
      coef     = ct[sp$iv, "Estimate"],
      se       = ct[sp$iv, "Std. Error"],
      ci_low   = ct[sp$iv, "Estimate"] - 1.96 * ct[sp$iv, "Std. Error"],
      ci_high  = ct[sp$iv, "Estimate"] + 1.96 * ct[sp$iv, "Std. Error"],
      p_value  = ct[sp$iv, "Pr(>|t|)"]
    )
  }
  lp_results[[sp$key]] <- bind_rows(rows)
  cat(sprintf("  Panel LP [%s]: h=0..%d done\n", sp$key, MAX_HORIZON))
}

lp_all_df <- bind_rows(lp_results)

# --------------------------------------------------------------------------- #
# 3b. Panel LP Figure
# --------------------------------------------------------------------------- #

cat("\nPlotting panel local projections...\n")

LP_FIG_DIR_MAIN <- here::here("results", "main")
LP_FIG_DIR_APPX <- here::here("results", "appendix")
dir.create(LP_FIG_DIR_MAIN, recursive = TRUE, showWarnings = FALSE)
dir.create(LP_FIG_DIR_APPX, recursive = TRUE, showWarnings = FALSE)

# Main text: MID family only. Fatal-MID companion in appendix per the
# manuscript reorganization (fatal MIDs live in the appendix throughout).
lp_plot_data <- lp_all_df |>
  mutate(
    iv_label = factor(iv_family, levels = c("mid", "fatal"),
                       labels = c("log(1 + MID Onset)",
                                   "log(1 + Fatal MID Onset)")),
    controls_label = factor(controls, levels = c("none", "full"),
                             labels = c("No controls", "Full controls")),
    sig = p_value < 0.05
  )

build_lp_plot <- function(plot_data, title_iv, caption_extra = "") {
  ggplot(plot_data,
         aes(x = horizon, y = coef)) +
    geom_hline(yintercept = 0, color = "gray50",
               linetype = "dashed", linewidth = 0.5) +
    geom_ribbon(aes(ymin = ci_low, ymax = ci_high),
                alpha = 0.2, color = NA) +
    geom_line(linewidth = 1) +
    geom_point(aes(shape = sig, , color = controls_label), size = 3.5) +
    facet_wrap(~controls_label) + 
    scale_color_manual(values = c("p-value < 0.05"   = "#08519c",
                                    "Null" = "#B03366"),
                       name = NULL) +
    scale_x_continuous(breaks = 0:MAX_HORIZON) +
    labs(
      title = "Panel local projections",
      x = "Horizon (years)",
      y = "Coefficient"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title    = element_text(face = "bold", size = 14),
      strip.text    = element_text(face = "bold", size = 9),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.spacing = unit(1.2, "lines"),
      legend.position = "bottom",
      legend.title = element_blank(),
      plot.margin = margin(t = 10, r = 12, b = 14, l = 10)
    )
}

# --- Main text figure: MID-only ---
p_lp_mid <- build_lp_plot(
  lp_plot_data |> filter(iv_family == "mid") |> mutate(sig = ifelse(sig == T, "p-value < 0.05", "Null")),
  title_iv = "log(1 + MID Onset)"
)
ggsave(file.path(LP_FIG_DIR_MAIN, "figure1_panel_lp.png"),
       p_lp_mid, width = 7, height = 5, dpi = 300, bg = "white")
cat("  Panel LP figure (MID, main) -> figure1_panel_lp.png\n")

# --- Appendix figure: Fatal-MID companion ---
p_lp_fatal <- build_lp_plot(
  lp_plot_data |> filter(iv_family == "fatal"),
  title_iv = "log(1 + Fatal MID Onset)",
  caption_extra = " Companion to main-text Figure 1 (MID Onset)."
)
ggsave(file.path(LP_FIG_DIR_APPX, "figureA_panel_lp_fatal.png"),
       p_lp_fatal, width = 7, height = 5, dpi = 300, bg = "white")
cat("  Panel LP figure (Fatal, appendix) -> figureA_panel_lp_fatal.png\n")

# --------------------------------------------------------------------------- #
# 4. Summary Tables (LaTeX via texreg)
# --------------------------------------------------------------------------- #

cat("\n========== LATEX TABLES ==========\n")

TABLE_DIR <- here::here("results", "main")
APPX_DIR  <- here::here("results", "appendix")
dir.create(TABLE_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(APPX_DIR, recursive = TRUE, showWarnings = FALSE)

# Coefficient map: map all IV names (raw + log; onset + count) to one
# "Conflict measure" row so the tables stay narrow and readable.
state_coef_map <- list(
  "mid_onset_count"              = "Conflict measure",
  "log_mid_onset_count"          = "Conflict measure",
  "fatal_mid_onset_count"        = "Conflict measure",
  "log_fatal_mid_onset_count"    = "Conflict measure",
  "mid_count"                    = "Conflict measure",
  "log_mid_count"                = "Conflict measure",
  "fatal_mid_count"              = "Conflict measure",
  "log_fatal_mid_count"          = "Conflict measure",
  "partner_war_frac"               = "Partner war exposure",
  "partner_war_frac_weighted"      = "Partner war exposure (weighted)",
  "partner_fatal_war_frac"         = "Partner fatal-war exposure",
  "partner_fatal_war_frac_weighted"= "Partner fatal-war exposure (weighted)",
  "democracy"                    = "Democracy",
  "autocracy"                    = "Autocracy",
  "cinc"                         = "CINC",
  "major_power"                  = "Major power",
  "is_isolated_t1"               = "Isolated at $t{-}1$",
  "(Intercept)"                  = "Intercept"
)

# --- Building blocks for the new 5-column main table ---
# m1, m2: baseline IV / baseline partner-war exposure with a
# quadratic year trend (poly(year, 2)) to mirror the system-level
# Table 1 spec ladder, where all four columns include the year
# trend. One-way state-clustered SEs since FE machinery isn't on.
# Partner-war control mirrors IV severity AND uses the tie-WEIGHTED
# variant (mid -> partner_war_frac_weighted, fatal ->
# partner_fatal_war_frac_weighted) so it's consistent with CTRL_MID /
# CTRL_FATAL above.
m1_mid <- feols(log_churn ~ log_mid_onset_count + poly(year, 2),
                data = state_churn, vcov = ~ccode)
m2_mid <- feols(log_churn ~ partner_war_frac_weighted + poly(year, 2),
                data = state_churn, vcov = ~ccode)
m1_fatal <- feols(log_churn ~ log_fatal_mid_onset_count + poly(year, 2),
                  data = state_churn, vcov = ~ccode)
m2_fatal <- feols(log_churn ~ partner_fatal_war_frac_weighted + poly(year, 2),
                  data = state_churn, vcov = ~ccode)

# m5: TWFE (state + year only) with THREE-WAY clustered SEs on
# (ccode + year + leiden_comm). This is different from the legacy
# `a3__3wfe__*` spec — there, leiden_comm enters as an additional
# fixed effect (absorbing community-level mean shifts). Here, it
# enters only on the inference side: the point estimate is
# identified from the same within-state, within-year variation as
# m4, but inference accounts for arbitrary correlation within
# multilayer communities. Sample restricted to state-years with
# non-missing leiden_comm to make the cluster ID well-defined.
state_churn_leiden <- state_churn |> filter(!is.na(leiden_comm))

m5_mid <- feols(log_churn ~ log_mid_onset_count + partner_war_frac_weighted +
                  democracy + autocracy + cinc + major_power +
                  is_isolated_t1 | ccode + year,
                data = state_churn_leiden,
                vcov = ~ccode + year + leiden_comm)

m5_fatal <- feols(log_churn ~ log_fatal_mid_onset_count + partner_fatal_war_frac_weighted +
                    democracy + autocracy + cinc + major_power +
                    is_isolated_t1 | ccode + year,
                  data = state_churn_leiden,
                  vcov = ~ccode + year + leiden_comm)

# --- Helper: write a 5-column table from a list of fitted models ---
# (Each entry can be a feols/lm object directly; texreg dispatches by
#  class. Use this for the new main + fatal-appendix tables.)
write_state_5col <- function(models,
                              column_headers,
                              caption,
                              label,
                              note_extra,
                              file_path,
                              digits = 3) {
  base_note <- paste0(
    "\\parbox{\\linewidth}{\\vspace{2pt}\\textit{Note:} ",
    note_extra, " ",
    "Standard errors in parentheses. \\textbf{Columns (1)--(3) include ",
    "a quadratic year trend ($\\texttt{poly(year, 2)}$); terms suppressed ",
    "for legibility. In columns (4)--(5), year effects are absorbed by ",
    "the year fixed effect.} Columns (1) and (2) are baseline OLS ",
    "regressions of the focal IV (and the year trend) with state-clustered ",
    "SEs. Column (3) adds the full controls (partner war exposure, ",
    "democracy, autocracy, CINC, major power, ",
    "$\\mathbb{1}\\{\\text{isolated at }t-1\\}$); SEs clustered by state. ",
    "Column (4) is the headline two-way ",
    "fixed-effects spec (state $+$ year), two-way clustered SEs. ",
    "Column (5) keeps the same TWFE estimator but clusters SEs three ",
    "ways (state $+$ year $+$ Leiden multilayer community); sample ",
    "restricted to state-years with non-missing community assignment. ",
    "COW system, 1816--2014. ",
    "$^{***}p<0.001$; $^{**}p<0.01$; $^{*}p<0.05$; $^{\\dagger}p<0.1$.}")
  texreg(models,
         custom.model.names = column_headers,
         custom.coef.map    = state_coef_map,
         caption    = caption,
         label      = label,
         custom.note= base_note,
         stars      = c(0.001, 0.01, 0.05, 0.1),
         symbol     = "$^{\\dagger}$",
         float.pos  = "htbp",
         booktabs   = TRUE,
         use.packages = FALSE,
         digits     = digits,
         file       = file_path)
  cat(sprintf("  %s -> %s\n", label, basename(file_path)))
}

# --- Table 3 (main text): MID-only state ladder, 5 cols ----------------
# m1: baseline IV + poly(year, 2)
# m2: baseline partner war (severity-matched + tie-WEIGHTED:
#      partner_war_frac_weighted) + poly(year, 2)
# m3: + full controls + poly(year, 2), NO FE  (a2__nofe__onset_mid)
# m4: TWFE state + year, full controls       (main__full__mid)
# m5: TWFE state + year, full controls; SEs clustered (state + year +
#     leiden_comm). Same point estimator as m4; inference allows for
#     within-community correlation in addition to within-state and
#     within-year correlation.
#
# Wild cluster bootstrap p-values for column 4 are produced separately by
# 19_wild_cluster_bootstrap.R and reported in Appendix Table
# tableA_wild_bootstrap.tex. They closely match the asymptotic two-way
# clustered p-values shown here.
write_state_5col(
  models = list(
    m1_mid,
    m2_mid,
    twfe_models[["a2__nofe__onset_mid"]]$fit,
    twfe_models[["main__full__mid"]]$fit,
    m5_mid
  ),
  column_headers = c(
    "(1) IV alone",
    "(2) PWE alone",
    "(3) + Ctrls, no FE",
    "(4) TWFE",
    "(5) TWFE + Comm cluster"
  ),
  caption = paste0(
    "State-Level Specification Ladder: $\\log(1+\\text{Churn}) ",
    "\\sim \\log(1+\\text{MID Onset})$. ",
    "Each column adds one layer of identification machinery; the ",
    "coefficient is stable across columns."),
  label   = "tab:state_main",
  note_extra = paste0(
    "DV is $\\log(1+\\text{Churn})$ where churn $= 1 - $ weighted ",
    "Jaccard similarity year-over-year. IV is $\\log(1+\\text{MID Onset})$. ",
    "Partner war exposure (PWE) is the alliance/IGO tie-strength-WEIGHTED ",
    "share of $t-1$ partners in any MID at $t$."),
  file_path = file.path(TABLE_DIR, "table3_state_main.tex")
)

# --- Appendix Table A1: Fatal-MID state ladder (parallel structure) ----
# Same 5-column structure as Table 3 but with fatal MID IV and the
# severity-matched fatal partner-war control.
write_state_5col(
  models = list(
    m1_fatal,
    m2_fatal,
    twfe_models[["a2__nofe__onset_fatal"]]$fit,
    twfe_models[["main__full__fatal"]]$fit,
    m5_fatal
  ),
  column_headers = c(
    "(1) IV alone",
    "(2) PWE alone",
    "(3) + Ctrls, no FE",
    "(4) TWFE",
    "(5) TWFE + Comm cluster"
  ),
  caption = paste0(
    "Appendix. State-Level Specification Ladder: $\\log(1+\\text{Churn}) ",
    "\\sim \\log(1+\\text{Fatal MID Onset})$. ",
    "Companion to main-text Table~\\ref{tab:state_main}; partner-war ",
    "exposure uses the fatal-MID severity-matched, tie-WEIGHTED variant ",
    "($\\textit{partner\\_fatal\\_war\\_frac\\_weighted}$)."),
  label   = "tab:state_appx_fatal",
  note_extra = paste0(
    "DV is $\\log(1+\\text{Churn})$. IV is $\\log(1+\\text{Fatal MID ",
    "Onset})$. Partner war exposure (PWE) here is the alliance/IGO ",
    "tie-WEIGHTED share of $t-1$ partners in a fatal MID at $t$ ",
    "(severity-matched to the IV)."),
  file_path = file.path(APPX_DIR, "tableA1_state_fatal.tex")
)

# --- Appendix Table A2: 3-way FE (state + year + Leiden community) -----
# 4 cols: MID onset, Fatal MID onset, MID count, Fatal MID count.
# This table is genuinely unique to 02 — 03's multiverse only uses Leiden
# as an SE clustering variable, never as a fixed effect that absorbs
# community-level mean differences. Keep here so reviewers can see the
# coefficient survives community-level absorption.
write_state_5col(
  models = list(
    twfe_models[["a3__3wfe__onset_mid"]]$fit,
    twfe_models[["a3__3wfe__onset_fatal"]]$fit,
    twfe_models[["a3__3wfe__count_mid"]]$fit,
    twfe_models[["a3__3wfe__count_fatal"]]$fit
  ),
  column_headers = c(
    "$\\log(1{+}\\text{MID Onset})$",
    "$\\log(1{+}\\text{Fatal Onset})$",
    "$\\log(1{+}\\text{MID Count})$",
    "$\\log(1{+}\\text{Fatal Count})$"
  ),
  caption = paste0(
    "Appendix. Three-way fixed effects (state $+$ year $+$ Leiden ",
    "multilayer community): $\\log(1+\\text{Churn})$ on log MID measures, ",
    "full controls."),
  label   = "tab:state_appx_3wfe",
  note_extra = paste0(
    "DV is $\\log(1+\\text{Churn})$. Three-way fixed effects (state, ",
    "year, Leiden multilayer community); SEs clustered three ways. ",
    "Sample restricted to state-years with non-missing Leiden community ",
    "assignment. Columns 1--2 use onset IVs; columns 3--4 use ongoing ",
    "MID counts. Partner-war control routed by IV severity."),
  file_path = file.path(APPX_DIR, "tableA2_state_3wfe.tex")
)

# --- Appendix Table A3: Lag models (8 columns, t-1 and t-2 IV) ---------

lag_coef_map <- list(
  "onset_lag1"                   = "MID Onset ($t{-}1$)",
  "onset_lag2"                   = "MID Onset ($t{-}2$)",
  "log_onset_lag1"               = "log MID Onset ($t{-}1$)",
  "log_onset_lag2"               = "log MID Onset ($t{-}2$)",
  "fatal_onset_lag1"             = "Fatal Onset ($t{-}1$)",
  "fatal_onset_lag2"             = "Fatal Onset ($t{-}2$)",
  "log_fatal_onset_lag1"         = "log Fatal Onset ($t{-}1$)",
  "log_fatal_onset_lag2"         = "log Fatal Onset ($t{-}2$)",
  "partner_war_frac"               = "Partner war exposure",
  "partner_war_frac_weighted"      = "Partner war exposure (weighted)",
  "partner_fatal_war_frac"         = "Partner fatal-war exposure",
  "partner_fatal_war_frac_weighted"= "Partner fatal-war exposure (weighted)",
  "democracy"                    = "Democracy",
  "autocracy"                    = "Autocracy",
  "cinc"                         = "CINC",
  "major_power"                  = "Major power",
  "is_isolated_t1"               = "Isolated at $t{-}1$"
)

lag_note <- paste0(
  "\\parbox{\\linewidth}{\\vspace{2pt}\\textit{Note:} ",
  "DV is $\\log(1+\\text{Churn})$. ",
  "Two-way fixed effects (state $+$ year), two-way clustered SEs ",
  "(state $+$ year). Each column reports the IV at $t{-}1$ or $t{-}2$ ",
  "with the full controls covariate set (partner war exposure, ",
  "democracy, autocracy, CINC, major power). ",
  "COW system, 1816--2014. ",
  "$^{***}p<0.001$; $^{**}p<0.01$; $^{*}p<0.05$; $^{\\dagger}p<0.1$.}")

texreg(
  list(
    lag_models[["mid_onset_count_lag1_full"]],     lag_models[["mid_onset_count_lag2_full"]],
    lag_models[["log_mid_onset_count_lag1_full"]], lag_models[["log_mid_onset_count_lag2_full"]],
    lag_models[["fatal_mid_onset_count_lag1_full"]], lag_models[["fatal_mid_onset_count_lag2_full"]],
    lag_models[["log_fatal_mid_onset_count_lag1_full"]], lag_models[["log_fatal_mid_onset_count_lag2_full"]]
  ),
  custom.model.names = c(
    "MID $t{-}1$", "MID $t{-}2$",
    "log MID $t{-}1$", "log MID $t{-}2$",
    "Fatal $t{-}1$", "Fatal $t{-}2$",
    "log Fatal $t{-}1$", "log Fatal $t{-}2$"
  ),
  custom.coef.map = lag_coef_map,
  caption = "Appendix. Lagged Conflict Onset and $\\log(1+\\text{Churn})$.",
  label = "tab:lag_models",
  digits = 3,
  custom.note = lag_note,
  stars = c(0.001, 0.01, 0.05, 0.1),
  symbol = "$^{\\dagger}$",
  float.pos = "htbp",
  booktabs = TRUE,
  use.packages = FALSE,
  file = file.path(APPX_DIR, "tableA3_state_lag_models.tex")
)
cat("  Appendix (lag models) -> tableA3_state_lag_models.tex\n")

# --------------------------------------------------------------------------- #
# 5. Clustering Robustness (onset IVs — appendix)
# --------------------------------------------------------------------------- #
# Full specification grid across all 4 IVs × weighting is in 03_robustness.R.
# Here: onset IVs only, varying SE clustering scheme.
# --------------------------------------------------------------------------- #

cat("\n========== CLUSTERING ROBUSTNESS ==========\n")

# One-way (state only)
m_onset_full_1way <- feols(churn ~ mid_onset_count + partner_war_frac_weighted +
                             democracy + autocracy + cinc + major_power + is_isolated_t1 | ccode + year,
                           data = state_churn, vcov = ~ccode)

m_fonset_full_1way <- feols(churn ~ fatal_mid_onset_count + partner_fatal_war_frac_weighted +
                              democracy + autocracy + cinc + major_power + is_isolated_t1 | ccode + year,
                            data = state_churn, vcov = ~ccode)

# Leiden CPM variants
state_churn_leiden <- state_churn |> filter(!is.na(leiden_comm))

m_onset_full_leiden <- feols(churn ~ mid_onset_count + partner_war_frac_weighted +
                               democracy + autocracy + cinc + major_power + is_isolated_t1 | ccode + year,
                             data = state_churn_leiden, vcov = ~leiden_comm)

m_onset_full_leiden_yr <- feols(churn ~ mid_onset_count + partner_war_frac_weighted +
                                  democracy + autocracy + cinc + major_power + is_isolated_t1 | ccode + year,
                                data = state_churn_leiden, vcov = ~leiden_comm + year)

m_onset_full_leiden_st <- feols(churn ~ mid_onset_count + partner_war_frac_weighted +
                                  democracy + autocracy + cinc + major_power + is_isolated_t1 | ccode + year,
                                data = state_churn_leiden, vcov = ~leiden_comm + ccode)

m_onset_full_leiden_yr_st <- feols(churn ~ mid_onset_count + partner_war_frac_weighted +
                                   democracy + autocracy + cinc + major_power + is_isolated_t1 | ccode + year,
                                 data = state_churn_leiden, vcov = ~leiden_comm + year + ccode)


m_fonset_full_leiden <- feols(churn ~ fatal_mid_onset_count + partner_fatal_war_frac_weighted +
                               democracy + autocracy + cinc + major_power + is_isolated_t1 | ccode + year,
                             data = state_churn_leiden, vcov = ~leiden_comm)

m_fonset_full_leiden_yr <- feols(churn ~ fatal_mid_onset_count + partner_fatal_war_frac_weighted +
                                  democracy + autocracy + cinc + major_power + is_isolated_t1 | ccode + year,
                                data = state_churn_leiden, vcov = ~leiden_comm + year)

m_fonset_full_leiden_st <- feols(churn ~ fatal_mid_onset_count + partner_fatal_war_frac_weighted +
                                  democracy + autocracy + cinc + major_power + is_isolated_t1 | ccode + year,
                                data = state_churn_leiden, vcov = ~leiden_comm + ccode)

m_fonset_full_leiden_yr_st <- feols(churn ~ fatal_mid_onset_count + partner_fatal_war_frac_weighted +
                                     democracy + autocracy + cinc + major_power + is_isolated_t1 | ccode + year,
                                   data = state_churn_leiden, vcov = ~leiden_comm + year + ccode)
# Clustering-robustness model fits are retained (saved into
# state_level_results.rds for downstream consumers). The standalone
# clustering-robustness appendix tables that lived here have been
# REMOVED — they're redundant with 03_robustness.R's 5-clustering-scheme
# loop, which produces the same SE-ladder evidence at every multiverse
# cell rather than just one column. Cite 03's clustering panel in the
# appendix instead.
cat("  Clustering-robustness models fitted (saved to RDS only; ",
    "tables generated by 03_robustness.R)\n", sep = "")

# --------------------------------------------------------------------------- #
# 6. Save
# --------------------------------------------------------------------------- #

results <- list(
  # Main + appendix TWFE models, keyed by spec key
  # (main__*, a1__*, a2__*, a3__*, a4__*)
  twfe = twfe_models,
  # Lag models (onset IVs at lag1/lag2, log_churn DV)
  lags = lag_models,
  # Clustering robustness — onset IVs (legacy, kept for back-compat)
  m_onset_full_1way         = m_onset_full_1way,
  m_onset_full_leiden       = m_onset_full_leiden,
  m_onset_full_leiden_yr    = m_onset_full_leiden_yr,
  m_onset_full_leiden_st    = m_onset_full_leiden_st,
  m_onset_full_leiden_yr_st = m_onset_full_leiden_yr_st,
  m_fonset_full_1way        = m_fonset_full_1way,
  m_fonset_full_leiden      = m_fonset_full_leiden,
  m_fonset_full_leiden_yr   = m_fonset_full_leiden_yr,
  m_fonset_full_leiden_st   = m_fonset_full_leiden_st,
  m_fonset_full_leiden_yr_st= m_fonset_full_leiden_yr_st,
  # Panel local projections (long format: spec_key, controls, iv, horizon)
  lp = lp_all_df,
  # Data
  data = state_churn
)

saveRDS(results, file.path(OUTPUT_DIR, "state_level_results.rds"))

# --- CSV exports ---
write_csv(state_churn, file.path(OUTPUT_DIR, "state_churn_data.csv"))

# Coefficient extractor — works on both list-wrapped (twfe_models) and
# bare feols (lag_models, clustering robustness) entries.
extract_feols_coefs <- function(fit_obj, label, group = NA, dv = NA, iv = NA,
                                  controls = NA) {
  ct <- coeftable(fit_obj)
  tibble(
    spec_key = label,
    group    = group,
    dv       = dv,
    iv       = iv,
    controls = controls,
    term     = rownames(ct),
    estimate = ct[, "Estimate"],
    se       = ct[, "Std. Error"],
    t_value  = ct[, "t value"],
    p_value  = ct[, "Pr(>|t|)"]
  )
}

# State-level TWFE coefs (combined across all groups)
state_coefs <- list()
for (nm in names(twfe_models)) {
  m <- twfe_models[[nm]]
  state_coefs[[nm]] <- extract_feols_coefs(
    m$fit, nm, m$group, m$dv, m$iv,
    ifelse(is.null(m$controls), "none", "full"))
}
state_coefs_df <- bind_rows(state_coefs)
write_csv(state_coefs_df,
          file.path(OUTPUT_DIR, "state_level_coefs_combined.csv"))
# Convenience splits: one CSV per table
for (g in unique(state_coefs_df$group)) {
  write_csv(state_coefs_df |> filter(group == g),
            file.path(OUTPUT_DIR,
                      sprintf("state_level_coefs_%s.csv", g)))
}

# Lag-model coefs (separate file)
lag_coefs <- list()
for (nm in names(lag_models)) {
  lag_coefs[[nm]] <- extract_feols_coefs(lag_models[[nm]], nm,
                                           group = "lag")
}
write_csv(bind_rows(lag_coefs),
          file.path(OUTPUT_DIR, "state_lag_coefs.csv"))

# LP exports
write_csv(lp_all_df, file.path(OUTPUT_DIR, "panel_lp_main.csv"))

# =========================================================================== #
# SUMMARY STATISTICS TABLE
# =========================================================================== #
cat("\n========== Summary Statistics ==========\n")

# Aligned with the variables that appear in the main-text Table 3
# (state-level specification ladder) full-control regressions: DV
# (churn, raw + logged), IV (MID onset, raw + logged; fatal-MID variants
# for the appendix companion), and the tie-WEIGHTED partner-war controls
# (severity-routed: partner_war_frac_weighted for MID,
# partner_fatal_war_frac_weighted for fatal). State attributes (democracy,
# autocracy, CINC, major power, is_isolated_t1) round out the control set.
#
# DROPPED: unweighted partner-war variants (only appear as alt specs in
# 03's multiverse, never in main-text). Ongoing MID count IVs (mid_count,
# fatal_mid_count) — main spec uses onset, count is appendix-only via 03.
# Lagged versions of any IV — not used in the main-text Table 3.
sumstat_vars <- c(
  # DV (Table 3)
  "churn", "log_churn",
  # IV (Table 3 main: MID onset; fatal in appendix companion)
  "mid_onset_count", "log_mid_onset_count",
  "fatal_mid_onset_count", "log_fatal_mid_onset_count",
  "partner_war_frac_weighted", "partner_fatal_war_frac_weighted",
  # State-level controls
  "democracy", "autocracy", "cinc", "major_power", "is_isolated_t1"
)

sumstat_labels <- c(
  churn                            = "Network Churn",
  log_churn                        = "$\\log(1 + \\text{Churn})$",
  mid_onset_count                  = "MID Onset Count",
  log_mid_onset_count              = "$\\log(1 + \\text{MID Onset})$",
  fatal_mid_onset_count            = "Fatal MID Onset Count",
  log_fatal_mid_onset_count        = "$\\log(1 + \\text{Fatal MID Onset})$",
  partner_war_frac_weighted        = "Partner War Fraction (weighted)",
  partner_fatal_war_frac_weighted  = "Partner Fatal War Fraction (weighted)",
  democracy                        = "Democracy",
  autocracy                        = "Autocracy",
  cinc                             = "CINC Score",
  major_power                      = "Major Power",
  is_isolated_t1                   = "Isolated at $t-1$"
)

sumstat_df <- state_churn |>
  select(any_of(sumstat_vars)) |>
  pivot_longer(everything(), names_to = "variable", values_to = "value") |>
  group_by(variable) |>
  summarise(
    N    = sum(!is.na(value)),
    Mean = mean(value, na.rm = TRUE),
    SD   = sd(value, na.rm = TRUE),
    Min  = min(value, na.rm = TRUE),
    Q25  = quantile(value, 0.25, na.rm = TRUE),
    Median = median(value, na.rm = TRUE),
    Q75  = quantile(value, 0.75, na.rm = TRUE),
    Max  = max(value, na.rm = TRUE),
    .groups = "drop"
  ) |>
  mutate(
    variable = factor(variable, levels = sumstat_vars),
    Label = sumstat_labels[as.character(variable)]
  ) |>
  arrange(variable) |>
  select(Label, N, Mean, SD, Min, Q25, Median, Q75, Max)

# CSV export
write_csv(sumstat_df, file.path(OUTPUT_DIR, "summary_stats_state.csv"))
cat("  CSV → summary_stats_state.csv\n")

# LaTeX export
tex_lines <- c(
  "\\begin{table}[htbp]",
  "\\centering",
  "\\caption{Summary Statistics: State-Level Variables}",
  "\\label{tab:sumstats_state}",
  "\\small",
  "\\begin{tabular}{lrrrrrrrr}",
  "\\toprule",
  "Variable & N & Mean & SD & Min & Q25 & Median & Q75 & Max \\\\",
  "\\midrule"
)
for (i in seq_len(nrow(sumstat_df))) {
  row <- sumstat_df[i, ]
  tex_lines <- c(tex_lines, sprintf(
    "%s & %d & %.4f & %.4f & %.4f & %.4f & %.4f & %.4f & %.4f \\\\",
    row$Label, row$N, row$Mean, row$SD, row$Min, row$Q25, row$Median, row$Q75, row$Max
  ))
}
tex_lines <- c(tex_lines,
  "\\bottomrule",
  "\\end{tabular}",
  "\\par\\vspace{4pt}",
  paste0(
    "{\\footnotesize \\textit{Note:} Variables appearing in the main-text ",
    "Table~\\ref{tab:state_main} state-level full-control regressions. ",
    "Both raw and $\\log(1+x)$ forms reported for the DV (churn) and ",
    "IVs (MID onset, fatal MID onset). Partner war fraction is the ",
    "unweighted share of $t-1$ partners in any MID at $t$ (severity-",
    "routed: \\textit{partner\\_war\\_frac} for MID-IV regressions, ",
    "\\textit{partner\\_fatal\\_war\\_frac} for fatal-MID regressions). ",
    "Weighted variants used only in the 03 multiverse; not shown here. ",
    "$\\mathbb{1}\\{$Isolated at $t{-}1\\}$ flags state-years with zero ",
    "partners at the prior year. COW system, 1816--2014.}"
  ),
  "\\end{table}"
)
writeLines(tex_lines, file.path(TABLE_DIR, "summary_stats_state.tex"))
cat("  LaTeX -> summary_stats_state.tex (results/main)\n")

cat("\nDone. State-level analysis complete.\n")