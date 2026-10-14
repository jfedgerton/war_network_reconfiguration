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
#   democracy, autocracy, cinc, major_power, is_isolated_t1
# [2026-10-12] Partner war exposure (PWE) is no longer a control in any model
# in this script (H3 dropped from the paper). The PWE variables are still
# computed (section 0b / 0b2) for the internal checks, but they are dropped
# from the saved panels (state_level_results.rds, state_churn_data.csv).
###############################################################################

library(tidyverse)
library(fixest)
library(peacesciencer)
library(igraph)
library(dynamicmultiplex)
library(texreg)

select <- dplyr::select
set.seed(123)

cat(sprintf("R %s | dynamicmultiplex %s | igraph %s | fixest %s\n",
            getRversion(), packageVersion("dynamicmultiplex"),
            packageVersion("igraph"), packageVersion("fixest")))
stopifnot(packageVersion("dynamicmultiplex") >= "1.4.0")   # fit_multilayer_jaccard(tracking, repair_chains)

OUTPUT_DIR <- here::here("data", "processed")

# --------------------------------------------------------------------------- #
# 0. Load data
# --------------------------------------------------------------------------- #

# Inputs must come from the same (most recent) Stage 1 run: print their mtimes
# and stop if they were written more than 6 hours apart.
INPUT_FILES <- file.path(OUTPUT_DIR, c("state_network_churn.rds", "dyad_year_network_layers.rds"))
stopifnot(all(file.exists(INPUT_FILES)))
input_mtime <- file.mtime(INPUT_FILES)
for (i in seq_along(INPUT_FILES)) cat(sprintf("  input %-32s mtime %s\n", basename(INPUT_FILES[i]), format(input_mtime[i], "%Y-%m-%d %H:%M:%S")))
if (diff(range(input_mtime)) > as.difftime(6, units = "hours"))
  stop("input files were written more than 6 hours apart; they are not from one Stage 1 run")

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
# 0b2. Peace-dyad churn and partner war exposure (Table 3, Panel B)
# --------------------------------------------------------------------------- #
# Partner j is dropped from state i's portfolio in BOTH years (t-1 and t) if i
# and j are in an ongoing MID (cowmidongoing) in t-1 or t. Churn is then the
# same weighted Jaccard as the canonical measure (02_compute_network_measures.R),
# computed on the remaining partners. If removing those partners empties a
# portfolio that was non-empty, the state-year is set to NA (no peace dyads to
# measure). Partner war exposure drops the same partners from numerator and
# denominator. Same rule as 18_robustness_excl_belligerent.R (validated against
# its per-state-year function on 2,500 random state-years, 2026-10-01).
cat("Computing peace-dyad churn and partner war exposure (MID in t-1 or t)...\n")
opp_pairs <- dyad_data |>
  filter(cowmidongoing == 1) |>
  distinct(ccode1, ccode2, year)
opp_pairs <- bind_rows(opp_pairs, opp_pairs |> rename(ccode1 = ccode2, ccode2 = ccode1)) |>
  distinct() |>
  transmute(i = ccode1, j = ccode2, year, opp = 1L)

# directed tie list, one row per (i, j, year); max over the two stored directions (as in the canonical churn)
peace_ties <- dyad_data |>
  select(ccode1, ccode2, year, w = !!sym(MAIN_WEIGHT_VAR)) |>
  bind_rows(dyad_data |> select(ccode1 = ccode2, ccode2 = ccode1, year, w = !!sym(MAIN_WEIGHT_VAR))) |>
  filter(w > 0) |>
  group_by(i = ccode1, j = ccode2, year) |>
  summarise(w = max(w), .groups = "drop")

peace_pairs <- full_join(peace_ties |> transmute(i, j, year = year + 1L, w1 = w),
                         peace_ties |> transmute(i, j, year, w2 = w),
                         by = c("i", "j", "year")) |>
  mutate(w1 = replace_na(w1, 0), w2 = replace_na(w2, 0)) |>
  left_join(opp_pairs |> rename(opp_t = opp), by = c("i", "j", "year")) |>
  left_join(opp_pairs |> transmute(i, j, year = year + 1L, opp_t1 = opp), by = c("i", "j", "year")) |>
  mutate(drop = !is.na(opp_t) | !is.na(opp_t1))

peace_churn <- peace_pairs |>
  group_by(ccode = i, year) |>
  summarise(
    nf1 = sum(w1 > 0), nf2 = sum(w2 > 0),
    n1 = sum(w1 > 0 & !drop), n2 = sum(w2 > 0 & !drop),
    all_min = sum(pmin(w1, w2)), all_max = sum(pmax(w1, w2)),
    pc_min = sum(pmin(w1, w2)[!drop]), pc_max = sum(pmax(w1, w2)[!drop]),
    .groups = "drop") |>
  mutate(
    churn_check = case_when(nf1 == 0 & nf2 == 0 ~ 0, nf1 == 0 | nf2 == 0 ~ 1, TRUE ~ 1 - all_min / all_max),
    churn_peace = case_when((n2 == 0 & nf2 > 0) | (n1 == 0 & nf1 > 0) ~ NA_real_,
                            n1 == 0 & n2 == 0 ~ 0, n1 == 0 | n2 == 0 ~ 1,
                            TRUE ~ 1 - pc_min / pc_max)) |>
  select(ccode, year, churn_check, churn_peace)

state_churn <- state_churn |>
  left_join(peace_churn, by = c("ccode", "year")) |>
  mutate(churn_check = ifelse(is.na(churn_check) & !is.na(churn), 0, churn_check),    # no ties in either year
         churn_peace = ifelse(!paste(ccode, year) %in% paste(peace_churn$ccode, peace_churn$year) & !is.na(churn), 0, churn_peace))
# The no-exclusion recomputation must reproduce the canonical churn; otherwise the
# dyad file and state_network_churn.rds come from different runs of Stage 1.
n_bad <- sum(abs(state_churn$churn_check - state_churn$churn) > 1e-6, na.rm = TRUE)
if (n_bad > 0) stop(n_bad, " state-years: churn recomputed from dyad_year_network_layers.rds != state_network_churn.rds; rerun Stage 1")
state_churn$churn_check <- NULL
cat(sprintf("  Peace-dyad churn: %d state-years undefined (every partner was an MID opponent)\n",
            sum(is.na(state_churn$churn_peace) & !is.na(state_churn$churn))))

peace_exposure <- partner_sets |>
  mutate(year_next = year + 1L) |>
  left_join(opp_pairs |> transmute(ccode1 = i, ccode2 = j, year, opp_t1 = opp), by = c("ccode1", "ccode2", "year")) |>
  left_join(opp_pairs |> transmute(ccode1 = i, ccode2 = j, year_next = year, opp_t = opp), by = c("ccode1", "ccode2", "year_next")) |>
  filter(is.na(opp_t1), is.na(opp_t)) |>
  left_join(all_state_war, by = c("ccode2" = "ccode", "year_next" = "year")) |>
  mutate(at_war = replace_na(at_war, 0L)) |>
  group_by(ccode = ccode1, year = year_next) |>
  summarise(partner_war_frac_weighted_peace = sum(at_war * w) / sum(w), .groups = "drop")

state_churn <- state_churn |>
  left_join(peace_exposure, by = c("ccode", "year")) |>
  mutate(partner_war_frac_weighted_peace = replace_na(partner_war_frac_weighted_peace, 0),
         log_churn_peace = log1p(churn_peace))
stopifnot(all(state_churn$partner_war_frac_weighted_peace >= 0 & state_churn$partner_war_frac_weighted_peace <= 1))

# [2026-10-12] HEADLINE peace rule (main-text Table 3, Panel B): partner j is dropped from state i's
# portfolio in BOTH years if i and j were opponents in a FATAL MID ongoing in t-1 (and only t-1). The
# treatment is fatal MID onset at t-1, so this removes exactly the treated dyads and fatal disputes still
# running in t-1. Same churn and isolation rules as churn_peace above. Result: churn_peace_ft1.
opp_pairs_f <- dyad_data |>
  filter(cowmidongoing == 1, cowmidfatal == 1) |>
  distinct(ccode1, ccode2, year)
opp_pairs_f <- bind_rows(opp_pairs_f, opp_pairs_f |> rename(ccode1 = ccode2, ccode2 = ccode1)) |>
  distinct() |>
  transmute(i = ccode1, j = ccode2, year, opp = 1L)
peace_pairs_f <- full_join(peace_ties |> transmute(i, j, year = year + 1L, w1 = w),
                           peace_ties |> transmute(i, j, year, w2 = w),
                           by = c("i", "j", "year")) |>
  mutate(w1 = replace_na(w1, 0), w2 = replace_na(w2, 0)) |>
  left_join(opp_pairs_f |> transmute(i, j, year = year + 1L, opp_t1 = opp), by = c("i", "j", "year")) |>
  mutate(drop = !is.na(opp_t1))
peace_churn_f <- peace_pairs_f |>
  group_by(ccode = i, year) |>
  summarise(
    nf1 = sum(w1 > 0), nf2 = sum(w2 > 0),
    n1 = sum(w1 > 0 & !drop), n2 = sum(w2 > 0 & !drop), n_drop = sum(drop & (w1 > 0 | w2 > 0)),
    pc_min = sum(pmin(w1, w2)[!drop]), pc_max = sum(pmax(w1, w2)[!drop]),
    .groups = "drop") |>
  mutate(churn_peace_ft1 = case_when((n2 == 0 & nf2 > 0) | (n1 == 0 & nf1 > 0) ~ NA_real_,
                                     n1 == 0 & n2 == 0 ~ 0, n1 == 0 | n2 == 0 ~ 1,
                                     TRUE ~ 1 - pc_min / pc_max)) |>
  select(ccode, year, churn_peace_ft1, n_drop_ft1 = n_drop)
state_churn <- state_churn |>
  left_join(peace_churn_f, by = c("ccode", "year")) |>
  mutate(churn_peace_ft1 = ifelse(!paste(ccode, year) %in% paste(peace_churn_f$ccode, peace_churn_f$year) & !is.na(churn),
                                  0, churn_peace_ft1),
         n_drop_ft1 = replace_na(n_drop_ft1, 0L),
         log_churn_peace_ft1 = log1p(churn_peace_ft1))
# Where no fatal opponent is dropped the peace measure must equal canonical churn.
n_bad_f <- sum(state_churn$n_drop_ft1 == 0 & abs(state_churn$churn_peace_ft1 - state_churn$churn) > 1e-9, na.rm = TRUE)
if (n_bad_f > 0) stop(n_bad_f, " state-years with no fatal opponent dropped have churn_peace_ft1 != churn")
cat(sprintf("  Fatal t-1 peace churn: %d state-years with fatal opponents dropped; %d undefined\n",
            sum(state_churn$n_drop_ft1 > 0), sum(is.na(state_churn$churn_peace_ft1) & !is.na(state_churn$churn))))

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
  )

# [2026-10-12] Onset lags from the full COW state-year universe in dyad_data (same construction as
# 01_data_cleaning/02_compute_network_measures.R). state_churn holds transition rows only, so lagging the
# previous ROW left the IV missing in every state's first transition year (217 state-years, mostly new
# states) although the state was in the system, and its conflict observed, in t-1. A lag is now NA only
# if the state was not in the system in t-k (in the current data every transition row has t-1 defined).
onset_univ <- bind_rows(dyad_data |> distinct(ccode = ccode1, year),
                        dyad_data |> distinct(ccode = ccode2, year)) |>
  distinct() |>
  left_join(
    dyad_data |>
      filter(cowmidongoing == 1, ccode1 < ccode2) |>
      select(ccode1, ccode2, year, cowmidfatal, cowmidonset) |>
      pivot_longer(c(ccode1, ccode2), values_to = "ccode") |>
      group_by(ccode, year) |>
      summarise(u_onset = sum(cowmidonset == 1, na.rm = TRUE),
                u_fatal_onset = sum(cowmidfatal == 1 & cowmidonset == 1, na.rm = TRUE),
                .groups = "drop"),
    by = c("ccode", "year")) |>
  mutate(u_onset = as.numeric(replace_na(u_onset, 0L)),
         u_fatal_onset = as.numeric(replace_na(u_fatal_onset, 0L)))
# the universe must reproduce the panel's own year-t counts
chk_univ <- inner_join(state_churn |> select(ccode, year, mid_onset_count, fatal_mid_onset_count),
                       onset_univ, by = c("ccode", "year"))
stopifnot(nrow(chk_univ) == nrow(state_churn), all(chk_univ$mid_onset_count == chk_univ$u_onset),
          all(chk_univ$fatal_mid_onset_count == chk_univ$u_fatal_onset))
rm(chk_univ)
state_churn <- state_churn |>
  left_join(onset_univ |> transmute(ccode, year = year + 1L,
                                    onset_lag1 = u_onset, fatal_onset_lag1 = u_fatal_onset,
                                    log_onset_lag1 = log1p(u_onset), log_fatal_onset_lag1 = log1p(u_fatal_onset)),
            by = c("ccode", "year")) |>
  left_join(onset_univ |> transmute(ccode, year = year + 2L,
                                    onset_lag2 = u_onset, fatal_onset_lag2 = u_fatal_onset,
                                    log_onset_lag2 = log1p(u_onset), log_fatal_onset_lag2 = log1p(u_fatal_onset)),
            by = c("ccode", "year")) |>
  arrange(ccode, year)
stopifnot(!anyDuplicated(state_churn[c("ccode", "year")]))
cat(sprintf("  t-1 onset defined for %d of %d state-years (NA = not in the system in t-1)\n",
            sum(!is.na(state_churn$log_fatal_onset_lag1)), nrow(state_churn)))

cat(sprintf("  MID count: mean=%.3f, sd=%.3f, max=%d\n",
            mean(state_churn$mid_count, na.rm = TRUE),
            sd(state_churn$mid_count, na.rm = TRUE),
            max(state_churn$mid_count, na.rm = TRUE)))
cat(sprintf("  MID onset: mean=%.3f, sd=%.3f, max=%d\n",
            mean(state_churn$mid_onset_count, na.rm = TRUE),
            sd(state_churn$mid_onset_count, na.rm = TRUE),
            max(state_churn$mid_onset_count, na.rm = TRUE)))

# --------------------------------------------------------------------------- #
# 0d. Two-stage multilayer community membership (DynMux, unweighted Jaccard)
# --------------------------------------------------------------------------- #
#
# Stage 1: Leiden CPM within each year on main_weight (unweighted Manhattan,
#   range01-normalized; all weights in [0, 1]). CPM with gamma > 0 splits
#   communities when internal density falls below gamma.
# Stage 2: communities in adjacent years are linked by UNWEIGHTED Jaccard
#   similarity of their node sets and tracked forward in time ("filtered":
#   the meta id of year t depends only on years <= t), with chain repair.
# leiden_comm below is the TRACKED META-COMMUNITY id, so one id means the
# same community across years. [2026-10-12] The previous code used the
# per-year labels (ml_fit$layer_communities[[i]]$membership), so community 3
# in 1850 and community 3 in 1950 were treated as one cluster. That was a bug.
# Unweighted (not weighted) Jaccard: in fit_multilayer_weighted_jaccard a
# zero-strength isolate has similarity 0/0 -> 0 with itself, so isolates never
# chain across years; unweighted Jaccard chains them.
# Layer construction is unchanged: it keeps weight-0 edges (!is.na, not > 0).
# --------------------------------------------------------------------------- #

cat("Computing multilayer communities (two-stage DynMux: Leiden CPM + Jaccard tracking)...\n")

# Build list of igraph objects using main_weight (positive edges only)
years_with_data <- sort(unique(state_churn$year))
layer_list <- list()
year_map <- c()

for (yr in years_with_data) {
  # [2026-10-12] dyad_data stores BOTH directions of every pair (i,j) and (j,i). Building an undirected
  # graph from all rows created two parallel edges per pair, doubling every tie's strength relative to
  # gamma. Collapse to one edge per unordered pair. 0.05% of pairs carry different weights in the two
  # directions; the stronger direction (max) is kept.
  edges_yr <- dyad_data |>
    filter(year == yr, !is.na(.data[[MAIN_WEIGHT_VAR]])) |>
    transmute(a = pmin(ccode1, ccode2), b = pmax(ccode1, ccode2), w = .data[[MAIN_WEIGHT_VAR]]) |>
    group_by(a, b) |>
    summarise(w = max(w), .groups = "drop")
  stopifnot(all(edges_yr$a != edges_yr$b), !anyDuplicated(edges_yr[c("a", "b")]))
  if (nrow(edges_yr) < 2) next
  g <- graph_from_data_frame(edges_yr |> select(a, b), directed = FALSE)
  E(g)$weight <- edges_yr$w
  stopifnot(!any(which_multiple(g)), !any(which_loop(g)), ecount(g) == nrow(edges_yr))
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
ml_fit <- fit_multilayer_jaccard(
  layer_list,
  algorithm = "leiden",
  layer_links = layer_links,           # adjacent years, weight 1 (unchanged)
  resolution_parameter = LEIDEN_GAMMA, # mean edge weight, unchanged (~0.1506)
  objective = "cpm",
  directed = FALSE,
  seed = 123,
  allow_unequal_nodes = TRUE,          # states enter and exit across years
  tracking = "filtered",
  repair_chains = TRUE,
  repair_min_similarity = 0.2
)
meta <- extract_meta_membership(ml_fit)   # tracked meta-communities, NOT layer_communities
stopifnot(length(meta) == length(layer_list))

# Extract the tracked meta-community of every node back to state-year format.
leiden_data <- vector("list", length(meta))
for (i in seq_along(meta)) {
  nm <- V(layer_list[[i]])$name
  stopifnot(length(nm) == length(meta[[i]]), !anyNA(meta[[i]]))
  leiden_data[[i]] <- data.frame(ccode = as.numeric(nm), year = year_map[i],
                                 leiden_comm = as.integer(meta[[i]]))
}
leiden_df <- bind_rows(leiden_data)
stopifnot(!anyDuplicated(leiden_df[, c("ccode", "year")]), !anyNA(leiden_df$ccode))

cat(sprintf("  Meta-communities: %d state-year assignments, %d distinct meta ids, median %d communities per year\n",
            nrow(leiden_df), n_distinct(leiden_df$leiden_comm),
            as.integer(median(tapply(leiden_df$leiden_comm, leiden_df$year, n_distinct)))))

state_churn <- state_churn |>
  left_join(leiden_df |> select(ccode, year, leiden_comm), by = c("ccode", "year"))
stopifnot(!anyNA(state_churn$leiden_comm))   # every state-year has a tracked community

# Community diagnostics on the analysis sample (one row per state-year).
# Isolate = state with zero tie strength at t (is_isolated_t). A spell is a
# run of consecutive isolate years for one state.
iso <- state_churn |> filter(is_isolated_t == 1) |> arrange(ccode, year) |>
  group_by(ccode) |> mutate(spell = cumsum(c(1, diff(year) != 1))) |> ungroup()
cat(sprintf("  Diagnostics: %d non-isolate meta-communities | %d isolate state-years (%.1f%%) | %d isolate spells | %d isolate meta ids\n",
            n_distinct(state_churn$leiden_comm[state_churn$is_isolated_t == 0]),
            nrow(iso), 100 * nrow(iso) / nrow(state_churn),
            n_distinct(paste(iso$ccode, iso$spell)), n_distinct(iso$leiden_comm)))

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
#   main: TWFE raw churn on log MID ONSET   (no/full ctrl x mid/fatal)
#   A1:   TWFE raw churn on log MID COUNT   (no/full ctrl x mid/fatal)
#   A2:   no-FE  raw churn on log {onset, count} + poly(year, 2) (full ctrl)
#   A3:   3-way FE (ccode + year + leiden_comm) raw churn on log
#         {onset, count} (full ctrl)
#   A4:   TWFE raw churn on raw MID ONSET      (no/full ctrl x mid/fatal)
# (A5 -- TWFE raw churn on raw MID COUNT -- dropped; see note below.)
#
# All TWFE models cluster two-way (ccode + year). The three-way FE model
# clusters three-way (ccode + year + leiden_comm). The no-FE model
# clusters one-way (ccode) since there are no within-state deviations to
# absorb.
# =========================================================================== #

# [2026-10-12] Partner war exposure is no longer a control (H3 dropped). The
# MID and fatal-MID control sets are now identical; ctrl_for_iv() is kept so
# the spec grid below is unchanged.
#
# is_isolated_t1 is a predetermined indicator for "state had no partners
# at t-1". It absorbs the level shift from `total_gain` (degree_t1 = 0,
# churn = 1) and `both_isolated` (degree_t1 = 0, churn = 0) cases that
# otherwise flatten the within-state slope. Predetermined ⇒ exogenous to
# year-t MID involvement.
CTRL_MID   <- "democracy + autocracy + cinc + major_power + is_isolated_t1"
CTRL_FATAL <- "democracy + autocracy + cinc + major_power + is_isolated_t1"
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
       dv = "churn", iv = "log_mid_onset_count",       controls = NULL),
  list(key = "main__full__mid",
       dv = "churn", iv = "log_mid_onset_count",       controls = CTRL_MID),
  list(key = "main__none__fatal",
       dv = "churn", iv = "log_fatal_mid_onset_count", controls = NULL),
  list(key = "main__full__fatal",
       dv = "churn", iv = "log_fatal_mid_onset_count", controls = CTRL_FATAL)
)

specs_a1 <- list(
  list(key = "a1__none__mid",
       dv = "churn", iv = "log_mid_count",        controls = NULL),
  list(key = "a1__full__mid",
       dv = "churn", iv = "log_mid_count",        controls = CTRL_MID),
  list(key = "a1__none__fatal",
       dv = "churn", iv = "log_fatal_mid_count",  controls = NULL),
  list(key = "a1__full__fatal",
       dv = "churn", iv = "log_fatal_mid_count",  controls = CTRL_FATAL)
)

# A2: no fixed effects + cubic year trend
specs_a2 <- list(
  list(key = "a2__nofe__onset_mid",
       dv = "churn", iv = "log_mid_onset_count",
       controls = paste(CTRL_MID, "+ poly(year, 2)"),
       fe = "", vcov_form = "ccode"),
  list(key = "a2__nofe__onset_fatal",
       dv = "churn", iv = "log_fatal_mid_onset_count",
       controls = paste(CTRL_FATAL, "+ poly(year, 2)"),
       fe = "", vcov_form = "ccode"),
  list(key = "a2__nofe__count_mid",
       dv = "churn", iv = "log_mid_count",
       controls = paste(CTRL_MID, "+ poly(year, 2)"),
       fe = "", vcov_form = "ccode"),
  list(key = "a2__nofe__count_fatal",
       dv = "churn", iv = "log_fatal_mid_count",
       controls = paste(CTRL_FATAL, "+ poly(year, 2)"),
       fe = "", vcov_form = "ccode")
)

# A3: three-way FE (state + year + Leiden community)
specs_a3 <- list(
  list(key = "a3__3wfe__onset_mid",
       dv = "churn", iv = "log_mid_onset_count",  controls = CTRL_MID,
       fe = "ccode + year + leiden_comm",
       vcov_form = "ccode + year + leiden_comm"),
  list(key = "a3__3wfe__onset_fatal",
       dv = "churn", iv = "log_fatal_mid_onset_count", controls = CTRL_FATAL,
       fe = "ccode + year + leiden_comm",
       vcov_form = "ccode + year + leiden_comm"),
  list(key = "a3__3wfe__count_mid",
       dv = "churn", iv = "log_mid_count",        controls = CTRL_MID,
       fe = "ccode + year + leiden_comm",
       vcov_form = "ccode + year + leiden_comm"),
  list(key = "a3__3wfe__count_fatal",
       dv = "churn", iv = "log_fatal_mid_count",  controls = CTRL_FATAL,
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
# main-text "lag robustness" appendix). These use raw churn on the
# LHS to stay consistent with the headline TWFE.
# --------------------------------------------------------------------------- #

cat("\n========== LAG MODELS (lagged onset IVs, raw churn DV) ==========\n")

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
      as.formula(paste("churn ~", lv, "| ccode + year")),
      data = state_churn, vcov = ~ccode + year)
    lag_models[[paste0(wv, "_", tag, "_full")]] <- feols(
      as.formula(paste("churn ~", lv, "+", ctrl_wv, "| ccode + year")),
      data = state_churn, vcov = ~ccode + year)
  }
  cat(sprintf("  %s lags: done\n", wv))
}

# --------------------------------------------------------------------------- #
# Panel Local Projections (Jordà 2005) — matches the new main spec.
# --------------------------------------------------------------------------- #
# 2 controls x 2 IVs x 6 horizons = 24 LP coefficients, logged IV, raw churn DV.
#
# For each horizon h:
#   churn_{i, t+h} = alpha_i + gamma_t + beta_h * log1p(IV)_{i,t}
#                            + (controls)_{i,t} + epsilon_{i,t+h}
# --------------------------------------------------------------------------- #

cat("\n========== PANEL LOCAL PROJECTIONS ==========\n")

MAX_HORIZON <- 5

# [2026-10-12] gap-safe leads: churn at t+h only if that row is exactly h calendar years later.
lp_data <- state_churn |>
  group_by(ccode) |>
  arrange(year, .by_group = TRUE) |>
  mutate(across(churn, list(
    h0 = ~.x,
    h1 = ~if_else(lead(year, 1) == year + 1, lead(.x, 1), NA_real_),
    h2 = ~if_else(lead(year, 2) == year + 2, lead(.x, 2), NA_real_),
    h3 = ~if_else(lead(year, 3) == year + 3, lead(.x, 3), NA_real_),
    h4 = ~if_else(lead(year, 4) == year + 4, lead(.x, 4), NA_real_),
    h5 = ~if_else(lead(year, 5) == year + 5, lead(.x, 5), NA_real_)
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

# [2026-10-12] Main text: fatal MID onset (headline). All-MID companion in the appendix.
# LP horizon h is churn at t+h on conflict at t; h = 1 corresponds to the headline t-1 timing.
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
    geom_point(aes(shape = sig, color = sig), size = 3.5) +   # [2026-10-12] colour keyed to sig, matching the scale
    facet_wrap(~controls_label) + 
    scale_color_manual(values = c("p-value < 0.05"   = "#08519c",
                                    "Null" = "#B03366"),
                       name = NULL) +
    scale_shape_manual(values = c("p-value < 0.05" = 17, "Null" = 16), name = NULL) +   # one merged legend
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

# --- Main text figure [2026-10-12]: fatal MID onset (the headline conflict measure) ---
p_lp_fatal <- build_lp_plot(
  lp_plot_data |> filter(iv_family == "fatal") |> mutate(sig = ifelse(sig == T, "p-value < 0.05", "Null")),
  title_iv = "log(1 + Fatal MID Onset)"
)
ggsave(file.path(LP_FIG_DIR_MAIN, "figure1_panel_lp.png"),
       p_lp_fatal, width = 7, height = 5, dpi = 300, bg = "white")
cat("  Panel LP figure (fatal MID, main) -> figure1_panel_lp.png\n")

# --- Appendix figure: all-MID companion ---
p_lp_mid <- build_lp_plot(
  lp_plot_data |> filter(iv_family == "mid") |> mutate(sig = ifelse(sig == T, "p-value < 0.05", "Null")),
  title_iv = "log(1 + MID Onset)",
  caption_extra = " Companion to main-text Figure 1 (fatal MID onset)."
)
ggsave(file.path(LP_FIG_DIR_APPX, "figureA_panel_lp_mid.png"),
       p_lp_mid, width = 7, height = 5, dpi = 300, bg = "white")
cat("  Panel LP figure (all MIDs, appendix) -> figureA_panel_lp_mid.png\n")

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

# [2026-10-12] The m1/m2/m5 building blocks (IV alone, PWE alone, TWFE with
# community-clustered SEs) were removed: the main table is built by t3_ladder()
# below and no longer has a PWE column.
state_churn_leiden <- state_churn |> filter(!is.na(leiden_comm))   # == state_churn (asserted above)

# --- Helper: rescale extracted coefficients and SEs by 100 for display ---
# Churn is a 0-1 index, so raw coefficients print as 0.00x. Multiplying the
# extracted coefficients and SEs by 100 (DV in percentage points) lets the
# appendix tables use two decimals. Fit statistics and p-values are unchanged.
scale_tr <- function(ms, k = 100) lapply(ms, function(m) {
  tr <- texreg::extract(m)
  tr@coef <- tr@coef * k
  tr@se   <- tr@se * k
  tr
})

# --- Helper: write a 5-column table from a list of fitted models ---
# (Each entry can be a feols/lm object directly; texreg dispatches by
#  class. Use this for the new main + fatal-appendix tables.)
write_state_5col <- function(models,
                              column_headers,
                              caption,
                              label,
                              note_extra,
                              file_path,
                              digits = 2) {
  base_note <- paste0(
    "\\parbox{\\linewidth}{\\vspace{2pt}\\textit{Note:} ",
    note_extra, " ",
    "Coefficients and standard errors are multiplied by 100. ",
    "Standard errors in parentheses. \\textbf{Columns (1)--(3) include ",
    "a quadratic year trend ($\\texttt{poly(year, 2)}$); terms suppressed ",
    "for legibility. In columns (4)--(5), year effects are absorbed by ",
    "the year fixed effect.} Columns (1) and (2) are baseline OLS ",
    "regressions of the focal IV (and the year trend) with state-clustered ",
    "SEs. Column (3) adds the full controls (",
    "democracy, autocracy, CINC, major power, ",
    "$\\mathbb{1}\\{\\text{isolated at }t-1\\}$); SEs clustered by state. ",
    "Column (4) is the headline two-way ",
    "fixed-effects spec (state $+$ year), two-way clustered SEs. ",
    "Column (5) keeps the same TWFE estimator but clusters SEs three ",
    "ways (state $+$ year $+$ two-stage multilayer community); sample ",
    "restricted to state-years with non-missing community assignment. ",
    "COW system, 1816--2014. ",
    "$^{***}p<0.001$; $^{**}p<0.01$; $^{*}p<0.05$; $^{\\dagger}p<0.1$.}")
  texreg(scale_tr(models),
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

# --- Table 3 (main text): four stacked panels, 5 cols each [2026-10-12] --
# DV is RAW churn (not logged). Churn is a 0-1 share, so log(1 + churn) is
# nearly linear in it; the logged version is the appendix companion below.
#   Panel A: log(1 + MID onset),       churn over all dyads
#   Panel B: log(1 + MID onset),       churn_peace (partners in an MID with i in t-1 or t dropped; section 0b2)
#   Panel C: log(1 + fatal MID onset), churn over all dyads
#   Panel D: log(1 + fatal MID onset), churn_peace
# Columns:
#   (1) IV + poly(year, 2), SEs clustered by state
#   (2) + controls + poly(year, 2), no FE, SEs clustered by state
#   (3) TWFE (state + year), SEs clustered by state and year
#   (4) TWFE (state + year), SEs clustered by state, year and meta-community
#   (5) three-way FE (state + year + meta-community), SEs clustered the same three ways
# Every state-year has a tracked meta-community (asserted in 0d), so column 4 uses the
# full sample; column 5 drops the fixed-effect singletons (its N is reported).
# Wild cluster bootstrap p-values for column 3 are produced by 23_wild_cluster_bootstrap.R.
state_churn_peace <- state_churn |> filter(!is.na(churn_peace))
CTRL_MID_PEACE <- "democracy + autocracy + cinc + major_power + is_isolated_t1"
t3_ladder <- function(dv, iv, ctrl, data) {
  list(
    feols(as.formula(paste(dv, "~", iv, "+ poly(year, 2)")), data = data, vcov = ~ccode),
    feols(as.formula(paste(dv, "~", iv, "+", ctrl, "+ poly(year, 2)")), data = data, vcov = ~ccode),
    feols(as.formula(paste(dv, "~", iv, "+", ctrl, "| ccode + year")), data = data, vcov = ~ccode + year),
    feols(as.formula(paste(dv, "~", iv, "+", ctrl, "| ccode + year")), data = data, vcov = ~ccode + year + leiden_comm),
    feols(as.formula(paste(dv, "~", iv, "+", ctrl, "| ccode + year + leiden_comm")), data = data,
          vcov = ~ccode + year + leiden_comm))
}
panel_A <- t3_ladder("churn",       "log_mid_onset_count",       CTRL_MID,       state_churn)
panel_B <- t3_ladder("churn_peace", "log_mid_onset_count",       CTRL_MID_PEACE, state_churn_peace)
panel_C <- t3_ladder("churn",       "log_fatal_mid_onset_count", CTRL_FATAL,     state_churn)
panel_D <- t3_ladder("churn_peace", "log_fatal_mid_onset_count", CTRL_MID_PEACE, state_churn_peace)
t3_panels <- list(A = panel_A, B = panel_B, C = panel_C, D = panel_D)
# column 2 must be a true no-FE model; columns 3-4 share the point estimates; column 4 uses the full sample
for (pn in names(t3_panels)) {
  pm <- t3_panels[[pn]]
  stopifnot(length(pm[[1]]$fixef_vars) == 0, length(pm[[2]]$fixef_vars) == 0,
            isTRUE(all.equal(coef(pm[[3]]), coef(pm[[4]]))), nobs(pm[[3]]) == nobs(pm[[4]]),
            identical(pm[[5]]$fixef_vars, c("ccode", "year", "leiden_comm")))
  cat(sprintf("  Panel %s: N = %s (cols 1-2), %s (cols 3-4), %s (col 5; %d FE singletons dropped)\n", pn,
              format(nobs(pm[[1]]), big.mark = ","), format(nobs(pm[[3]]), big.mark = ","),
              format(nobs(pm[[5]]), big.mark = ","), nobs(pm[[3]]) - nobs(pm[[5]])))
}

# Hand-built LaTeX (texreg cannot stack panels). Rows, in order, with the term name per IV.
t3_rows <- function(iv) list(
  c(paste0(if (grepl("fatal", iv)) "Log(1 + fatal MID onset)" else "Log(1 + MID onset)",
           if (grepl("lag1", iv)) " ($t{-}1$)" else ""), iv),
  c("Democracy",           "democracy"),
  c("Autocracy",           "autocracy"),
  c("CINC",                "cinc"),
  c("Major power",         "major_power"),
  c("Isolated at $t{-}1$", "is_isolated_t1"),
  c("Intercept",           "(Intercept)"))
# Coefficients and standard errors multiplied by 100 (percentage points of churn), two decimals,
# so no cell rounds to 0.00. Significance, R2 and N are unaffected by the rescaling.
T3_SCALE <- 100
t3_star <- function(p) ifelse(p < 0.001, "^{***}", ifelse(p < 0.01, "^{**}", ifelse(p < 0.05, "^{*}", ifelse(p < 0.1, "^{\\dagger}", ""))))
t3_panel <- function(models, iv) {
  out <- character(0)
  for (r in t3_rows(iv)) {
    est <- se <- character(length(models))
    for (m in seq_along(models)) {
      ct <- coeftable(models[[m]])
      if (r[2] %in% rownames(ct)) {
        est[m] <- sprintf("$%.2f%s$", T3_SCALE * ct[r[2], 1], t3_star(ct[r[2], 4]))
        se[m]  <- sprintf("$(%.2f)$", T3_SCALE * ct[r[2], 2])
      } else { est[m] <- "--"; se[m] <- "" }
    }
    out <- c(out, paste(r[1], "&", paste(est, collapse = " & "), "\\\\"), paste(" &", paste(se, collapse = " & "), "\\\\"))
  }
  g <- function(m, v) { n <- m$fixef_sizes; if (is.null(n) || !v %in% names(n)) "--" else as.character(n[[v]]) }
  w <- function(m, type) { x <- tryCatch(fixest::r2(m, type), error = function(e) NA); if (is.na(x)) "--" else if (x < 0.005) "$<$0.01" else sprintf("%.2f", x) }
  c(out, "\\midrule",
    paste("Num. obs. &", paste(sapply(models, function(m) format(nobs(m), big.mark = ",")), collapse = " & "), "\\\\"),
    paste("R$^2$ &", paste(sapply(models, w, "r2"), collapse = " & "), "\\\\"),
    paste("Within R$^2$ &", paste(sapply(models, w, "wr2"), collapse = " & "), "\\\\"),
    paste("Num. groups: States &", paste(sapply(models, g, "ccode"), collapse = " & "), "\\\\"),
    paste("Num. groups: Years &", paste(sapply(models, g, "year"), collapse = " & "), "\\\\"),
    paste("Num. groups: Communities &", paste(sapply(models, g, "leiden_comm"), collapse = " & "), "\\\\"))
}
t3_panel_titles <- c(
  A = "Panel A. $\\log(1+\\text{MID onset})$; churn over all dyads",
  B = "Panel B. $\\log(1+\\text{MID onset})$; churn over peace dyads only",
  C = "Panel C. $\\log(1+\\text{fatal MID onset})$; churn over all dyads",
  D = "Panel D. $\\log(1+\\text{fatal MID onset})$; churn over peace dyads only")
t3_ivs <- c(A = "log_mid_onset_count", B = "log_mid_onset_count", C = "log_fatal_mid_onset_count", D = "log_fatal_mid_onset_count")
t3_body <- function(panels) {
  out <- character(0)
  for (pn in names(panels)) {
    out <- c(out, if (pn != "A") "\\midrule",
             sprintf("\\multicolumn{6}{l}{\\textit{%s}} \\\\", t3_panel_titles[[pn]]), "\\midrule",
             t3_panel(panels[[pn]], t3_ivs[[pn]]))
  }
  out
}
t3_note <- paste0(
  "\\multicolumn{6}{l}{\\scriptsize{\\parbox{\\linewidth}{\\vspace{2pt}\\textit{Note:} ",
  "DV is churn, $1 - $ weighted Jaccard similarity of a state's ties in $t-1$ and $t$ (not logged). ",
  "Coefficients and standard errors are multiplied by 100. ",
  "The conflict measure is $\\log(1+\\text{MID onset})$ in Panels A--B and $\\log(1+\\text{fatal MID onset})$ in Panels C--D. ",
  "Panels B and D drop from churn every partner that was in an MID with the state in $t-1$ or $t$; ",
  "state-years whose partners were all MID opponents are undefined and excluded. ",
  "Columns (1)--(2) include a quadratic year trend and cluster SEs by state; column (2) adds the controls and has no fixed effects. Columns (2) onward use the state-years with non-missing controls (democracy and autocracy come from Polity), hence the smaller N. ",
  "Column (3) is the two-way fixed-effects spec (state $+$ year) with SEs clustered by state and year. ",
  "Column (4) keeps the TWFE estimator and clusters SEs by state, year and two-stage multilayer community ",
  "(Leiden CPM within year, Jaccard-tracked across years). ",
  "Column (5) adds the community as a third fixed effect with the same three-way clustering; ",
  "state-years that are singletons within a fixed effect are dropped. COW system, 1816--2014. ",
  "$^{***}p<0.001$; $^{**}p<0.01$; $^{*}p<0.05$; $^{\\dagger}p<0.1$.}}} \\\\")
t3_tex <- c(
  "\\begin{table}[p]", "\\scriptsize", "\\begin{center}", "\\begin{tabular}{l c c c c c}", "\\toprule",
  " & (1) IV + trend & (2) + Ctrls, no FE & (3) TWFE & (4) TWFE + Comm cluster & (5) 3-way FE \\\\", "\\midrule",
  t3_body(t3_panels), "\\bottomrule", t3_note,
  "\\end{tabular}", "\\end{center}",
  "\\caption{State-level regressions for conflict and churn. MID onset (Panels A--B) and fatal MID onset (Panels C--D), each over all dyads and over peace dyads only.}",
  "\\label{tab:state_main}", "\\end{table}")
stopifnot(!any(grepl("artner war|PWE", t3_tex)))

# [2026-10-12] MAIN TEXT = columns 1-3 only (4 panels x 3 columns = 12 specifications):
#   (1) IV + quadratic year trend, no controls; (2) + controls, still no FE; (3) TWFE (state + year),
#   SEs clustered by state and year. Columns 4-5 (community-clustered SEs; community FE) go to an
#   appendix table that repeats column 3 for comparison. t3_tex (all five columns) is still used below
#   for the log(1 + churn) appendix companion.
# [2026-10-12] Split a four-panel table into two floats (Panels A-B; Panels C-D) so neither exceeds a
# page. Part 2 repeats the column header and note; its caption is prefixed "(Continued)" and its label
# gets the suffix "_cd". \ContinuedFloat (caption package) keeps one table number for both parts.
# Fails loudly if the table does not have the expected layout.
split_ab_cd <- function(tex) {
  i_top <- which(tex == "\\toprule"); i_bot <- which(tex == "\\bottomrule")
  i_C   <- grep("textit\\{Panel C\\.", tex)
  i_cap <- grep("^\\\\caption\\{", tex); i_lab <- grep("^\\\\label\\{", tex)
  stopifnot(length(i_top) == 1, length(i_bot) == 1, length(i_C) == 1, length(i_cap) == 1, length(i_lab) == 1,
            tex[i_top + 2] == "\\midrule", tex[i_C - 1] == "\\midrule", i_top + 2 < i_C - 2, i_C < i_bot)
  pre   <- tex[1:(i_top + 2)]
  post  <- tex[i_bot:length(tex)]
  post2 <- post
  post2[post2 == tex[i_cap]] <- paste0("\\ContinuedFloat ", sub("^\\\\caption\\{", "\\\\caption{(Continued) ", tex[i_cap]))
  post2[post2 == tex[i_lab]] <- sub("\\}$", "_cd}", tex[i_lab])
  part1 <- c(pre, tex[(i_top + 3):(i_C - 2)], post)
  part2 <- c(pre, tex[i_C:(i_bot - 1)], post2)
  stopifnot(!any(grepl("Panel C\\.|Panel D\\.", part1)), !any(grepl("Panel A\\.|Panel B\\.", part2)),
            sum(grepl("Panel [AB]\\.", part1)) == 2, sum(grepl("Panel [CD]\\.", part2)) == 2,
            any(grepl("_cd\\}$", part2)), any(grepl("\\(Continued\\)", part2)), any(grepl("ContinuedFloat", part2)))
  c(part1, "", part2)
}
t3_body_k <- function(panels, cols, k) {
  out <- character(0)
  for (pn in names(panels)) {
    pan <- t3_panel(panels[[pn]][cols], t3_ivs[[pn]])
    if (!5 %in% cols) pan <- pan[!grepl("^Num. groups: Communities", pan)]
    if (!1 %in% cols && !2 %in% cols) {   # FE-only tables: the Intercept row is all "--"; drop it and its SE row
      ii <- which(grepl("^Intercept &", pan)); stopifnot(length(ii) == 1, grepl("^ &( +&)+ *\\\\\\\\$", pan[ii + 1]))
      pan <- pan[-c(ii, ii + 1)]
    }
    out <- c(out, if (pn != "A") "\\midrule",
             sprintf("\\multicolumn{%d}{l}{\\textit{%s}} \\\\", k, t3_panel_titles[[pn]]), "\\midrule", pan)
  }
  out
}
t3_note_common <- paste0(
  "DV is churn, $1 - $ weighted Jaccard similarity of a state's ties in $t-1$ and $t$ (not logged). ",
  "Coefficients and standard errors are multiplied by 100. ",
  "The conflict measure, measured in year $t$, is $\\log(1+\\text{MID onset})$ in Panels A--B and $\\log(1+\\text{fatal MID onset})$ in Panels C--D. ",
  "Panels B and D drop from churn every partner that was in an MID with the state in $t-1$ or $t$; ",
  "state-years whose partners were all MID opponents are undefined and excluded. ")
t3_stars <- "$^{***}p<0.001$; $^{**}p<0.01$; $^{*}p<0.05$; $^{\\dagger}p<0.1$.}}} \\\\"
t3_note_main <- paste0(
  "\\multicolumn{4}{l}{\\scriptsize{\\parbox{\\linewidth}{\\vspace{2pt}\\textit{Note:} ", t3_note_common,
  "Columns (1)--(2) include a quadratic year trend and cluster SEs by state; column (2) adds the controls and has no fixed effects. Columns (2) onward use the state-years with non-missing controls (democracy and autocracy come from Polity), hence the smaller N. ",
  "Column (3) is the two-way fixed-effects spec (state $+$ year) with SEs clustered by state and year. ",
  "Community-clustered SEs and community fixed effects are in Appendix Table~\\ref{tab:state_main_comm}. ",
  "COW system, 1816--2014. ", t3_stars)
t3_tex_main <- c(
  "\\begin{table}[p]", "\\scriptsize", "\\begin{center}", "\\begin{tabular}{l c c c}", "\\toprule",
  " & (1) IV + trend & (2) + Ctrls, no FE & (3) TWFE \\\\", "\\midrule",
  t3_body_k(t3_panels, 1:3, 4), "\\bottomrule", t3_note_main,
  "\\end{tabular}", "\\end{center}",
  "\\caption{State-level regressions for conflict and churn. MID onset (Panels A--B) and fatal MID onset (Panels C--D), each over all dyads and over peace dyads only.}",
  "\\label{tab:state_main}", "\\end{table}")
# Fail loudly if the main table is not exactly 4 panels x 3 value columns.
t3_iv_rows <- t3_tex_main[grepl("^Log\\(1 \\+ (fatal )?MID onset", t3_tex_main)]
stopifnot(length(t3_iv_rows) == 4,
          all(lengths(regmatches(t3_iv_rows, gregexpr("&", t3_iv_rows))) == 3),
          !any(grepl("(4)", t3_tex_main, fixed = TRUE)), !any(grepl("(5)", t3_tex_main, fixed = TRUE)),
          sum(grepl("^Num. groups: Communities", t3_tex_main)) == 0,
          !any(grepl("artner war|PWE", t3_tex_main)))
# [2026-10-12] superseded by the fatal-onset t-1 headline below; t3_tex_main is built for validation only.

t3_note_app <- paste0(
  "\\multicolumn{4}{l}{\\scriptsize{\\parbox{\\linewidth}{\\vspace{2pt}\\textit{Note:} ", t3_note_common,
  "All columns include the controls and state and year fixed effects; the conflict measure is in year $t$. ",
  "Column (3) repeats the main-text TWFE specification (SEs clustered by state and year). ",
  "Column (4) keeps the TWFE estimator and clusters SEs by state, year and two-stage multilayer community ",
  "(Leiden CPM within year, Jaccard-tracked across years). ",
  "Column (5) adds the community as a third fixed effect with the same three-way clustering; ",
  "state-years that are singletons within a fixed effect are dropped. COW system, 1816--2014. ", t3_stars)
t3_tex_app <- c(
  "\\begin{table}[p]", "\\scriptsize", "\\begin{center}", "\\begin{tabular}{l c c c}", "\\toprule",
  " & (3) TWFE & (4) TWFE + Comm cluster & (5) 3-way FE \\\\", "\\midrule",
  t3_body_k(t3_panels, 3:5, 4), "\\bottomrule", t3_note_app,
  "\\end{tabular}", "\\end{center}",
  "\\caption{Appendix. State-level regressions with community-clustered standard errors and community fixed effects, companion to main-text Table~\\ref{tab:state_main}.}",
  "\\label{tab:state_main_comm}", "\\end{table}")
t3_iv_rows_app <- t3_tex_app[grepl("^Log\\(1 \\+ (fatal )?MID onset", t3_tex_app)]
stopifnot(length(t3_iv_rows_app) == 4,
          all(lengths(regmatches(t3_iv_rows_app, gregexpr("&", t3_iv_rows_app))) == 3),
          sum(grepl("^Num. groups: Communities", t3_tex_app)) == 4,
          !any(grepl("artner war|PWE", t3_tex_app)))
# [2026-10-12] superseded: the community appendix table now accompanies the fatal-onset t-1 headline (below).

# Appendix companion: the same four-panel table with log(1 + churn) as the DV.
# The main text and every other appendix model use raw churn; this is the one logged-DV table.
t3_panels_log <- list(
  A = t3_ladder("log_churn",       "log_mid_onset_count",       CTRL_MID,       state_churn),
  B = t3_ladder("log_churn_peace", "log_mid_onset_count",       CTRL_MID_PEACE, state_churn_peace),
  C = t3_ladder("log_churn",       "log_fatal_mid_onset_count", CTRL_FATAL,     state_churn),
  D = t3_ladder("log_churn_peace", "log_fatal_mid_onset_count", CTRL_MID_PEACE, state_churn_peace))
for (pn in names(t3_panels_log)) stopifnot(length(t3_panels_log[[pn]][[2]]$fixef_vars) == 0,
                                           all(sapply(t3_panels_log[[pn]], nobs) == sapply(t3_panels[[pn]], nobs)))
t3_tex_log <- t3_tex
iB <- which(t3_tex == t3_body(t3_panels)[1]); stopifnot(length(iB) == 1)
nB <- length(t3_body(t3_panels)); stopifnot(identical(t3_tex[iB:(iB + nB - 1)], t3_body(t3_panels)))
t3_tex_log[iB:(iB + nB - 1)] <- t3_body(t3_panels_log)
n_sub <- sum(grepl("DV is churn, ", t3_tex_log, fixed = TRUE)) + sum(grepl("label{tab:state_main}", t3_tex_log, fixed = TRUE)) +
  sum(grepl("caption{State-level regressions for conflict and churn.", t3_tex_log, fixed = TRUE))
stopifnot(n_sub == 3)
t3_tex_log <- sub("DV is churn, $1 - $ weighted Jaccard similarity of a state's ties in $t-1$ and $t$ (not logged). ",
                  "DV is $\\log(1+\\text{churn})$, where churn is $1 - $ weighted Jaccard similarity of a state's ties in $t-1$ and $t$. ",
                  t3_tex_log, fixed = TRUE)
t3_tex_log <- sub("\\label{tab:state_main}", "\\label{tab:state_main_logchurn}", t3_tex_log, fixed = TRUE)
t3_tex_log <- sub("\\caption{State-level regressions for conflict and churn.",
                  "\\caption{Appendix. State-level regressions with $\\log(1+\\text{churn})$ as the dependent variable, companion to main-text Table~\\ref{tab:state_main}.",
                  t3_tex_log, fixed = TRUE)
stopifnot(any(grepl("state_main_logchurn", t3_tex_log, fixed = TRUE)), any(grepl("DV is $\\log(1+", t3_tex_log, fixed = TRUE)))
# [2026-10-12] superseded: the log(1+churn) companion now accompanies the fatal-onset t-1 headline (below).
write_csv(bind_rows(lapply(names(t3_panels), function(pn) bind_rows(lapply(seq_along(t3_panels[[pn]]), function(m)
  broom::tidy(t3_panels[[pn]][[m]]) |> mutate(panel = pn, column = m, iv = t3_ivs[[pn]],
                                             dv = if (pn %in% c("A", "C")) "churn" else "churn_peace",
                                             nobs = nobs(t3_panels[[pn]][[m]])))))),
  file.path(OUTPUT_DIR, "tableA_state_allmid_t_coefs.csv"))

# Appendix: the four-panel ladder with the IV at t (MID and fatal MID onset; any-MID peace rule in t-1 or t).
t3_tex_allmid <- t3_tex
stopifnot(sum(t3_tex_allmid == "\\label{tab:state_main}") == 1,
          sum(grepl("^\\\\caption\\{State-level regressions for conflict and churn\\.", t3_tex_allmid)) == 1)
t3_tex_allmid[t3_tex_allmid == "\\label{tab:state_main}"] <- "\\label{tab:state_allmid_t}"
t3_tex_allmid <- sub("\\caption{State-level regressions for conflict and churn.",
                     "\\caption{Appendix. State-level regressions with the conflict measure in year $t$ (MID onset and fatal MID onset), companion to main-text Table~\\ref{tab:state_main}.",
                     t3_tex_allmid, fixed = TRUE)
writeLines(split_ab_cd(t3_tex_allmid), file.path(APPX_DIR, "tableA_state_allmid_t.tex"))
cat("  tab:state_allmid_t -> results/appendix/tableA_state_allmid_t.tex (Panels A-D, IV at t, columns 1-5)\n")

# =========================================================================== #
# [2026-10-12] HEADLINE, main-text Table 3: log(1 + fatal MID onset) at t-1.
#   Panel A: churn over all dyads. Panel B: churn over peace dyads (partners that were opponents in a
#   fatal MID in t-1 dropped; churn_peace_ft1). Main text shows ladder columns 1, 2, 3 and 5 (5 = three-way FE,
#   printed as column 4) [2026-10-13]; the appendix shows columns 3-5 and the same table with log(1 + churn).
# =========================================================================== #
HL_IV <- "log_fatal_onset_lag1"
state_churn_peace_f <- state_churn |> filter(!is.na(churn_peace_ft1))
hl_panels <- list(
  A = t3_ladder("churn",           HL_IV, CTRL_FATAL, state_churn),
  B = t3_ladder("churn_peace_ft1", HL_IV, CTRL_FATAL, state_churn_peace_f))
hl_panels_log <- list(
  A = t3_ladder("log_churn",           HL_IV, CTRL_FATAL, state_churn),
  B = t3_ladder("log_churn_peace_ft1", HL_IV, CTRL_FATAL, state_churn_peace_f))
for (pn in names(hl_panels)) {
  pm <- hl_panels[[pn]]
  stopifnot(length(pm[[1]]$fixef_vars) == 0, length(pm[[2]]$fixef_vars) == 0,
            isTRUE(all.equal(coef(pm[[3]]), coef(pm[[4]]))), nobs(pm[[3]]) == nobs(pm[[4]]),
            identical(pm[[5]]$fixef_vars, c("ccode", "year", "leiden_comm")),
            all(sapply(hl_panels_log[[pn]], nobs) == sapply(pm, nobs)))
  cat(sprintf("  Headline panel %s: N = %s (cols 1-2), %s (cols 3-4), %s (col 5)\n", pn,
              format(nobs(pm[[1]]), big.mark = ","), format(nobs(pm[[3]]), big.mark = ","), format(nobs(pm[[5]]), big.mark = ",")))
}
hl_titles <- c(
  A = "Panel A. $\\log(1+\\text{fatal MID onset})_{t-1}$; churn over all dyads",
  B = "Panel B. $\\log(1+\\text{fatal MID onset})_{t-1}$; churn over peace dyads (fatal-MID opponents in $t-1$ dropped)")
hl_body <- function(panels, cols, k) {
  out <- character(0)
  for (pn in names(panels)) {
    pan <- t3_panel(panels[[pn]][cols], HL_IV)
    if (!5 %in% cols) pan <- pan[!grepl("^Num. groups: Communities", pan)]
    if (!1 %in% cols && !2 %in% cols) {
      ii <- which(grepl("^Intercept &", pan)); stopifnot(length(ii) == 1, grepl("^ &( +&)+ *\\\\\\\\$", pan[ii + 1]))
      pan <- pan[-c(ii, ii + 1)]
    }
    out <- c(out, if (pn != "A") "\\midrule",
             sprintf("\\multicolumn{%d}{l}{\\textit{%s}} \\\\", k, hl_titles[[pn]]), "\\midrule", pan)
  }
  out
}
hl_note_common <- paste0(
  "Coefficients and standard errors are multiplied by 100. ",
  "The conflict measure is $\\log(1+\\text{fatal MID onset})$ in year $t-1$. ",
  "Panel B drops from churn every partner that was the state's opponent in a fatal MID ongoing in $t-1$; ",
  "state-years whose partners were all such opponents are undefined and excluded. ")
hl_stars <- "$^{***}p<0.001$; $^{**}p<0.01$; $^{*}p<0.05$; $^{\\dagger}p<0.1$.}}} \\\\"
hl_note <- function(k, dv_text, cols_text) paste0(
  "\\multicolumn{", k, "}{l}{\\scriptsize{\\parbox{\\linewidth}{\\vspace{2pt}\\textit{Note:} ", dv_text, hl_note_common, cols_text,
  "COW system, 1816--2014. ", hl_stars)
HL_DV_RAW <- "DV is churn, $1 - $ weighted Jaccard similarity of a state's ties in $t-1$ and $t$ (not logged). "
HL_COLS_14 <- paste0(
  "Columns (1)--(2) include a quadratic year trend and cluster SEs by state; column (2) adds the controls and has no fixed effects. Columns (2) onward use the state-years with non-missing controls (democracy and autocracy come from Polity), hence the smaller N. ",
  "Column (3) adds state and year fixed effects with SEs clustered by state and year. ",
  "Column (4) adds the two-stage multilayer community (Leiden CPM within year, Jaccard-tracked across years) ",
  "as a third fixed effect and clusters SEs by state, year and community; state-years that are singletons within a fixed effect are dropped. ",
  "Community-clustered SEs without community fixed effects are in Appendix Table~\\ref{tab:state_main_comm}. ")
hl_tex <- c(
  "\\begin{table}[p]", "\\scriptsize", "\\begin{center}", "\\begin{tabular}{l c c c c}", "\\toprule",
  " & (1) IV + trend & (2) + Ctrls, no FE & (3) TWFE & (4) 3-way FE \\\\", "\\midrule",
  hl_body(hl_panels, c(1, 2, 3, 5), 5), "\\bottomrule", hl_note(5, HL_DV_RAW, HL_COLS_14),
  "\\end{tabular}", "\\end{center}",
  "\\caption{State-level regressions for fatal conflict and churn. Fatal MID onset in $t-1$, churn over all dyads (Panel A) and over peace dyads (Panel B).}",
  "\\label{tab:state_main}", "\\end{table}")
hl_iv_rows <- hl_tex[grepl("^Log\\(1 \\+ fatal MID onset", hl_tex)]
stopifnot(length(hl_iv_rows) == 2, all(lengths(regmatches(hl_iv_rows, gregexpr("&", hl_iv_rows))) == 4),
          sum(grepl("^Num. groups: Communities", hl_tex)) == 2, !any(grepl("artner war|PWE", hl_tex)))
writeLines(hl_tex, file.path(TABLE_DIR, "table3_state_main.tex"))
cat("  tab:state_main -> table3_state_main.tex (fatal MID onset t-1; Panels A-B; ladder columns 1, 2, 3, 5)\n")

hl_tex_comm <- c(
  "\\begin{table}[p]", "\\scriptsize", "\\begin{center}", "\\begin{tabular}{l c c c}", "\\toprule",
  " & (3) TWFE & (4) TWFE + Comm cluster & (5) 3-way FE \\\\", "\\midrule",
  hl_body(hl_panels, 3:5, 4), "\\bottomrule",
  hl_note(4, HL_DV_RAW, paste0("All columns include the controls and state and year fixed effects. ",
                               "Column (3) repeats main-text Table~\\ref{tab:state_main} column (3); column (4) keeps it and clusters SEs by state, year and community. ",
                               "Column (5) adds the community as a third fixed effect with three-way clustering (main-text column 4); ",
                               "state-years that are singletons within a fixed effect are dropped. ")),
  "\\end{tabular}", "\\end{center}",
  "\\caption{Appendix. State-level regressions with community fixed effects, companion to main-text Table~\\ref{tab:state_main}.}",
  "\\label{tab:state_main_comm}", "\\end{table}")
stopifnot(sum(grepl("^Log\\(1 \\+ fatal MID onset", hl_tex_comm)) == 2, sum(grepl("^Num. groups: Communities", hl_tex_comm)) == 2)
writeLines(hl_tex_comm, file.path(APPX_DIR, "tableA_state_main_comm.tex"))
cat("  tab:state_main_comm -> results/appendix/tableA_state_main_comm.tex (headline, columns 3-5)\n")

hl_tex_log <- c(
  "\\begin{table}[p]", "\\scriptsize", "\\begin{center}", "\\begin{tabular}{l c c c c c}", "\\toprule",
  " & (1) IV + trend & (2) + Ctrls, no FE & (3) TWFE & (4) TWFE + Comm cluster & (5) 3-way FE \\\\", "\\midrule",
  hl_body(hl_panels_log, 1:5, 6), "\\bottomrule",
  hl_note(6, "DV is $\\log(1+\\text{churn})$. ", paste0("Columns as in Table~\\ref{tab:state_main}; column (5) adds community fixed effects. ")),
  "\\end{tabular}", "\\end{center}",
  "\\caption{Appendix. State-level regressions with $\\log(1+\\text{churn})$ as the dependent variable, companion to main-text Table~\\ref{tab:state_main}.}",
  "\\label{tab:state_main_logchurn}", "\\end{table}")
stopifnot(sum(grepl("^Log\\(1 \\+ fatal MID onset", hl_tex_log)) == 2)
writeLines(hl_tex_log, file.path(APPX_DIR, "tableA_state_main_logchurn.tex"))
cat("  tab:state_main_logchurn -> results/appendix/tableA_state_main_logchurn.tex (headline, log churn)\n")

write_csv(bind_rows(lapply(names(hl_panels), function(pn) bind_rows(lapply(seq_along(hl_panels[[pn]]), function(m)
  broom::tidy(hl_panels[[pn]][[m]]) |> mutate(panel = pn, column = m, iv = HL_IV,
                                             dv = if (pn == "A") "churn" else "churn_peace_ft1",
                                             nobs = nobs(hl_panels[[pn]][[m]])))))),
  file.path(OUTPUT_DIR, "table3_state_main_coefs.csv"))

# [2026-10-12] Appendix Table A1 (fatal-MID ladder, tableA1_state_fatal.tex) is retired:
# the fatal results are Panels C-D of the main table.

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
    "Appendix. Three-way fixed effects (state $+$ year $+$ two-stage ",
    "multilayer community): churn (not logged) on log MID measures, ",
    "full controls."),
  label   = "tab:state_appx_3wfe",
  note_extra = paste0(
    "DV is churn (not logged). Three-way fixed effects (state, ",
    "year, two-stage multilayer community); SEs clustered three ways. ",
    "Sample restricted to state-years with non-missing Leiden community ",
    "assignment. Columns 1--2 use onset IVs; columns 3--4 use ongoing ",
    "MID counts."),
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
  "DV is churn (not logged). ",
  "Coefficients and standard errors are multiplied by 100. ",
  "Two-way fixed effects (state $+$ year), two-way clustered SEs ",
  "(state $+$ year). Each column reports the IV at $t{-}1$ or $t{-}2$ ",
  "with the full controls covariate set (",
  "democracy, autocracy, CINC, major power, isolated at $t{-}1$). ",
  "COW system, 1816--2014. ",
  "$^{***}p<0.001$; $^{**}p<0.01$; $^{*}p<0.05$; $^{\\dagger}p<0.1$.}")

texreg(
  scale_tr(list(
    lag_models[["mid_onset_count_lag1_full"]],     lag_models[["mid_onset_count_lag2_full"]],
    lag_models[["log_mid_onset_count_lag1_full"]], lag_models[["log_mid_onset_count_lag2_full"]],
    lag_models[["fatal_mid_onset_count_lag1_full"]], lag_models[["fatal_mid_onset_count_lag2_full"]],
    lag_models[["log_fatal_mid_onset_count_lag1_full"]], lag_models[["log_fatal_mid_onset_count_lag2_full"]]
  )),
  custom.model.names = c(
    "MID $t{-}1$", "MID $t{-}2$",
    "log MID $t{-}1$", "log MID $t{-}2$",
    "Fatal $t{-}1$", "Fatal $t{-}2$",
    "log Fatal $t{-}1$", "log Fatal $t{-}2$"
  ),
  custom.coef.map = lag_coef_map,
  caption = "Appendix. Lagged Conflict Onset and Churn.",
  label = "tab:lag_models",
  digits = 2,
  custom.note = lag_note,
  stars = c(0.001, 0.01, 0.05, 0.1),
  symbol = "$^{\\dagger}$",
  float.pos = "htbp",
  booktabs = TRUE,
  use.packages = FALSE,
  file = file.path(APPX_DIR, "tableA3_state_lag_models.tex")
)
cat("  Appendix (lag models) -> tableA3_state_lag_models.tex\n")

# [2026-10-08] Within (proj model) R^2 is tiny in the appendix TWFE tables and prints as 0.00 at two
# decimals. Rescaling cannot change R^2, so write those cells as "<0.01" instead of a misleading 0.00.
for (f in file.path(APPX_DIR, c("tableA2_state_3wfe.tex", "tableA3_state_lag_models.tex"))) {
  tx <- readLines(f)
  k <- grepl("R\\$\\^2\\$ \\(proj model\\)", tx)
  tx[k] <- gsub("$0.00$", "$<0.01$", tx[k], fixed = TRUE)
  writeLines(tx, f)
  stopifnot(!any(grepl("$0.00$", tx, fixed = TRUE)))
}

# --------------------------------------------------------------------------- #
# 5. Clustering Robustness (onset IVs — appendix)
# --------------------------------------------------------------------------- #
# Full specification grid across all 4 IVs × weighting is in 03_robustness.R.
# Here: onset IVs only, varying SE clustering scheme.
# --------------------------------------------------------------------------- #

cat("\n========== CLUSTERING ROBUSTNESS ==========\n")

# One-way (state only)
m_onset_full_1way <- feols(churn ~ mid_onset_count +
                             democracy + autocracy + cinc + major_power + is_isolated_t1 | ccode + year,
                           data = state_churn, vcov = ~ccode)

m_fonset_full_1way <- feols(churn ~ fatal_mid_onset_count +
                              democracy + autocracy + cinc + major_power + is_isolated_t1 | ccode + year,
                            data = state_churn, vcov = ~ccode)

# Leiden CPM variants
state_churn_leiden <- state_churn |> filter(!is.na(leiden_comm))

m_onset_full_leiden <- feols(churn ~ mid_onset_count +
                               democracy + autocracy + cinc + major_power + is_isolated_t1 | ccode + year,
                             data = state_churn_leiden, vcov = ~leiden_comm)

m_onset_full_leiden_yr <- feols(churn ~ mid_onset_count +
                                  democracy + autocracy + cinc + major_power + is_isolated_t1 | ccode + year,
                                data = state_churn_leiden, vcov = ~leiden_comm + year)

m_onset_full_leiden_st <- feols(churn ~ mid_onset_count +
                                  democracy + autocracy + cinc + major_power + is_isolated_t1 | ccode + year,
                                data = state_churn_leiden, vcov = ~leiden_comm + ccode)

m_onset_full_leiden_yr_st <- feols(churn ~ mid_onset_count +
                                   democracy + autocracy + cinc + major_power + is_isolated_t1 | ccode + year,
                                 data = state_churn_leiden, vcov = ~leiden_comm + year + ccode)


m_fonset_full_leiden <- feols(churn ~ fatal_mid_onset_count +
                               democracy + autocracy + cinc + major_power + is_isolated_t1 | ccode + year,
                             data = state_churn_leiden, vcov = ~leiden_comm)

m_fonset_full_leiden_yr <- feols(churn ~ fatal_mid_onset_count +
                                  democracy + autocracy + cinc + major_power + is_isolated_t1 | ccode + year,
                                data = state_churn_leiden, vcov = ~leiden_comm + year)

m_fonset_full_leiden_st <- feols(churn ~ fatal_mid_onset_count +
                                  democracy + autocracy + cinc + major_power + is_isolated_t1 | ccode + year,
                                data = state_churn_leiden, vcov = ~leiden_comm + ccode)

m_fonset_full_leiden_yr_st <- feols(churn ~ fatal_mid_onset_count +
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

# [2026-10-12] Partner war exposure (H3) is retired. Drop its columns from the saved panels so no downstream
# consumer can pick them up from state_level_results.rds or state_churn_data.csv.
PWE_COLS <- c("partner_war_frac", "partner_war_frac_weighted", "partner_war_count",
              "partner_fatal_war_frac", "partner_fatal_war_frac_weighted", "partner_fatal_war_count",
              "partner_war_frac_weighted_peace")
state_churn_export <- state_churn |> select(-any_of(PWE_COLS))
stopifnot(!any(PWE_COLS %in% names(state_churn_export)))

results <- list(
  # Main + appendix TWFE models, keyed by spec key
  # (main__*, a1__*, a2__*, a3__*, a4__*)
  twfe = twfe_models,
  # Lag models (onset IVs at lag1/lag2, raw churn DV)
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
  data = state_churn_export
)

saveRDS(results, file.path(OUTPUT_DIR, "state_level_results.rds"))

# --- CSV exports ---
write_csv(state_churn_export, file.path(OUTPUT_DIR, "state_churn_data.csv"))

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
# (churn, raw + logged), IV (MID onset and fatal MID onset, raw + logged),
# and the state attributes (democracy, autocracy, CINC, major power,
# is_isolated_t1) that make up the control set.
#
# DROPPED: partner war exposure [2026-10-12: no longer in any main-text
# model; the variables are still saved in state_churn_data.csv]. Ongoing
# MID count IVs (mid_count, fatal_mid_count) — main spec uses onset, count
# is appendix-only via 03. Lagged versions of any IV — not used in Table 3.
sumstat_vars <- c(
  # DV (Table 3)
  "churn", "churn_peace_ft1", "log_churn",
  # IV (Table 3 main [2026-10-12]: fatal MID onset at t-1; MID onset in the appendix)
  "fatal_onset_lag1", "log_fatal_onset_lag1",
  "mid_onset_count", "log_mid_onset_count",
  "fatal_mid_onset_count", "log_fatal_mid_onset_count",
  # State-level controls
  "democracy", "autocracy", "cinc", "major_power", "is_isolated_t1"
)

sumstat_labels <- c(
  churn                            = "Network Churn",
  churn_peace_ft1                  = "Network Churn, peace dyads",
  log_churn                        = "$\\log(1 + \\text{Churn})$",
  fatal_onset_lag1                 = "Fatal MID Onset Count$_{t-1}$",
  log_fatal_onset_lag1             = "$\\log(1 + \\text{Fatal MID Onset})_{t-1}$",
  mid_onset_count                  = "MID Onset Count",
  log_mid_onset_count              = "$\\log(1 + \\text{MID Onset})$",
  fatal_mid_onset_count            = "Fatal MID Onset Count",
  log_fatal_mid_onset_count        = "$\\log(1 + \\text{Fatal MID Onset})$",
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
    "IVs (MID onset, fatal MID onset). ",
    "$\\mathbb{1}\\{$Isolated at $t{-}1\\}$ flags state-years with zero ",
    "partners at the prior year. COW system, 1816--2014.}"
  ),
  "\\end{table}"
)
writeLines(tex_lines, file.path(TABLE_DIR, "summary_stats_state.tex"))
cat("  LaTeX -> summary_stats_state.tex (results/main)\n")

cat("\nDone. State-level analysis complete.\n")