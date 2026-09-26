###############################################################################
# 02_compute_network_measures.R
# Compute network divergence measures for system-level and state-level analysis
#
# MAIN SPECIFICATION (Option 3):
#   Uses only layers available for the full panel (alliances + IGO).
#   This avoids mechanical divergence spikes when layers enter/exit.
#
# APPENDIX SPECIFICATION (Option 2):
#   Computes layer-specific divergence measures for each layer individually,
#   then averages across available layers per year. Also computes measures
#   for the full composite (all layers) for comparison.
#
# System-level measures (year t vs year t-1):
#   - KL divergence (weighted)
#   - Jensen-Shannon divergence
#   - Edge turnover (count and %)
#   - Network density change
#   - Centralization change
#   - Modularity change
#   - Wasserstein distance on degree distributions
#
# State-level measures:
#   - Weighted Jaccard similarity of partner set
#   - Churn (1 - Jaccard)
#   - Binary Jaccard / churn
#   - Degree change
###############################################################################

library(tidyverse)
library(igraph)
library(data.table)
select <- dplyr::select

OUTPUT_DIR <- here::here("data", "processed")

# --------------------------------------------------------------------------- #
# 0. Load network data
# --------------------------------------------------------------------------- #

cat("Loading network data...\n")
dyad_data <- readRDS(file.path(OUTPUT_DIR, "dyad_year_network_layers.rds"))
metadata <- readRDS(file.path(OUTPUT_DIR, "network_metadata.rds"))

years <- sort(unique(dyad_data$year))

# --------------------------------------------------------------------------- #
# 0a. Define layer sets
# --------------------------------------------------------------------------- #

# Full-panel layers: available 1816-2018 (ATOP) / 1816-2014 (IGO)
# These are the MAIN specification layers.
MAIN_LAYERS <- c("layer_offensive", "layer_defensive", "layer_igo")
MAIN_WEIGHTED <- paste0(MAIN_LAYERS, "_weighted")

# All layers (for appendix composite)
ALL_LAYERS <- metadata$layers   # c("layer_offensive", "layer_defensive", "layer_igo", "layer_trade", "layer_dca")
ALL_WEIGHTED <- metadata$weighted_layers

# Individual layers for appendix layer-specific analysis
INDIVIDUAL_LAYERS <- ALL_LAYERS

cat("Main specification layers:", paste(MAIN_LAYERS, collapse = ", "), "\n")
cat("All layers:", paste(ALL_LAYERS, collapse = ", "), "\n\n")

# --------------------------------------------------------------------------- #
# 0b. Compute main-specification edge weight
#     (uses only alliance + IGO layers, normalized within year)
# --------------------------------------------------------------------------- #

range01 <- function(x) {
  rng <- max(x, na.rm = TRUE) - min(x, na.rm = TRUE)
  if (rng == 0) return(rep(0, length(x)))
  (x - min(x, na.rm = TRUE)) / rng
}

cat("Computing main-specification edge weight (alliance + IGO only)...\n")

# Raw sum of main layers
main_mat <- as.matrix(dyad_data[, MAIN_LAYERS, drop = FALSE])
main_mat[is.na(main_mat)] <- 0
dyad_data$main_raw <- rowSums(main_mat)

# Normalize within year to [0, 1]
dyad_data <- dyad_data |>
  group_by(year) |>
  mutate(
    main_weight = if (max(main_raw, na.rm = TRUE) > 0) {
      main_raw / max(main_raw, na.rm = TRUE)
    } else {
      0
    }
  ) |>
  ungroup() |>
  mutate(main_weight = round(main_weight, 4))

cat(sprintf("  main_weight: pct > 0 = %.1f%%\n",
            100 * mean(dyad_data$main_weight > 0, na.rm = TRUE)))

# Also do contiguity-weighted version for main layers
main_w_mat <- as.matrix(dyad_data[, MAIN_WEIGHTED, drop = FALSE])
main_w_mat[is.na(main_w_mat)] <- 0
n_main <- length(MAIN_LAYERS)
dyad_data$main_weighted_raw <- rowSums(main_w_mat) / n_main

# Shift residuals so all values are non-negative: tie_shifted = tie - min(tie)
# This preserves relative distances between dyads while making all values >= 0.
# Then range01 normalize within year.
dyad_data <- dyad_data |>
  group_by(year) |>
  mutate(
    main_weight_contig_shifted = main_weighted_raw - min(main_weighted_raw, na.rm = TRUE),
    main_weight_contig = range01(main_weight_contig_shifted)
  ) |>
  ungroup() |>
  mutate(main_weight_contig = round(main_weight_contig, 4))

cat(sprintf("  main_weight_contig: pct > 0 = %.1f%%, pct == 0 = %.1f%%\n",
            100 * mean(dyad_data$main_weight_contig > 0, na.rm = TRUE),
            100 * mean(dyad_data$main_weight_contig == 0, na.rm = TRUE)))
cat(sprintf("  main_weighted_raw range: [%.4f, %.4f]\n",
            min(dyad_data$main_weighted_raw, na.rm = TRUE),
            max(dyad_data$main_weighted_raw, na.rm = TRUE)))

# Save updated dyad_data with main_weight and main_weight_contig
# so downstream scripts (02_state_level_analysis.R) can access them
cat("Saving updated dyad_year_network_layers.rds (with main_weight_contig)...\n")
saveRDS(dyad_data, file.path(OUTPUT_DIR, "dyad_year_network_layers.rds"))

# --------------------------------------------------------------------------- #
# 1. Helper functions: divergence measures
# --------------------------------------------------------------------------- #

#' KL Divergence: D_KL(P || Q) = sum(P * log(P/Q))
#' Handles zeros by adding small epsilon
kl_divergence <- function(p, q, eps = 1e-10) {
  p <- p + eps
  q <- q + eps
  p <- p / sum(p)
  q <- q / sum(q)
  sum(p * log(p / q))
}

#' Jensen-Shannon Divergence: symmetric version of KL
js_divergence <- function(p, q, eps = 1e-10) {
  p <- p + eps
  q <- q + eps
  p <- p / sum(p)
  q <- q / sum(q)
  m <- 0.5 * (p + q)
  0.5 * kl_divergence(p, m, eps = 0) + 0.5 * kl_divergence(q, m, eps = 0)
}

#' Wasserstein distance (1D, for degree distributions)
wasserstein_dist <- function(p, q) {
  all_vals <- sort(unique(c(p, q)))
  cdf_p <- ecdf(p)
  cdf_q <- ecdf(q)
  sum(abs(cdf_p(all_vals) - cdf_q(all_vals))) * (max(all_vals) - min(all_vals)) / length(all_vals)
}

#' Compute all divergence measures between two years of edge data
#' Returns a one-row tibble
compute_divergence <- function(edges_t, edges_t1, weight_var) {

  # Build igraph objects
  make_graph <- function(edges, wvar) {
    if (nrow(edges) == 0) return(make_empty_graph(directed = FALSE))
    g <- graph_from_data_frame(edges |> select(ccode1, ccode2), directed = FALSE)
    E(g)$weight <- edges[[wvar]]
    g
  }

  g_t  <- make_graph(edges_t, weight_var)
  g_t1 <- make_graph(edges_t1, weight_var)

  n_t  <- vcount(g_t)
  n_t1 <- vcount(g_t1)

  # Network statistics
  density_t  <- if (n_t > 1) edge_density(g_t) else 0
  density_t1 <- if (n_t1 > 1) edge_density(g_t1) else 0

  centralization_t  <- if (n_t > 2) centr_degree(g_t)$centralization else NA_real_
  centralization_t1 <- if (n_t1 > 2) centr_degree(g_t1)$centralization else NA_real_

  modularity_t <- if (n_t > 2 & ecount(g_t) > 0) {
    modularity(cluster_louvain(g_t, weights = E(g_t)$weight))
  } else NA_real_
  modularity_t1 <- if (n_t1 > 2 & ecount(g_t1) > 0) {
    modularity(cluster_louvain(g_t1, weights = E(g_t1)$weight))
  } else NA_real_

  # KL and JS divergence on edge weight distributions
  merged <- full_join(
    edges_t |> select(ccode1, ccode2, w_t = !!sym(weight_var)),
    edges_t1 |> select(ccode1, ccode2, w_t1 = !!sym(weight_var)),
    by = c("ccode1", "ccode2")
  ) |> mutate(w_t = replace_na(w_t, 0), w_t1 = replace_na(w_t1, 0))

  if (min(merged$w_t) < 0) merged$w_t <- merged$w_t - min(merged$w_t)
  if (min(merged$w_t1) < 0) merged$w_t1 <- merged$w_t1 - min(merged$w_t1)

  kl_div <- kl_divergence(merged$w_t, merged$w_t1)
  js_div <- js_divergence(merged$w_t, merged$w_t1)

  # Edge turnover
  dyad_key <- function(df) {
    df |> mutate(dyad = paste(pmin(ccode1, ccode2), pmax(ccode1, ccode2), sep = "-")) |> pull(dyad)
  }
  edges_t_set  <- dyad_key(edges_t)
  edges_t1_set <- dyad_key(edges_t1)

  edges_added   <- length(setdiff(edges_t_set, edges_t1_set))
  edges_removed <- length(setdiff(edges_t1_set, edges_t_set))
  edge_turnover <- edges_added + edges_removed
  edge_turnover_pct <- edge_turnover / max(length(union(edges_t_set, edges_t1_set)), 1)

  # Degree distribution divergence
  deg_t  <- if (n_t > 0) degree(g_t) else numeric(0)
  deg_t1 <- if (n_t1 > 0) degree(g_t1) else numeric(0)

  if (length(deg_t) > 0 & length(deg_t1) > 0) {
    degree_wasserstein <- wasserstein_dist(deg_t, deg_t1)
    max_deg <- max(c(deg_t, deg_t1))
    hist_t  <- tabulate(deg_t + 1, nbins = max_deg + 1)
    hist_t1 <- tabulate(deg_t1 + 1, nbins = max_deg + 1)
    degree_kl <- kl_divergence(hist_t, hist_t1)
  } else {
    degree_wasserstein <- NA_real_
    degree_kl <- NA_real_
  }

  tibble(
    n_nodes          = n_t,
    n_edges          = ecount(g_t),
    density          = density_t,
    centralization   = centralization_t,
    modularity       = modularity_t,
    kl_divergence    = kl_div,
    js_divergence    = js_div,
    edges_added      = edges_added,
    edges_removed    = edges_removed,
    edge_turnover    = edge_turnover,
    edge_turnover_pct = edge_turnover_pct,
    density_change   = density_t - density_t1,
    centralization_change = centralization_t - centralization_t1,
    modularity_change     = modularity_t - modularity_t1,
    degree_wasserstein    = degree_wasserstein,
    degree_kl             = degree_kl
  )
}

# --------------------------------------------------------------------------- #
# 2. MAIN SPECIFICATION: System-level measures using fixed layer set
#    Edge weight: main_weight (unweighted Manhattan, alliance + IGO layers only)
# --------------------------------------------------------------------------- #

cat("=" |> strrep(70), "\n")
cat("MAIN SPECIFICATION: Alliance + IGO layers, unweighted Manhattan\n")
cat("=" |> strrep(70), "\n\n")

MAIN_WEIGHT_VAR <- "main_weight"

system_main <- list()

for (i in seq_along(years)) {
  yr <- years[i]
  if (i %% 20 == 0) cat(sprintf("  Processing %d...\n", yr))

  edges_t <- dyad_data |>
    filter(year == yr, .data[[MAIN_WEIGHT_VAR]] > 0) |>
    select(ccode1, ccode2, !!sym(MAIN_WEIGHT_VAR))

  if (i == 1) {
    # First year: store static measures only
    g_t <- if (nrow(edges_t) > 0) {
      g <- graph_from_data_frame(edges_t |> select(ccode1, ccode2), directed = FALSE)
      E(g)$weight <- edges_t[[MAIN_WEIGHT_VAR]]
      g
    } else make_empty_graph(directed = FALSE)

    system_main[[i]] <- tibble(
      year = yr,
      n_nodes = vcount(g_t), n_edges = ecount(g_t),
      density = if (vcount(g_t) > 1) edge_density(g_t) else 0,
      centralization = if (vcount(g_t) > 2) centr_degree(g_t)$centralization else NA_real_,
      modularity = if (vcount(g_t) > 2 & ecount(g_t) > 0) modularity(cluster_louvain(g_t, weights = E(g_t)$weight)) else NA_real_,
      kl_divergence = NA_real_, js_divergence = NA_real_,
      edges_added = NA_integer_, edges_removed = NA_integer_,
      edge_turnover = NA_integer_, edge_turnover_pct = NA_real_,
      density_change = NA_real_, centralization_change = NA_real_,
      modularity_change = NA_real_, degree_wasserstein = NA_real_,
      degree_kl = NA_real_
    )
    next
  }

  yr_prev <- years[i - 1]
  edges_t1 <- dyad_data |>
    filter(year == yr_prev, .data[[MAIN_WEIGHT_VAR]] > 0) |>
    select(ccode1, ccode2, !!sym(MAIN_WEIGHT_VAR))

  measures <- compute_divergence(edges_t, edges_t1, MAIN_WEIGHT_VAR)
  system_main[[i]] <- bind_cols(tibble(year = yr), measures)
}

system_main_ts <- bind_rows(system_main)
cat(sprintf("  Main specification: %d years computed\n\n", nrow(system_main_ts)))

# --------------------------------------------------------------------------- #
# 2b. ROBUSTNESS: System-level measures across all 5 weighting variants
# --------------------------------------------------------------------------- #
#
# Script 01 produces 5 composite edge weights (all 5 layers):
#   manhattan              — unweighted sum, range01
#   manhattan_weighted     — contiguity-residualized, Manhattan norm
#   manhattan_contig_euclid — contiguity-residualized, Euclidean norm
#   manhattan_dist         — distance-residualized, Manhattan norm
#   manhattan_dist_euclid  — distance-residualized, Euclidean norm
#
# For each: compute year-over-year KL, JS, edge turnover, density,
# centralization, modularity, Wasserstein on degree distributions.
# --------------------------------------------------------------------------- #

cat("=" |> strrep(70), "\n")
cat("ROBUSTNESS: System-level measures across 5 weighting variants\n")
cat("=" |> strrep(70), "\n\n")

WEIGHT_VARIANTS <- c(
  "manhattan",
  "manhattan_weighted",
  "manhattan_contig_euclid",
  "manhattan_dist",
  "manhattan_dist_euclid"
)

# Check which variants actually exist in the data
available_variants <- intersect(WEIGHT_VARIANTS, names(dyad_data))
missing_variants <- setdiff(WEIGHT_VARIANTS, names(dyad_data))
if (length(missing_variants) > 0) {
  cat(sprintf("  WARNING: Missing weight variants: %s\n", paste(missing_variants, collapse = ", ")))
}
cat(sprintf("  Computing for %d variants: %s\n\n", length(available_variants),
            paste(available_variants, collapse = ", ")))

system_robustness_list <- list()

for (wv in available_variants) {
  cat(sprintf("  Processing variant: %s\n", wv))

  variant_results <- list()

  for (i in seq_along(years)) {
    yr <- years[i]

    edges_t <- dyad_data |>
      filter(year == yr, .data[[wv]] > 0) |>
      select(ccode1, ccode2, !!sym(wv))

    if (i == 1) {
      # First year: static measures only
      g_t <- if (nrow(edges_t) > 0) {
        g <- graph_from_data_frame(edges_t |> select(ccode1, ccode2), directed = FALSE)
        E(g)$weight <- edges_t[[wv]]
        g
      } else make_empty_graph(directed = FALSE)

      variant_results[[i]] <- tibble(
        year = yr, weight_variant = wv,
        n_nodes = vcount(g_t), n_edges = ecount(g_t),
        density = if (vcount(g_t) > 1) edge_density(g_t) else 0,
        centralization = if (vcount(g_t) > 2) centr_degree(g_t)$centralization else NA_real_,
        modularity = if (vcount(g_t) > 2 & ecount(g_t) > 0) modularity(cluster_louvain(g_t, weights = E(g_t)$weight)) else NA_real_,
        kl_divergence = NA_real_, js_divergence = NA_real_,
        edges_added = NA_integer_, edges_removed = NA_integer_,
        edge_turnover = NA_integer_, edge_turnover_pct = NA_real_,
        density_change = NA_real_, centralization_change = NA_real_,
        modularity_change = NA_real_, degree_wasserstein = NA_real_,
        degree_kl = NA_real_
      )
      next
    }

    yr_prev <- years[i - 1]
    edges_t1 <- dyad_data |>
      filter(year == yr_prev, .data[[wv]] > 0) |>
      select(ccode1, ccode2, !!sym(wv))

    measures <- compute_divergence(edges_t, edges_t1, wv)
    variant_results[[i]] <- bind_cols(tibble(year = yr, weight_variant = wv), measures)
  }

  system_robustness_list[[wv]] <- bind_rows(variant_results)
  cat(sprintf("    %s: %d years, mean KL = %.4f\n", wv,
              nrow(system_robustness_list[[wv]]),
              mean(system_robustness_list[[wv]]$kl_divergence, na.rm = TRUE)))
}

system_robustness_ts <- bind_rows(system_robustness_list)
cat(sprintf("\n  Robustness total: %d variant-year observations\n\n",
            nrow(system_robustness_ts)))

# --------------------------------------------------------------------------- #
# 3. APPENDIX: Layer-specific divergence measures
# --------------------------------------------------------------------------- #

cat("=" |> strrep(70), "\n")
cat("APPENDIX: Layer-specific divergence measures\n")
cat("=" |> strrep(70), "\n\n")

# For each individual layer, compute KL and JS divergence year-over-year
layer_specific_results <- list()

for (lv in INDIVIDUAL_LAYERS) {
  cat(sprintf("  Computing measures for %s...\n", lv))

  # Check which years have nonzero values for this layer
  layer_years <- dyad_data |>
    group_by(year) |>
    summarise(has_data = any(.data[[lv]] > 0, na.rm = TRUE), .groups = "drop") |>
    filter(has_data) |>
    pull(year) |>
    sort()

  if (length(layer_years) == 0) {
    cat(sprintf("    Skipping %s: no nonzero data in any year\n", lv))
    next
  }

  cat(sprintf("    Coverage: %d-%d (%d years)\n",
              min(layer_years), max(layer_years), length(layer_years)))

  for (i in seq_along(layer_years)) {
    yr <- layer_years[i]

    if (i == 1) {
      layer_specific_results[[length(layer_specific_results) + 1]] <- tibble(
        layer = lv, year = yr,
        kl_divergence = NA_real_, js_divergence = NA_real_,
        edge_turnover_pct = NA_real_, n_edges = NA_integer_
      )
      next
    }

    yr_prev <- layer_years[i - 1]

    edges_t <- dyad_data |>
      filter(year == yr, .data[[lv]] > 0) |>
      select(ccode1, ccode2, !!sym(lv))
    edges_t1 <- dyad_data |>
      filter(year == yr_prev, .data[[lv]] > 0) |>
      select(ccode1, ccode2, !!sym(lv))

    # KL and JS on edge weight distributions
    merged <- full_join(
      edges_t |> rename(w_t = !!sym(lv)),
      edges_t1 |> rename(w_t1 = !!sym(lv)),
      by = c("ccode1", "ccode2")
    ) |> mutate(w_t = replace_na(w_t, 0), w_t1 = replace_na(w_t1, 0))

    if (min(merged$w_t) < 0) merged$w_t <- merged$w_t - min(merged$w_t)
    if (min(merged$w_t1) < 0) merged$w_t1 <- merged$w_t1 - min(merged$w_t1)

    kl_val <- kl_divergence(merged$w_t, merged$w_t1)
    js_val <- js_divergence(merged$w_t, merged$w_t1)

    # Edge turnover for this layer
    dyad_key <- function(df) {
      df |> mutate(dyad = paste(pmin(ccode1, ccode2), pmax(ccode1, ccode2), sep = "-")) |> pull(dyad)
    }
    et_set  <- dyad_key(edges_t)
    et1_set <- dyad_key(edges_t1)
    turnover <- length(setdiff(et_set, et1_set)) + length(setdiff(et1_set, et_set))
    turnover_pct <- turnover / max(length(union(et_set, et1_set)), 1)

    layer_specific_results[[length(layer_specific_results) + 1]] <- tibble(
      layer = lv, year = yr,
      kl_divergence = kl_val, js_divergence = js_val,
      edge_turnover_pct = turnover_pct, n_edges = nrow(edges_t)
    )
  }
}

layer_specific_ts <- bind_rows(layer_specific_results)
cat(sprintf("  Layer-specific: %d layer-year observations\n\n", nrow(layer_specific_ts)))

# Compute cross-layer average divergence per year (appendix composite)
appendix_avg <- layer_specific_ts |>
  filter(!is.na(kl_divergence)) |>
  group_by(year) |>
  summarise(
    n_layers_available   = n(),
    avg_kl_divergence    = mean(kl_divergence, na.rm = TRUE),
    avg_js_divergence    = mean(js_divergence, na.rm = TRUE),
    avg_edge_turnover_pct = mean(edge_turnover_pct, na.rm = TRUE),
    layers_present       = paste(sort(layer), collapse = "+"),
    .groups = "drop"
  )

cat("  Cross-layer average computed.\n")
cat(sprintf("  Layer availability transitions:\n"))
appendix_avg |>
  distinct(n_layers_available, layers_present) |>
  arrange(n_layers_available) |>
  print(n = Inf)
cat("\n")

# --------------------------------------------------------------------------- #
# 4. MAIN: State-level measures using fixed layer set (alliance + IGO)
# --------------------------------------------------------------------------- #

cat("=" |> strrep(70), "\n")
cat("MAIN: State-level churn (alliance + IGO layers only)\n")
cat("=" |> strrep(70), "\n\n")

# Undirected state ties using main_weight (unweighted Manhattan)
state_ties <- dyad_data |>
  select(ccode1, ccode2, year, w = !!sym(MAIN_WEIGHT_VAR)) |>
  bind_rows(
    dyad_data |>
      select(ccode1 = ccode2, ccode2 = ccode1, year, w = !!sym(MAIN_WEIGHT_VAR))
  ) |>
  group_by(ccode1, ccode2, year) |>
  summarise(w = max(w, na.rm = TRUE), .groups = "drop")

# --------------------------------------------------------------------------- #
# State-year universe: every (ccode, year) appearing in any dyad-row in
# dyad_data — equivalent to COW system membership for that year. This is
# the FULL state-year panel (not conditional on having positive ties), so
# isolated state-years (degree = 0) are present and get a defined churn
# value rather than being silently dropped.
#
# Four-case decision rule for churn:
#   (a) degree_t == 0 AND degree_t-1 == 0:  churn = 0  (vacuously identical
#       partner sets; Jaccard convention J(∅, ∅) = 1, so churn = 1 - 1 = 0).
#       Tagged case = "both_isolated" so robustness filters can drop them.
#   (b) degree_t == 0 XOR degree_t-1 == 0:  churn = 1  (total replacement
#       — one set is a strict subset of the other, intersection is empty).
#       Tagged case = "total_gain" or "total_loss" accordingly.
#   (c) Both nonzero:                        churn = 1 - weighted Jaccard.
#       Tagged case = "normal".
#
# Years must be CONSECUTIVE in the state's COW panel (yr_t - yr_t1 == 1).
# Non-consecutive transitions (state exits and re-enters COW) are dropped
# because the "year-over-year" interpretation breaks down.
# --------------------------------------------------------------------------- #

state_universe <- dyad_data |>
  select(ccode = ccode1, year) |>
  bind_rows(dyad_data |> select(ccode = ccode2, year)) |>
  distinct(ccode, year) |>
  arrange(ccode, year)

# Compute degree per state-year on the FULL universe (zero-rows allowed).
# Left join state_universe ← state_ties (filtered to w > 0) so missing
# joins yield degree = 0 and weighted_degree = 0.
state_partner_long <- state_ties |>
  filter(w > 0) |>
  rename(ccode = ccode1, partner = ccode2)

state_year_summary <- state_universe |>
  left_join(
    state_partner_long |>
      group_by(ccode, year) |>
      summarise(
        degree          = n_distinct(partner),
        weighted_degree = sum(w, na.rm = TRUE),
        .groups = "drop"
      ),
    by = c("ccode", "year")
  ) |>
  mutate(
    degree          = replace_na(degree, 0L),
    weighted_degree = replace_na(weighted_degree, 0)
  )

cat(sprintf("  Full state-year universe: %d obs across %d states, %d years\n",
            nrow(state_universe),
            n_distinct(state_universe$ccode),
            n_distinct(state_universe$year)))
cat(sprintf("  State-years with degree > 0: %d (%.1f%%)\n",
            sum(state_year_summary$degree > 0),
            100 * mean(state_year_summary$degree > 0)))
cat(sprintf("  State-years with degree == 0 (isolated): %d (%.1f%%)\n",
            sum(state_year_summary$degree == 0),
            100 * mean(state_year_summary$degree == 0)))

# Compute churn year-over-year for each state.
state_churn_list <- list()
states <- sort(unique(state_universe$ccode))

for (st in states) {
  st_years <- state_year_summary |>
    filter(ccode == st) |>
    arrange(year)

  if (nrow(st_years) < 2) next

  st_ties <- state_partner_long |> filter(ccode == st)

  for (j in 2:nrow(st_years)) {
    yr_t  <- st_years$year[j]
    yr_t1 <- st_years$year[j - 1]

    # Only consecutive COW state-years.
    if (yr_t - yr_t1 != 1) next

    deg_t  <- st_years$degree[j]
    deg_t1 <- st_years$degree[j - 1]

    if (deg_t == 0 && deg_t1 == 0) {
      # (a) Both isolated. Jaccard(∅, ∅) = 1 by convention. Churn = 0.
      jaccard <- 1; churn <- 0
      jaccard_binary <- 1; churn_binary <- 0
      partners_gained <- 0L; partners_lost <- 0L
      case_label <- "both_isolated"

    } else if (deg_t == 0 || deg_t1 == 0) {
      # (b) Total replacement. Jaccard = 0, Churn = 1.
      jaccard <- 0; churn <- 1
      jaccard_binary <- 0; churn_binary <- 1
      if (deg_t == 0) {
        partners_gained <- 0L
        partners_lost   <- deg_t1
        case_label      <- "total_loss"
      } else {
        partners_gained <- deg_t
        partners_lost   <- 0L
        case_label      <- "total_gain"
      }

    } else {
      # (c) Both nonzero — weighted Jaccard on the union of partner sets.
      w_t  <- st_ties |> filter(year == yr_t)  |> select(partner, w)
      w_t1 <- st_ties |> filter(year == yr_t1) |> select(partner, w)

      w_both <- full_join(w_t, w_t1, by = "partner",
                          suffix = c("_t", "_t1")) |>
        mutate(w_t = replace_na(w_t, 0), w_t1 = replace_na(w_t1, 0))

      sum_min <- sum(pmin(w_both$w_t, w_both$w_t1))
      sum_max <- sum(pmax(w_both$w_t, w_both$w_t1))
      jaccard <- if (sum_max > 0) sum_min / sum_max else NA_real_
      churn   <- if (!is.na(jaccard)) 1 - jaccard else NA_real_

      partners_t  <- w_both$partner[w_both$w_t  > 0]
      partners_t1 <- w_both$partner[w_both$w_t1 > 0]
      bin_inter <- length(intersect(partners_t, partners_t1))
      bin_union <- length(union(partners_t, partners_t1))
      jaccard_binary <- if (bin_union > 0) bin_inter / bin_union else NA_real_
      churn_binary   <- if (!is.na(jaccard_binary)) 1 - jaccard_binary else NA_real_

      partners_gained <- length(setdiff(partners_t,  partners_t1))
      partners_lost   <- length(setdiff(partners_t1, partners_t))
      case_label      <- "normal"
    }

    state_churn_list[[length(state_churn_list) + 1]] <- tibble(
      ccode = st,
      year  = yr_t,
      jaccard = jaccard,
      churn   = churn,
      jaccard_binary = jaccard_binary,
      churn_binary   = churn_binary,
      partners_gained = partners_gained,
      partners_lost   = partners_lost,
      degree_t  = deg_t,
      degree_t1 = deg_t1,
      degree_change = deg_t - deg_t1,
      weighted_degree_t  = st_years$weighted_degree[j],
      weighted_degree_t1 = st_years$weighted_degree[j - 1],
      case = case_label
    )
  }
}

if (length(state_churn_list) == 0) {
  cat("  WARNING: No consecutive-year state churn pairs found.\n")
  state_churn <- tibble(
    ccode = integer(0), year = integer(0), jaccard = numeric(0),
    churn = numeric(0), jaccard_binary = numeric(0), churn_binary = numeric(0),
    partners_gained = integer(0), partners_lost = integer(0),
    degree_t = integer(0), degree_t1 = integer(0), degree_change = integer(0),
    weighted_degree_t = numeric(0), weighted_degree_t1 = numeric(0),
    case = character(0)
  )
} else {
  state_churn <- bind_rows(state_churn_list)
}

# Predetermined isolation indicators for use as regression controls.
# is_isolated_t1 == 1 when the state had no partners at t-1 — captures
# both `total_gain` (degree_t1 = 0, degree_t > 0, churn = 1) and
# `both_isolated` (churn = 0). The mean shift between these cases and
# `normal` would otherwise flatten the within-state IV->churn slope.
# Including is_isolated_t1 in the regression absorbs this level shift
# without conditioning on the outcome (degree_t1 is fixed before t).
#
# is_isolated_t exists only for diagnostic use; it's a function of the
# year-t network and so is endogenous to the DV — DON'T put it on the
# RHS of the main regression.
state_churn <- state_churn |>
  mutate(
    is_isolated_t1 = as.integer(degree_t1 == 0),
    is_isolated_t  = as.integer(degree_t  == 0)
  )

# Diagnostic: distribution of cases — useful for the appendix and to make
# the panel construction transparent. A state-year is in `state_churn` iff
# the state was in COW at both yr and yr-1.
case_counts <- state_churn |> count(case) |>
  mutate(pct = 100 * n / sum(n))
cat("  Churn case distribution:\n")
for (i in seq_len(nrow(case_counts))) {
  cat(sprintf("    %-15s : %6d (%.2f%%)\n",
              case_counts$case[i], case_counts$n[i], case_counts$pct[i]))
}

cat(sprintf("  State-level churn observations (full panel): %d\n\n",
            nrow(state_churn)))

# --------------------------------------------------------------------------- #
# 5. Merge war participation into state-level data
# --------------------------------------------------------------------------- #

state_war <- dyad_data |>
  filter(cowmidongoing == 1) |>
  select(ccode1, ccode2, year, cowmidfatal) |>
  pivot_longer(cols = c(ccode1, ccode2), values_to = "ccode") |>
  group_by(ccode, year) |>
  summarise(
    war_participation = 1,
    fatal_mid_participation = as.numeric(any(cowmidfatal == 1, na.rm = TRUE)),
    mid_count = n_distinct(paste(ccode, year)),  # count of MID dyads per state-year
    .groups = "drop"
  )

# State-level MID counts: number of MID involvements per state-year
# mid_count: all ongoing; fatal_mid_count: fatal only; mid_onset_count: new MIDs;
# fatal_mid_onset_count: new MIDs that are also fatal
state_mid_count <- dyad_data |>
  filter(cowmidongoing == 1) |>
  select(ccode1, ccode2, year, cowmidfatal, cowmidonset) |>
  pivot_longer(cols = c(ccode1, ccode2), values_to = "ccode") |>
  group_by(ccode, year) |>
  summarise(
    mid_count = n(),
    fatal_mid_count = sum(cowmidfatal == 1, na.rm = TRUE),
    mid_onset_count = sum(cowmidonset == 1, na.rm = TRUE),
    fatal_mid_onset_count = sum(cowmidfatal == 1 & cowmidonset == 1, na.rm = TRUE),
    .groups = "drop"
  )

state_churn <- state_churn |>
  left_join(state_war |> select(ccode, year, war_participation, fatal_mid_participation),
            by = c("ccode", "year")) |>
  left_join(state_mid_count, by = c("ccode", "year")) |>
  mutate(
    war_participation = replace_na(war_participation, 0),
    fatal_mid_participation = replace_na(fatal_mid_participation, 0),
    mid_count = replace_na(mid_count, 0L),
    fatal_mid_count = replace_na(fatal_mid_count, 0L),
    mid_onset_count = replace_na(mid_onset_count, 0L),
    fatal_mid_onset_count = replace_na(fatal_mid_onset_count, 0L)
  )

# --------------------------------------------------------------------------- #
# 6. Save outputs
# --------------------------------------------------------------------------- #

cat("=" |> strrep(70), "\n")
cat("Saving outputs...\n")
cat("=" |> strrep(70), "\n\n")

# --- Main specification ---
saveRDS(system_main_ts, file.path(OUTPUT_DIR, "system_network_measures.rds"))
write_csv(system_main_ts, file.path(OUTPUT_DIR, "system_network_measures.csv"))

saveRDS(state_churn, file.path(OUTPUT_DIR, "state_network_churn.rds"))

# --- Robustness: weighting variants ---
saveRDS(system_robustness_ts, file.path(OUTPUT_DIR, "system_robustness_weighting.rds"))
write_csv(system_robustness_ts, file.path(OUTPUT_DIR, "system_robustness_weighting.csv"))

# --- Appendix: layer-specific ---
saveRDS(layer_specific_ts, file.path(OUTPUT_DIR, "layer_specific_measures.rds"))
write_csv(layer_specific_ts, file.path(OUTPUT_DIR, "layer_specific_measures.csv"))

saveRDS(appendix_avg, file.path(OUTPUT_DIR, "appendix_avg_divergence.rds"))
write_csv(appendix_avg, file.path(OUTPUT_DIR, "appendix_avg_divergence.csv"))

# --- System-level analysis data (main + war counts + controls) ---
system_war_counts <- readRDS(file.path(OUTPUT_DIR, "system_war_counts.rds"))

# System-level controls
cat("Computing system-level controls (state count, democracy count)...\n")

system_controls <- dyad_data |>
  group_by(year) |>
  summarise(n_states = n_distinct(c(ccode1, ccode2)), .groups = "drop")

polity_added <- tryCatch({
  library(peacesciencer)
  state_years <- create_stateyears(subset_years = min(dyad_data$year):max(dyad_data$year))
  attr(state_years, "ps_data_type") <- "state_year"
  attr(state_years, "ps_system") <- "cow"
  state_years <- state_years |> add_democracy()

  if ("polity2" %in% names(state_years)) {
    dem_counts <- state_years |>
      group_by(year) |>
      summarise(
        n_democracies = sum(polity2 >= 6, na.rm = TRUE),
        n_autocracies = sum(polity2 <= -6, na.rm = TRUE),
        .groups = "drop"
      )
    system_controls <- system_controls |> left_join(dem_counts, by = "year")
    TRUE
  } else FALSE
}, error = function(e) {
  cat(sprintf("  Warning: Could not load Polity data: %s\n", e$message))
  FALSE
})

if (!polity_added) {
  system_controls$n_democracies <- NA_integer_
  system_controls$n_autocracies <- NA_integer_
}

system_analysis <- system_main_ts |>
  left_join(system_war_counts, by = "year") |>
  left_join(system_controls, by = "year")

saveRDS(system_analysis, file.path(OUTPUT_DIR, "system_analysis_data.rds"))

# --- Metadata update ---
metadata$main_layers <- MAIN_LAYERS
metadata$layer_specific <- TRUE
metadata$note <- "Main spec uses alliance+IGO only (fixed layer set). Appendix has layer-specific measures."
saveRDS(metadata, file.path(OUTPUT_DIR, "network_metadata.rds"))

cat("\nDone.\n")
cat(sprintf("  Main system-level: %d years\n", nrow(system_main_ts)))
cat(sprintf("  Main state-level: %s state-years\n", format(nrow(state_churn), big.mark = ",")))
cat(sprintf("  Robustness weighting variants: %d variant-year obs (%d variants)\n",
            nrow(system_robustness_ts), n_distinct(system_robustness_ts$weight_variant)))
cat(sprintf("  Appendix layer-specific: %d layer-year obs\n", nrow(layer_specific_ts)))
cat(sprintf("  Appendix avg divergence: %d years\n", nrow(appendix_avg)))
