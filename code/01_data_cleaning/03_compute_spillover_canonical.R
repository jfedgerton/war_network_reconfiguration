###############################################################################
# 03_compute_spillover_canonical.R
# One-time precompute of state-year spillover variables on the canonical
# alliance + IGO main_weight network. Output is consumed by:
#   - 02_analysis/03_robustness.R (left-joined into state_data so spillover
#     specs run in every multiverse cell)
#   - 02_analysis/21_spillover_nonfighters.R (for the focused tests)
#
# Two variables produced per state-year:
#   inv_dist_to_war              = 1 / (1 + network_distance_to_nearest_other_warring_state)
#                                  (self-distance excluded; 0 = no path or
#                                  no other warring state reachable;
#                                  1 = direct neighbor of a warring state)
#   partner_war_frac_weighted_2hop = sum(A^2 %*% war_indicator) / sum(rowSums(A^2 - diag))
#                                    (weighted share of partners-of-partners at war)
#
# Timing: predetermined network — t-1 partners, t war indicator. Matches
# the timing of partner_war_frac_weighted in 02_state_level_analysis.R so
# 1-hop, 2-hop, and distance variables are all on the same clock.
#
# Vectorization: uses Matrix::sparseMatrix instead of an inner edge loop.
# Speedup ~50-100x relative to naive for-loop construction.
#
# Drop-in: from project root,
#   Rscript code/01_data_cleaning/03_compute_spillover_canonical.R
#
# Output:
#   data/processed/spillover_canonical.rds
###############################################################################

suppressPackageStartupMessages({
  library(tidyverse)
  library(igraph)
  library(Matrix)
})
select <- dplyr::select

set.seed(123)

OUTPUT_DIR <- here::here("data", "processed")

cat("Loading dyad-year network layers...\n")
dyad_data <- readRDS(file.path(OUTPUT_DIR, "dyad_year_network_layers.rds"))

# Build canonical main_weight if not already present (alliance + IGO).
if (!"main_weight" %in% names(dyad_data)) {
  cat("Building canonical main_weight (offensive + defensive + IGO)...\n")
  ML <- c("layer_offensive", "layer_defensive", "layer_igo")
  mat <- as.matrix(dyad_data[, ML, drop = FALSE])
  mat[is.na(mat)] <- 0
  dyad_data$main_raw <- rowSums(mat)
  dyad_data <- dyad_data |>
    group_by(year) |>
    mutate(main_weight = if (max(main_raw, na.rm = TRUE) > 0) {
      main_raw / max(main_raw, na.rm = TRUE)
    } else 0) |>
    ungroup()
  dyad_data$main_raw <- NULL
}

# State-year war indicator: 1 if state has any ongoing MID this year.
at_war_state <- dyad_data |>
  filter(cowmidongoing == 1) |>
  select(ccode1, ccode2, year) |>
  pivot_longer(c(ccode1, ccode2), values_to = "ccode") |>
  distinct(ccode, year) |>
  mutate(at_war = 1L)

years_loop <- sort(unique(dyad_data$year))
cat(sprintf("Computing spillover variables for %d years (%d-%d)...\n",
              length(years_loop), min(years_loop), max(years_loop)))

# ---------------------------------------------------------------------------- #
# Per-year computation: distance to nearest warring state + 2-hop partner war.
# Vectorized sparse-matrix construction; single distances() call per year.
# ---------------------------------------------------------------------------- #

dist_results   <- list()
twohop_results <- list()

t0 <- Sys.time()

for (yr in years_loop) {
  # Predetermined network: use t-1 edges with year-t war indicator. This
  # matches the timing of `partner_war_frac_weighted` (1-hop spillover
  # built in 02_state_level_analysis.R), so 1-hop, 2-hop, and distance
  # are all on the same (t-1 partners, t war) clock. The earliest year
  # in years_loop yields no edges and is skipped via the `next`.
  edges_yr <- dyad_data |>
    filter(year == yr - 1L, main_weight > 0) |>
    select(ccode1, ccode2, w = main_weight)
  if (nrow(edges_yr) == 0) next

  ccodes   <- sort(unique(c(edges_yr$ccode1, edges_yr$ccode2)))
  n_ccodes <- length(ccodes)
  idx      <- setNames(seq_along(ccodes), as.character(ccodes))

  i_idx <- idx[as.character(edges_yr$ccode1)]
  j_idx <- idx[as.character(edges_yr$ccode2)]

  # ---- (1) Symmetric weighted adjacency via vectorized sparseMatrix ----
  # Each undirected edge contributes one entry in (i, j) and one in (j, i).
  A <- Matrix::sparseMatrix(
    i    = c(i_idx, j_idx),
    j    = c(j_idx, i_idx),
    x    = c(edges_yr$w, edges_yr$w),
    dims = c(n_ccodes, n_ccodes)
  )

  # ---- (2) 2-hop reachability with weights = sum of products of edge weights ----
  # Two variants:
  #   partner_war_frac_weighted_2hop:        all length-2 paths (loose)
  #   partner_war_frac_weighted_2hop_strict: zero out (i,j) cells where i and
  #     j are also DIRECT partners (i.e., A[i,j] > 0). This isolates exposure
  #     to "partners of partners I'm not directly connected to" rather than
  #     direct partners reached via triangle paths. The loose version
  #     conflates these two; the strict version cleanly separates them.
  A2 <- A %*% A
  diag(A2) <- 0
  A2_strict <- A2 * (A == 0)   # zero where direct edge exists
  diag(A2_strict) <- 0          # safety: also zero self

  war_yr  <- at_war_state |> filter(year == yr) |> pull(ccode) |> as.integer()
  war_vec <- as.integer(ccodes %in% war_yr)

  num        <- as.numeric(A2 %*% war_vec)
  den        <- as.numeric(Matrix::rowSums(A2))
  pwf2       <- ifelse(den > 0, num / den, 0)

  num_strict <- as.numeric(A2_strict %*% war_vec)
  den_strict <- as.numeric(Matrix::rowSums(A2_strict))
  pwf2_strict <- ifelse(den_strict > 0, num_strict / den_strict, 0)

  twohop_results[[length(twohop_results) + 1]] <- tibble(
    ccode = ccodes,
    year  = yr,
    partner_war_frac_weighted_2hop        = pwf2,
    partner_war_frac_weighted_2hop_strict = pwf2_strict
  )

  # ---- (3) Network distance to nearest warring state ----
  # igraph BFS on the unweighted version of the graph (presence/absence of tie).
  g <- graph_from_data_frame(
    d = edges_yr |> select(ccode1, ccode2),
    directed = FALSE,
    vertices = data.frame(name = ccodes)
  )
  V(g)$name <- as.character(V(g)$name)

  warring_in_g <- intersect(as.character(war_yr), V(g)$name)

  if (length(warring_in_g) == 0) {
    dist_results[[length(dist_results) + 1]] <- tibble(
      ccode = ccodes, year = yr, dist_to_war = NA_real_
    )
    next
  }

  d_mat <- distances(g, v = V(g), to = warring_in_g, mode = "all")
  d_mat[is.infinite(d_mat)] <- NA_real_
  # Exclude self-distance: if focal i is itself in warring_in_g, the column
  # for i contains 0 (self), which mechanically forces min_d_i = 0 and
  # inv_dist_to_war = 1 — producing perfect collinearity with the focal's
  # own war involvement in downstream regressions once TWFE absorbs
  # ccode/year. Set self-cells to NA so each state's distance is to OTHER
  # warring states only.
  v_idx_of_warring <- match(warring_in_g, V(g)$name)
  self_pairs <- cbind(v_idx_of_warring, seq_along(warring_in_g))
  d_mat[self_pairs] <- NA_real_
  min_d <- suppressWarnings(apply(d_mat, 1, min, na.rm = TRUE))
  min_d[is.infinite(min_d)] <- NA_real_

  dist_results[[length(dist_results) + 1]] <- tibble(
    ccode       = as.integer(V(g)$name),
    year        = yr,
    dist_to_war = min_d
  )

  if (yr %% 25 == 0) {
    cat(sprintf("  ...year %d processed\n", yr))
  }
}

elapsed <- difftime(Sys.time(), t0, units = "secs")
cat(sprintf("Per-year computation finished (%.1f sec).\n", as.numeric(elapsed)))

# ---------------------------------------------------------------------------- #
# Combine + transform.
# ---------------------------------------------------------------------------- #

dist_df <- bind_rows(dist_results) |>
  mutate(inv_dist_to_war = ifelse(is.na(dist_to_war), 0,
                                    1 / (1 + dist_to_war)))

twohop_df <- bind_rows(twohop_results) |>
  mutate(
    partner_war_frac_weighted_2hop        = replace_na(partner_war_frac_weighted_2hop, 0),
    partner_war_frac_weighted_2hop_strict = replace_na(partner_war_frac_weighted_2hop_strict, 0)
  )

spillover_canonical <- full_join(dist_df, twohop_df, by = c("ccode", "year")) |>
  mutate(
    dist_to_war                           = dist_to_war,
    inv_dist_to_war                       = replace_na(inv_dist_to_war, 0),
    partner_war_frac_weighted_2hop        = replace_na(partner_war_frac_weighted_2hop, 0),
    partner_war_frac_weighted_2hop_strict = replace_na(partner_war_frac_weighted_2hop_strict, 0)
  ) |>
  arrange(ccode, year)

cat(sprintf("\nspillover_canonical.rds: %s state-year rows\n",
              format(nrow(spillover_canonical), big.mark = ",")))
cat(sprintf("  range inv_dist_to_war:                 [%.3f, %.3f] (mean %.3f)\n",
              min(spillover_canonical$inv_dist_to_war, na.rm = TRUE),
              max(spillover_canonical$inv_dist_to_war, na.rm = TRUE),
              mean(spillover_canonical$inv_dist_to_war, na.rm = TRUE)))
cat(sprintf("  range partner_war_frac_weighted_2hop:  [%.3f, %.3f] (mean %.3f)\n",
              min(spillover_canonical$partner_war_frac_weighted_2hop, na.rm = TRUE),
              max(spillover_canonical$partner_war_frac_weighted_2hop, na.rm = TRUE),
              mean(spillover_canonical$partner_war_frac_weighted_2hop, na.rm = TRUE)))

saveRDS(spillover_canonical,
         file.path(OUTPUT_DIR, "spillover_canonical.rds"))
cat(sprintf("\nSaved: %s\n",
              file.path("data", "processed", "spillover_canonical.rds")))
cat("Done.\n")
