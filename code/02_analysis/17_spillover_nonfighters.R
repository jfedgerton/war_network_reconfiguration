###############################################################################
# 21_spillover_nonfighters.R
# Tests whether war affects the network configuration of CONNECTED but
# NON-BELLIGERENT states — five complementary tests motivated by the
# observation that the main-spec partner_war_frac_weighted control is
# largely null. Null on the coarse measure may reflect the measure's
# weakness rather than absence of spillovers.
#
# The five tests (all 1-hop and 2-hop variants are weighted, on the
# canonical alliance + IGO main_weight network):
#   1. Non-fighter subsample (fatal_mid_count == 0; threats / verbal MIDs
#      remain in the pool, only use-of-force/battle-death cases drop):
#      1a. System-level IV (system_mid_count, system_fatal_mid_count) with
#          state FE + poly(year, 2) — no year FE because IV is year-level.
#      1b. Partner-war IV with state + year FE — heterogeneity within
#          non-fighter year cells.
#   2. Network-geodesic distance to nearest warring state — converts the
#      partner-share measure into a topological-proximity measure.
#   3. Higher-order (2-hop) partner war exposure — partners-of-partners:
#      3.  loose:   length-2 paths via A^2 (conflates direct + indirect
#                   neighbors reached through triangles).
#      3b. strict:  A^2 * (A == 0) — zeroes cells where a direct edge
#                   exists, isolating partners-of-partners NOT directly
#                   connected to the focal state.
#   4. Lead-lag analysis on partner_war_frac_weighted (t-2 ... t+2),
#      proper-timing construction: NETWORK held fixed at t-1, WAR YEAR
#      varied. Leads are placebos; lags test for delayed reshaping.
#      Avoids the built-in-lag confound of `lag(partner_war_frac_weighted, k)`,
#      which would shift both the network and the war year together.
#   5. Heterogeneous effects on partner_war_frac_weighted (the only
#      spillover IV with a positive significant main effect):
#      a. x mid_onset_count        (does focal MID involvement amplify?)
#      b. x fatal_mid_onset_count  (use-of-force version of a)
#      c. x major_power            (great-power moderation)
#      d. x democracy              (regime-type moderation)
#      e. x major_war_year         (1914-18, 1939-45, 1950-53, 1990-91, 2003)
#      Each runs on the full panel and the non-fighter subsample.
#
# Outputs:
#   data/processed/spillover_nonfighters_results.csv
#   results/appendix/tableA_spillover_nonfighters.tex
#   code/03_figures/output/manuscript/fig_spillover_nonfighters.png
#
# Convention: each model object is keyed `s24__<test>__<iv>` so downstream
# tables can filter cleanly.
###############################################################################

library(tidyverse)
library(fixest)
library(igraph)
library(peacesciencer)
library(Matrix)
select <- dplyr::select

set.seed(123)

OUTPUT_DIR <- here::here("data", "processed")
APPX_DIR   <- here::here("results", "appendix")
FIG_DIR    <- here::here("code", "03_figures", "output", "manuscript")
dir.create(APPX_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG_DIR,  recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------
# Manuscript theme
# ------------------------------------------------------------
ACCENT  <- "#2C6CB0"
ACCENT2 <- "#B03366"
ACCENT3 <- "#5A8F3D"

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

# ============================================================
# Load data
# ============================================================

cat("Loading data...\n")

state_churn <- readRDS(file.path(OUTPUT_DIR, "state_network_churn.rds"))
system_data <- readRDS(file.path(OUTPUT_DIR, "system_analysis_data.rds")) |>
  select(year, system_mid_count = mid_count,
         system_fatal_mid_count = fatal_mid_count,
         n_democracies, n_autocracies, n_states)
dyad_data <- readRDS(file.path(OUTPUT_DIR, "dyad_year_network_layers.rds"))

# Compute main_weight if missing.
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

# State attributes (cinc, major_power)
state_attrs <- dyad_data |>
  select(ccode = ccode1, year, cinc = cinc1, major_power = cowmaj1) |>
  distinct(ccode, year, .keep_all = TRUE) |>
  mutate(major_power = as.integer(coalesce(major_power, 0L) == 1))

# Polity (democracy / autocracy)
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
  left_join(state_attrs, by = c("ccode", "year")) |>
  left_join(polity,      by = c("ccode", "year")) |>
  left_join(system_data, by = "year")

# ============================================================
# Inline 1-hop partner-war exposure (standard harmonized block)
# ============================================================

partner_sets <- dyad_data |>
  filter(main_weight > 0) |>
  select(ccode1, ccode2, year, w = main_weight) |>
  bind_rows(dyad_data |> filter(main_weight > 0) |>
              select(ccode1 = ccode2, ccode2 = ccode1, year, w = main_weight))

at_war_dyad <- dyad_data |>
  filter(cowmidongoing == 1) |>
  select(ccode1, ccode2, year) |>
  pivot_longer(c(ccode1, ccode2), values_to = "ccode") |>
  distinct(ccode, year) |>
  mutate(at_war = 1L)

at_war_dyad_fatal <- dyad_data |>
  filter(cowmidfatal == 1) |>
  select(ccode1, ccode2, year) |>
  pivot_longer(c(ccode1, ccode2), values_to = "ccode") |>
  distinct(ccode, year) |>
  mutate(at_war_fatal = 1L)

pwe <- partner_sets |>
  mutate(year_next = year + 1L) |>
  left_join(at_war_dyad,
            by = c("ccode2" = "ccode", "year_next" = "year")) |>
  left_join(at_war_dyad_fatal,
            by = c("ccode2" = "ccode", "year_next" = "year")) |>
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
         major_power                     = replace_na(major_power, 0))

# --- Partner-war control routing by IV severity ---
# Non-fatal IVs (any MID exposure) -> partner_war_frac
# Fatal IVs (battle-death-bearing MIDs) -> partner_fatal_war_frac
# Weighted variants are only used as explicit alternative specifications
# (Tests 1b / 4 / 5 use partner_war_frac_weighted as the IV being tested).
CTRL_MID   <- "partner_war_frac + democracy + autocracy + cinc + major_power + is_isolated_t1"
CTRL_FATAL <- "partner_fatal_war_frac + democracy + autocracy + cinc + major_power + is_isolated_t1"
ctrl_for_iv <- function(iv) if (grepl("fatal", iv)) CTRL_FATAL else CTRL_MID

# ============================================================
# Test 2 PREP: network-geodesic distance to nearest warring state
# ============================================================

cat("Computing network distance to nearest warring state per year...\n")

at_war_state <- dyad_data |>
  filter(cowmidongoing == 1) |>
  select(ccode1, ccode2, year) |>
  pivot_longer(c(ccode1, ccode2), values_to = "ccode") |>
  distinct(ccode, year)

# Use precomputed canonical spillover variables if available
# (01_data_cleaning/03_compute_spillover_canonical.R produces these via vectorized sparseMatrix
# construction; ~50x faster than inline). Fall back to inline vectorized
# computation otherwise — also fast, but only run if the RDS is missing.

spillover_path <- file.path(OUTPUT_DIR, "spillover_canonical.rds")
years_loop <- sort(unique(dyad_data$year))

if (file.exists(spillover_path)) {
  cat("Loading precomputed spillover variables from spillover_canonical.rds...\n")
  spillover_df <- readRDS(spillover_path)
  state_churn <- state_churn |>
    left_join(spillover_df |>
                select(ccode, year,
                        dist_to_war, inv_dist_to_war,
                        partner_war_frac_weighted_2hop,
                        any_of("partner_war_frac_weighted_2hop_strict")),
              by = c("ccode", "year")) |>
    mutate(
      inv_dist_to_war                = replace_na(inv_dist_to_war, 0),
      partner_war_frac_weighted_2hop = replace_na(partner_war_frac_weighted_2hop, 0)
    )
  if ("partner_war_frac_weighted_2hop_strict" %in% names(state_churn)) {
    state_churn <- state_churn |>
      mutate(partner_war_frac_weighted_2hop_strict =
               replace_na(partner_war_frac_weighted_2hop_strict, 0))
  }
} else {
  cat("spillover_canonical.rds not found — computing inline (vectorized).\n")

  dist_results   <- list()
  twohop_results <- list()

  for (yr in years_loop) {
    # Predetermined network: use t-1 edges with year-t war indicator
    # (matches partner_war_frac_weighted timing in 02). Mirrors the
    # canonical version in 01_data_cleaning/03_compute_spillover_canonical.R.
    edges_yr <- dyad_data |>
      filter(year == yr - 1L, main_weight > 0) |>
      select(ccode1, ccode2, w = main_weight)
    if (nrow(edges_yr) == 0) next

    ccodes   <- sort(unique(c(edges_yr$ccode1, edges_yr$ccode2)))
    n_ccodes <- length(ccodes)
    idx      <- setNames(seq_along(ccodes), as.character(ccodes))
    i_idx    <- idx[as.character(edges_yr$ccode1)]
    j_idx    <- idx[as.character(edges_yr$ccode2)]

    # Vectorized symmetric weighted adjacency (single sparseMatrix call).
    A <- Matrix::sparseMatrix(
      i = c(i_idx, j_idx),
      j = c(j_idx, i_idx),
      x = c(edges_yr$w, edges_yr$w),
      dims = c(n_ccodes, n_ccodes)
    )

    # ---- 2-hop (loose + strict) ----
    A2 <- A %*% A
    diag(A2) <- 0
    A2_strict <- A2 * (A == 0)   # zero out cells where direct edge exists
    diag(A2_strict) <- 0

    war_yr  <- at_war_state |> filter(year == yr) |> pull(ccode) |> as.integer()
    war_vec <- as.integer(ccodes %in% war_yr)

    num <- as.numeric(A2 %*% war_vec)
    den <- as.numeric(Matrix::rowSums(A2))
    pwf2 <- ifelse(den > 0, num / den, 0)

    num_strict <- as.numeric(A2_strict %*% war_vec)
    den_strict <- as.numeric(Matrix::rowSums(A2_strict))
    pwf2_strict <- ifelse(den_strict > 0, num_strict / den_strict, 0)

    twohop_results[[length(twohop_results) + 1]] <- tibble(
      ccode = ccodes, year = yr,
      partner_war_frac_weighted_2hop        = pwf2,
      partner_war_frac_weighted_2hop_strict = pwf2_strict
    )

    # ---- distance ----
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
    # Exclude self-distance: warring focal states would otherwise have
    # min_d = 0 (self) and inv_dist_to_war = 1 mechanically, contaminating
    # the second-order spillover regression with the focal's own war.
    v_idx_of_warring <- match(warring_in_g, V(g)$name)
    self_pairs <- cbind(v_idx_of_warring, seq_along(warring_in_g))
    d_mat[self_pairs] <- NA_real_
    min_d <- suppressWarnings(apply(d_mat, 1, min, na.rm = TRUE))
    min_d[is.infinite(min_d)] <- NA_real_

    dist_results[[length(dist_results) + 1]] <- tibble(
      ccode = as.integer(V(g)$name), year = yr, dist_to_war = min_d
    )
  }

  dist_df <- bind_rows(dist_results) |>
    mutate(inv_dist_to_war = ifelse(is.na(dist_to_war), 0,
                                      1 / (1 + dist_to_war)))
  twohop_df <- bind_rows(twohop_results)

  state_churn <- state_churn |>
    left_join(dist_df,   by = c("ccode", "year")) |>
    left_join(twohop_df, by = c("ccode", "year")) |>
    mutate(
      inv_dist_to_war                       = replace_na(inv_dist_to_war, 0),
      partner_war_frac_weighted_2hop        = replace_na(partner_war_frac_weighted_2hop, 0),
      partner_war_frac_weighted_2hop_strict = replace_na(partner_war_frac_weighted_2hop_strict, 0)
    )
}

# ============================================================
# Final harmonization steps (DV + IV transforms)
# ============================================================

state_churn <- state_churn |>
  mutate(
    log_churn                     = log1p(churn),
    log_mid_onset_count           = log1p(mid_onset_count),
    log_fatal_mid_onset_count     = log1p(fatal_mid_onset_count),
    log_system_mid_count          = log1p(system_mid_count),
    log_system_fatal_mid_count    = log1p(system_fatal_mid_count),
    is_fighter                    = as.integer(mid_count > 0),
    # Major-war years (system-wide great-power conflicts)
    major_war_year = as.integer(
      year %in% c(1914:1918, 1939:1945, 1950:1953, 1990:1991, 2003))
  )

# --- Properly-timed lead-lag of weighted partner-war exposure ---
# Build pwf_war_K for K in {-2, -1, +1, +2} where K is the war-year offset
# relative to the focal year. ALL variants use the SAME network (partners
# at t-1, the predetermined network); only the WAR YEAR varies. This
# avoids the built-in-lag confound of `lag(partner_war_frac_weighted)`,
# which would shift BOTH the network AND the war year together.
#
#   pwf_war_lag2  = (1/sum(w)) sum_j w * 1{partner j at war at t-2}
#   pwf_war_lag1  = ...                              at war at t-1
#   partner_war_frac_weighted (existing) =                at war at t   <-- current
#   pwf_war_lead1 = ...                              at war at t+1   <-- placebo
#   pwf_war_lead2 = ...                              at war at t+2   <-- placebo
#
# Each is built by joining partner_sets (year = t-1) with at_war_dyad
# at year_war = (year + 1) + offset, then aggregating weighted means.
build_pwf_offset <- function(offset_focal, label) {
  partner_sets |>
    mutate(year_focal = year + 1L,
           year_war   = year_focal + offset_focal) |>
    left_join(at_war_dyad,
              by = c("ccode2" = "ccode", "year_war" = "year")) |>
    mutate(at_war = replace_na(at_war, 0L)) |>
    group_by(ccode1, year_focal) |>
    summarise(!!label := sum(at_war * w) / sum(w), .groups = "drop") |>
    rename(ccode = ccode1, year = year_focal)
}
pwf_offsets <- list(
  list(off = -2L, lab = "pwf_war_lag2"),
  list(off = -1L, lab = "pwf_war_lag1"),
  list(off =  1L, lab = "pwf_war_lead1"),
  list(off =  2L, lab = "pwf_war_lead2")
)
for (po in pwf_offsets) {
  state_churn <- state_churn |>
    left_join(build_pwf_offset(po$off, po$lab), by = c("ccode", "year")) |>
    mutate(!!po$lab := replace_na(.data[[po$lab]], 0))
}

# Two non-fighter subsamples. Definition uses FATAL MID rather than ANY
# MID: a state filing a low-intensity threat (HostLev 1-2) shouldn't be
# excluded from the "non-fighter" sample because the spillover question
# is whether STATES NOT EXPERIENCING USE-OF-FORCE/WAR have their networks
# reshaped by other states' wars. Using `fatal_mid_count == 0` recovers
# ~30% more sample size (state-years with only verbal protests / threats
# come back into the spillover pool) without compromising the
# "non-belligerent" interpretation.
state_churn <- state_churn |>
  group_by(ccode) |>
  arrange(year, .by_group = TRUE) |>
  mutate(prev_fatal_mid_count = dplyr::lag(fatal_mid_count, 1)) |>
  ungroup()

state_nonfighter        <- state_churn |> filter(fatal_mid_count == 0)
state_nonfighter_strict <- state_churn |> filter(fatal_mid_count == 0,
                                                  prev_fatal_mid_count == 0)

cat(sprintf("  Non-fighter sample: %s state-years (%.1f%% of full)\n",
            format(nrow(state_nonfighter), big.mark = ","),
            100 * nrow(state_nonfighter) / nrow(state_churn)))
cat(sprintf("  Strict non-fighter sample: %s state-years (%.1f%% of full)\n",
            format(nrow(state_nonfighter_strict), big.mark = ","),
            100 * nrow(state_nonfighter_strict) / nrow(state_churn)))

# ============================================================
# Common controls (state level)
# ============================================================

ctrl_full_mid    <- CTRL_MID    # partner_war_frac     + democracy + ... + is_isolated_t1
ctrl_full_fatal  <- CTRL_FATAL  # partner_fatal_war_frac + democracy + ... + is_isolated_t1
ctrl_full_no_pwf <- "democracy + autocracy + cinc + major_power + is_isolated_t1"

# ============================================================
# Reference: harmonized main spec on full sample (anchor row)
# ============================================================

m_main <- feols(
  as.formula(paste0("log_churn ~ mid_onset_count + ",
                      ctrl_full_mid, " | ccode + year")),
  data = state_churn, cluster = ~ ccode + year)
m_main_fatal <- feols(
  as.formula(paste0("log_churn ~ fatal_mid_onset_count + ",
                      ctrl_full_fatal, " | ccode + year")),
  data = state_churn, cluster = ~ ccode + year)

# ============================================================
# TEST 1a — system-level IV on non-fighters
#   State FE only (no year FE because IV is year-level)
#   poly(year, 2) absorbs secular trend
# ============================================================

cat("\n=== TEST 1a: system-level IV on non-fighter subsample ===\n")

m_1a_mid <- feols(
  log_churn ~ system_mid_count + partner_war_frac +
    democracy + autocracy + cinc + major_power + is_isolated_t1 + poly(year, 2) | ccode,
  data    = state_nonfighter, cluster = ~ ccode + year)
m_1a_fatal <- feols(
  log_churn ~ system_fatal_mid_count + partner_fatal_war_frac +
    democracy + autocracy + cinc + major_power + is_isolated_t1 + poly(year, 2) | ccode,
  data    = state_nonfighter, cluster = ~ ccode + year)
m_1a_log_mid <- feols(
  log_churn ~ log_system_mid_count + partner_war_frac +
    democracy + autocracy + cinc + major_power + is_isolated_t1 + poly(year, 2) | ccode,
  data    = state_nonfighter, cluster = ~ ccode + year)
m_1a_log_fatal <- feols(
  log_churn ~ log_system_fatal_mid_count + partner_fatal_war_frac +
    democracy + autocracy + cinc + major_power + is_isolated_t1 + poly(year, 2) | ccode,
  data    = state_nonfighter, cluster = ~ ccode + year)

# Strict non-fighters (also peaceful at t-1)
m_1a_mid_strict <- feols(
  log_churn ~ system_mid_count + partner_war_frac +
    democracy + autocracy + cinc + major_power + is_isolated_t1 + poly(year, 2) | ccode,
  data    = state_nonfighter_strict, cluster = ~ ccode + year)
m_1a_fatal_strict <- feols(
  log_churn ~ system_fatal_mid_count + partner_fatal_war_frac +
    democracy + autocracy + cinc + major_power + is_isolated_t1 + poly(year, 2) | ccode,
  data    = state_nonfighter_strict, cluster = ~ ccode + year)

# ============================================================
# TEST 1b — partner-war IV on non-fighters with state + year FE
# ============================================================

cat("=== TEST 1b: partner-war IV on non-fighter subsample ===\n")

m_1b <- feols(
  as.formula(paste0("log_churn ~ partner_war_frac_weighted + ",
                      ctrl_full_no_pwf, " | ccode + year")),
  data    = state_nonfighter, cluster = ~ ccode + year)
m_1b_strict <- feols(
  as.formula(paste0("log_churn ~ partner_war_frac_weighted + ",
                      ctrl_full_no_pwf, " | ccode + year")),
  data    = state_nonfighter_strict, cluster = ~ ccode + year)

# ============================================================
# TEST 2 — network-geodesic distance to nearest warring state
#   Pair with partner_war_frac_weighted (1-hop weighted) as control to
#   keep the spillover ladder on a consistent weighted network.
# ============================================================

cat("=== TEST 2: network-distance-to-war ===\n")

m_2_full <- feols(
  log_churn ~ inv_dist_to_war + partner_war_frac_weighted +
    democracy + autocracy + cinc + major_power + is_isolated_t1 | ccode + year,
  data    = state_churn, cluster = ~ ccode + year)
m_2_nonfighter <- feols(
  log_churn ~ inv_dist_to_war + partner_war_frac_weighted +
    democracy + autocracy + cinc + major_power + is_isolated_t1 | ccode + year,
  data    = state_nonfighter, cluster = ~ ccode + year)

# ============================================================
# TEST 3 — 2-hop partner war exposure
#   Two variants:
#     3 (loose):  partner_war_frac_weighted_2hop  (length-2 paths; conflates
#                                                   direct + indirect partners)
#     3b (strict): partner_war_frac_weighted_2hop_strict (zeroes cells where
#                  a direct edge exists; isolates "partners-of-partners I
#                  am NOT directly connected to")
#   Both pair with partner_war_frac_weighted (1-hop weighted) so the IV /
#   control are on the same weighted network.
# ============================================================

cat("=== TEST 3: 2-hop partner war exposure (loose) ===\n")

m_3_full <- feols(
  log_churn ~ partner_war_frac_weighted_2hop + partner_war_frac_weighted +
    democracy + autocracy + cinc + major_power + is_isolated_t1 | ccode + year,
  data    = state_churn, cluster = ~ ccode + year)
m_3_nonfighter <- feols(
  log_churn ~ partner_war_frac_weighted_2hop + partner_war_frac_weighted +
    democracy + autocracy + cinc + major_power + is_isolated_t1 | ccode + year,
  data    = state_nonfighter, cluster = ~ ccode + year)

cat("=== TEST 3b: 2-hop partner war exposure (strict, no direct edges) ===\n")

m_3b_full <- feols(
  log_churn ~ partner_war_frac_weighted_2hop_strict + partner_war_frac_weighted +
    democracy + autocracy + cinc + major_power + is_isolated_t1 | ccode + year,
  data    = state_churn, cluster = ~ ccode + year)
m_3b_nonfighter <- feols(
  log_churn ~ partner_war_frac_weighted_2hop_strict + partner_war_frac_weighted +
    democracy + autocracy + cinc + major_power + is_isolated_t1 | ccode + year,
  data    = state_nonfighter, cluster = ~ ccode + year)

# ============================================================
# TEST 4 — lead-lag of partner_war_frac_weighted
#   Uses the properly-timed pwf_war_{lag2, lag1, lead1, lead2}: network
#   held fixed at t-1, war year varied. Avoids the built-in-lag confound
#   of `lag(partner_war_frac_weighted, k)`, which would shift BOTH the
#   network AND the war year together.
# ============================================================

cat("=== TEST 4: lead-lag dynamics (proper-timing) ===\n")

m_4 <- feols(
  log_churn ~ pwf_war_lag2 + pwf_war_lag1 + partner_war_frac_weighted +
    pwf_war_lead1 + pwf_war_lead2 +
    democracy + autocracy + cinc + major_power + is_isolated_t1 | ccode + year,
  data    = state_churn, cluster = ~ ccode + year)
m_4_nonfighter <- feols(
  log_churn ~ pwf_war_lag2 + pwf_war_lag1 + partner_war_frac_weighted +
    pwf_war_lead1 + pwf_war_lead2 +
    democracy + autocracy + cinc + major_power + is_isolated_t1 | ccode + year,
  data    = state_nonfighter, cluster = ~ ccode + year)

# ============================================================
# TEST 5 — heterogeneous effects (interactions)
#   Five interactions on partner_war_frac_weighted (the spillover IV that
#   is positive and significant in the harmonized main spec):
#     a. PWF x mid_onset_count          (does focal MID involvement
#                                         amplify spillover sensitivity?)
#     b. PWF x fatal_mid_onset_count    (same, for use-of-force MIDs)
#     c. PWF x major_power              (great-power moderation)
#     d. PWF x democracy                (regime-type moderation)
#     e. PWF x major_war_year           (system-shock years amplification)
#   Each runs on the full panel and the non-fighter subsample (a, b, c, d
#   only — interaction with onset is mechanically zero in non-fighter
#   sample, so we still fit it to confirm null and document attenuation).
# ============================================================

cat("=== TEST 5: heterogeneous effects ===\n")

# --- a. PWF x onset (any-MID) ---
m_5_onset_mid <- feols(
  log_churn ~ partner_war_frac_weighted * mid_onset_count +
    democracy + autocracy + cinc + major_power + is_isolated_t1 | ccode + year,
  data    = state_churn, cluster = ~ ccode + year)
m_5_onset_mid_nf <- feols(
  log_churn ~ partner_war_frac_weighted * mid_onset_count +
    democracy + autocracy + cinc + major_power + is_isolated_t1 | ccode + year,
  data    = state_nonfighter, cluster = ~ ccode + year)

# --- b. PWF x fatal onset ---
m_5_onset_fatal <- feols(
  log_churn ~ partner_war_frac_weighted * fatal_mid_onset_count +
    democracy + autocracy + cinc + major_power + is_isolated_t1 | ccode + year,
  data    = state_churn, cluster = ~ ccode + year)
m_5_onset_fatal_nf <- feols(
  log_churn ~ partner_war_frac_weighted * fatal_mid_onset_count +
    democracy + autocracy + cinc + major_power + is_isolated_t1 | ccode + year,
  data    = state_nonfighter, cluster = ~ ccode + year)

# --- c. PWF x major power ---
m_5_majorpower <- feols(
  log_churn ~ partner_war_frac_weighted * major_power +
    democracy + autocracy + cinc + is_isolated_t1 | ccode + year,
  data    = state_churn, cluster = ~ ccode + year)
m_5_majorpower_nf <- feols(
  log_churn ~ partner_war_frac_weighted * major_power +
    democracy + autocracy + cinc + is_isolated_t1 | ccode + year,
  data    = state_nonfighter, cluster = ~ ccode + year)

# --- d. PWF x democracy ---
m_5_democracy <- feols(
  log_churn ~ partner_war_frac_weighted * democracy +
    autocracy + cinc + major_power + is_isolated_t1 | ccode + year,
  data    = state_churn, cluster = ~ ccode + year)
m_5_democracy_nf <- feols(
  log_churn ~ partner_war_frac_weighted * democracy +
    autocracy + cinc + major_power + is_isolated_t1 | ccode + year,
  data    = state_nonfighter, cluster = ~ ccode + year)

# --- e. PWF x major-war year ---
m_5_majorwar <- feols(
  log_churn ~ partner_war_frac_weighted * major_war_year +
    democracy + autocracy + cinc + major_power + is_isolated_t1 | ccode + year,
  data    = state_churn, cluster = ~ ccode + year)
m_5_majorwar_nf <- feols(
  log_churn ~ partner_war_frac_weighted * major_war_year +
    democracy + autocracy + cinc + major_power + is_isolated_t1 | ccode + year,
  data    = state_nonfighter, cluster = ~ ccode + year)

# ============================================================
# Pull a single coefficient out of each fit
# ============================================================

extract_one <- function(fit, term, test, label, sample) {
  ct <- coeftable(fit)
  if (!term %in% rownames(ct)) {
    return(tibble(test = test, label = label, sample = sample,
                  term = term, estimate = NA_real_, se = NA_real_,
                  pval = NA_real_, n_obs = nobs(fit)))
  }
  tibble(test = test, label = label, sample = sample,
         term      = term,
         estimate  = ct[term, "Estimate"],
         se        = ct[term, "Std. Error"],
         pval      = ct[term, "Pr(>|t|)"],
         n_obs     = nobs(fit))
}

results_rows <- bind_rows(
  # Anchor — main spec on full sample
  extract_one(m_main,           "mid_onset_count",
              "main_anchor",    "MID onset (full)",          "all"),
  extract_one(m_main_fatal,     "fatal_mid_onset_count",
              "main_anchor",    "Fatal MID onset (full)",    "all"),

  # Test 1a — system-level IV
  extract_one(m_1a_mid,          "system_mid_count",
              "1a_system_iv",   "System MID count",          "nonfighter"),
  extract_one(m_1a_fatal,        "system_fatal_mid_count",
              "1a_system_iv",   "System fatal MID count",    "nonfighter"),
  extract_one(m_1a_log_mid,      "log_system_mid_count",
              "1a_system_iv",   "log(1+System MID)",         "nonfighter"),
  extract_one(m_1a_log_fatal,    "log_system_fatal_mid_count",
              "1a_system_iv",   "log(1+Sys. fatal)",         "nonfighter"),
  extract_one(m_1a_mid_strict,   "system_mid_count",
              "1a_system_iv",   "System MID count",          "nonfighter_strict"),
  extract_one(m_1a_fatal_strict, "system_fatal_mid_count",
              "1a_system_iv",   "System fatal MID count",    "nonfighter_strict"),

  # Test 1b — partner-war IV
  extract_one(m_1b,         "partner_war_frac_weighted",
              "1b_partner_iv", "Partner war (weighted)",     "nonfighter"),
  extract_one(m_1b_strict,  "partner_war_frac_weighted",
              "1b_partner_iv", "Partner war (weighted)",     "nonfighter_strict"),

  # Test 2 — network distance
  extract_one(m_2_full,        "inv_dist_to_war",
              "2_distance",    "1/(1+dist to war)",          "all"),
  extract_one(m_2_nonfighter,  "inv_dist_to_war",
              "2_distance",    "1/(1+dist to war)",          "nonfighter"),

  # Test 3 — 2-hop partner war (loose)
  extract_one(m_3_full,        "partner_war_frac_weighted_2hop",
              "3_twohop",      "Partner war (2-hop, wt.)",   "all"),
  extract_one(m_3_nonfighter,  "partner_war_frac_weighted_2hop",
              "3_twohop",      "Partner war (2-hop, wt.)",   "nonfighter"),

  # Test 3b — 2-hop partner war (strict; excludes direct partners)
  extract_one(m_3b_full,       "partner_war_frac_weighted_2hop_strict",
              "3b_twohop_strict", "Partner war (2-hop strict, wt.)", "all"),
  extract_one(m_3b_nonfighter, "partner_war_frac_weighted_2hop_strict",
              "3b_twohop_strict", "Partner war (2-hop strict, wt.)", "nonfighter"),

  # Test 4 — lead-lag (proper-timing: fixed t-1 network, varying war year)
  extract_one(m_4, "pwf_war_lag2",               "4_leadlag",
              "PWF war t-2",                "all"),
  extract_one(m_4, "pwf_war_lag1",               "4_leadlag",
              "PWF war t-1",                "all"),
  extract_one(m_4, "partner_war_frac_weighted",  "4_leadlag",
              "PWF war t",                  "all"),
  extract_one(m_4, "pwf_war_lead1",              "4_leadlag",
              "PWF war t+1 (placebo)",      "all"),
  extract_one(m_4, "pwf_war_lead2",              "4_leadlag",
              "PWF war t+2 (placebo)",      "all"),
  extract_one(m_4_nonfighter, "pwf_war_lag2",    "4_leadlag",
              "PWF war t-2",                "nonfighter"),
  extract_one(m_4_nonfighter, "pwf_war_lag1",    "4_leadlag",
              "PWF war t-1",                "nonfighter"),
  extract_one(m_4_nonfighter, "partner_war_frac_weighted", "4_leadlag",
              "PWF war t",                  "nonfighter"),
  extract_one(m_4_nonfighter, "pwf_war_lead1",   "4_leadlag",
              "PWF war t+1 (placebo)",      "nonfighter"),
  extract_one(m_4_nonfighter, "pwf_war_lead2",   "4_leadlag",
              "PWF war t+2 (placebo)",      "nonfighter"),

  # Test 5 — heterogeneous
  extract_one(m_5_onset_mid,        "partner_war_frac_weighted:mid_onset_count",
              "5_hetero", "PWF x MID onset",            "all"),
  extract_one(m_5_onset_mid_nf,     "partner_war_frac_weighted:mid_onset_count",
              "5_hetero", "PWF x MID onset",            "nonfighter"),
  extract_one(m_5_onset_fatal,      "partner_war_frac_weighted:fatal_mid_onset_count",
              "5_hetero", "PWF x fatal MID onset",      "all"),
  extract_one(m_5_onset_fatal_nf,   "partner_war_frac_weighted:fatal_mid_onset_count",
              "5_hetero", "PWF x fatal MID onset",      "nonfighter"),
  extract_one(m_5_majorpower,       "partner_war_frac_weighted:major_power",
              "5_hetero", "PWF x major power",          "all"),
  extract_one(m_5_majorpower_nf,    "partner_war_frac_weighted:major_power",
              "5_hetero", "PWF x major power",          "nonfighter"),
  extract_one(m_5_democracy,        "partner_war_frac_weighted:democracy",
              "5_hetero", "PWF x democracy",            "all"),
  extract_one(m_5_democracy_nf,     "partner_war_frac_weighted:democracy",
              "5_hetero", "PWF x democracy",            "nonfighter"),
  extract_one(m_5_majorwar,         "partner_war_frac_weighted:major_war_year",
              "5_hetero", "PWF x major-war year",       "all"),
  extract_one(m_5_majorwar_nf,      "partner_war_frac_weighted:major_war_year",
              "5_hetero", "PWF x major-war year",       "nonfighter")
)

results_rows <- results_rows |>
  mutate(ci_low  = estimate - 1.96 * se,
         ci_high = estimate + 1.96 * se,
         sig05   = !is.na(pval) & pval < 0.05)

print(results_rows, n = Inf, digits = 4)

write_csv(results_rows,
          file.path(OUTPUT_DIR, "spillover_nonfighters_results.csv"))
cat(sprintf("\nSaved: spillover_nonfighters_results.csv (%d rows)\n",
            nrow(results_rows)))

# ============================================================
# LaTeX summary table → results/appendix
# ============================================================

tex <- c(
  "\\begin{table}[htbp]", "\\centering",
  "\\caption{Network spillover tests: does war affect non-belligerents' network configuration?}",
  "\\label{tab:spillover_nonfighters}",
  "\\small",
  "\\begin{tabular}{llrrrr}", "\\toprule",
  "Test & Coefficient & Sample & Estimate & SE & $p$ \\\\",
  "\\midrule"
)
test_lab <- c(
  main_anchor       = "Main anchor (full)",
  `1a_system_iv`    = "1a. System IV",
  `1b_partner_iv`   = "1b. Partner IV",
  `2_distance`      = "2. Network distance",
  `3_twohop`        = "3. 2-hop partner",
  `3b_twohop_strict`= "3b. 2-hop (strict)",
  `4_leadlag`       = "4. Lead-lag",
  `5_hetero`        = "5. Heterogeneous"
)
for (i in seq_len(nrow(results_rows))) {
  r <- results_rows[i, ]
  tex <- c(tex, sprintf(
    "%s & %s & %s & %.4f & %.4f & %.3f \\\\",
    test_lab[r$test],
    gsub("&", "\\&", r$label, fixed = TRUE),
    r$sample,
    ifelse(is.na(r$estimate), NA_real_, r$estimate),
    ifelse(is.na(r$se),       NA_real_, r$se),
    ifelse(is.na(r$pval),     NA_real_, r$pval)
  ))
}
tex <- c(tex, "\\bottomrule", "\\end{tabular}",
  "\\begin{tablenotes}\\small",
  "\\item DV: $\\log(1+\\text{churn})$. State and (where applicable) year FE; SEs two-way clustered (state + year). \"Non-fighter\" subsample = state-years with \\texttt{fatal\\_mid\\_count} = 0 (use-of-force / battle-death MIDs); \"strict\" also requires \\texttt{fatal\\_mid\\_count} = 0 at $t-1$. Test 1a uses state FE + quadratic year trend (year-level IV precludes year FE). Tests 1b, 2-5 use state and year FE. Test 3 = loose 2-hop (length-2 paths); Test 3b = strict 2-hop (zeroes cells where a direct edge already exists, isolating partners-of-partners not directly connected). Test 4 lead-lag uses a fixed $t-1$ network with varying war year (avoids the built-in-lag confound of lagging an already $t-1$-indexed variable).",
  "\\end{tablenotes}",
  "\\end{table}")
writeLines(tex, file.path(APPX_DIR, "tableA_spillover_nonfighters.tex"))
cat("  LaTeX → results/appendix/tableA_spillover_nonfighters.tex\n")

# ============================================================
# Plot — coefficient grid
# ============================================================

cat("\n--- Plotting ---\n")

plot_df <- results_rows |>
  filter(test != "main_anchor") |>
  mutate(
    test_label = recode(test, !!!test_lab),
    sample_label = recode(sample,
      all                = "Full sample",
      nonfighter         = "Non-fighter",
      nonfighter_strict  = "Non-fighter (strict)"))

p_spill <- ggplot(plot_df,
  aes(x = estimate,
      y = forcats::fct_rev(label),
      color = sig05)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray60") +
  geom_errorbarh(aes(xmin = ci_low, xmax = ci_high),
                 height = 0.18, linewidth = 0.6) +
  geom_point(size = 2.4) +
  scale_color_manual(values = c(`TRUE` = ACCENT2, `FALSE` = "gray45"),
                     labels = c("p ≥ 0.05", "p < 0.05"),
                     name = NULL) +
  facet_grid(test_label ~ sample_label, scales = "free_y", space = "free_y",
             labeller = label_wrap_gen(width = 22)) +
  labs(x = "Coefficient (95% CI)",
       y = NULL,
       title = "Spillover tests") +
  theme_pres(base_size = 11) +
  theme(strip.text.y = element_text(angle = 0, hjust = 0, size = 9),
        strip.text.x = element_text(size = 9),
        panel.spacing = unit(1.2, "lines"),
        panel.spacing.x = unit(1.5, "lines"),
        plot.margin = margin(t = 10, r = 12, b = 14, l = 10))

ggsave(file.path(FIG_DIR, "fig_spillover_nonfighters.png"),
       p_spill, width = 10, height = 13, dpi = 400, bg = "white")
cat("  -> fig_spillover_nonfighters.png\n")

cat("\nDone. 21_spillover_nonfighters.R complete.\n")
