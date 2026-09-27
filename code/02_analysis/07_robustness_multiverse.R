###############################################################################
# 03_robustness.R
# Robustness: Replicate main models across layer combos, weighting schemes,
# IVs, DVs, and transformations.
#
# SYSTEM-LEVEL GRID:
#   IVs: mid_count, log_mid_count, fatal_mid_count, log_fatal_mid_count
#         (no onset — ongoing wars reshape system topology, onset is state-level)
#   Layer combos: all subsets of {offensive, defensive, IGO, trade, DCA}
#     filtered to require IGO or defensive alliances (24 combos)
#   Weighting: contig Manhattan, contig Euclidean, dist Manhattan,
#     dist Euclidean, unweighted (5 methods)
#   Transforms: raw, log(1+x), z-score (3; zscore+euclidean skipped)
#   DVs: 5 system-level (KL, JS, edge turnover, degree KL, degree Wasserstein)
#   Models: OLS full with lagged IVs (t-1, all 5 DVs),
#     VAR/Granger fwd+rev (4 divergence/degree DVs), LP h=5 (all 5 DVs)
#
# STATE-LEVEL GRID:
#   IVs: 8 total — 4 onset (mid_onset_count, log_mid_onset_count,
#     fatal_mid_onset_count, log_fatal_mid_onset_count) + 4 count
#     (mid_count, log_mid_count, fatal_mid_count, log_fatal_mid_count)
#   Layer combos: same 24
#   Weighting: same 5 (raw transform only)
#   DV: churn (1 - weighted Jaccard)
#   Clustering: 3 schemes (ccode+year, year, ccode). Leiden community-
#   clustered SE specs were removed — they affect only inference and
#   would require per-cell Leiden recomputation across the multiverse,
#   which is computationally expensive without changing point estimates.
#   Models: FE_full (TWFE with controls), baseline (no FE, no controls,
#     focal IV + poly(year, 2) — mirrors Table 3 cols 1-2), full_noFE
#     (controls but no FE) → 7 models per IV (5 clustering × FE_full +
#     baseline + full_noFE)
#
# TEMPORAL SUBSETTING:
#   Layer combos with trade/DCA restricted to years where those layers
#   have data. Applied to system-level analyses. State-level TWFE with
#   year FE absorbs compositional level shifts.
#
# Also: meta-analysis, leave-one-layer-out, layer-specific
# Reverse Granger causality (Network → War) for all VAR specifications
#
# PLOTS: tile (layers×weighting by IV), Simonsohn spec curves (KL, churn),
#   forest (standardized effects), FEVD horizon evolution,
#   forward vs reverse Granger comparison, meta-regression coefplot
#
# Approx specs: system ~22.5k (OLS + Granger fwd/rev + LP) + state ~6.7k ≈ ~29.2k total
###############################################################################

library(tidyverse)
library(fixest)
library(igraph)
# Matrix package no longer required (spillover removed); keep available
# in case downstream meta-regression needs sparse handling.
suppressPackageStartupMessages(library(Matrix))
library(vars)
library(sandwich)
library(patchwork)
library(texreg)
library(furrr)        # parallel apply over (tfm, method, combo) triples
library(future)

select <- dplyr::select

# Parallel plan: detect cores and spin up workers. Override with the
# environment variable ROBUSTNESS_WORKERS if you want to cap manually
# (e.g. on a shared machine: ROBUSTNESS_WORKERS=4 Rscript 03_robustness.R).
.n_workers <- as.integer(Sys.getenv("ROBUSTNESS_WORKERS",
                                      max(1L, parallelly::availableCores() - 14L)))
plan(multisession, workers = .n_workers)
cat(sprintf("Parallel plan: multisession with %d workers\n", .n_workers))

# Lift the future.globals size cap. Default is 500 MiB. The helper
# functions (compute_manhattan_vec, compute_system_dvs, etc.) carry
# their enclosing environments — which include the global dyad_data
# (~1.9 GB) — so each function "weighs" ~2 GB when serialized. Lifting
# this cap to 100 GB is safe on a workstation with adequate RAM and
# silences future's protective check. With N workers each receiving
# its own copy of the globals, peak memory is .n_workers * ~2 GB.
options(future.globals.maxSize = 100 * 1024^3)   # 100 GB

# Checkpointing: save partial results to a per-batch RDS after each
# (tfm, method) dispatch. If the script crashes mid-run, restart and the
# completed batches will be loaded from disk and skipped. Set the env var
# ROBUSTNESS_RESET=1 to delete checkpoints and start fresh.
CHECKPOINT_DIR <- here::here("data", "processed", ".robustness_checkpoints")
dir.create(CHECKPOINT_DIR, recursive = TRUE, showWarnings = FALSE)
# Guarded checkpoint reset: only fires under Rscript / batch mode (where
# `interactive()` returns FALSE), never when the script is sourced into an
# interactive R session. This protects against the failure mode where a
# stale ROBUSTNESS_RESET env var (set by run_all.R Stage 0 in a previous
# run) silently wipes checkpoints when the user re-sources 03 to recover
# post-processing outputs.
if (Sys.getenv("ROBUSTNESS_RESET", "0") == "1" && !interactive()) {
  unlink(list.files(CHECKPOINT_DIR, full.names = TRUE))
  cat("Checkpoints reset (ROBUSTNESS_RESET=1).\n")
} else if (Sys.getenv("ROBUSTNESS_RESET", "0") == "1" && interactive()) {
  cat("ROBUSTNESS_RESET=1 detected but interactive() — skipping unlink ",
      "(safety guard).\n", sep = "")
}
cat(sprintf("Checkpoint dir: %s\n", CHECKPOINT_DIR))

OUTPUT_DIR <- here::here("data", "processed")
FIG_DIR <- here::here("code", "03_figures", "output")
APPX_DIR <- here::here("results", "appendix")
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(APPX_DIR, recursive = TRUE, showWarnings = FALSE)

# --------------------------------------------------------------------------- #
# Helper functions
# --------------------------------------------------------------------------- #

kl_divergence <- function(p, q, eps = 1e-10) {
  p <- p + eps; q <- q + eps
  p <- p / sum(p); q <- q / sum(q)
  sum(p * log(p / q))
}

js_divergence <- function(p, q, eps = 1e-10) {
  p <- p + eps; q <- q + eps
  p <- p / sum(p); q <- q / sum(q)
  m <- 0.5 * (p + q)
  0.5 * kl_divergence(p, m, eps = 0) + 0.5 * kl_divergence(q, m, eps = 0)
}

range01 <- function(x) {
  rng <- max(x, na.rm = TRUE) - min(x, na.rm = TRUE)
  if (rng == 0) return(rep(0, length(x)))
  (x - min(x, na.rm = TRUE)) / rng
}

# SE-ladder helper: given an lm fit and a focal term, compute the iid /
# HC1 / HC3 / NW(5) SEs, p-values, and 95% CIs in one call. Returns a
# named list ready to splat into a tibble row. Used by the multiverse
# OLS specs (system Table 1 ladder + lagged OLS + state full_noFE_year_poly)
# to keep the SE-ladder code path DRY.
compute_lm_se_ladder <- function(fit, term) {
  sm <- summary(fit)$coefficients
  if (!term %in% rownames(sm)) return(NULL)
  est_val  <- sm[term, "Estimate"]
  se_iid   <- sm[term, "Std. Error"]
  p_iid    <- sm[term, "Pr(>|t|)"]
  df_resid <- fit$df.residual
  t_crit   <- qt(0.975, df = df_resid)

  ci_lo_iid <- est_val - t_crit * se_iid
  ci_hi_iid <- est_val + t_crit * se_iid

  hc1_vv <- tryCatch(sandwich::vcovHC(fit, type = "HC1"),
                     error = function(e) NULL)
  hc1_se <- if (!is.null(hc1_vv)) sqrt(diag(hc1_vv))[term] else NA_real_
  hc1_p  <- if (!is.na(hc1_se)) {
    2 * pt(abs(est_val / hc1_se), df = df_resid, lower.tail = FALSE)
  } else NA_real_
  hc1_lo <- if (!is.na(hc1_se)) est_val - t_crit * hc1_se else NA_real_
  hc1_hi <- if (!is.na(hc1_se)) est_val + t_crit * hc1_se else NA_real_

  hc3_vv <- tryCatch(sandwich::vcovHC(fit, type = "HC3"),
                     error = function(e) NULL)
  hc3_se <- if (!is.null(hc3_vv)) sqrt(diag(hc3_vv))[term] else NA_real_
  hc3_p  <- if (!is.na(hc3_se)) {
    2 * pt(abs(est_val / hc3_se), df = df_resid, lower.tail = FALSE)
  } else NA_real_
  hc3_lo <- if (!is.na(hc3_se)) est_val - t_crit * hc3_se else NA_real_
  hc3_hi <- if (!is.na(hc3_se)) est_val + t_crit * hc3_se else NA_real_

  nw_vv <- tryCatch(sandwich::NeweyWest(fit, lag = 5,
                                         prewhite = FALSE,
                                         adjust = TRUE),
                    error = function(e) NULL)
  nw_se <- if (!is.null(nw_vv)) sqrt(diag(nw_vv))[term] else NA_real_
  nw_p  <- if (!is.na(nw_se)) {
    2 * pt(abs(est_val / nw_se), df = df_resid, lower.tail = FALSE)
  } else NA_real_
  nw_lo <- if (!is.na(nw_se)) est_val - t_crit * nw_se else NA_real_
  nw_hi <- if (!is.na(nw_se)) est_val + t_crit * nw_se else NA_real_

  list(
    coef = est_val, se = se_iid, pval = p_iid,
    ci_low_iid = ci_lo_iid, ci_high_iid = ci_hi_iid,
    hc1_se = unname(hc1_se), hc1_pval = hc1_p,
    hc1_ci_low = unname(hc1_lo), hc1_ci_high = unname(hc1_hi),
    hc3_se = unname(hc3_se), hc3_pval = hc3_p,
    hc3_ci_low = unname(hc3_lo), hc3_ci_high = unname(hc3_hi),
    nw_se = unname(nw_se), nw_pval = nw_p,
    nw_ci_low = unname(nw_lo), nw_ci_high = unname(nw_hi)
  )
}

# --------------------------------------------------------------------------- #
# 0. Load data
# --------------------------------------------------------------------------- #

dyad_data <- readRDS(file.path(OUTPUT_DIR, "dyad_year_network_layers.rds"))
metadata <- readRDS(file.path(OUTPUT_DIR, "network_metadata.rds"))
system_war <- readRDS(file.path(OUTPUT_DIR, "system_war_counts.rds"))

MAIN_LAYERS <- c("layer_offensive", "layer_defensive", "layer_igo")

all_layer_vars <- metadata$layers
all_layer_vars <- all_layer_vars[sapply(all_layer_vars, function(lv) {
  sd(dyad_data[[lv]], na.rm = TRUE) > 0
})]

# --------------------------------------------------------------------------- #
# DCA IS EXCLUDED FROM THE MULTIVERSE  [2026-07-14]
#
# Kinne's DCA data begins in 1980. Because every spec containing a restricted
# layer is estimated on that layer's coverage window (see get_valid_year_range
# below), DCA combos collapse to ~30 system-years. You cannot estimate a lagged
# OLS with controls, a VAR, or a local projection on 30 observations -- those
# specs were not evidence against the hypothesis, they were noise being counted
# as failures. They accounted for 33% of the system grid and dragged the
# significant-positive share from 50.4% to 40.3%.
#
# DCA is instead reported separately as a period-matched robustness check
# (1980-2014, DCA-augmented measure vs. the same measure without DCA on the
# same years), which is the only comparison that identifies anything.
#
# TRADE IS RETAINED. Its coverage is long enough (~144 system-years) to
# estimate on, and the year-range restriction is now applied at BOTH the
# system and state levels -- see the state loop. Previously the state loop
# computed yr_range_st but never filtered on it, so trade/DCA layers
# contributed STRUCTURAL ZEROS before their coverage began rather than
# restricting the sample. That made the churn DV a different construct before
# and after the coverage boundary, and left the "robustness" specs ~70%
# duplicates of the baseline.
# --------------------------------------------------------------------------- #
EXCLUDED_LAYERS <- c("layer_dca")
all_layer_vars  <- setdiff(all_layer_vars, EXCLUDED_LAYERS)

cat(sprintf("Main layers: %s\n", paste(MAIN_LAYERS, collapse = ", ")))
cat(sprintf("All layers: %s\n", paste(all_layer_vars, collapse = ", ")))
cat(sprintf("EXCLUDED: %s (coverage too short -- see note above)\n",
            paste(EXCLUDED_LAYERS, collapse = ", ")))

# --------------------------------------------------------------------------- #
# 0a. Determine temporal coverage for each layer
# --------------------------------------------------------------------------- #

# COVERAGE = where the layer's DATA EXISTS, not where it is positive. These
# differ, and the difference matters [FIXED 2026-07-15]:
#
#   - Offensive alliances have their last POSITIVE tie in 1990, but the data
#     span the full 1816-2014 panel. The post-1990 zeros are REAL (no offensive
#     alliances formed), not missing. Inferring coverage from `layer > 0` would
#     wrongly cap offensive combos at 1990 and discard 24 years of valid zeros.
#   - Trade data genuinely BEGIN in 1870; DCA data in 1980. Years before those
#     are true absence of measurement and must be excluded.
#
# So coverage is hard-coded from the data sources, not inferred:
PANEL_END <- 2014L
LAYER_COVERAGE <- list(
  layer_offensive = c(1816L, PANEL_END),  # ATOP; real zeros post-1990
  layer_defensive = c(1816L, PANEL_END),  # ATOP
  layer_igo       = c(1816L, PANEL_END),  # Pevehouse COW IGO
  layer_trade     = c(1870L, PANEL_END),  # COW trade begins 1870
  layer_dca       = c(1980L, 2010L)       # Kinne DCA, 1980-2010 (excluded from 07)
)

layer_year_ranges <- list()
for (lv in all_layer_vars) {
  layer_year_ranges[[lv]] <- if (!is.null(LAYER_COVERAGE[[lv]])) {
    LAYER_COVERAGE[[lv]]
  } else {
    # Fallback for any unlisted layer: infer from non-missing, non-zero data.
    ay <- dyad_data |> filter(.data[[lv]] > 0) |> pull(year)
    if (length(ay)) c(min(ay), max(ay)) else c(NA, NA)
  }
  cat(sprintf("  %s: coverage %d-%d\n", lv,
              layer_year_ranges[[lv]][1], layer_year_ranges[[lv]][2]))
}

get_valid_year_range <- function(layers) {
  start_year <- max(sapply(layers, function(lv) layer_year_ranges[[lv]][1]))
  end_year <- min(sapply(layers, function(lv) layer_year_ranges[[lv]][2]))
  c(start_year, end_year)
}

# Load system controls (log-transformed to match main analysis).
# Harmonized control set matches 01_system_level_analysis.R: includes
# n_states alongside n_democracies and n_autocracies. poly(year, 2) is
# applied inline in the regression formulas (year is already on every row).
system_analysis <- readRDS(file.path(OUTPUT_DIR, "system_analysis_data.rds"))
sys_controls <- system_analysis |>
  select(year, any_of(c("n_democracies", "n_autocracies", "n_states"))) |>
  distinct() |>
  mutate(order_period = factor(case_when(
    year <= 1913 ~ "Concert of Europe",
    year <= 1945 ~ "WWI / Interwar / WWII",
    year <= 1991 ~ "Cold War",
    TRUE         ~ "Post-Cold War"
  ),
  levels = c("Concert of Europe", "WWI / Interwar / WWII",
             "Cold War", "Post-Cold War")))

# Load state data (with controls already merged)
state_data <- readRDS(file.path(OUTPUT_DIR, "state_level_results.rds"))$data

# Build at-war lookups ONCE for combo-specific weighted PWF computation
# inside compute_state_churn(). The any-MID and fatal-MID variants mirror
# script 22 / 25's convention: PWF at year y = weighted fraction of year
# (y-1)'s partners that are at war at year y.
# These are captured by furrr::future_map_dfr workers via global capture.
at_war_lookup_any <- dyad_data |>
  filter(cowmidongoing == 1L) |>
  select(ccode1, ccode2, year) |>
  pivot_longer(c(ccode1, ccode2), values_to = "ccode") |>
  distinct(ccode, year) |>
  mutate(at_war = 1L)

at_war_lookup_fatal <- dyad_data |>
  filter(cowmidongoing == 1L, cowmidfatal == 1L) |>
  select(ccode1, ccode2, year) |>
  pivot_longer(c(ccode1, ccode2), values_to = "ccode") |>
  distinct(ccode, year) |>
  mutate(at_war_fatal = 1L)

# Add order_period factor — used by the state_period_fe model class
# (folded in from former 12_order_period_fe.R). Same 4-regime
# periodization as 01's Table 1 col 4.
state_data <- state_data |>
  mutate(order_period = factor(case_when(
    year <= 1913 ~ "Concert of Europe",
    year <= 1945 ~ "WWI / Interwar / WWII",
    year <= 1991 ~ "Cold War",
    TRUE         ~ "Post-Cold War"
  ),
  levels = c("Concert of Europe", "WWI / Interwar / WWII",
             "Cold War", "Post-Cold War")))

# NOTE: spillover specifications (inv_dist_to_war, 2-hop partner war) used
# to be part of this multiverse but were removed when the manuscript scope
# narrowed. Spillover analysis now lives entirely in
# 21_spillover_nonfighters.R, which consumes the canonical RDS produced by
# 01_data_cleaning/03_compute_spillover_canonical.R.

# --------------------------------------------------------------------------- #
# 1. Compute edge weight for a given layer combo + weighting method
# --------------------------------------------------------------------------- #

compute_manhattan_vec <- function(data, layers, method) {
  mat <- as.matrix(data[, layers, drop = FALSE])
  mat[is.na(mat)] <- 0

  if (method == "contig_manhattan") {
    weighted_cols <- intersect(paste0(layers, "_weighted"), names(data))
    if (length(weighted_cols) == 0) return(rep(0, nrow(data)))
    mat_w <- as.matrix(data[, weighted_cols, drop = FALSE])
    mat_w[is.na(mat_w)] <- 0
    result <- rowSums(mat_w) / length(weighted_cols)
  } else if (method == "contig_euclidean") {
    weighted_cols <- intersect(paste0(layers, "_weighted"), names(data))
    if (length(weighted_cols) == 0) return(rep(0, nrow(data)))
    mat_w <- as.matrix(data[, weighted_cols, drop = FALSE])
    mat_w[is.na(mat_w)] <- 0
    mat_w[mat_w < 0] <- 0
    result <- sqrt(rowSums(mat_w^2)) / length(weighted_cols)
  } else if (method == "dist_manhattan") {
    dist_cols <- intersect(paste0(layers, "_dist_weighted"), names(data))
    if (length(dist_cols) == 0) return(rep(0, nrow(data)))
    mat_d <- as.matrix(data[, dist_cols, drop = FALSE])
    mat_d[is.na(mat_d)] <- 0
    result <- rowSums(mat_d) / length(dist_cols)
  } else if (method == "dist_euclidean") {
    dist_cols <- intersect(paste0(layers, "_dist_weighted"), names(data))
    if (length(dist_cols) == 0) return(rep(0, nrow(data)))
    mat_d <- as.matrix(data[, dist_cols, drop = FALSE])
    mat_d[is.na(mat_d)] <- 0
    mat_d[mat_d < 0] <- 0
    result <- sqrt(rowSums(mat_d^2)) / length(dist_cols)
  } else if (method == "unweighted") {
    result <- rowSums(mat)
  }

  # Normalize by year
  years <- data$year
  for (yr in unique(years)) {
    idx <- which(years == yr)
    result[idx] <- range01(result[idx])
  }
  result
}

# --------------------------------------------------------------------------- #
# 2. Compute system-level divergence series
# --------------------------------------------------------------------------- #

compute_system_dvs <- function(data, manhattan_vec) {
  data$w_alt <- manhattan_vec
  years <- sort(unique(data$year))

  out <- tibble(
    year = years,
    kl_divergence = NA_real_,
    js_divergence = NA_real_,
    edge_turnover_pct = NA_real_,
    degree_kl = NA_real_,
    degree_wasserstein = NA_real_
  )

  wasserstein_dist <- function(a, b) {
    a_s <- sort(a); b_s <- sort(b)
    n <- max(length(a_s), length(b_s))
    a_q <- quantile(a_s, probs = seq(0, 1, length.out = n), type = 1)
    b_q <- quantile(b_s, probs = seq(0, 1, length.out = n), type = 1)
    mean(abs(a_q - b_q))
  }

  for (i in 2:length(years)) {
    e_t  <- data[data$year == years[i] & data$w_alt > 0,
                 c("ccode1", "ccode2", "w_alt")]
    e_t1 <- data[data$year == years[i-1] & data$w_alt > 0,
                 c("ccode1", "ccode2", "w_alt")]

    if (nrow(e_t) == 0 | nrow(e_t1) == 0) next

    merged <- merge(e_t, e_t1, by = c("ccode1", "ccode2"), all = TRUE,
                    suffixes = c("_t", "_t1"))
    merged$w_alt_t[is.na(merged$w_alt_t)] <- 0
    merged$w_alt_t1[is.na(merged$w_alt_t1)] <- 0
    if (min(merged$w_alt_t) < 0) merged$w_alt_t <- merged$w_alt_t - min(merged$w_alt_t)
    if (min(merged$w_alt_t1) < 0) merged$w_alt_t1 <- merged$w_alt_t1 - min(merged$w_alt_t1)

    out$kl_divergence[i] <- kl_divergence(merged$w_alt_t, merged$w_alt_t1)
    out$js_divergence[i] <- js_divergence(merged$w_alt_t, merged$w_alt_t1)

    # Edge turnover
    dk <- function(df) paste(pmin(df$ccode1, df$ccode2), pmax(df$ccode1, df$ccode2), sep = "-")
    s_t <- dk(e_t); s_t1 <- dk(e_t1)
    turnover <- length(setdiff(s_t, s_t1)) + length(setdiff(s_t1, s_t))
    out$edge_turnover_pct[i] <- turnover / max(length(union(s_t, s_t1)), 1)

    # Graph-level measures
    g_t <- tryCatch({
      g <- graph_from_data_frame(e_t[, c("ccode1", "ccode2")], directed = FALSE)
      E(g)$weight <- e_t$w_alt; g
    }, error = function(e) NULL)
    g_t1 <- tryCatch({
      g <- graph_from_data_frame(e_t1[, c("ccode1", "ccode2")], directed = FALSE)
      E(g)$weight <- e_t1$w_alt; g
    }, error = function(e) NULL)

    if (!is.null(g_t) && !is.null(g_t1)) {
      n_t <- vcount(g_t); n_t1 <- vcount(g_t1)

      # Degree distributions
      deg_t <- degree(g_t); deg_t1 <- degree(g_t1)
      if (length(deg_t) > 1 && length(deg_t1) > 1) {
        out$degree_wasserstein[i] <- wasserstein_dist(deg_t, deg_t1)
        max_deg <- max(c(deg_t, deg_t1)) + 1
        hist_t <- tabulate(deg_t + 1, nbins = max_deg)
        hist_t1 <- tabulate(deg_t1 + 1, nbins = max_deg)
        out$degree_kl[i] <- kl_divergence(hist_t, hist_t1)
      }

      # (modularity, centralization, density removed — not theoretically predicted)
    }
  }

  # Logged variants of every system DV. These are the canonical form used in
  # the harmonized 01 main spec; the spec curve scans both raw and logged.
  out$log_kl_divergence       <- log1p(out$kl_divergence)
  out$log_js_divergence       <- log1p(out$js_divergence)
  out$log_edge_turnover_pct   <- log1p(out$edge_turnover_pct)
  out$log_degree_kl           <- log1p(out$degree_kl)
  out$log_degree_wasserstein  <- log1p(out$degree_wasserstein)

  out
}

# --------------------------------------------------------------------------- #
# 3. Compute state-level churn for a given layer combo
# --------------------------------------------------------------------------- #

compute_state_churn <- function(data, manhattan_vec,
                                  at_war_any   = at_war_lookup_any,
                                  at_war_fatal = at_war_lookup_fatal) {
  data$w_alt <- manhattan_vec

  # Build undirected ties
  ties <- data |>
    select(ccode1, ccode2, year, w_alt) |>
    filter(w_alt > 0) |>
    bind_rows(
      data |> filter(w_alt > 0) |>
        select(ccode1 = ccode2, ccode2 = ccode1, year, w_alt)
    ) |>
    group_by(ccode1, ccode2, year) |>
    summarise(w = max(w_alt, na.rm = TRUE), .groups = "drop")

  # --- Combo-specific weighted PWF (any-MID and fatal-MID) ---
  # PWF at year y = weighted fraction of year (y-1)'s partners that are
  # at war at year y. Computed using THIS combo's edge weights so the
  # control is internally consistent with the cell's churn outcome.
  pwf_any <- ties |>
    mutate(year_next = year + 1L) |>
    left_join(at_war_any,
              by = c("ccode2" = "ccode", "year_next" = "year")) |>
    mutate(at_war = replace_na(at_war, 0L)) |>
    group_by(ccode1, year_next) |>
    summarise(pwf_weighted_combo = sum(at_war * w) / sum(w),
              .groups = "drop") |>
    rename(ccode = ccode1, year = year_next)

  pwf_fatal <- ties |>
    mutate(year_next = year + 1L) |>
    left_join(at_war_fatal,
              by = c("ccode2" = "ccode", "year_next" = "year")) |>
    mutate(at_war_fatal = replace_na(at_war_fatal, 0L)) |>
    group_by(ccode1, year_next) |>
    summarise(pwf_fatal_weighted_combo = sum(at_war_fatal * w) / sum(w),
              .groups = "drop") |>
    rename(ccode = ccode1, year = year_next)

  # --- Churn per (state, year) ---
  states <- unique(ties$ccode1)
  result_list <- list()

  for (st in states) {
    st_ties <- ties |> filter(ccode1 == st)
    st_years <- sort(unique(st_ties$year))
    if (length(st_years) < 2) next

    for (j in 2:length(st_years)) {
      yr_t  <- st_years[j]
      yr_t1 <- st_years[j - 1]
      if (yr_t - yr_t1 != 1) next

      w_t  <- st_ties |> filter(year == yr_t)  |> select(ccode2, w)
      w_t1 <- st_ties |> filter(year == yr_t1) |> select(ccode2, w)

      # Churn (1 - weighted Jaccard)
      w_both <- full_join(w_t, w_t1, by = "ccode2", suffix = c("_t", "_t1")) |>
        mutate(w_t = replace_na(w_t, 0), w_t1 = replace_na(w_t1, 0))
      sum_min <- sum(pmin(w_both$w_t, w_both$w_t1))
      sum_max <- sum(pmax(w_both$w_t, w_both$w_t1))
      churn <- if (sum_max > 0) 1 - sum_min / sum_max else NA_real_

      result_list[[length(result_list) + 1]] <- tibble(
        ccode = st, year = yr_t,
        churn_alt     = churn,
        log_churn_alt = log1p(churn)
      )
    }
  }
  churn_df <- bind_rows(result_list)

  # Merge in combo-specific PWF columns
  churn_df |>
    left_join(pwf_any,   by = c("ccode", "year")) |>
    left_join(pwf_fatal, by = c("ccode", "year")) |>
    mutate(pwf_weighted_combo       = replace_na(pwf_weighted_combo, 0),
           pwf_fatal_weighted_combo = replace_na(pwf_fatal_weighted_combo, 0))
}

# --------------------------------------------------------------------------- #
# 3b. (REMOVED) Per-combo spillover helper
# --------------------------------------------------------------------------- #
# `compute_spillover_for_combo()` and `at_war_state_global` were removed when
# the spillover analyses were pulled out of the multiverse. Spillover now
# lives only in 21_spillover_nonfighters.R, which uses the canonical RDS
# produced by 01_data_cleaning/03_compute_spillover_canonical.R.

# --------------------------------------------------------------------------- #
# 4. Generate all specifications
# --------------------------------------------------------------------------- #

cat("=== Generating specifications ===\n")

weighting_methods <- c("contig_manhattan", "contig_euclidean",
                       "dist_manhattan", "dist_euclidean", "unweighted")

# Layer combos: all subsets, filtered to require IGO or defensive
all_combos_raw <- list()
for (k in 1:length(all_layer_vars)) {
  combos <- combn(all_layer_vars, k, simplify = FALSE)
  all_combos_raw <- c(all_combos_raw, combos)
}

# Also add main + extras combos
extra_layers <- setdiff(all_layer_vars, MAIN_LAYERS)
if (length(extra_layers) > 0) {
  for (k in 1:length(extra_layers)) {
    extras <- combn(extra_layers, k, simplify = FALSE)
    for (ex in extras) {
      all_combos_raw <- c(all_combos_raw, list(c(MAIN_LAYERS, ex)))
    }
  }
}

all_combos_raw <- unique(lapply(all_combos_raw, sort))

all_combos <- Filter(function(combo) {
  "layer_igo" %in% combo | "layer_defensive" %in% combo
}, all_combos_raw)

transform_methods <- c("raw", "log1p", "zscore")

# System IVs: counts only (no onset), raw + log
system_war_vars <- c("mid_count", "log_mid_count",
                     "fatal_mid_count", "log_fatal_mid_count")

# State IVs: onset + count, raw + log, AND partner-war-fraction variants.
# FULL SWEEP retained for the baseline (no extra control). The GDP/strength
# control dimension is what tripled the state grid, so it is restricted -- see
# GDP_STRENGTH_IVS / GDP_STRENGTH_DVS below and the guard in the FE_full loop.
state_war_vars <- c("mid_onset_count", "log_mid_onset_count",
                    "fatal_mid_onset_count", "log_fatal_mid_onset_count",
                    "mid_count", "log_mid_count",
                    "fatal_mid_count", "log_fatal_mid_count",
                    "partner_war_frac", "partner_war_frac_weighted",
                    "partner_fatal_war_frac", "partner_fatal_war_frac_weighted")

# GDP and strength controls run ONLY on the theoretically-motivated
# operationalizations [2026-07-15]: logged MID onset/count (H2) and both
# partner-war fractions (H3), on the logged churn DVs. Rationale: fatal-MID
# severity is a separate appendix section; raw counts are dominated by the
# logged form. This keeps the main sweep comprehensive while making the
# GDP/strength robustness cheap (it does not need every mechanical variant).
GDP_STRENGTH_IVS <- c("log_mid_onset_count", "log_mid_count",
                      "partner_war_frac", "partner_war_frac_weighted")
GDP_STRENGTH_DVS <- c("log_churn_alt", "log_churn")

# System DVs: 5 raw + 5 log1p variants. The harmonized 01 main spec uses
# log_kl_divergence; 03 scans both forms so the spec curve shows the
# distribution under each DV transform.
system_dvs <- c(
  # raw
  "kl_divergence", "js_divergence", "edge_turnover_pct",
  "degree_kl", "degree_wasserstein",
  # log1p
  "log_kl_divergence", "log_js_divergence", "log_edge_turnover_pct",
  "log_degree_kl", "log_degree_wasserstein"
)

# DVs eligible for VAR/Granger (distribution-based measures, raw + log)
granger_dvs <- c(
  "kl_divergence", "js_divergence", "degree_kl", "degree_wasserstein",
  "log_kl_divergence", "log_js_divergence", "log_degree_kl", "log_degree_wasserstein"
)

# State DVs: 4 = {churn_alt, log_churn_alt, churn, log_churn}. FULL SET retained
# for the baseline sweep. The GDP/strength specs use only the logged pair
# (GDP_STRENGTH_DVS) via the guard in the FE_full loop.
#   churn_alt / log_churn_alt: per-combo recomputed (varies with layer+weighting).
#   churn / log_churn: canonical from 02 (constant across layer combos).
state_dvs <- c("churn_alt", "log_churn_alt", "churn", "log_churn")

# State clustering schemes (Leiden community-clustered SE removed —
# would require per-cell Leiden recomputation across the multiverse;
# the three remaining schemes cover the canonical inference angles).
clustering_schemes <- list(
  "ccode_year" = ~ccode + year,
  "year_only"  = ~year,
  "ccode_only" = ~ccode
)

# Harmonization: main spec uses the UNWEIGHTED partner war fraction,
# routed by IV severity — mid-count IVs use partner_war_frac (any MID),
# fatal-MID IVs use partner_fatal_war_frac (only fatal MIDs). The
# tie-strength-weighted versions are kept as alternative specs (looped
# over for the robustness multiverse). Both severities × both weighting
# variants are stored on every row of state_data.
#
# get_control_sets(war_var) returns a 2-element list keyed by
# `unweighted` (main, anchors headline) and `weighted` (alt, appendix).
# Each entry uses the partner-war variable that mirrors the IV's
# severity level.
#
# is_isolated_t1: predetermined indicator for "state had no partners at
# t-1". Absorbs the level shift from `total_gain` and `both_isolated`
# state-years (degree_t1 = 0) introduced by the complete-COW-panel
# refactor in 02_compute_network_measures.R, so the within-state IV
# coefficient is identified primarily from `normal` and `total_loss` rows.
get_control_sets <- function(war_var) {
  base       <- "+ democracy + autocracy + cinc + major_power + is_isolated_t1"
  base_clean <- "democracy + autocracy + cinc + major_power + is_isolated_t1"
  # When PWF is itself the IV, drop the PWF control (perfectly collinear).
  pwf_ivs <- c("partner_war_frac", "partner_war_frac_weighted",
               "partner_fatal_war_frac", "partner_fatal_war_frac_weighted")
  if (war_var %in% pwf_ivs) {
    return(list(self = base_clean))
  }
  if (grepl("fatal", war_var)) {
    list(
      unweighted = paste0("partner_fatal_war_frac ", base),
      weighted   = paste0("partner_fatal_war_frac_weighted ", base)
    )
  } else {
    list(
      unweighted = paste0("partner_war_frac ", base),
      weighted   = paste0("partner_war_frac_weighted ", base)
    )
  }
}
# Static reference for any code that still expects a top-level binding;
# the IV-routed version above is the source of truth.
control_sets <- get_control_sets("mid_count")

# --------------------------------------------------------------------------- #
# EXTRA CONTROL SETS  [ADDED 2026-07-14]
#
# Folds in the two scripts that used to run these as standalone robustness:
#   - former 26_robustness_gdp_control.R      (log GDP)
#   - former 25_robustness_strength_control.R (log network strength at t-1)
#
# Applied to the FE_full* model class ONLY. The point of these controls is
# whether the POINT ESTIMATE survives them; running them across the no-FE and
# bivariate classes would triple the grid to no purpose.
#
# `strength` DROPS is_isolated_t1. log1p(strength_t1) already encodes isolation
# (strength = 0 -> log = 0), so the two are collinear by construction. This
# mirrors what 25 did.
#
# log_strength_t1 is COMBO-SPECIFIC -- it is the sum of THIS cell's edge
# weights at t-1, so it is recomputed inside the combo loop, not merged in
# from a canonical panel. log_gdp is a state-year attribute and is merged once.
# --------------------------------------------------------------------------- #
EXTRA_CTRLS <- c("none", "gdp", "strength")

augment_controls <- function(rhs, extra) {
  if (extra == "none")     return(rhs)
  if (extra == "gdp")      return(paste0(rhs, " + log_gdp"))
  if (extra == "strength") {
    rhs2 <- gsub("\\s*\\+\\s*is_isolated_t1", "", rhs)   # collinear -- drop
    return(paste0(rhs2, " + log_strength_t1"))
  }
  stop("unknown extra control: ", extra)
}

# Combo-specific network strength at t-1 (mirrors former 25's construction).
compute_state_strength <- function(data, manhattan_vec) {
  d <- tibble::tibble(ccode1 = data$ccode1, ccode2 = data$ccode2,
                      year = data$year, w = manhattan_vec)
  dplyr::bind_rows(
    d |> dplyr::select(ccode = ccode1, year, w),
    d |> dplyr::select(ccode = ccode2, year, w)
  ) |>
    dplyr::group_by(ccode, year) |>
    dplyr::summarise(strength = sum(w, na.rm = TRUE), .groups = "drop") |>
    dplyr::mutate(year_next = year + 1L) |>
    dplyr::select(ccode, year = year_next, strength_t1 = strength) |>
    dplyr::mutate(log_strength_t1 = log1p(strength_t1))
}

# --- log GDP (state-year; merged once) ------------------------------------- #
# Markowitz/Fariss simulated real GDP via peacesciencer. Coverage is
# incomplete, so `gdp` specs will drop observations -- that is expected and is
# recorded in n_obs.
cat("Adding log GDP via peacesciencer::add_sim_gdp_pop()...\n")
gdp_panel <- tryCatch({
  sy <- peacesciencer::create_stateyears(system = "cow", subset_years = 1816:2014)
  attr(sy, "ps_data_type") <- "state_year"
  attr(sy, "ps_system")    <- "cow"
  sy |>
    peacesciencer::add_sim_gdp_pop() |>
    dplyr::select(ccode, year, pwtrgdp) |>
    dplyr::mutate(log_gdp = log1p(pwtrgdp)) |>
    dplyr::select(ccode, year, log_gdp)
}, error = function(e) {
  cat("  WARNING: GDP merge failed (", conditionMessage(e), ")\n")
  cat("  -> `gdp` specs will be skipped.\n")
  NULL
})
if (!is.null(gdp_panel)) {
  state_data <- state_data |> dplyr::left_join(gdp_panel, by = c("ccode", "year"))
  cat(sprintf("  log_gdp non-missing: %.1f%% of state-years\n",
              100 * mean(!is.na(state_data$log_gdp))))
} else {
  state_data$log_gdp <- NA_real_
  EXTRA_CTRLS <- setdiff(EXTRA_CTRLS, "gdp")
}

cat(sprintf("  Layer combos: %d\n", length(all_combos)))
cat(sprintf("  Weighting: %d methods\n", length(weighting_methods)))
cat(sprintf("  Transforms: %s\n", paste(transform_methods, collapse = ", ")))
cat(sprintf("  System IVs: %s\n", paste(system_war_vars, collapse = ", ")))
cat(sprintf("  State IVs: %s\n", paste(state_war_vars, collapse = ", ")))
cat(sprintf("  System DVs: %d\n", length(system_dvs)))
cat(sprintf("  State DVs: %s\n", paste(state_dvs, collapse = ", ")))
cat(sprintf("  Clustering schemes: %d\n", length(clustering_schemes)))

# Ensure log war vars exist in system_war
if (!"log_mid_count" %in% names(system_war)) {
  system_war <- system_war |>
    mutate(
      log_mid_count = log1p(mid_count),
      log_fatal_mid_count = log1p(fatal_mid_count)
    )
}

# Ensure log war vars exist in state_data
if (!"log_mid_onset_count" %in% names(state_data)) {
  state_data <- state_data |>
    mutate(
      log_mid_onset_count = log1p(mid_onset_count),
      log_fatal_mid_onset_count = log1p(fatal_mid_onset_count)
    )
}
if (!"log_mid_count" %in% names(state_data)) {
  state_data <- state_data |>
    mutate(
      log_mid_count = log1p(mid_count),
      log_fatal_mid_count = log1p(fatal_mid_count)
    )
}
# Ensure log_churn exists for logged DV spec curve
if (!"log_churn" %in% names(state_data)) {
  state_data <- state_data |>
    mutate(log_churn = log1p(churn))
}

# =========================================================================== #
# 5. SYSTEM-LEVEL ROBUSTNESS (Lagged OLS + VAR/Granger + LP)
# =========================================================================== #

cat("\n=== SYSTEM-LEVEL ROBUSTNESS ===\n")

spec_results <- list()
spec_idx <- 0
analysis_data_main <- NULL  # will be populated when unweighted/raw/defensive+igo is hit

for (tfm in transform_methods) {
  # Apply layer transform
  if (tfm == "log1p") {
    continuous_layers <- intersect(c("layer_igo", "layer_trade"), all_layer_vars)
    dyad_work <- dyad_data
    for (cl in continuous_layers) dyad_work[[cl]] <- log1p(dyad_work[[cl]])
  } else if (tfm == "zscore") {
    dyad_work <- dyad_data |>
      group_by(year) |>
      mutate(across(all_of(all_layer_vars), ~ {
        s <- sd(., na.rm = TRUE)
        if (is.na(s) || s == 0) 0 else (. - mean(., na.rm = TRUE)) / s
      })) |>
      ungroup()
  } else {
    dyad_work <- dyad_data
  }

  for (method in weighting_methods) {
    if (tfm == "zscore" && grepl("euclidean", method)) next

    # Checkpoint: skip this batch if it's already been computed.
    # Memory hygiene — only read the row count, don't load the full tibble
    # into the parent process. The full tibble is reloaded at the
    # post-processing block from disk.
    ckpt_path <- file.path(CHECKPOINT_DIR,
                            sprintf("system__%s__%s.rds", tfm, method))
    if (file.exists(ckpt_path)) {
      n_loaded <- nrow(readRDS(ckpt_path))
      cat(sprintf("  [%s/%s] LOADED CHECKPOINT (%d specs)\n",
                  tfm, method, n_loaded))
      spec_idx <- spec_idx + n_loaded
      gc(verbose = FALSE)
      next
    }

    cat(sprintf("  [%s/%s] dispatching %d combos in parallel (%d workers)...\n",
                tfm, method, length(all_combos), .n_workers))

    triple_dfs <- furrr::future_map_dfr(seq_along(all_combos),
                                          function(combo_idx) {
      out_list <- list()
      layers_subset <- all_combos[[combo_idx]]
      combo_name <- paste(gsub("layer_", "", layers_subset), collapse = "+")

      manhattan_vec <- compute_manhattan_vec(dyad_work, layers_subset, method)
      dv_df <- compute_system_dvs(dyad_work, manhattan_vec)

      # Temporal subsetting
      yr_range <- get_valid_year_range(layers_subset)

      analysis_data <- dv_df |>
        left_join(system_war, by = "year") |>
        left_join(sys_controls, by = "year") |>
        filter(year >= yr_range[1], year <= yr_range[2]) |>
        filter(!is.na(kl_divergence)) |>
        mutate(
          log_kl_divergence      = log1p(kl_divergence),
          log_js_divergence      = log1p(js_divergence),
          log_edge_turnover_pct  = log1p(edge_turnover_pct),
          log_degree_kl          = log1p(degree_kl),
          log_degree_wasserstein = log1p(degree_wasserstein)
        )

      # Skip combo if too few obs — return empty bound rows to future_map_dfr
      if (nrow(analysis_data) < 20) return(dplyr::bind_rows(out_list))

      has_offensive <- "layer_offensive" %in% layers_subset
      has_defensive <- "layer_defensive" %in% layers_subset
      has_igo <- "layer_igo" %in% layers_subset
      has_trade <- "layer_trade" %in% layers_subset
      has_dca <- "layer_dca" %in% layers_subset
      year_start <- yr_range[1]
      year_end <- yr_range[2]

      # NOTE: analysis_data_main (used by FEVD horizon plot below) is now
      # re-computed sequentially after the parallel system loop finishes —
      # the parallel workers can't propagate side-effects back to the parent.

      # Create lagged IVs for OLS (t-1, matching main analysis)
      analysis_data <- analysis_data |>
        arrange(year) |>
        mutate(
          mid_count_lag1 = lag(mid_count),
          log_mid_count_lag1 = lag(log_mid_count),
          fatal_mid_count_lag1 = lag(fatal_mid_count),
          log_fatal_mid_count_lag1 = lag(log_fatal_mid_count)
        )

      ols_lag_vars <- c("mid_count_lag1", "log_mid_count_lag1",
                        "fatal_mid_count_lag1", "log_fatal_mid_count_lag1")

      # Map lagged IV names back to their war_var name for consistent labeling
      ols_lag_to_war <- c(
        mid_count_lag1 = "mid_count",
        log_mid_count_lag1 = "log_mid_count",
        fatal_mid_count_lag1 = "fatal_mid_count",
        log_fatal_mid_count_lag1 = "log_fatal_mid_count"
      )

      # ── Lagged OLS full (all 8 DVs) ──
      for (ols_iv in ols_lag_vars) {
        if (all(is.na(analysis_data[[ols_iv]]))) next

        for (dv in system_dvs) {
          if (all(is.na(analysis_data[[dv]]))) next

          m_ols <- tryCatch(
            lm(as.formula(paste(dv, "~", ols_iv,
                                "+ log1p(n_democracies) + log1p(n_autocracies)",
                                "+ log1p(n_states) + poly(year, 2)")),
               data = analysis_data),
            error = function(e) NULL
          )
          if (!is.null(m_ols)) {
            sm <- summary(m_ols)$coefficients
            if (ols_iv %in% rownames(sm)) {
              # SE ladder: iid (lm default) + HC1 (White) + HC3
              # (jackknife) + Newey-West HAC(5). HC1 is the standard
              # White-style robust SE; HC3 weights residuals by
              # 1/(1 - h_ii)^2 (Long & Ervin 2000 recommend for N<250);
              # NW(5) handles serial correlation in the time-series
              # residuals. All via sandwich, no refitting. The four-way
              # SE ladder gives the spec curve a complete inference
              # picture for each cell of the multiverse.
              df_resid <- m_ols$df.residual
              t_crit   <- qt(0.975, df = df_resid)
              est_val  <- sm[ols_iv, "Estimate"]

              hc1_vv <- tryCatch(sandwich::vcovHC(m_ols, type = "HC1"),
                                  error = function(e) NULL)
              hc1_se_val <- if (!is.null(hc1_vv)) sqrt(diag(hc1_vv))[ols_iv] else NA_real_
              hc1_p_val  <- if (!is.na(hc1_se_val)) {
                2 * pt(abs(est_val / hc1_se_val), df = df_resid, lower.tail = FALSE)
              } else NA_real_
              hc1_ci_lo  <- if (!is.na(hc1_se_val)) est_val - t_crit * hc1_se_val else NA_real_
              hc1_ci_hi  <- if (!is.na(hc1_se_val)) est_val + t_crit * hc1_se_val else NA_real_

              hc3_vv <- tryCatch(sandwich::vcovHC(m_ols, type = "HC3"),
                                  error = function(e) NULL)
              hc3_se_val <- if (!is.null(hc3_vv)) sqrt(diag(hc3_vv))[ols_iv] else NA_real_
              hc3_p_val  <- if (!is.na(hc3_se_val)) {
                2 * pt(abs(est_val / hc3_se_val), df = df_resid, lower.tail = FALSE)
              } else NA_real_
              hc3_ci_lo  <- if (!is.na(hc3_se_val)) est_val - t_crit * hc3_se_val else NA_real_
              hc3_ci_hi  <- if (!is.na(hc3_se_val)) est_val + t_crit * hc3_se_val else NA_real_

              nw_vv <- tryCatch(sandwich::NeweyWest(m_ols, lag = 5,
                                                      prewhite = FALSE,
                                                      adjust = TRUE),
                                  error = function(e) NULL)
              nw_se_val <- if (!is.null(nw_vv)) sqrt(diag(nw_vv))[ols_iv] else NA_real_
              nw_p_val  <- if (!is.na(nw_se_val)) {
                2 * pt(abs(est_val / nw_se_val), df = df_resid, lower.tail = FALSE)
              } else NA_real_
              nw_ci_lo  <- if (!is.na(nw_se_val)) est_val - t_crit * nw_se_val else NA_real_
              nw_ci_hi  <- if (!is.na(nw_se_val)) est_val + t_crit * nw_se_val else NA_real_

              ci_lo_iid <- est_val - t_crit * sm[ols_iv, "Std. Error"]
              ci_hi_iid <- est_val + t_crit * sm[ols_iv, "Std. Error"]

              out_list[[length(out_list) + 1]] <- tibble(
                level = "system", dv = dv,
                war_variable = ols_lag_to_war[[ols_iv]],
                weighting = method, transform = tfm, layers = combo_name,
                n_layers = length(layers_subset), model = "OLS_lagged",
                coef = est_val,
                se = sm[ols_iv, "Std. Error"],
                pval = sm[ols_iv, "Pr(>|t|)"],
                ci_low_iid = ci_lo_iid, ci_high_iid = ci_hi_iid,
                hc1_se = unname(hc1_se_val), hc1_pval = hc1_p_val,
                hc1_ci_low = unname(hc1_ci_lo), hc1_ci_high = unname(hc1_ci_hi),
                hc3_se = unname(hc3_se_val), hc3_pval = hc3_p_val,
                hc3_ci_low = unname(hc3_ci_lo), hc3_ci_high = unname(hc3_ci_hi),
                nw_se = unname(nw_se_val), nw_pval = nw_p_val,
                nw_ci_low = unname(nw_ci_lo), nw_ci_high = unname(nw_ci_hi),
                n_obs = sum(!is.na(analysis_data[[ols_iv]]) &
                            !is.na(analysis_data[[dv]])),
                has_offensive = has_offensive, has_defensive = has_defensive,
                has_igo = has_igo, has_trade = has_trade, has_dca = has_dca,
                is_onset = FALSE, year_start = year_start, year_end = year_end
              )
            }
          }
        }
      }

      # ── Manuscript Table 1 (system) ladder ──
      # Folds the four-column main system table from
      # 01_system_level_analysis.R into the multiverse so each (layers,
      # weighting, transform) cell yields the full ladder. Col 2 of
      # Table 1 (lagged IV + full controls + poly(year, 2)) is already
      # captured by the OLS_lagged block above; this block adds:
      #   - sys_T1_iv_only      = IV + poly(year, 2)
      #   - sys_T1_ar1          = OLS_lagged spec + log_kl_lag1 (AR(1))
      #   - sys_T1_order_period = OLS_lagged spec but with order_period
      #                            FE replacing the demographic ctrls'
      #                            secular trend (order_period absorbs
      #                            between-period mean shifts).
      # All three use the same NW(5) HAC SE ladder (HC1 + HC3 + NW)
      # via compute_lm_se_ladder() so they appear cleanly in the spec
      # curve and meta-regression alongside OLS_lagged.

      # Build lagged DV columns (AR(1) input). Only computed once per
      # combo since they depend on (analysis_data, dv) — recompute
      # inside the dv loop below.
      sys_t1_specs <- list(
        list(key = "sys_T1_iv_only",
             rhs = function(iv) paste(iv, "+ poly(year, 2)")),
        list(key = "sys_T1_order_period",
             rhs = function(iv) paste(iv,
                                       "+ log1p(n_democracies) + log1p(n_autocracies)",
                                       "+ log1p(n_states) + order_period"))
      )

      for (ols_iv in ols_lag_vars) {
        if (all(is.na(analysis_data[[ols_iv]]))) next

        for (dv in system_dvs) {
          if (all(is.na(analysis_data[[dv]]))) next

          # Cols 1 and 4: standard ladder via helper.
          for (sp in sys_t1_specs) {
            m_t1 <- tryCatch(
              lm(as.formula(paste(dv, "~", sp$rhs(ols_iv))),
                 data = analysis_data),
              error = function(e) NULL
            )
            if (is.null(m_t1)) next
            ladder <- compute_lm_se_ladder(m_t1, ols_iv)
            if (is.null(ladder)) next

            out_list[[length(out_list) + 1]] <- tibble(
              level = "system", dv = dv,
              war_variable = ols_lag_to_war[[ols_iv]],
              weighting = method, transform = tfm, layers = combo_name,
              n_layers = length(layers_subset), model = sp$key,
              coef       = ladder$coef,
              se         = ladder$se,
              pval       = ladder$pval,
              ci_low_iid = ladder$ci_low_iid, ci_high_iid = ladder$ci_high_iid,
              hc1_se     = ladder$hc1_se,     hc1_pval     = ladder$hc1_pval,
              hc1_ci_low = ladder$hc1_ci_low, hc1_ci_high  = ladder$hc1_ci_high,
              hc3_se     = ladder$hc3_se,     hc3_pval     = ladder$hc3_pval,
              hc3_ci_low = ladder$hc3_ci_low, hc3_ci_high  = ladder$hc3_ci_high,
              nw_se      = ladder$nw_se,      nw_pval      = ladder$nw_pval,
              nw_ci_low  = ladder$nw_ci_low,  nw_ci_high   = ladder$nw_ci_high,
              n_obs = sum(!is.na(analysis_data[[ols_iv]]) &
                          !is.na(analysis_data[[dv]])),
              has_offensive = has_offensive, has_defensive = has_defensive,
              has_igo = has_igo, has_trade = has_trade, has_dca = has_dca,
              is_onset = FALSE, year_start = year_start, year_end = year_end
            )
          }

          # Col 3: + AR(1) lagged DV. Lag is per-DV so build it inline.
          dv_lag_col <- paste0(dv, "_lag1_internal")
          ar_data <- analysis_data |>
            arrange(year) |>
            mutate(!!dv_lag_col := dplyr::lag(.data[[dv]]))
          ar_rhs <- paste(ols_iv,
                          "+ log1p(n_democracies) + log1p(n_autocracies)",
                          "+ log1p(n_states) + poly(year, 2) +",
                          dv_lag_col)
          m_t1_ar <- tryCatch(
            lm(as.formula(paste(dv, "~", ar_rhs)), data = ar_data),
            error = function(e) NULL
          )
          if (!is.null(m_t1_ar)) {
            ladder_ar <- compute_lm_se_ladder(m_t1_ar, ols_iv)
            if (!is.null(ladder_ar)) {
              out_list[[length(out_list) + 1]] <- tibble(
                level = "system", dv = dv,
                war_variable = ols_lag_to_war[[ols_iv]],
                weighting = method, transform = tfm, layers = combo_name,
                n_layers = length(layers_subset), model = "sys_T1_ar1",
                coef       = ladder_ar$coef,
                se         = ladder_ar$se,
                pval       = ladder_ar$pval,
                ci_low_iid = ladder_ar$ci_low_iid, ci_high_iid = ladder_ar$ci_high_iid,
                hc1_se     = ladder_ar$hc1_se,     hc1_pval     = ladder_ar$hc1_pval,
                hc1_ci_low = ladder_ar$hc1_ci_low, hc1_ci_high  = ladder_ar$hc1_ci_high,
                hc3_se     = ladder_ar$hc3_se,     hc3_pval     = ladder_ar$hc3_pval,
                hc3_ci_low = ladder_ar$hc3_ci_low, hc3_ci_high  = ladder_ar$hc3_ci_high,
                nw_se      = ladder_ar$nw_se,      nw_pval      = ladder_ar$nw_pval,
                nw_ci_low  = ladder_ar$nw_ci_low,  nw_ci_high   = ladder_ar$nw_ci_high,
                n_obs = sum(!is.na(ar_data[[ols_iv]]) &
                            !is.na(ar_data[[dv]]) &
                            !is.na(ar_data[[dv_lag_col]])),
                has_offensive = has_offensive, has_defensive = has_defensive,
                has_igo = has_igo, has_trade = has_trade, has_dca = has_dca,
                is_onset = FALSE, year_start = year_start, year_end = year_end
              )
            }
          }
        }
      }

      for (war_var in system_war_vars) {
        if (all(is.na(analysis_data[[war_var]]))) next

        # ── VAR / Granger (divergence/degree DVs only) ──
        for (dv in granger_dvs) {
          if (all(is.na(analysis_data[[dv]]))) next

          m_var <- tryCatch({
            var_block <- analysis_data |>
              select(all_of(c(war_var, dv,
                              "n_democracies", "n_autocracies", "n_states"))) |>
              na.omit()
            if (nrow(var_block) < 30) stop("too few obs for VAR")

            var_endo <- var_block |>
              select(all_of(c(war_var, dv))) |>
              as.ts()
            # Harmonization rule 5: VAR exogenous controls are
            # log1p(n_democracies, n_autocracies, n_states).
            # poly(year, 2) is intentionally not added — type = "both" already
            # supplies the constant + linear trend, and the VAR's own lag
            # structure handles longer-run dynamics.
            var_exog <- var_block |>
              select(n_democracies, n_autocracies, n_states) |>
              mutate_all(log1p) |>
              as.ts()

            opt_lag <- tryCatch({
              sel <- VARselect(var_endo, lag.max = 8, type = "both",
                               exogen = var_exog)
              max(1L, as.integer(sel$selection["AIC(n)"]))
            }, error = function(e) 2L)

            vf <- VAR(var_endo, p = opt_lag, type = "both", exogen = var_exog)
            gc <- causality(vf, cause = war_var, boot.runs = 100)
            list(f = as.numeric(gc$Granger$statistic),
                 p = as.numeric(gc$Granger$p.value),
                 lag = opt_lag, vf = vf)
          }, error = function(e) NULL)

          if (!is.null(m_var)) {
            out_list[[length(out_list) + 1]] <- tibble(
              level = "system", dv = dv, war_variable = war_var,
              weighting = method, transform = tfm, layers = combo_name,
              n_layers = length(layers_subset), model = "VAR_Granger",
              coef = m_var$f, se = NA_real_, pval = m_var$p,
              n_obs = nrow(analysis_data),
              has_offensive = has_offensive, has_defensive = has_defensive,
              has_igo = has_igo, has_trade = has_trade, has_dca = has_dca,
              is_onset = FALSE, year_start = year_start, year_end = year_end
            )

            # ── Reverse Granger: DV → war_var ──
            m_var_rev <- tryCatch({
              gc_rev <- causality(m_var$vf, cause = dv, boot.runs = 100)
              list(f = as.numeric(gc_rev$Granger$statistic),
                   p = as.numeric(gc_rev$Granger$p.value))
            }, error = function(e) NULL)
            if (!is.null(m_var_rev)) {
              out_list[[length(out_list) + 1]] <- tibble(
                level = "system", dv = dv, war_variable = war_var,
                weighting = method, transform = tfm, layers = combo_name,
                n_layers = length(layers_subset), model = "VAR_Granger_reverse",
                coef = m_var_rev$f, se = NA_real_, pval = m_var_rev$p,
                n_obs = nrow(analysis_data),
                has_offensive = has_offensive, has_defensive = has_defensive,
                has_igo = has_igo, has_trade = has_trade, has_dca = has_dca,
                is_onset = FALSE, year_start = year_start, year_end = year_end
              )
            }

            # ── FEVD per cell: shares at h ∈ {1, 3, 5, 10} ──
            # Reuses the already-fitted VAR (m_var$vf). Emits 4 rows tagged
            # `model = "VAR_FEVD_h{H}"`, with `coef` set to the share of dv
            # forecast variance attributable to war_var shocks at horizon H.
            # `se` and `pval` are NA — FEVD shares are variance decompositions,
            # not test statistics, so they're excluded from the p-value
            # meta-regression downstream (filtered out before pval_long_meta).
            # They get their own appendix CSV (fevd_multiverse_results.csv)
            # and figure for the H1 magnitude evidence.
            fevd_shares <- tryCatch({
              fd <- vars::fevd(m_var$vf, n.ahead = 10)
              fd_dv <- fd[[dv]]
              c(h1  = fd_dv[1,  war_var],
                h3  = fd_dv[3,  war_var],
                h5  = fd_dv[5,  war_var],
                h10 = fd_dv[10, war_var])
            }, error = function(e) NULL)
            if (!is.null(fevd_shares)) {
              for (hname in names(fevd_shares)) {
                out_list[[length(out_list) + 1]] <- tibble(
                  level = "system", dv = dv, war_variable = war_var,
                  weighting = method, transform = tfm, layers = combo_name,
                  n_layers = length(layers_subset),
                  model = paste0("VAR_FEVD_", hname),
                  coef = unname(fevd_shares[[hname]]),
                  se = NA_real_, pval = NA_real_,
                  n_obs = nrow(analysis_data),
                  has_offensive = has_offensive, has_defensive = has_defensive,
                  has_igo = has_igo, has_trade = has_trade, has_dca = has_dca,
                  is_onset = FALSE, year_start = year_start, year_end = year_end
                )
              }
            }
          }
        }

        # ── LP h=5 (all 8 DVs) ──
        for (dv in system_dvs) {
          if (all(is.na(analysis_data[[dv]]))) next

          lp_data <- analysis_data |>
            arrange(year) |>
            mutate(dv_val = .data[[dv]],
                   dv_lead5 = lead(dv_val, 5),
                   dv_lag1 = lag(dv_val, 1),
                   war_lag1 = lag(.data[[war_var]], 1)) |>
            filter(!is.na(dv_lead5), !is.na(dv_lag1), !is.na(war_lag1))

          m_lp <- tryCatch(
            lm(dv_lead5 ~ war_lag1 + dv_lag1 +
                 log1p(n_democracies) + log1p(n_autocracies) +
                 log1p(n_states) + poly(year, 2),
               data = lp_data),
            error = function(e) NULL
          )
          if (!is.null(m_lp)) {
            sm_lp <- summary(m_lp)$coefficients
            if ("war_lag1" %in% rownames(sm_lp)) {
              # HC3 + Newey-West HAC SEs / 95% CIs for the LP coefficient.
              # The lead-DV setup induces serial correlation in residuals
              # by construction, so HAC matters here even more than in OLS.
              df_resid_lp <- m_lp$df.residual
              t_crit_lp   <- qt(0.975, df = df_resid_lp)
              est_lp      <- sm_lp["war_lag1", "Estimate"]
              se_iid_lp   <- sm_lp["war_lag1", "Std. Error"]

              hc3_vv_lp <- tryCatch(sandwich::vcovHC(m_lp, type = "HC3"),
                                     error = function(e) NULL)
              hc3_se_lp <- if (!is.null(hc3_vv_lp)) sqrt(diag(hc3_vv_lp))["war_lag1"] else NA_real_
              hc3_p_lp  <- if (!is.na(hc3_se_lp)) {
                2 * pt(abs(est_lp / hc3_se_lp), df = df_resid_lp, lower.tail = FALSE)
              } else NA_real_
              hc3_ci_lo_lp <- if (!is.na(hc3_se_lp)) est_lp - t_crit_lp * hc3_se_lp else NA_real_
              hc3_ci_hi_lp <- if (!is.na(hc3_se_lp)) est_lp + t_crit_lp * hc3_se_lp else NA_real_

              nw_vv_lp <- tryCatch(sandwich::NeweyWest(m_lp, lag = 5,
                                                         prewhite = FALSE,
                                                         adjust = TRUE),
                                     error = function(e) NULL)
              nw_se_lp <- if (!is.null(nw_vv_lp)) sqrt(diag(nw_vv_lp))["war_lag1"] else NA_real_
              nw_p_lp  <- if (!is.na(nw_se_lp)) {
                2 * pt(abs(est_lp / nw_se_lp), df = df_resid_lp, lower.tail = FALSE)
              } else NA_real_
              nw_ci_lo_lp <- if (!is.na(nw_se_lp)) est_lp - t_crit_lp * nw_se_lp else NA_real_
              nw_ci_hi_lp <- if (!is.na(nw_se_lp)) est_lp + t_crit_lp * nw_se_lp else NA_real_

              ci_lo_iid_lp <- est_lp - t_crit_lp * se_iid_lp
              ci_hi_iid_lp <- est_lp + t_crit_lp * se_iid_lp

              out_list[[length(out_list) + 1]] <- tibble(
                level = "system", dv = dv, war_variable = war_var,
                weighting = method, transform = tfm, layers = combo_name,
                n_layers = length(layers_subset), model = "LP_h5",
                coef = est_lp,
                se = se_iid_lp,
                pval = sm_lp["war_lag1", "Pr(>|t|)"],
                ci_low_iid = ci_lo_iid_lp, ci_high_iid = ci_hi_iid_lp,
                hc3_se = unname(hc3_se_lp), hc3_pval = hc3_p_lp,
                hc3_ci_low = unname(hc3_ci_lo_lp), hc3_ci_high = unname(hc3_ci_hi_lp),
                nw_se = unname(nw_se_lp), nw_pval = nw_p_lp,
                nw_ci_low = unname(nw_ci_lo_lp), nw_ci_high = unname(nw_ci_hi_lp),
                n_obs = nrow(lp_data),
                has_offensive = has_offensive, has_defensive = has_defensive,
                has_igo = has_igo, has_trade = has_trade, has_dca = has_dca,
                is_onset = FALSE, year_start = year_start, year_end = year_end
              )
            }
          }
        }
      }  # end war_var loop

      result <- dplyr::bind_rows(out_list)
      rm(out_list)
      gc(verbose = FALSE)
      result
    }, .options = furrr::furrr_options(seed = TRUE))  # end combo future_map

    # Persist this batch to disk and FREE from memory. Reloaded later in the
    # post-processing block so the parent process doesn't accumulate
    # 100k+ tibble rows during long runs.
    saveRDS(triple_dfs, ckpt_path)
    n_saved <- nrow(triple_dfs)
    rm(triple_dfs)
    gc(verbose = FALSE)

    if (n_saved > 0) {
      spec_idx <- spec_idx + n_saved
      cat(sprintf("    -> %d specs (running total: %d) [checkpoint saved]\n",
                   n_saved, spec_idx))
    }
  }  # end method loop
}  # end transform loop

n_system_specs <- spec_idx
cat(sprintf("\n  System-level specs: %d\n", n_system_specs))

# Re-compute analysis_data_main for the FEVD horizon plot.
# Previously this was captured as a side-effect inside the inner combo loop
# when (tfm, method, layers) = (raw, unweighted, defensive+igo). Since the
# parallel workers can't write to the parent frame, we re-compute explicitly
# here. Same logic as the inner-loop body, single combo.
{
  layers_main <- c("layer_defensive", "layer_igo")
  manhattan_main <- compute_manhattan_vec(dyad_data, layers_main, "unweighted")
  dv_df_main <- compute_system_dvs(dyad_data, manhattan_main)
  yr_main <- get_valid_year_range(layers_main)
  analysis_data_main <- dv_df_main |>
    left_join(system_war,   by = "year") |>
    left_join(sys_controls, by = "year") |>
    filter(year >= yr_main[1], year <= yr_main[2]) |>
    filter(!is.na(kl_divergence)) |>
    mutate(
      log_kl_divergence      = log1p(kl_divergence),
      log_js_divergence      = log1p(js_divergence),
      log_edge_turnover_pct  = log1p(edge_turnover_pct),
      log_degree_kl          = log1p(degree_kl),
      log_degree_wasserstein = log1p(degree_wasserstein)
    )
}

# =========================================================================== #
# 6. STATE-LEVEL ROBUSTNESS
# =========================================================================== #

cat("\n=== STATE-LEVEL ROBUSTNESS ===\n")

for (tfm in transform_methods) {
  # Apply layer transform (same as system level)
  if (tfm == "log1p") {
    continuous_layers <- intersect(c("layer_igo", "layer_trade"), all_layer_vars)
    dyad_work_state <- dyad_data
    for (cl in continuous_layers) dyad_work_state[[cl]] <- log1p(dyad_work_state[[cl]])
  } else if (tfm == "zscore") {
    dyad_work_state <- dyad_data |>
      group_by(year) |>
      mutate(across(all_of(all_layer_vars), ~ {
        s <- sd(., na.rm = TRUE)
        if (is.na(s) || s == 0) 0 else (. - mean(., na.rm = TRUE)) / s
      })) |>
      ungroup()
  } else {
    dyad_work_state <- dyad_data
  }

  for (method in weighting_methods) {
    if (tfm == "zscore" && grepl("euclidean", method)) next

    # Checkpoint: skip this batch if already computed (count only, no load).
    ckpt_path_state <- file.path(CHECKPOINT_DIR,
                                   sprintf("state__%s__%s.rds", tfm, method))
    if (file.exists(ckpt_path_state)) {
      n_loaded <- nrow(readRDS(ckpt_path_state))
      cat(sprintf("  [state %s/%s] LOADED CHECKPOINT (%d specs)\n",
                  tfm, method, n_loaded))
      spec_idx <- spec_idx + n_loaded
      gc(verbose = FALSE)
      next
    }

    cat(sprintf("  [state %s/%s] dispatching %d combos in parallel (%d workers)...\n",
                tfm, method, length(all_combos), .n_workers))

    triple_dfs_state <- furrr::future_map_dfr(seq_along(all_combos),
                                                function(combo_idx) {
      out_list <- list()
      layers_subset <- all_combos[[combo_idx]]
      combo_name <- paste(gsub("layer_", "", layers_subset), collapse = "+")

      manhattan_vec <- compute_manhattan_vec(dyad_work_state, layers_subset, method)
      state_dvs_df <- compute_state_churn(dyad_work_state, manhattan_vec)

      # Combo-specific network strength at t-1 (for the `strength` control set).
      strength_df <- compute_state_strength(dyad_work_state, manhattan_vec)

    if (nrow(state_dvs_df) < 100) return(dplyr::bind_rows(out_list))

    # Merge with state controls. Canonical churn / log_churn from 02 come
    # along for the ride so the state spec curve has both the per-combo
    # recomputed DV (churn_alt / log_churn_alt) and the canonical 02 DV
    # (churn / log_churn) on the same row.
    #
    # IMPORTANT — combo-specific weighted PWF: state_dvs_df now carries
    # `pwf_weighted_combo` and `pwf_fatal_weighted_combo` (computed with
    # this cell's own edge weights). We OVERWRITE the canonical columns
    # `partner_war_frac_weighted` and `partner_fatal_war_frac_weighted`
    # with these combo-specific values so the existing regression
    # formulas (which reference `partner_war_frac_weighted`) pick up
    # the internally-consistent control. The unweighted variants
    # (`partner_war_frac`, `partner_fatal_war_frac`) don't depend on
    # network weights — they count partners — so they stay canonical.
    churn_merged <- state_dvs_df |>
      left_join(
        state_data |>
          select(ccode, year,
                 churn, log_churn,
                 mid_count, fatal_mid_count,
                 mid_onset_count, fatal_mid_onset_count,
                 log_mid_count, log_fatal_mid_count,
                 log_mid_onset_count, log_fatal_mid_onset_count,
                 partner_war_frac, partner_war_frac_weighted,
                 partner_fatal_war_frac, partner_fatal_war_frac_weighted,
                 is_isolated_t1,
                 democracy, autocracy,
                 cinc, major_power,
                 log_gdp,
                 order_period) |>
          distinct(),
        by = c("ccode", "year")
      ) |>
      # Combo-specific network strength at t-1 (for the `strength` control set).
      left_join(strength_df, by = c("ccode", "year")) |>
      mutate(
        strength_t1     = tidyr::replace_na(strength_t1, 0),
        log_strength_t1 = tidyr::replace_na(log_strength_t1, 0)
      ) |>
      # Overwrite canonical weighted PWF with combo-specific version.
      mutate(
        partner_war_frac_weighted       = pwf_weighted_combo,
        partner_fatal_war_frac_weighted = pwf_fatal_weighted_combo
      )

    # (Spillover variables removed from the multiverse — see top of file.)

    has_offensive <- "layer_offensive" %in% layers_subset
    has_defensive <- "layer_defensive" %in% layers_subset
    has_igo <- "layer_igo" %in% layers_subset
    has_trade <- "layer_trade" %in% layers_subset
    has_dca <- "layer_dca" %in% layers_subset
    yr_range_st <- get_valid_year_range(layers_subset)
    year_start <- yr_range_st[1]
    year_end <- yr_range_st[2]

    # ------------------------------------------------------------------ #
    # TEMPORAL RESTRICTION  [FIXED 2026-07-14]
    #
    # Previously yr_range_st was computed here and then only WRITTEN OUT as
    # metadata (year_start / year_end columns). The state data was never
    # actually filtered by it -- unlike the system loop, which does filter.
    #
    # Consequence of the bug: a layer whose coverage begins late (trade ~1870,
    # DCA 1980) contributed STRUCTURAL ZEROS to the composite tie for every
    # year before its coverage began. Adding zero changes neither the sum nor
    # the within-year max-normalization, so for the pre-coverage era the tie --
    # and therefore churn, the DV -- was IDENTICAL to the combo without that
    # layer. Post-coverage it differed. The DV was thus a different construct
    # before and after the coverage boundary, year FE absorb only the level
    # shift (not a construct change), and roughly 70% of each "robustness" spec
    # was a literal duplicate of the baseline spec.
    #
    # That is why median N was identical (15,309 vs 15,318) with and without
    # trade/DCA, and why those arms produced near-identical significance rates:
    # they were largely the same regression.
    #
    # Fix: mirror the system loop. Restrict to the coverage window of the
    # layer set, so each spec is estimated only where all of its layers exist.
    # ------------------------------------------------------------------ #
    churn_combo <- churn_merged |>
      filter(year >= year_start, year <= year_end)

    for (war_var in state_war_vars) {
      is_onset <- grepl("onset", war_var)

      for (dv in state_dvs) {
        dv_col <- dv

        analysis_sub <- churn_combo |>
          filter(!is.na(.data[[dv_col]]), !is.na(.data[[war_var]]))

        if (nrow(analysis_sub) < 100) next

        # Loop over control sets — IV-routed: severity of partner-war
        # control mirrors the IV (mid IVs use partner_war_frac;
        # fatal IVs use partner_fatal_war_frac). Two weighting variants
        # per IV: unweighted (main) and weighted (appendix robustness).
        iv_control_sets <- get_control_sets(war_var)
        for (ctrl_name in names(iv_control_sets)) {
          controls_rhs <- iv_control_sets[[ctrl_name]]
          exposure_label <- ctrl_name  # "unweighted" or "weighted"

        # ── A. TWFE full × clustering schemes × EXTRA CONTROL SETS ──
        # extra_ctrl in {none, gdp, strength}. Applied to FE_full ONLY --
        # see the EXTRA_CTRLS note above. `strength` drops is_isolated_t1.
        for (extra_ctrl in EXTRA_CTRLS) {

          # GDP/strength run ONLY on the theory-motivated IV x DV subset
          # [2026-07-15]. The baseline (extra_ctrl == "none") keeps the full
          # sweep. This is what makes the GDP/strength robustness cheap.
          if (extra_ctrl != "none" &&
              (!war_var %in% GDP_STRENGTH_IVS || !dv %in% GDP_STRENGTH_DVS)) next

          rhs_x <- augment_controls(controls_rhs, extra_ctrl)

          # Skip if the augmenting variable is absent or degenerate in this cell.
          if (extra_ctrl == "gdp" &&
              (!"log_gdp" %in% names(analysis_sub) ||
               all(is.na(analysis_sub$log_gdp)))) next
          if (extra_ctrl == "strength" &&
              (!"log_strength_t1" %in% names(analysis_sub) ||
               sd(analysis_sub$log_strength_t1, na.rm = TRUE) == 0)) next

          fml_fe <- as.formula(paste(dv_col, "~", war_var, "+",
                                     rhs_x, "| ccode + year"))

          for (clust_name in names(clustering_schemes)) {
            sub_data <- analysis_sub

            m_fe <- tryCatch(
              feols(fml_fe, data = sub_data,
                    vcov = clustering_schemes[[clust_name]]),
              error = function(e) NULL
            )
            if (!is.null(m_fe)) {
              ct <- coeftable(m_fe)
              if (war_var %in% rownames(ct)) {
                out_list[[length(out_list) + 1]] <- tibble(
                  level = "state", dv = dv, war_variable = war_var,
                  weighting = method, transform = tfm, layers = combo_name,
                  n_layers = length(layers_subset),
                  model = paste0("FE_full_", clust_name),
                  extra_control = extra_ctrl,
                  exposure = exposure_label,
                  coef = ct[war_var, "Estimate"],
                  se = ct[war_var, "Std. Error"],
                  pval = ct[war_var, "Pr(>|t|)"],
                  n_obs = nobs(m_fe),
                  has_offensive = has_offensive, has_defensive = has_defensive,
                  has_igo = has_igo, has_trade = has_trade, has_dca = has_dca,
                  is_onset = is_onset,
                  year_start = year_start, year_end = year_end
                )
              }
            }
          }
        }  # end extra_ctrl loop

        # ── C. Full model, no FE ──
        fml_full_nofe <- as.formula(paste(dv_col, "~", war_var, "+",
                                          controls_rhs))

        m_full_nofe <- tryCatch(
          lm(fml_full_nofe, data = analysis_sub),
          error = function(e) NULL
        )
        if (!is.null(m_full_nofe)) {
          sm <- summary(m_full_nofe)$coefficients
          if (war_var %in% rownames(sm)) {
            out_list[[length(out_list) + 1]] <- tibble(
              level = "state", dv = dv, war_variable = war_var,
              weighting = method, transform = tfm, layers = combo_name,
              n_layers = length(layers_subset), model = "full_noFE",
              exposure = exposure_label,
              coef = sm[war_var, "Estimate"],
              se = sm[war_var, "Std. Error"],
              pval = sm[war_var, "Pr(>|t|)"],
              n_obs = nrow(analysis_sub),
              has_offensive = has_offensive, has_defensive = has_defensive,
              has_igo = has_igo, has_trade = has_trade, has_dca = has_dca,
              is_onset = is_onset,
              year_start = year_start, year_end = year_end
            )
          }
        }

        # ── D. Full model, no FE, with poly(year, 2) trend ──
        # Mirror of (C) but adds a quadratic year trend in lieu of state /
        # year fixed effects. Coefficients and p-values pulled with the same
        # SE ladder (iid, HC3, NW) as the system-level OLS specs above so
        # they appear cleanly in the spec curve and meta-regression.
        fml_full_nofe_y <- as.formula(paste(dv_col, "~", war_var, "+",
                                              controls_rhs, "+ poly(year, 2)"))

        m_full_nofe_y <- tryCatch(
          lm(fml_full_nofe_y, data = analysis_sub),
          error = function(e) NULL
        )
        if (!is.null(m_full_nofe_y)) {
          sm_y <- summary(m_full_nofe_y)$coefficients
          if (war_var %in% rownames(sm_y)) {
            df_resid <- m_full_nofe_y$df.residual
            t_crit   <- qt(0.975, df = df_resid)
            est_val  <- sm_y[war_var, "Estimate"]
            se_iid   <- sm_y[war_var, "Std. Error"]

            hc3_vv <- tryCatch(sandwich::vcovHC(m_full_nofe_y, type = "HC3"),
                                error = function(e) NULL)
            hc3_se_val <- if (!is.null(hc3_vv)) sqrt(diag(hc3_vv))[war_var] else NA_real_
            hc3_p_val  <- if (!is.na(hc3_se_val)) {
              2 * pt(abs(est_val / hc3_se_val), df = df_resid, lower.tail = FALSE)
            } else NA_real_
            hc3_ci_lo  <- if (!is.na(hc3_se_val)) est_val - t_crit * hc3_se_val else NA_real_
            hc3_ci_hi  <- if (!is.na(hc3_se_val)) est_val + t_crit * hc3_se_val else NA_real_

            nw_vv <- tryCatch(sandwich::NeweyWest(m_full_nofe_y, lag = 5,
                                                    prewhite = FALSE,
                                                    adjust = TRUE),
                                error = function(e) NULL)
            nw_se_val <- if (!is.null(nw_vv)) sqrt(diag(nw_vv))[war_var] else NA_real_
            nw_p_val  <- if (!is.na(nw_se_val)) {
              2 * pt(abs(est_val / nw_se_val), df = df_resid, lower.tail = FALSE)
            } else NA_real_
            nw_ci_lo  <- if (!is.na(nw_se_val)) est_val - t_crit * nw_se_val else NA_real_
            nw_ci_hi  <- if (!is.na(nw_se_val)) est_val + t_crit * nw_se_val else NA_real_

            out_list[[length(out_list) + 1]] <- tibble(
              level = "state", dv = dv, war_variable = war_var,
              weighting = method, transform = tfm, layers = combo_name,
              n_layers = length(layers_subset), model = "full_noFE_year_poly",
              exposure = exposure_label,
              coef = est_val,
              se = se_iid,
              pval = sm_y[war_var, "Pr(>|t|)"],
              ci_low_iid  = est_val - t_crit * se_iid,
              ci_high_iid = est_val + t_crit * se_iid,
              hc3_se = unname(hc3_se_val), hc3_pval = hc3_p_val,
              hc3_ci_low = unname(hc3_ci_lo), hc3_ci_high = unname(hc3_ci_hi),
              nw_se = unname(nw_se_val), nw_pval = nw_p_val,
              nw_ci_low = unname(nw_ci_lo), nw_ci_high = unname(nw_ci_hi),
              n_obs = nrow(analysis_sub),
              has_offensive = has_offensive, has_defensive = has_defensive,
              has_igo = has_igo, has_trade = has_trade, has_dca = has_dca,
              is_onset = is_onset,
              year_start = year_start, year_end = year_end
            )
          }
        }

        # ── D2. State + order_period FE  [REMOVED 2026-07-14] ──
        #
        # This model class clustered on `order_period`, which has FOUR levels.
        # fixest takes the MINIMUM cluster count for the degrees of freedom, so
        # every one of these 12,480 specs was tested against a t-distribution
        # with 3 df (critical value ~3.18 rather than 1.96). Cluster-robust
        # inference with G = 4 is not conservative, it is invalid
        # (Cameron & Miller 2015).
        #
        # The diagnostic: matched against FE_full_ccode_year on identical spec
        # cells (n = 12,480), the point estimates were TWICE AS LARGE
        # (median 0.0048 vs 0.0025) and agreed in sign 94.4% of the time, but
        # the SEs were 2.4x larger. The class reported only 29% significant-
        # positive -- not because the effect failed under coarse time FE, but
        # because the standard errors were computed off four clusters.
        #
        # Reporting it would have made the result look FRAGILE when the
        # correctly-estimated version is STRONGER than the headline model.
        # Dropped rather than re-clustered: the year-FE specs already cover
        # the "does it survive time controls" question.
        if (FALSE) {
        fml_state_period <- as.formula(paste(dv_col, "~", war_var, "+",
                                               controls_rhs,
                                               "| ccode + order_period"))
        m_state_period <- tryCatch(
          feols(fml_state_period, data = analysis_sub,
                vcov = ~ ccode + order_period),
          error = function(e) NULL
        )
        if (!is.null(m_state_period)) {
          ct_sp <- coeftable(m_state_period)
          if (war_var %in% rownames(ct_sp)) {
            out_list[[length(out_list) + 1]] <- tibble(
              level = "state", dv = dv, war_variable = war_var,
              weighting = method, transform = tfm, layers = combo_name,
              n_layers = length(layers_subset),
              model = "state_period_fe",
              exposure = exposure_label,
              coef = ct_sp[war_var, "Estimate"],
              se = ct_sp[war_var, "Std. Error"],
              pval = ct_sp[war_var, "Pr(>|t|)"],
              n_obs = nrow(analysis_sub),
              has_offensive = has_offensive, has_defensive = has_defensive,
              has_igo = has_igo, has_trade = has_trade, has_dca = has_dca,
              is_onset = is_onset,
              year_start = year_start, year_end = year_end
            )
          }
        }
        }  # end if (FALSE) -- state_period_fe disabled, see note above

        # (E. Spillover specs removed — spillover analysis is now confined
        # to 21_spillover_nonfighters.R; see notes at top of file.)

        }  # end control_sets loop

        # ── B. Baseline (no FE, no controls) — run once, no exposure variant.
        # Includes a quadratic year trend (poly(year, 2)) to mirror the
        # main-text Table 3 cols 1-2 spec ladder. ──
        fml_biv <- as.formula(paste(dv_col, "~", war_var, "+ poly(year, 2)"))

        m_biv <- tryCatch(
          lm(fml_biv, data = analysis_sub),
          error = function(e) NULL
        )
        if (!is.null(m_biv)) {
          sm <- summary(m_biv)$coefficients
          if (war_var %in% rownames(sm)) {
            out_list[[length(out_list) + 1]] <- tibble(
              level = "state", dv = dv, war_variable = war_var,
              weighting = method, transform = tfm, layers = combo_name,
              n_layers = length(layers_subset), model = "baseline_noFE",
              exposure = "none",
              coef = sm[war_var, "Estimate"],
              se = sm[war_var, "Std. Error"],
              pval = sm[war_var, "Pr(>|t|)"],
              n_obs = nrow(analysis_sub),
              has_offensive = has_offensive, has_defensive = has_defensive,
              has_igo = has_igo, has_trade = has_trade, has_dca = has_dca,
              is_onset = is_onset,
              year_start = year_start, year_end = year_end
            )
          }
        }

      }  # end dv loop
    }  # end war_var loop

      result <- dplyr::bind_rows(out_list)
      rm(out_list)
      gc(verbose = FALSE)
      result
    }, .options = furrr::furrr_options(seed = TRUE))  # end combo future_map (state)

    saveRDS(triple_dfs_state, ckpt_path_state)
    n_saved <- nrow(triple_dfs_state)
    rm(triple_dfs_state)
    gc(verbose = FALSE)

    if (n_saved > 0) {
      spec_idx <- spec_idx + n_saved
      cat(sprintf("    -> %d state specs (running total: %d) [checkpoint saved]\n",
                   n_saved, spec_idx))
    }
  }  # end method loop
}  # end transform loop (state)

n_state_specs <- spec_idx - n_system_specs
cat(sprintf("\n  State-level specs: %d\n", n_state_specs))

# =========================================================================== #
# 7. Combine and Summarize
# =========================================================================== #
# Reload all per-batch checkpoints from disk and bind into a single
# tibble. The parallel loops above intentionally don't accumulate
# spec_results in the parent process to bound memory; we reconstruct
# the full multiverse here from the per-(tfm, method) RDS files.

cat("\n  Reloading checkpoints from disk for aggregation...\n")
ckpt_files <- list.files(CHECKPOINT_DIR, pattern = "\\.rds$", full.names = TRUE)
spec_results <- lapply(ckpt_files, readRDS)
specs_df <- bind_rows(spec_results)
rm(spec_results)
gc(verbose = FALSE)
cat(sprintf("  Loaded %d checkpoint files, %d total spec rows.\n",
              length(ckpt_files), nrow(specs_df)))

# `extra_control` only exists on FE_full rows (the GDP / strength dimension).
# Every other row -- system level, no-FE classes, bivariate -- is the baseline
# control set. bind_rows fills the column with NA; label it explicitly so
# downstream filters and facets can key on it without NA handling.
if (!"extra_control" %in% names(specs_df)) specs_df$extra_control <- "none"
specs_df <- specs_df |>
  dplyr::mutate(extra_control = tidyr::replace_na(extra_control, "none"))

cat("\n  extra_control breakdown:\n")
print(specs_df |> dplyr::count(extra_control))

# Split FEVD-share rows out of specs_df. They have NA pval / se by
# construction (variance decomposition, not test statistic) and would
# pollute the meta-regression on -log10(p). Persist them to their own
# CSV for the appendix figure / table downstream.
# Strict pattern match: only keep rows whose model is exactly
# "VAR_FEVD_h{N}" with N integer. The earlier looser regex pulled in any
# row starting "VAR_FEVD_" and threw a coercion warning when the suffix
# wasn't a number; the tightened version side-steps that.
fevd_pattern <- "^VAR_FEVD_h([0-9]+)$"
fevd_specs_df <- specs_df |>
  dplyr::filter(grepl(fevd_pattern, model)) |>
  dplyr::mutate(horizon = as.integer(sub(fevd_pattern, "\\1", model)))
specs_df <- specs_df |> dplyr::filter(!grepl(fevd_pattern, model))

if (nrow(fevd_specs_df) > 0) {
  readr::write_csv(fevd_specs_df,
                    file.path(OUTPUT_DIR, "fevd_multiverse_results.csv"))
  cat(sprintf("  FEVD multiverse: %d rows -> %s\n",
              nrow(fevd_specs_df),
              file.path("data", "processed",
                          "fevd_multiverse_results.csv")))
}

cat(sprintf("\n=== TOTAL SPECIFICATIONS: %d ===\n", nrow(specs_df)))

cat(sprintf("  System: %d\n", sum(specs_df$level == "system")))
# One-sided positive significance threshold used everywhere downstream.
# Spec-curve / forest / tile / console summaries highlight a spec only if
# its coefficient is in the predicted direction (positive) AND the one-sided
# p-value clears the threshold. Equivalently: pval / 2 < ALPHA_ONESIDE & coef > 0.
ALPHA_ONESIDE <- 0.10

is_sig_pos <- function(coef, p) coef > 0 & p / 2 < ALPHA_ONESIDE

cat(sprintf("  State: %d\n", sum(specs_df$level == "state")))
cat(sprintf("  Significant positive (one-sided, alpha=%.2f): %d / %d (%.1f%%)\n",
            ALPHA_ONESIDE,
            sum(is_sig_pos(specs_df$coef, specs_df$pval), na.rm = TRUE),
            nrow(specs_df),
            100 * mean(is_sig_pos(specs_df$coef, specs_df$pval), na.rm = TRUE)))

cat("\nBy level and DV:\n")
specs_df |>
  group_by(level, dv) |>
  summarise(n = n(),
            pct_sig_pos = 100 * mean(is_sig_pos(coef, pval), na.rm = TRUE),
            median_pval = median(pval, na.rm = TRUE), .groups = "drop") |>
  print()

cat("\nBy war variable:\n")
specs_df |>
  group_by(war_variable) |>
  summarise(n = n(),
            pct_sig_pos = 100 * mean(is_sig_pos(coef, pval), na.rm = TRUE),
            .groups = "drop") |>
  print()

cat("\nBy model type:\n")
specs_df |>
  group_by(model) |>
  summarise(n = n(),
            pct_sig_pos = 100 * mean(is_sig_pos(coef, pval), na.rm = TRUE),
            .groups = "drop") |>
  print()

cat("\nBy weighting:\n")
specs_df |>
  group_by(weighting) |>
  summarise(n = n(),
            pct_sig_pos = 100 * mean(is_sig_pos(coef, pval), na.rm = TRUE),
            .groups = "drop") |>
  print()

# =========================================================================== #
# 8. Meta-Analysis
# =========================================================================== #

cat("\n=== Meta-Analysis ===\n")

# Derive perturbation indicators that decompose the spec-curve mass.
# These join the layer / weighting / transform / model / level / DV / exposure
# regressors already in the original meta-regression.
#
#   dv_logged           — DV is log1p-transformed (vs. raw scale)
#   model_no_fe         — spec uses no fixed effects (OLS / pooled)
#   model_year_poly     — spec includes a poly(year, 2) trend either alongside
#                          FE or in lieu of FE
#   state_dv_canonical  — for state specs, DV is the canonical 02-output churn
#                          (vs. the per-combo recomputed churn_alt)
#
# These cross-cut the existing `model` and `dv` factors but make the
# perturbations explicit so the meta-regression coefficient on, say,
# `dv_logged` directly tells you "log1p the DV, holding everything else
# constant — what does that do to the estimate?"
specs_df <- specs_df |>
  mutate(
    dv_logged          = grepl("^log_", dv),
    model_no_fe        = model %in% c("OLS_lagged", "LP_h5",
                                         "full_noFE", "full_noFE_year_poly",
                                         "baseline_noFE"),
    model_year_poly    = model %in% c("OLS_lagged", "LP_h5",
                                         "full_noFE_year_poly"),
    state_dv_canonical = level == "state" & dv %in% c("churn", "log_churn"),
    # ── Explicit specification-choice indicators (added so the
    # meta-regression spells out each design decision separately
    # rather than burying it in a multi-level factor). ──
    iv_is_fatal           = grepl("fatal", war_variable),
    iv_is_log             = grepl("^log_", war_variable),
    partner_ctrl_weighted = exposure == "weighted",
    # The IV-routed control mirrors IV severity by construction in this
    # multiverse, so partner_ctrl_severity_matched is constant TRUE and
    # is recorded for documentation but excluded from the regression.
    partner_ctrl_severity_matched = TRUE
  )

# ── Benjamini-Hochberg FDR adjustment, applied within each (level, dv)
# cell, separately for each p-value flavor. Relocated to here from below
# so the FDR columns are available to every downstream consumer (spec
# curve plot, summary tables, meta-regs, etc.). The standard spec-curve
# adjustment for dependent multiple tests within a single DV's multiverse.
specs_df <- specs_df |>
  group_by(level, dv) |>
  mutate(
    pval_fdr_iid = if (sum(!is.na(pval))     > 0) p.adjust(pval,     method = "BH") else NA_real_,
    pval_fdr_hc1 = if ("hc1_pval" %in% names(cur_data()) && sum(!is.na(hc1_pval)) > 0)
                     p.adjust(hc1_pval, method = "BH") else NA_real_,
    pval_fdr_hc3 = if ("hc3_pval" %in% names(cur_data()) && sum(!is.na(hc3_pval)) > 0)
                     p.adjust(hc3_pval, method = "BH") else NA_real_,
    pval_fdr_nw  = if ("nw_pval"  %in% names(cur_data()) && sum(!is.na(nw_pval))  > 0)
                     p.adjust(nw_pval,  method = "BH") else NA_real_
  ) |>
  ungroup()

specs_for_coef <- specs_df |> filter(!model %in% c("VAR_Granger",
                                                     "VAR_Granger_reverse"))

# Helper: drop RHS variables with fewer than 2 unique non-NA values.
# Required because system and state subsets each have a different set of
# variables that are constant within them (e.g. `level` is constant within
# each subset; `state_dv_canonical` is always FALSE on system rows;
# `exposure` is NA on system rows under the original schema).
drop_singletons <- function(df, rhs_vars) {
  Filter(function(v) {
    if (!v %in% names(df)) return(FALSE)
    length(unique(stats::na.omit(df[[v]]))) >= 2
  }, rhs_vars)
}

# Helper: prep a data frame for the meta-regression. Three transforms:
#
#   (1) Convert logical columns to 0/1 integers. Under the no-intercept
#       `~ 0 + ...` formula, R would otherwise encode each logical with
#       two dummies (e.g. has_offensiveTRUE *and* has_offensiveFALSE),
#       bloating the coefficient table with redundant rows. Converting
#       to integer gives one slope per binary while still letting
#       multi-level factors keep all their levels.
#
#   (2) Drop unused factor levels. When a subset contains zero rows for
#       a particular level, R still tries to fit a coefficient for that
#       level and returns NA. Calling droplevels() removes empty levels
#       so only populated cells get coefficients.
#
#   (3) Set canonical baseline reference levels for the multi-level
#       factors so coefficients read as effects relative to the
#       manuscript's headline specifications:
#         - weighting    : "unweighted"  (headline weighting scheme)
#         - transform    : "raw"         (headline IV transform)
#         - model        : "FE_full_ccode_year" (headline state TWFE) for
#                            state subsets, "OLS_lagged" for system.
#         - exposure     : "weighted"    (headline state PWE control)
binarize_logicals <- function(df) {
  df <- df |>
    dplyr::mutate(dplyr::across(dplyr::where(is.logical), as.integer)) |>
    dplyr::mutate(dplyr::across(dplyr::where(is.factor), forcats::fct_drop)) |>
    droplevels()

  # Promote relevant character columns to factor with explicit references
  # so summary(lm(...)) coefficients read as deviations from the headline.
  set_ref <- function(df, col, ref) {
    if (col %in% names(df)) {
      vec <- df[[col]]
      lv  <- if (is.factor(vec)) levels(vec) else unique(stats::na.omit(vec))
      if (ref %in% lv) {
        df[[col]] <- stats::relevel(factor(vec), ref = ref)
      } else if (length(lv) >= 2) {
        df[[col]] <- factor(vec)  # promote without releveling if ref absent
      }
    }
    df
  }

  df <- set_ref(df, "weighting", "unweighted")
  df <- set_ref(df, "transform", "raw")
  df <- set_ref(df, "exposure",  "weighted")
  # model reference depends on level — pick FE_full_ccode_year if the
  # subset has it (state), otherwise OLS_lagged (system).
  if ("model" %in% names(df)) {
    mvals <- if (is.factor(df$model)) levels(df$model) else unique(stats::na.omit(df$model))
    if ("FE_full_ccode_year" %in% mvals) {
      df <- set_ref(df, "model", "FE_full_ccode_year")
    } else if ("OLS_lagged" %in% mvals) {
      df <- set_ref(df, "model", "OLS_lagged")
    }
  }
  df
}

# Long-form pval frame builder (mirrors the original pval_long_meta logic
# but on a passed-in subset).
build_pval_long <- function(df) {
  df |>
    dplyr::select(coef, pval,
                  dplyr::any_of(c("hc1_pval", "hc3_pval", "nw_pval")),
                  has_offensive, has_defensive, has_igo,
                  has_trade, has_dca,
                  is_onset, iv_is_fatal, iv_is_log,
                  weighting, transform, model, war_variable,
                  level, dv, exposure, partner_ctrl_weighted, year_start,
                  dv_logged, model_no_fe, model_year_poly,
                  state_dv_canonical) |>
    tidyr::pivot_longer(
      cols = dplyr::any_of(c("pval", "hc1_pval", "hc3_pval", "nw_pval")),
      names_to  = "se_method",
      values_to = "p_value") |>
    dplyr::filter(!is.na(p_value)) |>
    dplyr::mutate(
      se_method = factor(se_method,
                         levels = c("pval", "hc1_pval", "hc3_pval", "nw_pval"),
                         labels = c("iid", "HC1", "HC3", "NW(5)")),
      neg_log_p = -log10(pmax(p_value, 1e-16)),
      sig_05    = as.integer(p_value < 0.05))
}

# Helper: run the trio (coef LM, -log10(p) LM, sig_05 logit) on one subset
# and persist coefficient tables to disk. Drop-singleton RHS handling makes
# the helper robust to whatever subset is passed in.
run_meta_split <- function(label, coef_df, pval_df,
                            rhs_vars_full, out_dir, out_tag) {
  cat(sprintf("\n=== Meta-Analysis: %s (n_coef=%d, n_pval_long=%d) ===\n",
              label, nrow(coef_df), nrow(pval_df)))

  # Convert logicals to 0/1 so they don't double-encode under ~ 0 +
  coef_df <- binarize_logicals(coef_df)
  pval_df <- binarize_logicals(pval_df)

  rhs_coef <- drop_singletons(coef_df, rhs_vars_full)
  m_coef   <- NULL
  if (length(rhs_coef) > 0 && nrow(coef_df) > 0) {
    # Standard intercept-coded formula. Avoids rank deficiency that
    # `~ 0 + factor1 + factor2 + ...` introduces when multiple factors'
    # level-dummies all sum to 1 (causing R to alias out one level per
    # factor as NA — e.g. exposureself, modelbaseline_noFE).
    fml_coef <- as.formula(paste("coef ~",
                                  paste(rhs_coef, collapse = " + ")))
    m_coef <- tryCatch(lm(fml_coef, data = coef_df),
                       error = function(e) {
                         cat("  coef LM failed:", conditionMessage(e), "\n")
                         NULL
                       })
    if (!is.null(m_coef)) { cat("  --- coef as DV ---\n"); print(summary(m_coef)) }
  }

  rhs_pval <- drop_singletons(pval_df, c(rhs_vars_full, "se_method"))
  m_pval   <- NULL
  if (length(rhs_pval) > 0 && nrow(pval_df) > 0) {
    fml_pval <- as.formula(paste("neg_log_p ~",
                                  paste(rhs_pval, collapse = " + ")))
    m_pval <- tryCatch(lm(fml_pval, data = pval_df),
                       error = function(e) {
                         cat("  pval LM failed:", conditionMessage(e), "\n")
                         NULL
                       })
    if (!is.null(m_pval)) { cat("  --- -log10(p) as DV ---\n"); print(summary(m_pval)) }
  }

  m_sig <- NULL
  if (length(rhs_pval) > 0 && nrow(pval_df) > 0) {
    fml_sig <- as.formula(paste("sig_05 ~",
                                 paste(rhs_pval, collapse = " + ")))
    m_sig <- tryCatch(glm(fml_sig, data = pval_df, family = binomial),
                      error = function(e) {
                        cat("  sig logit failed:", conditionMessage(e), "\n")
                        NULL
                      })
    if (!is.null(m_sig)) { cat("  --- sig_05 logit ---\n"); print(summary(m_sig)) }
  }

  for (mm in list(list(nm = "coef",     fit = m_coef),
                  list(nm = "negLogP",  fit = m_pval),
                  list(nm = "sigLogit", fit = m_sig))) {
    if (is.null(mm$fit)) next
    sm <- tryCatch(summary(mm$fit)$coefficients, error = function(e) NULL)
    if (is.null(sm)) next
    tib <- tibble::as_tibble(sm, rownames = "term")
    nms <- c("term", "estimate", "se", "tstat", "pval")
    names(tib) <- nms[seq_len(ncol(tib))]
    readr::write_csv(tib,
                     file.path(out_dir,
                                sprintf("meta_%s_%s.csv", out_tag, mm$nm)))
  }

  invisible(list(coef = m_coef, pval = m_pval, sig = m_sig))
}

# RHS menu — keeps the substantive multi-level factors and independent
# layer indicators; OMITS:
#   - Derived binary indicators (is_onset, iv_is_fatal, iv_is_log,
#     dv_logged, state_dv_canonical, model_no_fe, model_year_poly,
#     partner_ctrl_weighted) because they are perfect linear functions
#     of factors already in the model.
#   - `exposure`, because exposurenone is a perfect alias for
#     modelbaseline_noFE, exposureself is captured by
#     (war_variable ∈ PWF) × (model ≠ baseline), and exposureweighted
#     vs exposureunweighted is the only independent comparison (handled
#     in a focused side regression in §A8 if needed).
# The factor-level interpretation is substantively richer — e.g.,
# `modelFE_full_ccode_year` tells you what the headline TWFE does, while
# `model_no_feTRUE` would collapse across multiple FE schemes.
# Singletons get dropped per subset by drop_singletons().
RHS_VARS_FULL <- c("has_offensive", "has_defensive", "has_igo",
                   "has_trade", "has_dca",
                   "weighting", "transform", "model", "war_variable",
                   "dv", "year_start")

OMNI_OUT_DIR <- APPX_DIR
dir.create(OMNI_OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# Split coef df by level
coef_sys <- specs_for_coef |> dplyr::filter(level == "system")
coef_st  <- specs_for_coef |> dplyr::filter(level == "state")

# Long-form pval frames per level
pval_sys <- build_pval_long(coef_sys)
pval_st  <- build_pval_long(coef_st)

# Run separate omnibus regressions per level
m_meta_sys <- run_meta_split(
  label   = "Omnibus (system level)",
  coef_df = coef_sys,
  pval_df = pval_sys,
  rhs_vars_full = RHS_VARS_FULL,
  out_dir = OMNI_OUT_DIR,
  out_tag = "omnibus_system"
)
m_meta_st <- run_meta_split(
  label   = "Omnibus (state level)",
  coef_df = coef_st,
  pval_df = pval_st,
  rhs_vars_full = RHS_VARS_FULL,
  out_dir = OMNI_OUT_DIR,
  out_tag = "omnibus_state"
)

# ── Backwards-compatibility shims for downstream blocks ──
# Some downstream code (extract_meta_coefs, conditional-effect summary,
# p-value distribution plot, texreg combined output) expects
# `m_meta_coef`, `m_meta_pval`, `m_meta_sig`, and `pval_long_meta` to
# exist as single objects. Keep those names alive by aliasing to whichever
# subset is present (preferring state since it's larger and carries the
# partner-war-exposure factor). The split-level CSVs in
# results/appendix/meta_omnibus_{system,state}_*.csv are the new
# authoritative omnibus outputs.
or_null <- function(x, y) if (!is.null(x)) x else y
m_meta_coef    <- or_null(m_meta_st$coef, m_meta_sys$coef)
m_meta_pval    <- or_null(m_meta_st$pval, m_meta_sys$pval)
m_meta_sig     <- or_null(m_meta_st$sig,  m_meta_sys$sig)
pval_long_meta <- dplyr::bind_rows(pval_sys, pval_st)

# =========================================================================== #
# 8b. H1 / H2 / H3 filtered meta-regressions
# =========================================================================== #
# Three filtered meta-regressions tied directly to the manuscript's three
# hypotheses. Each runs the same trio (coef LM, -log10(p) LM, sig logit)
# but on a hypothesis-specific slice of the multiverse:
#
#   H1: war reshapes the SYSTEM network
#       level == "system", IV ∈ {MID count family}, model ∈ {sys_T1_*,
#       OLS_lagged, VAR_Granger}.
#   H2: war reshapes a state's OWN connections
#       level == "state", IV ∈ {own-MID onset/count family}, model ∈ Table 3
#       column classes (baseline, full_noFE_year_poly, FE_full_*).
#   H3: war reshapes a PARTNER's network (1-hop spillover)
#       level == "state", IV ∈ {partner_war_frac variants}, same model set
#       as H2.
#
# Each H-meta-reg writes its own coef table to results/appendix/.

OWN_MID_IVS <- c("mid_onset_count", "log_mid_onset_count",
                 "fatal_mid_onset_count", "log_fatal_mid_onset_count",
                 "mid_count", "log_mid_count",
                 "fatal_mid_count", "log_fatal_mid_count")

PWF_IVS <- c("partner_war_frac", "partner_war_frac_weighted",
             "partner_fatal_war_frac", "partner_fatal_war_frac_weighted")

H1_MODELS <- c("sys_T1_iv_only", "OLS_lagged",
               "sys_T1_ar1", "sys_T1_order_period",
               "VAR_Granger")

# State-level model classes that mirror Table 3 columns (Leiden-clustered
# SE class names removed; clustering schemes are now ccode+year, year,
# ccode only):
#   col 1 (IV alone + poly(year, 2)):     baseline_noFE
#   col 2 (PWE alone + poly(year, 2)):    baseline_noFE filtered to PWF IVs (H3)
#   col 3 (+ ctrls, no FE):               full_noFE / full_noFE_year_poly
#   col 4 (TWFE):                         FE_full_ccode_year
H_STATE_MODELS <- c("baseline_noFE",
                    "full_noFE", "full_noFE_year_poly",
                    "FE_full_ccode_year", "FE_full_year",
                    "FE_full_year_only", "FE_full_ccode_only",
                    "state_period_fe")

run_h_meta <- function(label, coef_df, pval_df, h_dir, h_tag) {
  cat(sprintf("\n=== H-meta: %s (n_coef=%d, n_pval_long=%d) ===\n",
              label, nrow(coef_df), nrow(pval_df)))

  # Drop factor levels with only one value to avoid singular fits.
  drop_singletons <- function(df, vars) {
    keep <- vars
    for (v in vars) {
      if (v %in% names(df) && dplyr::n_distinct(df[[v]]) <= 1) {
        keep <- setdiff(keep, v)
      }
    }
    keep
  }

  # Reuse the omnibus prep helper: logical→integer, drop empty levels,
  # set canonical baseline reference levels (weighting="unweighted",
  # transform="raw", exposure="weighted", model=headline TWFE/OLS).
  coef_df <- binarize_logicals(coef_df)
  pval_df <- binarize_logicals(pval_df)

  # See note on RHS_VARS_FULL above: derived binary indicators and
  # `exposure` factor omitted because they alias with parent factors
  # already in the model and force NA coefficients.
  rhs_full <- c("has_offensive", "has_defensive", "has_igo",
                "has_trade", "has_dca",
                "weighting", "transform", "model", "war_variable",
                "dv", "year_start")

  rhs_coef <- drop_singletons(coef_df, rhs_full)
  # Intercept-coded formula: reference level for each factor is absorbed
  # into the intercept; binaries (post binarize) get one slope each.
  fml_coef <- as.formula(paste("coef ~",
                                 paste(rhs_coef, collapse = " + ")))

  m_coef <- tryCatch(lm(fml_coef, data = coef_df),
                     error = function(e) {
                       cat("  coef LM failed:", conditionMessage(e), "\n")
                       NULL
                     })
  if (!is.null(m_coef)) {
    cat("  --- coef LM ---\n")
    print(summary(m_coef))
  }

  rhs_pval <- drop_singletons(pval_df, c(rhs_full, "se_method"))
  fml_pval <- as.formula(paste("neg_log_p ~",
                                 paste(rhs_pval, collapse = " + ")))
  m_pval <- tryCatch(lm(fml_pval, data = pval_df),
                     error = function(e) {
                       cat("  pval LM failed:", conditionMessage(e), "\n")
                       NULL
                     })
  if (!is.null(m_pval)) {
    cat("  --- -log10(p) LM ---\n")
    print(summary(m_pval))
  }

  fml_sig <- as.formula(paste("sig_05 ~",
                                paste(rhs_pval, collapse = " + ")))
  m_sig <- tryCatch(glm(fml_sig, data = pval_df, family = binomial),
                    error = function(e) {
                      cat("  sig logit failed:", conditionMessage(e), "\n")
                      NULL
                    })
  if (!is.null(m_sig)) {
    cat("  --- sig_05 logit ---\n")
    print(summary(m_sig))
  }

  # Persist coef tables for the appendix.
  for (mm in list(list(name = "coef", fit = m_coef),
                  list(name = "negLogP", fit = m_pval),
                  list(name = "sigLogit", fit = m_sig))) {
    if (is.null(mm$fit)) next
    sm <- tryCatch(summary(mm$fit)$coefficients, error = function(e) NULL)
    if (is.null(sm)) next
    tib <- tibble::as_tibble(sm, rownames = "term")
    names(tib) <- c("term", "estimate", "se", "tstat", "pval")[
      seq_len(ncol(tib))]
    readr::write_csv(
      tib,
      file.path(h_dir, sprintf("meta_%s_%s.csv", h_tag, mm$name))
    )
  }

  invisible(list(coef = m_coef, p = m_pval, sig = m_sig))
}

H_DIR <- file.path(APPX_DIR)

# H1 ------------------------------------------------------------
# NOTE: VAR_Granger F-statistics are excluded from the coef-as-DV regression
# (different statistic semantics from OLS coefficients) and are also absent
# from pval_long_meta because that long form is derived from specs_for_coef
# (which filters them out). The Granger evidence for H1 is summarized
# separately via the spec-curve Granger plot and the per-cell hit rate; it
# does not feed the H1 meta-regression. This is a deliberate scope choice.
h1_coef <- specs_for_coef |>
  dplyr::filter(level == "system",
                model %in% setdiff(H1_MODELS, "VAR_Granger"))
h1_pval <- pval_long_meta |>
  dplyr::filter(level == "system",
                model %in% setdiff(H1_MODELS, "VAR_Granger"))
h1_meta <- run_h_meta("H1: war reshapes the system network",
                       h1_coef, h1_pval, H_DIR, "H1_system")

# H2 ------------------------------------------------------------
h2_coef <- specs_for_coef |>
  dplyr::filter(level == "state",
                war_variable %in% OWN_MID_IVS,
                model %in% H_STATE_MODELS)
h2_pval <- pval_long_meta |>
  dplyr::filter(level == "state",
                war_variable %in% OWN_MID_IVS,
                model %in% H_STATE_MODELS)
h2_meta <- run_h_meta("H2: war reshapes a state's own connections",
                       h2_coef, h2_pval, H_DIR, "H2_state_own")

# H3 ------------------------------------------------------------
h3_coef <- specs_for_coef |>
  dplyr::filter(level == "state",
                war_variable %in% PWF_IVS,
                model %in% H_STATE_MODELS)
h3_pval <- pval_long_meta |>
  dplyr::filter(level == "state",
                war_variable %in% PWF_IVS,
                model %in% H_STATE_MODELS)
h3_meta <- run_h_meta("H3: war reshapes partners' networks (1-hop spillover)",
                       h3_coef, h3_pval, H_DIR, "H3_state_partner")

# Quick console summary: how much do new indicators move the conditional
# mean of -log10(p) relative to the reference cell?
cat("\n--- Conditional effect of new perturbation indicators on -log10(p) ---\n")
for (term in c("dv_loggedTRUE", "model_no_feTRUE", "model_year_polyTRUE",
                "state_dv_canonicalTRUE",
                "se_methodHC3", "se_methodNW(5)")) {
  if (term %in% names(coef(m_meta_pval))) {
    cf <- coef(m_meta_pval)[term]
    se <- sqrt(diag(vcov(m_meta_pval)))[term]
    cat(sprintf("  %-30s %+.4f  (SE=%.4f, p=%.4f)\n",
                term, unname(cf), unname(se),
                2 * (1 - pnorm(abs(cf / se)))))
  }
}

# Meta-analysis coefficient extraction
extract_meta_coefs <- function(model, model_name) {
  if (is.null(model)) return(tibble())
  sm <- summary(model)$coefficients
  tibble(
    term = rownames(sm), estimate = sm[, "Estimate"],
    se = sm[, "Std. Error"], pval = sm[, ncol(sm)],
    model = model_name
  ) |>
    mutate(
      ci_lower = estimate - 1.96 * se,
      ci_upper = estimate + 1.96 * se,
      significant = pval < 0.05,
      term_clean = case_when(
        term == "has_offensiveTRUE" ~ "Incl. Offensive Alliances",
        term == "has_offensiveFALSE" ~ "No Offensive Alliances",
        term == "has_defensiveTRUE" ~ "Incl. Defensive Alliances",
        term == "has_defensiveFALSE" ~ "No Defensive Alliances",
        term == "has_igoTRUE" ~ "Incl. IGO Co-membership",
        term == "has_igoFALSE" ~ "No IGO Co-membership",
        term == "has_tradeTRUE" ~ "Incl. Bilateral Trade",
        term == "has_tradeFALSE" ~ "No Bilateral Trade",
        term == "has_dcaTRUE" ~ "Incl. Defense Coop. (DCA)",
        term == "has_dcaFALSE" ~ "No Defense Coop. (DCA)",
        term == "is_onsetTRUE" ~ "Onset (vs Count)",
        term == "is_onsetFALSE" ~ "Count (vs Onset)",
        term == "weightingcontig_euclidean" ~ "Contiguity + Euclidean",
        term == "weightingcontig_manhattan" ~ "Contiguity + Manhattan",
        term == "weightingdist_euclidean" ~ "Distance + Euclidean",
        term == "weightingdist_manhattan" ~ "Distance + Manhattan",
        term == "weightingunweighted" ~ "Unweighted",
        term == "transformraw" ~ "Transform: Raw",
        term == "transformlog1p" ~ "Transform: Log(1+x)",
        term == "transformzscore" ~ "Transform: Z-Score",
        term == "levelstate" ~ "Level: State",
        term == "levelsystem" ~ "Level: System",
        term == "modelOLS_lagged" ~ "OLS (Lagged)",
        term == "modelFE_full_ccode_year" ~ "FE (State + Year)",
        term == "modelFE_full_year" ~ "FE (Year Only)",
        term == "modelFE_simple_ccode_year" ~ "FE Simple (State + Year)",
        term == "modelFE_simple_year" ~ "FE Simple (Year Only)",
        term == "war_variablemid_count" ~ "MID Count",
        term == "war_variablelog_mid_count" ~ "Log(1+MID Count)",
        term == "war_variablefatal_mid_count" ~ "Fatal MID Count",
        term == "war_variablelog_fatal_mid_count" ~ "Log(1+Fatal MID Count)",
        term == "war_variablemid_onset_count" ~ "MID Onset Count",
        term == "war_variablelog_mid_onset_count" ~ "Log(1+MID Onset Count)",
        term == "war_variablefatal_mid_onset_count" ~ "Fatal MID Onset Count",
        term == "war_variablelog_fatal_mid_onset_count" ~ "Log(1+Fatal MID Onset Count)",
        term == "dvkl_divergence" ~ "KL Divergence",
        term == "dvjs_divergence" ~ "JS Divergence",
        term == "dvedge_turnover_pct" ~ "Edge Turnover (%)",
        term == "dvdegree_kl" ~ "Degree KL",
        term == "dvdegree_wasserstein" ~ "Degree Wasserstein",
        term == "dvchurn_alt" ~ "Churn (State)",
        # Newly added perturbation indicators
        term == "dv_loggedTRUE"           ~ "DV: log1p-transformed",
        term == "dv_loggedFALSE"          ~ "DV: raw scale",
        term == "model_no_feTRUE"         ~ "Estimator: No FE",
        term == "model_no_feFALSE"        ~ "Estimator: includes FE",
        term == "model_year_polyTRUE"     ~ "Estimator: incl. poly(year, 2)",
        term == "model_year_polyFALSE"    ~ "Estimator: no year polynomial",
        term == "state_dv_canonicalTRUE"  ~ "DV: canonical churn (vs. per-combo)",
        term == "state_dv_canonicalFALSE" ~ "DV: per-combo churn",
        # SE method dummies (only present in p-value meta-regressions)
        term == "se_methodHC1"            ~ "SE method: HC1",
        term == "se_methodHC3"            ~ "SE method: HC3",
        term == "se_methodNW(5)"          ~ "SE method: NW(5)",
        TRUE ~ term |> str_remove("^(model|war_variable|dv)") |>
          str_replace_all("_", " ") |> str_to_title()
      )
    )
}

meta_coefs_all <- bind_rows(
  extract_meta_coefs(m_meta_coef, "DV: Coefficient"),
  extract_meta_coefs(m_meta_pval, "DV: -log10(p)"),
  extract_meta_coefs(m_meta_sig, "DV: Significant (logit)")
)

# =========================================================================== #
# 9. Plots (PNG only)
# =========================================================================== #

cat("\n=== Generating plots ===\n")

# ── 9a. Meta-analysis coefficient plot ──
p_meta <- meta_coefs_all |>
  ggplot(aes(x = estimate, y = reorder(term_clean, estimate), color = significant)) +
  geom_point(size = 2.8) +
  geom_errorbarh(aes(xmin = ci_lower, xmax = ci_upper), height = 0.25) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  facet_wrap(~ model, scales = "free_x", ncol = 3,
             labeller = label_wrap_gen(width = 22)) +
  scale_color_manual(values = c("FALSE" = "gray50", "TRUE" = "black")) +
  scale_y_discrete(labels = function(x) stringr::str_wrap(x, width = 20)) +
  labs(x = "Coefficient", y = NULL,
       title = "Meta-analysis",
       color = "p < 0.05") +
  theme_minimal() +
  theme(panel.grid.minor = element_blank(), legend.position = "bottom",
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12, lineheight = 0.9),
        axis.title.x = element_text(size = 14),
        strip.text = element_text(size = 11, face = "bold"),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12),
        panel.spacing = unit(1.4, "lines"),
        panel.spacing.x = unit(1.8, "lines"),
        plot.margin = margin(t = 10, r = 12, b = 14, l = 10))

ggsave(file.path(FIG_DIR, "fig_meta_coefplot.png"), p_meta, width = 15, height = 14, dpi = 300)

# ── 9b. Tile plot: layers × weighting/transform, faceted by IV ──
cat("  Tile plot...\n")

tile_data <- specs_df |>
  filter(model != "VAR_Granger", model != "VAR_Granger_reverse") |>
  mutate(
    scheme = paste0(
      recode(weighting,
        unweighted = "Unwt.", contig_manhattan = "Contig.\nManh.",
        contig_euclidean = "Contig.\nEucl.", dist_manhattan = "Dist.\nManh.",
        dist_euclidean = "Dist.\nEucl."),
      " / ",
      recode(transform, raw = "Raw", log1p = "Log", zscore = "Z")
    ),
    layers_short = layers |>
      str_replace_all("offensive", "Off") |>
      str_replace_all("defensive", "Def") |>
      str_replace_all("trade", "Trd") |>
      str_replace_all("igo", "IGO") |>
      str_replace_all("dca", "DCA")
  ) |>
  group_by(war_variable, layers_short, scheme) |>
  summarise(
    # One-sided positive at ALPHA_ONESIDE; matches s-curve / forest / console.
    sig_pos_share = mean(is_sig_pos(coef, pval), na.rm = TRUE),
    n_specs = n(), .groups = "drop"
  )

iv_labels_tile <- c(
  mid_count = "MID Count",
  log_mid_count = "Log(1+MID Count)",
  fatal_mid_count = "Fatal MID Count",
  log_fatal_mid_count = "Log(1+Fatal MID Count)",
  mid_onset_count = "MID Onset Count",
  log_mid_onset_count = "Log(1+MID Onset Count)",
  fatal_mid_onset_count = "Fatal MID Onset Count",
  log_fatal_mid_onset_count = "Log(1+Fatal MID Onset)"
)

tile_data <- tile_data |>
  mutate(
    iv_label = recode(war_variable, !!!iv_labels_tile)
  )

# Split into 4 separate plots, each with count + log(1+count) pair
tile_groups <- list(
  list(vars = c("mid_count", "log_mid_count"),
       title = "System: MID count", file = "fig_tile_mid_count.png"),
  list(vars = c("fatal_mid_count", "log_fatal_mid_count"),
       title = "System: fatal MID", file = "fig_tile_fatal_mid_count.png"),
  list(vars = c("mid_onset_count", "log_mid_onset_count"),
       title = "State: MID onset", file = "fig_tile_mid_onset_count.png"),
  list(vars = c("fatal_mid_onset_count", "log_fatal_mid_onset_count"),
       title = "State: fatal MID onset", file = "fig_tile_fatal_mid_onset_count.png")
)

tile_theme <- theme_minimal(base_size = 10) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 10),
        axis.text.y = element_text(size = 10),
        strip.text = element_text(size = 9, face = "bold"),
        plot.title = element_text(size = 15, face = "bold"),
        legend.position = "bottom",
        legend.key.height = unit(0.6, "cm"),
        legend.key.width = unit(1.5, "cm"),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 11, face = "bold"),
        panel.spacing = unit(1.2, "lines"),
        plot.margin = margin(t = 10, r = 12, b = 14, l = 10))

for (tg in tile_groups) {
  td <- tile_data |> filter(war_variable %in% tg$vars) |>
    mutate(iv_label = factor(iv_label, levels = iv_labels_tile[tg$vars]))
  if (nrow(td) == 0) next

  p_tile <- ggplot(td, aes(x = scheme, y = layers_short, fill = sig_pos_share)) +
    geom_tile(color = "white", linewidth = 0.3) +
    scale_fill_distiller(palette = "YlOrRd", direction = 1,
                         limits = c(0, 1), name = "Share Sig+") +
    facet_wrap(~ iv_label, ncol = 2,
               labeller = label_wrap_gen(width = 22)) +
    labs(x = NULL, y = NULL, title = tg$title) +
    tile_theme

  ggsave(file.path(FIG_DIR, tg$file), p_tile,
         width = 16, height = 10, dpi = 200)
}

# ── 9c. Specification curve (Simonsohn-style) ──
# Plot 1: Side-by-side S-curves (system KL + state churn)
# Plot 2: Side-by-side specification indicator panels
cat("  Specification curve...\n")

curve_panels <- list()
indicator_panels <- list()

wv_readable <- c(
  mid_count = "MID\nCount", log_mid_count = "Log(1+MID\nCount)",
  fatal_mid_count = "Fatal\nMID Count", log_fatal_mid_count = "Log(1+Fatal\nMID Count)",
  mid_onset_count = "MID Onset\nCount", log_mid_onset_count = "Log(1+MID\nOnset Count)",
  fatal_mid_onset_count = "Fatal MID\nOnset Count",
  log_fatal_mid_onset_count = "Log(1+Fatal\nMID Onset)"
)

weighting_readable <- c(
  unweighted = "Unwt.",
  contig_manhattan = "Contig.\nManh.",
  contig_euclidean = "Contig.\nEucl.",
  dist_manhattan = "Dist.\nManh.",
  dist_euclidean = "Dist.\nEucl."
)

model_readable <- c(
  OLS_lagged = "OLS\n(Lagged)",
  FE_full_ccode_year = "FE Full\n(St+Yr)",
  FE_full_year = "FE Full\n(Yr)",
  FE_simple_ccode_year = "FE Simp.\n(St+Yr)",
  FE_simple_year = "FE Simp.\n(Yr)",
  FE_full_year_only = "FE\n(Year FE)",
  FE_full_ccode_only = "FE\n(State FE)",
  baseline_noFE = "Baseline",
  baseline_pwf  = "PWF Baseline",
  full_noFE = "OLS Full\n(No FE)",
  full_noFE_year_poly = "OLS Full\n(No FE +\npoly(year,2))",
  state_period_fe = "State + Period FE",
  # Main-text Table 1 (system) ladder
  sys_T1_iv_only       = "Sys T1: IV +\npoly(yr,2)",
  sys_T1_full_ctrl     = "Sys T1: + ctrls",
  sys_T1_ar1           = "Sys T1: + AR(1)",
  sys_T1_order_period  = "Sys T1: + order\nperiod FE"
)

transform_readable <- c(raw = "Raw", log1p = "Log(1+x)", zscore = "Z-Score")

# Abbreviate layer combo labels for indicator plot
layers_readable <- function(x) {
  x |>
    str_replace_all("offensive", "Off") |>
    str_replace_all("defensive", "Def") |>
    str_replace_all("trade", "Trd") |>
    str_replace_all("igo", "IGO") |>
    str_replace_all("dca", "DCA")
}

for (focus in list(
  list(lev = "system", dv_name = "kl_divergence", label = "System: KL Divergence"),
  list(lev = "state", dv_name = "churn_alt", label = "State: Churn")
)) {
  curve_data <- specs_df |>
    filter(level == focus$lev, dv == focus$dv_name,
           !model %in% c("VAR_Granger", "VAR_Granger_reverse")) |>
    arrange(coef) |>
    mutate(
      spec_rank = row_number(),
      # One-sided test: significant if coef > 0 and p/2 < 0.10
      sig = case_when(
        coef > 0 & pval / 2 < 0.10 ~ "sig_pos",
        TRUE ~ "other"
      ),
      # Indicator binaries for each p-value flavor (one-sided, p/2 < 0.10)
      sig_iid     = as.integer(coef > 0 & pval         / 2 < 0.10),
      sig_hc3     = as.integer(coef > 0 & hc3_pval     / 2 < 0.10),
      sig_nw      = as.integer(coef > 0 & nw_pval      / 2 < 0.10),
      sig_iid_fdr = as.integer(coef > 0 & pval_fdr_iid / 2 < 0.10),
      sig_hc3_fdr = as.integer(coef > 0 & pval_fdr_hc3 / 2 < 0.10),
      sig_nw_fdr  = as.integer(coef > 0 & pval_fdr_nw  / 2 < 0.10)
    )

  if (nrow(curve_data) == 0) next

  pct_sig_pos <- round(100 * mean(curve_data$sig == "sig_pos"), 1)

  # Console summary of % significant under each p-value flavor (system OLS
  # specs only; state and VAR rows have NA for hc3/nw and contribute via
  # na.rm = TRUE).
  cat(sprintf("\n  [%s / %s] %% positive + significant (one-sided, p/2 < 0.10):\n",
              focus$lev, focus$dv_name))
  cat(sprintf("    iid       %5.1f%%   iid (FDR)  %5.1f%%\n",
              100 * mean(curve_data$sig_iid,     na.rm = TRUE),
              100 * mean(curve_data$sig_iid_fdr, na.rm = TRUE)))
  cat(sprintf("    HC3       %5.1f%%   HC3 (FDR)  %5.1f%%\n",
              100 * mean(curve_data$sig_hc3,     na.rm = TRUE),
              100 * mean(curve_data$sig_hc3_fdr, na.rm = TRUE)))
  cat(sprintf("    NW(5)     %5.1f%%   NW (FDR)   %5.1f%%\n",
              100 * mean(curve_data$sig_nw,      na.rm = TRUE),
              100 * mean(curve_data$sig_nw_fdr,  na.rm = TRUE)))

  # S-curve panel
  curve_panels[[focus$lev]] <- ggplot(curve_data, aes(x = spec_rank, y = coef, color = sig)) +
    geom_point(size = 0.6, alpha = 0.7) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
    scale_color_manual(values = c(sig_pos = "forestgreen", other = "gray60"),
                       labels = c(sig_pos = "Sig. Positive (one-sided)",
                                  other = "Not Significant / Negative"),
                       name = NULL) +
    labs(y = "Coefficient", x = "Specification (sorted by coefficient)",
         title = focus$label) +
    theme_minimal(base_size = 13) +
    theme(axis.text.x = element_blank(),
          axis.text.y = element_text(size = 11),
          axis.title.y = element_text(size = 12),
          axis.title.x = element_text(size = 11),
          plot.title = element_text(size = 14, face = "bold"),
          legend.text = element_text(size = 10),
          legend.position = "top",
          panel.grid.minor = element_blank())

  # Indicator panel
  choice_labels <- c(
    war_variable = "Conflict\nMeasure",
    weighting    = "Edge\nWeighting",
    transform    = "Transform",
    model        = "Estimator",
    exposure     = "Partner War\nExposure",
    layers       = "Layer\nCombination",
    n_layers     = "Number of\nLayers"
  )

  exposure_readable <- c(
    unweighted = "Unweighted",
    weighted   = "Tie-Weighted",
    none       = "None (Bivariate)"
  )

  indicator_data <- curve_data |>
    mutate(
      war_variable = recode(war_variable, !!!wv_readable),
      weighting = recode(weighting, !!!weighting_readable),
      model = recode(model, !!!model_readable),
      transform = recode(transform, !!!transform_readable),
      exposure = recode(replace_na(exposure, "none"), !!!exposure_readable),
      layers = layers_readable(layers)
    ) |>
    select(spec_rank, war_variable, weighting, transform, model, exposure, layers) |>
    mutate(n_layers = as.character(str_count(layers, "\\+") + 1)) |>
    pivot_longer(-spec_rank, names_to = "choice", values_to = "value") |>
    mutate(
      choice = recode(choice, !!!choice_labels),
      value = ifelse(
        str_count(value, "\\+") >= 2,
        {
          parts <- str_split(value, "\\+")
          sapply(parts, function(p) {
            mid <- ceiling(length(p) / 2)
            paste0(paste(p[1:mid], collapse = "+"), "\n", paste(p[(mid+1):length(p)], collapse = "+"))
          })
        },
        value
      )
    )

  indicator_panels[[focus$lev]] <- ggplot(indicator_data,
                                          aes(x = spec_rank, y = choice, fill = value)) +
    geom_tile(linewidth = 0) +
    scale_fill_viridis_d(option = "turbo", name = "Specification\nChoice") +
    labs(x = "Specification (sorted by coefficient)", y = NULL,
         title = focus$label) +
    theme_minimal(base_size = 13) +
    theme(axis.text.x = element_blank(),
          axis.text.y = element_text(size = 11),
          axis.title.x = element_text(size = 11),
          plot.title = element_text(size = 14, face = "bold"),
          legend.position = "bottom",
          legend.title = element_text(size = 11, face = "bold"),
          legend.text = element_text(size = 8),
          legend.box = "horizontal",
          panel.grid = element_blank()) +
    guides(fill = guide_legend(nrow = 8, byrow = TRUE))
}

# Plot 1: S-curves side by side
if (length(curve_panels) == 2) {
  p_scurves <- curve_panels[["system"]] + curve_panels[["state"]] +
    plot_layout(ncol = 2) +
    plot_annotation(
      title = "Specification curves",
      theme = theme(plot.title = element_text(face = "bold", size = 16))
    )
  ggsave(file.path(FIG_DIR, "fig_spec_curves_combined.png"), p_scurves,
         width = 16, height = 7, dpi = 300)
}

# Plot 2: Indicator panels side by side
if (length(indicator_panels) == 2) {
  p_indicators <- indicator_panels[["system"]] + indicator_panels[["state"]] +
    plot_layout(ncol = 2, guides = "collect") &
    theme(legend.position = "bottom") &
    guides(fill = guide_legend(nrow = 8, byrow = TRUE))
  ggsave(file.path(FIG_DIR, "fig_spec_indicators_combined.png"), p_indicators,
         width = 18, height = 10, dpi = 300)
}

# ── 9c'. P-value distribution plots: ECDF + histogram across SE flavors ──
# Pivot the per-spec p-values to long form so a single ggplot can compare
# the full distribution under iid / HC3 / NW (with and without BH-FDR
# adjustment). State-level TWFE specs and VAR/Granger specs only carry the
# iid p-value; HC3 / NW columns are NA for those rows and are dropped here.
cat("  P-value distribution plots...\n")

pval_long <- specs_df |>
  filter(!model %in% c("VAR_Granger", "VAR_Granger_reverse")) |>
  select(level, dv, model, coef,
         iid          = pval,
         any_of(c(HC1 = "hc1_pval")),
         HC3          = hc3_pval,
         `NW(5)`      = nw_pval,
         `iid (FDR)`  = pval_fdr_iid,
         any_of(c(`HC1 (FDR)` = "pval_fdr_hc1")),
         `HC3 (FDR)`  = pval_fdr_hc3,
         `NW (FDR)`   = pval_fdr_nw) |>
  pivot_longer(cols = any_of(c("iid", "HC1", "HC3", "NW(5)",
                                 "iid (FDR)", "HC1 (FDR)",
                                 "HC3 (FDR)", "NW (FDR)")),
                names_to  = "se_method",
                values_to = "p_value") |>
  filter(!is.na(p_value)) |>
  mutate(
    se_method = factor(se_method,
      levels = c("iid", "HC1", "HC3", "NW(5)",
                  "iid (FDR)", "HC1 (FDR)", "HC3 (FDR)", "NW (FDR)")),
    adj = ifelse(grepl("FDR", se_method), "FDR-adjusted (BH)", "Unadjusted"),
    adj = factor(adj, levels = c("Unadjusted", "FDR-adjusted (BH)")),
    se_base = factor(gsub(" \\(FDR\\)", "", as.character(se_method)),
                       levels = c("iid", "HC1", "HC3", "NW(5)"))
  )

# Restrict to the main four DVs people will look at; the full multiverse
# has 14 DVs and the panel grid would be unreadable.
# Log DVs only (non-log KL / churn dropped from the box / p-value plots).
focus_dvs <- c("log_kl_divergence", "log_churn")
pval_focus <- pval_long |>
  filter(dv %in% focus_dvs) |>
  mutate(dv = factor(dv, levels = focus_dvs,
                       labels = c("log(1+KL)", "log(1+Churn) (canon.)")))

# (i) ECDF — cleanest for comparing distribution shape across SE flavors.
p_pval_ecdf <- ggplot(pval_focus,
    aes(x = p_value, color = se_base, linetype = adj)) +
  stat_ecdf(geom = "step", linewidth = 0.8, alpha = 0.85) +
  geom_vline(xintercept = 0.05, linetype = "dashed",
             color = "gray40", linewidth = 0.4) +
  geom_vline(xintercept = 0.10, linetype = "dotted",
             color = "gray60", linewidth = 0.4) +
  geom_abline(slope = 1, intercept = 0,
              linetype = "dotted", color = "gray70", linewidth = 0.3) +
  facet_wrap(~ dv, scales = "fixed", ncol = 2,
             labeller = label_wrap_gen(width = 22)) +
  scale_color_manual(values = c("iid"   = "gray40",
                                  "HC1"   = "#5A8F3D",
                                  "HC3"   = "#2C6CB0",
                                  "NW(5)" = "#B03366"),
                     name = "SE method") +
  scale_linetype_manual(values = c("Unadjusted" = "solid",
                                     "FDR-adjusted (BH)" = "dashed"),
                        name = NULL) +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.25, 0.5, 0.75, 1)) +
  labs(
    title    = "P-value ECDF",
    x        = "p-value",
    y        = "Cumulative share of specifications"
  ) +
  theme_minimal(base_size = 12) +
  theme(strip.text = element_text(face = "bold", size = 9),
        legend.position = "bottom",
        legend.box      = "horizontal",
        panel.spacing = unit(1.2, "lines"),
        plot.margin = margin(t = 10, r = 12, b = 14, l = 10))

ggsave(file.path(FIG_DIR, "fig_pvalue_ecdf.png"), p_pval_ecdf,
        width = 11, height = 8, dpi = 300, bg = "white")

# (ii) Boxplot — p-value distribution by SE method, faceted by DV.
# x-axis is horizontal text (no angled labels per Jared's preference).
p_pval_box <- ggplot(pval_focus,
    aes(x = se_base, y = p_value, fill = adj)) +
  geom_hline(yintercept = 0.05, linetype = "dashed",
             color = "gray30", linewidth = 0.4) +
  geom_hline(yintercept = 0.10, linetype = "dotted",
             color = "gray60", linewidth = 0.4) +
  geom_boxplot(outlier.size = 0.4, outlier.alpha = 0.4,
               linewidth = 0.4, position = position_dodge(width = 0.75)) +
  facet_wrap(~ dv, scales = "fixed", ncol = 2,
             labeller = label_wrap_gen(width = 22)) +
  scale_fill_manual(values = c("Unadjusted" = "gray85",
                                  "FDR-adjusted (BH)" = "#5A8F3D"),
                     name = NULL) +
  scale_y_continuous(breaks = c(0, 0.05, 0.1, 0.25, 0.5, 0.75, 1)) +
  labs(
    title    = "P-value distribution",
    x        = "Standard error method",
    y        = "p-value"
  ) +
  theme_minimal(base_size = 12) +
  theme(strip.text      = element_text(face = "bold", size = 9),
        legend.position = "bottom",
        axis.text       = element_text(size = 11),     # no angle
        axis.title      = element_text(size = 12),
        panel.spacing = unit(1.2, "lines"),
        plot.margin = margin(t = 10, r = 12, b = 14, l = 10))

ggsave(file.path(FIG_DIR, "fig_pvalue_boxplot.png"), p_pval_box,
        width = 11, height = 8, dpi = 300, bg = "white")

# (iii) Coefficient distribution — one box per DV, dropping VAR/Granger specs,
# paired with the share of specifications that are positive AND significant
# at the one-sided 0.05 level (using is_sig_pos defined above).
coef_focus <- specs_df |>
  filter(!model %in% c("VAR_Granger", "VAR_Granger_reverse"),
         dv %in% focus_dvs) |>
  mutate(dv = factor(dv, levels = focus_dvs,
                       labels = c("log(1+KL)", "log(1+Churn) (canon.)")))

# (iii-a) Coefficient boxplot — faceted with free y so the two DVs (very
# different coefficient scales) each render cleanly.
p_coef_box <- ggplot(coef_focus, aes(x = dv, y = coef)) +
  geom_hline(yintercept = 0, linetype = "dashed",
             color = "gray40", linewidth = 0.4) +
  geom_boxplot(fill = "#2C6CB0", alpha = 0.6,
               outlier.size = 0.4, outlier.alpha = 0.4,
               linewidth = 0.4, width = 0.55) +
  facet_wrap(~ dv, scales = "free", ncol = 2) +
  labs(
    title    = "(a) Coefficient distribution",
    x        = NULL,
    y        = "Estimated coefficient on conflict measure"
  ) +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y  = element_text(size = 11),
        axis.title   = element_text(size = 12),
        strip.text   = element_text(face = "bold", size = 10),
        plot.title   = element_text(face = "bold"),
        plot.margin  = margin(t = 10, r = 12, b = 14, l = 10))

# (iii-b) Share positive & significant (one-sided p < 0.05): blue = pos+sig,
# red = null. Uses is_sig_pos(coef, pval) = coef > 0 & pval/2 < ALPHA_ONESIDE.
coef_sig_df <- coef_focus |>
  filter(!is.na(coef), !is.na(pval)) |>
  mutate(sig_cat = factor(
    if_else(is_sig_pos(coef, pval),
            "Positive & significant (one-sided p<0.05)", "Null"),
    levels = c("Positive & significant (one-sided p<0.05)", "Null"))) |>
  count(dv, sig_cat, .drop = FALSE) |>
  group_by(dv) |>
  mutate(share = n / sum(n)) |>
  ungroup()

p_coef_bar <- ggplot(coef_sig_df, aes(x = dv, y = share, fill = sig_cat)) +
  geom_col(width = 0.55) +
  geom_text(aes(label = ifelse(
              sig_cat == "Positive & significant (one-sided p<0.05)" & share > 0.06,
              scales::percent(share, accuracy = 1), "")),
            position = position_stack(vjust = 0.5),
            color = "white", size = 3, fontface = "bold") +
  scale_fill_manual(values = c(
      "Positive & significant (one-sided p<0.05)" = "#2C6CB0",
      "Null"                                      = "#B03366"),
      name = NULL) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                     expand = expansion(mult = c(0, 0.02))) +
  labs(
    title    = "(b) % positive & significant (one-sided p < 0.05)",
    x        = "Dependent variable",
    y        = "Share of specifications"
  ) +
  theme_minimal(base_size = 12) +
  theme(axis.text       = element_text(size = 11),     # no angle
        axis.title      = element_text(size = 12),
        plot.title      = element_text(face = "bold"),
        legend.position = "bottom",
        plot.margin = margin(t = 10, r = 12, b = 14, l = 10))

# Combined: coefficient boxplot | %-significant bar (single row).
p_coef_combined <- p_coef_box + p_coef_bar + patchwork::plot_layout(widths = c(1, 1))
ggsave(file.path(FIG_DIR, "fig_coef_boxplot.png"), p_coef_combined,
        width = 13, height = 6, dpi = 300, bg = "white")

# (iv) Combined panel: coefficient boxplot | p-value boxplot | s-curve
# Built with patchwork. S-curves require pre-built panels (curve_panels list
# from above), so we stack the system + state s-curves vertically and pair
# them with the boxplots on the left.
if (length(curve_panels) >= 1) {
  scurve_panel <- if (length(curve_panels) == 2) {
    curve_panels[["system"]] / curve_panels[["state"]] +
      plot_layout(heights = c(1, 1))
  } else {
    curve_panels[[1]]
  }

  p_combined <- (p_coef_box / p_pval_box) | scurve_panel +
    plot_layout(widths = c(1, 1)) +
    plot_annotation(
      title = "Multiverse summary",
      theme = theme(plot.title = element_text(face = "bold", size = 15))
    )
  ggsave(file.path(FIG_DIR, "fig_multiverse_combined.png"), p_combined,
          width = 18, height = 10, dpi = 300, bg = "white")
}

# Save the long-form pvalue tibble for downstream / reproducible plotting.
write_csv(pval_long,
            file.path(OUTPUT_DIR, "robustness_pvalue_long.csv"))
cat("  -> fig_pvalue_ecdf.png\n")
cat("  -> fig_pvalue_boxplot.png\n")
cat("  -> fig_coef_boxplot.png\n")
cat("  -> fig_multiverse_combined.png\n")
cat("  -> robustness_pvalue_long.csv\n")

# ── 9d. Forest plot (2 columns: System | State; 2 rows: MID Count | Fatal MID Count) ──
cat("  Forest plot...\n")

forest_all <- specs_df |>
  filter(model %in% c("OLS_lagged", "FE_full_ccode_year"),
         weighting == "unweighted", transform == "raw",
         layers == "defensive+igo") |>
  mutate(
    ci_lo = coef - 1.96 * se,
    ci_hi = coef + 1.96 * se,
    dv_label = case_when(
      dv == "kl_divergence"          ~ "KL Divergence",
      dv == "js_divergence"          ~ "JS Divergence",
      dv == "edge_turnover_pct"      ~ "Edge Turnover %",
      dv == "degree_kl"              ~ "Degree KL",
      dv == "degree_wasserstein"     ~ "Degree Wasserstein",
      dv == "log_kl_divergence"      ~ "log(1+KL)",
      dv == "log_js_divergence"      ~ "log(1+JS)",
      dv == "log_edge_turnover_pct"  ~ "log(1+Edge Turn.)",
      dv == "log_degree_kl"          ~ "log(1+Deg. KL)",
      dv == "log_degree_wasserstein" ~ "log(1+Deg. Wass.)",
      dv == "churn"                  ~ "Churn",
      dv == "log_churn"              ~ "log(1+Churn)",
      dv == "churn_alt"              ~ "Churn (alt)",
      dv == "log_churn_alt"          ~ "log(1+Churn alt)",
      TRUE ~ dv
    ),
    row_label = paste0(dv_label, "\n(",
                       ifelse(model == "OLS_lagged", "OLS", "FE"), ")"),
    level = str_to_title(level),
    # One-sided positive at ALPHA_ONESIDE — same convention as the s-curve,
    # tile, and console summary.
    sig = is_sig_pos(coef, pval)
  )

forest_theme <- theme_minimal(base_size = 13) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank(),
        strip.text = element_text(face = "bold", size = 12),
        plot.title = element_text(size = 13, face = "bold"),
        axis.text.y = element_text(size = 15),
        axis.text.x = element_text(size = 11),
        axis.title.x = element_text(size = 13))

make_forest_panel <- function(data, wv, title_text) {
  d <- data |> filter(war_variable == wv)
  if (nrow(d) == 0) return(NULL)
  ggplot(d, aes(x = coef, y = reorder(row_label, coef))) +
    geom_point(aes(color = sig), size = 3) +
    geom_errorbarh(aes(xmin = ci_lo, xmax = ci_hi, color = sig), height = 0.2) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    scale_color_manual(values = c("FALSE" = "gray50", "TRUE" = "forestgreen"),
                       name = sprintf("Sig. positive (one-sided, alpha = %.2f)",
                                       ALPHA_ONESIDE)) +
    labs(x = "Coefficient (95% CI)", y = NULL, title = title_text) +
    forest_theme
}

p_sys_mid <- make_forest_panel(forest_all |> filter(level == "System"),
                                "mid_count", "System: MID Count")
p_state_mid <- make_forest_panel(forest_all |> filter(level == "State"),
                                  "mid_count", "State: MID Count")

p_sys_fatal <- make_forest_panel(forest_all |> filter(level == "System"),
                                  "fatal_mid_count", "System: Fatal MID Count")
p_state_fatal <- make_forest_panel(forest_all |> filter(level == "State"),
                                    "fatal_mid_count", "State: Fatal MID Count")

if (!is.null(p_sys_mid)) {
  p_forest <- (p_sys_mid | (p_state_mid %||% plot_spacer())) /
              (p_sys_fatal | (p_state_fatal %||% plot_spacer())) +
    plot_layout(guides = "collect") +
    plot_annotation(
      title = "Effect sizes across DVs",
      theme = theme(plot.title = element_text(face = "bold", size = 16))
    ) &
    theme(legend.position = "bottom")

  ggsave(file.path(FIG_DIR, "fig_forest_dvs.png"), p_forest,
         width = 14, height = 14, dpi = 300)
}

# ── 9e. FEVD horizon evolution (computed from VAR fits in main spec) ──
cat("  FEVD horizon plot...\n")

# Re-fit the main-spec VAR and compute FEVD at horizons 1:15
fevd_horizon_data <- list()

# Defensive: rebuild analysis_data_main inline if it's not in the
# environment (interactive sourcing of post-processing chunks can skip
# the rebuild block above). Headline spec: unweighted, raw, defensive+IGO.
if (!exists("analysis_data_main") || is.null(analysis_data_main)) {
  cat("  Rebuilding analysis_data_main (headline spec) inline...\n")
  layers_main_   <- c("layer_defensive", "layer_igo")
  manhattan_main_ <- compute_manhattan_vec(dyad_data, layers_main_, "unweighted")
  dv_df_main_    <- compute_system_dvs(dyad_data, manhattan_main_)
  yr_main_       <- get_valid_year_range(layers_main_)
  analysis_data_main <- dv_df_main_ |>
    dplyr::left_join(system_war,   by = "year") |>
    dplyr::left_join(sys_controls, by = "year") |>
    dplyr::filter(year >= yr_main_[1], year <= yr_main_[2]) |>
    dplyr::filter(!is.na(kl_divergence)) |>
    dplyr::mutate(
      log_kl_divergence      = log1p(kl_divergence),
      log_js_divergence      = log1p(js_divergence),
      log_edge_turnover_pct  = log1p(edge_turnover_pct),
      log_degree_kl          = log1p(degree_kl),
      log_degree_wasserstein = log1p(degree_wasserstein)
    )
}

if (!is.null(analysis_data_main)) {
for (wv in system_war_vars) {
  for (dv in granger_dvs) {
    fevd_fit <- tryCatch({
      var_block <- analysis_data_main |>
        select(all_of(c(wv, dv, "n_democracies", "n_autocracies"))) |>
        na.omit()
      if (nrow(var_block) < 30) stop("too few")

      var_endo <- var_block |> select(all_of(c(wv, dv))) |> as.ts()
      var_exog <- var_block |>
        select(n_democracies, n_autocracies) |> mutate_all(log1p) |> as.ts()

      opt_lag <- tryCatch({
        sel <- VARselect(var_endo, lag.max = 8, type = "both", exogen = var_exog)
        max(1L, as.integer(sel$selection["AIC(n)"]))
      }, error = function(e) 2L)

      vf <- VAR(var_endo, p = opt_lag, type = "both", exogen = var_exog)
      fd <- fevd(vf, n.ahead = 15, boot = TRUE, nboot = 100)
      # Extract share of DV variance explained by war_var at each horizon
      fd_dv <- fd[[dv]]
      tibble(
        horizon = 1:15,
        war_share = fd_dv[, wv],
        dv_share = fd_dv[, dv]
      )
    }, error = function(e) NULL)

    if (!is.null(fevd_fit)) {
      fevd_fit$war_variable <- wv
      fevd_fit$dv <- dv
      fevd_horizon_data[[length(fevd_horizon_data) + 1]] <- fevd_fit
    }
  }
}
}  # end if (!is.null(analysis_data_main))

if (length(fevd_horizon_data) > 0) {
  fevd_all <- bind_rows(fevd_horizon_data)

  wv_labels <- c(
    mid_count = "MID Count", log_mid_count = "Log(1+MIDs)",
    fatal_mid_count = "Fatal MIDs", log_fatal_mid_count = "Log(1+Fatal MIDs)"
  )
  dv_labels <- c(
    kl_divergence = "KL Divergence", js_divergence = "JS Divergence",
    degree_kl = "Degree KL", degree_wasserstein = "Degree Wasserstein"
  )

  fevd_all <- fevd_all |>
    mutate(
      wv_label = wv_labels[war_variable],
      dv_label = dv_labels[dv]
    ) |>
    # The VAR sweep includes log-transformed DVs, but this horizon figure
    # shows the four raw network divergence measures. Drop variants without a
    # display label so they don't render as an empty "NA" facet.
    filter(!is.na(dv_label), !is.na(wv_label))

  p_fevd <- ggplot(fevd_all, aes(x = horizon, y = war_share * 100,
                                  color = wv_label, linetype = wv_label)) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 1.8) +
    facet_wrap(~ dv_label, scales = "free_y", ncol = 2,
               labeller = label_wrap_gen(width = 22)) +
    scale_color_brewer(palette = "Dark2") +
    scale_linetype_manual(values = c("solid", "dashed", "dotted", "dotdash")) +
    labs(x = "Forecast Horizon", y = "% of Variance Explained by Conflict",
         title = "FEVD horizon",
         color = "Conflict Measure", linetype = "Conflict Measure") +
    theme_minimal(base_size = 13) +
    theme(legend.position = "bottom", panel.grid.minor = element_blank(),
          strip.text = element_text(face = "bold", size = 9),
          panel.spacing = unit(1.2, "lines"),
          panel.spacing.x = unit(1.5, "lines"),
          plot.margin = margin(t = 10, r = 12, b = 14, l = 10))

  ggsave(file.path(FIG_DIR, "fig_fevd_horizon.png"), p_fevd, width = 12, height = 8, dpi = 300)
}

# ── 9f. Reverse vs Forward Granger comparison ──
cat("  Granger comparison plot...\n")

granger_compare <- specs_df |>
  filter(model %in% c("VAR_Granger", "VAR_Granger_reverse")) |>
  mutate(
    direction = ifelse(model == "VAR_Granger", "War → Network", "Network → War"),
    sig = pval < 0.05
  )

if (nrow(granger_compare) > 0) {
  granger_wv_labels <- c(
    mid_count = "MID\nCount",
    log_mid_count = "Log(1+MID\nCount)",
    fatal_mid_count = "Fatal MID\nCount",
    log_fatal_mid_count = "Log(1+Fatal\nMID Count)",
    mid_onset_count = "MID Onset\nCount",
    log_mid_onset_count = "Log(1+MID\nOnset Count)",
    fatal_mid_onset_count = "Fatal MID\nOnset Count",
    log_fatal_mid_onset_count = "Log(1+Fatal\nMID Onset Count)"
  )
  granger_dv_labels <- c(
    kl_divergence = "KL\nDivergence",
    js_divergence = "JS\nDivergence",
    edge_turnover_pct = "Edge\nTurnover (%)",
    degree_kl = "Degree\nKL",
    degree_wasserstein = "Degree\nWasserstein"
  )

  granger_summary <- granger_compare |>
    group_by(direction, dv, war_variable) |>
    summarise(
      pct_sig = mean(sig) * 100,
      median_pval = median(pval, na.rm = TRUE),
      n = n(), .groups = "drop"
    ) |>
    mutate(
      wv_label = recode(war_variable, !!!granger_wv_labels),
      dv_label = recode(dv, !!!granger_dv_labels),
      wv_label = factor(wv_label, levels = c(
        "MID\nCount", "Log(1+MID\nCount)",
        "Fatal MID\nCount", "Log(1+Fatal\nMID Count)",
        "MID Onset\nCount", "Log(1+MID\nOnset Count)",
        "Fatal MID\nOnset Count", "Log(1+Fatal\nMID Onset Count)"
      ))
    )

  p_granger <- ggplot(granger_summary, aes(x = dv_label, y = pct_sig, fill = direction)) +
    geom_col(position = position_dodge(width = 0.7), width = 0.6) +
    facet_wrap(~ wv_label, ncol = 2,
               labeller = label_wrap_gen(width = 22)) +
    scale_fill_manual(values = c("War \u2192 Network" = "steelblue",
                                 "Network \u2192 War" = "firebrick")) +
    geom_hline(yintercept = 5, linetype = "dashed", color = "gray40") +
    labs(x = NULL, y = "% Significant (p < 0.05)",
         title = "Forward vs. reverse Granger",
         fill = "Direction") +
    theme_minimal(base_size = 13) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 10),
          strip.text = element_text(face = "bold", size = 9),
          legend.position = "top",
          legend.title = element_text(size = 13, face = "bold"),
          legend.text = element_text(size = 12),
          legend.key.size = unit(1, "cm"),
          panel.grid.minor = element_blank(),
          panel.spacing = unit(1.2, "lines"),
          plot.margin = margin(t = 10, r = 12, b = 14, l = 10))

  ggsave(file.path(FIG_DIR, "fig_granger_forward_vs_reverse.png"), p_granger,
         width = 14, height = 8, dpi = 300)
}

cat("  Plots saved.\n")

# =========================================================================== #
# 10. Leave-One-Layer-Out
# =========================================================================== #

cat("\n=== Leave-One-Layer-Out ===\n")

loo_results <- list()
for (drop_layer in all_layer_vars) {
  remaining <- setdiff(all_layer_vars, drop_layer)
  if (length(remaining) == 0) next

  loo_manhattan <- compute_manhattan_vec(dyad_data, remaining, "unweighted")
  loo_dvs <- compute_system_dvs(dyad_data, loo_manhattan)
  loo_analysis <- loo_dvs |>
    inner_join(system_war, by = "year") |>
    inner_join(sys_controls, by = "year") |>
    filter(!is.na(kl_divergence))
  if (nrow(loo_analysis) < 20) next

  lc <- gsub("layer_", "", drop_layer)

  for (wv in system_war_vars) {
    # Granger with exogenous controls
    g <- tryCatch({
      var_block <- loo_analysis |>
        select(all_of(c(wv, "kl_divergence", "n_democracies", "n_autocracies"))) |>
        na.omit()
      if (nrow(var_block) < 30) stop("too few")

      var_endo <- var_block |> select(all_of(c(wv, "kl_divergence"))) |> as.ts()
      var_exog <- var_block |> select(n_democracies, n_autocracies) |>
        mutate_all(log1p) |> as.ts()

      opt_lag <- tryCatch({
        sel <- VARselect(var_endo, lag.max = 8, type = "both", exogen = var_exog)
        max(1L, as.integer(sel$selection["AIC(n)"]))
      }, error = function(e) 2L)

      vf <- VAR(var_endo, p = opt_lag, type = "both", exogen = var_exog)
      gc <- causality(vf, cause = wv, boot.runs = 100)
      list(f = as.numeric(gc$Granger$statistic),
           p = as.numeric(gc$Granger$p.value),
           vf = vf)
    }, error = function(e) NULL)

    if (!is.null(g)) {
      loo_results[[length(loo_results) + 1]] <- tibble(
        dropped_layer = lc, war_var = wv, model = "Granger",
        coef = g$f, pval = g$p
      )
      # Reverse Granger
      g_rev <- tryCatch({
        gc_rev <- causality(g$vf, cause = "kl_divergence", boot.runs = 100)
        list(f = as.numeric(gc_rev$Granger$statistic),
             p = as.numeric(gc_rev$Granger$p.value))
      }, error = function(e) NULL)
      if (!is.null(g_rev)) {
        loo_results[[length(loo_results) + 1]] <- tibble(
          dropped_layer = lc, war_var = wv, model = "Granger_reverse",
          coef = g_rev$f, pval = g_rev$p
        )
      }
    }
  }

  cat(sprintf("  Drop %s: done\n", lc))
}

loo_summary <- bind_rows(loo_results)
cat("\nLOO summary:\n")
print(loo_summary |> arrange(model, war_var))

# =========================================================================== #
# 11. Layer-Specific Analyses
# =========================================================================== #

cat("\n=== Layer-Specific Analyses ===\n")

layer_specific_ts <- readRDS(file.path(OUTPUT_DIR, "layer_specific_measures.rds"))
layer_results <- list()

for (lv in unique(layer_specific_ts$layer)) {
  lc <- gsub("layer_", "", lv)
  ld <- layer_specific_ts |>
    filter(layer == lv, !is.na(kl_divergence)) |>
    inner_join(system_war, by = "year") |>
    inner_join(sys_controls, by = "year")
  if (nrow(ld) < 20) next

  for (wv in system_war_vars) {
    g <- tryCatch({
      var_block <- ld |>
        select(all_of(c(wv, "kl_divergence", "n_democracies", "n_autocracies"))) |>
        na.omit()
      if (nrow(var_block) < 30) stop("too few")

      var_endo <- var_block |> select(all_of(c(wv, "kl_divergence"))) |> as.ts()
      var_exog <- var_block |> select(n_democracies, n_autocracies) |>
        mutate_all(log1p) |> as.ts()

      opt_lag <- tryCatch({
        sel <- VARselect(var_endo, lag.max = 8, type = "both", exogen = var_exog)
        max(1L, as.integer(sel$selection["AIC(n)"]))
      }, error = function(e) 2L)

      vf <- VAR(var_endo, p = opt_lag, type = "both", exogen = var_exog)
      gc <- causality(vf, cause = wv, boot.runs = 100)
      list(f = as.numeric(gc$Granger$statistic),
           p = as.numeric(gc$Granger$p.value),
           vf = vf)
    }, error = function(e) NULL)

    if (!is.null(g)) {
      layer_results[[length(layer_results) + 1]] <- tibble(
        layer = lc, model = "Granger", war_var = wv,
        coef = g$f, pval = g$p, n_obs = nrow(ld)
      )
      # Reverse Granger
      g_rev <- tryCatch({
        gc_rev <- causality(g$vf, cause = "kl_divergence", boot.runs = 100)
        list(f = as.numeric(gc_rev$Granger$statistic),
             p = as.numeric(gc_rev$Granger$p.value))
      }, error = function(e) NULL)
      if (!is.null(g_rev)) {
        layer_results[[length(layer_results) + 1]] <- tibble(
          layer = lc, model = "Granger_reverse", war_var = wv,
          coef = g_rev$f, pval = g_rev$p, n_obs = nrow(ld)
        )
      }
    }
  }
  cat(sprintf("  %s: %d obs\n", lc, nrow(ld)))
}

layer_analysis_summary <- bind_rows(layer_results)

# =========================================================================== #
# 12. Save
# =========================================================================== #

cat("\nSaving...\n")

robustness <- list(
  specifications = specs_df,
  meta_coef = m_meta_coef,
  meta_pval = m_meta_pval,
  meta_sig = m_meta_sig,
  meta_coefs_plot_data = meta_coefs_all,
  leave_one_out = loo_summary,
  layer_specific = layer_analysis_summary,
  fevd_horizon = if (length(fevd_horizon_data) > 0) bind_rows(fevd_horizon_data) else NULL
)

saveRDS(robustness, file.path(OUTPUT_DIR, "robustness_results.rds"))
write_csv(specs_df, file.path(OUTPUT_DIR, "robustness_specifications.csv"))
write_csv(meta_coefs_all, file.path(OUTPUT_DIR, "robustness_meta_coefs.csv"))

# =========================================================================== #
# 13. LaTeX Table Exports
# =========================================================================== #

cat("\n=== LaTeX table exports ===\n")

# ── 13a. Meta-regression table (texreg for lm/glm) ──
texreg(list(m_meta_coef, m_meta_pval, m_meta_sig),
       file = file.path(APPX_DIR, "robustness_meta_regression.tex"),
       custom.model.names = c("Coefficient", "$-\\log_{10}(p)$", "Significant (Logit)"),
       caption = "Meta-Regression: What Drives Specification-Level Results?",
       label = "tab:meta_regression",
       float.pos = "htbp",
       booktabs = TRUE, use.packages = FALSE,
       stars = c(0.001, 0.01, 0.05, 0.1),
       symbol = "$^{\\dagger}$",
       custom.note = "Each observation is one specification across the full multiverse (system + state, raw + logged DV, raw + logged IV, all layer combos and weighting schemes). All models suppress the intercept. $^{***}p<0.001$; $^{**}p<0.01$; $^{*}p<0.05$; $^{\\dagger}p<0.1$.")
cat("  Meta-regression → results/appendix/robustness_meta_regression.tex\n")

# ── 13b. Specification curve summary table ──
#
# Three p-value flavors per spec, all already saved in specs_df by the
# system-level OLS_lagged / LP_h5 blocks above:
#   pval      — iid OLS (unadjusted)
#   hc3_pval  — HC3 White-style heteroskedasticity-robust
#   nw_pval   — Newey-West HAC (lag 5)
#
# State-level specs and VAR/Granger specs don't carry hc3_pval / nw_pval
# (they use different SE machinery: feols TWFE clustering and bootstrap-based
# Granger F-tests respectively). Those rows show NA for the system-only
# columns and the % calculations use na.rm = TRUE so they reflect only the
# system-level subsamples.
#
# Multiple-comparison adjustment: Benjamini-Hochberg FDR applied within
# each (level, dv) cell, separately for each p-value flavor. This block
# was relocated upstream (right after the specs_df enrichment, around
# the dv_logged / model_no_fe mutate) so the FDR columns are available
# to every downstream consumer including the spec-curve plot. Kept here
# as a no-op safety net: if pval_fdr_iid is somehow missing, recompute.
if (!"pval_fdr_iid" %in% names(specs_df)) {
  specs_df <- specs_df |>
    group_by(level, dv) |>
    mutate(
      pval_fdr_iid = if (sum(!is.na(pval))     > 0) p.adjust(pval,     method = "BH") else NA_real_,
      pval_fdr_hc1 = if ("hc1_pval" %in% names(cur_data()) && sum(!is.na(hc1_pval)) > 0)
                       p.adjust(hc1_pval, method = "BH") else NA_real_,
      pval_fdr_hc3 = if ("hc3_pval" %in% names(cur_data()) && sum(!is.na(hc3_pval)) > 0)
                       p.adjust(hc3_pval, method = "BH") else NA_real_,
      pval_fdr_nw  = if ("nw_pval"  %in% names(cur_data()) && sum(!is.na(nw_pval))  > 0)
                       p.adjust(nw_pval,  method = "BH") else NA_real_
    ) |>
    ungroup()
}

spec_summary <- specs_df |>
  filter(!model %in% c("VAR_Granger", "VAR_Granger_reverse")) |>
  group_by(level, dv) |>
  summarise(
    N = n(),
    Mean_Coef = mean(coef, na.rm = TRUE),
    Median_Coef = median(coef, na.rm = TRUE),
    SD_Coef = sd(coef, na.rm = TRUE),
    Pct_Positive = 100 * mean(coef > 0, na.rm = TRUE),
    # One-sided test at p/2 < 0.10 (directional hypothesis), unadjusted
    Pct_Sig_Pos_iid = 100 * mean(coef > 0 & pval     / 2 < 0.10, na.rm = TRUE),
    Pct_Sig_Neg_iid = 100 * mean(coef < 0 & pval     / 2 < 0.10, na.rm = TRUE),
    Pct_Sig_Pos_hc1 = if ("hc1_pval" %in% names(cur_data())) 100 * mean(coef > 0 & hc1_pval / 2 < 0.10, na.rm = TRUE) else NA_real_,
    Pct_Sig_Neg_hc1 = if ("hc1_pval" %in% names(cur_data())) 100 * mean(coef < 0 & hc1_pval / 2 < 0.10, na.rm = TRUE) else NA_real_,
    Pct_Sig_Pos_hc3 = 100 * mean(coef > 0 & hc3_pval / 2 < 0.10, na.rm = TRUE),
    Pct_Sig_Neg_hc3 = 100 * mean(coef < 0 & hc3_pval / 2 < 0.10, na.rm = TRUE),
    Pct_Sig_Pos_nw  = 100 * mean(coef > 0 & nw_pval  / 2 < 0.10, na.rm = TRUE),
    Pct_Sig_Neg_nw  = 100 * mean(coef < 0 & nw_pval  / 2 < 0.10, na.rm = TRUE),
    # FDR-adjusted (BH) within (level, dv)
    Pct_Sig_Pos_iid_fdr = 100 * mean(coef > 0 & pval_fdr_iid / 2 < 0.10, na.rm = TRUE),
    Pct_Sig_Neg_iid_fdr = 100 * mean(coef < 0 & pval_fdr_iid / 2 < 0.10, na.rm = TRUE),
    Pct_Sig_Pos_hc1_fdr = if ("pval_fdr_hc1" %in% names(cur_data())) 100 * mean(coef > 0 & pval_fdr_hc1 / 2 < 0.10, na.rm = TRUE) else NA_real_,
    Pct_Sig_Neg_hc1_fdr = if ("pval_fdr_hc1" %in% names(cur_data())) 100 * mean(coef < 0 & pval_fdr_hc1 / 2 < 0.10, na.rm = TRUE) else NA_real_,
    Pct_Sig_Pos_hc3_fdr = 100 * mean(coef > 0 & pval_fdr_hc3 / 2 < 0.10, na.rm = TRUE),
    Pct_Sig_Neg_hc3_fdr = 100 * mean(coef < 0 & pval_fdr_hc3 / 2 < 0.10, na.rm = TRUE),
    Pct_Sig_Pos_nw_fdr  = 100 * mean(coef > 0 & pval_fdr_nw  / 2 < 0.10, na.rm = TRUE),
    Pct_Sig_Neg_nw_fdr  = 100 * mean(coef < 0 & pval_fdr_nw  / 2 < 0.10, na.rm = TRUE),
    .groups = "drop"
  )
write_csv(spec_summary, file.path(OUTPUT_DIR, "robustness_spec_curve_summary.csv"))

dv_labels <- c(
  # raw system DVs
  kl_divergence = "KL Divergence", js_divergence = "JS Divergence",
  edge_turnover_pct = "Edge Turnover (\\%)", degree_kl = "Degree KL",
  degree_wasserstein = "Degree Wasserstein",
  # log1p system DVs
  log_kl_divergence = "log(1+KL)", log_js_divergence = "log(1+JS)",
  log_edge_turnover_pct = "log(1+Edge Turn.)",
  log_degree_kl = "log(1+Deg. KL)",
  log_degree_wasserstein = "log(1+Deg. Wass.)",
  # state DVs (raw + logged, per-combo + canonical)
  churn_alt     = "Churn (State, per-combo)",
  log_churn_alt = "log(1+Churn) (State, per-combo)",
  churn         = "Churn (State, canonical)",
  log_churn     = "log(1+Churn) (State, canonical)"
)
tex_sc <- c(
  "\\begin{table}[htbp]", "\\centering",
  "\\caption{Specification Curve Summary (One-Sided Tests, $p/2 < 0.10$). \\% Sig.$+$ columns report the share of specifications with positive coefficient and $p/2 < 0.10$ under each SE method (iid, HC3, Newey-West) and after Benjamini-Hochberg FDR adjustment within each (level, DV) cell.}",
  "\\label{tab:spec_curve}",
  "\\footnotesize",
  "\\begin{tabular}{llrrrrrrrrrrr}",
  "\\toprule",
  paste("& & & & & & &",
         "\\multicolumn{3}{c}{\\% Sig.$+$ unadjusted} &",
         "\\multicolumn{3}{c}{\\% Sig.$+$ FDR adj.} \\\\"),
  "\\cmidrule(lr){8-10} \\cmidrule(lr){11-13}",
  paste("Level & DV & N & Mean & Median & SD & \\% Pos. &",
         "iid & HC3 & NW & iid & HC3 & NW \\\\"),
  "\\midrule"
)
for (i in seq_len(nrow(spec_summary))) {
  r <- spec_summary[i, ]
  dv_lab <- ifelse(r$dv %in% names(dv_labels), dv_labels[r$dv], r$dv)
  tex_sc <- c(tex_sc, sprintf(
    "%s & %s & %s & %.4f & %.4f & %.4f & %.1f & %.1f & %.1f & %.1f & %.1f & %.1f & %.1f \\\\",
    str_to_title(r$level), dv_lab, format(r$N, big.mark = ","),
    r$Mean_Coef, r$Median_Coef, r$SD_Coef,
    r$Pct_Positive,
    r$Pct_Sig_Pos_iid,     r$Pct_Sig_Pos_hc3,     r$Pct_Sig_Pos_nw,
    r$Pct_Sig_Pos_iid_fdr, r$Pct_Sig_Pos_hc3_fdr, r$Pct_Sig_Pos_nw_fdr
  ))
}
tex_sc <- c(tex_sc, "\\bottomrule", "\\end{tabular}",
  "\\begin{tablenotes}\\small",
  "\\item Significance assessed using one-sided tests ($p/2 < 0.10$) given the directional hypothesis that conflict increases network reconfiguration. iid: conventional OLS standard errors. HC3: heteroskedasticity-robust (jackknife-style; Long-Ervin 2000 recommend for $N < 250$). NW: Newey-West HAC with bandwidth 5. FDR adj.: Benjamini-Hochberg adjustment applied within each (level, DV) cell. State-level and VAR/Granger specifications use cluster-robust or bootstrap-based inference and contribute only the iid column; HC3 and NW columns reflect system-level OLS / LP specifications only.",
  "\\end{tablenotes}",
  "\\end{table}")
writeLines(tex_sc, file.path(APPX_DIR, "robustness_spec_curve_summary.tex"))
cat("  Spec curve summary → results/appendix/robustness_spec_curve_summary.tex\n")

# ── 13c. Granger forward vs reverse table ──
granger_specs <- specs_df |>
  filter(model %in% c("VAR_Granger", "VAR_Granger_reverse")) |>
  group_by(model, war_variable, dv) |>
  summarise(
    N = n(),
    Mean_Fstat = mean(coef, na.rm = TRUE),
    Mean_Pval  = mean(pval, na.rm = TRUE),
    Pct_Sig    = 100 * mean(pval < 0.05, na.rm = TRUE),
    .groups = "drop"
  )
write_csv(granger_specs, file.path(OUTPUT_DIR, "robustness_granger_comparison.csv"))

tex_gr <- c(
  "\\begin{table}[htbp]", "\\centering",
  "\\caption{Granger Causality: Forward vs.\\ Reverse}", "\\label{tab:granger_compare}",
  "\\small",
  "\\begin{tabular}{llrrrr}", "\\toprule",
  "Direction & Conflict Measure & N & Mean $F$ & Mean $p$ & \\% Sig. \\\\",
  "\\midrule"
)
dir_labels <- c(VAR_Granger = "War $\\rightarrow$ Network",
                VAR_Granger_reverse = "Network $\\rightarrow$ War")
for (i in seq_len(nrow(granger_specs))) {
  r <- granger_specs[i, ]
  tex_gr <- c(tex_gr, sprintf(
    "%s & %s & %d & %.3f & %.4f & %.1f \\\\",
    dir_labels[r$model], r$war_variable, r$N,
    r$Mean_Fstat, r$Mean_Pval, r$Pct_Sig
  ))
}
tex_gr <- c(tex_gr, "\\bottomrule", "\\end{tabular}", "\\end{table}")
writeLines(tex_gr, file.path(APPX_DIR, "robustness_granger_comparison.tex"))
cat("  Granger comparison → results/appendix/robustness_granger_comparison.tex\n")

# ── 13d. FEVD horizon table ──
if (length(fevd_horizon_data) > 0) {
  fevd_h <- bind_rows(fevd_horizon_data) |>
    filter(horizon %in% c(1, 5, 10, 15))
  write_csv(fevd_h, file.path(OUTPUT_DIR, "robustness_fevd_horizon.csv"))

  tex_fh <- c(
    "\\begin{table}[htbp]", "\\centering",
    "\\caption{Forecast Error Variance Decomposition by Horizon}",
    "\\label{tab:fevd_horizon}",
    "\\small",
    "\\begin{tabular}{llrr}", "\\toprule",
    "Conflict Measure & DV & Horizon & \\% Variance Explained \\\\",
    "\\midrule"
  )
  for (i in seq_len(nrow(fevd_h))) {
    r <- fevd_h[i, ]
    tex_fh <- c(tex_fh, sprintf(
      "%s & %s & %d & %.2f \\\\",
      r$war_variable, r$dv, r$horizon, r$war_share * 100
    ))
  }
  tex_fh <- c(tex_fh, "\\bottomrule", "\\end{tabular}", "\\end{table}")
  writeLines(tex_fh, file.path(APPX_DIR, "robustness_fevd_horizon.tex"))
  cat("  FEVD horizon → results/appendix/robustness_fevd_horizon.tex\n")
}

# =========================================================================== #
# 14. System OLS SE-ladder appendix table (pulled from spec curve)
# =========================================================================== #
# The 4-method SE ladder (iid / HC1 / HC3 / NW HAC) is now part of the
# OLS_lagged multiverse cells — every cell carries all four SEs and they
# enter the spec-curve coefficient and -log10(p) distributions through
# the `se_method` factor. This block extracts the canonical headline
# cell of each of the 4 system IVs and writes the appendix HAC table
# (formerly produced by 19_appendix_hac.R) directly from specs_df.

cat("\n========== System OLS SE-ladder appendix table ==========\n")

# Headline cell: log_kl_divergence DV, raw transform, unweighted method,
# layers = the 3 main-spec layers (offensive + defensive + igo).
HEADLINE_LAYERS <- "offensive+defensive+igo"

hac_appendix_rows <- specs_df |>
  filter(level == "system",
         model == "OLS_lagged",
         dv == "log_kl_divergence",
         transform == "raw",
         weighting == "unweighted",
         layers == HEADLINE_LAYERS) |>
  mutate(iv_label = recode(war_variable,
                             "mid_count"            = "MID count (t-1)",
                             "log_mid_count"        = "log(1 + MID count) (t-1)",
                             "fatal_mid_count"      = "Fatal MID count (t-1)",
                             "log_fatal_mid_count"  = "log(1 + Fatal MID count) (t-1)"),
         se_ratio_hac_over_ols = nw_se / se) |>
  select(war_variable, iv_label, estimate = coef,
         ols_se = se, ols_p = pval,
         hc1_se = any_of("hc1_se"), hc1_pval = any_of("hc1_pval"),
         hc3_se, hc3_pval, nw_se, nw_pval,
         se_ratio_hac_over_ols, n_obs)

if (nrow(hac_appendix_rows) > 0) {
  write_csv(hac_appendix_rows,
            file.path(OUTPUT_DIR, "appendix_hac_comparison.csv"))
  cat(sprintf("  Saved: appendix_hac_comparison.csv (%d rows)\n",
              nrow(hac_appendix_rows)))

  fmt_hac <- function(x, d = 4) formatC(x, format = "f", digits = d)
  fmt_p_hac <- function(p) {
    ifelse(p < 0.001, "$<\\!0.001$",
           paste0("$", formatC(p, format = "f", digits = 3), "$"))
  }
  hac_tex_rows <- hac_appendix_rows |>
    transmute(
      IV = iv_label,
      `Estimate` = fmt_hac(estimate),
      `OLS SE`   = fmt_hac(ols_se),
      `OLS p`    = fmt_p_hac(ols_p),
      `HAC SE`   = fmt_hac(nw_se),
      `HAC p`    = fmt_p_hac(nw_pval),
      `HAC/OLS ratio` = fmt_hac(se_ratio_hac_over_ols, 2),
      `N` = n_obs
    )
  hac_body <- apply(hac_tex_rows, 1, function(r) paste(r, collapse = " & "))
  hac_tex_lines <- c(
    "\\begin{table}[t]", "\\centering", "\\small",
    "\\begin{tabular}{lrrrrrrr}", "\\toprule",
    paste("IV (t-1) & Estimate & OLS SE & OLS $p$ & HAC SE & HAC $p$ & HAC/OLS ratio & $N$ \\\\"),
    "\\midrule",
    paste0(hac_body, " \\\\"),
    "\\bottomrule", "\\end{tabular}",
    paste0("\\caption{Appendix: system-level baseline with Newey-West HAC ",
           "standard errors (lag 5). Headline cell of the spec curve ",
           "(log(1+KL) DV, raw transform, unweighted aggregation, ",
           "alliance + IGO layers). Point estimates identical to the ",
           "main-text OLS baseline; only standard errors differ. Ratio ",
           "$> 1$ indicates serial correlation that i.i.d. SEs ",
           "understate. The full 4-way SE ladder (iid / HC1 / HC3 / ",
           "NW HAC) populates the spec-curve $-\\log_{10}(p)$ ",
           "distribution. Two-sided $p$-values.}"),
    "\\label{tab:appendix-hac}", "\\end{table}"
  )
  writeLines(hac_tex_lines, file.path(APPX_DIR, "appendix_hac_table.tex"))
  cat("  Saved: results/appendix/appendix_hac_table.tex\n")
} else {
  cat("  WARNING: no headline-cell rows found; appendix HAC table not written.\n")
  cat("  Check that HEADLINE_LAYERS matches the layers naming used in specs_df.\n")
}

###############################################################################
# tmp_multiverse_appendix_fig.R   [TEMPORARY / standalone]
#
# Regenerates the multiverse appendix figures from the SAVED specification data
# (data/processed/robustness_specifications.csv) WITHOUT rerunning script 07.
#
# Produces (to code/03_figures/output/):
#   fig_multiverse_hypotheses.png   -- (left) heatmap: share significant-positive
#                                      by decision x hypothesis; (right) overall
#                                      share supporting each hypothesis.
#   fig_multiverse_meta_by_hyp.png  -- coefficient plot: how each design choice
#                                      shifts the downstream estimate, faceted
#                                      across H1 / H2 / H3 (free scales).
#
# Conventions:
#   * Significance = ONE-SIDED POSITIVE at p < 0.10  (coef > 0 & pval/2 < 0.10).
#   * DCA is dropped: there are NO DCA specs in the 07 multiverse (has_dca is
#     always FALSE), which was the source of the NA cell in the earlier figure.
#   * H1 = system; H2 = state, own conflict; H3 = state, partner exposure.
###############################################################################

library(tidyverse)
library(patchwork)
set.seed(123)

OUTPUT_DIR <- here::here("data", "processed")
FIG_DIR    <- here::here("code", "03_figures", "output")
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

specs <- readr::read_csv(file.path(OUTPUT_DIR, "robustness_specifications.csv"),
                         show_col_types = FALSE)

# Python wrote "True"/"False"; coerce the design booleans to logical to be safe.
bool_cols <- c("has_offensive", "has_defensive", "has_igo", "has_trade",
               "has_dca", "iv_is_fatal", "iv_is_log", "dv_logged")
specs <- specs |> mutate(across(any_of(bool_cols), ~ as.logical(.x)))

# ===========================================================================
# 1. Classify specs into hypotheses; flag one-sided positive significance
# ===========================================================================
PARTNER <- c("partner_war_frac", "partner_war_frac_weighted",
             "partner_fatal_war_frac", "partner_fatal_war_frac_weighted")

specs <- specs |>
  filter(!str_starts(model, "VAR_")) |>           # drop VAR / Granger / FEVD cells
  filter(!is.na(coef), !is.na(pval)) |>
  mutate(
    hypothesis = case_when(
      level == "system"                            ~ "H1",
      level == "state" & war_variable %in% PARTNER ~ "H3",
      level == "state"                             ~ "H2"
    ),
    sig_pos = coef > 0 & (pval / 2) < 0.10          # one-sided positive, p < 0.10
  ) |>
  filter(!is.na(hypothesis))

HYP <- c("H1", "H2", "H3")

# ===========================================================================
# 2. Readable labels for each decision value
# ===========================================================================
dv_measure <- c(
  churn = "Churn", log_churn = "Churn",
  churn_alt = "Churn (alternate)", log_churn_alt = "Churn (alternate)",
  degree_kl = "Degree KL", log_degree_kl = "Degree KL",
  degree_wasserstein = "Degree Wasserstein", log_degree_wasserstein = "Degree Wasserstein",
  edge_turnover_pct = "Edge turnover", log_edge_turnover_pct = "Edge turnover",
  js_divergence = "JS divergence", log_js_divergence = "JS divergence",
  kl_divergence = "KL divergence", log_kl_divergence = "KL divergence"
)
weight_lab <- c(
  contig_euclidean = "Contiguity (Euclidean)", contig_manhattan = "Contiguity (Manhattan)",
  dist_euclidean = "Distance (Euclidean)", dist_manhattan = "Distance (Manhattan)",
  unweighted = "Unweighted"
)
transform_lab <- c(log1p = "log(1 + x)", raw = "Raw", zscore = "z-score")

# Display order (top -> bottom) within each category.
ITEM_ORDER <- c(
  # (a) reconfiguration measure
  "Churn", "Churn (alternate)", "Degree KL", "Degree Wasserstein",
  "Edge turnover", "JS divergence", "KL divergence",
  # (b) network layer
  "Defensive alliances", "Offensive alliances", "Shared IGOs", "Trade",
  # (c) transform
  "log(1 + x)", "Raw", "z-score",
  # (d) conflict measure
  "All MIDs", "Fatal MIDs",
  # (e) edge weighting
  "Contiguity (Euclidean)", "Contiguity (Manhattan)",
  "Distance (Euclidean)", "Distance (Manhattan)", "Unweighted"
)
CAT_A <- "(a) Reconfig-\nuration measure"
CAT_B <- "(b) Network layer\nincluded"
CAT_C <- "(c) Variable\ntransform"
CAT_D <- "(d) Conflict\nmeasure"
CAT_E <- "(e) Edge\nweighting"

pct <- function(x) 100 * mean(x)   # share TRUE over a non-empty logical vector

# ===========================================================================
# 3. Heatmap data: % significant-positive by (decision value, hypothesis)
# ===========================================================================
d_recon <- specs |>
  mutate(item = dv_measure[dv]) |>
  group_by(hypothesis, item) |>
  summarise(val = pct(sig_pos), .groups = "drop") |>
  mutate(category = CAT_A)

layer_map <- tibble(
  col  = c("has_offensive", "has_defensive", "has_igo", "has_trade"),
  item = c("Offensive alliances", "Defensive alliances", "Shared IGOs", "Trade")
)
d_layer <- purrr::pmap_dfr(layer_map, function(col, item) {
  specs |> filter(.data[[col]]) |>
    group_by(hypothesis) |>
    summarise(val = pct(sig_pos), .groups = "drop") |>
    mutate(item = item)
}) |> mutate(category = CAT_B)

d_trans <- specs |>
  mutate(item = transform_lab[transform]) |>
  group_by(hypothesis, item) |>
  summarise(val = pct(sig_pos), .groups = "drop") |>
  mutate(category = CAT_C)

d_conf <- specs |>
  mutate(item = if_else(iv_is_fatal, "Fatal MIDs", "All MIDs")) |>
  group_by(hypothesis, item) |>
  summarise(val = pct(sig_pos), .groups = "drop") |>
  mutate(category = CAT_D)

d_weight <- specs |>
  mutate(item = weight_lab[weighting]) |>
  group_by(hypothesis, item) |>
  summarise(val = pct(sig_pos), .groups = "drop") |>
  mutate(category = CAT_E)

heat <- bind_rows(d_recon, d_layer, d_trans, d_conf, d_weight) |>
  # Fill the (item x hypothesis) grid so undefined combinations become NA (dash),
  # NOT a mislabeled numeric. This is what removes the stray "NA" cell.
  complete(nesting(category, item), hypothesis = HYP) |>
  mutate(
    category   = factor(category, levels = c(CAT_A, CAT_B, CAT_C, CAT_D, CAT_E)),
    item       = factor(item, levels = rev(ITEM_ORDER)),   # first listed on top
    hypothesis = factor(hypothesis, levels = HYP),
    label      = ifelse(is.na(val), "—", sprintf("%.0f", val))
  )

p_heat <- ggplot(heat, aes(x = hypothesis, y = item, fill = val)) +
  geom_tile(color = "white", linewidth = 0.6) +
  geom_text(aes(label = label), size = 3.1, color = "black") +
  facet_grid(category ~ ., scales = "free_y", space = "free_y", switch = "y") +
  scale_fill_gradient2(low = "#C0392B", mid = "white", high = "#2C6CB0",
                       midpoint = 50, limits = c(0, 100),
                       na.value = "grey90", name = "% significant positive") +
  scale_x_discrete(position = "top") +
  labs(title = "Which specification choices matter?",
       subtitle = "Share of specifications significant-positive, by decision. Dashes = measure not defined at that level.",
       x = NULL, y = NULL) +
  guides(fill = guide_colorbar(barwidth = 12, barheight = 0.5, title.position = "top")) +
  theme_minimal(base_size = 11) +
  theme(panel.grid = element_blank(),
        strip.text.y.left = element_text(angle = 90, face = "bold", size = 8),
        strip.placement = "outside",
        legend.position = "bottom",
        axis.text = element_text(color = "black"),
        plot.title = element_text(face = "bold"),
        plot.subtitle = element_text(color = "gray30", size = 8.5))

# ===========================================================================
# 4. Overall bars: share of all specs supporting each hypothesis
# ===========================================================================
overall <- specs |>
  group_by(hypothesis) |>
  summarise(sig = 100 * mean(sig_pos), .groups = "drop") |>
  mutate(null = 100 - sig) |>
  pivot_longer(c(sig, null), names_to = "outcome", values_to = "share") |>
  mutate(
    outcome = factor(outcome, levels = c("sig", "null"),
                     labels = c("Significant positive (one-sided p < 0.1)",
                                "Null or significant negative")),
    hypothesis = factor(hypothesis, levels = rev(HYP))   # H1 at top
  )

p_bar <- ggplot(overall, aes(x = share, y = hypothesis, fill = outcome)) +
  geom_col(width = 0.62) +
  geom_text(aes(label = sprintf("%.1f%%", share)),
            position = position_stack(vjust = 0.5), color = "white",
            fontface = "bold", size = 3.4) +
  scale_fill_manual(values = c("Significant positive (one-sided p < 0.1)" = "#1F2C74",
                               "Null or significant negative" = "#8B1E1E"),
                    name = NULL) +
  scale_x_continuous(labels = function(x) paste0(x, "%"),
                     expand = expansion(mult = c(0, 0))) +
  labs(title = "Overall multiverse result",
       subtitle = "Share of all specifications supporting each hypothesis.",
       x = "Share of specifications", y = NULL) +
  theme_minimal(base_size = 11) +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "bottom",
        legend.direction = "vertical",
        axis.text = element_text(color = "black"),
        plot.title = element_text(face = "bold"),
        plot.subtitle = element_text(color = "gray30", size = 8.5))

p_fig1 <- p_heat + p_bar + plot_layout(widths = c(1.35, 1))
ggsave(file.path(FIG_DIR, "fig_multiverse_hypotheses.png"), p_fig1,
       width = 13, height = 8.5, dpi = 400, bg = "white")
cat("  -> fig_multiverse_hypotheses.png\n")

# ===========================================================================
# 5. Coefficient plot: how each decision shifts the downstream estimate,
#    faceted across H1 / H2 / H3 (free scales).
#    Outcome = the spec's coefficient, z-scored WITHIN hypothesis, so the
#    meta-regression coefficients read as SD-unit shifts in the downstream
#    estimate and are comparable within each facet.
# ===========================================================================
meta_df <- specs |>
  group_by(hypothesis) |>
  mutate(coef_z = as.numeric(scale(coef))) |>
  ungroup() |>
  mutate(conflict = if_else(iv_is_fatal, "Fatal MIDs", "All MIDs"),
         recon    = dv_measure[dv],
         # Layer choices recoded as DEVIATIONS from the main spec, which includes
         # offensive + defensive + IGO and excludes trade. Each indicator is
         # FALSE in the main spec, so the main layer set is the baseline.
         drop_offensive = !has_offensive,
         drop_defensive = !has_defensive,
         drop_igo       = !has_igo,
         add_trade      = has_trade,
         neglog10p      = -log10(pmax(pval, 1e-300)),   # significance, continuous
         sig_num        = as.integer(sig_pos))          # 0/1 for the LPM

# Reference = the MAIN specification for each hypothesis, so the intercept IS the
# main model and every coefficient reads as a deviation from it:
#   H1 (system) : OLS_lagged, JS divergence, raw, unweighted, main layers, All MIDs
#   H2/H3 (state): FE_full_ccode_year (state+year TWFE = "model 4"), Churn, raw,
#                  unweighted, main layers, All MIDs
relevel_to <- function(x, ref) {
  x <- droplevels(factor(x))
  if (ref %in% levels(x)) relevel(x, ref = ref) else x
}

# Wald tidy: estimate +/- 1.96*SE (fast; avoids profiling confint() on 100k rows).
tidy_wald <- function(m) {
  ct <- summary(m)$coefficients
  tibble(term = rownames(ct), estimate = ct[, 1], se = ct[, 2]) |>
    filter(term != "(Intercept)") |>
    mutate(conf.low = estimate - 1.96 * se, conf.high = estimate + 1.96 * se) |>
    select(term, estimate, conf.low, conf.high)
}

run_meta <- function(df, hyp) {
  ref_model   <- if (hyp == "H1") "OLS_lagged" else "FE_full_ccode_year"
  ref_measure <- if (hyp == "H1") "JS divergence" else "Churn"  # JS = main-text DV (2026-09-26)
  
  df <- df |> mutate(
    model_f   = relevel_to(model, ref_model),
    transform = relevel_to(transform, "raw"),
    weighting = relevel_to(weighting, "unweighted"),
    conflict  = relevel_to(conflict, "All MIDs"),
    recon     = relevel_to(recon, ref_measure)
  )
  # Keep only predictors that actually vary within this hypothesis.
  keep <- function(nm, x) if (length(unique(x[!is.na(x)])) > 1) nm else NULL
  preds <- c(keep("model_f", df$model_f),
             keep("transform", df$transform), keep("weighting", df$weighting),
             keep("conflict", df$conflict), keep("recon", df$recon),
             keep("drop_offensive", df$drop_offensive),
             keep("drop_defensive", df$drop_defensive),
             keep("drop_igo", df$drop_igo), keep("add_trade", df$add_trade))
  rhs <- paste(preds, collapse = " + ")
  
  # THREE downstream outcomes, same predictors, same main-spec baseline:
  #   (1) the estimate (z-scored within hypothesis),
  #   (2) statistical significance as -log10(p),
  #   (3) P(significant-positive) via an LPM (avoids logit separation; reads as a
  #       change in the share of significant-positive specs).
  m_coef <- lm(as.formula(paste("coef_z ~",    rhs)), data = df)
  m_p    <- lm(as.formula(paste("neglog10p ~", rhs)), data = df)
  m_sig  <- lm(as.formula(paste("sig_num ~",   rhs)), data = df)
  
  bind_rows(
    tidy_wald(m_coef) |> mutate(outcome = "Coefficient (SD)"),
    tidy_wald(m_p)    |> mutate(outcome = "-log10(p-value)"),
    tidy_wald(m_sig)  |> mutate(outcome = "P(significant-positive)")
  )
}

model_labels <- c(
  sys_T1_ar1          = "Model: AR(1)",
  sys_T1_iv_only      = "Model: IV only",
  sys_T1_order_period = "Model: order/period",
  LP_h5               = "Model: Local projection (h5)",
  FE_full_year_only   = "Model: Year FE only",
  FE_full_ccode_only  = "Model: Unit FE only",
  full_noFE           = "Model: No FE (full controls)",
  full_noFE_year_poly = "Model: No FE + year poly",
  baseline_noFE       = "Model: Baseline (no FE)"
)

clean_term <- function(t) dplyr::case_when(
  str_starts(t, "model_f")         ~ unname(model_labels[str_remove(t, "^model_f")]),
  t == "transformlog1p"            ~ "Transform: log(1+x)",
  t == "transformzscore"           ~ "Transform: z-score",
  t == "weightingcontig_manhattan" ~ "Weight: Contiguity (Manhattan)",
  t == "weightingcontig_euclidean" ~ "Weight: Contiguity (Euclidean)",
  t == "weightingdist_manhattan"   ~ "Weight: Distance (Manhattan)",
  t == "weightingdist_euclidean"   ~ "Weight: Distance (Euclidean)",
  t == "conflictFatal MIDs"        ~ "Conflict: Fatal MIDs",
  t == "drop_offensiveTRUE"        ~ "Drop offensive alliances",
  t == "drop_defensiveTRUE"        ~ "Drop defensive alliances",
  t == "drop_igoTRUE"              ~ "Drop shared IGOs",
  t == "add_tradeTRUE"             ~ "Add trade",
  str_starts(t, "recon")           ~ paste0("Measure: ", str_remove(t, "^recon")),
  TRUE ~ t
)

meta_coef <- purrr::imap_dfr(split(meta_df, meta_df$hypothesis),
                             ~ run_meta(.x, .y) |> mutate(hypothesis = .y)) |>
  mutate(term_clean = clean_term(term))

# Order y-axis terms by decision category (model -> transform -> weighting ->
# conflict -> measure -> layers), consistent across all facets.
cat_rank <- function(tc) dplyr::case_when(
  str_starts(tc, "Model:")     ~ 1,
  str_starts(tc, "Transform:") ~ 2,
  str_starts(tc, "Weight:")    ~ 3,
  str_starts(tc, "Conflict:")  ~ 4,
  str_starts(tc, "Measure:")   ~ 5,
  TRUE                         ~ 6            # layer deviations (Drop.../Add...)
)
term_levels <- meta_coef |> distinct(term_clean) |>
  mutate(r = cat_rank(term_clean)) |> arrange(r, term_clean) |> pull(term_clean)

meta_coef <- meta_coef |>
  mutate(hypothesis = factor(hypothesis, levels = HYP),
         term_clean = factor(term_clean, levels = rev(term_levels)),
         outcome    = factor(outcome, levels = c("Coefficient (SD)",
                                                 "-log10(p-value)",
                                                 "P(significant-positive)")))

# facet_grid free scales: x free per column (the 3 outcomes have very different
# scales), y free per row (each hypothesis shows only its own terms).
p_fig2 <- ggplot(meta_coef, aes(x = estimate, y = term_clean)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray55") +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), height = 0.2,
                 color = "#2C6CB0") +
  geom_point(size = 1.9, color = "#2C6CB0") +
  facet_grid(hypothesis ~ outcome, scales = "free", space = "free_y") +
  labs(title = "How each measurement decision affects the downstream analyses",
       subtitle = paste0("Meta-regressions of the estimate, its significance (-log10 p), ",
                         "and P(significant-positive) on design choices,\nrelative to the ",
                         "MAIN specification (H1: OLS lagged-DV; H2/H3: state+year TWFE). ",
                         "Rows = hypotheses; columns = outcomes (free x)."),
       x = "Effect vs. main specification", y = NULL) +
  theme_bw(base_size = 10) +
  theme(panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "gray92", color = "gray70"),
        strip.text = element_text(face = "bold"),
        axis.text.y = element_text(size = 7.5),
        axis.text.x = element_text(size = 8),
        plot.title = element_text(face = "bold"),
        plot.subtitle = element_text(color = "gray30", size = 8.5))

ggsave(file.path(FIG_DIR, "fig_multiverse_meta_by_hyp.png"), p_fig2,
       width = 13, height = 11, dpi = 400, bg = "white")
cat("  -> fig_multiverse_meta_by_hyp.png\n")

cat("\nDone. Overall support by hypothesis (one-sided positive, p<0.10):\n")
print(specs |> group_by(hypothesis) |>
        summarise(n = n(), pct_sig_pos = round(100 * mean(sig_pos), 1)))


cat("\nDone. 03_robustness.R complete.\n")

