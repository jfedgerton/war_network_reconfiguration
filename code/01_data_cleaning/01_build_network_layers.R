###############################################################################
# 01_build_network_layers.R
# Construct yearly network layers for Part 1: War → Network Reconfiguration
#
# Layers:
#   1. Offensive alliances (ATOP v5.1, 1815-2018)
#   2. Defensive alliances (ATOP v5.1, 1815-2018)
#   3. Shared IGO membership (COW via peacesciencer)
#   4. Bilateral trade volume (COW via peacesciencer)
#   5. Defense Cooperation Agreements (Kinne — external)
#
# Weighting: Contiguity-adjusted Manhattan distance (Edgerton 2024, JOP)
# Alternative weighting schemes as robustness
###############################################################################

library(tidyverse)
library(peacesciencer)
library(data.table)
select <- dplyr::select

# --------------------------------------------------------------------------- #
# 0. Parameters
# --------------------------------------------------------------------------- #

YEAR_MIN <- 1816
YEAR_MAX <- 2014
OUTPUT_DIR <- here::here("data", "processed")
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

# --------------------------------------------------------------------------- #
# 1. Build dyad-year base from peacesciencer
# --------------------------------------------------------------------------- #

cat("Building dyad-year base...\n")

dyad_base <- create_dyadyears(subset_years = YEAR_MIN:YEAR_MAX) |>
  # Add state attributes
  add_cow_majors() |>
  add_contiguity() |>
  add_capital_distance() |>
  add_nmc()

cat(sprintf("Base dyad-years: %s rows\n", format(nrow(dyad_base), big.mark = ",")))

# --------------------------------------------------------------------------- #
# 2. Add network layers
# --------------------------------------------------------------------------- #

cat("Adding alliance data (ATOP v5.1, 1815-2018)...\n")
dyad_data <- dyad_base |>
  add_atop_alliance()

cat("Adding IGO co-membership...\n")
dyad_data <- dyad_data |>
  add_igos()

cat("Adding bilateral trade...\n")
dyad_data <- dyad_data |>
  add_cow_trade()

# --------------------------------------------------------------------------- #
# 3. Construct layer-specific edge variables
# --------------------------------------------------------------------------- #

cat("Constructing network layers...\n")

# ATOP alliance types (v5.1):
#   atop_offense  = 1 if offense pledge
#   atop_defense  = 1 if defense pledge
#   atop_neutrality = 1 if neutrality pledge
#   atop_nonagg  = 1 if nonaggression pledge
#   atop_consul  = 1 if consultation pledge
#
# Mapping to layers:
#   layer_offensive = atop_offense (attack obligations)
#   layer_defensive = atop_defense (defense obligations)

dyad_data <- dyad_data |>
  mutate(
    # Alliance layers (binary from ATOP)
    layer_offensive = as.numeric(atop_offense == 1),
    layer_defensive = as.numeric(atop_defense == 1),

    # IGO co-membership (continuous — number of shared IGOs, normalized 0-1 by year)
    layer_igo = dyadigos,

    # Trade (continuous — total bilateral trade flow, smoothed for coverage)
    layer_trade = smoothflow1 + smoothflow2,
    layer_trade = ifelse(is.na(layer_trade) | layer_trade <= 0, 0, layer_trade),

    # Contiguity (for weighting)
    contiguous = as.numeric(conttype %in% 1:2)
  )

# Normalize IGO within year (0-1)
dyad_data <- dyad_data |>
  group_by(year) |>
  mutate(
    layer_igo = if (max(layer_igo, na.rm = TRUE) > 0) {
      (layer_igo - min(layer_igo, na.rm = TRUE)) /
        (max(layer_igo, na.rm = TRUE) - min(layer_igo, na.rm = TRUE))
    } else {
      0
    },
    # Normalize trade within year (0-1)
    layer_trade = if (max(layer_trade, na.rm = TRUE) > 0) {
      (layer_trade - min(layer_trade, na.rm = TRUE)) /
        (max(layer_trade, na.rm = TRUE) - min(layer_trade, na.rm = TRUE))
    } else {
      0
    }
  ) |>
  ungroup()

# Replace NAs with 0 in layer variables
layer_vars <- c("layer_offensive", "layer_defensive", "layer_igo", "layer_trade")
dyad_data <- dyad_data |>
  mutate(across(all_of(layer_vars), ~ replace_na(.x, 0)))

# --------------------------------------------------------------------------- #
# 4. Merge Kinne DCA data
# --------------------------------------------------------------------------- #

cat("Loading Kinne DCA data...\n")

# Kinne's DCA data: defense cooperation agreements
# Source: Brandon Kinne (https://brandonkinne.com)
# File should be placed in data/raw/kinne_dca/
DCA_PATH <- here::here("data", "raw", "kinne_dca")
dir.create(DCA_PATH, recursive = TRUE, showWarnings = FALSE)

cat(sprintf("  DCA_PATH resolved to: %s\n", DCA_PATH))
cat(sprintf("  DCA_PATH exists: %s\n", dir.exists(DCA_PATH)))

dca_files <- list.files(DCA_PATH, pattern = "\\.csv$|\\.dta$|\\.rds$", full.names = TRUE)

# Prefer the dyadic file if multiple CSVs exist
dyadic_file <- grep("dyadic", dca_files, value = TRUE, ignore.case = TRUE)
if (length(dyadic_file) > 0) {
  dca_file <- dyadic_file[1]
} else if (length(dca_files) > 0) {
  dca_file <- dca_files[1]
} else {
  dca_file <- NULL
}

if (!is.null(dca_file)) {
  cat("  Found DCA file(s):", basename(dca_files), "\n")
  cat("  Using:", basename(dca_file), "\n")

  # Read DCA data
  if (grepl("\\.csv$", dca_file)) {
    dca_raw <- read_csv(dca_file, show_col_types = FALSE)
  } else if (grepl("\\.dta$", dca_file)) {
    dca_raw <- haven::read_dta(dca_file)
  } else if (grepl("\\.rds$", dca_file)) {
    dca_raw <- readRDS(dca_file)
  }

  # Standardize DCA to dyad-year format
  # Kinne's data uses COW codes; merge on ccode1, ccode2, year
  # Use dcaAnyV2 (any DCA, version 2) per JOP reference
  cat("  DCA columns: ", paste(names(dca_raw), collapse = ", "), "\n")

  if ("ccode1" %in% names(dca_raw) & "ccode2" %in% names(dca_raw)) {
    # Explicitly select dcaAnyV2 (matches JOP reference which renames to dcaV2)
    dca_col <- if ("dcaAnyV2" %in% names(dca_raw)) {
      "dcaAnyV2"
    } else if ("dcaanyv2" %in% tolower(names(dca_raw))) {
      names(dca_raw)[tolower(names(dca_raw)) == "dcaanyv2"][1]
    } else {
      # Fallback: first dca column
      grep("dca", names(dca_raw), value = TRUE, ignore.case = TRUE)[1]
    }
    cat(sprintf("  Using DCA column: %s\n", dca_col))

    dca_clean <- dca_raw |>
      select(ccode1, ccode2, year, all_of(dca_col)) |>
      rename(layer_dca = !!sym(dca_col)) |>
      mutate(
        layer_dca = as.numeric(layer_dca),
        # Ensure directed dyad ordering matches peacesciencer (ccode1 < ccode2)
        low = pmin(ccode1, ccode2),
        high = pmax(ccode1, ccode2)
      ) |>
      distinct(low, high, year, .keep_all = TRUE)

    # Merge into main data
    dyad_data <- dyad_data |>
      mutate(
        low = pmin(ccode1, ccode2),
        high = pmax(ccode1, ccode2)
      ) |>
      left_join(
        dca_clean |> select(low, high, year, layer_dca),
        by = c("low", "high", "year")
      ) |>
      mutate(layer_dca = replace_na(layer_dca, 0)) |>
      select(-low, -high)

    layer_vars <- c(layer_vars, "layer_dca")
    cat("  DCA layer merged successfully.\n")
    cat(sprintf("  DCA diagnostics: %d nonzero rows (%.2f%%), years %d-%d\n",
                sum(dyad_data$layer_dca > 0, na.rm = TRUE),
                100 * mean(dyad_data$layer_dca > 0, na.rm = TRUE),
                min(dyad_data$year[dyad_data$layer_dca > 0], na.rm = TRUE),
                max(dyad_data$year[dyad_data$layer_dca > 0], na.rm = TRUE)))
    cat(sprintf("  DCA clean had %d nonzero rows before merge\n",
                sum(dca_clean$layer_dca > 0, na.rm = TRUE)))
  } else {
    warning("DCA file structure not recognized. Columns found: ",
            paste(names(dca_raw), collapse = ", "))
    dyad_data$layer_dca <- 0
    layer_vars <- c(layer_vars, "layer_dca")
  }
} else {
  cat("  WARNING: No DCA file found in", DCA_PATH, "\n")
  cat("  Looked for .csv/.dta/.rds files\n")
  cat("  Download Kinne's DCA data from https://brandonkinne.com/data\n")
  cat("  Place the dyadic file (DCAD-v1.0-dyadic.csv) in data/raw/kinne_dca/\n")
  cat("  Proceeding without DCA layer.\n")
  dyad_data$layer_dca <- 0
  layer_vars <- c(layer_vars, "layer_dca")
}

# --------------------------------------------------------------------------- #
# 5. Contiguity weighting (Edgerton 2024, JOP method)
# --------------------------------------------------------------------------- #
#
# For each network layer, regress the layer on contiguity with sender FE,
# receiver FE, and year FE. Take residuals. This strips out geographic
# proximity so distant interactions carry more weight.
#
# Formula: y_{ijk} = β₁ contiguity + η_i + γ_j + ζ_k
#   where η_i = sender FE, γ_j = receiver FE, ζ_k = year FE
#
# The residual (observed - predicted) captures the portion of the tie NOT
# explained by contiguity. Positive residuals → more interaction than
# geographic proximity alone would predict.
# --------------------------------------------------------------------------- #

# --------------------------------------------------------------------------- #
# Verify multiplex alignment: all layers share the same dyad-year universe
# --------------------------------------------------------------------------- #

cat("Verifying multiplex alignment...\n")

alignment_check <- dyad_data |>
  summarise(
    total_rows = n(),
    n_dyads = n_distinct(paste(ccode1, ccode2)),
    n_years = n_distinct(year),
    year_range = paste(min(year), max(year), sep = "-"),
    across(all_of(layer_vars), ~ sum(!is.na(.x)), .names = "nonmissing_{.col}"),
    across(all_of(layer_vars), ~ mean(.x > 0, na.rm = TRUE), .names = "pct_active_{.col}")
  )

cat("  Alignment summary:\n")
cat(sprintf("    Total dyad-years: %s\n", format(alignment_check$total_rows, big.mark = ",")))
cat(sprintf("    Unique dyads: %s\n", format(alignment_check$n_dyads, big.mark = ",")))
cat(sprintf("    Year range: %s\n", alignment_check$year_range))
cat("    Layer coverage:\n")
for (lv in layer_vars) {
  nm <- paste0("nonmissing_", lv)
  pct <- paste0("pct_active_", lv)
  cat(sprintf("      %-20s: %s non-missing (%.1f%% active)\n",
              lv, format(alignment_check[[nm]], big.mark = ","),
              alignment_check[[pct]] * 100))
}

# Check for any year gaps
year_seq <- min(dyad_data$year):max(dyad_data$year)
year_counts <- dyad_data |> count(year) |> pull(year)
missing_years <- setdiff(year_seq, year_counts)
if (length(missing_years) > 0) {
  cat(sprintf("  WARNING: Missing years in panel: %s\n", paste(missing_years, collapse = ", ")))
} else {
  cat("  No missing years in panel.\n")
}

cat("Applying contiguity weighting (JOP predict method)...\n")

# range01: normalize to [0,1] (from JOP reference code)
range01 <- function(x) {
  rng <- max(x, na.rm = TRUE) - min(x, na.rm = TRUE)
  if (rng == 0) return(rep(0, length(x)))
  (x - min(x, na.rm = TRUE)) / rng
}

# For each layer: feols(value ~ contiguous | year + ccode1 + ccode2)
# Then: tie = value - predict(model)  [JOP method, avoids residuals() naming bug]
for (lv in layer_vars) {
  weighted_name <- paste0(lv, "_weighted")
  cat(sprintf("  Weighting %s...\n", lv))

  layer_vals <- dyad_data[[lv]]

  # Skip constant layers
  if (all(is.na(layer_vals)) || sd(layer_vals, na.rm = TRUE) == 0) {
    cat(sprintf("    Skipping %s (constant/all-zero)\n", lv))
    dyad_data[[weighted_name]] <- 0
    next
  }

  # Filter to years with variation (JOP lines 218-222)
  year_var <- dyad_data |>
    group_by(year) |>
    summarise(n_unique = n_distinct(.data[[lv]]), .groups = "drop")
  min_year_with_var <- min(year_var$year[year_var$n_unique > 1])

  temp_data <- dyad_data |>
    filter(year >= min_year_with_var) |>
    select(ccode1, ccode2, year, contiguous, all_of(lv))

  # Run FE regression
  formula_str <- paste0(lv, " ~ contiguous | year + ccode1 + ccode2")
  fit <- tryCatch(
    fixest::feols(as.formula(formula_str), data = temp_data,
                  warn = FALSE, notes = FALSE),
    error = function(e) {
      cat(sprintf("    FE regression failed for %s: %s\n", lv, e$message))
      NULL
    }
  )

  if (is.null(fit)) {
    dyad_data[[weighted_name]] <- 0
    next
  }

  # Use predict() — JOP method (avoids residuals naming bug)
  temp_data$pred_values <- predict(fit)
  temp_data$tie <- temp_data[[lv]] - temp_data$pred_values

  cat(sprintf("    %s: n_obs=%d, tie range=[%.4f, %.4f], pct>0=%.1f%%\n",
              lv, fit$nobs, min(temp_data$tie, na.rm = TRUE),
              max(temp_data$tie, na.rm = TRUE),
              100 * mean(temp_data$tie > 0, na.rm = TRUE)))

  # Merge back to full data
  tie_data <- temp_data |> select(ccode1, ccode2, year, tie)
  dyad_data <- dyad_data |>
    left_join(tie_data, by = c("ccode1", "ccode2", "year"), suffix = c("", ".tmp"))
  dyad_data[[weighted_name]] <- replace_na(dyad_data$tie, 0)
  dyad_data$tie <- NULL
  if ("tie.tmp" %in% names(dyad_data)) dyad_data$tie.tmp <- NULL
}

weighted_layer_vars <- paste0(layer_vars, "_weighted")

# --------------------------------------------------------------------------- #
# 5b. Distance weighting (capital distance, log(1+d))
# --------------------------------------------------------------------------- #
#
# Same FE regression approach as contiguity, but using log(1 + capital distance)
# as the geographic predictor instead of binary contiguity.
#
# Formula: y_{ijk} = β₁ log(1 + capdist) + η_i + γ_j + ζ_k
# --------------------------------------------------------------------------- #

cat("Applying distance weighting (log1p capital distance)...\n")

# Check capdist availability
capdist_na <- sum(is.na(dyad_data$capdist))
cat(sprintf("  capdist: %d NAs out of %d (%.1f%%)\n",
            capdist_na, nrow(dyad_data), 100 * capdist_na / nrow(dyad_data)))

for (lv in layer_vars) {
  dist_weighted_name <- paste0(lv, "_dist_weighted")
  cat(sprintf("  Distance-weighting %s...\n", lv))

  layer_vals <- dyad_data[[lv]]

  # Skip constant layers
  if (all(is.na(layer_vals)) || sd(layer_vals, na.rm = TRUE) == 0) {
    cat(sprintf("    Skipping %s (constant/all-zero)\n", lv))
    dyad_data[[dist_weighted_name]] <- 0
    next
  }

  # Filter to years with variation
  year_var <- dyad_data |>
    group_by(year) |>
    summarise(n_unique = n_distinct(.data[[lv]]), .groups = "drop")
  min_year_with_var <- min(year_var$year[year_var$n_unique > 1])

  temp_data <- dyad_data |>
    filter(year >= min_year_with_var, !is.na(capdist)) |>
    mutate(log_capdist = log1p(capdist)) |>
    select(ccode1, ccode2, year, log_capdist, all_of(lv))

  # Run FE regression with log(1 + capdist)
  formula_str <- paste0(lv, " ~ log_capdist | year + ccode1 + ccode2")
  fit <- tryCatch(
    fixest::feols(as.formula(formula_str), data = temp_data,
                  warn = FALSE, notes = FALSE),
    error = function(e) {
      cat(sprintf("    FE regression failed for %s: %s\n", lv, e$message))
      NULL
    }
  )

  if (is.null(fit)) {
    dyad_data[[dist_weighted_name]] <- 0
    next
  }

  # Use predict() — same approach as contiguity weighting
  temp_data$pred_values <- predict(fit)
  temp_data$tie <- temp_data[[lv]] - temp_data$pred_values

  cat(sprintf("    %s: n_obs=%d, tie range=[%.4f, %.4f], pct>0=%.1f%%\n",
              lv, fit$nobs, min(temp_data$tie, na.rm = TRUE),
              max(temp_data$tie, na.rm = TRUE),
              100 * mean(temp_data$tie > 0, na.rm = TRUE)))

  # Merge back to full data
  tie_data <- temp_data |> select(ccode1, ccode2, year, tie)
  dyad_data <- dyad_data |>
    left_join(tie_data, by = c("ccode1", "ccode2", "year"), suffix = c("", ".tmp"))
  dyad_data[[dist_weighted_name]] <- replace_na(dyad_data$tie, 0)
  dyad_data$tie <- NULL
  if ("tie.tmp" %in% names(dyad_data)) dyad_data$tie.tmp <- NULL
}

dist_weighted_layer_vars <- paste0(layer_vars, "_dist_weighted")

# --------------------------------------------------------------------------- #
# 6. Compute aggregate distance measures (2×2 + unweighted)
# --------------------------------------------------------------------------- #
#
# 5 aggregation variants:
#   1. manhattan (unweighted, Manhattan norm) — baseline
#   2. manhattan_weighted (contiguity, Manhattan norm) — JOP method
#   3. manhattan_contig_euclid (contiguity, Euclidean norm)
#   4. manhattan_dist (distance, Manhattan norm)
#   5. manhattan_dist_euclid (distance, Euclidean norm)
# --------------------------------------------------------------------------- #

cat("Computing aggregate distance measures (2x2 + unweighted)...\n")

n_layers <- length(layer_vars)

# --- 6a. Unweighted Manhattan (baseline) ---
raw_mat <- as.matrix(dyad_data[, layer_vars, drop = FALSE])
raw_mat[is.na(raw_mat)] <- 0
dyad_data$manhattan_unweighted_raw <- rowSums(raw_mat)

dyad_data <- dyad_data |>
  group_by(year) |>
  mutate(
    manhattan = if (max(manhattan_unweighted_raw, na.rm = TRUE) > 0) {
      manhattan_unweighted_raw / max(manhattan_unweighted_raw, na.rm = TRUE)
    } else {
      0
    }
  ) |>
  ungroup() |>
  mutate(manhattan = round(manhattan, 2))

cat(sprintf("  Unweighted manhattan: pct > 0 = %.1f%%\n",
            100 * mean(dyad_data$manhattan > 0, na.rm = TRUE)))

# --- 6b. Contiguity-weighted Manhattan (JOP method) ---
manhattan_mat <- as.matrix(dyad_data[, weighted_layer_vars, drop = FALSE])
manhattan_mat[is.na(manhattan_mat)] <- 0
dyad_data$manhattan_weighted_raw <- rowSums(manhattan_mat) / n_layers

dyad_data <- dyad_data |>
  group_by(year) |>
  mutate(manhattan_weighted = range01(manhattan_weighted_raw)) |>
  ungroup() |>
  mutate(manhattan_weighted = round(manhattan_weighted, 2))

cat(sprintf("  Contiguity Manhattan: pct > 0 = %.1f%%\n",
            100 * mean(dyad_data$manhattan_weighted > 0, na.rm = TRUE)))

# --- 6c. Contiguity-weighted Euclidean ---
# Zero out negative residuals so only positive deviations (more interaction
# than geography predicts) contribute. Without this, squaring turns large
# negative residuals into large positive values, making dissimilar dyads
# appear similar.
manhattan_mat_pos <- manhattan_mat
manhattan_mat_pos[manhattan_mat_pos < 0] <- 0
dyad_data$contig_euclid_raw <- sqrt(rowSums(manhattan_mat_pos^2)) / n_layers

dyad_data <- dyad_data |>
  group_by(year) |>
  mutate(manhattan_contig_euclid = range01(contig_euclid_raw)) |>
  ungroup() |>
  mutate(manhattan_contig_euclid = round(manhattan_contig_euclid, 2))

cat(sprintf("  Contiguity Euclidean: pct > 0 = %.1f%%\n",
            100 * mean(dyad_data$manhattan_contig_euclid > 0, na.rm = TRUE)))

# --- 6d. Distance-weighted Manhattan ---
dist_mat <- as.matrix(dyad_data[, dist_weighted_layer_vars, drop = FALSE])
dist_mat[is.na(dist_mat)] <- 0
dyad_data$dist_manhattan_raw <- rowSums(dist_mat) / n_layers

dyad_data <- dyad_data |>
  group_by(year) |>
  mutate(manhattan_dist = range01(dist_manhattan_raw)) |>
  ungroup() |>
  mutate(manhattan_dist = round(manhattan_dist, 2))

cat(sprintf("  Distance Manhattan: pct > 0 = %.1f%%\n",
            100 * mean(dyad_data$manhattan_dist > 0, na.rm = TRUE)))

# --- 6e. Distance-weighted Euclidean ---
# Same fix: zero out negative residuals before squaring
dist_mat_pos <- dist_mat
dist_mat_pos[dist_mat_pos < 0] <- 0
dyad_data$dist_euclid_raw <- sqrt(rowSums(dist_mat_pos^2)) / n_layers

dyad_data <- dyad_data |>
  group_by(year) |>
  mutate(manhattan_dist_euclid = range01(dist_euclid_raw)) |>
  ungroup() |>
  mutate(manhattan_dist_euclid = round(manhattan_dist_euclid, 2))

cat(sprintf("  Distance Euclidean: pct > 0 = %.1f%%\n",
            100 * mean(dyad_data$manhattan_dist_euclid > 0, na.rm = TRUE)))

# --------------------------------------------------------------------------- #
# 7. Add war / MID data
# --------------------------------------------------------------------------- #

cat("Adding MID data...\n")

# Restore peacesciencer attribute (gets stripped by mutate/data.table operations)
attr(dyad_data, "ps_data_type") <- "dyad_year"
attr(dyad_data, "ps_system") <- "cow"

# add_cow_mids() adds cowmidongoing and cowmidonset
# Use keep= to retain fatality columns from cow_mid_dirdisps
dyad_data <- dyad_data |>
  add_cow_mids(keep = c("fatality1", "fatality2", "hostlev1", "hostlev2"))

# Create fatal MID indicator: fatality1 or fatality2 > 0
dyad_data <- dyad_data |>
  mutate(
    cowmidfatal = as.integer(
      cowmidongoing == 1 &
      (coalesce(fatality1, 0L) > 0 | coalesce(fatality2, 0L) > 0)
    )
  )

# System-level war counts per year (two DVs for analysis)
system_war_counts <- dyad_data |>
  group_by(year) |>
  summarise(
    mid_count = sum(cowmidongoing == 1, na.rm = TRUE),
    fatal_mid_count = sum(cowmidfatal == 1, na.rm = TRUE),
    mid_onset_count = sum(cowmidonset == 1, na.rm = TRUE),
    fatal_mid_onset_count = sum(cowmidfatal == 1 & cowmidonset == 1, na.rm = TRUE),
    .groups = "drop"
  )

# --------------------------------------------------------------------------- #
# 8. Save outputs
# --------------------------------------------------------------------------- #

cat("Saving outputs...\n")

# Full dyad-year data with all layers
saveRDS(dyad_data, file.path(OUTPUT_DIR, "dyad_year_network_layers.rds"))

# System-level war counts
saveRDS(system_war_counts, file.path(OUTPUT_DIR, "system_war_counts.rds"))

# Metadata
metadata <- list(
  layers = layer_vars,
  weighted_layers = weighted_layer_vars,
  dist_weighted_layers = dist_weighted_layer_vars,
  year_range = c(YEAR_MIN, YEAR_MAX),
  weighting = "contiguity_and_distance_fe",
  n_dyad_years = nrow(dyad_data),
  created = Sys.time()
)
saveRDS(metadata, file.path(OUTPUT_DIR, "network_metadata.rds"))

cat("Done. Outputs saved to", OUTPUT_DIR, "\n")
cat(sprintf("  Layers: %s\n", paste(layer_vars, collapse = ", ")))
cat(sprintf("  Dyad-years: %s\n", format(nrow(dyad_data), big.mark = ",")))
