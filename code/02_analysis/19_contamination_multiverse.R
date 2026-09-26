###############################################################################
# 20_contamination_multiverse.R
#
# MULTIVERSE FOR THE CONTAMINATION RESULT (appendix)
# --------------------------------------------------
# 03_contamination_demo.R produces the main-text figure on ONE measurement
# choice. This script re-runs that entire design across the same measurement
# grid used by 07_robustness_multiverse.R, plus the two degrees of freedom that
# are specific to the contamination design.
#
# WHY THIS IS A SEPARATE SCRIPT
#   It is slow (~2-4 h). Folding it into 03 would make every run of the
#   pipeline pay that cost just to redraw Figure 8. Same reasoning as
#   19_robustness_excl_belligerent.R.
#
# THE OUTCOME IS NOT THE COEFFICIENT
#   A coefficient distribution would only say "the estimate moves around."
#   The paper's claim is sharper: THE SIGN OF THE TIE EQUATION PREDICTS THE
#   DIRECTION THE LADDER MOVES. If prior conflict weakens a tie, conditioning
#   on it should pull the network coefficient TOWARD zero (attenuation). If
#   prior conflict strengthens a tie, conditioning should push the coefficient
#   AWAY from zero (amplification).
#
#   So every cell records `prediction_held`, and the headline statistic is the
#   share of specifications in which the prediction was correct. That converts
#   a 6-for-6 result in the main text into an N-for-N result here.
#
# THE GRID
#   layer combos     24  (all subsets of offensive/defensive/IGO/trade/DCA
#                         requiring IGO or defensive)
#   weighting         5  (contig/dist x manhattan/euclidean, unweighted)
#   transform         3  (raw, log1p, zscore; zscore x euclidean skipped)
#   history window    3  (t-2..t-4, t-2..t-6 [main], t-2..t-11)
#   monadic spec      3  (sum, max, min+max)
#
# DCA IS DROPPED [2026-07-16, author request] for consistency with the rest of
#   the paper, which excludes DCA throughout. Layer inputs are therefore
#   offensive, defensive, IGO, and trade (4 singles) plus one all-layers
#   combined tie = 5 inputs x 3 weightings = 15 tie-cells (30 dyadic+monadic).
#
# YEAR-RANGE RESTRICTION
#   Every spec is estimated only on years where ALL of its layers exist. This
#   mirrors 07 and prevents the structural-zero problem: a late-coverage layer
#   contributing zeros before its start date would leave the tie IDENTICAL to
#   the combo without it, making the "robustness" spec a duplicate.
#
# NO GDP OR STRENGTH CONTROLS HERE.
#   GDP is a state-year attribute; the dyadic literature's answer to relative
#   material capacity is the capability ratio, which is already in CTRL.
#   Network strength is worse: it is the sum of a state's tie weights at t-1,
#   and TIE STRENGTH AT t-1 IS THE TREATMENT. Conditioning on it would absorb
#   the treatment's own variation and make the attenuation meaningless.
#
# Outputs:
#   data/processed/contamination_multiverse.csv
#   code/03_figures/output/manuscript/fig_ms_contamination_multiverse.png
###############################################################################

library(tidyverse)
library(fixest)
library(peacesciencer)
library(furrr)
library(future)
select <- dplyr::select
set.seed(123)

OUTPUT_DIR <- here::here("data", "processed")
FIG_DIR    <- here::here("code", "03_figures", "output", "manuscript")
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

YR_MIN <- 1816L
YR_MAX <- 2014L

# --------------------------------------------------------------------------- #
# CONFIG -- dial the grid down here, do not edit the loop
# --------------------------------------------------------------------------- #
# REDUCED GRID [2026-07-15, author request]: this is no longer a full multiverse.
# We hold the transform, history window, and monadic parameterisation at their
# main-text values and vary only (i) the layer input -- each layer on its own
# plus ONE all-layers combined tie -- and (ii) three tie measures: unweighted
# Manhattan, contiguity-weighted Manhattan, and distance-weighted Manhattan.
MV_WEIGHTING  <- c("unweighted", "contig_manhattan", "dist_manhattan")
MV_TRANSFORMS <- c("raw")
MV_HIST <- list("t2_t6" = 2:6)   # main-text window only
MV_MONADIC <- c("sum")           # main-text monadic (third-party sum, see sec 4)

MAIN_HIST    <- "t2_t6"      # the window reported in the main text
MAIN_MONADIC <- "sum"

# Workers. EACH worker serialises dyad_data (~1.9 GB). Peak RAM ~ workers x 2GB.
.n_workers <- as.integer(Sys.getenv("CONTAM_WORKERS", "4"))
plan(multisession, workers = .n_workers)
options(future.globals.maxSize = 100 * 1024^3)
cat(sprintf("Parallel plan: %d workers\n", .n_workers))

# ===========================================================================
# 1. Dyad-year panel (same construction as 03)
# ===========================================================================
cat("Building dyad-year panel...\n")

d <- create_dyadyears(subset_years = YR_MIN:YR_MAX, directed = FALSE) |>
  add_cow_mids()          |>
  add_democracy()         |>
  add_nmc()               |>
  add_cow_majors()        |>
  add_peace_years()

py_col <- intersect(c("cowmidspell", "gmlmidspell", "midspell"), names(d))[1]
stopifnot(!is.na(py_col))
d$peace_years <- d[[py_col]]

d <- d |>
  mutate(dyad_id = paste0(pmin(ccode1, ccode2), "_", pmax(ccode1, ccode2)))

# Dyadic covariates. NOTE: contiguity and log distance are DROPPED -- they are
# time-invariant within a dyad (absorbed by dyad FE), and whatever within-dyad
# variation they do have is driven by war-driven territorial change, which makes
# them post-treatment relative to the conflict history being conditioned on.
d <- d |>
  mutate(
    # Dyadic capability asymmetry, bounded [0,1]; matches 03/20/21.
    cap_ratio       = abs(cinc1 - cinc2) / (cinc1 + cinc2),
    joint_democracy = as.integer(pmin(polity21, polity22, na.rm = TRUE) >= 6),
    major_power     = as.integer(cowmaj1 == 1 | cowmaj2 == 1)
  )

CTRL <- "cap_ratio + joint_democracy + major_power +
         peace_years + I(peace_years^2) + I(peace_years^3)"

# ===========================================================================
# 2. Network layers + the helpers borrowed from 07
# ===========================================================================
dyad_data <- readRDS(file.path(OUTPUT_DIR, "dyad_year_network_layers.rds"))
metadata  <- readRDS(file.path(OUTPUT_DIR, "network_metadata.rds"))

all_layer_vars <- metadata$layers
all_layer_vars <- all_layer_vars[sapply(all_layer_vars, function(lv) {
  sd(dyad_data[[lv]], na.rm = TRUE) > 0
})]
# DCA dropped for consistency across the paper [2026-07-16, author request].
all_layer_vars <- setdiff(all_layer_vars, "layer_dca")
cat(sprintf("Layers: %s\n", paste(all_layer_vars, collapse = ", ")))

# Layer coverage, for the year-range restriction.
# COVERAGE = where the layer's DATA EXISTS, not where it is positive
# [FIXED 2026-07-15; see 07 for the full rationale]. Offensive alliances have
# their last positive tie in 1990 but the data span the full panel -- the
# post-1990 zeros are REAL, not missing -- so inferring coverage from
# `layer > 0` would wrongly discard 24 years of valid zeros. Trade begins 1870
# and DCA 1980 as genuine data-start boundaries. Hard-coded from the sources:
LAYER_COVERAGE <- list(
  layer_offensive = c(1816L, YR_MAX),   # ATOP; real zeros post-1990
  layer_defensive = c(1816L, YR_MAX),   # ATOP
  layer_igo       = c(1816L, YR_MAX),   # Pevehouse COW IGO
  layer_trade     = c(1870L, YR_MAX),   # COW trade begins 1870
  layer_dca       = c(1980L, 2010L)     # Kinne DCA, 1980-2010
)

layer_year_ranges <- list()
for (lv in all_layer_vars) {
  layer_year_ranges[[lv]] <- if (!is.null(LAYER_COVERAGE[[lv]])) {
    LAYER_COVERAGE[[lv]]
  } else {
    ay <- dyad_data |> filter(.data[[lv]] > 0) |> pull(year)
    if (length(ay)) c(min(ay), max(ay)) else c(NA, NA)
  }
  cat(sprintf("  %-18s coverage %d-%d\n", lv,
              layer_year_ranges[[lv]][1], layer_year_ranges[[lv]][2]))
}
get_valid_year_range <- function(layers) {
  c(max(sapply(layers, function(lv) layer_year_ranges[[lv]][1])),
    min(sapply(layers, function(lv) layer_year_ranges[[lv]][2])))
}

range01 <- function(x) {
  rng <- range(x, na.rm = TRUE)
  if (!is.finite(rng[2]) || rng[2] <= 0) return(rep(0, length(x)))
  x / rng[2]
}

# Identical to 07's compute_manhattan_vec(), copied rather than sourced so this
# script does not inherit 07's ~3,400 lines of side effects.
compute_manhattan_vec <- function(data, layers, method) {
  mat <- as.matrix(data[, layers, drop = FALSE]); mat[is.na(mat)] <- 0
  if (method == "contig_manhattan") {
    wc <- intersect(paste0(layers, "_weighted"), names(data))
    if (!length(wc)) return(rep(0, nrow(data)))
    m <- as.matrix(data[, wc, drop = FALSE]); m[is.na(m)] <- 0
    result <- rowSums(m) / length(wc)
  } else if (method == "contig_euclidean") {
    wc <- intersect(paste0(layers, "_weighted"), names(data))
    if (!length(wc)) return(rep(0, nrow(data)))
    m <- as.matrix(data[, wc, drop = FALSE]); m[is.na(m)] <- 0; m[m < 0] <- 0
    result <- sqrt(rowSums(m^2)) / length(wc)
  } else if (method == "dist_manhattan") {
    dc <- intersect(paste0(layers, "_dist_weighted"), names(data))
    if (!length(dc)) return(rep(0, nrow(data)))
    m <- as.matrix(data[, dc, drop = FALSE]); m[is.na(m)] <- 0
    result <- rowSums(m) / length(dc)
  } else if (method == "dist_euclidean") {
    dc <- intersect(paste0(layers, "_dist_weighted"), names(data))
    if (!length(dc)) return(rep(0, nrow(data)))
    m <- as.matrix(data[, dc, drop = FALSE]); m[is.na(m)] <- 0; m[m < 0] <- 0
    result <- sqrt(rowSums(m^2)) / length(dc)
  } else if (method == "unweighted") {
    result <- rowSums(mat)
  }
  years <- data$year
  for (yr in unique(years)) {
    idx <- which(years == yr)
    result[idx] <- range01(result[idx])
  }
  result
}

# Layer inputs [REDUCED 2026-07-15]: each layer on its own, plus ONE combined
# tie built from ALL layers together. DCA retained (this is the one script it
# belongs in -- see header). No IGO/defensive filter: every single layer runs.
all_combos <- c(
  lapply(all_layer_vars, function(lv) lv),   # each layer alone
  list(sort(all_layer_vars))                 # one all-layers combined tie
)
cat(sprintf("Layer inputs: %d singles + 1 combined = %d\n",
            length(all_layer_vars), length(all_combos)))

n_cells <- length(all_combos) * length(MV_WEIGHTING) * length(MV_TRANSFORMS) *
           length(MV_HIST) * length(MV_MONADIC)
cat(sprintf("Grid: %d layer-inputs x %d weighting x %d transform x %d window x %d monadic\n",
            length(all_combos), length(MV_WEIGHTING), length(MV_TRANSFORMS),
            length(MV_HIST), length(MV_MONADIC)))
cat(sprintf("      = %s cells\n", format(n_cells, big.mark = ",")))

# ===========================================================================
# 3. Monadic (own-conflict) history, per state-year
# ===========================================================================
# Built once per history window and reused across every layer cell -- it does
# not depend on the network measure.
state_churn <- readRDS(file.path(OUTPUT_DIR, "state_network_churn.rds"))

# Matrix-free windowed lag sum. rowSums(sapply(...)) collapses to a plain
# vector -- and then errors -- for single-row groups (a state or dyad present
# for one year only, e.g. ccode 511). Reduce over a list is size-safe.
lag_sum <- function(x, ks) {
  Reduce(`+`, lapply(ks, function(k) dplyr::lag(x, k, default = 0)))
}

build_state_hist <- function(lags) {
  state_churn |>
    select(ccode, year, mid_onset_count) |>
    arrange(ccode, year) |>
    group_by(ccode) |>
    mutate(mid_hist = lag_sum(mid_onset_count, lags)) |>
    ungroup() |>
    select(ccode, year, mid_hist)
}
STATE_HIST <- lapply(MV_HIST, build_state_hist)

# Dyadic history: did i and j fight EACH OTHER during the window?
build_dyad_hist <- function(dat, lags) {
  dat |>
    arrange(dyad_id, year) |>
    group_by(dyad_id) |>
    mutate(dyad_hist = lag_sum(cowmidonset, lags)) |>
    ungroup()
}
DYAD_HIST <- lapply(MV_HIST, function(lags) {
  build_dyad_hist(d |> select(dyad_id, ccode1, ccode2, year, cowmidonset), lags) |>
    select(dyad_id, year, dyad_hist)
})

# ===========================================================================
# 4. One cell = (combo, weighting, transform, window, monadic)
# ===========================================================================
# All lags use default = 0, so the estimation sample does NOT change with the
# history window -- which lets rung (a) be estimated once per tie-cell and
# reused across windows. Rungs (b) and (c) must be refit per window.

run_tie_cell <- function(layers_subset, method, tfm, dyad_work) {

  combo_name <- paste(gsub("layer_", "", layers_subset), collapse = "+")
  yr <- get_valid_year_range(layers_subset)

  tie_vec <- compute_manhattan_vec(dyad_work, layers_subset, method)

  ties <- tibble(
    ccode1 = pmin(dyad_work$ccode1, dyad_work$ccode2),
    ccode2 = pmax(dyad_work$ccode1, dyad_work$ccode2),
    year   = dyad_work$year,
    tie    = tie_vec
  ) |>
    group_by(ccode1, ccode2, year) |>
    summarise(tie = max(tie, na.rm = TRUE), .groups = "drop")

  base <- d |>
    filter(year >= yr[1], year <= yr[2]) |>          # <- coverage restriction
    left_join(ties, by = c("ccode1", "ccode2", "year")) |>
    mutate(tie = replace_na(tie, 0)) |>
    arrange(dyad_id, year) |>
    group_by(dyad_id) |>
    mutate(tie_lag = dplyr::lag(tie, 1)) |>
    ungroup()

  out <- list()

  for (hw in names(MV_HIST)) {

    dat <- base |>
      left_join(DYAD_HIST[[hw]], by = c("dyad_id", "year")) |>
      left_join(STATE_HIST[[hw]] |> rename(mh_i = mid_hist),
                by = c("ccode1" = "ccode", "year")) |>
      left_join(STATE_HIST[[hw]] |> rename(mh_j = mid_hist),
                by = c("ccode2" = "ccode", "year")) |>
      mutate(
        dyad_hist = replace_na(dyad_hist, 0),
        mh_i = replace_na(mh_i, 0), mh_j = replace_na(mh_j, 0),
        # THIRD-PARTY monadic sum (mirrors 03's mid5_sum_tp): each state's own
        # onsets over the window MINUS its i-j onsets, so the monadic term is
        # conflict with THIRD parties only and does not double-count the dyadic
        # process. Clamp at 0 in case the monadic/dyadic counts differ at the margin.
        mon_sum = pmax(mh_i - dyad_hist, 0) + pmax(mh_j - dyad_hist, 0),
        mon_max = pmax(mh_i, mh_j),
        mon_min = pmin(mh_i, mh_j)
      )

    need <- c("cowmidonset", "tie_lag", "cap_ratio", "joint_democracy",
              "major_power", "peace_years", "dyad_hist")
    est <- dat |> filter(if_all(all_of(need), ~ !is.na(.x)))

    # At-risk sample: a dyad already in an ongoing MID cannot have a NEW onset.
    est <- est |> filter(coalesce(cowmidongoing, 0L) == 0 | cowmidonset == 1)
    if (nrow(est) < 5000) { rm(dat, est); gc(verbose = FALSE); next }

    sd_safe <- function(x) { s <- sd(x, na.rm = TRUE); if (is.na(s) || s == 0) 1 else s }
    est <- est |>
      mutate(
        tie_sd       = tie_lag  / sd_safe(tie_lag),
        # Conflict-history onset counts: log1p FIRST (mirrors the main analyses'
        # log1p(onset count)) THEN per-SD, for both the dyadic and monadic terms.
        dyad_hist_sd = { z <- log1p(dyad_hist); z / sd_safe(z) },
        mon_sum_sd   = { z <- log1p(mon_sum);   z / sd_safe(z) },
        mon_max_sd   = mon_max  / sd_safe(mon_max),
        mon_min_sd   = mon_min  / sd_safe(mon_min)
      )
    if (sd(est$tie_sd, na.rm = TRUE) == 0) { rm(dat, est); gc(verbose = FALSE); next }

    # ---- rung (a): baseline, no conflict history. Same across windows, but
    #      the sample can differ by coverage so it is refit inside the loop.
    m_a <- tryCatch(
      feols(as.formula(paste("cowmidonset ~ tie_sd +", CTRL, "| dyad_id + year")),
            data = est, vcov = ~dyad_id),
      error = function(e) NULL)
    if (is.null(m_a) || !"tie_sd" %in% rownames(coeftable(m_a))) { rm(dat, est); next }
    beta_base <- coeftable(m_a)["tie_sd", "Estimate"]
    p_base    <- coeftable(m_a)["tie_sd", "Pr(>|t|)"]
    rm(m_a)   # only the scalars are needed downstream; free the model frame now
    if (!is.finite(beta_base) || beta_base == 0) { rm(dat, est); next }

    for (ms in MV_MONADIC) {

      mon_term <- switch(ms,
        "sum"     = "mon_sum_sd",
        "max"     = "mon_max_sd",
        "min_max" = "mon_max_sd")   # only "sum" runs in the reduced grid

      # TWO SEPARATE PROCESSES, each conditioned on its OWN history term and
      # compared to the same baseline. Emitted as separate rows (process col),
      # so the figure can show the dyadic and monadic ladders side by side.
      #   Dyadic  = did i and j fight EACH OTHER?      (i-j MIDs)
      #   Monadic = did i or j fight THIRD parties?    (own onsets minus i-j)
      PROC <- c(Dyadic = "dyad_hist_sd", Monadic = mon_term)

      for (pnm in names(PROC)) {
        hist_term <- PROC[[pnm]]

        # ---- Step 1: does THIS process's history reshape the tie?
        m_tie <- tryCatch(
          feols(as.formula(paste("tie_lag ~", hist_term, "+", CTRL,
                                 "| dyad_id + year")),
                data = est, vcov = ~dyad_id),
          error = function(e) NULL)
        if (is.null(m_tie) || !hist_term %in% rownames(coeftable(m_tie))) {
          if (!is.null(m_tie)) rm(m_tie); next
        }
        beta_tie <- coeftable(m_tie)[hist_term, "Estimate"]
        p_tie    <- coeftable(m_tie)[hist_term, "Pr(>|t|)"]

        # ---- Step 2: conditioned ladder for this process only.
        m_c <- tryCatch(
          feols(as.formula(paste("cowmidonset ~ tie_sd +", CTRL, "+",
                                 hist_term, "| dyad_id + year")),
                data = est, vcov = ~dyad_id),
          error = function(e) NULL)
        if (is.null(m_c) || !"tie_sd" %in% rownames(coeftable(m_c))) { rm(m_tie); next }
        beta_cond <- coeftable(m_c)["tie_sd", "Estimate"]
        p_cond    <- coeftable(m_c)["tie_sd", "Pr(>|t|)"]

        # ---- leg 2: does this history predict onset, net of the tie?
        beta_leg2 <- coeftable(m_c)[hist_term, "Estimate"]
        p_leg2    <- coeftable(m_c)[hist_term, "Pr(>|t|)"]
        rm(m_tie, m_c)   # scalars extracted; free both model frames

        # ---- THE PREDICTION.
        # bias ~ (history -> onset) x (history <-> tie). With history -> onset
        # positive (leg 2), the SIGN of the tie equation determines the sign of
        # the omitted-variable bias, and therefore which way conditioning moves
        # the coefficient.
        #   tie equation NEGATIVE (war severs the tie)  -> attenuate (toward 0)
        #   tie equation POSITIVE (war builds the tie)  -> amplify  (away from 0)
        attenuation <- 100 * (1 - abs(beta_cond) / abs(beta_base))
        predicted   <- if (beta_tie < 0) "attenuate" else "amplify"
        observed    <- if (abs(beta_cond) < abs(beta_base)) "attenuate" else "amplify"

        out[[length(out) + 1]] <- tibble(
          layers = combo_name, n_layers = length(layers_subset),
          weighting = method, transform = tfm,
          history_window = hw, monadic_spec = ms, process = pnm,
          year_start = yr[1], year_end = yr[2], n_obs = nrow(est),
          has_offensive = "layer_offensive" %in% layers_subset,
          has_defensive = "layer_defensive" %in% layers_subset,
          has_igo       = "layer_igo"       %in% layers_subset,
          has_trade     = "layer_trade"     %in% layers_subset,
          has_dca       = "layer_dca"       %in% layers_subset,
          beta_tie_eq = beta_tie, p_tie_eq = p_tie,
          sign_tie_eq = if (beta_tie < 0) "negative" else "positive",
          beta_leg2 = beta_leg2, p_leg2 = p_leg2,
          beta_base = beta_base, p_base = p_base,
          beta_cond = beta_cond, p_cond = p_cond,
          attenuation_pct = attenuation,
          predicted_dir = predicted, observed_dir = observed,
          prediction_held = predicted == observed
        )
      }  # end process loop
    }  # end monadic loop

    # Free the per-window frames before the next window. Without this the
    # worker's heap climbs across combos until it is OOM-killed.
    rm(dat, est)
    gc(verbose = FALSE)
  }  # end window loop

  rm(base, ties, tie_vec)
  gc(verbose = FALSE)
  dplyr::bind_rows(out)
}

# ===========================================================================
# 5. Run the grid
# ===========================================================================
results <- list()
t0 <- proc.time()

for (tfm in MV_TRANSFORMS) {

  # Layer transform, applied before weighting (mirrors 07).
  if (tfm == "log1p") {
    cont <- intersect(c("layer_igo", "layer_trade"), all_layer_vars)
    dyad_work <- dyad_data
    for (cl in cont) dyad_work[[cl]] <- log1p(dyad_work[[cl]])
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

  for (method in MV_WEIGHTING) {
    if (tfm == "zscore" && grepl("euclidean", method)) next   # mirrors 07

    cat(sprintf("  [%s / %s] %d combos on %d workers...\n",
                tfm, method, length(all_combos), .n_workers))

    chunk <- furrr::future_map_dfr(
      all_combos,
      function(cb) run_tie_cell(cb, method, tfm, dyad_work),
      .options = furrr_options(seed = 123)
    )
    results[[length(results) + 1]] <- chunk
    cat(sprintf("      -> %d cells (%.1f min elapsed)\n",
                nrow(chunk), (proc.time() - t0)["elapsed"] / 60))
    rm(chunk)
    gc(verbose = FALSE)
  }
  # dyad_work is a full copy of the ~1M-row panel; drop it before the next
  # transform builds a fresh one.
  rm(dyad_work)
  gc(verbose = FALSE)

  # Respawn the parallel workers. gc() reclaims memory INSIDE each R worker, but
  # multisession workers do not reliably return freed memory to the OS, so the
  # heap creeps across batches until one worker is OOM-killed (the "Future
  # interrupted" error). Tearing down and rebuilding the plan resets the workers
  # to a clean slate between transforms.
  plan(sequential)
  plan(multisession, workers = .n_workers)
}

mv <- bind_rows(results)
cat(sprintf("\nTOTAL CELLS: %s   (%.1f min)\n",
            format(nrow(mv), big.mark = ","),
            (proc.time() - t0)["elapsed"] / 60))

write_csv(mv, file.path(OUTPUT_DIR, "contamination_multiverse.csv"))

# ===========================================================================
# 6. THE HEADLINE
# ===========================================================================
cat("\n=====================================================================\n")
cat("HEADLINE: does the tie equation's sign predict the ladder's movement?\n")
cat("=====================================================================\n")
cat(sprintf("\n  Prediction held in %.1f%% of %s specifications.\n",
            100 * mean(mv$prediction_held), format(nrow(mv), big.mark = ",")))

cat("\n--- by process x tie-equation sign ---\n")
print(mv |> group_by(process, sign_tie_eq) |>
        summarise(n = n(), held = 100 * mean(prediction_held),
                  median_atten = median(attenuation_pct), .groups = "drop"), n = Inf)

cat("\n--- leg 2 by process: is this history a significant predictor of onset? ---\n")
print(mv |> group_by(process) |>
        summarise(pct_sig05 = 100 * mean(p_leg2 < 0.05),
                  pct_pos   = 100 * mean(beta_leg2 > 0), .groups = "drop"), n = Inf)

# Readable label for each layer input, used in the console breakdown and figure.
label_layer <- function(x) dplyr::case_when(
  x == "offensive" ~ "Offensive\nalliance",
  x == "defensive" ~ "Defensive\nalliance",
  x == "igo"       ~ "IGO",
  x == "trade"     ~ "Trade",
  x == "dca"       ~ "DCA",
  grepl("\\+", x)  ~ "Combined\n(all layers)",
  TRUE             ~ x
)
weight_label <- c("unweighted"      = "Manhattan (unweighted)",
                  "contig_manhattan" = "Contiguity-weighted Manhattan",
                  "dist_manhattan"   = "Distance-weighted Manhattan")

cat("\n--- attenuation by process x layer input (across the 3 Manhattan measures) ---\n")
print(mv |>
        mutate(layer_type = gsub("\n", " ", label_layer(layers))) |>
        group_by(process, layer_type) |>
        summarise(n = n(),
                  median_atten = median(attenuation_pct),
                  pct_amplify  = 100 * mean(attenuation_pct < 0),
                  held         = 100 * mean(prediction_held), .groups = "drop"), n = Inf)

cat("\n--- every cell, for reference ---\n")
print(mv |>
        mutate(input = gsub("\n", " ", label_layer(layers))) |>
        select(process, input, weighting, beta_tie_eq, beta_base, beta_cond,
               attenuation_pct, prediction_held), n = Inf)

# ===========================================================================
# 7. FIGURE
# ===========================================================================
cat("\n--- Plotting ---\n")

ACCENT <- "#2C6CB0"

# Order the layer inputs: singles first, combined last (top of the plot).
input_order <- c("Offensive\nalliance", "Defensive\nalliance", "IGO",
                 "Trade", "Combined\n(all layers)")
# One row, two columns: dyadic (i-j MIDs) on the left, monadic (third-party) on
# the right -- mirroring the two separate conflict-history rungs in the main text.
proc_label <- c(Dyadic  = "Dyadic conflict history\n(i-j MIDs)",
                Monadic = "Monadic conflict history\n(third-party, i-j removed)")
plot_df <- mv |>
  mutate(layer_type = factor(label_layer(layers), levels = rev(input_order)),
         measure    = factor(weight_label[weighting],
                             levels = unname(weight_label)),
         process    = factor(proc_label[process], levels = unname(proc_label))) |> 
  filter(!is.na(layer_type), layer_type != "DCA")

p <- ggplot(plot_df, aes(x = attenuation_pct, y = layer_type,
                         color = measure, shape = measure)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray55") +
  geom_point(size = 3, stroke = 0.9,
             position = position_dodge(width = 0.5)) +
  facet_wrap(~ process, nrow = 1, ncol = 2, scales = "free_x") +
  scale_color_manual(name = NULL,
                     values = c(ACCENT, "#C2A25A", "#6BA368")) +
  scale_shape_manual(name = NULL, values = c(16, 17, 15)) +
  scale_x_continuous(breaks = scales::pretty_breaks(4)) +
  labs(
    title = "Contamination is not an artifact of one measurement choice",
    subtitle = paste0(
      "Attenuation of the tie-conflict coefficient by network input, under three\n",
      "Manhattan tie measures, for the dyadic and monadic conflict-history ",
      "processes separately.\nPositive = attenuation (toward zero); negative = ",
      "amplification (away from zero)."),
    x = "Attenuation (%)", y = NULL
  ) +
  theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(),
        strip.background = element_rect(fill = "gray92", color = "gray70"),
        strip.text = element_text(face = "bold"),
        axis.text.x = element_text(angle = 0, hjust = 0.5),
        axis.text.y = element_text(lineheight = 0.9),
        legend.position = "bottom",
        plot.title = element_text(face = "bold"),
        plot.subtitle = element_text(color = "gray30", size = 9))

ggsave(file.path(FIG_DIR, "fig_ms_contamination_multiverse.png"), p,
       width = 11, height = 5.5, dpi = 400, bg = "white")
cat("  -> fig_ms_contamination_multiverse.png\n")

cat("\nSaved -> data/processed/contamination_multiverse.csv\n")
plan(sequential)
