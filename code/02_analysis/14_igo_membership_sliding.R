# [2026-09-26] Headline system-level DV switched from KL to Jensen-Shannon (JS)
# divergence. KL results are retained in results/appendix/kl_canonical/.
###############################################################################
# 11_igo_membership_sliding.R
# Robustness: rebuild the multiplex network dropping IGOs by per-year
# membership share, then refit the preliminary analyses.
#
# Motivation: the main IGO layer lumps all IGOs together (a dyad-year count
# of shared memberships). A reviewer might worry that results are driven by
# (a) niche / regional IGOs with very few members, or (b) near-universal
# IGOs like the UN that almost everyone joins. Both are plausible
# nuisances. This script lets us drop either tail and see whether the
# main-text findings survive.
#
# Requires raw Pevehouse IGO V3 state-unit data (IGO x state x year panel).
# Obtained from CoW or the igoR package. See `load_igo_state_year()` below.
#
# Three sliding grids:
#   (1) Lower tail    : drop IGOs with < L % membership share;
#                       L in {0, 5, 10, 15}
#   (2) Upper tail    : drop IGOs with > U % membership share;
#                       U in {85, 90, 95, 100}
#   (3) Both tails    : lower L AND upper U; L x U grid
#
# At each threshold combination:
#   - Rebuild dyad-year shared IGO count from the filtered IGO set
#   - Normalize layer_igo within year (0-1)
#   - Combine with alliance layer -> new main_weight
#   - Recompute year-over-year KL divergence (system) and
#     weighted Jaccard churn (state)
#   - Refit the main system-level OLS and state-level TWFE
#   - Store fatal / non-fatal MID coefficients and SEs
#
# Outputs:
#   data/processed/robustness_igo_share_results.csv
#   code/03_figures/output/presentation/
#       pres_fig15_igo_share_lower.png
#       pres_fig15_igo_share_upper.png
#       pres_fig15_igo_share_both.png
###############################################################################

library(tidyverse)
library(fixest)
library(igraph)
library(peacesciencer)
select <- dplyr::select

set.seed(123)

OUTPUT_DIR <- here::here("data", "processed")
RAW_DIR    <- here::here("data", "raw", "igo")
FIG_DIR    <- here::here("code", "03_figures", "output", "manuscript")
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(RAW_DIR, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# Manuscript theme — white background, publication-suitable
# ============================================================

ACCENT  <- "#2C6CB0"
ACCENT2 <- "#B03366"

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

save_pres <- function(filename, plot, width = 7, height = 5, dpi = 400) {
  ggsave(file.path(FIG_DIR, filename), plot,
         width = width, height = height, dpi = dpi, bg = "white")
  cat(sprintf("  -> %s\n", filename))
}

# ============================================================
# IGO membership loader — tries igoR, falls back to raw CSV
# ============================================================
#
# Expected return: tibble with columns
#   igo_name  (char or factor)  -- unique IGO identifier
#   ccode     (int)             -- CoW country code
#   year      (int)
#   member    (0/1)             -- 1 if full member
#
# igoR convention: the long-form state-year panel is exposed as an object
# named something like `state_year_format3` or accessed via a function.
# The loader below handles the most common two access patterns; if neither
# works it reads the Pevehouse CSV from data/raw/igo/.
#
# If either path produces a tibble with different column names, edit the
# final `select()` call to map onto (igo_name, ccode, year, member).

load_igo_state_year <- function() {
  # --- Option 1: igoR package ---
  if (requireNamespace("igoR", quietly = TRUE)) {
    cat("Using igoR package...\n")
    # igoR exposes the Pevehouse V3 state-unit panel as `state_year_format3`
    # (lazy-loaded; needs an explicit data() call to materialize).
    # If igoR's actual object name differs, edit this block.
    e <- new.env()
    utils::data(list = "state_year_format3", package = "igoR",
                envir = e)
    df <- as_tibble(e$state_year_format3)

    # Pevehouse V3 state-unit format: one row per state-year, columns include
    # state identifiers + one column per IGO whose entries are status codes
    # (1 = full member, 0 = no member, 2/3 = associate/observer, -9 = NA).
    # Anything that isn't a known id column is treated as an IGO column.
    
    io_pct <- data.frame()
    data_years <- unique(df$year)
    for (i in 4:ncol(df))
    {
      temp_active <- df[,c(2, i)]
      active_years <- c(min(df$year[temp_active[,2] == 1]), max(df$year[temp_active[,2] == 1]))
      io_pct_temp <- data.frame(
            io = colnames(df)[i],
            year = seq(active_years[1], active_years[2]),
            pct = NA
        )
      for (z in 1:nrow(io_pct_temp))
      {
        temp_df <- subset(df, year == io_pct_temp$year[z])
        io_pct_temp$pct[z] <- sum(temp_df[,i] == 1)/nrow(temp_df)
      }
      io_pct <- io_pct |> 
        bind_rows(
          io_pct_temp
        )
    }
    
    return(list(io_pct, df))
  }
}

# ============================================================
# Load inputs
# ============================================================

cat("Loading inputs...\n")

dyad_data   <- readRDS(file.path(OUTPUT_DIR, "dyad_year_network_layers.rds"))
system_data <- readRDS(file.path(OUTPUT_DIR, "system_analysis_data.rds")) |>
  arrange(year) |>
  # [2026-10-12] conflict lags are built BEFORE dropping years without JS, so 1817 keeps its 1816 conflict
  (\(d) { stopifnot(all(diff(d$year) == 1)); d })() |>   # consecutive years
  mutate(
    log_js_divergence    = log1p(js_divergence),
    log_mid_count        = log1p(mid_count),
    log_fatal_mid_count  = log1p(fatal_mid_count),
    mid_count_lag1           = dplyr::lag(mid_count,           1),
    log_mid_count_lag1       = dplyr::lag(log_mid_count,       1),
    fatal_mid_count_lag1     = dplyr::lag(fatal_mid_count,     1),
    log_fatal_mid_count_lag1 = dplyr::lag(log_fatal_mid_count, 1),
    unit = 1L
  ) |>
  filter(!is.na(mid_count_lag1), !is.na(js_divergence))
state_churn_base <- readRDS(file.path(OUTPUT_DIR, "state_network_churn.rds"))
# [2026-10-12] Headline state IVs are onset counts at t-1, taken from the full COW state-year universe in
# dyad_data (same construction as 01_data_cleaning/02_compute_network_measures.R and 02). The churn panel
# holds transition rows only, so lagging the previous ROW would drop every state's first transition year.
# A t-1 value is NA only if the state was not in the system in t-1 (after a COW membership gap).
onset_univ <- dplyr::bind_rows(dyad_data |> dplyr::distinct(ccode = ccode1, year),
                               dyad_data |> dplyr::distinct(ccode = ccode2, year)) |>
  dplyr::distinct() |>
  dplyr::left_join(
    dyad_data |>
      dplyr::filter(cowmidongoing == 1, ccode1 < ccode2) |>
      dplyr::select(ccode1, ccode2, year, cowmidfatal, cowmidonset) |>
      tidyr::pivot_longer(c(ccode1, ccode2), values_to = "ccode") |>
      dplyr::group_by(ccode, year) |>
      dplyr::summarise(u_onset = sum(cowmidonset == 1, na.rm = TRUE),
                       u_fatal_onset = sum(cowmidfatal == 1 & cowmidonset == 1, na.rm = TRUE),
                       .groups = "drop"),
    by = c("ccode", "year")) |>
  dplyr::mutate(u_onset = as.numeric(tidyr::replace_na(u_onset, 0L)),
                u_fatal_onset = as.numeric(tidyr::replace_na(u_fatal_onset, 0L)))
chk_univ <- dplyr::inner_join(state_churn_base |> dplyr::select(ccode, year, mid_onset_count, fatal_mid_onset_count),
                              onset_univ, by = c("ccode", "year"))
stopifnot(nrow(chk_univ) == nrow(state_churn_base), all(chk_univ$mid_onset_count == chk_univ$u_onset),
          all(chk_univ$fatal_mid_onset_count == chk_univ$u_fatal_onset))
rm(chk_univ)
state_churn_base <- state_churn_base |>
  dplyr::select(-dplyr::any_of(c("onset_lag1", "fatal_onset_lag1", "log_onset_lag1", "log_fatal_onset_lag1"))) |>
  dplyr::left_join(onset_univ |>
                     dplyr::transmute(ccode, year = year + 1L,
                                      onset_lag1 = u_onset, fatal_onset_lag1 = u_fatal_onset,
                                      log_onset_lag1 = log1p(u_onset), log_fatal_onset_lag1 = log1p(u_fatal_onset)),
                   by = c("ccode", "year"))
stopifnot(!anyDuplicated(state_churn_base[c("ccode", "year")]))

igo_input <- load_igo_state_year()
io_pct    <- igo_input[[1]]
igo_wide  <- igo_input[[2]]

# Standard id columns in Pevehouse v3 state-year format; everything else
# is treated as an IGO column.
.IGO_ID_COLS <- intersect(names(igo_wide),
                           c("ccode","year","state","country","stateabb",
                             "ioname","version","longname"))
.IGO_COLS    <- setdiff(names(igo_wide), .IGO_ID_COLS)

cat(sprintf("  IGO panel: %d state-years, %d IGO columns, years %d-%d\n",
            nrow(igo_wide),
            length(.IGO_COLS),
            min(igo_wide$year, na.rm = TRUE),
            max(igo_wide$year, na.rm = TRUE)))
cat(sprintf("  IGO-year share rows: %d (median pct %.3f)\n",
            nrow(io_pct), median(io_pct$pct, na.rm = TRUE)))

# ============================================================
# Core rebuild: given a filtered IGO list for each year, compute
# the new dyad-year shared IGO count, build the composite main_weight
# with the (untouched) alliance layer, and recompute KL / churn.
# ============================================================

# Pre-compute dyad pairs once per year from the current dyad_data panel.
# Use ccode1 < ccode2 to match peacesciencer ordering.

# Save alliance-only contribution (layer_offensive + layer_defensive),
# which does NOT change across thresholds.
alliance_layer <- dyad_data |>
  select(ccode1, ccode2, year,
         layer_offensive, layer_defensive) |>
  mutate(alliance_raw = coalesce(layer_offensive, 0) +
                       coalesce(layer_defensive, 0),
         ccode1 = as.integer(ccode1),
         ccode2 = as.integer(ccode2),
         year   = as.integer(year)) |>
  # peacesciencer::create_dyadyears() defaults to directed = TRUE, so
  # dyad_data ships both (i, j) and (j, i) for every undirected dyad-year.
  # build_shared_counts_wide only emits (i, j) with i < j (upper triangle
  # of M %*% t(M)). Without filtering here, the (j, i) rows in
  # alliance_layer would left_join to NA layer_igo and become alliance-
  # only main_weight while the (i, j) row carries alliance + IGO --
  # asymmetric weights for the same undirected tie. Then compute_state_churn's
  # bind_rows + swap double-counts via the duplicated directed rows and
  # blows up sum(pmin)/sum(pmax). Symmetric dyadigos in the main pipeline
  # masks this; the rebuild can't.
  filter(ccode1 < ccode2) |>
  select(ccode1, ccode2, year, alliance_raw)

# Baseline: what fraction of the current layer_igo coverage falls in each year?
baseline_igo_edges <- dyad_data |>
  select(ccode1, ccode2, year, layer_igo_base = layer_igo)

# Build function (wide format / matrix product): for each year, build a
# state x IGO membership matrix from `igo_wide` restricted to IGOs whose
# per-year share lies in [L, U] (per `io_pct`), then `tcrossprod` gives
# the dyad shared-count matrix in one shot. Faster and uses io_pct
# directly as the share table.
build_shared_counts_wide <- function(igo_wide, io_pct, L = 0, U = 1) {
  id_cols  <- intersect(names(igo_wide),
                         c("ccode","year","state","country","stateabb",
                           "ioname","version","longname"))
  igo_cols <- setdiff(names(igo_wide), id_cols)
  years    <- sort(unique(igo_wide$year))

  purrr::map_dfr(years, function(yr) {
    df_yr <- igo_wide |> filter(year == yr)
    if (nrow(df_yr) < 2) return(NULL)

    # IGOs to KEEP this year: pct in [L, U]
    keep <- io_pct |>
      filter(year == yr, pct >= L, pct <= U) |>
      pull(io)
    keep <- intersect(keep, igo_cols)
    if (length(keep) == 0) return(NULL)

    # state x IGO binary membership (1 = full member, else 0)
    M <- as.matrix(df_yr[, keep, drop = FALSE])
    M[] <- suppressWarnings(as.integer(M))
    M[is.na(M)] <- 0L
    M <- (M == 1L) * 1L

    # Sort rows by ccode so upper.tri yields ccode1 < ccode2
    ord    <- order(as.integer(df_yr$ccode))
    M      <- M[ord, , drop = FALSE]
    ccodes <- as.integer(df_yr$ccode)[ord]

    # Off-diagonal shared full memberships per dyad
    S  <- tcrossprod(M)
    ix <- which(upper.tri(S) & S > 0, arr.ind = TRUE)
    if (nrow(ix) == 0) return(NULL)

    tibble(
      ccode1   = ccodes[ix[, 1]],
      ccode2   = ccodes[ix[, 2]],
      year     = as.integer(yr),
      n_shared = as.integer(S[ix])
    )
  })
}

# Combine shared_counts with alliance_layer -> main_weight panel.
# Matches 02_compute_network_measures.R:
#   layer_igo   = n_shared / max(n_shared) within year
#   main_raw    = alliance_raw + layer_igo
#   main_weight = main_raw / max(main_raw) within year
rebuild_panel_from_shared <- function(shared_counts, alliance_layer) {
  shared_counts <- shared_counts |>
    group_by(year) |>
    mutate(layer_igo = if (max(n_shared, na.rm = TRUE) > 0) {
      n_shared / max(n_shared, na.rm = TRUE)
    } else 0) |>
    ungroup() |>
    select(ccode1, ccode2, year, layer_igo)

  alliance_layer |>
    left_join(shared_counts, by = c("ccode1", "ccode2", "year")) |>
    mutate(layer_igo = replace_na(layer_igo, 0),
           main_raw  = alliance_raw + layer_igo) |>
    group_by(year) |>
    mutate(main_weight = if (max(main_raw, na.rm = TRUE) > 0) {
      main_raw / max(main_raw, na.rm = TRUE)
    } else 0) |>
    ungroup() |>
    select(-main_raw)
}

# System-level KL divergence and state-level churn from the rebuilt panel.
compute_system_kl <- function(panel) {
  years_sorted <- sort(unique(panel$year))
  out <- tibble(year = integer(), js_divergence = double())
  prev <- NULL
  for (yr in years_sorted) {
    d <- panel |> filter(year == yr, main_weight > 0)
    if (nrow(d) == 0) { prev <- NULL; next }
    dyad_key <- paste0(d$ccode1, "_", d$ccode2)
    w <- d$main_weight
    p_curr <- setNames(w / sum(w), dyad_key)
    if (!is.null(prev)) {
      all_k <- union(names(prev), names(p_curr))
      eps <- 1e-10
      p <- (p_curr[all_k] %||% 0); p[is.na(p)] <- 0; p <- p + eps; p <- p/sum(p)
      q <- (prev[all_k]  %||% 0);   q[is.na(q)] <- 0; q <- q + eps; q <- q/sum(q)
      m <- 0.5 * (p + q)
      kl <- 0.5 * sum(p * log(p / m)) + 0.5 * sum(q * log(q / m))  # JS divergence
      out <- bind_rows(out, tibble(year = yr, js_divergence = kl))
    }
    prev <- p_curr
  }
  out
}

compute_state_churn <- function(panel) {
  # [2026-10-12] Canonical four-case rule over the full membership universe, as in
  # 01_data_cleaning/02_compute_network_measures.R and 07: for every state in the system in both t-1 and t,
  # both isolated -> 0; exactly one side isolated (first partner gained / last partner lost) -> 1;
  # otherwise 1 - weighted Jaccard. The old gate kept only states with a positive tie in BOTH years, so
  # isolation transitions were dropped and the dropped set moved with the IGO threshold.
  # `panel` has one row per undirected dyad (ccode1 < ccode2) for every system member pair.
  universe <- bind_rows(panel |> distinct(ccode = ccode1, year),
                        panel |> distinct(ccode = ccode2, year)) |>
    distinct()
  trans <- universe |> semi_join(universe |> mutate(year = year + 1L), by = c("ccode", "year"))
  ties1 <- panel |> filter(main_weight > 0) |> transmute(i = ccode1, j = ccode2, year, w = main_weight)
  ties <- bind_rows(ties1, ties1 |> select(i = j, j = i, year, w)) |>
    group_by(i, j, year) |>
    summarise(w = max(w), .groups = "drop")
  pairs <- full_join(ties |> transmute(i, j, year = year + 1L, w1 = w),
                     ties |> transmute(i, j, year, w2 = w),
                     by = c("i", "j", "year")) |>
    mutate(w1 = replace_na(w1, 0), w2 = replace_na(w2, 0))
  agg <- pairs |>
    group_by(ccode = i, year) |>
    summarise(n1 = sum(w1 > 0), n2 = sum(w2 > 0),
              num = sum(pmin(w1, w2)), den = sum(pmax(w1, w2)), .groups = "drop")
  trans |>
    left_join(agg, by = c("ccode", "year")) |>
    mutate(n1 = replace_na(n1, 0L), n2 = replace_na(n2, 0L),
           num = replace_na(num, 0), den = replace_na(den, 0),
           churn = case_when(n1 == 0 & n2 == 0 ~ 0,
                             n1 == 0 | n2 == 0 ~ 1,
                             den > 0 ~ 1 - num / den,
                             TRUE ~ NA_real_)) |>
    arrange(ccode, year) |>
    select(ccode, year, churn)
}

# ============================================================
# Sweep across thresholds
# ============================================================
#
# For each (L, U) combination:
#   - per year, keep IGOs with share in (L, U)
#   - rebuild panel, KL, churn
#   - refit system OLS (HAC NW lag 5) and state TWFE
#   - store fatal / non-fatal coefs

# IVs assessed separately (not jointly). Each threshold spec refits 4
# system-level and 4 state-level models.
sys_ivs <- list(
  list(col = "mid_count_lag1",           label = "MID count (t-1)"),
  list(col = "log_mid_count_lag1",       label = "log(1 + MID count) (t-1)"),
  list(col = "fatal_mid_count_lag1",     label = "Fatal MID count (t-1)"),
  list(col = "log_fatal_mid_count_lag1", label = "log(1 + Fatal MID count) (t-1)")
)
# [2026-10-12] state onset IVs at t-1 (headline timing)
state_ivs <- list(
  list(col = "onset_lag1",           label = "MID onset (t-1)"),
  list(col = "log_onset_lag1",       label = "log(1 + MID onset) (t-1)"),
  list(col = "fatal_onset_lag1",     label = "Fatal MID onset (t-1)"),
  list(col = "log_fatal_onset_lag1", label = "log(1 + Fatal MID onset) (t-1)")
)

# Re-use state_churn controls if already present; otherwise build them
# inline (same pattern as 11 / 14 / 15 / 16 / 18).
state_controls_needed <- c("partner_war_frac", "democracy", "autocracy",
                            "cinc", "major_power")
if (!all(state_controls_needed %in% names(state_churn_base))) {
  cat("Preparing state-level controls for IGO robustness pipeline...\n")

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
  state_churn_base <- state_churn_base |>
    left_join(state_attrs, by = c("ccode", "year"))

  state_years <- create_stateyears(
    subset_years = min(state_churn_base$year):max(state_churn_base$year))
  attr(state_years, "ps_data_type") <- "state_year"
  attr(state_years, "ps_system")    <- "cow"
  state_years <- state_years |> add_democracy()
  polity <- state_years |>
    select(ccode, year, polity2) |>
    mutate(democracy = as.integer(polity2 >= 6),
           autocracy = as.integer(polity2 <= -6))
  state_churn_base <- state_churn_base |>
    left_join(polity, by = c("ccode", "year"))

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
  at_war_fatal <- dyad_data |>
    filter(cowmidfatal == 1) |>
    select(ccode1, ccode2, year) |>
    pivot_longer(c(ccode1, ccode2), values_to = "ccode") |>
    distinct(ccode, year) |>
    mutate(at_war_fatal = 1L)
  pwe <- partner_sets |>
    mutate(year_next = year + 1L) |>
    left_join(at_war,       by = c("ccode2" = "ccode", "year_next" = "year")) |>
    left_join(at_war_fatal, by = c("ccode2" = "ccode", "year_next" = "year")) |>
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
  state_churn_base <- state_churn_base |>
    left_join(pwe, by = c("ccode", "year")) |>
    mutate(partner_war_frac                = replace_na(partner_war_frac, 0),
           partner_war_frac_weighted       = replace_na(partner_war_frac_weighted, 0),
           partner_fatal_war_frac          = replace_na(partner_fatal_war_frac, 0),
           partner_fatal_war_frac_weighted = replace_na(partner_fatal_war_frac_weighted, 0),
           major_power = replace_na(major_power, 0))
}

# --- State-level control set ---
# [2026-10-12] partner war exposure removed from the control set (H3 dropped from the paper).
# CTRL_MID and CTRL_FATAL are now identical (democracy + autocracy + cinc +
# major_power + is_isolated_t1); both names and ctrl_for_iv() are kept so the
# call sites are unchanged.
CTRL_MID   <- "democracy + autocracy + cinc + major_power + is_isolated_t1"
CTRL_FATAL <- "democracy + autocracy + cinc + major_power + is_isolated_t1"
ctrl_for_iv <- function(iv) if (grepl("fatal", iv)) CTRL_FATAL else CTRL_MID

sys_ctrl_str   <- "log1p(n_democracies) + log1p(n_autocracies) + log1p(n_states) + poly(year, 2)"
# state_ctrl_str is now routed per IV via ctrl_for_iv() inside run_spec().

# Thresholds are FRACTIONS now (matching io_pct$pct in [0, 1]).
#   Lower: drop IGOs with share <  L
#   Upper: drop IGOs with share >  U
lower_grid <- seq(0.01, 0.15, length.out = 4)   # 1% .. 15%
upper_grid <- seq(0.85, 0.99, length.out = 4)   # 85% .. 99%

run_spec <- function(L, U) {
  cat(sprintf("  L=%.2f%%, U=%.2f%%\n", L * 100, U * 100))

  # Build shared-IGO counts directly from wide df + io_pct
  shared <- build_shared_counts_wide(igo_wide, io_pct, L = L, U = U)
  panel  <- rebuild_panel_from_shared(shared, alliance_layer)
  sys_kl <- compute_system_kl(panel)
  st_ch  <- compute_state_churn(panel)

  # Merge new KL onto system regression data; log-transform for headline DV.
  sys_df <- system_data |>
    select(-js_divergence, -log_js_divergence) |>
    left_join(sys_kl, by = "year") |>
    filter(!is.na(js_divergence)) |>
    mutate(log_js_divergence = log1p(js_divergence))

  # Merge new churn onto state regression data; add log_churn + log IVs.
  st_df <- state_churn_base |>
    left_join(st_ch |> rename(churn_recomp = churn),
              by = c("ccode", "year")) |>
    mutate(log_churn_recomp          = log1p(churn_recomp),
           log_mid_onset_count       = log1p(mid_onset_count),
           log_fatal_mid_onset_count = log1p(fatal_mid_onset_count)) |>
    filter(!is.na(churn_recomp))

  # Fit each IV separately at each level
  sys_rows <- purrr::map_dfr(sys_ivs, function(iv) {
    iv_col <- iv$col; iv_lab <- iv$label
    f <- as.formula(paste0("js_divergence ~ ", iv_col, " + ",
                            sys_ctrl_str))
    m <- tryCatch(lm(f, data = sys_df), error = function(e) NULL)
    if (is.null(m)) return(NULL)
    cf <- coef(m); vv <- vcov(m)
    tibble(level = "system", iv = iv_col, iv_label = iv_lab,
           estimate = unname(cf[iv_col]),
           se       = sqrt(diag(vv)[iv_col]),
           n_obs    = nobs(m))
  })

  st_rows <- purrr::map_dfr(state_ivs, function(iv) {
    iv_col <- iv$col; iv_lab <- iv$label
    f <- as.formula(paste0("churn_recomp ~ ", iv_col, " + ",
                            ctrl_for_iv(iv_col), " | ccode + year"))
    m <- tryCatch(feols(f, data = st_df, cluster = ~ ccode + year),
                   error = function(e) NULL)
    if (is.null(m)) return(NULL)
    cf <- coef(m); vv <- vcov(m)
    tibble(level = "state", iv = iv_col, iv_label = iv_lab,
           estimate = unname(cf[iv_col]),
           se       = sqrt(diag(vv)[iv_col]),
           n_obs    = nobs(m))
  })

  bind_rows(sys_rows, st_rows) |>
    mutate(L = L, U = U,
           ci_low  = estimate - 1.96 * se,
           ci_high = estimate + 1.96 * se)
}

# ------------------------------------------------------------
# Grid 1: lower-tail sweep (U = 100 fixed; drop IGOs below L%)
# ------------------------------------------------------------
cat("\n--- Grid 1: Lower-tail (drop niche IGOs) ---\n")
grid_lower <- purrr::map_dfr(lower_grid, ~ run_spec(L = .x, U = 1)) |>
  mutate(grid = "lower")

# ------------------------------------------------------------
# Grid 2: upper-tail sweep (L = 0 fixed; drop IGOs above U)
# ------------------------------------------------------------
cat("\n--- Grid 2: Upper-tail (drop near-universal IGOs) ---\n")
grid_upper <- purrr::map_dfr(upper_grid, ~ run_spec(L = 0, U = .x)) |>
  mutate(grid = "upper")

# ------------------------------------------------------------
# Grid 3: both tails sweep (full L x U Cartesian product)
# ------------------------------------------------------------
# Drop IGOs with share < L% OR share > U% simultaneously, across the
# full 4 x 4 grid of (L, U) pairs.
cat("\n--- Grid 3: Both tails (Cartesian L \u00d7 U) ---\n")
grid_both <- expand_grid(L = lower_grid, U = upper_grid) |>
  pmap_dfr(function(L, U) run_spec(L = L, U = U)) |>
  mutate(grid = "both")

iv_levels_all <- c(sapply(sys_ivs, `[[`, "label"),
                    sapply(state_ivs, `[[`, "label"))

all_results <- bind_rows(grid_lower, grid_upper, grid_both) |>
  mutate(iv_label = factor(iv_label, levels = iv_levels_all),
         level_label = factor(level, levels = c("system", "state"),
                               labels = c("System (JS)", "State (churn)")))

write_csv(all_results,
          file.path(OUTPUT_DIR, "robustness_igo_share_results.csv"))
cat(sprintf("\nSaved: robustness_igo_share_results.csv (%d rows)\n",
            nrow(all_results)))

# ============================================================
# Plots
# ============================================================

cat("\n--- Plotting ---\n")

# Round threshold labels for readability. L and U are stored as fractions
# (e.g. 0.10 -> "10.00%"), so multiply by 100 for the axis text.
fmt_pct <- function(x) sprintf("%.2f%%", x * 100)

# --- Lower tail ---
p_lower_system <- all_results |>
  filter(grid == "lower", 
         level == "system") |>
  mutate(L_label = factor(fmt_pct(L), levels = fmt_pct(sort(unique(L))))) |>
  ggplot(aes(x = L_label, y = estimate)) +
  geom_hline(yintercept = 0, color = "gray55", linewidth = 0.3) +
  geom_errorbar(aes(ymin = ci_low, ymax = ci_high),
                width = 0.15, linewidth = 0.8, color = ACCENT) +
  geom_point(size = 2.5, color = ACCENT) +
  facet_wrap( ~ iv_label, scales = "free_y",
              labeller = label_wrap_gen(width = 22)) +
  labs(
    title = "System: lower-tail IGO drop",
    x = "Threshold L", y = "Coefficient (95% CI)"
  ) +
  theme_pres() +
  theme(panel.spacing = unit(1.2, "lines"),
        panel.spacing.x = unit(1.5, "lines"),
        strip.placement = "outside",
        strip.text = element_text(size = 9, face = "bold",
                                    margin = margin(4, 4, 4, 4)),
        plot.margin = margin(t = 10, r = 12, b = 14, l = 10))

save_pres("fig_ms_igo_share_lower_system.png", p_lower_system,
          width = 9, height = 5.5)


p_lower_state <- all_results |>
  filter(grid == "lower", 
         level != "system") |>
  mutate(L_label = factor(fmt_pct(L), levels = fmt_pct(sort(unique(L))))) |>
  ggplot(aes(x = L_label, y = estimate)) +
  geom_hline(yintercept = 0, color = "gray55", linewidth = 0.3) +
  geom_errorbar(aes(ymin = ci_low, ymax = ci_high),
                width = 0.15, linewidth = 0.8, color = ACCENT) +
  geom_point(size = 2.5, color = ACCENT) +
  facet_wrap( ~ iv_label, scales = "free_y",
              labeller = label_wrap_gen(width = 22)) +
  labs(
    title = "State: lower-tail IGO drop",
    x = "Threshold L", y = "Coefficient (95% CI)"
  ) +
  theme_pres() +
  theme(panel.spacing = unit(1.2, "lines"),
        panel.spacing.x = unit(1.5, "lines"),
        strip.placement = "outside",
        strip.text = element_text(size = 9, face = "bold",
                                  margin = margin(4, 4, 4, 4)),
        plot.margin = margin(t = 10, r = 12, b = 14, l = 10))

save_pres("fig_ms_igo_share_lower_state.png", p_lower_state,
          width = 9, height = 5.5)

# --- Upper tail: SYSTEM ---
p_upper_system <- all_results |>
  filter(grid == "upper",
         level == "system") |>
  mutate(U_label = factor(fmt_pct(U), levels = fmt_pct(sort(unique(U))))) |>
  ggplot(aes(x = U_label, y = estimate)) +
  geom_hline(yintercept = 0, color = "gray55", linewidth = 0.3) +
  geom_errorbar(aes(ymin = ci_low, ymax = ci_high),
                width = 0.15, linewidth = 0.8, color = ACCENT) +
  geom_point(size = 2.5, color = ACCENT) +
  facet_wrap(~ iv_label, scales = "free_y",
             labeller = label_wrap_gen(width = 22)) +
  labs(title = "System: upper-tail IGO drop",
       x = "Threshold U", y = "Coefficient (95% CI)") +
  theme_pres() +
  theme(panel.spacing = unit(1.2, "lines"),
        panel.spacing.x = unit(1.5, "lines"),
        strip.placement = "outside",
        strip.text = element_text(size = 9, face = "bold",
                                  margin = margin(4, 4, 4, 4)),
        plot.margin = margin(t = 10, r = 12, b = 14, l = 10))

save_pres("fig_ms_igo_share_upper_system.png", p_upper_system,
          width = 9, height = 5.5)

# --- Upper tail: STATE ---
p_upper_state <- all_results |>
  filter(grid == "upper",
         level != "system") |>
  mutate(U_label = factor(fmt_pct(U), levels = fmt_pct(sort(unique(U))))) |>
  ggplot(aes(x = U_label, y = estimate)) +
  geom_hline(yintercept = 0, color = "gray55", linewidth = 0.3) +
  geom_errorbar(aes(ymin = ci_low, ymax = ci_high),
                width = 0.15, linewidth = 0.8, color = ACCENT) +
  geom_point(size = 2.5, color = ACCENT) +
  facet_wrap(~ iv_label, scales = "free_y",
             labeller = label_wrap_gen(width = 22)) +
  labs(title = "State: upper-tail IGO drop",
       x = "Threshold U", y = "Coefficient (95% CI)") +
  theme_pres() +
  theme(panel.spacing = unit(1.2, "lines"),
        panel.spacing.x = unit(1.5, "lines"),
        strip.placement = "outside",
        strip.text = element_text(size = 9, face = "bold",
                                  margin = margin(4, 4, 4, 4)),
        plot.margin = margin(t = 10, r = 12, b = 14, l = 10))

save_pres("fig_ms_igo_share_upper_state.png", p_upper_state,
          width = 9, height = 5.5)

# --- Both tails (L x U): SYSTEM ---
p_both_system <- all_results |>
  filter(grid == "both",
         level == "system") |>
  mutate(L_label = factor(fmt_pct(L), levels = fmt_pct(sort(unique(L)))),
         U_label = factor(fmt_pct(U), levels = fmt_pct(sort(unique(U))))) |>
  ggplot(aes(x = L_label, y = estimate,
             color = U_label, group = U_label)) +
  geom_hline(yintercept = 0, color = "gray55", linewidth = 0.3) +
  geom_errorbar(aes(ymin = ci_low, ymax = ci_high),
                width = 0.15, linewidth = 0.6,
                position = position_dodge(width = 0.55)) +
  geom_point(size = 2.2, position = position_dodge(width = 0.55)) +
  scale_color_brewer(palette = "Dark2", name = "Upper U") +
  facet_wrap(~ iv_label, scales = "free_y",
             labeller = label_wrap_gen(width = 22)) +
  labs(title = "System: both-tail IGO drop",
       x = "Lower threshold L", y = "Coefficient (95% CI)") +
  theme_pres() +
  theme(panel.spacing = unit(1.2, "lines"),
        panel.spacing.x = unit(1.5, "lines"),
        strip.placement = "outside",
        strip.text = element_text(size = 9, face = "bold",
                                  margin = margin(4, 4, 4, 4)),
        plot.margin = margin(t = 10, r = 12, b = 14, l = 10))

save_pres("fig_ms_igo_share_both_system.png", p_both_system,
          width = 10, height = 5.5)

# --- Both tails (L x U): STATE ---
p_both_state <- all_results |>
  filter(grid == "both",
         level != "system") |>
  mutate(L_label = factor(fmt_pct(L), levels = fmt_pct(sort(unique(L)))),
         U_label = factor(fmt_pct(U), levels = fmt_pct(sort(unique(U))))) |>
  ggplot(aes(x = L_label, y = estimate,
             color = U_label, group = U_label)) +
  geom_hline(yintercept = 0, color = "gray55", linewidth = 0.3) +
  geom_errorbar(aes(ymin = ci_low, ymax = ci_high),
                width = 0.15, linewidth = 0.6,
                position = position_dodge(width = 0.55)) +
  geom_point(size = 2.2, position = position_dodge(width = 0.55)) +
  scale_color_brewer(palette = "Dark2", name = "Upper U") +
  facet_wrap(~ iv_label, scales = "free_y",
             labeller = label_wrap_gen(width = 22)) +
  labs(title = "State: both-tail IGO drop",
       x = "Lower threshold L", y = "Coefficient (95% CI)") +
  theme_pres() +
  theme(panel.spacing = unit(1.2, "lines"),
        panel.spacing.x = unit(1.5, "lines"),
        strip.placement = "outside",
        strip.text = element_text(size = 9, face = "bold",
                                  margin = margin(4, 4, 4, 4)),
        plot.margin = margin(t = 10, r = 12, b = 14, l = 10))

save_pres("fig_ms_igo_share_both_state.png", p_both_state,
          width = 10, height = 5.5)

cat("\nDone.\n")
