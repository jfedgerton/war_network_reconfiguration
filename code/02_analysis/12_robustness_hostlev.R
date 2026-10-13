###############################################################################
# 13_robustness_hostlev.R
# Robustness: replace fatal/non-fatal binary with COW MID hostility level
# (HostLev 1-5). Fatal/non-fatal is a coarse proxy for "informativeness"
# of a conflict event; HostLev is a sharper ordinal measure of how
# coercive the dispute was.
#
# HostLev coding (COW MID v5):
#   1 = no militarized action
#   2 = threat to use force
#   3 = display of force
#   4 = use of force
#   5 = war
#
# At each level we try three specifications:
#   (a) max_hostlev as a linear continuous IV (1-5)
#   (b) max_hostlev as ordered factor (dummies for levels 3, 4, 5
#       with <=2 as reference) — shows dose-response shape
#   (c) sum_hostlev: total severity exposure (sum across MIDs)
#
# State-level: for each (ccode, year), find all MIDs the state was in.
#   max_hostlev_{i,t} = max HostLev across its MIDs that year
#   sum_hostlev_{i,t} = sum across its MIDs that year
#
# System-level: over all MIDs active system-wide that year.
#
# Data source: peacesciencer's directed dyad-year MID panel (if available)
# OR CoW MIDB v5 CSV. Loader below tries both.
###############################################################################

library(tidyverse)
library(fixest)
library(peacesciencer)
select <- dplyr::select

set.seed(123)

OUTPUT_DIR <- here::here("data", "processed")
RAW_DIR    <- here::here("data", "raw", "mid")
FIG_DIR    <- here::here("code", "03_figures", "output", "manuscript")
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(RAW_DIR, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------
# Manuscript theme — white background, publication-suitable
# ------------------------------------------------------------
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

# ------------------------------------------------------------
# Load MID hostility data
# ------------------------------------------------------------
# peacesciencer bundles `cow_mid_dirdisps` (directed dispute-level) with
# HostLev; fallback is the raw COW MIDB v5 CSV at data/raw/mid/MIDB_5.0.csv.

load_mid_hostlev <- function() {
  if (requireNamespace("peacesciencer", quietly = TRUE)) {

    # --- (1) Try directed-dyadic dispute-level tables (have ccode1 +
    #         hostlev). Different peacesciencer versions use different
    #         names: cow_mid_ddydisps, cow_mid_dirdisps, etc.
    candidates <- c("cow_mid_ddydisps", "cow_mid_dirdisps",
                    "cow_mid_dirdisp",  "cow_mid_dyads",
                    "cow_mid_ddyads")
    for (nm in candidates) {
      e <- new.env()
      ok <- tryCatch({
        utils::data(list = nm, package = "peacesciencer", envir = e)
        TRUE
      }, warning = function(w) FALSE, error = function(e) FALSE)
      if (!(ok && exists(nm, envir = e))) next

      df <- as_tibble(get(nm, envir = e))
      if (!"ccode1" %in% names(df)) next

      # Find a hostility column (case-insensitive)
      hl_pat <- grepl("(?i)host[._-]?lev", names(df), perl = TRUE)
      hl_cols <- names(df)[hl_pat]
      if (length(hl_cols) == 0) next

      # Prefer hihostlev > hostlev > anything matching
      hl_priority <- c("hihostlev", "hostlev", "HostLev", "host_lev")
      hl_col <- intersect(hl_priority, hl_cols)
      hl_col <- if (length(hl_col) > 0) hl_col[1] else hl_cols[1]

      cat(sprintf("Using peacesciencer::%s (hostlev col: %s)\n", nm, hl_col))
      out <- df |>
        rename(hostlev = !!sym(hl_col))

      # Find fatality column (case-insensitive regex on "fatal")
      fat_pat <- grepl("(?i)fatal", names(out), perl = TRUE)
      fat_cols <- names(out)[fat_pat]
      fat_priority <- c("fatality", "fatal", "cowmidfatal", "fatalmid")
      fat_col <- intersect(fat_priority, fat_cols)
      fat_col <- if (length(fat_col) > 0) fat_col[1]
                 else if (length(fat_cols) > 0) fat_cols[1]
                 else NA_character_
      if (!is.na(fat_col)) {
        cat(sprintf("    Fatality col: %s\n", fat_col))
        out$fatal <- as.integer(suppressWarnings(
          as.numeric(out[[fat_col]])) > 0)
      } else {
        cat("    No fatality column found; proxying fatal = HostLev >= 4 (use of force or war).\n")
        out$fatal <- as.integer(out$hostlev >= 4)
      }

      # Look for a dispute-id column so we can identify onsets cleanly
      # (year == min(year) per (state, dispnum)). If absent, the loader
      # treats each row as an onset (most peacesciencer dyad-dispute
      # panels are already keyed at the dispute-start year).
      disp_pat <- grepl("(?i)disp(num|id)?", names(out), perl = TRUE)
      disp_cols <- names(out)[disp_pat]
      disp_priority <- c("dispnum3", "dispnum", "dispid")
      disp_col <- intersect(disp_priority, disp_cols)
      disp_col <- if (length(disp_col) > 0) disp_col[1]
                  else if (length(disp_cols) > 0) disp_cols[1]
                  else NA_character_
      if (!is.na(disp_col)) {
        out$dispnum <- as.integer(out[[disp_col]])
      } else {
        out$dispnum <- NA_integer_
      }

      keep <- intersect(c("ccode1", "ccode2", "year",
                           "hostlev", "fatal", "dispnum"),
                         names(out))
      out <- out |> select(all_of(keep))
      if (!"ccode2"  %in% names(out)) out$ccode2  <- NA_integer_
      if (!"dispnum" %in% names(out)) out$dispnum <- NA_integer_
      return(out)
    }

    # --- (2) State-year aggregation via add_cow_mids(). This gives one
    #         row per (ccode, year) with the maximum hostility level of
    #         any MID that state was in. We collapse it into the same
    #         (ccode1, ccode2, year, hostlev) shape the downstream
    #         pipeline expects (with ccode2 = NA).
    cat("Trying peacesciencer::add_cow_mids() state-year aggregation...\n")
    sy <- create_stateyears(
      subset_years = min(1816, min(state_churn$year)):max(2014, max(state_churn$year))
    )
    attr(sy, "ps_data_type") <- "state_year"
    attr(sy, "ps_system")    <- "cow"
    sy_mids <- tryCatch(sy |> add_cow_mids(), error = function(e) NULL)
    if (!is.null(sy_mids)) {
      hl_pat <- grepl("(?i)host[._-]?lev", names(sy_mids), perl = TRUE)
      hl_cols <- names(sy_mids)[hl_pat]
      # [2026-10-12] The IV is the hostility of the state's MID ONSETS, so only an onset-level hostility column
      # is acceptable here. Ongoing-max hostility would give a continuing war without a new onset a positive
      # value. If no onset column exists, fall through to the MIDB CSV fallback below.
      hl_cols <- intersect(c("cowmidonsethihostlev", "midonsethihostlev"), hl_cols)
      if (length(hl_cols) > 0) {
        hl_col <- hl_cols[1]
        cat(sprintf("  Using state-year col: %s\n", hl_col))
        # In state-year aggregation, fatality is approximated by HostLev>=4
        # and dispute-id is unavailable; we treat each state-year row as
        # an onset (loses true onset granularity, but is a usable proxy).
        return(sy_mids |>
          transmute(ccode1 = ccode,
                    ccode2 = NA_integer_,
                    year   = year,
                    hostlev = .data[[hl_col]],
                    fatal  = as.integer(.data[[hl_col]] >= 4),
                    dispnum = NA_integer_) |>
          filter(!is.na(hostlev), hostlev >= 1, hostlev <= 5))
      }
    }
  }

  # --- (3) Raw MIDB v5 CSV fallback ---
  csv_path <- file.path(RAW_DIR, "MIDB_5.0.csv")
  if (!file.exists(csv_path)) {
    stop("MID hostility data not found. Tried peacesciencer dyadic tables ",
         "and add_cow_mids() but neither produced a usable hostlev column. ",
         "Place COW MIDB v5 CSV at ", csv_path,
         " as a final fallback.")
  }
  cat(sprintf("Reading MIDB from %s ...\n", csv_path))
  raw <- readr::read_csv(csv_path, show_col_types = FALSE) |>
    rename_with(tolower)
  need <- c("ccode", "styear", "endyear", "hostlev")
  if (!all(need %in% names(raw))) {
    stop("MIDB CSV missing expected columns: ",
         paste(setdiff(need, names(raw)), collapse = ", "))
  }
  fatal_col <- intersect(c("fatality", "fatal"), names(raw))
  fatal_col <- if (length(fatal_col) > 0) fatal_col[1] else NA_character_
  disp_col_csv <- intersect(c("dispnum3", "dispnum", "dispid"), names(raw))
  raw <- raw |>
    mutate(styear = pmax(styear, 1816L),
           endyear = pmin(endyear, 2014L),
           fatal = if (!is.na(fatal_col))
             as.integer(suppressWarnings(as.numeric(.data[[fatal_col]])) > 0)
           else as.integer(hostlev >= 4),
           dispnum = if (length(disp_col_csv) > 0)
             as.integer(.data[[disp_col_csv[1]]])
           else NA_integer_) |>
    filter(!is.na(styear), !is.na(endyear), endyear >= styear) |>
    rowwise() |>
    mutate(years = list(styear:endyear)) |>
    ungroup() |>
    tidyr::unnest(years)
  raw |>
    transmute(ccode1 = ccode, ccode2 = NA_integer_,
              year = years, hostlev, fatal, dispnum)
}

# ------------------------------------------------------------
# Load other inputs
# ------------------------------------------------------------

cat("Loading inputs...\n")

# State-level only — system-level KL hostility analysis intentionally
# omitted (max/mean HostLev aggregates are uninformative at the system
# level; see note further down in the aggregation block).
state_churn <- readRDS(file.path(OUTPUT_DIR, "state_network_churn.rds"))
dyad_data   <- readRDS(file.path(OUTPUT_DIR, "dyad_year_network_layers.rds"))

# Build state-level controls inline if not already attached.
needs_prep <- !all(c("cinc", "major_power", "democracy", "autocracy",
                      "partner_war_frac") %in% names(state_churn))
if (needs_prep) {
  cat("Preparing state-level controls inline (cinc / democracy / partner_war_frac / ...)...\n")

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
  state_churn <- state_churn |> left_join(state_attrs, by = c("ccode", "year"))

  state_years <- create_stateyears(
    subset_years = min(state_churn$year):max(state_churn$year))
  attr(state_years, "ps_data_type") <- "state_year"
  attr(state_years, "ps_system")    <- "cow"
  state_years <- state_years |> add_democracy()
  polity <- state_years |>
    select(ccode, year, polity2) |>
    mutate(democracy = as.integer(polity2 >= 6),
           autocracy = as.integer(polity2 <= -6))
  state_churn <- state_churn |> left_join(polity, by = c("ccode", "year"))

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
  state_churn <- state_churn |>
    left_join(pwe, by = c("ccode", "year")) |>
    mutate(partner_war_frac                = replace_na(partner_war_frac, 0),
           partner_war_frac_weighted       = replace_na(partner_war_frac_weighted, 0),
           partner_fatal_war_frac          = replace_na(partner_fatal_war_frac, 0),
           partner_fatal_war_frac_weighted = replace_na(partner_fatal_war_frac_weighted, 0),
           major_power = replace_na(major_power, 0))
}

# --- State-level control set ---
# Non-fatal and fatal IVs share the same state-level controls:
# democracy + autocracy + cinc + major_power + is_isolated_t1.
# HostLev IVs (mean_hostlev_onset / log_mean_hostlev_onset) are not
# named "fatal", so they route to CTRL_MID (any-MID baseline).
# [2026-10-12] partner war exposure removed from the control set (H3 dropped from the paper).
CTRL_MID   <- "democracy + autocracy + cinc + major_power + is_isolated_t1"
CTRL_FATAL <- "democracy + autocracy + cinc + major_power + is_isolated_t1"
ctrl_for_iv <- function(iv) if (grepl("fatal", iv)) CTRL_FATAL else CTRL_MID

mid_hl <- load_mid_hostlev()
cat(sprintf("  MID panel: %d rows\n", nrow(mid_hl)))

# ------------------------------------------------------------
# Aggregate hostility to state-year and system-year
# ------------------------------------------------------------

# State-year aggregation. Mean HostLev across the state's MID ONSETS:
#   mean_hostlev_onset = average HostLev across the state's MID onsets
#                        in that year. An "onset" is the first year a
#                        state appears in a given dispute (year ==
#                        min(year) per (ccode, dispnum)).
#
# Mean isolates *intensity* from quantity: a state with one war (HL=5)
# and a state with five wars both score 5; a state with one threat
# (HL=2) and a state with five threats both score 2. State-years with
# no onsets are coded 0 (no MID at all = no violence).
# This is parallel to the main-text contemporaneous onset IV but tests
# "does the average violence intensity of a state's conflicts predict
# churn", not "does the total count predict churn".

# Identify onsets. If dispnum is missing for all rows, fall back to
# treating every row as an onset (peacesciencer's dyad-dispute panels
# typically already key at the dispute-start year).
if (all(is.na(mid_hl$dispnum))) {
  cat("  No dispute identifier in loader output; treating every MID-state-year row as an onset.\n")
  mid_hl_onset <- mid_hl
} else {
  cat("  Identifying onsets via min(year) per (ccode, dispnum).\n")
  mid_hl_onset <- mid_hl |>
    filter(!is.na(hostlev), hostlev >= 1, hostlev <= 5,
           !is.na(ccode1), !is.na(year)) |>
    group_by(ccode1, dispnum) |>
    filter(year == min(year, na.rm = TRUE)) |>
    ungroup()
}

state_hl <- mid_hl_onset |>
  filter(!is.na(hostlev), hostlev >= 1, hostlev <= 5) |>
  group_by(ccode = ccode1, year) |>
  summarise(
    mean_hostlev_onset = mean(hostlev, na.rm = TRUE),
    n_onsets           = n(),
    .groups = "drop"
  )

# NOTE: system-level (KL divergence) hostility analysis is intentionally
# omitted. Aggregating HostLev across all MIDs in a year (max or mean)
# is degenerate at the system level: max is almost always 5 in any year
# with an ongoing war, and mean is dominated by the relative frequency
# of low-hostility events rather than the substantive coercive load on
# the system. Hostility-level scaling is meaningful at the STATE level,
# where each state has a small number of MIDs whose hostility we can
# meaningfully summarise.

state_churn <- state_churn |>
  left_join(state_hl, by = c("ccode", "year")) |>
  mutate(mean_hostlev_onset     = replace_na(mean_hostlev_onset, 0),
         n_onsets               = replace_na(n_onsets, 0L),
         log_mean_hostlev_onset = log1p(mean_hostlev_onset),
         log_churn              = log1p(churn))

# [2026-10-12] Hostility of the state's MID onsets in t-1 (headline timing), from the full COW state-year
# universe (as in 02): 0 if the state was in the system in t-1 without an onset, NA if it was not in the system.
hl_univ <- dplyr::bind_rows(dyad_data |> dplyr::distinct(ccode = ccode1, year),
                            dyad_data |> dplyr::distinct(ccode = ccode2, year)) |>
  dplyr::distinct() |>
  dplyr::left_join(state_hl |> dplyr::select(ccode, year, mean_hostlev_onset), by = c("ccode", "year")) |>
  dplyr::mutate(mean_hostlev_onset = as.numeric(tidyr::replace_na(mean_hostlev_onset, 0)))
state_churn <- state_churn |>
  dplyr::left_join(hl_univ |> dplyr::transmute(ccode, year = year + 1L,
                                                mean_hostlev_onset_lag1 = mean_hostlev_onset,
                                                log_mean_hostlev_onset_lag1 = log1p(mean_hostlev_onset)),
                   by = c("ccode", "year"))
stopifnot(!anyDuplicated(state_churn[c("ccode", "year")]))

# ------------------------------------------------------------
# Helpers (reused pattern)
# ------------------------------------------------------------

coef_se      <- function(m) sqrt(diag(vcov(m)))
coef_pvalues <- function(m) 2 * (1 - pnorm(abs(coef(m) / coef_se(m))))

extract_row <- function(m, terms, label, level) {
  if (is.null(m)) return(NULL)
  cf <- coef(m); vv <- vcov(m)
  present <- intersect(terms, names(cf))
  if (length(present) == 0) return(NULL)
  tibble(
    spec  = label,
    level = level,
    term  = present,
    estimate = unname(cf[present]),
    se       = sqrt(diag(vv)[present]),
    n_obs    = nobs(m)
  )
}

# ============================================================
# IV LIST — 2 state-level IVs (mean HostLev + log version)
# ============================================================
# Hostility-level intensity measure: average violence level across the
# state's MID onsets that year. Bounded [0, 5] (0 = no onset that year).
# The "fatal / non-fatal" axis is not used because HostLev already
# encodes violence intensity (fatal events occupy the high end by
# construction). Two IVs is the honest cardinality for hostility
# scaling at the state level.

state_ivs <- list(
  list(col = "mean_hostlev_onset_lag1",
       label = "Mean HostLev (onsets, t-1)"),
  list(col = "log_mean_hostlev_onset_lag1",
       label = "log(1 + Mean HostLev) (onsets, t-1)")
)

# ctrl_state is now routed per IV via ctrl_for_iv() inside the iteration.

cat("\n=== State level ===\n")
state_rows <- purrr::map_dfr(state_ivs, function(iv) {
  iv_col <- iv$col; iv_lab <- iv$label
  f <- as.formula(paste0("churn ~ ", iv_col, " + ", ctrl_for_iv(iv_col),
                          " | ccode + year"))
  m <- feols(f, data = state_churn, cluster = ~ ccode + year)
  cf <- coef(m); vv <- vcov(m)
  tibble(level = "state", iv = iv_col, iv_label = iv_lab,
         estimate = unname(cf[iv_col]),
         se       = sqrt(diag(vv)[iv_col]),
         n_obs    = nobs(m))
})

results <- state_rows |>
  mutate(ci_low  = estimate - 1.96 * se,
         ci_high = estimate + 1.96 * se)

print(results, digits = 4, n = Inf)

write_csv(results, file.path(OUTPUT_DIR, "robustness_hostlev_results.csv"))
cat(sprintf("\nSaved: robustness_hostlev_results.csv (%d rows)\n",
            nrow(results)))

# ============================================================
# PLOT — coefficient grid (8 IVs per level)
# ============================================================

cat("\n--- Plotting ---\n")

plot_df <- results |>
  mutate(iv_label = factor(iv_label,
                            levels = c("Mean HostLev (onsets, t-1)",
                                        "log(1 + Mean HostLev) (onsets, t-1)"),
                           labels = c("Mean Hostility\nLevel (onsets, t-1)",
                                      "log(1 + Mean\nHostility\nLevel) (onsets, t-1)")))

p_hl <- ggplot(plot_df,
               aes(x = estimate, y = forcats::fct_rev(iv_label))) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray55") +
  geom_errorbarh(aes(xmin = ci_low, xmax = ci_high),
                 height = 0.15, linewidth = 0.8, color = ACCENT) +
  geom_point(size = 3, color = ACCENT) +
  labs(
    title = "Hostility-level robustness",
    x = "Coefficient (95% CI)",
    y = NULL
  ) +
  theme_pres() +
  theme(plot.margin = margin(t = 10, r = 12, b = 14, l = 10))

save_pres("fig_ms_hostlev.png", p_hl, width = 6, height = 4)

# ============================================================
# Combined figure: state-level war-sample (full sample) + hostility-level
# specifications, each shown as a separate coefficient effect in one panel.
# Reads the war-sample results produced by 12_robustness_war_samples.R.
# ============================================================
ws_path <- file.path(OUTPUT_DIR, "robustness_war_samples.csv")
if (file.exists(ws_path)) {
  ws_state_full <- readr::read_csv(ws_path, show_col_types = FALSE) |>
    filter(level == "state", sample == "(1) Full sample") |>
    transmute(iv_label, estimate, se, ci_low, ci_high,
              group = "MID sample (full)")

  hl_state <- results |>                       # state HostLev rows from above
    transmute(iv_label, estimate, se, ci_low, ci_high,
              group = "Hostility level")

  # labels written by 11_robustness_war_samples.R (state IVs at t-1)
  lev_order <- c("MID onset (t-1)", "log(1 + MID onset) (t-1)",
                 "Fatal MID onset (t-1)", "log(1 + Fatal MID onset) (t-1)",
                 "Mean HostLev (onsets, t-1)", "log(1 + Mean HostLev) (onsets, t-1)")

  combined <- bind_rows(ws_state_full, hl_state) |>
    mutate(iv_label = factor(iv_label, levels = lev_order),
           group    = factor(group,
                             levels = c("MID sample (full)", "Hostility level")))

  p_combined <- ggplot(combined,
                       aes(x = estimate, y = forcats::fct_rev(iv_label),
                           color = group)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray55") +
    geom_errorbarh(aes(xmin = ci_low, xmax = ci_high),
                   height = 0.18, linewidth = 0.8) +
    geom_point(size = 3) +
    scale_color_manual(values = c("MID sample (full)" = ACCENT,
                                  "Hostility level"    = ACCENT2),
                       name = NULL) +
    labs(
      title = "State-level robustness: MID-sample and hostility-level specifications",
      x = "Coefficient (95% CI)", y = NULL
    ) +
    theme_pres() +
    theme(plot.margin = margin(t = 10, r = 12, b = 14, l = 10))

  save_pres("fig_ms_hostlev_war_samples_state.png", p_combined,
            width = 8.75, height = 5.5)
} else {
  cat("  (Combined figure skipped: robustness_war_samples.csv not found;",
      "run 12_robustness_war_samples.R first.)\n")
}

cat("\nDone.\n")
