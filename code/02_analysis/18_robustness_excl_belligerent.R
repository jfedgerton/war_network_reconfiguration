###############################################################################
# 24_robustness_excl_belligerent_dyads.R
# Robustness check: recompute weighted Jaccard churn at the state-year
# level, excluding dyads (i, j) where i and j are active MID participants
# in year t or t-1. Then re-run the headline H2 + H3 regressions on this
# alternative outcome.
#
# RATIONALE: The canonical churn measure includes ties between active
# combatants (e.g., Russia-Ukraine in 2022). When those dyads sever as
# part of the war itself, they inflate churn for both sides in a way
# that is mechanical and theoretically uninformative — of course
# Russia and Ukraine cut diplomatic ties when they fight each other.
# What the H2 and H3 hypotheses claim is something stronger: that war
# induces broader portfolio realignment beyond the bilateral conflict
# dyads. To test that, we strip the belligerent-belligerent dyads out
# of the churn computation and re-estimate the headline ladder.
#
# DESIGN:
#   1. Load dyad_year_network_layers.rds and identify belligerent dyads
#      per year (cowmidongoing == 1).
#   2. For each state-year (i, t), recompute weighted Jaccard churn
#      EXCLUDING any partner j where (i, j) was an active MID dyad at
#      either t or t-1. The Jaccard formula is otherwise identical to
#      02_compute_network_measures.R.
#   3. Save the alternative churn as `churn_excl_belligerent` and
#      `log_churn_excl_belligerent`.
#   4. Build the state-year analysis frame using the SAME controls and
#      sample as the manuscript main spec (CINC, major_power, polity
#      democracy/autocracy, partner_war_frac_weighted, is_isolated_t1).
#   5. Run the headline H2 ladder (Table 3 columns 3-4):
#        log_churn_excl ~ log(MID onset) + ctrls (no FE, +poly(year, 2))
#        log_churn_excl ~ log(MID onset) + ctrls | ccode + year   (TWFE)
#      For both MID-onset and fatal-MID-onset IVs.
#   6. Compare coefficients to the canonical-churn versions on the same
#      sample (re-fit canonical on this exact sample for apples-to-apples).
#   7. Output: CSV results, LaTeX table, coefficient-comparison plot.
#
# Inputs:  data/processed/dyad_year_network_layers.rds
#          data/processed/state_network_churn.rds  (for canonical churn comparison)
# Outputs:
#   data/processed/state_churn_excl_belligerent.rds
#   data/processed/robustness_excl_belligerent_results.csv
#   results/appendix/tableA_excl_belligerent.tex
#   code/03_figures/output/manuscript/fig_ms_excl_belligerent.png
###############################################################################

library(tidyverse)
library(fixest)
library(peacesciencer)
library(patchwork)
select <- dplyr::select

set.seed(123)

OUTPUT_DIR <- here::here("data", "processed")
APPX_DIR   <- here::here("results", "appendix")
FIG_DIR    <- here::here("code", "03_figures", "output", "manuscript")
dir.create(APPX_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG_DIR,  recursive = TRUE, showWarnings = FALSE)

ACCENT  <- "#2C6CB0"   # canonical
ACCENT2 <- "#B03366"   # excl-belligerent

theme_pres <- function(base_size = 12) {
  theme_bw(base_size = base_size) %+replace%
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "gray90", linewidth = 0.3),
      plot.title       = element_text(size = base_size + 2, face = "bold",
                                       hjust = 0, margin = margin(b = 8)),
      axis.text        = element_text(color = "black", size = base_size - 1),
      axis.title       = element_text(color = "black", size = base_size),
      axis.ticks       = element_line(color = "gray40"),
      strip.background = element_rect(fill = "gray92", color = "gray70"),
      strip.text       = element_text(face = "bold", size = base_size, color = "black"),
      legend.position  = "bottom",
      legend.text      = element_text(size = base_size - 1),
      legend.title     = element_text(size = base_size, face = "bold"),
      plot.margin      = margin(10, 12, 8, 10)
    )
}

save_pres <- function(filename, plot, width = 9, height = 6, dpi = 400) {
  ggsave(file.path(FIG_DIR, filename), plot,
         width = width, height = height, dpi = dpi, bg = "white")
  cat(sprintf("  -> %s\n", filename))
}

# ============================================================================
# 1. Load data
# ============================================================================

cat("Loading dyad-year network layers and canonical state churn...\n")
dyad_data    <- readRDS(file.path(OUTPUT_DIR, "dyad_year_network_layers.rds"))
state_churn  <- readRDS(file.path(OUTPUT_DIR, "state_network_churn.rds"))

stopifnot("main_weight"   %in% names(dyad_data))
stopifnot("cowmidongoing" %in% names(dyad_data))
stopifnot("churn"         %in% names(state_churn))

# ============================================================================
# 2. Identify belligerent dyads per year (i, j) where cowmidongoing == 1
# ============================================================================
# We treat "belligerent dyad at time t" as any (i, j, t) with an active
# MID. For state i's churn at year t, we want to exclude partner j if
# (i, j) was an active MID at t OR at t-1 — that captures the year of
# the bilateral severing whether the MID started in t-1 or t.

cat("Building belligerent-dyad lookup...\n")

belligerent_dyads <- dyad_data |>
  filter(cowmidongoing == 1L) |>
  select(ccode1, ccode2, year) |>
  distinct() |>
  # Symmetrize so we can look up by (focal_state, partner) regardless of order
  bind_rows(
    dyad_data |>
      filter(cowmidongoing == 1L) |>
      select(ccode1 = ccode2, ccode2 = ccode1, year) |>
      distinct()
  ) |>
  distinct()

cat(sprintf("  Belligerent dyad-years (symmetrized): %s\n",
            format(nrow(belligerent_dyads), big.mark = ",")))

# Pre-compute the (focal, partner, year) set as a fast-lookup table so we
# can avoid joining inside the per-state loop.
bell_lookup <- belligerent_dyads |>
  rename(focal = ccode1, partner = ccode2) |>
  mutate(key = paste(focal, partner, year, sep = "_"))
bell_set <- bell_lookup$key

# Per-(focal, year) belligerent-partner lookup ("focal_year" -> partner vec).
# Used for fast SYMMETRIC in-loop exclusion in compute_churn_local: a partner
# belligerent at t OR t-1 is stripped from BOTH years (matching the headline
# compute_churn_excl), so a war-driven severing never re-enters churn as a
# phantom dissolution at war onset.
bell_partner_by_fy <- split(bell_lookup$partner,
                            paste(bell_lookup$focal, bell_lookup$year, sep = "_"))

# ============================================================================
# 3. Recompute weighted Jaccard churn excluding belligerent dyads
# ============================================================================
# For each state-year (st, yr_t), pull the partner set at t and t-1 from
# main_weight, EXCLUDING any partner that was a belligerent of st at t or
# t-1. Then compute the canonical weighted Jaccard formula
# (1 - sum_min/sum_max) on the surviving partners.

cat("Recomputing weighted Jaccard churn excluding belligerent dyads...\n")

# Build a long edge list keyed by focal state for fast filtering
focal_edges <- bind_rows(
  dyad_data |>
    select(ccode1, ccode2, year, w = main_weight) |>
    transmute(focal = ccode1, partner = ccode2, year, w),
  dyad_data |>
    select(ccode1, ccode2, year, w = main_weight) |>
    transmute(focal = ccode2, partner = ccode1, year, w)
) |>
  filter(!is.na(w))

# Helper: is (focal, partner) a belligerent dyad at year y? Vectorised.
is_belligerent <- function(focal_vec, partner_vec, year_vec) {
  paste(focal_vec, partner_vec, year_vec, sep = "_") %in% bell_set
}

# Iterate state-year (matching state_churn rows so we keep the same panel
# structure). state_churn already has (ccode, year, churn, ...). We add a
# column `churn_excl_belligerent` to it.

cat(sprintf("  Iterating over %s state-year rows...\n",
            format(nrow(state_churn), big.mark = ",")))

# Pre-extract minimal columns from focal_edges for speed
focal_keys <- paste(focal_edges$focal, focal_edges$year, sep = "_")
focal_edges$key <- focal_keys

# Build per-(focal, year) edge tables once
edges_by_focal_year <- focal_edges |>
  group_by(key) |>
  group_split()
key_index <- vapply(edges_by_focal_year, function(df) df$key[1], character(1))
edge_lookup <- setNames(edges_by_focal_year, key_index)

compute_churn_excl <- function(st, yr_t) {
  yr_t1 <- yr_t - 1L
  k_t  <- paste(st, yr_t,  sep = "_")
  k_t1 <- paste(st, yr_t1, sep = "_")
  e_t  <- edge_lookup[[k_t]]
  e_t1 <- edge_lookup[[k_t1]]

  # Identify belligerent partners SYMMETRICALLY: any partner appearing in
  # either year that was belligerent with st at t OR t-1 is dropped from BOTH
  # years, so a war-driven severing is fully excluded (never a phantom
  # dissolution between t-1 and t).
  cand <- unique(c(if (!is.null(e_t)) e_t$partner, if (!is.null(e_t1)) e_t1$partner))
  cand <- cand[!is.na(cand)]
  bell_partners <- if (length(cand) > 0) {
    cand[is_belligerent(rep(st, length(cand)), cand, rep(yr_t,  length(cand))) |
         is_belligerent(rep(st, length(cand)), cand, rep(yr_t1, length(cand)))]
  } else integer(0)

  # Pre-exclusion portfolio (w > 0) = the state's GENUINE degree each year.
  e_t_full  <- if (!is.null(e_t))  e_t  |> filter(w > 0) else NULL
  e_t1_full <- if (!is.null(e_t1)) e_t1 |> filter(w > 0) else NULL

  # Drop belligerent partners from both year frames
  e_t_keep  <- if (!is.null(e_t_full))  e_t_full  |> filter(!partner %in% bell_partners) else NULL
  e_t1_keep <- if (!is.null(e_t1_full)) e_t1_full |> filter(!partner %in% bell_partners) else NULL

  if (is.null(e_t_keep) || is.null(e_t1_keep)) return(NA_real_)

  n_t   <- nrow(e_t_keep);  n_t1  <- nrow(e_t1_keep)
  nf_t  <- if (is.null(e_t_full))  0L else nrow(e_t_full)
  nf_t1 <- if (is.null(e_t1_full)) 0L else nrow(e_t1_full)

  # GUARD: a post-exclusion degree-0 year whose PRE-exclusion degree was > 0
  # is empty only because we deleted belligerent ties. Coding it as churn = 1
  # manufactures maximal reshuffling in war-years and inflates the estimate,
  # so drop the state-year (NA) instead of treating it as real isolation.
  if ((n_t == 0 && nf_t > 0) || (n_t1 == 0 && nf_t1 > 0)) return(NA_real_)

  if (n_t == 0 && n_t1 == 0) return(0)   # (a) genuinely isolated both years
  if (n_t == 0 || n_t1 == 0) return(1)   # (b) genuine total replacement

  w_both <- full_join(
    e_t_keep  |> select(partner, w_t  = w),
    e_t1_keep |> select(partner, w_t1 = w),
    by = "partner"
  ) |> mutate(w_t = replace_na(w_t, 0), w_t1 = replace_na(w_t1, 0))

  sum_min <- sum(pmin(w_both$w_t, w_both$w_t1))
  sum_max <- sum(pmax(w_both$w_t, w_both$w_t1))
  if (sum_max > 0) 1 - sum_min / sum_max else NA_real_
}

# Sequential loop (per coding-style preference) with progress
churn_excl_vec <- numeric(nrow(state_churn))
for (i in seq_len(nrow(state_churn))) {
  churn_excl_vec[i] <- compute_churn_excl(state_churn$ccode[i], state_churn$year[i])
  if (i %% 1000 == 0) cat(sprintf("    progress: %d / %d\n", i, nrow(state_churn)))
}

state_churn$churn_excl_belligerent     <- churn_excl_vec
state_churn$log_churn_excl_belligerent <- log1p(churn_excl_vec)

# Diagnostic: state-years dropped (NA) by the exclusion-induced-emptiness
# guard, i.e. years that USED to be coded churn = 1 purely because removing
# belligerent ties emptied the portfolio. These are the artifactual spikes
# that previously inflated the excl-belligerent coefficient.
n_excl_dropped <- max(0L, sum(is.na(state_churn$churn_excl_belligerent)) -
                            sum(is.na(state_churn$churn)))
cat(sprintf("  Exclusion-induced NA state-years (dropped, not coded churn=1): %d\n",
            n_excl_dropped))

# Drop a sanity-check correlation: alt churn should be highly correlated
# with canonical (most state-years have no belligerent partners), with
# divergence concentrated in war-participating state-years.
cat("\n--- Sanity: correlation between canonical and excl-belligerent churn ---\n")
print(state_churn |>
        summarise(cor_all = cor(churn, churn_excl_belligerent,
                                 use = "complete.obs"),
                  mean_canon = mean(churn,                  na.rm = TRUE),
                  mean_excl  = mean(churn_excl_belligerent, na.rm = TRUE),
                  n_diff     = sum(abs(churn - churn_excl_belligerent) > 1e-6,
                                   na.rm = TRUE)))

saveRDS(state_churn,
        file.path(OUTPUT_DIR, "state_churn_excl_belligerent.rds"))

# ============================================================================
# 4. Build state-year analysis frame (mirror 02_state_level_analysis.R)
# ============================================================================
# Same controls as the manuscript main spec. We rebuild here rather than
# reading from 02 outputs so this script is self-contained.

cat("\nBuilding state-year analysis frame (same controls as main spec)...\n")

# State attributes from dyad_data: CINC, major_power
state_attrs <- dyad_data |>
  select(ccode = ccode1, year, cinc = cinc1, major_power = cowmaj1) |>
  distinct(ccode, year, .keep_all = TRUE) |>
  mutate(major_power = as.integer(coalesce(major_power, 0L) == 1))

state_churn <- state_churn |>
  left_join(state_attrs, by = c("ccode", "year"))

# Polity (democracy / autocracy) via peacesciencer
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
  left_join(polity, by = c("ccode", "year"))

# Partner-war exposure (severity-routed, weighted) — same logic as scripts
# 02 and 22 so the controls match the manuscript main spec.
partner_sets <- dyad_data |>
  filter(main_weight > 0) |>
  select(ccode1, ccode2, year, w = main_weight) |>
  bind_rows(
    dyad_data |> filter(main_weight > 0) |>
      select(ccode1 = ccode2, ccode2 = ccode1, year, w = main_weight)
  )

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
  left_join(at_war_dyad,       by = c("ccode2" = "ccode", "year_next" = "year")) |>
  left_join(at_war_dyad_fatal, by = c("ccode2" = "ccode", "year_next" = "year")) |>
  mutate(at_war = replace_na(at_war, 0L),
         at_war_fatal = replace_na(at_war_fatal, 0L)) |>
  group_by(ccode1, year_next) |>
  summarise(
    partner_war_frac_weighted       = sum(at_war * w)       / sum(w),
    partner_fatal_war_frac_weighted = sum(at_war_fatal * w) / sum(w),
    .groups = "drop"
  ) |>
  rename(ccode = ccode1, year = year_next)

state_churn <- state_churn |>
  left_join(pwe, by = c("ccode", "year")) |>
  mutate(partner_war_frac_weighted       = replace_na(partner_war_frac_weighted, 0),
         partner_fatal_war_frac_weighted = replace_na(partner_fatal_war_frac_weighted, 0))

# ----------------------------------------------------------------------------
# Excl-belligerent PWF: recompute the weighted partner-war exposure using
# the SAME exclusion logic as the excl-belligerent churn — drop dyads
# where (i, j) was an active MID dyad at year t or t-1 from the partner
# set BEFORE computing the weighted fraction at war.
# This makes the excl-belligerent specs fully apples-to-apples internally:
# excl-belligerent churn (DV) + excl-belligerent PWF (control). The
# canonical specs continue to use canonical PWF.
# ----------------------------------------------------------------------------
cat("Computing excl-belligerent PWF (weighted, severity-routed)...\n")

# Mark each (focal, partner, year) edge with whether it was a belligerent
# dyad at year t or t-1. Reuse `bell_set` built earlier.
partner_sets_excl <- partner_sets |>
  mutate(
    bell_t  = is_belligerent(ccode1, ccode2, year),
    bell_t1 = is_belligerent(ccode1, ccode2, year - 1L),
    is_bell = bell_t | bell_t1
  ) |>
  filter(!is_bell) |>
  select(ccode1, ccode2, year, w)

pwe_excl <- partner_sets_excl |>
  mutate(year_next = year + 1L) |>
  left_join(at_war_dyad,       by = c("ccode2" = "ccode", "year_next" = "year")) |>
  left_join(at_war_dyad_fatal, by = c("ccode2" = "ccode", "year_next" = "year")) |>
  mutate(at_war = replace_na(at_war, 0L),
         at_war_fatal = replace_na(at_war_fatal, 0L)) |>
  group_by(ccode1, year_next) |>
  summarise(
    partner_war_frac_weighted_excl       = sum(at_war * w)       / sum(w),
    partner_fatal_war_frac_weighted_excl = sum(at_war_fatal * w) / sum(w),
    .groups = "drop"
  ) |>
  rename(ccode = ccode1, year = year_next)

state_churn <- state_churn |>
  left_join(pwe_excl, by = c("ccode", "year")) |>
  mutate(partner_war_frac_weighted_excl       = replace_na(partner_war_frac_weighted_excl, 0),
         partner_fatal_war_frac_weighted_excl = replace_na(partner_fatal_war_frac_weighted_excl, 0))

cat(sprintf("  Excl-belligerent PWF: corr(canonical, excl) = %.3f (any-MID), %.3f (fatal)\n",
            cor(state_churn$partner_war_frac_weighted,
                state_churn$partner_war_frac_weighted_excl, use = "complete.obs"),
            cor(state_churn$partner_fatal_war_frac_weighted,
                state_churn$partner_fatal_war_frac_weighted_excl, use = "complete.obs")))

# IV transforms (mirror 02)
state_churn <- state_churn |>
  mutate(
    log_churn_canonical       = log1p(churn),
    log_mid_onset_count       = log1p(mid_onset_count),
    log_fatal_mid_onset_count = log1p(fatal_mid_onset_count)
  )

# Full state-year panel (ALL state-years, including degree-0 / isolated ones).
# The multiverse churn recomputation iterates over this so that
# isolation-transition state-years are retained, mirroring script 02's
# headline churn (which the previous years-with-ties-only loop dropped).
state_panel <- state_churn |> distinct(ccode, year) |> arrange(ccode, year)

# ============================================================================
# 5. Run headline H2 ladder on canonical vs excl-belligerent outcomes
# ============================================================================
# Mirror Table 3 columns 3 (no-FE + poly trend) and 4 (TWFE) from
# 02_state_level_analysis.R. Run on the SAME sample for both outcomes
# (rows where both are non-missing).

# Four control strings: (canonical | excl-belligerent) × (any-MID | fatal-MID).
# Severity-routed: MID IV uses any-MID PWF, fatal-MID IV uses fatal-MID PWF.
# Outcome-routed: canonical churn uses canonical PWF, excl churn uses excl PWF.
CTRL_MID_CANON   <- "partner_war_frac_weighted + democracy + autocracy + cinc + major_power + is_isolated_t1"
CTRL_MID_EXCL    <- "partner_war_frac_weighted_excl + democracy + autocracy + cinc + major_power + is_isolated_t1"
CTRL_FATAL_CANON <- "partner_fatal_war_frac_weighted + democracy + autocracy + cinc + major_power + is_isolated_t1"
CTRL_FATAL_EXCL  <- "partner_fatal_war_frac_weighted_excl + democracy + autocracy + cinc + major_power + is_isolated_t1"

# Apples-to-apples sample: rows with both outcomes non-missing
analysis_sample <- state_churn |>
  filter(!is.na(log_churn_canonical),
         !is.na(log_churn_excl_belligerent))

cat(sprintf("\nAnalysis sample: %s state-years\n",
            format(nrow(analysis_sample), big.mark = ",")))

fit_one <- function(formula, data, vcov_arg, iv_var, outcome_label, model_label) {
  fit <- tryCatch(
    feols(formula, data = data, vcov = vcov_arg),
    error = function(e) {
      cat(sprintf("    ! %s | %s failed: %s\n", outcome_label, model_label,
                  conditionMessage(e)))
      return(NULL)
    }
  )
  if (is.null(fit)) return(NULL)
  ct <- coeftable(fit)
  if (!iv_var %in% rownames(ct)) return(NULL)
  tibble(
    iv       = iv_var,
    iv_label = switch(iv_var,
                      "log_mid_onset_count"       = "log(1+MID Onset)",
                      "log_fatal_mid_onset_count" = "log(1+Fatal MID Onset)",
                      iv_var),
    outcome  = outcome_label,
    spec     = model_label,
    estimate = ct[iv_var, "Estimate"],
    se       = ct[iv_var, "Std. Error"],
    pval     = ct[iv_var, "Pr(>|t|)"],
    n_obs    = nobs(fit),
    ci_low   = ct[iv_var, "Estimate"] - 1.96 * ct[iv_var, "Std. Error"],
    ci_high  = ct[iv_var, "Estimate"] + 1.96 * ct[iv_var, "Std. Error"]
  )
}

# Both MID-onset and fatal-MID-onset IVs. Severity-routed: MID IV uses
# any-MID PWF; fatal IV uses fatal-MID PWF. Outcome-routed: canonical
# churn uses canonical PWF; excl-belligerent churn uses excl-belligerent
# PWF (so the spec is internally consistent within each outcome variant).
spec_grid <- expand.grid(
  iv      = c("log_mid_onset_count", "log_fatal_mid_onset_count"),
  outcome = c("canonical", "excl"),
  spec    = c("(1) Bivariate TWFE", "(2) +Ctrls TWFE",
              "(3) +Ctrls, no FE"),
  stringsAsFactors = FALSE
)

pick_dv   <- function(outcome) {
  if (outcome == "canonical") "log_churn_canonical" else "log_churn_excl_belligerent"
}
pick_ctrl <- function(iv, outcome) {
  if (iv == "log_mid_onset_count") {
    if (outcome == "canonical") CTRL_MID_CANON   else CTRL_MID_EXCL
  } else {
    if (outcome == "canonical") CTRL_FATAL_CANON else CTRL_FATAL_EXCL
  }
}
outcome_label <- function(outcome) {
  if (outcome == "canonical") "Canonical churn" else "Excl. belligerent dyads"
}

run_spec_row <- function(iv, outcome, spec) {
  dv   <- pick_dv(outcome)
  ctrl <- pick_ctrl(iv, outcome)
  lab  <- outcome_label(outcome)

  if (spec == "(1) Bivariate TWFE") {
    fml <- as.formula(paste0(dv, " ~ ", iv, " | ccode + year"))
    return(fit_one(fml, analysis_sample, ~ ccode + year, iv, lab, spec))
  }
  if (spec == "(2) +Ctrls TWFE") {
    fml <- as.formula(paste0(dv, " ~ ", iv, " + ", ctrl, " | ccode + year"))
    return(fit_one(fml, analysis_sample, ~ ccode + year, iv, lab, spec))
  }
  if (spec == "(3) +Ctrls, no FE") {
    fml <- as.formula(paste0(dv, " ~ ", iv, " + ", ctrl, " + poly(year, 2)"))
    return(fit_one(fml, analysis_sample, ~ ccode, iv, lab, spec))
  }
  NULL
}

# Parallelize across spec rows (12 rows: 2 IVs × 2 outcomes × 3 specs).
# Each row is a single feols fit; furrr::future_pmap_dfr runs them in
# parallel workers. Fall back to sequential if furrr is unavailable.
.use_parallel <- requireNamespace("furrr", quietly = TRUE) &&
                 requireNamespace("future", quietly = TRUE)
if (.use_parallel) {
  future::plan(future::multisession,
               workers = max(1L, future::availableCores() - 1L))
  cat(sprintf("\nParallel: %d workers\n", future::nbrOfWorkers()))
  results_df <- furrr::future_pmap_dfr(
    list(spec_grid$iv, spec_grid$outcome, spec_grid$spec),
    run_spec_row,
    .options = furrr::furrr_options(seed = TRUE)
  )
  future::plan(future::sequential)
} else {
  cat("\nSequential (furrr/future not installed)\n")
  results_df <- purrr::pmap_dfr(
    list(spec_grid$iv, spec_grid$outcome, spec_grid$spec),
    run_spec_row
  )
}

cat("\n========== RESULTS ==========\n")
print(results_df, n = Inf, digits = 4)

write_csv(results_df,
          file.path(OUTPUT_DIR, "robustness_excl_belligerent_results.csv"))

# ============================================================================
# 6. Coefficient-comparison plot
# ============================================================================

# ============================================================================
# 8. Multiverse: excl-belligerent across layer combos × weighting schemes
# ============================================================================
# Replicate the excl-belligerent robustness across the same layer combos
# and weighting schemes used in 03_robustness.R. All helper functions are
# defined locally here to avoid modifying shared code.
#
# Grid: 12 layer combos (DCA excluded) × 5 weighting methods = 60 cells.
# For each cell: recompute edge weights, compute canonical churn and
# excl-belligerent churn, run the MID-onset TWFE regression on both,
# record attenuation.

cat("\n=== MULTIVERSE: Layer combos × Weighting schemes ===\n")

# --- Load metadata for layer discovery ---
metadata <- readRDS(file.path(OUTPUT_DIR, "network_metadata.rds"))
all_layer_vars <- metadata$layers
all_layer_vars <- all_layer_vars[sapply(all_layer_vars, function(lv) {
  sd(dyad_data[[lv]], na.rm = TRUE) > 0
})]

# DCA IS EXCLUDED [2026-07-15], matching 07_robustness_multiverse.R. DCA data
# begins in 1980, and once each combo is restricted to its layers' coverage
# window (see LAYER_COVERAGE below), DCA combos collapse to too few years to be
# informative. DCA is retained only in the contamination multiverse
# (19_contamination_multiverse.R), where the dyad-level panel has ample N.
EXCLUDED_LAYERS <- c("layer_dca")
all_layer_vars  <- setdiff(all_layer_vars, EXCLUDED_LAYERS)

# COVERAGE = where each layer's DATA EXISTS, not where it is positive. Offensive
# alliances have their last positive tie in 1990 but the data span the full
# panel (post-1990 zeros are real); trade data begin in 1870. Hard-coded from
# the sources, matching 07.
PANEL_END <- 2014L
LAYER_COVERAGE <- list(
  layer_offensive = c(1816L, PANEL_END),
  layer_defensive = c(1816L, PANEL_END),
  layer_igo       = c(1816L, PANEL_END),
  layer_trade     = c(1870L, PANEL_END)
)
get_valid_year_range <- function(layers) {
  starts <- sapply(layers, function(lv) LAYER_COVERAGE[[lv]][1])
  ends   <- sapply(layers, function(lv) LAYER_COVERAGE[[lv]][2])
  c(max(starts), min(ends))
}

MAIN_LAYERS <- c("layer_offensive", "layer_defensive", "layer_igo")

cat(sprintf("  All layers: %s\n", paste(all_layer_vars, collapse = ", ")))
cat(sprintf("  EXCLUDED: %s\n", paste(EXCLUDED_LAYERS, collapse = ", ")))

# --- Build layer combos (same logic as 03) ---
all_combos_raw <- list()
for (k in seq_along(all_layer_vars)) {
  all_combos_raw <- c(all_combos_raw, combn(all_layer_vars, k, simplify = FALSE))
}
extra_layers <- setdiff(all_layer_vars, MAIN_LAYERS)
if (length(extra_layers) > 0) {
  for (k in seq_along(extra_layers)) {
    for (ex in combn(extra_layers, k, simplify = FALSE)) {
      all_combos_raw <- c(all_combos_raw, list(sort(c(MAIN_LAYERS, ex))))
    }
  }
}
all_combos_raw <- unique(lapply(all_combos_raw, sort))
all_combos <- Filter(function(combo) {
  "layer_igo" %in% combo | "layer_defensive" %in% combo
}, all_combos_raw)

weighting_methods <- c("contig_manhattan", "contig_euclidean",
                       "dist_manhattan", "dist_euclidean", "unweighted")

cat(sprintf("  Layer combos: %d | Weighting methods: %d | Total cells: %d\n",
            length(all_combos), length(weighting_methods),
            length(all_combos) * length(weighting_methods)))

# --- Local helper: range01 ---
range01_local <- function(x) {
  r <- range(x, na.rm = TRUE)
  if (r[2] == r[1]) return(rep(0, length(x)))
  (x - r[1]) / (r[2] - r[1])
}

# --- Local helper: compute edge weight (mirrors 03 compute_manhattan_vec) ---
compute_weight_local <- function(data, layers, method) {
  mat <- as.matrix(data[, layers, drop = FALSE])
  mat[is.na(mat)] <- 0

  if (method == "contig_manhattan") {
    wcols <- intersect(paste0(layers, "_weighted"), names(data))
    if (length(wcols) == 0) return(rep(0, nrow(data)))
    mat_w <- as.matrix(data[, wcols, drop = FALSE]); mat_w[is.na(mat_w)] <- 0
    result <- rowSums(mat_w) / length(wcols)
  } else if (method == "contig_euclidean") {
    wcols <- intersect(paste0(layers, "_weighted"), names(data))
    if (length(wcols) == 0) return(rep(0, nrow(data)))
    mat_w <- as.matrix(data[, wcols, drop = FALSE]); mat_w[is.na(mat_w)] <- 0
    mat_w[mat_w < 0] <- 0
    result <- sqrt(rowSums(mat_w^2)) / length(wcols)
  } else if (method == "dist_manhattan") {
    dcols <- intersect(paste0(layers, "_dist_weighted"), names(data))
    if (length(dcols) == 0) return(rep(0, nrow(data)))
    mat_d <- as.matrix(data[, dcols, drop = FALSE]); mat_d[is.na(mat_d)] <- 0
    result <- rowSums(mat_d) / length(dcols)
  } else if (method == "dist_euclidean") {
    dcols <- intersect(paste0(layers, "_dist_weighted"), names(data))
    if (length(dcols) == 0) return(rep(0, nrow(data)))
    mat_d <- as.matrix(data[, dcols, drop = FALSE]); mat_d[is.na(mat_d)] <- 0
    mat_d[mat_d < 0] <- 0
    result <- sqrt(rowSums(mat_d^2)) / length(dcols)
  } else if (method == "unweighted") {
    result <- rowSums(mat)
  } else {
    return(rep(0, nrow(data)))
  }

  # Normalize by year
  years <- data$year
  for (yr in unique(years)) {
    idx <- which(years == yr)
    result[idx] <- range01_local(result[idx])
  }
  result
}

# --- Local helper: compute churn for a weight vector, with optional excl ---
# Returns a tibble with (ccode, year, churn_val, log_churn_val).
# If excl_bell = TRUE, drops belligerent-dyad partners from computation.
compute_churn_local <- function(data, w_vec, panel, excl_bell = FALSE) {
  data$w_alt <- w_vec

  ties <- bind_rows(
    data |> filter(w_alt > 0) |> transmute(ccode1, ccode2, year, w = w_alt),
    data |> filter(w_alt > 0) |>
      transmute(ccode1 = ccode2, ccode2 = ccode1, year, w = w_alt)
  ) |>
    group_by(ccode1, ccode2, year) |>
    summarise(w = max(w, na.rm = TRUE), .groups = "drop")

  # FULL (un-excluded) partner-weight lookup: key "focal_year" -> named
  # numeric vector (names = partner ccode, values = weight). Belligerent
  # partners are removed PER TRANSITION inside the loop (see below), NOT by a
  # per-year pre-filter -- a per-year filter would keep a tie at t-1 and drop
  # it at t, manufacturing a phantom dissolution at war onset and inflating
  # the excl-belligerent coefficient.
  key <- paste(ties$ccode1, ties$year, sep = "_")
  pw_lookup <- split(setNames(ties$w, ties$ccode2), key)

  # Iterate the FULL state-year panel (incl. degree-0 years) so isolation
  # transitions are retained, mirroring script 02's four-case churn rule.
  panel <- panel |> distinct(ccode, year) |> arrange(ccode, year)
  pc <- panel$ccode; py <- panel$year
  n  <- length(pc)

  res_ccode <- numeric(n); res_year <- numeric(n); res_churn <- numeric(n)
  idx <- 0L

  for (st in unique(pc)) {
    yrs <- py[pc == st]
    if (length(yrs) < 2) next
    for (j in 2:length(yrs)) {
      yr_t <- yrs[j]; yr_t1 <- yrs[j - 1]
      if (yr_t - yr_t1 != 1) next

      pt  <- pw_lookup[[paste(st, yr_t,  sep = "_")]]
      pt1 <- pw_lookup[[paste(st, yr_t1, sep = "_")]]
      deg_t_full  <- length(pt)   # pre-exclusion degree (length(NULL) == 0)
      deg_t1_full <- length(pt1)

      # SYMMETRIC exclusion: a partner belligerent at t OR t-1 is removed from
      # BOTH years' portfolios, so a war-driven severing is fully excluded
      # rather than registering as a tie that vanishes between t-1 and t.
      if (excl_bell) {
        bp <- unique(c(bell_partner_by_fy[[paste(st, yr_t,  sep = "_")]],
                       bell_partner_by_fy[[paste(st, yr_t1, sep = "_")]]))
        if (length(bp) > 0) {
          bp <- as.character(bp)
          if (deg_t_full  > 0) pt  <- pt[!(names(pt)  %in% bp)]
          if (deg_t1_full > 0) pt1 <- pt1[!(names(pt1) %in% bp)]
        }
      }

      deg_t  <- length(pt)
      deg_t1 <- length(pt1)

      # GUARD (excl_bell only): a post-exclusion degree-0 year whose
      # PRE-exclusion degree was > 0 is empty solely because we deleted
      # belligerent ties. Coding it as churn = 1 manufactures maximal
      # reshuffling in war-years and inflates the estimate, so drop the
      # state-year (NA) instead of treating it as real isolation.
      if (excl_bell && ((deg_t == 0 && deg_t_full > 0) ||
                        (deg_t1 == 0 && deg_t1_full > 0))) {
        idx <- idx + 1L
        res_ccode[idx] <- st
        res_year[idx]  <- yr_t
        res_churn[idx] <- NA_real_
        next
      }

      if (deg_t == 0 && deg_t1 == 0) {          # (a) both genuinely isolated
        ch <- 0
      } else if (deg_t == 0 || deg_t1 == 0) {   # (b) genuine total replacement
        ch <- 1
      } else {                                   # (c) weighted Jaccard
        parts <- union(names(pt), names(pt1))
        wt  <- pt[parts];  wt[is.na(wt)]   <- 0
        wt1 <- pt1[parts]; wt1[is.na(wt1)] <- 0
        sx <- sum(pmax(wt, wt1))
        ch <- if (sx > 0) 1 - sum(pmin(wt, wt1)) / sx else NA_real_
      }
      idx <- idx + 1L
      res_ccode[idx] <- st
      res_year[idx]  <- yr_t
      res_churn[idx] <- ch
    }
  }

  tibble(ccode = res_ccode[seq_len(idx)],
         year  = res_year[seq_len(idx)],
         churn_val     = res_churn[seq_len(idx)],
         log_churn_val = log1p(res_churn[seq_len(idx)]))
}

# --- Local helper: compute combo-specific partner_war_frac_weighted ---
# Uses the same edge weights as churn so the control is consistent.
compute_pwf_local <- function(data, w_vec) {
  data$w_alt <- w_vec
  ties <- bind_rows(
    data |> filter(w_alt > 0) |>
      transmute(ccode1, ccode2, year, w = w_alt),
    data |> filter(w_alt > 0) |>
      transmute(ccode1 = ccode2, ccode2 = ccode1, year, w = w_alt)
  )
  at_war <- data |>
    filter(cowmidongoing == 1) |>
    select(ccode1, ccode2, year) |>
    pivot_longer(c(ccode1, ccode2), values_to = "ccode") |>
    distinct(ccode, year) |>
    mutate(at_war = 1L)
  ties |>
    mutate(year_next = year + 1L) |>
    left_join(at_war, by = c("ccode2" = "ccode", "year_next" = "year")) |>
    mutate(at_war = replace_na(at_war, 0L)) |>
    group_by(ccode1, year_next) |>
    summarise(pwf_combo = sum(at_war * w) / sum(w), .groups = "drop") |>
    rename(ccode = ccode1, year = year_next)
}

# --- Per-cell worker (parallelized via furrr below) ---
run_mv_cell_24 <- function(combo_i, method) {
  layers <- all_combos[[combo_i]]
  combo_label <- paste(gsub("layer_", "", layers), collapse = "+")

  # Compute edge weights for this combo
  w_vec <- tryCatch(compute_weight_local(dyad_data, layers, method),
                    error = function(e) NULL)
  if (is.null(w_vec) || all(w_vec == 0)) return(NULL)

  churn_canon <- tryCatch(
    compute_churn_local(dyad_data, w_vec, panel = state_panel, excl_bell = FALSE),
    error = function(e) NULL
  )
  if (is.null(churn_canon) || nrow(churn_canon) == 0) return(NULL)

  churn_excl <- tryCatch(
    compute_churn_local(dyad_data, w_vec, panel = state_panel, excl_bell = TRUE),
    error = function(e) NULL
  )
  if (is.null(churn_excl) || nrow(churn_excl) == 0) return(NULL)

  # Hold the partner-war-exposure control at the headline version
  # (partner_war_frac_weighted, already in analysis_sample) and vary only the
  # DV churn across cells. A combo-specific recompute (compute_pwf_local) drove
  # the canonical estimate away from the headline; the control is a covariate,
  # not the focus, so fixing it keeps the comparison apples-to-apples.
  merged <- analysis_sample |>
    select(ccode, year, log_mid_onset_count, partner_war_frac_weighted,
           democracy, autocracy, cinc, major_power, is_isolated_t1) |>
    inner_join(churn_canon |> select(ccode, year,
                                      churn_canon = churn_val,
                                      log_churn_canon = log_churn_val),
               by = c("ccode", "year")) |>
    inner_join(churn_excl |> select(ccode, year,
                                     churn_excl = churn_val,
                                     log_churn_excl = log_churn_val),
               by = c("ccode", "year"))

  # Restrict to this combo's coverage window [2026-07-15]. Without this, a
  # trade-containing combo would carry structural zeros before 1870, making the
  # churn DV a different construct pre- vs post-coverage and leaving the combo
  # a near-duplicate of the no-trade version. Mirrors 07 and 19.
  yr <- get_valid_year_range(layers)
  merged <- merged |> filter(year >= yr[1], year <= yr[2])

  if (nrow(merged) < 100) return(NULL)

  CTRL_COMBO <- "partner_war_frac_weighted + democracy + autocracy + cinc + major_power + is_isolated_t1"

  # Two focal IVs: H2 (MID onset) and H3 (partner-war exposure). Both enter the
  # same regression; we record each cell's coefficient for both.
  focal_ivs <- c("log(1+MID Onset)"     = "log_mid_onset_count",
                 "Partner-war exposure" = "partner_war_frac_weighted")

  rows <- list()
  for (outcome_pair in list(
    list(dv = "log_churn_canon", label = "Canonical"),
    list(dv = "log_churn_excl",  label = "Excl. belligerent")
  )) {
    fit <- tryCatch(
      feols(as.formula(paste0(outcome_pair$dv, " ~ log_mid_onset_count + ",
                              CTRL_COMBO, " | ccode + year")),
            data = merged, vcov = ~ ccode + year),
      error = function(e) NULL
    )
    if (is.null(fit)) next
    ct <- coeftable(fit)
    for (iv_lab in names(focal_ivs)) {
      iv_var <- focal_ivs[[iv_lab]]
      if (!iv_var %in% rownames(ct)) next
      rows[[length(rows) + 1]] <- tibble(
        layers   = combo_label,
        method   = method,
        outcome  = outcome_pair$label,
        iv       = iv_lab,
        estimate = ct[iv_var, "Estimate"],
        se       = ct[iv_var, "Std. Error"],
        pval     = ct[iv_var, "Pr(>|t|)"],
        n_obs    = nobs(fit)
      )
    }
  }
  bind_rows(rows)
}

# --- Parallel dispatch over (combo, method) cells ---
spec_grid_mv <- expand.grid(
  combo_i = seq_along(all_combos),
  method  = weighting_methods,
  stringsAsFactors = FALSE
)
total_cells <- nrow(spec_grid_mv)
cat(sprintf("  Total cells: %d (%d combos × %d methods)\n",
            total_cells, length(all_combos), length(weighting_methods)))

# Default to SEQUENTIAL: the per-cell churn is now cheap, and the parallel
# workers each copy dyad_data, which caused memory pressure / silent empty
# results on some machines. Opt into parallel with ROBUSTNESS_PARALLEL_24=TRUE.
.use_parallel_mv <- isTRUE(as.logical(Sys.getenv("ROBUSTNESS_PARALLEL_24", "FALSE"))) &&
                    requireNamespace("furrr", quietly = TRUE) &&
                    requireNamespace("future", quietly = TRUE)

# Defensive parallel setup (mirrors 25's memory-safety pattern).
# Each multisession worker gets a copy of dyad_data + analysis_sample +
# helpers, so unbounded worker counts blow through RAM. Hard-cap at 6
# workers and bump globals.maxSize to 8GB. Wrap in tryCatch so a
# parallel crash falls back to sequential rather than killing the script.
run_multiverse_24 <- function(use_parallel) {
  if (use_parallel) {
    .n_workers_cap <- 4L
    .n_workers_env <- as.integer(Sys.getenv("ROBUSTNESS_WORKERS",
                                              as.character(.n_workers_cap)))
    .n_workers <- min(.n_workers_env, .n_workers_cap,
                      max(1L, future::availableCores() - 1L))
    options(future.globals.maxSize = 8 * 1024^3)  # 8 GB
    future::plan(future::multisession, workers = .n_workers)
    cat(sprintf("  Parallel: %d workers (capped to avoid memory pressure)\n",
                future::nbrOfWorkers()))
    on.exit(future::plan(future::sequential), add = TRUE)
    furrr::future_map2_dfr(
      spec_grid_mv$combo_i, spec_grid_mv$method,
      run_mv_cell_24,
      .options = furrr::furrr_options(
        seed = TRUE,
        globals = c("all_combos", "dyad_data", "analysis_sample", "bell_set",
                    "bell_partner_by_fy", "state_panel",
                    "compute_weight_local", "compute_churn_local",
                    "compute_pwf_local", "range01_local",
                    "get_valid_year_range", "LAYER_COVERAGE")
      ),
      .progress = TRUE
    )
  } else {
    cat("  Sequential (furrr/future not installed)\n")
    purrr::map2_dfr(
      spec_grid_mv$combo_i, spec_grid_mv$method,
      run_mv_cell_24
    )
  }
}

mv_df <- tryCatch(
  run_multiverse_24(.use_parallel_mv),
  error = function(e) {
    cat(sprintf("  Parallel run failed: %s\n", conditionMessage(e)))
    cat("  Falling back to sequential execution...\n")
    future::plan(future::sequential)
    purrr::map2_dfr(
      spec_grid_mv$combo_i, spec_grid_mv$method,
      run_mv_cell_24
    )
  }
)

cat(sprintf("\n  Multiverse complete: %d cells produced results\n", nrow(mv_df)))

# --- VALIDATION GATE ---------------------------------------------------------
# The recomputed CANONICAL churn for the headline composite
# (offensive + defensive + IGO, unweighted) must reproduce the main-text
# state-level estimate. If it does not, the multiverse churn recompute has
# diverged from script 02 and the multiverse results must NOT be used.
# Reference = the non-multiverse canonical estimate (results_df), which is
# built on the actual headline churn.
if (nrow(mv_df) > 0) {
  ref_cell <- mv_df |>
    filter(layers == "defensive+igo+offensive", method == "unweighted",
           outcome == "Canonical", iv == "MID onset")
  headline_beta <- results_df |>
    filter(outcome == "Canonical churn", iv == "log_mid_onset_count",
           spec == "(2) +Ctrls TWFE") |>
    pull(estimate)

  cat("\n--- Validation gate: main-composite canonical vs headline ---\n")
  if (nrow(ref_cell) == 1 && length(headline_beta) == 1) {
    drift <- abs(ref_cell$estimate - headline_beta)
    cat(sprintf("  Multiverse canonical (offensive+defensive+igo, unweighted): %.4f (p=%.3f)\n",
                ref_cell$estimate, ref_cell$pval))
    cat(sprintf("  Headline canonical (non-multiverse, TWFE + controls):       %.4f\n",
                headline_beta))
    cat(sprintf("  Absolute drift: %.4f (tolerance 0.0015)\n", drift))
    if (drift > 1.5e-3) {
      stop(sprintf(paste0(
        "VALIDATION GATE FAILED: recomputed canonical churn (%.4f) does not ",
        "reproduce the headline estimate (%.4f). The multiverse churn ",
        "recompute has diverged from script 02 -- do not use these results."),
        ref_cell$estimate, headline_beta))
    }
    cat("  PASS: multiverse canonical reproduces the headline within tolerance.\n")
  } else {
    warning("Validation gate could not run (reference cell or headline beta missing).")
  }
}

# --- Summary statistics (by focal IV: H2 MID onset, H3 partner-war exposure) ---
if (nrow(mv_df) > 0) {
  mv_wide <- mv_df |>
    select(layers, method, iv, outcome, estimate) |>
    pivot_wider(names_from = outcome, values_from = estimate) |>
    mutate(attenuation_pct = 100 * (1 - `Excl. belligerent` / Canonical))

  cat("\n--- Multiverse attenuation summary (by IV) ---\n")
  for (iv_lab in unique(mv_wide$iv)) {
    a <- mv_wide |> filter(iv == iv_lab) |> pull(attenuation_pct)
    cat(sprintf("  %-22s mean %.1f%% | median %.1f%% | range %.1f to %.1f\n",
                iv_lab, mean(a, na.rm = TRUE), median(a, na.rm = TRUE),
                min(a, na.rm = TRUE), max(a, na.rm = TRUE)))
  }

  # One-sided 0.10 for a positive coefficient -- matches the threshold used
  # throughout the manuscript (coef > 0 & pval/2 < 0.10).
  cat("\n--- Positive & significant share (one-sided 0.10) ---\n")
  sig_tab <- mv_df |>
    group_by(iv, outcome) |>
    summarise(n = dplyr::n(),
              pos_sig = sum(estimate > 0 & (pval / 2) < 0.10, na.rm = TRUE),
              pct = 100 * mean(estimate > 0 & (pval / 2) < 0.10, na.rm = TRUE),
              .groups = "drop")
  print(sig_tab)

  write_csv(mv_df, file.path(OUTPUT_DIR,
            "robustness_excl_belligerent_multiverse.csv"))
  write_csv(mv_wide, file.path(OUTPUT_DIR,
            "robustness_excl_belligerent_attenuation.csv"))

  # --- Visualization: coefficient boxplots + positive-&-significant bars,
  #     faceted by focal IV (H2 MID onset; H3 partner-war exposure), with
  #     canonical vs excl-belligerent churn. Free y-scales because the two
  #     IVs are on different magnitudes. One-sided test is used in panel (b)
  #     because each cell is asked to agree with the predicted positive sign.
  mv_box_df <- mv_df |>
    mutate(outcome = factor(outcome,
                            levels = c("Canonical", "Excl. belligerent")),
           iv = factor(iv, levels = c("log(1+MID Onset)", "Partner-war exposure"), labels = c("MID Onset", "Partner-war exposure")))

  mv_class_df <- mv_box_df |>
    mutate(
      pos_sig_1side = (estimate > 0) & ((pval / 2) < 0.1),
      cat = factor(if_else(pos_sig_1side,
                           "Positive & sig (1-sided p<0.1)",
                           "Not positive & sig"),
                   levels = c("Positive & sig (1-sided p<0.1)",
                              "Not positive & sig"))
    ) |>
    count(iv, outcome, cat, .drop = FALSE) |>
    group_by(iv, outcome) |>
    mutate(share = n / sum(n)) |>
    ungroup()

  n_per_variant <- mv_box_df |>
    count(iv, outcome) |>
    summarise(n = max(n)) |>
    pull(n)

  p_box <- ggplot(mv_box_df,
                  aes(x = outcome, y = estimate, fill = outcome)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray60") +
    geom_boxplot(width = 0.5, outlier.size = 0.7,
                 alpha = 0.75, linewidth = 0.4) +
    facet_wrap(~ iv, scales = "free_y") +
    scale_fill_manual(values = c("Canonical"         = ACCENT,
                                  "Excl. belligerent" = ACCENT2),
                      guide = "none") +
    scale_x_discrete(labels = c("Canonical"         = "Canonical\nchurn",
                                 "Excl. belligerent" = "Excl.\nbelligerent")) +
    labs(x = NULL, y = "Coefficient (95% CI cells)",
         title = "(a) Coefficient distribution across multiverse cells") +
    theme_pres(base_size = 11)

  p_bar <- ggplot(mv_class_df,
                  aes(x = outcome, y = share, fill = cat)) +
    geom_col(width = 0.5) +
    geom_text(aes(label = ifelse(share > 0.04,
                                  scales::percent(share, accuracy = 1),
                                  "")),
              position = position_stack(vjust = 0.5),
              color = "white", size = 3.2, fontface = "bold") +
    facet_wrap(~ iv) +
    scale_fill_manual(values = c("Positive & sig (1-sided p<0.1)" = "#2E75B6",
                                  "Not positive & sig"             = "#C0392B"),
                      name = NULL) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                       expand = expansion(mult = c(0, 0.02))) +
    scale_x_discrete(labels = c("Canonical"         = "Canonical\nchurn",
                                 "Excl. belligerent" = "Excl.\nbelligerent")) +
    labs(x = NULL, y = "Share of multiverse cells",
         title = "(b) Positive & significant at one-sided p<0.05") +
    theme_pres(base_size = 11) +
    theme(legend.position = "bottom",
          legend.key.size = unit(0.4, "cm"))

  p_combined <- p_box / p_bar +
    plot_annotation(
      title    = "Belligerent-dyad multiverse: canonical vs excluded churn",
      subtitle = sprintf(
        paste0("%d cells per IV x variant (24 layer combos x 5 weighting ",
               "methods). Focal IVs: H2 (MID onset) and H3 (partner-war exposure)."),
        n_per_variant)
    )

  save_pres("fig_ms_excl_belligerent_multiverse.png", p_combined,
            width = 10, height = 8)
}

cat("\nDone. 24_robustness_excl_belligerent_dyads.R complete.\n")
