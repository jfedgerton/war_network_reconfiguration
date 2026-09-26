###############################################################################
# 21_contamination_falsification.R
#
# TWO FALSIFICATION TESTS FOR THE CONTAMINATION RESULT
# ----------------------------------------------------
# 03_contamination_demo.R shows that conditioning on conflict history moves the
# tie-conflict coefficient. This script runs the two falsifications that do NOT
# appear in the main text. (The non-conflict PLACEBO test -- capability/regime
# shifts -- was moved to the main-text benchmark ladder in 03 and REMOVED here
# [2026-07-15].) Both remaining tests are run for the dyadic and monadic
# processes SEPARATELY, mirroring the 03 ladder.
#
#   (B) PERMUTED CONFLICT HISTORY
#       Objection: "Adding regressors mechanically jostles coefficients."
#       Test: randomly permute conflict history ACROSS dyads within year,
#             preserving the year marginals but breaking the dyad-specific
#             link. Re-run the conditioned model N_PERM times, once for the
#             dyadic process and once for the monadic process. This yields a
#             NULL DISTRIBUTION for each attenuation statistic, so the observed
#             attenuation becomes a p-value rather than a number.
#       PASS = observed attenuation lies outside the permutation null.
#
#   (C) FUTURE-CONFLICT CONDITIONING
#       Objection: "This is rivalry, not war -> ties. Conflict-prone dyads
#                   simply have different ties."
#       Test: condition on conflict at t+2..t+6 INSTEAD of t-2..t-6. Future
#             conflict cannot have caused the tie at t-1. If it attenuates the
#             coefficient AS MUCH as past conflict does, then what is being
#             absorbed is a time-invariant dyadic propensity -- which the dyad
#             fixed effects were supposed to absorb already -- and the temporal
#             ordering is doing no work.
#       PASS = past-conflict attenuation is MATERIALLY LARGER than
#              future-conflict attenuation.
#
#       NOTE: (C) is a DIAGNOSTIC, not a control. Future conflict is
#       post-treatment (the tie helps cause it). It is not a specification I
#       would ever report as an estimate -- only as a comparison.
#
#       THIS IS THE TEST THAT CAN HURT. If past and future attenuate equally,
#       the temporal-ordering defence in the Discussion is weaker than the
#       paper claims, and the footnote needs rewriting. Run it before
#       finalising the text, not after.
#
# Main-text specification only -- unweighted, raw, alliance / IGO / multiplex.
# No multiverse indexing here; that lives in 20_contamination_multiverse.R.
#
# Outputs:
#   data/processed/contamination_falsification.csv
#   data/processed/contamination_permutation_null.csv
#   code/03_figures/output/manuscript/fig_ms_contamination_falsification.png
###############################################################################

library(tidyverse)
library(fixest)
library(peacesciencer)
select <- dplyr::select
set.seed(123)

OUTPUT_DIR <- here::here("data", "processed")
FIG_DIR    <- here::here("code", "03_figures", "output", "manuscript")
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

YR_MIN <- 1816L
YR_MAX <- 2014L
N_PERM <- 200L      # permutation reps for test (B). Now 2 fits/layer (dyadic +
                    # monadic) per rep, so ~2x the earlier runtime (~80 min).

# ===========================================================================
# 1. Dyad-year panel  (identical construction to 03)
# ===========================================================================
cat("Building dyad-year panel...\n")

d <- create_dyadyears(subset_years = YR_MIN:YR_MAX, directed = FALSE) |>
  add_cow_mids()   |>
  add_democracy()  |>
  add_nmc()        |>
  add_cow_majors() |>
  add_peace_years()

py_col <- intersect(c("cowmidspell", "gmlmidspell", "midspell"), names(d))[1]
d$peace_years <- d[[py_col]]

d <- d |>
  mutate(
    dyad_id         = paste0(pmin(ccode1, ccode2), "_", pmax(ccode1, ccode2)),
    # Dyadic capability asymmetry, bounded [0,1]; matches 03/19/21.
    cap_ratio       = abs(cinc1 - cinc2) / (cinc1 + cinc2),
    joint_democracy = as.integer(pmin(polity21, polity22, na.rm = TRUE) >= 6),
    major_power     = as.integer(cowmaj1 == 1 | cowmaj2 == 1),
    min_polity      = pmin(polity21, polity22, na.rm = TRUE)
  )

CTRL <- "cap_ratio + joint_democracy + major_power +
         peace_years + I(peace_years^2) + I(peace_years^3)"

# ===========================================================================
# 2. Tie measures (main-text: unweighted, raw)
# ===========================================================================
cat("Building alliance / IGO / multiplex ties...\n")

dyad_data <- readRDS(file.path(OUTPUT_DIR, "dyad_year_network_layers.rds"))

norm_year <- function(x, yr) {
  tibble(x = x, yr = yr) |>
    group_by(yr) |>
    mutate(out = if (max(x, na.rm = TRUE) > 0) x / max(x, na.rm = TRUE) else 0) |>
    ungroup() |> pull(out)
}

L <- dyad_data |>
  mutate(across(c(layer_offensive, layer_defensive, layer_igo),
                ~ replace_na(.x, 0))) |>
  mutate(
    tie_alliance  = norm_year(layer_offensive + layer_defensive, year),
    tie_igo       = norm_year(layer_igo, year),
    tie_multiplex = norm_year(layer_offensive + layer_defensive + layer_igo, year)
  )

ties <- L |>
  transmute(ccode1 = pmin(ccode1, ccode2), ccode2 = pmax(ccode1, ccode2), year,
            tie_alliance, tie_igo, tie_multiplex) |>
  group_by(ccode1, ccode2, year) |>
  summarise(across(starts_with("tie_"), ~ max(.x, na.rm = TRUE)),
            .groups = "drop")

d <- d |>
  left_join(ties, by = c("ccode1", "ccode2", "year")) |>
  mutate(across(starts_with("tie_"), ~ replace_na(.x, 0))) |>
  arrange(dyad_id, year) |>
  group_by(dyad_id) |>
  mutate(across(starts_with("tie_"), ~ dplyr::lag(.x, 1), .names = "{.col}_lag")) |>
  ungroup()

# ===========================================================================
# 3. Conflict history  (t-2 .. t-6, the main-text window)
# ===========================================================================
cat("Building conflict history (t-2..t-6) and its counterfactuals...\n")

# Matrix-free windowed sums. rowSums(sapply(...)) errors on single-row groups
# (a state/dyad present for one year), where sapply returns a vector not a
# matrix. Reduce over a list is size-safe.
lag_sum  <- function(x, ks) Reduce(`+`, lapply(ks, function(k) dplyr::lag(x, k, default = 0)))
lead_sum <- function(x, ks) Reduce(`+`, lapply(ks, function(k) dplyr::lead(x, k, default = 0)))

# --- dyadic: did i and j fight EACH OTHER?
d <- d |>
  arrange(dyad_id, year) |>
  group_by(dyad_id) |>
  mutate(
    dyad_hist_past   = lag_sum(cowmidonset, 2:6),    # (main)
    dyad_hist_future = lead_sum(cowmidonset, 2:6)    # (C) diagnostic
  ) |>
  ungroup()

# --- monadic: did i or j fight ANYONE?
state_churn <- readRDS(file.path(OUTPUT_DIR, "state_network_churn.rds"))
state_hist <- state_churn |>
  select(ccode, year, mid_onset_count) |>
  arrange(ccode, year) |>
  group_by(ccode) |>
  mutate(mon_past   = lag_sum(mid_onset_count, 2:6),
         mon_future = lead_sum(mid_onset_count, 2:6)) |>
  ungroup() |>
  select(ccode, year, mon_past, mon_future)

d <- d |>
  left_join(state_hist |> rename(mp_i = mon_past, mf_i = mon_future),
            by = c("ccode1" = "ccode", "year")) |>
  left_join(state_hist |> rename(mp_j = mon_past, mf_j = mon_future),
            by = c("ccode2" = "ccode", "year")) |>
  mutate(across(c(mp_i, mp_j, mf_i, mf_j), ~ replace_na(.x, 0)),
         # THIRD-PARTY monadic (mirrors 03's mid5_sum_tp): each state's own
         # onsets over the window MINUS its i-j onsets, so the monadic term is
         # conflict with THIRD parties only. dyad_hist_past/future count the i-j
         # onsets over the matching window. Clamp at 0.
         mon_hist_past   = pmax(mp_i - dyad_hist_past,   0) + pmax(mp_j - dyad_hist_past,   0),
         mon_hist_future = pmax(mf_i - dyad_hist_future, 0) + pmax(mf_j - dyad_hist_future, 0))

# ===========================================================================
# 4. Estimation sample (main-text: at-risk)
# ===========================================================================
# NOTE: the non-conflict PLACEBO histories (capability/regime shifts) were
# removed here [2026-07-15]. That falsification now lives in the main-text
# benchmark ladder (script 03). This script retains the two falsifications that
# do NOT appear in the main text: the permutation null (B) and the future-
# conflict diagnostic (C).
need <- c("cowmidonset", "tie_alliance_lag", "tie_igo_lag", "tie_multiplex_lag",
          "cap_ratio", "joint_democracy", "major_power", "peace_years",
          "dyad_hist_past", "mon_hist_past",
          "dyad_hist_future", "mon_hist_future")

est <- d |> filter(if_all(all_of(need), ~ !is.na(.x)))
n_pre <- nrow(est)
est <- est |> filter(coalesce(cowmidongoing, 0L) == 0 | cowmidonset == 1)

sd_safe <- function(x) { s <- sd(x, na.rm = TRUE); if (is.na(s) || s == 0) 1 else s }
est <- est |>
  mutate(
    across(c(tie_alliance_lag, tie_igo_lag, tie_multiplex_lag),
           ~ .x / sd_safe(.x), .names = "{.col}_sd"),
    # Conflict-history onset counts: log1p (mirrors the main analyses) then
    # per-SD, for the dyadic and monadic terms, past and future.
    across(c(dyad_hist_past, mon_hist_past, dyad_hist_future, mon_hist_future),
           ~ { z <- log1p(.x); z / sd_safe(z) }, .names = "{.col}_sd"))

cat(sprintf("  at-risk sample: %s -> %s dyad-years | onset rate %.3f%%\n",
            format(n_pre, big.mark = ","), format(nrow(est), big.mark = ","),
            100 * mean(est$cowmidonset)))

LAYERS <- c("Alliance"  = "tie_alliance_lag_sd",
            "IGO"       = "tie_igo_lag_sd",
            "Multiplex" = "tie_multiplex_lag_sd")

# Helper: fit the ladder and return the tie coefficient.
fit_tie <- function(iv, add = "") {
  f <- as.formula(paste("cowmidonset ~", iv, "+", CTRL, add, "| dyad_id + year"))
  m <- tryCatch(feols(f, data = est, vcov = ~dyad_id), error = function(e) NULL)
  if (is.null(m) || !iv %in% rownames(coeftable(m))) return(c(NA, NA, NA))
  ct <- coeftable(m)
  c(ct[iv, "Estimate"], ct[iv, "Std. Error"], ct[iv, "Pr(>|t|)"])
}

# ===========================================================================
# 6. TEST C -- FUTURE-CONFLICT DIAGNOSTIC (dyadic and monadic, separately)
# ===========================================================================
cat("\n=====================================================================\n")
cat("TEST C: future-conflict diagnostic (dyadic and monadic processes)\n")
cat("=====================================================================\n")

# Dyadic and monadic are estimated SEPARATELY, mirroring the 03 ladder. Each
# past process has a matching future counterfactual (t+2..t+6).
CONDITION_SETS <- list(
  "Baseline (no history)"                         = "",
  "Dyadic conflict, past t-2..t-6"                = "+ dyad_hist_past_sd",
  "Monadic conflict, past t-2..t-6"               = "+ mon_hist_past_sd",
  "DIAGNOSTIC: dyadic conflict, future t+2..t+6"  = "+ dyad_hist_future_sd",
  "DIAGNOSTIC: monadic conflict, future t+2..t+6" = "+ mon_hist_future_sd"
)

fals <- purrr::imap_dfr(LAYERS, function(iv, lyr) {
  purrr::imap_dfr(CONDITION_SETS, function(add, nm) {
    r <- fit_tie(iv, add)
    tibble(layer = lyr, condition = nm,
           beta = r[1], se = r[2], p = r[3])
  })
}) |>
  group_by(layer) |>
  mutate(beta_base       = beta[condition == "Baseline (no history)"],
         attenuation_pct = 100 * (1 - abs(beta) / abs(beta_base))) |>
  ungroup()

print(fals |> select(layer, condition, beta, se, p, attenuation_pct),
      digits = 3, n = Inf)

write_csv(fals, file.path(OUTPUT_DIR, "contamination_falsification.csv"))

cat("\n--- VERDICT (C: does FUTURE conflict absorb as much as PAST?) ---\n")
for (l in names(LAYERS)) {
  s   <- fals |> filter(layer == l)
  get <- function(cond) s$attenuation_pct[s$condition == cond]
  d_past <- get("Dyadic conflict, past t-2..t-6")
  d_fut  <- get("DIAGNOSTIC: dyadic conflict, future t+2..t+6")
  m_past <- get("Monadic conflict, past t-2..t-6")
  m_fut  <- get("DIAGNOSTIC: monadic conflict, future t+2..t+6")
  cat(sprintf("\n  %s\n", l))
  cat(sprintf("    dyadic  past %+6.1f%% | future %+6.1f%%  -> C %s\n", d_past, d_fut,
              if (abs(d_fut) < 0.7 * abs(d_past)) "PASS"
              else "FAIL (future absorbs as much -- rivalry story, not war->ties)"))
  cat(sprintf("    monadic past %+6.1f%% | future %+6.1f%%  -> C %s\n", m_past, m_fut,
              if (abs(m_fut) < 0.7 * abs(m_past)) "PASS"
              else "FAIL (future absorbs as much -- rivalry story, not war->ties)"))
}

# ===========================================================================
# 7. TEST B -- PERMUTED CONFLICT HISTORY
# ===========================================================================
# Permute the conflict-history columns ACROSS DYADS WITHIN YEAR. This preserves
# the year marginals (how much conflict history exists in each year) but breaks
# the link between a dyad's OWN history and its OWN tie. Under the null that the
# attenuation is a mechanical artifact of adding regressors, permuted history
# should attenuate just as much.
#
# Only rung (c) is refit -- rung (a) contains no history and does not change.
# ===========================================================================
cat("\n=====================================================================\n")
cat(sprintf("TEST B: permutation null (%d replications)\n", N_PERM))
cat("=====================================================================\n")

# The dyadic and monadic processes each get their OWN permutation null, against
# their OWN observed (single-process) attenuation -- mirroring the 03 ladder.
obs_atten <- fals |>
  filter(condition %in% c("Dyadic conflict, past t-2..t-6",
                          "Monadic conflict, past t-2..t-6")) |>
  transmute(layer,
            measure   = if_else(grepl("Dyadic", condition), "Dyadic", "Monadic"),
            obs_atten = attenuation_pct)

PERM_MEASURES <- c(Dyadic = "+ dyad_hist_perm", Monadic = "+ mon_hist_perm")

perm_rows <- list()
t0 <- proc.time()

for (b in seq_len(N_PERM)) {
  gc()
  est_p <- est |>
    group_by(year) |>
    mutate(dyad_hist_perm = sample(dyad_hist_past_sd),
           mon_hist_perm  = sample(mon_hist_past_sd)) |>
    ungroup()

  for (l in names(LAYERS)) {
    iv    <- LAYERS[[l]]
    bbase <- fals$beta[fals$layer == l & fals$condition == "Baseline (no history)"]
    for (mz in names(PERM_MEASURES)) {
      f <- as.formula(paste("cowmidonset ~", iv, "+", CTRL,
                            PERM_MEASURES[[mz]], "| dyad_id + year"))
      m <- tryCatch(feols(f, data = est_p, vcov = ~dyad_id), error = function(e) NULL)
      if (is.null(m) || !iv %in% rownames(coeftable(m))) next
      bcond <- coeftable(m)[iv, "Estimate"]
      perm_rows[[length(perm_rows) + 1]] <- tibble(
        rep = b, layer = l, measure = mz,
        attenuation_pct = 100 * (1 - abs(bcond) / abs(bbase))
      )
    }
  }
  if (b %% 25 == 0) {
    cat(sprintf("  %d / %d  (%.1f min)\n", b, N_PERM,
                (proc.time() - t0)["elapsed"] / 60))
  }
}

perm <- bind_rows(perm_rows)
write_csv(perm, file.path(OUTPUT_DIR, "contamination_permutation_null.csv"))

perm_summary <- perm |>
  left_join(obs_atten, by = c("layer", "measure")) |>
  group_by(layer, measure) |>
  summarise(
    obs_atten   = first(obs_atten),
    null_median = median(attenuation_pct),
    null_q025   = quantile(attenuation_pct, 0.025),
    null_q975   = quantile(attenuation_pct, 0.975),
    # one-sided p: share of permutations at least as extreme, in the OBSERVED
    # direction (attenuation for alliance/multiplex, amplification for IGO)
    p_perm      = if (first(obs_atten) > 0)
                    mean(attenuation_pct >= first(obs_atten))
                  else
                    mean(attenuation_pct <= first(obs_atten)),
    .groups = "drop"
  )

cat("\n--- PERMUTATION NULL vs OBSERVED ---\n")
print(perm_summary, digits = 3, n = Inf)
cat("\n  p_perm = share of permutations that reproduced the observed movement.\n")
cat("  Small p => the attenuation is NOT a mechanical artifact of adding regressors.\n")

write_csv(perm_summary, file.path(OUTPUT_DIR, "contamination_permutation_summary.csv"))

# ===========================================================================
# 8. FIGURE
# ===========================================================================
cat("\n--- Plotting ---\n")

ACCENT <- "#2C6CB0"; ACCENT2 <- "#B03366"

pa <- fals |>
  filter(condition != "Baseline (no history)") |>
  mutate(
    # kind must be computed from the ORIGINAL labels (before the line break).
    kind = if_else(grepl("DIAGNOSTIC", condition),
                   "Future conflict (diagnostic)", "Conflict history (past)"),
    # Break the y-axis label after the colon so long labels render cleanly.
    condition = factor(gsub(":\\s*|,\\s*", "\n", condition),
                       levels = rev(gsub(":\\s*|,\\s*", "\n", names(CONDITION_SETS)[-1])))
  ) |>
  ggplot(aes(x = attenuation_pct, y = condition, fill = kind)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray55") +
  geom_col(width = 0.6) +
  facet_wrap(~ layer, ncol = 3, scales = "free_x") +
  scale_fill_manual(values = c("Conflict history (past)"        = ACCENT,
                               "Future conflict (diagnostic)"   = ACCENT2)) +
  labs(title = "Falsification: PAST conflict moves the estimate, FUTURE conflict does not",
       subtitle = paste0("Attenuation of the tie-conflict coefficient, dyadic and monadic ",
                         "processes separately.\nPositive = toward zero, negative = away from ",
                         "zero. Future conflict cannot have shaped a past tie."),
       x = "Attenuation (%)", y = NULL, fill = NULL) +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(),
        strip.background = element_rect(fill = "gray92", color = "gray70"),
        strip.text = element_text(face = "bold"),
        legend.position = "bottom",
        axis.text.x = element_text(angle = 0, hjust = 0.5),
        plot.title = element_text(face = "bold"),
        plot.subtitle = element_text(color = "gray30", size = 9))

pb <- perm |>
  left_join(obs_atten, by = c("layer", "measure")) |>
  ggplot(aes(x = attenuation_pct)) +
  geom_histogram(aes(fill = "Permuted draws (null)"), bins = 40, color = "white") +
  geom_vline(aes(xintercept = obs_atten, color = "Observed"), linewidth = 1) +
  facet_grid(measure ~ layer, scales = "free") +
  scale_fill_manual(name = NULL, values = c("Permuted draws (null)" = "gray75")) +
  scale_color_manual(name = NULL, values = c("Observed" = ACCENT2)) +
  labs(title = sprintf("Permutation null (%d draws): observed movement is not mechanical",
                       N_PERM),
       subtitle = paste0("Attenuation when conflict history is randomly reassigned across ",
                         "dyads within year,\ndyadic and monadic processes separately."),
       x = "Attenuation (%)", y = "Permutations") +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "gray92", color = "gray70"),
        strip.text = element_text(face = "bold"),
        axis.text.x = element_text(angle = 0, hjust = 0.5),
        legend.position = "bottom",
        plot.title = element_text(face = "bold"),
        plot.subtitle = element_text(color = "gray30", size = 9))

library(patchwork)
p_full <- pa / pb + plot_layout(heights = c(1.1, 1))

ggsave(file.path(FIG_DIR, "fig_ms_contamination_falsification.png"), p_full,
       width = 11, height = 9, dpi = 400, bg = "white")
cat("  -> fig_ms_contamination_falsification.png\n")

cat("\nSaved -> contamination_falsification.csv\n")
cat("Saved -> contamination_permutation_null.csv\n")
cat("Saved -> contamination_permutation_summary.csv\n")
