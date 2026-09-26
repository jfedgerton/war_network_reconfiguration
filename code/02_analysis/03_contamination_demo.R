###############################################################################
# 28_contamination_demo.R
#
# THE CONTAMINATION RESULT
# ------------------------
# CLAIM (explicitly NOT causal): the network->conflict association is not
# identified independently of war history. Conditioning on the conflict history
# that PRODUCES network ties attenuates the estimated pacifying association.
#
# The argument is that conflict alters the MULTIPLEX structure. So the ladder is
# run separately on each layer:
#     (a) ALLIANCE ties   -- the alliance-and-conflict literature
#     (b) IGO ties        -- the institutionalist / IGO literature
#     (c) MULTIPLEX ties  -- the composite this paper uses
# If conflict contaminates the multiplex network, it must contaminate its
# components. Each layer gets its own baseline, its own guard, and its own
# attenuation.
#
# THE DAG
#   W(t-2..t-6) --> N(t-1) --> W(t)     war shapes the tie (this paper)
#   W(t-2..t-6) ------------> W(t)      recurrence
#   U ----------> N, W(t), W(t-1)       unobserved time-invariant rivalry
#
#   => W is a CONFOUNDER of N -> W. The backdoor N <- W <- ... -> W(t) is left
#      OPEN by the literature. Conditioning on it closes it, and the MOVEMENT IN
#      THE TIE COEFFICIENT IS THE FINDING.
#   => Conditioning on W would open a collider path through U. DYAD FIXED EFFECTS
#      absorb U -- which is why this, not the pooled model, is the baseline.
#
# FOUR DESIGN CHOICES THAT MAKE IT CREDIBLE
#   (1) ACYCLIC. War history is measured t-2..t-6, strictly BEFORE the tie at
#       t-1. No simultaneity between war and the tie. Conservative => the
#       attenuation is a LOWER BOUND.
#   (2) NON-MONOTONE. Cumulative counts rise monotonically within a dyad and
#       would act as a dyad-specific TREND inside a FE model. 5-year windows only.
#   (3) AT-RISK SAMPLE. A dyad already in an ongoing MID cannot have a new onset.
#   (4) NO year-t leakage. Every history regressor is built from lag(onset, >=2).
#
# TWO THINGS DELIBERATELY EXCLUDED
#   * PARTNER-WAR EXPOSURE. PWE = sum_k w_ik * 1[k at war] / sum_k w_ik -- its
#     denominator IS the sum of tie weights, and j is one of those k. So PWE is a
#     function of the very tie it would be controlling for. Circular. (It also
#     moved the estimate by <1pp, so nothing is lost.)
#   * PRIOR CHURN. War -> churn -> tie, so churn is a MEDIATOR, not a confounder.
#     Conditioning on it is post-treatment adjustment.
#
# Output: data/processed/contamination_demo.csv
###############################################################################

library(tidyverse)
library(fixest)
library(peacesciencer)
select <- dplyr::select
set.seed(123)

OUTPUT_DIR <- here::here("data", "processed")
YR_MIN <- 1816L
YR_MAX <- 2014L

# ===========================================================================
# 1. Dyad-year panel with the literature's own covariates
# ===========================================================================
cat("Building dyad-year panel...\n")

d <- create_dyadyears(subset_years = YR_MIN:YR_MAX, directed = FALSE) |>
  add_cow_mids()          |>   # cowmidonset (DV), cowmidongoing
  add_contiguity()        |>
  add_minimum_distance()  |>
  add_democracy()         |>
  add_nmc()               |>
  add_cow_majors()        |>
  add_peace_years()

py_col <- intersect(c("cowmidspell", "gmlmidspell", "midspell"), names(d))[1]
stopifnot(!is.na(py_col))
d$peace_years <- d[[py_col]]

# ===========================================================================
# 2. THE THREE TIE MEASURES (the layers of the multiplex network)
# ===========================================================================
cat("Building alliance / IGO / multiplex tie measures...\n")

dyad_data <- readRDS(file.path(OUTPUT_DIR, "dyad_year_network_layers.rds"))

# Normalize each layer set to [0,1] by the within-year maximum, exactly as the
# main paper's multiplex weight is constructed.
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
    raw_alliance  = layer_offensive + layer_defensive,
    raw_igo       = layer_igo,
    raw_multiplex = layer_offensive + layer_defensive + layer_igo
  ) |>
  mutate(
    tie_alliance  = norm_year(raw_alliance,  year),
    tie_igo       = norm_year(raw_igo,       year),
    tie_multiplex = norm_year(raw_multiplex, year)
  )

# Canonicalize to undirected (ccode1 < ccode2) to match create_dyadyears()
ties <- L |>
  transmute(a = pmin(ccode1, ccode2), b = pmax(ccode1, ccode2), year,
            tie_alliance, tie_igo, tie_multiplex) |>
  group_by(a, b, year) |>
  summarise(across(starts_with("tie_"), ~ max(.x, na.rm = TRUE)), .groups = "drop") |>
  rename(ccode1 = a, ccode2 = b)

d <- d |>
  mutate(dyad_id = paste0(pmin(ccode1, ccode2), "_", pmax(ccode1, ccode2))) |>
  left_join(ties, by = c("ccode1", "ccode2", "year")) |>
  # absent from the network file = NO tie (not missing). Leaving these NA would
  # silently restrict the sample to tie-having dyads.
  mutate(across(starts_with("tie_"), ~ replace_na(.x, 0))) |>
  arrange(dyad_id, year) |>
  group_by(dyad_id) |>
  mutate(across(starts_with("tie_"), ~ lag(.x, 1), .names = "{.col}_lag")) |>
  ungroup()

cat(sprintf("  tie > 0:  alliance %.1f%% | IGO %.1f%% | multiplex %.1f%%\n",
            100 * mean(d$tie_alliance > 0),
            100 * mean(d$tie_igo > 0),
            100 * mean(d$tie_multiplex > 0)))

# ===========================================================================
# 3. WAR HISTORY -- non-monotone, and STRICTLY PRIOR TO THE TIE (t-2 .. t-6)
#
#    W(t-2..t-6) --> N(t-1) --> W(t).  War precedes the tie; the tie precedes
#    the onset. Acyclic. Conservative => attenuation is a LOWER BOUND.
# ===========================================================================
cat("Building war-history regressors (t-2 .. t-6)...\n")

# (a) DYADIC: i vs j
d <- d |>
  arrange(dyad_id, year) |>
  group_by(dyad_id) |>
  mutate(dyad_mid_5yr_lag =
           lag(cowmidonset, 2, default = 0) + lag(cowmidonset, 3, default = 0) +
           lag(cowmidonset, 4, default = 0) + lag(cowmidonset, 5, default = 0) +
           lag(cowmidonset, 6, default = 0)) |>
  ungroup()

# (b) MONADIC: each state's own wars
state_churn <- readRDS(file.path(OUTPUT_DIR, "state_network_churn.rds"))

state_hist <- state_churn |>
  select(ccode, year, mid_onset_count) |>
  arrange(ccode, year) |>
  group_by(ccode) |>
  mutate(mid_5yr_lag =
           lag(mid_onset_count, 2, default = 0) + lag(mid_onset_count, 3, default = 0) +
           lag(mid_onset_count, 4, default = 0) + lag(mid_onset_count, 5, default = 0) +
           lag(mid_onset_count, 6, default = 0)) |>
  ungroup() |>
  select(ccode, year, mid_5yr_lag)

# Symmetrize: in an UNDIRECTED dyad, i and j are arbitrary (ccode1 < ccode2), so
# entering _i and _j separately would make coefficients depend on the sort order.
#
# Idempotent: drop any mid5_* left over from a previous run of this block. Without
# this, re-running turns mid5_i into mid5_i.x / mid5_i.y and the mutate below fails.
d <- d |> select(-any_of(c("mid5_i", "mid5_j", "mid5_sum", "mid5_max", "mid5_min")))

stopifnot(
  "state_hist has duplicate ccode-year rows" =
    !anyDuplicated(state_hist[, c("ccode", "year")])
)

d <- d |>
  left_join(state_hist |> rename(mid5_i = mid_5yr_lag), by = c("ccode1" = "ccode", "year")) |>
  left_join(state_hist |> rename(mid5_j = mid_5yr_lag), by = c("ccode2" = "ccode", "year"))

cat("monadic war history unmatched: ",
    sum(is.na(d$mid5_i) | is.na(d$mid5_j)), " of ", nrow(d), " dyad-years\n", sep = "")

# Unmatched = state-year not present in state_churn. Treat as zero prior wars.
# NOTE ON PARAMETERIZATION. The dyad is UNDIRECTED (ccode1 < ccode2 is a sorting
# rule, not a role), so only symmetric functions of (mid5_i, mid5_j) are
# admissible. sum, max, and min all qualify; mid5_i alone does not.
#
# But sum = min + max identically. So entering (sum, max) together is an exact
# reparameterization of (min, max), and the "sum" coefficient is then really the
# coefficient on MIN -- the quieter state's record. That is not the quantity we
# want to label "monadic conflict history."
#
# MAIN TEXT  -> mid5_sum ONLY. Interpretable: total recent conflict in the dyad.
# APPENDIX   -> max-only and (min + max). See section 9.
d <- d |>
  mutate(mid5_i = coalesce(mid5_i, 0),
         mid5_j = coalesce(mid5_j, 0),
         mid5_sum = mid5_i + mid5_j,
         mid5_max = pmax(mid5_i, mid5_j),
         mid5_min = pmin(mid5_i, mid5_j),
         # THIRD-PARTY monadic history: each state's onsets with ANYONE minus its
         # onsets against the dyad partner, so the monadic term captures conflict
         # with THIRD parties only. This cleanly separates the two processes:
         #   dyadic  = i vs j;  monadic = i or j vs others (no i-j double count).
         # dyad_mid_5yr_lag counts the i-j onsets over t-2..t-6. Clamp at 0 in
         # case the monadic and dyadic counts differ marginally.
         mid5_sum_tp = pmax(mid5_i - dyad_mid_5yr_lag, 0) +
                       pmax(mid5_j - dyad_mid_5yr_lag, 0))

# ===========================================================================
# 3b. BENCHMARK HISTORIES -- mirror the dyadic/monadic conflict-history structure
#     with NON-conflict shocks over the same t-2..t-6 window.
#
#     CAVEAT (stated in a paper footnote): conflict itself moves capabilities
#     and regimes, so these are conflict-ADJACENT, not clean benchmarks. They are
#     a MECHANISM reference -- if power/regime shifts do NOT reproduce the
#     attenuation, the contamination runs through conflict specifically (which
#     also aligns with the mediation finding that power change does not mediate).
#
#       Capability, between-unit (dyadic) : |change in capability asymmetry|
#       Capability, within-unit (monadic) : |change in own CINC|, summed
#       Regime,     between-unit (dyadic) : |change in polity gap|
#       Regime,     within-unit (monadic) : |change in own polity|, summed
# ===========================================================================
cat("Building benchmark histories (capability + regime shifts, t-2 .. t-6)...\n")

lag5 <- function(x) {
  x <- tidyr::replace_na(x, 0)
  Reduce(`+`, lapply(2:6, function(k) dplyr::lag(x, k, default = 0)))
}

# -- dyadic (between-unit) benchmarks, on the dyad panel --
d <- d |>
  arrange(dyad_id, year) |>
  group_by(dyad_id) |>
  mutate(
    .cap_asym      = abs(cinc1 - cinc2) / (cinc1 + cinc2),
    cap_between    = lag5(abs(.cap_asym - dplyr::lag(.cap_asym, 1))),
    .polity_gap    = abs(polity21 - polity22),
    regime_between = lag5(abs(.polity_gap - dplyr::lag(.polity_gap, 1)))
  ) |>
  ungroup() |>
  select(-.cap_asym, -.polity_gap)

# -- within-unit (monadic) benchmarks: built per state, summed across the dyad --
state_cov <- dplyr::bind_rows(
  d |> select(ccode = ccode1, year, cinc = cinc1, polity = polity21),
  d |> select(ccode = ccode2, year, cinc = cinc2, polity = polity22)
) |>
  group_by(ccode, year) |>
  summarise(cinc = first(cinc), polity = first(polity), .groups = "drop") |>
  arrange(ccode, year) |>
  group_by(ccode) |>
  mutate(cinc_shift   = lag5(abs(cinc   - dplyr::lag(cinc, 1))),
         polity_shift = lag5(abs(polity - dplyr::lag(polity, 1)))) |>
  ungroup() |>
  select(ccode, year, cinc_shift, polity_shift)

d <- d |> select(-any_of(c("cs_i","cs_j","ps_i","ps_j","cap_within","regime_within")))
d <- d |>
  left_join(state_cov |> rename(cs_i = cinc_shift, ps_i = polity_shift),
            by = c("ccode1" = "ccode", "year")) |>
  left_join(state_cov |> rename(cs_j = cinc_shift, ps_j = polity_shift),
            by = c("ccode2" = "ccode", "year")) |>
  mutate(across(c(cs_i, cs_j, ps_i, ps_j), ~ coalesce(.x, 0)),
         cap_within    = cs_i + cs_j,
         regime_within = ps_i + ps_j)

# ===========================================================================
# 4. Controls + estimation sample
# ===========================================================================
d <- d |>
  mutate(
    contig          = as.integer(!is.na(conttype) & conttype > 0),
    log_distance    = log1p(replace_na(mindist, max(mindist, na.rm = TRUE))),
    # Dyadic capability asymmetry |CINC_i - CINC_j| / (CINC_i + CINC_j), bounded
    # [0,1]: 0 = parity, ->1 = one state dwarfs the other. Standard, well-behaved
    # dyadic power control; order-invariant for undirected dyads. IDENTICAL
    # construction to 19/20/21 so the baseline matches the benchmark models.
    cap_ratio       = abs(cinc1 - cinc2) / (cinc1 + cinc2),
    joint_democracy = as.integer(pmin(polity21, polity22, na.rm = TRUE) >= 6),
    major_power     = as.integer(cowmaj1 == 1 | cowmaj2 == 1)
  )

# NO separate alliance control -- alliances are INSIDE the tie measure.
# Temporal dependence: peace years + cubic polynomial (Carter & Signorino 2010).
#
# GEOGRAPHY IS DROPPED (contiguity and minimum distance), for two reasons:
#   (1) Both are essentially time-invariant within a dyad, so the DYAD FIXED
#       EFFECTS absorb them. They contribute nothing to identification.
#   (2) More importantly, whatever within-dyad variation they DO have is driven
#       by TERRITORIAL CHANGE -- which is war-driven. That makes contiguity
#       POST-TREATMENT relative to the war history being conditioned on:
#         W(t-2) -> territory -> contiguity -> Y(t)
#       Conditioning on it would absorb part of the war effect we are measuring
#       and could induce collider bias. Same objection as prior churn.
# Referee answer: "geography is absorbed by the dyad fixed effects."
#
# NOTE: `peace_years` IS ITSELF A CONFLICT-HISTORY CONTROL, by construction.
# That is a feature: the baseline already contains the literature's own way of
# handling conflict history, and the tie coefficient STILL attenuates when the
# war-history terms are added. Temporal dependence in the OUTCOME cannot fix an
# endogenous REGRESSOR.
CTRL <- "cap_ratio + joint_democracy + major_power +
         peace_years + I(peace_years^2) + I(peace_years^3)"

need <- c("cowmidonset","tie_alliance_lag","tie_igo_lag","tie_multiplex_lag",
          "contig","log_distance","cap_ratio","joint_democracy","major_power",
          "peace_years","dyad_mid_5yr_lag","mid5_sum","mid5_max","mid5_min")
est <- d |> filter(if_all(all_of(need), ~ !is.na(.x)))
n_pre <- nrow(est)

# AT-RISK SAMPLE: a dyad already in an ongoing MID cannot have a NEW onset.
est <- est |> filter(coalesce(cowmidongoing, 0L) == 0 | cowmidonset == 1)

# Standardize each tie AFTER the restriction, so coefficients are per-SD and
# comparable ACROSS layers (alliance and IGO ties have very different scales).
# sd_safe matches 20's convention (returns 1 for a degenerate column) so the
# baseline here is identical to the benchmark models.
sd_safe <- function(x) { s <- sd(x, na.rm = TRUE); if (is.na(s) || s == 0) 1 else s }
est <- est |>
  mutate(across(c(tie_alliance_lag, tie_igo_lag, tie_multiplex_lag),
                ~ .x / sd_safe(.x), .names = "{.col}_sd"))

# STANDARDIZE THE CONFLICT-HISTORY TERMS TOO -- required for Step 1 to be
# readable. The two terms live on very different scales:
#   dyad_mid_5yr_lag : did THESE TWO fight each other, t-2..t-6   (range ~0-3)
#   mid5_sum         : disputes either fought with ANYONE, summed (range ~0-20+)
# Comparing raw slopes across them is meaningless: a small coefficient on a
# wide-ranging variable can move the tie as much as a large one on a narrow
# variable. Per-SD units make the two rows of Step 1 directly comparable, so a
# monadic null is a real null and not an artifact of the units.
#
# NOTE: this CANNOT change Step 2. Rescaling a control leaves the coefficient on
# the tie untouched. Attenuation percentages are identical either way.
est <- est |>
  mutate(
    # LADDER conflict-history terms: log1p FIRST -- mirrors the main analyses,
    # which enter conflict as log1p(onset count) -- THEN per-SD. Applies to both
    # the dyadic (i-j) and monadic (third-party) onset counts.
    across(c(dyad_mid_5yr_lag, mid5_sum_tp),
           ~ { z <- log1p(.x); z / sd_safe(z) }, .names = "{.col}_sd"),
    # Section-9 monadic carve-up stays in RAW counts so the sum = min + max
    # identity holds. Benchmarks are shifts (not counts), so also raw.
    across(c(mid5_sum, mid5_max, mid5_min,
             cap_between, cap_within, regime_between, regime_within),
           ~ .x / sd_safe(.x), .names = "{.col}_sd"))

HIST_SDS <- est |>
  summarise(across(c(dyad_mid_5yr_lag, mid5_sum, mid5_max, mid5_min), sd))

cat(sprintf("\n  at-risk restriction: %s -> %s dyad-years | onset rate %.3f%%\n",
            format(n_pre, big.mark = ","), format(nrow(est), big.mark = ","),
            100 * mean(est$cowmidonset)))

cat("\n  SDs of the conflict-history terms (for the footnote):\n")
print(HIST_SDS, digits = 4)

# ===========================================================================
# 5. THE LADDER, RUN SEPARATELY ON EACH LAYER
#    dyad + year FE, LPM, clustered by dyad.
#    LPM (not FE-logit): FE-logit drops every dyad that never has an onset.
# ===========================================================================
cat("\nEstimating (dyad + year FE, LPM, clustered by dyad)...\n")

# NB: these are the two war-history BLOCKS in the regression. They are NOT the
# paper's hypotheses (H1 = systemic conflict -> systemic reconfiguration;
# H2 = conflict onset -> belligerent churn; H3 = partner war exposure -> churn).
# Named HIST_* to avoid that collision.
# Per-SD units (see section 4). Step 2's tie coefficients are unaffected.
HIST_DYAD <- "+ dyad_mid_5yr_lag_sd"  # did i and j fight EACH OTHER?  (t-2..t-6)
# THIRD-PARTY monadic: did i or j fight OTHER states (not each other)? Excludes
# the i-j interstate onsets so the dyadic and monadic processes do not overlap.
HIST_MON  <- "+ mid5_sum_tp_sd"       # did i or j fight THIRD parties? (t-2..t-6)

# BENCHMARK conditioning sets, mirroring the dyadic + monadic conflict structure
# with NON-conflict shocks (see section 3b).
PLAC_CAP    <- "+ cap_between_sd + cap_within_sd"      # capability shift
PLAC_REGIME <- "+ regime_between_sd + regime_within_sd" # regime change

adds <- c("", HIST_DYAD, HIST_MON,
          "+ cap_between_sd", "+ cap_within_sd",
          "+ regime_between_sd", "+ regime_within_sd")

# "Conflict", not "war": the DV is MID onset (disputes), not COW war. Keep the
# vocabulary identical in the labels, the panel titles, and the prose.
labs <- c("(a) Baseline",
          "(b) + dyadic\nconflict history",
          "(c) + monadic\nconflict history",
          "Benchmark:\ncapability between",
          "Benchmark:\ncapability within",
          "Benchmark:\nregime change between",
          "Benchmark:\nregime change within")
# Family drives the fill/legend so the benchmark reference rungs stand out.
fam  <- c("Baseline", "Conflict history", "Conflict history",
          "Benchmark", "Benchmark","Benchmark", "Benchmark")

LAYERS <- c("Alliance"  = "tie_alliance_lag_sd",
            "IGO"       = "tie_igo_lag_sd",
            "Multiplex" = "tie_multiplex_lag_sd")

run_layer <- function(lyr_name, iv) {
  purrr::pmap_dfr(list(adds, labs, fam), function(a, l, fm) {
    f <- as.formula(paste("cowmidonset ~", iv, "+", CTRL, a, "| dyad_id + year"))
    m <- feols(f, data = est, vcov = ~dyad_id)
    ct <- coeftable(m)
    tibble(layer = lyr_name, model = l, family = fm,
           beta = ct[iv, "Estimate"], se = ct[iv, "Std. Error"],
           p    = ct[iv, "Pr(>|t|)"], n = nobs(m))
  }) |>
    mutate(
      # 95% (outer, thin) and 90% (inner, thick) bands. The 90% band is what
      # lets the IGO baseline (p = .056) be read honestly: it clears 90% but not
      # 95%, and the figure shows that rather than arguing about it.
      ci95_low  = beta - 1.960 * se, ci95_high = beta + 1.960 * se,
      ci90_low  = beta - 1.645 * se, ci90_high = beta + 1.645 * se,
      # pct_retained < 100 => coefficient moved TOWARD zero (ATTENUATION)
      pct_retained = 100 * abs(beta) / abs(beta[1]),
      attenuation  = 100 - pct_retained)
}

res <- purrr::imap_dfr(LAYERS, ~ run_layer(.y, .x))

cat("\n=====================================================================\n")
cat("TIE COEFFICIENT (per SD) BY LAYER -- dyad+year FE, war history t-2..t-6\n")
cat("=====================================================================\n")
for (l in names(LAYERS)) {
  cat("\n", l, " ties\n", sep = "")
  print(res |> filter(layer == l) |>
          select(model, beta, se, p, pct_retained, attenuation),
        digits = 4, n = Inf)
}

# labs[2] = dyadic-history rung; labs[3] = monadic (third-party) rung. These are
# SEPARATE conditioning sets (each added to baseline alone), NOT cumulative.
cat("\n--- HEADLINE: attenuation from each conflict-history process (vs baseline) ---\n")
res |> filter(model %in% c(labs[2], labs[3])) |>
  mutate(txt = sprintf("  %-10s %-24s %5.1f%%   (baseline p = %.1e)",
                       layer, gsub("\n"," ",model), attenuation, p)) |>
  pull(txt) |> cat(sep = "\n")
cat("\n")
cat("--- BENCHMARK rungs (should be ~0 attenuation) ---\n")
res |> filter(family == "Benchmark") |>
  mutate(txt = sprintf("  %-10s %-28s %5.1f%%", layer, gsub("\n"," ",model), attenuation)) |>
  pull(txt) |> cat(sep = "\n")
cat("\n")

# --- DYADIC vs MONADIC: which process moves the estimate? -------------------
# `attenuation` in res is ALWAYS measured against baseline (a). The dyadic and
# monadic rungs are estimated SEPARATELY (each conditioning set added to the
# baseline on its own), so the two columns are directly comparable movements
# from the same anchor -- not a nested decomposition.
cat("\n--- DYADIC vs MONADIC: which process moves the estimate? ---\n")
decomp <- res |>
  filter(family %in% c("Baseline", "Conflict history")) |>   # conflict ladder only, not benchmarks
  select(layer, model, attenuation) |>
  tidyr::pivot_wider(names_from = model, values_from = attenuation) |>
  rename(from_dyadic = 3, from_monadic = 4) |>
  select(layer, from_dyadic, from_monadic)
print(decomp, digits = 3, n = Inf)
cat("  (positive = attenuation, negative = amplification; each vs. baseline)\n")
cat("  from_monadic uses THIRD-PARTY onsets only (i-j onsets removed)\n")
write_csv(decomp, file.path(OUTPUT_DIR, "contamination_decomposition.csv"))

# ---------------------------------------------------------------------------
# 5b. LEG 2 OF THE CONFOUNDING ARGUMENT: does prior conflict predict conflict
#     at t, NET of the tie?
#
#     These coefficients are already estimated inside model (3) -- they were
#     simply never reported. Confounding needs TWO arrows:
#       (i)  prior conflict -> tie at t-1   [the tie equation, section 7]
#       (ii) prior conflict -> conflict at t [HERE]
#     Showing only (i) licenses a MEDIATOR story, not a CONFOUNDING story.
#
#     Note the strong form: peace_years + its cubic are ALREADY in CTRL, so
#     conflict recurrence is partly conditioned on before the ladder begins.
#     Anything these terms pick up is recurrence structure that peace-years
#     does NOT absorb.
# ---------------------------------------------------------------------------
# The ladder estimates the dyadic and monadic processes in SEPARATE models, so
# no single ladder rung carries both prior-conflict terms. Leg 2 is a claim
# about the confounding arrows, not the ladder, so estimate a dedicated model
# that conditions on BOTH histories at once and report each coefficient.
HIST_TERMS <- c("dyad_mid_5yr_lag_sd", "mid5_sum_tp_sd")

hist_out <- purrr::imap_dfr(LAYERS, function(iv, lyr) {
  f <- as.formula(paste("cowmidonset ~", iv, "+", CTRL,
                        HIST_DYAD, HIST_MON, "| dyad_id + year"))
  m <- feols(f, data = est, vcov = ~dyad_id)
  ct <- coeftable(m)
  tibble(layer = lyr, term = HIST_TERMS,
         beta = ct[HIST_TERMS, "Estimate"],
         se   = ct[HIST_TERMS, "Std. Error"],
         p    = ct[HIST_TERMS, "Pr(>|t|)"])
})

cat("\n=====================================================================\n")
cat("LEG 2: prior conflict -> conflict at t, NET of the tie and peace-years\n")
cat("       (coefficients from model (3), by layer)\n")
cat("=====================================================================\n")
print(hist_out, digits = 4, n = Inf)
write_csv(hist_out, file.path(OUTPUT_DIR, "contamination_leg2.csv"))
cat("\n  If these are positive and significant, prior conflict independently\n")
cat("  predicts onset. Combined with the tie equation (section 7), that makes\n")
cat("  it a CONFOUNDER rather than a mediator.\n")

# ===========================================================================
# 6. GUARDS -- one per layer
# ===========================================================================
cat("\n--- Guard 1: does column (1) reproduce the pacifying finding? ---\n")
base <- res |> filter(model == labs[1])
for (i in seq_len(nrow(base))) {
  ok <- base$beta[i] < 0 && base$p[i] < 0.05
  cat(sprintf("  %-10s beta = %9.5f  p = %.2e  -> %s\n",
              base$layer[i], base$beta[i], base$p[i],
              if (ok) "PASS" else "FAIL (cannot anchor attenuation on this layer)"))
}

cat("\n--- Guard 2: no year-t leakage into the history regressors ---\n")
print(est |> summarise(
  cor_dyad5_onset = cor(dyad_mid_5yr_lag_sd, cowmidonset),
  cor_mid5tp_onset = cor(mid5_sum_tp_sd,      cowmidonset)), digits = 3)

# ===========================================================================
# 7. THE TIE EQUATION -- is each layer war-shaped?
#
#    This is NOT a 2SLS first stage. Nothing is instrumented. It establishes the
#    PREMISE: network ties are a product of war history, and therefore are not
#    exogenous predictors of conflict.
#
#    DYAD + YEAR FE -- matching the ladder. The ladder identifies the tie
#    coefficient off WITHIN-dyad variation, so the premise must be demonstrated
#    in that same variation. Showing "war shapes ties" between dyads while the
#    main model discards between-dyad variation would be incoherent.
#
#    NOTE: no PWE here. PWE's denominator IS the sum of tie weights, so
#    regressing a tie on PWE would put the tie on both sides. Circular.
# ===========================================================================
cat("\n=====================================================================\n")
cat("TIE EQUATION: does war history shape each layer? (dyad+year FE)\n")
cat("=====================================================================\n")
# Tie equation, run per FAMILY (conflict / capability / regime), each with its
# dyadic (between-unit) and monadic (within-unit) term, mirroring the ladder.
TIE_FAMILIES <- list(
  "Conflict history"           = c("dyad_mid_5yr_lag_sd", "mid5_sum_tp_sd"),
  "Capability shift (benchmark)" = c("cap_between_sd", "cap_within_sd"),
  "Regime change (benchmark)"    = c("regime_between_sd", "regime_within_sd")
)
TERM_LABELS <- c(
  dyad_mid_5yr_lag_sd = "Dyadic / between", mid5_sum_tp_sd    = "Monadic / within",
  cap_between_sd      = "Dyadic / between", cap_within_sd     = "Monadic / within",
  regime_between_sd   = "Dyadic / between", regime_within_sd  = "Monadic / within"
)

fs <- purrr::imap_dfr(
  c("Alliance" = "tie_alliance", "IGO" = "tie_igo", "Multiplex" = "tie_multiplex"),
  function(v, nm) {
    purrr::imap_dfr(TIE_FAMILIES, function(terms, famname) {
      rhs <- paste(terms, collapse = " + ")
      m <- feols(as.formula(paste(v, "~", rhs, "+", CTRL, "| dyad_id + year")),
                 data = est, vcov = ~dyad_id)
      ct <- coeftable(m)
      tibble(layer = nm, family = famname, term = terms,
             term_label = TERM_LABELS[terms],
             beta = ct[terms, "Estimate"], se = ct[terms, "Std. Error"],
             p    = ct[terms, "Pr(>|t|)"])
    })
  }) |>
  mutate(ci95_low = beta - 1.960 * se, ci95_high = beta + 1.960 * se,
         ci90_low = beta - 1.645 * se, ci90_high = beta + 1.645 * se)

print(fs, digits = 4, n = Inf)
cat("\n  Read: WITHIN a dyad, does prior war change the tie? If yes, the tie is\n")
cat("  not exogenous -- which is the premise the attenuation result rests on.\n")
cat("  The SIGN here predicts the DIRECTION the ladder moves: war that severs a\n")
cat("  tie (negative) should attenuate; war that builds one (positive) amplifies.\n")

# ===========================================================================
# 8. FIGURE: the two legs of the confounding argument, stacked.
#
#    PANEL A (top)    prior conflict -> TIE STRENGTH        [leg 1]
#    PANEL B (bottom) tie strength   -> MID ONSET, laddered [leg 3]
#
#    The rhetorical work is done by reading DOWN a column. The sign in A
#    PREDICTS the direction of movement in B: war that severs a tie (negative
#    in A) should attenuate the estimate in B; war that builds one (positive
#    in A) should amplify it. A reader who sees both rows does not need to be
#    told the argument.
#
#    Two plots + patchwork, NOT facet_grid: the panels are in different units
#    (tie strength vs. probability of onset), and facet_grid's free_x frees
#    scales by column only, which would force a shared x across the two rows.
# ===========================================================================
cat("\n--- Plotting ---\n")

library(patchwork)

FIG_DIR <- here::here("code", "03_figures", "output", "manuscript")
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

ACCENT  <- "#2C6CB0"
LAYER_LEVELS  <- c("Alliance", "IGO", "Multiplex")
# Step 1 (tie equation) keeps the conflict/capability/regime split; Step 2
# (ladder) collapses the two benchmarks into a single "Benchmark" family.
STEP1_FAMILIES <- c("Conflict history",
                    "Capability shift (benchmark)",
                    "Regime change (benchmark)")
STEP2_FAMILIES <- c("Baseline", "Conflict history", "Benchmark")
# One colour vector covering every family label used in either panel. Conflict =
# prominent blue; benchmarks/baseline = muted greys so they read as reference.
FAMILY_COLORS <- c("Baseline"                     = "#C2A25A",
                   "Conflict history"             = ACCENT,
                   "Benchmark"                    = "#8C8C8C",
                   "Capability shift (benchmark)" = "#8C8C8C",
                   "Regime change (benchmark)"    = "#8C8C8C")
fam_short <- c("Conflict history"           = "Conflict",
               "Capability shift (benchmark)" = "Capability",
               "Regime change (benchmark)"    = "Regime")

# --- Panel A: the tie equation, all families -----------------------------
row_lev <- unlist(lapply(STEP1_FAMILIES, function(f)
  paste0(fam_short[f], ":\n", c("Dyadic / between", "Monadic / within"))))
fs_plot <- fs |>
  mutate(layer  = factor(layer, levels = LAYER_LEVELS),
         family = factor(family, levels = STEP1_FAMILIES),
         row    = factor(paste0(fam_short[as.character(family)], ":\n", term_label),
                         levels = rev(row_lev)))

p_tie <- ggplot(fs_plot, aes(x = beta, y = row, color = family)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray55") +
  geom_linerange(aes(xmin = ci95_low, xmax = ci95_high), linewidth = 0.5) +
  geom_linerange(aes(xmin = ci90_low, xmax = ci90_high), linewidth = 1.6) +
  geom_point(size = 2.3) +
  facet_wrap(~ layer, ncol = 3, scales = "free_x") +
  scale_color_manual(values = FAMILY_COLORS, name = NULL) +
  scale_x_continuous(expand = expansion(mult = 0.12)) +
  scale_y_discrete(expand = expansion(add = 0.5)) +
  labs(subtitle = "Step 1: Prior conflict reshapes the tie — benchmarks do not",
       x = "Coefficient on tie strength (per SD)", y = NULL) +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(),
        strip.background = element_rect(fill = "gray92", color = "gray70"),
        strip.text = element_text(face = "bold"),
        axis.text.x = element_text(angle = 0, hjust = 0.5),
        axis.text.y = element_text(size = 8, lineheight = 0.9),
        legend.position = "none",
        plot.subtitle = element_text(face = "bold", size = 11))

# --- Panel B: the attenuation ladder + benchmark reference rungs ------------
plot_df <- res |>
  mutate(layer  = factor(layer, levels = LAYER_LEVELS),
         family = factor(family, levels = STEP2_FAMILIES),
         model  = factor(model, levels = rev(labs)))   # baseline at the top

# Annotate every non-baseline rung (the two conflict-history processes and the
# four benchmark rungs) with its attenuation/amplification, so the contrast is
# explicit: conflict moves the estimate, benchmarks barely do.
lab_df <- plot_df |>
  filter(model != labs[1]) |>
  mutate(txt = sprintf("%.2f%% %s", round(attenuation, 2),
                       if_else(attenuation >= 0, "attenuation", "amplification")))

# Per-facet reference line at the BASELINE coefficient (labs[1]). Each rung's
# horizontal distance from this line is the attenuation/amplification -- a far
# more informative anchor for Step 2 than zero. Zero is kept as a faint hairline
# so the reader can still see whether a coefficient crosses it.
base_df <- plot_df |>
  filter(model == labs[1]) |>
  select(layer, base_beta = beta)

p_att <- ggplot(plot_df, aes(x = beta, y = model, color = family)) +
  geom_vline(xintercept = 0, linetype = "dotted", color = "gray80", linewidth = 0.4) +
  geom_vline(data = base_df, aes(xintercept = base_beta),
             linetype = "dashed", color = "gray45", linewidth = 0.6) +
  geom_linerange(aes(xmin = ci95_low, xmax = ci95_high), linewidth = 0.5) +
  geom_linerange(aes(xmin = ci90_low, xmax = ci90_high), linewidth = 1.6) +
  geom_point(size = 2.4) +
  geom_text(data = lab_df, aes(x = beta, y = model, label = txt),
            hjust = 0.5, vjust = -1.6, size = 2.9, color = "gray25",
            inherit.aes = FALSE, show.legend = FALSE) +
  facet_wrap(~ layer, ncol = 3, scales = "free_x") +
  scale_color_manual(values = FAMILY_COLORS, name = NULL) +
  scale_x_continuous(labels = scales::label_number(accuracy = 1e-4,
                                                   drop0trailing = TRUE),
                     expand = expansion(mult = 0.12)) +
  scale_y_discrete(expand = expansion(add = 0.5)) +
  labs(subtitle = "Step 2: So the tie-conflict estimate moves",
       x = "Coefficient on network tie (per SD)", y = NULL) +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(),
        strip.background = element_rect(fill = "gray92", color = "gray70"),
        strip.text = element_text(face = "bold"),
        axis.text.x = element_text(angle = 0, hjust = 0.5),
        axis.text.y = element_text(size = 8, lineheight = 0.9),
        legend.position = "bottom",
        plot.subtitle = element_text(face = "bold", size = 11))

# --- Assemble and export --------------------------------------------------
p_full <- (p_tie / p_att) +
  plot_annotation(
    title = "Conflict history confounds the network-conflict relationship",
    subtitle = paste0(
      "Dyad + year fixed effects. Conflict history measured t-2 to t-6, ",
      "strictly prior to the tie at t-1.\n",
      "Capability and regime shifts are non-conflict benchmarks over the ",
      "same window.\nThick bar = 90% CI, thin bar = 95% CI."),
    theme = theme(plot.title = element_text(face = "bold", size = 14),
                  plot.subtitle = element_text(color = "gray30", size = 10))
  )

ggsave(file.path(FIG_DIR, "fig_contamination_ladder.png"), p_full,
       width = 11, height = 9, dpi = 400, bg = "white")
cat("  -> fig_contamination_ladder.png\n")

# ===========================================================================
# 9. APPENDIX: alternative monadic parameterizations (sum vs max vs min+max)
#    sum = min + max identically, so these are three views of the same
#    information. The dyadic term is held fixed across all three.
#    NB: these use the FULL monadic sum (all onsets, incl. i-j) so the
#    min/max identity holds. The MAIN ladder (section 5) uses the THIRD-PARTY
#    monadic sum instead; this appendix is a separate carve-up robustness.
# ===========================================================================
cat("\n=====================================================================\n")
cat("APPENDIX: monadic parameterization -- sum vs. max vs. (min + max)\n")
cat("=====================================================================\n")

MON_SPECS <- list(
  "A. sum (all onsets)" = "+ mid5_sum_sd",
  "B. max"              = "+ mid5_max_sd",
  "C. min + max"        = "+ mid5_min_sd + mid5_max_sd"
)

# --- 9a. Tie equation (Step 1) under each parameterization ----------------
app_tie <- purrr::imap_dfr(MON_SPECS, function(mon, spec_nm) {
  purrr::imap_dfr(
    c("Alliance" = "tie_alliance", "IGO" = "tie_igo", "Multiplex" = "tie_multiplex"),
    function(v, lyr) {
      f <- as.formula(paste(v, "~ dyad_mid_5yr_lag_sd", mon, "+", CTRL,
                            "| dyad_id + year"))
      m  <- feols(f, data = est, vcov = ~dyad_id)
      ct <- coeftable(m)
      tms <- rownames(ct)[rownames(ct) %in%
                            c("dyad_mid_5yr_lag_sd", "mid5_sum_sd",
                              "mid5_max_sd", "mid5_min_sd")]
      tibble(spec = spec_nm, layer = lyr, term = tms,
             beta = ct[tms, "Estimate"], se = ct[tms, "Std. Error"],
             p    = ct[tms, "Pr(>|t|)"])
    })
})
cat("\n--- 9a. TIE EQUATION (Step 1) ---\n")
print(app_tie, digits = 4, n = Inf)

# --- 9b. Attenuation ladder (Step 2) under each parameterization ----------
app_lad <- purrr::imap_dfr(MON_SPECS, function(mon, spec_nm) {
  purrr::imap_dfr(LAYERS, function(iv, lyr) {
    a_full <- paste(HIST_DYAD, mon)
    b0 <- feols(as.formula(paste("cowmidonset ~", iv, "+", CTRL,
                                 "| dyad_id + year")),
                data = est, vcov = ~dyad_id)
    b3 <- feols(as.formula(paste("cowmidonset ~", iv, "+", CTRL, a_full,
                                 "| dyad_id + year")),
                data = est, vcov = ~dyad_id)
    beta0 <- coeftable(b0)[iv, "Estimate"]
    beta3 <- coeftable(b3)[iv, "Estimate"]
    tibble(spec = spec_nm, layer = lyr,
           beta_base = beta0, beta_cond = beta3,
           p_cond    = coeftable(b3)[iv, "Pr(>|t|)"],
           attenuation = 100 - 100 * abs(beta3) / abs(beta0))
  })
})
cat("\n--- 9b. ATTENUATION LADDER (Step 2), fully conditioned model ---\n")
cat("    Positive = attenuation (toward zero). Negative = amplification.\n\n")
print(app_lad, digits = 4, n = Inf)

write_csv(app_tie, file.path(OUTPUT_DIR, "contamination_appendix_tie.csv"))
write_csv(app_lad, file.path(OUTPUT_DIR, "contamination_appendix_ladder.csv"))

# ===========================================================================
# 10. Save
# ===========================================================================
write_csv(res, file.path(OUTPUT_DIR, "contamination_demo.csv"))
cat("\nSaved -> data/processed/contamination_demo.csv\n")
cat("Saved -> data/processed/contamination_appendix_tie.csv\n")
cat("Saved -> data/processed/contamination_appendix_ladder.csv\n")

cat("\n---------------------------------------------------------------------\n")
cat("FRAMING (do not overclaim):\n")
cat("  I do not estimate the causal effect of ties on conflict. I show that the\n")
cat("  network-conflict association is sensitive to conditioning on prior\n")
cat("  conflict histories, which is precisely what we should expect if conflict\n")
cat("  helps produce the ties later used as predictors.\n")
cat("---------------------------------------------------------------------\n")
