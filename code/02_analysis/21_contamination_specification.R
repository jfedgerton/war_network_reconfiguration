###############################################################################
# 22_contamination_specification.R
#
# SPECIFICATION CHECKS FOR THE CONTAMINATION RESULT
# -------------------------------------------------
# Sibling to 21 (falsification). Where 21 asks "could this be an artifact?",
# this script asks "does it survive reasonable alternative choices?".
#
#   (D) ESTIMATOR -- DELIBERATELY NOT TESTED. Reasoning, so nobody re-adds it:
#
#       FE-LOGIT does not test this claim. It drops every dyad whose onset
#       outcome never varies; with a ~0.3% onset rate that is most of the panel,
#       leaving a non-random subset of dyads that fought at least once. It is
#       therefore not estimating the same quantity on the same population -- a
#       "pass" would be weak evidence and a "fail" uninterpretable (link
#       function, or sample?). Its coefficients are log-odds and the LPM's are
#       probability differences, so only the attenuation RATIO is comparable.
#
#       POOLED LOGIT (the literature's own workhorse) would fix the sample
#       problem but requires dropping the dyad fixed effects -- which
#       reintroduces exactly the unobserved time-invariant rivalry (U) that the
#       dyad FE exist to absorb. That trades a sample problem for a confounding
#       problem.
#
#       THE ACTUAL DEFENCE: the contamination claim is a claim about OMITTED
#       VARIABLE BIAS, and OVB (beta_hat = beta + gamma*delta) is a property of
#       the LINEAR PROJECTION. Whether the underlying DGP is logistic is beside
#       the point: the question is whether THE ESTIMATE SCHOLARS REPORT absorbs
#       conflict history, and the at-risk LPM with dyad + year FE is the
#       workhorse specification of the dyadic-conflict literature. Changing the
#       link function does not test the claim; it changes the subject.
#
#   (E) SAMPLE
#       (1) at-risk        : drop dyad-years already in an ongoing MID (MAIN)
#       (2) all dyads      : no at-risk restriction
#       (3) politically relevant dyads (PRD): contiguous OR at least one major
#           power. The dyadic-conflict literature's standard sample. If the
#           contamination result only appears in the full 1.9M-dyad panel of
#           mostly-irrelevant pairs, it is not a claim about anybody's
#           published estimates.
#
#   (F) OVER-CONTROL CHECK: TIE AT t-2
#       The Discussion footnote concedes that because ties are persistent, the
#       tie at t-1 is close to the tie at t-2, so conditioning on conflict at
#       t-2..t-6 may block part of a path running from the EARLIER tie, through
#       conflict, to onset -- i.e. some of the attenuation could be over-control
#       rather than bias correction.
#
#       Test: move the tie to t-2, so there is a one-year gap between the tie
#       and the start of the conflict-history window... no: at t-2 the tie is
#       CONTEMPORANEOUS with the first history lag, which is worse. So instead
#       shift the HISTORY WINDOW BACK, to t-3..t-7, opening a two-year gap
#       between the tie (t-1) and the nearest history year (t-3).
#
#       If the attenuation survives with a wider gap, the persistence story is
#       weaker: the further the history is from the tie, the less scope there
#       is for the tie itself to have produced it.
#
# Main-text specification otherwise (unweighted, raw, alliance/IGO/multiplex).
#
# Outputs:
#   data/processed/contamination_specification.csv
#   code/03_figures/output/manuscript/fig_ms_contamination_specification.png
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

# ===========================================================================
# 1. Dyad-year panel
# ===========================================================================
cat("Building dyad-year panel...\n")

d <- create_dyadyears(subset_years = YR_MIN:YR_MAX, directed = FALSE) |>
  add_cow_mids()   |>
  add_contiguity() |>          # needed for the PRD sample (E3)
  add_democracy()  |>
  add_nmc()        |>
  add_cow_majors() |>
  add_peace_years()

py_col <- intersect(c("cowmidspell", "gmlmidspell", "midspell"), names(d))[1]
d$peace_years <- d[[py_col]]

contig_col <- intersect(c("conttype", "contig"), names(d))[1]

d <- d |>
  mutate(
    dyad_id         = paste0(pmin(ccode1, ccode2), "_", pmax(ccode1, ccode2)),
    # Dyadic capability asymmetry, bounded [0,1]; matches 03/19/20.
    cap_ratio       = abs(cinc1 - cinc2) / (cinc1 + cinc2),
    joint_democracy = as.integer(pmin(polity21, polity22, na.rm = TRUE) >= 6),
    major_power     = as.integer(cowmaj1 == 1 | cowmaj2 == 1),
    is_contig       = as.integer(coalesce(.data[[contig_col]], 0) > 0),
    # Politically relevant dyad: contiguous OR at least one major power.
    is_prd          = as.integer(is_contig == 1 | major_power == 1)
  )

CTRL <- "cap_ratio + joint_democracy + major_power +
         peace_years + I(peace_years^2) + I(peace_years^3)"

# ===========================================================================
# 2. Ties (main-text: unweighted, raw)
# ===========================================================================
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
  summarise(across(starts_with("tie_"), ~ max(.x, na.rm = TRUE)), .groups = "drop")

d <- d |>
  left_join(ties, by = c("ccode1", "ccode2", "year")) |>
  mutate(across(starts_with("tie_"), ~ replace_na(.x, 0))) |>
  arrange(dyad_id, year) |>
  group_by(dyad_id) |>
  mutate(across(c(tie_alliance, tie_igo, tie_multiplex),
                ~ dplyr::lag(.x, 1), .names = "{.col}_lag")) |>
  ungroup()

# ===========================================================================
# 3. Conflict history: main window (t-2..t-6) and the gapped window (t-3..t-7)
# ===========================================================================
# Matrix-free windowed lag sum. rowSums(sapply(...)) errors on single-row groups
# (a state/dyad present for one year); Reduce over a list is size-safe.
lag_sum <- function(x, ks) Reduce(`+`, lapply(ks, function(k) dplyr::lag(x, k, default = 0)))

d <- d |>
  arrange(dyad_id, year) |>
  group_by(dyad_id) |>
  mutate(
    dyad_hist      = lag_sum(cowmidonset, 2:6),   # MAIN
    dyad_hist_gap  = lag_sum(cowmidonset, 3:7)    # (F) wider gap from the tie
  ) |>
  ungroup()

state_churn <- readRDS(file.path(OUTPUT_DIR, "state_network_churn.rds"))
state_hist <- state_churn |>
  select(ccode, year, mid_onset_count) |>
  arrange(ccode, year) |>
  group_by(ccode) |>
  mutate(mon      = lag_sum(mid_onset_count, 2:6),
         mon_gap  = lag_sum(mid_onset_count, 3:7)) |>
  ungroup() |>
  select(ccode, year, mon, mon_gap)

d <- d |>
  left_join(state_hist |> rename(m_i = mon, mg_i = mon_gap),
            by = c("ccode1" = "ccode", "year")) |>
  left_join(state_hist |> rename(m_j = mon, mg_j = mon_gap),
            by = c("ccode2" = "ccode", "year")) |>
  mutate(across(c(m_i, m_j, mg_i, mg_j), ~ replace_na(.x, 0)),
         # THIRD-PARTY monadic (mirrors 03): own onsets minus i-j onsets, over
         # the matching window. dyad_hist / dyad_hist_gap count the i-j onsets.
         mon_hist     = pmax(m_i  - dyad_hist,     0) + pmax(m_j  - dyad_hist,     0),
         mon_hist_gap = pmax(mg_i - dyad_hist_gap, 0) + pmax(mg_j - dyad_hist_gap, 0))

# ===========================================================================
# 3b. Benchmark shift histories (for the no-controls appendix ladder, sec 9)
#     Ported from 03: capability and regime SHIFTS over t-2..t-6, dyadic
#     (between) and monadic (within). These are shifts, not counts, so they are
#     NOT log1p'd -- only per-SD standardized (see section 4).
# ===========================================================================
lag5 <- function(x) {
  x <- tidyr::replace_na(x, 0)
  Reduce(`+`, lapply(2:6, function(k) dplyr::lag(x, k, default = 0)))
}

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
# 4. Samples (E)
# ===========================================================================
need <- c("cowmidonset", "tie_alliance_lag", "tie_igo_lag", "tie_multiplex_lag",
          "cap_ratio", "joint_democracy", "major_power", "peace_years",
          "dyad_hist", "mon_hist", "dyad_hist_gap", "mon_hist_gap")

base <- d |> filter(if_all(all_of(need), ~ !is.na(.x)))

sd_safe <- function(x) { s <- sd(x, na.rm = TRUE); if (is.na(s) || s == 0) 1 else s }
standardize <- function(dat) {
  dat |> mutate(
    # Ties + benchmark shifts: plain per-SD.
    across(c(tie_alliance_lag, tie_igo_lag, tie_multiplex_lag,
             cap_between, cap_within, regime_between, regime_within),
           ~ .x / sd_safe(.x), .names = "{.col}_sd"),
    # Conflict-history onset counts: log1p (mirrors main analyses) then per-SD.
    across(c(dyad_hist, mon_hist, dyad_hist_gap, mon_hist_gap),
           ~ { z <- log1p(.x); z / sd_safe(z) }, .names = "{.col}_sd"))
}

# PRD sample removed [2026-07-15] per author request. The two retained samples
# are the at-risk main sample and the unrestricted all-dyads sample.
SAMPLES <- list(
  "At-risk (MAIN)" = base |>
    filter(coalesce(cowmidongoing, 0L) == 0 | cowmidonset == 1) |> standardize(),
  "All dyads"      = base |> standardize()
)

for (nm in names(SAMPLES)) {
  s <- SAMPLES[[nm]]
  cat(sprintf("  %-28s N = %10s   onset rate %.3f%%   dyads %s\n",
              nm, format(nrow(s), big.mark = ","),
              100 * mean(s$cowmidonset),
              format(dplyr::n_distinct(s$dyad_id), big.mark = ",")))
}

LAYERS <- c("Alliance"  = "tie_alliance_lag_sd",
            "IGO"       = "tie_igo_lag_sd",
            "Multiplex" = "tie_multiplex_lag_sd")

# ===========================================================================
# 5. Ladder runner -- LPM or FE-logit, main window or gapped window
# ===========================================================================
run_ladder <- function(dat, iv, window) {

  # Dyadic and monadic conflict-history processes, conditioned SEPARATELY (each
  # added to the same baseline on its own), mirroring the main-text 03 ladder.
  hist_terms <- if (window == "main")
    c(Dyadic = "dyad_hist_sd",     Monadic = "mon_hist_sd")
  else
    c(Dyadic = "dyad_hist_gap_sd", Monadic = "mon_hist_gap_sd")

  f0 <- as.formula(paste("cowmidonset ~", iv, "+", CTRL, "| dyad_id + year"))
  m0 <- tryCatch(feols(f0, data = dat, vcov = ~dyad_id), error = function(e) NULL)
  if (is.null(m0) || !iv %in% rownames(coeftable(m0))) return(NULL)
  b0 <- coeftable(m0)[iv, "Estimate"]; p0 <- coeftable(m0)[iv, "Pr(>|t|)"]

  purrr::imap_dfr(hist_terms, function(term, proc) {
    f1 <- as.formula(paste("cowmidonset ~", iv, "+", CTRL, "+", term,
                           "| dyad_id + year"))
    m1 <- tryCatch(feols(f1, data = dat, vcov = ~dyad_id), error = function(e) NULL)
    if (is.null(m1) || !iv %in% rownames(coeftable(m1)))
      return(tibble(process = proc, window = window, n_input = nrow(dat),
                    n_used = NA, beta_base = b0, p_base = p0,
                    beta_cond = NA, p_cond = NA, attenuation_pct = NA))
    b1 <- coeftable(m1)[iv, "Estimate"]; p1 <- coeftable(m1)[iv, "Pr(>|t|)"]
    tibble(process = proc, window = window, n_input = nrow(dat),
           n_used = nobs(m1), beta_base = b0, p_base = p0,
           beta_cond = b1, p_cond = p1,
           attenuation_pct = 100 * (1 - abs(b1) / abs(b0)))
  })
}

# ===========================================================================
# 6. Run the grid  (samples x estimators x windows x layers)
# ===========================================================================
cat("\n=====================================================================\n")
cat("E: samples  |  F: gapped history window     (D dropped -- see header)\n")
cat("=====================================================================\n")

GRID <- tidyr::expand_grid(
  sample_name = names(SAMPLES),
  window      = c("main", "gap"),
  layer       = names(LAYERS)
)

spec <- purrr::pmap_dfr(GRID, function(sample_name, window, layer) {
  cat(sprintf("  %-28s %-5s %-10s ...\n", sample_name, window, layer))
  r <- run_ladder(SAMPLES[[sample_name]], LAYERS[[layer]], window)
  if (is.null(r)) {
    return(tibble(sample_name = sample_name, layer = layer, window = window,
                  process = c("Dyadic", "Monadic"),
                  n_input = NA, n_used = NA,
                  beta_base = NA, p_base = NA, beta_cond = NA, p_cond = NA,
                  attenuation_pct = NA))
  }
  r |> mutate(sample_name = sample_name, layer = layer, .before = 1)
})

print(spec |> select(sample_name, layer, window, process,
                     n_used, attenuation_pct, p_base, p_cond),
      digits = 3, n = Inf)

write_csv(spec, file.path(OUTPUT_DIR, "contamination_specification.csv"))

# ===========================================================================
# 7. What to read
# ===========================================================================
cat("\n--- (E) SAMPLE SENSITIVITY (main window) ---\n")
print(spec |> filter(window == "main") |>
        select(sample_name, layer, process, n_used, attenuation_pct, p_base),
      digits = 3, n = Inf)
cat("\n  Alliance and multiplex should attenuate; IGO should amplify (dyadic),\n")
cat("  across both the at-risk and all-dyads samples.\n")

cat("\n--- (F) OVER-CONTROL: history moved back to t-3..t-7 ---\n")
gapcmp <- spec |>
  filter(sample_name == "At-risk (MAIN)") |>
  select(layer, process, window, attenuation_pct) |>
  tidyr::pivot_wider(names_from = window, values_from = attenuation_pct) |>
  mutate(retained_pct = 100 * gap / main)
print(gapcmp, digits = 3, n = Inf)
cat("\n  If `retained_pct` is high, the attenuation survives when the conflict\n")
cat("  history is pushed further from the tie -- which weakens the tie-\n")
cat("  persistence / over-control objection in the Discussion footnote.\n")

# ===========================================================================
# 8. FIGURE
# ===========================================================================
cat("\n--- Plotting ---\n")
ACCENT <- "#2C6CB0"

p <- spec |>
  filter(!is.na(attenuation_pct)) |>
  mutate(
    spec_label = paste0(sample_name, "\n",
                        ifelse(window == "main", "history t-2..t-6",
                                                 "history t-3..t-7")),
    layer   = factor(layer, levels = c("Alliance", "IGO", "Multiplex")),
    process = factor(process, levels = c("Dyadic", "Monadic"),
                     labels = c("Dyadic conflict\nhistory (i-j MIDs)",
                                "Monadic conflict\nhistory (third-party)"))
  ) |>
  ggplot(aes(x = attenuation_pct, y = spec_label)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray55") +
  geom_col(fill = ACCENT, width = 0.62) +
  facet_grid(process ~ layer, scales = "free_x") +
  labs(
    title = "Contamination survives the sample and the history window",
    subtitle = paste0("Attenuation of the tie-conflict coefficient, dyadic and monadic ",
                      "conflict-history processes separately.\nPositive = toward zero, ",
                      "negative = away from zero. The t-3..t-7 rows open a two-year gap ",
                      "between the\ntie (t-1) and the nearest history year, testing the ",
                      "tie-persistence / over-control objection."),
    x = "Attenuation (%)", y = NULL
  ) +
  theme_bw(base_size = 10) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(),
        strip.background = element_rect(fill = "gray92", color = "gray70"),
        strip.text = element_text(face = "bold"),
        #strip.text.y = element_text(angle = 0),
        axis.text.x = element_text(angle = 0, hjust = 0.5),
        axis.text.y = element_text(size = 8, lineheight = 0.95),
        plot.title = element_text(face = "bold"),
        plot.subtitle = element_text(color = "gray30", size = 8.5))

ggsave(file.path(FIG_DIR, "fig_ms_contamination_specification.png"), p,
       width = 11, height = 8.5, dpi = 400, bg = "white")
cat("  -> fig_ms_contamination_specification.png\n")

# ===========================================================================
# 9. APPENDIX: no-controls ladder (dyad + year FE only)
#
#    Repeats the MAIN 03 contamination ladder with the CTRL block DROPPED -- only
#    the dyad + year fixed effects and the conditioning term of interest remain.
#    This directly answers the pre-control-asymmetry objection: with cap_ratio /
#    joint_democracy no longer in the baseline, the capability and regime
#    benchmarks can no longer be pre-absorbed. Same three tie inputs as 03 --
#    combined alliance, IGO, and multiplex. Each rung conditions on ONE term
#    (added to the tie + FE only); attenuation is measured vs the no-controls
#    baseline. Estimated on the at-risk main sample.
# ===========================================================================
cat("\n=====================================================================\n")
cat("APPENDIX: no-controls ladder (dyad + year FE only) -- alliance/IGO/multiplex\n")
cat("=====================================================================\n")

# Same layers as the main 03 ladder (combined alliance, IGO, multiplex).
NC_LAYERS <- LAYERS

# Six conditioning terms, each added on its own (mirrors 03's rung structure).
NC_TERMS <- list(
  "Dyadic conflict history"        = "dyad_hist_sd",
  "Monadic conflict history"       = "mon_hist_sd",
  "Between polity (regime gap)"    = "regime_between_sd",
  "Monadic polity (own regime)"    = "regime_within_sd",
  "Between CINC (asymmetry)"       = "cap_between_sd",
  "Monadic CINC (own capability)"  = "cap_within_sd"
)

nc_dat <- SAMPLES[["At-risk (MAIN)"]]

run_nc <- function(dat, iv) {
  b0 <- tryCatch(feols(as.formula(paste("cowmidonset ~", iv, "| dyad_id + year")),
                       data = dat, vcov = ~dyad_id), error = function(e) NULL)
  if (is.null(b0) || !iv %in% rownames(coeftable(b0))) return(NULL)
  beta0 <- coeftable(b0)[iv, "Estimate"]; p0 <- coeftable(b0)[iv, "Pr(>|t|)"]
  purrr::imap_dfr(NC_TERMS, function(term, nm) {
    m <- tryCatch(feols(as.formula(paste("cowmidonset ~", iv, "+", term,
                                        "| dyad_id + year")),
                        data = dat, vcov = ~dyad_id), error = function(e) NULL)
    if (is.null(m) || !iv %in% rownames(coeftable(m)))
      return(tibble(term = nm, beta_base = beta0, p_base = p0,
                    beta_cond = NA, p_cond = NA, attenuation_pct = NA))
    bc <- coeftable(m)[iv, "Estimate"]; pc <- coeftable(m)[iv, "Pr(>|t|)"]
    tibble(term = nm, beta_base = beta0, p_base = p0,
           beta_cond = bc, p_cond = pc,
           attenuation_pct = 100 * (1 - abs(bc) / abs(beta0)))
  })
}

nc <- purrr::imap_dfr(NC_LAYERS, function(iv, lyr)
  run_nc(nc_dat, iv) |> mutate(layer = lyr, .before = 1))

print(nc |> select(layer, term, beta_base, beta_cond, attenuation_pct, p_cond),
      digits = 3, n = Inf)
write_csv(nc, file.path(OUTPUT_DIR, "contamination_nocontrols_ladder.csv"))

# --- Figure: no-controls ladder, one facet per layer (alliance/IGO/multiplex) ---
nc_plot <- nc |>
  filter(!is.na(attenuation_pct)) |>
  mutate(
    kind  = if_else(grepl("conflict", term), "Conflict history", "Benchmark"),
    term  = factor(term, levels = rev(names(NC_TERMS))),
    layer = factor(layer, levels = names(NC_LAYERS)))

p_nc <- ggplot(nc_plot, aes(x = attenuation_pct, y = term, fill = kind)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray55") +
  geom_col(width = 0.62) +
  facet_wrap(~ layer, ncol = 3, scales = "free_x") +
  scale_fill_manual(name = NULL,
                    values = c("Conflict history" = ACCENT, "Benchmark" = "#8C8C8C")) +
  labs(
    title = "No-controls ladder: contamination without the pre-controlled levels",
    subtitle = paste0("Attenuation of the tie-conflict coefficient with dyad + year fixed ",
                      "effects only (no controls).\nEach rung conditions on ONE term. ",
                      "Positive = toward zero, negative = away from zero."),
    x = "Attenuation (%)", y = NULL) +
  theme_bw(base_size = 10) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(),
        strip.background = element_rect(fill = "gray92", color = "gray70"),
        strip.text = element_text(face = "bold"),
        axis.text.x = element_text(angle = 0, hjust = 0.5),
        axis.text.y = element_text(size = 8, lineheight = 0.95),
        legend.position = "bottom",
        plot.title = element_text(face = "bold"),
        plot.subtitle = element_text(color = "gray30", size = 8.5))

ggsave(file.path(FIG_DIR, "fig_app_contamination_nocontrols.png"), p_nc,
       width = 11, height = 6, dpi = 400, bg = "white")
cat("  -> fig_app_contamination_nocontrols.png\n")

cat("\nSaved -> contamination_specification.csv\n")
cat("Saved -> contamination_nocontrols_ladder.csv\n")
