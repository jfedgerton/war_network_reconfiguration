###############################################################################
# 09_mediation_battle_deaths.R (formerly 09_mediation_power_delta.R)
#
# Mediation / informativeness analysis: War -> Battle Deaths -> log(Churn).
#
# Battle deaths is treated here as an *informativeness* / intensity measure
# of the conflict signal, not as an alternative mechanism. Under the
# information framing of the paper, more lethal wars reveal more
# capability information; the question this script asks is whether the
# magnitude of information revealed (proxied by log(1 + battle deaths))
# accounts for the war -> network-reconfiguration association.
#
# NOTE on scope: this script previously also ran ΔCINC mediation as an
# alternative-mechanism test. That analysis has been moved to
# 16_robustness_alt_mechanisms.R, which collects ΔCINC alongside the
# other downstream-channel mediators (|ΔPolity2|, major-power flip,
# territorial change). Battle deaths is *not* a downstream channel --
# it is a granular measure of the IV itself -- so it lives here.
#
# Structure:
#   1. Baron-Kenny: battle deaths mediator, full panel (4 IV specs)
#   2. Imai-Keele-Tingley causal mediation (mediation package)
#   3. Sensitivity analysis for sequential ignorability
#   4. Save results + LaTeX tables
#
# DV: log(1 + Churn) (matches 02 main spec)
# Mediator: log(1 + battle deaths)
# Controls: partner_war_frac + democracy + autocracy + cinc + major_power + is_isolated_t1
# All panel models: two-way FE (ccode + year), clustered SEs (~ccode + year)
###############################################################################

library(tidyverse)
library(fixest)
library(mediation)
library(texreg)
select <- dplyr::select

set.seed(123)

OUTPUT_DIR <- here::here("data", "processed")
FIG_DIR    <- here::here("code", "03_figures", "output")

cat("=== Mediation Analysis: War → Power Delta → Churn ===\n\n")

# ============================================================
# 0. Load and prepare data
# ============================================================

cat("Loading data...\n")
state_churn <- readRDS(file.path(OUTPUT_DIR, "state_network_churn.rds"))
dyad_data   <- readRDS(file.path(OUTPUT_DIR, "dyad_year_network_layers.rds"))

# --- State-level CINC ---
# CINC is kept as a CONTROL (matches 02 main spec). ΔCINC mediation
# was moved to 16_robustness_alt_mechanisms.R; the lagged / pct / abs
# Δ-CINC variables that powered that block are no longer computed here.
state_cinc <- dyad_data |>
  select(ccode = ccode1, year, cinc = cinc1) |>
  distinct(ccode, year, .keep_all = TRUE)

state_churn <- state_churn |>
  left_join(state_cinc, by = c("ccode", "year"))

# --- State controls ---
state_attrs <- dyad_data |>
  select(ccode = ccode1, year, cowmaj1) |>
  distinct(ccode, year, .keep_all = TRUE) |>
  rename(major_power = cowmaj1)

state_churn <- state_churn |>
  left_join(state_attrs, by = c("ccode", "year")) |>
  mutate(major_power = replace_na(major_power, 0))

# Democracy/autocracy from peacesciencer
library(peacesciencer)
state_years <- create_stateyears(
  subset_years = min(state_churn$year):max(state_churn$year)
)
attr(state_years, "ps_data_type") <- "state_year"
attr(state_years, "ps_system") <- "cow"
state_years <- state_years |> add_democracy()

polity_data <- state_years |>
  select(ccode, year, polity2) |>
  mutate(
    democracy = as.integer(polity2 >= 6),
    autocracy = as.integer(polity2 <= -6)
  )

state_churn <- state_churn |>
  left_join(polity_data, by = c("ccode", "year"))

# --- Log-transformed IVs and DV (matches 02 main spec) ---
state_churn <- state_churn |>
  mutate(
    log_churn                 = log1p(churn),
    log_mid_onset_count       = log1p(mid_onset_count),
    log_fatal_mid_onset_count = log1p(fatal_mid_onset_count)
  )

# --- Weighted partner war exposure (matches 02 main spec) ---
cat("Computing weighted partner war exposure...\n")
.MAIN_W_VAR <- "main_weight"
if (!.MAIN_W_VAR %in% names(dyad_data)) {
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

partner_sets_09 <- dyad_data |>
  filter(.data[[.MAIN_W_VAR]] > 0) |>
  select(ccode1, ccode2, year, w = !!sym(.MAIN_W_VAR)) |>
  bind_rows(
    dyad_data |> filter(.data[[.MAIN_W_VAR]] > 0) |>
      select(ccode1 = ccode2, ccode2 = ccode1, year, w = !!sym(.MAIN_W_VAR))
  )

all_state_war_09 <- dyad_data |>
  filter(cowmidongoing == 1) |>
  select(ccode1, ccode2, year) |>
  pivot_longer(cols = c(ccode1, ccode2), values_to = "ccode") |>
  distinct(ccode, year) |>
  mutate(at_war = 1L)

all_state_fatal_war_09 <- dyad_data |>
  filter(cowmidfatal == 1) |>
  select(ccode1, ccode2, year) |>
  pivot_longer(cols = c(ccode1, ccode2), values_to = "ccode") |>
  distinct(ccode, year) |>
  mutate(at_war_fatal = 1L)

pwe_09 <- partner_sets_09 |>
  mutate(year_next = year + 1L) |>
  left_join(all_state_war_09,       by = c("ccode2" = "ccode", "year_next" = "year")) |>
  left_join(all_state_fatal_war_09, by = c("ccode2" = "ccode", "year_next" = "year")) |>
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
  left_join(pwe_09, by = c("ccode", "year")) |>
  mutate(partner_war_frac                = replace_na(partner_war_frac, 0),
         partner_war_frac_weighted       = replace_na(partner_war_frac_weighted, 0),
         partner_fatal_war_frac          = replace_na(partner_fatal_war_frac, 0),
         partner_fatal_war_frac_weighted = replace_na(partner_fatal_war_frac_weighted, 0))

# --- Partner-war control routing by IV severity ---
# Non-fatal IVs (any MID exposure) -> partner_war_frac_weighted
# Fatal IVs (battle-death-bearing MIDs) -> partner_fatal_war_frac_weighted
# State-level churn regressions in this script use the tie-WEIGHTED partner-
# war exposure (alliance-strength-weighted share of partners at war), which is
# the harmonized robustness-script default. The unweighted variants are
# reserved for the main 02_state_level_analysis.R headline tables.
CTRL_MID   <- "partner_war_frac_weighted + democracy + autocracy + cinc + major_power + is_isolated_t1"
CTRL_FATAL <- "partner_fatal_war_frac_weighted + democracy + autocracy + cinc + major_power + is_isolated_t1"
ctrl_for_iv <- function(iv) if (grepl("fatal", iv)) CTRL_FATAL else CTRL_MID

# --- Lagged IVs / lead DV (kept for downstream temporal robustness) ---
state_churn <- state_churn |>
  group_by(ccode) |>
  arrange(year) |>
  mutate(
    mid_onset_lag1       = lag(mid_onset_count, 1),
    fatal_onset_lag1     = lag(fatal_mid_onset_count, 1),
    log_onset_lag1       = lag(log_mid_onset_count, 1),
    log_fatal_onset_lag1 = lag(log_fatal_mid_onset_count, 1),
    log_churn_lead1      = lead(log_churn, 1)
  ) |>
  ungroup()

# --- COW interstate wars: battle deaths ---
cat("Loading COW interstate war data...\n")
data(cow_war_inter, package = "peacesciencer")

war_state_year <- cow_war_inter |>
  select(warnum, ccode = ccode1, year, outcome = outcome1, batdeath = batdeath1) |>
  distinct(warnum, ccode, year, .keep_all = TRUE) |>
  group_by(ccode, year) |>
  summarise(
    in_cow_war     = 1L,
    n_wars         = n_distinct(warnum),
    total_batdeath = sum(batdeath, na.rm = TRUE),
    .groups = "drop"
  ) |>
  mutate(log_batdeath = log1p(total_batdeath))

state_churn <- state_churn |>
  left_join(war_state_year, by = c("ccode", "year")) |>
  mutate(
    in_cow_war     = replace_na(in_cow_war, 0L),
    n_wars         = replace_na(n_wars, 0L),
    total_batdeath = replace_na(total_batdeath, 0),
    log_batdeath   = replace_na(log_batdeath, 0)
  )

# Complete cases for mediation
med_data <- state_churn |>
  filter(!is.na(churn), !is.na(mid_onset_count), !is.na(log_batdeath)) |>
  mutate(
    ccode_f = factor(ccode),
    year_f  = factor(year)
  )

cat(sprintf("  Mediation panel: %d state-years, %d states, %d years\n",
            nrow(med_data), n_distinct(med_data$ccode), n_distinct(med_data$year)))

# ============================================================
# Helper: run Baron-Kenny for a given dataset + IV
# ============================================================

run_bk <- function(data, iv, mediator, iv_label, tag_prefix = "") {
  tag  <- paste0(tag_prefix, iv, "__log_churn")
  ctrl <- ctrl_for_iv(iv)

  f_total <- as.formula(paste("log_churn ~", iv, "+", ctrl,
                                "| ccode + year"))
  m_total <- tryCatch(feols(f_total, data = data, vcov = ~ccode + year),
                      error = function(e) NULL)

  f_a <- as.formula(paste(mediator, "~", iv, "+", ctrl, "| ccode + year"))
  m_a <- tryCatch(feols(f_a, data = data, vcov = ~ccode + year),
                  error = function(e) NULL)

  f_bc <- as.formula(paste("log_churn ~", iv, "+", mediator, "+", ctrl,
                            "| ccode + year"))
  m_bc <- tryCatch(feols(f_bc, data = data, vcov = ~ccode + year),
                   error = function(e) NULL)

  if (is.null(m_total) || is.null(m_a) || is.null(m_bc)) return(NULL)

  ct_total <- coeftable(m_total)
  ct_a     <- coeftable(m_a)
  ct_bc    <- coeftable(m_bc)

  c_total  <- ct_total[iv, "Estimate"]
  a_path   <- ct_a[iv, "Estimate"]
  b_path   <- ct_bc[mediator, "Estimate"]
  c_prime  <- ct_bc[iv, "Estimate"]
  indirect <- a_path * b_path
  pct_mediated <- if (abs(c_total) > 1e-12) indirect / c_total else NA_real_

  cat(sprintf("  %s -> %s -> log_churn: a=%.6f (p=%.3f), b=%.6f (p=%.3f), indirect=%.6f, %%med=%.1f%%\n",
              iv_label, mediator,
              a_path, ct_a[iv, "Pr(>|t|)"],
              b_path, ct_bc[mediator, "Pr(>|t|)"],
              indirect, ifelse(is.na(pct_mediated), NA, pct_mediated * 100)))

  tibble(
    tag = tag,
    iv = iv, iv_label = iv_label,
    mediator = mediator,
    c_total       = c_total,
    c_total_se    = ct_total[iv, "Std. Error"],
    c_total_p     = ct_total[iv, "Pr(>|t|)"],
    a_path        = a_path,
    a_path_se     = ct_a[iv, "Std. Error"],
    a_path_p      = ct_a[iv, "Pr(>|t|)"],
    b_path        = b_path,
    b_path_se     = ct_bc[mediator, "Std. Error"],
    b_path_p      = ct_bc[mediator, "Pr(>|t|)"],
    c_prime       = c_prime,
    c_prime_se    = ct_bc[iv, "Std. Error"],
    c_prime_p     = ct_bc[iv, "Pr(>|t|)"],
    indirect_ab   = indirect,
    pct_mediated  = pct_mediated
  )
}

# IV specifications
iv_specs <- list(
  list(iv = "mid_onset_count",            label = "MID Onset"),
  list(iv = "log_mid_onset_count",        label = "Log MID Onset"),
  list(iv = "fatal_mid_onset_count",      label = "Fatal MID Onset"),
  list(iv = "log_fatal_mid_onset_count",  label = "Log Fatal MID Onset")
)

# ============================================================
# 1. Baron-Kenny: battle deaths mediator, full panel
# ============================================================

cat("\n--- 1. Baron-Kenny: log(1 + Battle Deaths) Mediator (Full Panel) ---\n")

bk_batdeath <- list()
for (ivs in iv_specs) {
  res <- run_bk(med_data, ivs$iv, "log_batdeath", ivs$label, tag_prefix = "batdeath_")
  if (!is.null(res)) bk_batdeath[[res$tag]] <- res
}

# Also: COW war indicator → battle deaths → churn
res <- run_bk(med_data, "in_cow_war", "log_batdeath", "COW War", tag_prefix = "batdeath_")
if (!is.null(res)) bk_batdeath[[res$tag]] <- res

bk_batdeath_df <- bind_rows(bk_batdeath)
cat(sprintf("  Battle death mediation: %d specs tested\n", nrow(bk_batdeath_df)))

# ============================================================
# 2. Imai-Keele-Tingley causal mediation: battle deaths only
# ============================================================

cat("\n--- 2. Causal Mediation Analysis (IKT, battle deaths) ---\n")

# Subset to states with 20+ years for stable FE estimation
state_mid_counts <- med_data |>
  group_by(ccode) |>
  summarise(total_mids = sum(mid_onset_count, na.rm = TRUE), n_yrs = n()) |>
  filter(n_yrs >= 20)

med_sub <- med_data |>
  filter(ccode %in% state_mid_counts$ccode) |>
  filter(!is.na(churn), !is.na(log_batdeath), !is.na(mid_onset_count))

cat(sprintf("  IKT subsample: %d obs (%d states with 20+ years)\n",
            nrow(med_sub), nrow(state_mid_counts)))

ikt_results <- list()

# ΔCINC combos removed; live in 16_robustness_alt_mechanisms.R now.
# Battle-deaths combos retained, with logged DV and full controls
# matching 02 main spec.
ikt_combos <- list(
  list(iv = "mid_onset_count",       med = "log_batdeath",
       label = "MID -> log_batdeath -> log_churn"),
  list(iv = "fatal_mid_onset_count", med = "log_batdeath",
       label = "FatalMID -> log_batdeath -> log_churn"),
  list(iv = "log_mid_onset_count",   med = "log_batdeath",
       label = "LogMID -> log_batdeath -> log_churn"),
  list(iv = "log_fatal_mid_onset_count", med = "log_batdeath",
       label = "LogFatalMID -> log_batdeath -> log_churn")
)

# Routed per-IV inside the loop (mid vs fatal partner-war control).
ikt_ctrl_terms_for_iv <- function(iv) paste(ctrl_for_iv(iv), "+ ccode_f + year_f")

for (combo in ikt_combos) {
  iv  <- combo$iv
  med <- combo$med
  ikt_ctrl_terms <- ikt_ctrl_terms_for_iv(iv)

  d_complete <- med_sub |>
    filter(!is.na(.data[[iv]]), !is.na(log_churn),
           !is.na(.data[[med]])) |>
    droplevels()

  if (nrow(d_complete) < 100) {
    cat(sprintf("  SKIP %s: n=%d < 100\n", combo$label, nrow(d_complete)))
    next
  }

  cat(sprintf("  Running %s (n=%d)...\n", combo$label, nrow(d_complete)))

  f_med <- as.formula(paste(med, "~", iv, "+", ikt_ctrl_terms))
  f_out <- as.formula(paste("log_churn ~", iv, "+", med, "+", ikt_ctrl_terms))

  m_med <- tryCatch(lm(f_med, data = d_complete), error = function(e) NULL)
  m_out <- tryCatch(lm(f_out, data = d_complete), error = function(e) NULL)
  if (is.null(m_med) || is.null(m_out)) {
    cat(sprintf("  FAIL %s: model fitting error\n", combo$label)); next
  }

  med_fit <- tryCatch(
    mediate(m_med, m_out,
            treat = iv, mediator = med,
            boot = TRUE, sims = 1000, boot.ci.type = "perc"),
    error = function(e) { cat(sprintf("  ERROR %s: %s\n", combo$label, e$message)); NULL }
  )

  if (!is.null(med_fit)) {
    ikt_results[[combo$label]] <- tibble(
      iv = iv, mediator = med, label = combo$label,
      acme       = med_fit$d0,
      acme_ci_lo = med_fit$d0.ci[1],
      acme_ci_hi = med_fit$d0.ci[2],
      acme_p     = med_fit$d0.p,
      ade        = med_fit$z0,
      ade_ci_lo  = med_fit$z0.ci[1],
      ade_ci_hi  = med_fit$z0.ci[2],
      ade_p      = med_fit$z0.p,
      total      = med_fit$tau.coef,
      total_p    = med_fit$tau.p,
      prop_med   = med_fit$n0
    )

    cat(sprintf("    ACME=%.6f [%.6f, %.6f] p=%.3f | ADE=%.6f p=%.3f | %%med=%.1f%%\n",
                med_fit$d0, med_fit$d0.ci[1], med_fit$d0.ci[2], med_fit$d0.p,
                med_fit$z0, med_fit$z0.p,
                ifelse(is.na(med_fit$n0), NA, med_fit$n0 * 100)))
  }
}

ikt_df <- bind_rows(ikt_results)

# ============================================================
# 3. Sensitivity analysis (battle deaths only)
# ============================================================

cat("\n--- 3. Sensitivity Analysis (battle deaths) ---\n")

sens_results <- list()

d_sens_bd <- med_sub |>
  filter(!is.na(fatal_mid_onset_count), !is.na(log_churn),
         !is.na(log_batdeath)) |>
  droplevels()

if (nrow(d_sens_bd) >= 100) {
  ikt_ctrl_terms_bd <- ikt_ctrl_terms_for_iv("fatal_mid_onset_count")
  f_med_bd <- as.formula(paste("log_batdeath ~ fatal_mid_onset_count +",
                                ikt_ctrl_terms_bd))
  f_out_bd <- as.formula(paste("log_churn ~ fatal_mid_onset_count + log_batdeath +",
                                ikt_ctrl_terms_bd))

  m_med_bd <- tryCatch(lm(f_med_bd, data = d_sens_bd), error = function(e) NULL)
  m_out_bd <- tryCatch(lm(f_out_bd, data = d_sens_bd), error = function(e) NULL)

  if (!is.null(m_med_bd) && !is.null(m_out_bd)) {
    med_fit_bd <- tryCatch(
      mediate(m_med_bd, m_out_bd,
              treat = "fatal_mid_onset_count", mediator = "log_batdeath",
              boot = TRUE, sims = 1000, boot.ci.type = "perc"),
      error = function(e) NULL
    )

    if (!is.null(med_fit_bd)) {
      sens_fit_bd <- tryCatch(
        medsens(med_fit_bd, rho.by = 0.05, effect.type = "indirect", sims = 1000),
        error = function(e) { cat(sprintf("  Sensitivity error (BD): %s\n", e$message)); NULL }
      )

      if (!is.null(sens_fit_bd)) {
        sens_results$fatalmid_batdeath_churn <- list(
          rho_at_zero    = sens_fit_bd$rho.at.sign.change,
          r2_med_at_zero = sens_fit_bd$R2tilde.sigma.at.sign.change,
          r2_out_at_zero = sens_fit_bd$R2tilde.tau.at.sign.change
        )

        cat(sprintf("  FatalMID -> log_batdeath -> log_churn: rho at ACME=0: %.3f\n",
                    ifelse(length(sens_fit_bd$rho.at.sign.change) > 0,
                           sens_fit_bd$rho.at.sign.change, NA)))

        png(file.path(FIG_DIR, "fig_mediation_sensitivity_batdeath.png"),
            width = 8, height = 6, units = "in", res = 300)
        plot(sens_fit_bd, main = "Sensitivity: Fatal MID -> Battle Deaths -> log(Churn)",
             xlab = expression(rho), ylab = "ACME")
        dev.off()
        cat("  Saved sensitivity plot (battle deaths)\n")
      }
    }
  }
}

# ============================================================
# 4. Save results + LaTeX tables
# ============================================================

cat("\n--- 4. Saving Results ---\n")

mediation_results <- list(
  batdeath_panel = bk_batdeath_df,
  ikt_mediation  = ikt_df,
  sensitivity    = sens_results,
  panel_info = tibble(
    n_obs    = nrow(med_data),
    n_states = n_distinct(med_data$ccode),
    n_years  = n_distinct(med_data$year),
    year_min = min(med_data$year),
    year_max = max(med_data$year)
  )
)

saveRDS(mediation_results, file.path(OUTPUT_DIR, "mediation_results.rds"))
cat("  Saved mediation_results.rds\n")

APPX_DIR <- here::here("results", "appendix")
dir.create(APPX_DIR, recursive = TRUE, showWarnings = FALSE)

# CSV exports (battle-deaths only; ΔCINC moved to script 18)
if (nrow(bk_batdeath_df) > 0) {
  write_csv(bk_batdeath_df, file.path(OUTPUT_DIR, "mediation_batdeath_panel.csv"))
}
if (nrow(ikt_df) > 0) {
  write_csv(ikt_df, file.path(OUTPUT_DIR, "mediation_ikt_batdeath.csv"))
}

# LaTeX: Baron-Kenny battle-deaths table
if (nrow(bk_batdeath_df) > 0) {
  bk_disp <- bk_batdeath_df |>
    mutate(
      sig_a  = ifelse(a_path_p < 0.05, "*", ""),
      sig_b  = ifelse(b_path_p < 0.05, "*", ""),
      sig_c  = ifelse(c_total_p < 0.05, "*", ""),
      sig_cp = ifelse(c_prime_p < 0.05, "*", "")
    )

  tex_lines <- c(
    "\\begin{table}[htbp]",
    "\\centering",
    "\\caption{Battle Deaths as Informativeness Mediator: Conflict $\\to$ Battle Deaths $\\to$ $\\log(1+\\text{Churn})$.}",
    "\\label{tab:mediation_batdeath}",
    "\\footnotesize",
    "\\begin{tabular}{lrrrrrr}",
    "\\toprule",
    "IV & $c$ & $a$ (X$\\to$M) & $b$ (M$\\to$Y) & $c'$ & $a \\times b$ & \\% Med. \\\\",
    "\\midrule"
  )

  for (i in seq_len(nrow(bk_disp))) {
    r <- bk_disp[i, ]
    tex_lines <- c(tex_lines, sprintf(
      "%s & %.4f%s & %.6f%s & %.4f%s & %.4f%s & %.6f & %.1f \\\\",
      r$iv_label,
      r$c_total, r$sig_c,
      r$a_path,  r$sig_a,
      r$b_path,  r$sig_b,
      r$c_prime, r$sig_cp,
      r$indirect_ab,
      ifelse(is.na(r$pct_mediated), 0, r$pct_mediated * 100)
    ))
  }

  tex_lines <- c(tex_lines,
    "\\bottomrule",
    "\\multicolumn{7}{l}{\\footnotesize $^*p<0.05$. Two-way FE (state + year), two-way clustered SEs.} \\\\",
    "\\multicolumn{7}{l}{\\footnotesize DV: $\\log(1+\\text{Churn})$. Mediator: $\\log(1+\\text{Battle Deaths})$.} \\\\",
    "\\multicolumn{7}{l}{\\footnotesize Controls: \\texttt{partner\\_war\\_frac\\_weighted}, democracy, autocracy, CINC, major power.} \\\\",
    "\\end{tabular}",
    "\\end{table}"
  )

  writeLines(tex_lines, file.path(APPX_DIR, "tableA_mediation_batdeath.tex"))
  cat("  LaTeX -> tableA_mediation_batdeath.tex\n")
}

# LaTeX: IKT causal mediation (battle deaths only)
if (nrow(ikt_df) > 0) {
  tex_ikt <- c(
    "\\begin{table}[htbp]",
    "\\centering",
    "\\caption{Causal Mediation (Imai-Keele-Tingley): Battle Deaths as Informativeness Channel.}",
    "\\label{tab:mediation_ikt_batdeath}",
    "\\small",
    "\\begin{tabular}{lrrrrrr}",
    "\\toprule",
    "Specification & ACME & 95\\% CI & $p$ & ADE & $p$ & \\% Med. \\\\",
    "\\midrule"
  )

  for (i in seq_len(nrow(ikt_df))) {
    r <- ikt_df[i, ]
    tex_ikt <- c(tex_ikt, sprintf(
      "%s & %.5f & [%.5f, %.5f] & %.3f & %.5f & %.3f & %.1f \\\\",
      gsub("->", "$\\\\to$", r$label),
      r$acme, r$acme_ci_lo, r$acme_ci_hi, r$acme_p,
      r$ade, r$ade_p,
      ifelse(is.na(r$prop_med), 0, r$prop_med * 100)
    ))
  }

  tex_ikt <- c(tex_ikt,
    "\\bottomrule",
    "\\multicolumn{7}{l}{\\footnotesize Bootstrap (1000 draws). State + year FE.} \\\\",
    "\\multicolumn{7}{l}{\\footnotesize DV: $\\log(1+\\text{Churn})$. Mediator: $\\log(1+\\text{Battle Deaths})$.} \\\\",
    "\\end{tabular}",
    "\\end{table}"
  )

  writeLines(tex_ikt, file.path(APPX_DIR, "tableA_mediation_ikt_batdeath.tex"))
  cat("  LaTeX -> tableA_mediation_ikt_batdeath.tex\n")
}

cat("\nDone. Battle-deaths informativeness analysis complete.\n")
cat("Note: ΔCINC mediation lives in 16_robustness_alt_mechanisms.R.\n")
