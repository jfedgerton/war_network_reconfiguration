# [2026-09-26] Headline system-level DV switched from KL to Jensen-Shannon (JS)
# divergence. KL results are retained in results/appendix/kl_canonical/.
###############################################################################
# 26_robustness_igo_cadence.R
# Appendix robustness: IGO measurement cadence and layer responsiveness
#
# Motivation:
#   COW IGO co-membership is recorded at 5-year intervals (1815-1965) and
#   annually from 1966 onward, with values carried forward between updates.
#   A reader may worry that the IGO layer's low responsiveness is an artifact
#   of the recording schedule. This script shows it is structural.
#
# Produces (written to results/appendix/):
#   tab_igo_layers.tex          -- per-layer turnover / KL / active edges
#   tab_igo_loo.tex             -- main results WITH vs WITHOUT the IGO layer
#   tab_igo_annual_altspec.tex  -- system & state models re-estimated on the
#                                  1966+ annual-IGO subsample
# and prints the cadence decomposition (pre/post 1966, grid vs off-grid) plus
# the annual-era per-layer turnover to the console.
#
# Inputs (all already produced by the pipeline):
#   data/processed/layer_specific_measures.rds
#   data/processed/robustness_specifications.csv
#   data/processed/system_analysis_data.rds
#   data/processed/state_churn_data.csv
###############################################################################

library(tidyverse)
library(fixest)

OUTPUT_DIR <- here::here("data", "processed")
TAB_DIR    <- here::here("results", "appendix")
dir.create(TAB_DIR, recursive = TRUE, showWarnings = FALSE)

NW_LAG <- 5L

# Small formatting helpers ----------------------------------------------------
fmt  <- function(x, d = 3) formatC(x, format = "f", digits = d)
icom <- function(x) formatC(round(x), format = "d", big.mark = ",")
stars <- function(p) {
  ifelse(p < 0.001, "$^{***}$",
  ifelse(p < 0.01,  "$^{**}$",
  ifelse(p < 0.05,  "$^{*}$",
  ifelse(p < 0.10,  "$^{\\dagger}$", ""))))
}

# ===========================================================================
# 1. Layer responsiveness table
# ===========================================================================
cat("\n== Layer responsiveness ==\n")
layer_meas <- readRDS(file.path(OUTPUT_DIR, "layer_specific_measures.rds"))

layer_labels <- c(
  layer_offensive = "Offensive alliance",
  layer_defensive = "Defensive alliance",
  layer_igo       = "IGO",
  layer_trade     = "Trade",
  layer_dca       = "DCA"
)

layer_tab <- layer_meas |>
  filter(!is.na(edge_turnover_pct)) |>
  filter(layer != "layer_dca") |>   # DCA excluded outside the contamination analyses [2026-07-15]
  group_by(layer) |>
  summarise(turnover = mean(edge_turnover_pct, na.rm = TRUE),
            kl       = mean(js_divergence,     na.rm = TRUE),
            edges    = mean(n_edges,           na.rm = TRUE),
            .groups  = "drop") |>
  mutate(Layer = layer_labels[layer]) |>
  arrange(desc(turnover))

print(layer_tab)

layer_rows <- layer_tab |>
  transmute(line = sprintf("%s & %s & %s & %s \\\\",
                           Layer, fmt(turnover, 3), fmt(kl, 3), icom(edges)))

layer_tex <- c(
  "\\begin{table}[t]",
  "\\centering",
  "\\caption{Layer responsiveness. Mean annual edge turnover, JS divergence,",
  "and active-edge count by network layer. IGO is the densest and least",
  "volatile layer.}",
  "\\label{tab:igo_layers}",
  "\\begin{tabular}{lccc}",
  "\\toprule",
  "Layer & Edge turnover & JS divergence & Active edges \\\\",
  "\\midrule",
  layer_rows$line,
  "\\bottomrule",
  "\\end{tabular}",
  "\\end{table}"
)
writeLines(layer_tex, file.path(TAB_DIR, "tab_igo_layers.tex"))
cat("  -> tab_igo_layers.tex\n")

# ===========================================================================
# 2. Cadence decomposition (console + drives the annual-era note)
# ===========================================================================
cat("\n== IGO cadence decomposition ==\n")
igo <- layer_meas |>
  filter(layer == "layer_igo", !is.na(edge_turnover_pct)) |>
  mutate(era  = if_else(year <= 1965, "1816-1965 (5-yr COW)",
                        "1966-2014 (annual COW)"),
         grid = (year %% 5 == 0))

cad <- igo |>
  group_by(era, grid) |>
  summarise(turnover = mean(edge_turnover_pct), .groups = "drop")
print(cad)

# Annual-era per-layer turnover (the headline robustness number)
annual_turn <- layer_meas |>
  filter(year >= 1966, !is.na(edge_turnover_pct)) |>
  group_by(layer) |>
  summarise(turnover = mean(edge_turnover_pct), .groups = "drop") |>
  mutate(Layer = layer_labels[layer]) |>
  arrange(desc(turnover))
cat("\nAnnual-era (1966+) mean turnover by layer:\n")
print(annual_turn)

# ===========================================================================
# 3. Leave-IGO-out across the specification multiverse
# ===========================================================================
cat("\n== Main results WITH vs WITHOUT the IGO layer ==\n")
specs <- readr::read_csv(file.path(OUTPUT_DIR, "robustness_specifications.csv"),
                         show_col_types = FALSE) |>
  filter(!is.na(coef), !is.na(pval))

loo <- specs |>
  mutate(igo = if_else(has_igo, "Included", "Excluded")) |>
  group_by(level, igo) |>
  summarise(pct_pos     = mean(coef > 0) * 100,
            pct_pos_sig = mean(coef > 0 & pval < 0.05) * 100,
            mean_coef   = mean(coef),
            .groups = "drop") |>
  mutate(level = factor(level, levels = c("system", "state"),
                        labels = c("System", "State")),
         igo   = factor(igo, levels = c("Excluded", "Included"))) |>
  arrange(level, igo)
print(loo)

loo_rows <- loo |>
  group_by(level) |>
  mutate(lvl_lab = if_else(row_number() == 1,
                           sprintf("\\multirow{2}{*}{%s}", level), "")) |>
  ungroup() |>
  transmute(line = sprintf("%s & %s & %s & %s & %s \\\\",
                           lvl_lab, igo,
                           fmt(pct_pos, 1), fmt(pct_pos_sig, 1),
                           fmt(mean_coef, 2)))
# add a small vertical gap after the System block
loo_lines <- loo_rows$line
loo_lines[2] <- paste0(loo_lines[2], "[2pt]")

loo_tex <- c(
  "\\begin{table}[t]",
  "\\centering",
  "\\caption{Main results with and without the IGO layer, across the",
  "specification multiverse. Excluding IGO leaves the share of positive (and",
  "positive-and-significant) coefficients essentially unchanged.}",
  "\\label{tab:igo_loo}",
  "\\begin{tabular}{llccc}",
  "\\toprule",
  "Level & IGO layer & \\% positive & \\% pos.\\ \\& sig. & Mean coef. \\\\",
  "\\midrule",
  loo_lines,
  "\\bottomrule",
  "\\end{tabular}",
  "\\end{table}"
)
writeLines(loo_tex, file.path(TAB_DIR, "tab_igo_loo.tex"))
cat("  -> tab_igo_loo.tex\n")

# ===========================================================================
# 4. Alternative specification: re-estimate the MAIN models on 1966+ only
#    (annual-IGO era), for both the systemic and state-level models.
# ===========================================================================
cat("\n== 1966+ annual-IGO alternative specification ==\n")

## --- System: main spec (full controls + lagged DV), feols + NW(5) --------
## NB: the .rds is the pre-transform input; the log/lag columns the main
## spec uses are written by 01_system_level_analysis.R to the .csv.
sys <- readr::read_csv(file.path(OUTPUT_DIR, "system_analysis_data.csv"),
                       show_col_types = FALSE) |>
  filter(!is.na(log_js_divergence_lag1), !is.na(log_mid_count_lag1))

sys_fml <- log_js_divergence ~ log_mid_count_lag1 +
  log1p(n_democracies) + log1p(n_autocracies) + log1p(n_states) +
  poly(year, 2) + log_js_divergence_lag1

fit_sys <- function(d) feols(sys_fml, data = d,
                             panel.id = ~ unit + year, vcov = NW(lag = NW_LAG))

sys_full <- fit_sys(sys)
sys_66   <- fit_sys(filter(sys, year >= 1966))

## --- State: main TWFE spec, two-way clustered ----------------------------
st <- readr::read_csv(file.path(OUTPUT_DIR, "state_churn_data.csv"),
                      show_col_types = FALSE)

st_fml <- log_churn ~ log_mid_onset_count + partner_war_frac_weighted +
  democracy + autocracy + cinc + major_power + is_isolated_t1 |
  ccode + year

fit_st <- function(d) feols(st_fml, data = d, vcov = ~ ccode + year)

st_full <- fit_st(st)
st_66   <- fit_st(filter(st, year >= 1966))

## --- Pull one coefficient row from a feols fit ---------------------------
grab <- function(fit, term) {
  ct <- coeftable(fit)
  if (!term %in% rownames(ct)) return(c(b = NA, se = NA, p = NA, n = NA))
  c(b = ct[term, "Estimate"], se = ct[term, "Std. Error"],
    p = ct[term, "Pr(>|t|)"], n = nobs(fit))
}

cell <- function(v) {
  if (is.na(v["b"])) return("--")
  sprintf("%s%s (%s)", fmt(v["b"], 3), stars(v["p"]), fmt(v["se"], 3))
}

alt_specs <- list(
  list(label = "System: MID count$_{t-1}$",
       full = grab(sys_full, "log_mid_count_lag1"),
       yr66 = grab(sys_66,   "log_mid_count_lag1")),
  list(label = "State: MID onset",
       full = grab(st_full,  "log_mid_onset_count"),
       yr66 = grab(st_66,    "log_mid_onset_count")),
  list(label = "State: Partner-war exposure",
       full = grab(st_full,  "partner_war_frac_weighted"),
       yr66 = grab(st_66,    "partner_war_frac_weighted"))
)

cat("\nAlt-spec coefficients (b, se, p, n):\n")
for (r in alt_specs) {
  cat(sprintf("  %-30s full: b=%.4f p=%.4g n=%d | 1966+: b=%.4f p=%.4g n=%d\n",
              r$label,
              r$full["b"], r$full["p"], r$full["n"],
              r$yr66["b"], r$yr66["p"], r$yr66["n"]))
}

alt_rows <- vapply(alt_specs, function(r) {
  sprintf("%s & %s & %s & %s & %s \\\\",
          r$label,
          cell(r$full), icom(r$full["n"]),
          cell(r$yr66), icom(r$yr66["n"]))
}, character(1))

alt_tex <- c(
  "\\begin{table}[t]",
  "\\centering",
  "\\caption{Main systemic and state-level specifications re-estimated on the",
  "1966+ subsample, where COW records IGO co-membership annually (no five-year",
  "carry-forward). Estimates and significance are consistent with the full",
  "sample. Standard errors in parentheses (system: Newey--West HAC, lag 5;",
  "state: two-way clustered on state and year).}",
  "\\label{tab:igo_annual_altspec}",
  "\\begin{tabular}{lcccc}",
  "\\toprule",
  " & \\multicolumn{2}{c}{Full sample} & \\multicolumn{2}{c}{1966--2014} \\\\",
  "\\cmidrule(lr){2-3} \\cmidrule(lr){4-5}",
  "Coefficient & Estimate & $N$ & Estimate & $N$ \\\\",
  "\\midrule",
  alt_rows,
  "\\bottomrule",
  "\\end{tabular}",
  "\\vspace{2pt}",
  paste0("{\\footnotesize Significance: $^{***}p<0.001$; $^{**}p<0.01$; ",
         "$^{*}p<0.05$; $^{\\dagger}p<0.1$.}"),
  "\\end{table}"
)
writeLines(alt_tex, file.path(TAB_DIR, "tab_igo_annual_altspec.tex"))
cat("  -> tab_igo_annual_altspec.tex\n")

cat("\nDone: 26_robustness_igo_cadence.R\n")
