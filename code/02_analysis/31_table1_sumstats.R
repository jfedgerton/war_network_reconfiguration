# 31_table1_sumstats.R -- main-text Table 1 (summary statistics, systemic and state level).
# [2026-10-27] Adds the peace-dyad DVs (js_divergence_peace_ft1, churn_peace_ft1: dyads in an ongoing fatal MID in t-1 dropped)
#              next to the all-dyad DVs, matching the two panels of the main-text tables.
# Input:  data/processed/system_analysis_data.csv (written by 01), data/processed/state_churn_data.csv (written by 02).
# Output: results/main/table1_sumstats.tex (\label{tab:sumstats_system}), results/main/table1_sumstats.csv.
# Run after 01 and 02. No estimation; reads the processed data only.
suppressPackageStartupMessages({library(dplyr); library(readr); library(tidyr)})
set.seed(123)
IN  <- here::here("data", "processed")
OUT <- here::here("results", "main")

sys <- read_csv(file.path(IN, "system_analysis_data.csv"), show_col_types = FALSE)
st  <- read_csv(file.path(IN, "state_churn_data.csv"), show_col_types = FALSE)

sys_vars <- c(js_divergence = "JS Divergence", js_divergence_peace_ft1 = "JS Divergence, peace dyads",
              fatal_mid_count = "Fatal MID Count",
              n_democracies = "Count of Democracies", n_autocracies = "Count of Autocracies", n_states = "Count of States")
st_vars  <- c(churn = "Network Churn", churn_peace_ft1 = "Network Churn, peace dyads",
              fatal_onset_lag1 = "Fatal MID Onset Count$_{t-1}$",
              democracy = "Democracy", autocracy = "Autocracy", cinc = "CINC", major_power = "Major Power",
              is_isolated_t1 = "Isolated at $t-1$")
stopifnot(all(names(sys_vars) %in% names(sys)), all(names(st_vars) %in% names(st)))

summ <- function(df, vars) {
  df |> select(all_of(names(vars))) |>
    pivot_longer(everything(), names_to = "variable", values_to = "value") |>
    group_by(variable) |>
    summarise(N = sum(!is.na(value)), Mean = mean(value, na.rm = TRUE), SD = sd(value, na.rm = TRUE),
              Min = min(value, na.rm = TRUE), Q25 = quantile(value, 0.25, na.rm = TRUE),
              Median = median(value, na.rm = TRUE), Q75 = quantile(value, 0.75, na.rm = TRUE),
              Max = max(value, na.rm = TRUE), .groups = "drop") |>
    mutate(Label = vars[variable]) |>
    arrange(match(variable, names(vars)))
}
s_sys <- summ(sys, sys_vars)
s_st  <- summ(st, st_vars)

# ---- checks against the previous Table 1 (all-dyad rows must be unchanged) ----
stopifnot(s_sys$N[s_sys$variable == "js_divergence"] == 198,
          s_sys$N[s_sys$variable == "js_divergence_peace_ft1"] == 198,
          s_st$N[s_st$variable == "churn"] == 15318,
          abs(s_sys$Mean[s_sys$variable == "fatal_mid_count"] - 9.84) < 0.005,
          abs(s_st$Mean[s_st$variable == "churn"] - 0.074) < 0.0005)
# peace-dyad DVs can only be weakly below the all-dyad counts in N
stopifnot(s_st$N[s_st$variable == "churn_peace_ft1"] <= s_st$N[s_st$variable == "churn"])

write_csv(bind_rows(mutate(s_sys, level = "system"), mutate(s_st, level = "state")), file.path(OUT, "table1_sumstats.csv"))

# ---- LaTeX in the main-text format: 0 for exact zero, <0.01 for small positives, else two decimals ----
f <- function(x) ifelse(x == 0, "0", ifelse(abs(x) < 0.005, "$<$0.01", sprintf("%.2f", x)))
row_tex <- function(r) sprintf("\\hspace{1pt} %s & %d & %s & %s & %s & %s & %s & %s & %s \\\\",
                               r$Label, r$N, f(r$Mean), f(r$SD), f(r$Min), f(r$Q25), f(r$Median), f(r$Q75), f(r$Max))
rows <- function(s, v) vapply(which(s$variable %in% v), function(i) row_tex(s[i, ]), character(1))

tex <- c("\\begin{table}[t]", "\\scriptsize", "\\centering",
  "\\caption{Summary statistics for systemic and state level variables}", "\\label{tab:sumstats_system}",
  "\\begin{tabular}{lrrrrrrrr}", "\\toprule",
  "Name & Obs & Mean & Std. Dev. & Min. & 25$^\\text{th}$ Quant. & Median & 75$^\\text{th}$ Quant. & Max \\\\",
  "\\midrule", "\\textbf{Systemic variables} \\\\", "\\addlinespace",
  "\\textit{Dependent variables} \\\\", rows(s_sys, c("js_divergence", "js_divergence_peace_ft1")), "\\addlinespace",
  "\\textit{Independent variables} \\\\", rows(s_sys, "fatal_mid_count"), "\\addlinespace",
  "\\textit{Control variables} \\\\", rows(s_sys, c("n_democracies", "n_autocracies", "n_states")),
  "\\midrule", "\\textbf{State Level variables} \\\\", "\\addlinespace",
  "\\textit{Dependent variables} \\\\", rows(s_st, c("churn", "churn_peace_ft1")), "\\addlinespace",
  "\\textit{Independent variables} \\\\", rows(s_st, "fatal_onset_lag1"), "\\addlinespace",
  "\\textit{Control variables} \\\\", rows(s_st, c("democracy", "autocracy", "cinc", "major_power", "is_isolated_t1")),
  "\\bottomrule", "\\end{tabular}", "\\end{table}")
writeLines(tex, file.path(OUT, "table1_sumstats.tex"))
cat("wrote", file.path(OUT, "table1_sumstats.tex"), "\n")
