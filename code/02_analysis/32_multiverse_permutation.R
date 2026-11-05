###############################################################################
# 32_multiverse_permutation.R   [2026-11-05]
#
# Specification-curve inference for the multiverse (Simonsohn, Simmons and Nelson 2020).
# Question: how often would a multiverse built on data with NO conflict-reconfiguration link produce support
# shares as high as the observed ones? Specifications are highly correlated, so the 10% chance benchmark has
# no sampling distribution of its own; this script builds one.
#
# Design
#   * Each draw reruns 07_robustness_multiverse.R in permutation mode (env MV_PERM_DRAW = k):
#       - the same fixed random subset of layer combos in every draw (MV_PERM_N_COMBOS, seed 123),
#         crossed with every weighting, transform, IV, DV, model class and extra control (none, log GDP);
#       - draw 0 = observed data (no permutation) on that subset;
#       - draw k > 0 = conflict IVs permuted with seed 123 + k. System: years shuffled jointly across the four
#         system IVs. State: t-1 IVs shuffled jointly across years within each state (as in 04's permutation
#         placebo). Networks, controls and samples are untouched.
#   * Per draw and hypothesis (H1 system, H2 state) it records the share of specs significant-positive at
#     one-sided p < 0.10 (the heatmap rule), the share positive, and the median and mean z.
#   * p-value = (1 + #draws with statistic >= observed) / (1 + #draws).
#   * Check: the observed draw must reproduce the main 07 run exactly on the same subset.
#
# Draws run as separate R processes (each sequential, ~5-8 GB RAM). Resumable: finished draws are skipped.
#   MV_PERM_DRAWS     number of permuted draws            (default 500)
#   MV_PERM_PROCS     concurrent processes                (default 10; ~8 GB each)
#   MV_PERM_N_COMBOS  layer combos in the fixed subset    (default 3 of 12)
#
# Inputs : 07 and its inputs; data/processed/robustness_specifications.csv (main 07 run, for the check)
# Outputs: data/processed/mv_perm/draw_XXXX.csv (one per draw), logs/mv_perm/draw_XXXX.log
#          results/appendix/multiverse_permutation_summary.csv
#          code/03_figures/output/fig_multiverse_permutation.png
###############################################################################
suppressPackageStartupMessages({library(dplyr); library(readr); library(tidyr); library(ggplot2)})
setwd(here::here())
set.seed(123)

N_DRAWS  <- as.integer(Sys.getenv("MV_PERM_DRAWS", "500"))
N_PROCS  <- as.integer(Sys.getenv("MV_PERM_PROCS", "10"))
N_COMBOS <- as.integer(Sys.getenv("MV_PERM_N_COMBOS", "3"))
stopifnot(N_DRAWS >= 1, N_PROCS >= 1, N_COMBOS >= 1)

SCRIPT_07 <- here::here("code", "02_analysis", "07_robustness_multiverse.R")
PERM_DIR  <- here::here("data", "processed", "mv_perm")
LOG_DIR   <- here::here("logs", "mv_perm")
FIG_DIR   <- here::here("code", "03_figures", "output")
APPX_DIR  <- here::here("results", "appendix")
for (d in c(PERM_DIR, LOG_DIR, FIG_DIR, APPX_DIR)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
RS <- file.path(R.home("bin"), "Rscript")

# A subset change invalidates earlier draws: record the setting and refuse to mix.
cfg_file <- file.path(PERM_DIR, "config.txt")
cfg_now  <- sprintf("n_combos=%d", N_COMBOS)
if (file.exists(cfg_file) && readLines(cfg_file)[1] != cfg_now)
  stop("mv_perm/ holds draws from a different configuration (", readLines(cfg_file)[1], "); move or delete them first.")
writeLines(cfg_now, cfg_file)

run_draw <- function(k) {
  out <- file.path(PERM_DIR, sprintf("draw_%04d.csv", k))
  if (file.exists(out)) return("skipped")
  logf <- file.path(LOG_DIR, sprintf("draw_%04d.log", k))
  rc <- system2(RS, shQuote(SCRIPT_07), stdout = logf, stderr = logf,
                env = c(sprintf("MV_PERM_DRAW=%d", k), sprintf("MV_PERM_N_COMBOS=%d", N_COMBOS),
                        "ROBUSTNESS_WORKERS=1", "ROBUSTNESS_RESET=0"))
  if (rc != 0 || !file.exists(out)) return(sprintf("FAILED (rc=%s); see %s", rc, logf))
  "ok"
}

# Observed draw first (it is also the timing check), then the permuted draws in parallel.
t0 <- Sys.time()
st0 <- run_draw(0L)
cat(sprintf("draw 0: %s (%.1f min)\n", st0, as.numeric(difftime(Sys.time(), t0, units = "mins"))))
stopifnot(st0 %in% c("ok", "skipped"))

t1 <- Sys.time()
status <- unlist(parallel::mclapply(seq_len(N_DRAWS), run_draw, mc.cores = N_PROCS, mc.preschedule = FALSE))
cat(sprintf("permuted draws: %s (%.1f h)\n", paste(names(table(status)), table(status), collapse = ", "),
            as.numeric(difftime(Sys.time(), t1, units = "hours"))))
if (any(grepl("^FAILED", status))) stop(sum(grepl("^FAILED", status)), " draws failed: ",
                                        paste(head(status[grepl("^FAILED", status)], 5), collapse = " | "))

# ---- check: observed draw reproduces the main 07 run on the same subset ------------------------------------
KEYS <- c("level", "dv", "war_variable", "weighting", "transform", "layers", "model", "extra_control")
obs_specs  <- read_csv(file.path(PERM_DIR, "draw_0000_specs.csv"), show_col_types = FALSE)
main_specs <- read_csv(here::here("data", "processed", "robustness_specifications.csv"), show_col_types = FALSE) |>
  mutate(extra_control = coalesce(extra_control, "none")) |>
  filter(!startsWith(model, "VAR_")) |>
  select(all_of(KEYS), coef_main = coef, pval_main = pval)
stopifnot(!anyDuplicated(obs_specs[, KEYS]), !anyDuplicated(main_specs[, KEYS]))
chk <- inner_join(obs_specs, main_specs, by = KEYS)
cat(sprintf("observed-draw check: %d of %d specs matched to the main run; max |coef diff| = %.2e; max |p diff| = %.2e\n",
            nrow(chk), nrow(obs_specs), max(abs(chk$coef - chk$coef_main)), max(abs(chk$pval - chk$pval_main))))
stopifnot(nrow(chk) == nrow(obs_specs),
          max(abs(chk$coef - chk$coef_main)) < 1e-8, max(abs(chk$pval - chk$pval_main)) < 1e-8)

# ---- inference ------------------------------------------------------------------------------------------------
draws <- list.files(PERM_DIR, pattern = "^draw_[0-9]{4}\\.csv$", full.names = TRUE) |>
  lapply(read_csv, show_col_types = FALSE) |> bind_rows()
obs  <- draws |> filter(draw == 0)
null <- draws |> filter(draw > 0, draw <= N_DRAWS)
stopifnot(nrow(obs) == 2, n_distinct(null$draw) == N_DRAWS, all(table(null$hypothesis) == N_DRAWS))
# Spec counts can differ slightly across draws only if a permuted IV makes a fit fail; report, do not stop.
n_range <- draws |> group_by(hypothesis) |> summarise(n_min = min(n), n_max = max(n), .groups = "drop")
print(n_range)
if (any(n_range$n_max - n_range$n_min > 0.01 * n_range$n_max)) warning("spec counts vary by more than 1% across draws")

STATS <- c(share_sig_pos = "Share significant-positive", median_z = "Median z")
summ <- lapply(names(STATS), function(st) {
  obs |> select(hypothesis, n, observed = all_of(st)) |>
    left_join(null |> group_by(hypothesis) |>
                summarise(null_mean = mean(.data[[st]]), null_p95 = quantile(.data[[st]], 0.95),
                          null_max = max(.data[[st]]), .groups = "drop"), by = "hypothesis") |>
    left_join(null |> select(hypothesis, v = all_of(st)) |>
                left_join(obs |> select(hypothesis, o = all_of(st)), by = "hypothesis") |>
                group_by(hypothesis) |>
                summarise(p_value = (1 + sum(v >= o)) / (1 + n()), .groups = "drop"), by = "hypothesis") |>
    mutate(statistic = st, n_draws = N_DRAWS)
}) |> bind_rows() |> select(statistic, hypothesis, n_specs = n, observed, null_mean, null_p95, null_max, p_value, n_draws)
stopifnot(all(summ$p_value > 0), all(summ$p_value <= 1))
write_csv(summ, file.path(APPX_DIR, "multiverse_permutation_summary.csv"))
print(summ)

# ---- figure: null distributions with the observed value; p-values sit in the panel strips ------------------
lab <- summ |> mutate(panel = sprintf("%s, %s\n(observed %.2f; p %s)", hypothesis, STATS[statistic], observed,
                                      ifelse(p_value <= 1 / (1 + N_DRAWS), sprintf("< %.3f", 1 / N_DRAWS),
                                             sprintf("= %.3f", p_value))))
plot_null <- null |> pivot_longer(all_of(names(STATS)), names_to = "statistic", values_to = "value") |>
  left_join(lab |> select(statistic, hypothesis, panel), by = c("statistic", "hypothesis"))
plot_obs  <- lab |> select(panel, observed)
p <- ggplot(plot_null, aes(value)) +
  geom_histogram(bins = 40, fill = "grey70", colour = "grey35") +
  geom_vline(data = plot_obs, aes(xintercept = observed), colour = "#2C6CB0", linewidth = 1.1) +
  facet_wrap(~ panel, scales = "free", ncol = 2) +
  labs(title = "Multiverse support under permuted conflict",
       subtitle = sprintf("Grey: %d permuted draws. Blue: observed. Same %d layer combos in every draw, crossed with all weightings, transforms, IVs, DVs, models and controls.",
                          N_DRAWS, N_COMBOS),
       x = NULL, y = "Draws") +
  theme_minimal(base_size = 14) +
  theme(panel.grid.minor = element_blank(), strip.text = element_text(face = "bold", colour = "black"),
        axis.text = element_text(colour = "black"), plot.title = element_text(face = "bold"),
        plot.subtitle = element_text(size = 10, colour = "grey25"))
ggsave(file.path(FIG_DIR, "fig_multiverse_permutation.png"), p, width = 12, height = 8, dpi = 300, bg = "white")
cat("fig_multiverse_permutation.png written\n")
