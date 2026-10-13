###############################################################################
# 04_bayesian_structural_ts.R
# Bayesian Structural Time Series (CausalImpact)
#
# Counterfactual: what would network structure look like absent major wars?
# Uses pre-war period to fit a structural time series model, then compares
# observed post-war network divergence to the synthetic counterfactual.
#
# Mirrors the 2 transforms x 2 controls x 2 IVs design from
# 01_system_level_analysis.R, adapted to BSTS / intervention analysis:
#
#   transform = c("logged", "raw")   # log1p applied to DV (and covariates) or not
#   controls  = c("none", "full")    # full = log1p(n_democracies) +
#                                      log1p(n_autocracies) + log1p(n_states)
#
# DVs: KL divergence, JS divergence, edge turnover %, degree KL,
#      degree Wasserstein
###############################################################################

library(tidyverse)
library(CausalImpact)
select <- dplyr::select

OUTPUT_DIR <- here::here("data", "processed")
FIG_DIR <- here::here("code", "03_figures", "output")
TABLE_DIR <- here::here("results", "main")
APPX_DIR  <- here::here("results", "appendix")
dir.create(FIG_DIR,   recursive = TRUE, showWarnings = FALSE)
dir.create(TABLE_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(APPX_DIR,  recursive = TRUE, showWarnings = FALSE)

# --------------------------------------------------------------------------- #
# 0. Load data
# --------------------------------------------------------------------------- #

system_data <- readRDS(file.path(OUTPUT_DIR, "system_analysis_data.rds")) |>
  filter(!is.na(kl_divergence))

# --------------------------------------------------------------------------- #
# 1. CausalImpact: war x DV x transform x controls
# --------------------------------------------------------------------------- #

WARS <- list(
  WWI          = list(pre = c(1816, 1913), post = c(1914, 1920)),
  WWII         = list(pre = c(1920, 1938), post = c(1939, 1946)),
  Cold_War_End = list(pre = c(1970, 1990), post = c(1991, 2000))
)

CI_DVS <- c("kl_divergence", "js_divergence", "edge_turnover_pct",
             "degree_kl", "degree_wasserstein")

DV_LABELS <- c(
  kl_divergence       = "KL Divergence",
  js_divergence       = "JS Divergence",
  edge_turnover_pct   = "Edge Turnover (%)",
  degree_kl           = "Degree KL",
  degree_wasserstein  = "Degree Wasserstein"
)

# Covariate cols — full set requested for the new framework
COV_FULL <- intersect(c("n_democracies", "n_autocracies", "n_states"),
                       names(system_data))

run_one_ci <- function(war_name, dv, transform, controls) {
  war_def <- WARS[[war_name]]
  if (!dv %in% names(system_data)) {
    cat(sprintf("  Skipping %s/%s: column not found\n", war_name, dv))
    return(NULL)
  }

  ts_subset <- system_data |>
    filter(year >= war_def$pre[1], year <= war_def$post[2],
           !is.na(.data[[dv]])) |>
    arrange(year)

  if (nrow(ts_subset) < 10) {
    cat(sprintf("  Skipping %s/%s/%s/%s: too few obs (%d)\n",
                war_name, dv, transform, controls, nrow(ts_subset)))
    return(NULL)
  }

  # DV transform
  y <- ts_subset[[dv]]
  if (transform == "logged") y <- log1p(y)

  # Covariates (only if controls = "full"). Always log1p covariates --
  # they are non-negative, right-skewed, and the log keeps them on a
  # comparable scale to the logged DV.
  x <- if (controls == "full" && length(COV_FULL) > 0) {
    ts_subset |>
      select(all_of(COV_FULL)) |>
      mutate(across(everything(), log1p)) |>
      as.matrix()
  } else {
    NULL
  }

  pre_period  <- c(war_def$pre[1],  war_def$pre[2])
  post_period <- c(war_def$post[1], war_def$post[2])

  ci_data <- if (!is.null(x)) {
    zoo::zoo(cbind(y = y, x), order.by = ts_subset$year)
  } else {
    zoo::zoo(cbind(y = y), order.by = ts_subset$year)
  }

  cat(sprintf("  CI [%s | %s | %s | %s] n=%d, pre=%d-%d, post=%d-%d\n",
              war_name, dv, transform, controls, nrow(ts_subset),
              pre_period[1], pre_period[2], post_period[1], post_period[2]))

  impact <- tryCatch(
    CausalImpact(ci_data, pre.period = pre_period, post.period = post_period,
                  model.args = list(niter = 5000, nseasons = 1)),
    error = function(e) {
      cat(sprintf("  *** ERROR: %s\n", e$message))
      NULL
    }
  )
  if (is.null(impact)) return(NULL)

  list(
    war = war_name, dv = dv,
    transform = transform, controls = controls,
    impact = impact
  )
}

# Spec grid
spec_grid <- expand.grid(
  war       = names(WARS),
  dv        = CI_DVS,
  transform = c("logged", "raw"),
  controls  = c("none", "full"),
  stringsAsFactors = FALSE
)

ci_results <- list()
ci_summary_rows <- list()

for (i in seq_len(nrow(spec_grid))) {
  sp <- spec_grid[i, ]
  res <- run_one_ci(sp$war, sp$dv, sp$transform, sp$controls)
  if (is.null(res)) next

  key <- sprintf("%s__%s__%s__%s", sp$war, sp$dv, sp$transform, sp$controls)
  ci_results[[key]] <- res

  # Per-spec PNG
  war_label <- gsub("_", " ", sp$war) |> tools::toTitleCase()
  dv_clean <- ifelse(sp$dv %in% names(DV_LABELS), DV_LABELS[sp$dv],
                     gsub("_", " ", sp$dv) |> tools::toTitleCase())
  dv_disp <- if (sp$transform == "logged") {
    sprintf("log(1 + %s)", dv_clean)
  } else dv_clean
  ctrl_disp <- ifelse(sp$controls == "full",
                       "with controls", "no controls")
  p_ci <- plot(res$impact) +
    ggplot2::ggtitle(sprintf("Causal impact: %s", war_label)) +
    ggplot2::theme(
      strip.text = ggplot2::element_text(face = "bold", size = 9),
      panel.spacing = ggplot2::unit(1.2, "lines"),
      plot.margin = ggplot2::margin(t = 10, r = 12, b = 14, l = 10)
    )
  p_ci$facet$params$labeller <-
    ggplot2::as_labeller(function(x) tools::toTitleCase(x))
  ggplot2::ggsave(
    file.path(FIG_DIR,
              sprintf("causal_impact_%s_%s_%s_%s.png",
                       sp$war, sp$dv, sp$transform, sp$controls)),
    p_ci, width = 10, height = 8, dpi = 200)

  # Summary row
  s <- res$impact$summary
  cnames <- colnames(s)
  abs_col <- which(cnames == "AbsEffect")
  rel_col <- which(cnames == "RelEffect")
  p_col   <- which(cnames == "p")
  lo_col  <- which(cnames == "AbsEffect.lower")
  hi_col  <- which(cnames == "AbsEffect.upper")

  ci_summary_rows[[length(ci_summary_rows) + 1]] <- data.frame(
    war = sp$war, dv = sp$dv,
    transform = sp$transform, controls = sp$controls,
    spec_key  = key,
    avg_effect     = if (length(abs_col) == 1) as.numeric(s[1, abs_col]) else NA_real_,
    cum_effect     = if (length(abs_col) == 1) as.numeric(s[2, abs_col]) else NA_real_,
    avg_effect_rel = if (length(rel_col) == 1) as.numeric(s[1, rel_col]) else NA_real_,
    p_value        = if (length(p_col)   == 1) as.numeric(s[1, p_col])   else NA_real_,
    ci_lower       = if (length(lo_col)  == 1) as.numeric(s[1, lo_col])  else NA_real_,
    ci_upper       = if (length(hi_col)  == 1) as.numeric(s[1, hi_col])  else NA_real_,
    stringsAsFactors = FALSE
  )
}

# --------------------------------------------------------------------------- #
# 2. Summary table (combined; keyed by transform x controls)
# --------------------------------------------------------------------------- #

cat("\n--- Summary Table ---\n")
ci_summary <- as.data.frame(bind_rows(ci_summary_rows))
tryCatch(print(ci_summary), error = function(e) {
  cat("  (print failed:", e$message, ")\n")
})

# --------------------------------------------------------------------------- #
# 3. Combined panel plot — headline cell only
#    (logged DV + full controls), one facet per war, KL divergence.
# --------------------------------------------------------------------------- #

build_panel_plot <- function(transform, controls, file_name,
                              ylab_prefix = "KL Divergence") {
  panel_data <- list()
  for (war_name in names(WARS)) {
    key <- sprintf("%s__kl_divergence__%s__%s", war_name, transform, controls)
    if (!key %in% names(ci_results)) next
    impact <- ci_results[[key]]$impact
    series <- impact$series
    panel_data[[war_name]] <- tibble(
      year      = as.numeric(time(series)),
      observed  = as.numeric(series[, "response"]),
      predicted = as.numeric(series[, "point.pred"]),
      ci_lower  = as.numeric(series[, "point.pred.lower"]),
      ci_upper  = as.numeric(series[, "point.pred.upper"]),
      war       = war_name
    )
  }
  if (length(panel_data) == 0) return(invisible(NULL))

  war_order <- c(WWI = "WWI", WWII = "WWII",
                  "Cold War End" = "Cold War End")
  panel_df <- bind_rows(panel_data) |>
    mutate(war = gsub("_", " ", war),
           war = factor(war, levels = names(war_order),
                         labels = war_order))

  ylab_label <- if (transform == "logged") {
    sprintf("log(1 + %s)", ylab_prefix)
  } else ylab_prefix
  ctrl_disp <- ifelse(controls == "full",
                       "with full controls", "no controls")

  p_panel <- ggplot(panel_df, aes(x = year)) +
    geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper),
                fill = "lightblue", alpha = 0.5) +
    geom_line(aes(y = predicted), color = "blue", linetype = "dashed") +
    geom_line(aes(y = observed), color = "black") +
    facet_wrap(~ war, scales = "free", nrow = 1,
               labeller = label_wrap_gen(width = 22)) +
    labs(x = "Year",
         y = ylab_label,
         title = "Counterfactual vs observed") +
    theme_minimal(base_size = 14) +
    theme(strip.text   = element_text(face = "bold", size = 9),
          axis.text    = element_text(size = 11),
          axis.title   = element_text(size = 13),
          plot.title   = element_text(face = "bold", size = 16),
          panel.spacing = unit(1.2, "lines"),
          panel.spacing.x = unit(1.5, "lines"),
          plot.margin = margin(t = 10, r = 12, b = 14, l = 10))

  ggsave(file.path(FIG_DIR, file_name), p_panel,
         width = 16, height = 6, dpi = 300, bg = "white")
  cat(sprintf("  Panel plot [%s/%s] -> %s\n",
              transform, controls, file_name))
}

build_panel_plot("logged", "full", "fig_causal_impact_panel_logged_full.png")
build_panel_plot("logged", "none", "fig_causal_impact_panel_logged_none.png")
build_panel_plot("raw",    "full", "fig_causal_impact_panel_raw_full.png")
build_panel_plot("raw",    "none", "fig_causal_impact_panel_raw_none.png")

# --------------------------------------------------------------------------- #
# 4. Save results
# --------------------------------------------------------------------------- #

saveRDS(ci_results, file.path(OUTPUT_DIR, "causal_impact_results.rds"))

# Combined CSV (set / controls / spec_key columns)
write_csv(ci_summary,
          file.path(OUTPUT_DIR, "causal_impact_summary_combined.csv"))
# Convenience splits
for (tr in c("logged", "raw")) {
  for (ct in c("none", "full")) {
    sub <- ci_summary |> filter(transform == tr, controls == ct)
    if (nrow(sub) == 0) next
    write_csv(sub,
              file.path(OUTPUT_DIR,
                        sprintf("causal_impact_summary_%s_%s.csv", tr, ct)))
  }
}

# --- LaTeX tables: one per (transform, controls) cell -----------------------
write_ci_table <- function(transform, controls, file_path, caption, label) {
  sub <- ci_summary |> filter(transform == !!transform,
                                controls  == !!controls) |>
    arrange(war, dv)
  if (nrow(sub) == 0) {
    cat(sprintf("  No CI rows for %s/%s -- skipping table\n",
                transform, controls))
    return(invisible(NULL))
  }

  war_labels <- c(WWI = "World War I", WWII = "World War II",
                   Cold_War_End = "Cold War End")
  dv_labels <- c(
    kl_divergence       = "KL Divergence",
    js_divergence       = "JS Divergence",
    edge_turnover_pct   = "Edge Turnover (\\%)",
    degree_kl           = "Degree KL",
    degree_wasserstein  = "Degree Wasserstein"
  )
  tex_ci <- c(
    "\\begin{table}[htbp]", "\\centering",
    sprintf("\\caption{%s}", caption),
    sprintf("\\label{%s}", label),
    "\\small",
    "\\begin{tabular}{llrrrr}", "\\toprule",
    "Event & DV & Avg.\\ Effect & Rel.\\ Effect (\\%) & CI$_{95}$ & $p$ \\\\",
    "\\midrule"
  )
  for (i in seq_len(nrow(sub))) {
    r <- sub[i, ]
    w_lab <- ifelse(r$war %in% names(war_labels), war_labels[r$war], r$war)
    d_lab <- ifelse(r$dv %in% names(dv_labels),  dv_labels[r$dv],  r$dv)
    ci_str <- sprintf("[%.4f, %.4f]",
                      ifelse(is.na(r$ci_lower), NA_real_, r$ci_lower),
                      ifelse(is.na(r$ci_upper), NA_real_, r$ci_upper))
    tex_ci <- c(tex_ci, sprintf(
      "%s & %s & %.4f & %.1f & %s & %.3f \\\\",
      w_lab, d_lab,
      ifelse(is.na(r$avg_effect),     NA_real_, r$avg_effect),
      ifelse(is.na(r$avg_effect_rel), NA_real_, r$avg_effect_rel * 100),
      ci_str,
      ifelse(is.na(r$p_value),        NA_real_, r$p_value)
    ))
  }
  tex_ci <- c(tex_ci, "\\bottomrule", "\\end{tabular}",
    "\\begin{tablenotes}\\small",
    sprintf("\\item DV %s. Controls: %s.",
            ifelse(transform == "logged", "is log1p-transformed",
                                          "is on the raw scale"),
            ifelse(controls == "full",
                   "log1p(n\\_democracies), log1p(n\\_autocracies), log1p(n\\_states)",
                   "none (BSTS state-space level only)")),
    "\\end{tablenotes}",
    "\\end{table}")

  writeLines(tex_ci, file_path)
  cat(sprintf("  CI table [%s/%s] -> %s\n",
              transform, controls, basename(file_path)))
}

# Headline cell: logged + full controls -> main; others -> appendix.
write_ci_table("logged", "full",
               file.path(TABLE_DIR, "table_causal_impact_logged_full.tex"),
               "Bayesian Structural Time Series: Causal Impact of Major Wars on log(1 + Network Divergence), with Full Controls.",
               "tab:ci_logged_full")
write_ci_table("logged", "none",
               file.path(APPX_DIR,  "tableA_causal_impact_logged_none.tex"),
               "Appendix. Causal Impact on log(1 + Network Divergence), No Controls.",
               "tab:ci_logged_none")
write_ci_table("raw", "full",
               file.path(APPX_DIR,  "tableA_causal_impact_raw_full.tex"),
               "Appendix. Causal Impact on Raw Network Divergence, with Full Controls.",
               "tab:ci_raw_full")
write_ci_table("raw", "none",
               file.path(APPX_DIR,  "tableA_causal_impact_raw_none.tex"),
               "Appendix. Causal Impact on Raw Network Divergence, No Controls.",
               "tab:ci_raw_none")

cat("\nDone. Bayesian structural time series analysis complete.\n")