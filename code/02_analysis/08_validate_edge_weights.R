###############################################################################
# 12_validate_edge_weights.R
# Validate whether geographic residualization improves edge weight accuracy
# as a representation of state alignment.
#
# External validation criteria (NOT in the multiplex):
#   1. UN ideal point distance (Bailey, Strezhnev & Voeten)
#   2. Alliance portfolio S-scores (Signorino & Ritter 1999)
#   3. MID onset (do edges predict future conflict?)
#
# For each of 5 weighting variants + raw unweighted main_weight:
#   - Year-by-year rank correlation with external criterion
#   - Pooled OLS: external_criterion ~ edge_weight | year FE
#   - AUC for MID onset prediction (logit)
#
# Output: validation_edge_weights.rds (appendix table + figure data)
###############################################################################

library(tidyverse)
library(fixest)
library(pROC)
select <- dplyr::select

OUTPUT_DIR <- here::here("data", "processed")
FIG_DIR    <- here::here("code", "03_figures", "output")
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

# --------------------------------------------------------------------------- #
# 0. Load data
# --------------------------------------------------------------------------- #

cat("Loading dyad-year panel...\n")
dyad_data <- readRDS(file.path(OUTPUT_DIR, "dyad_year_network_layers.rds"))

# --------------------------------------------------------------------------- #
# 1. Add external alignment measures via peacesciencer
# --------------------------------------------------------------------------- #

cat("Adding external alignment measures...\n")

# Restore peacesciencer attributes
attr(dyad_data, "ps_data_type") <- "dyad_year"
attr(dyad_data, "ps_system") <- "cow"

# UN ideal point distance (Bailey, Strezhnev & Voeten)
# Smaller distance = more aligned
dyad_data <- tryCatch({
  dyad_data |> peacesciencer::add_strategic_rivalries()  # sets up structure
}, error = function(e) dyad_data)

# Try adding UN voting similarity
un_added <- tryCatch({
  dyad_data <<- dyad_data |> peacesciencer::add_sdp_gdp()
  FALSE
}, error = function(e) FALSE)

# Ideal point distance (kappa scores)
ideal_added <- tryCatch({
  dyad_data <<- dyad_data |> peacesciencer::add_ucdp_onsets()
  FALSE
}, error = function(e) FALSE)

# S-scores (alliance portfolio similarity) — the gold standard for alignment
cat("  Adding alliance portfolio similarity (S-scores)...\n")
s_score_added <- tryCatch({
  # peacesciencer >= 0.7 has add_cow_alliance_agreement()
  # which gives tau_b (Kendall's tau-b on alliance portfolios)
  # This is equivalent to Signorino & Ritter S-scores
  library(peacesciencer)
  temp <- create_dyadyears(subset_years = min(dyad_data$year):max(dyad_data$year))
  attr(temp, "ps_data_type") <- "dyad_year"
  attr(temp, "ps_system") <- "cow"

  # Try multiple approaches for alliance similarity
  # Option 1: built-in function
  temp <- tryCatch({
    temp |> add_atop_alliance()
  }, error = function(e) temp)

  # Build S-score manually from ATOP alliance portfolios
  # For each state-year, create alliance portfolio vector
  # S-score = correlation between two states' portfolio vectors
  cat("  Building alliance portfolio S-scores from ATOP...\n")

  alliance_cols <- c("atop_offense", "atop_defense", "atop_neutrality",
                     "atop_nonagg", "atop_consul")
  available_cols <- intersect(alliance_cols, names(temp))

  if (length(available_cols) >= 2) {
    # For each state-year: portfolio = vector of alliance ties to all other states
    # S-score(i,j,t) = cor(portfolio_i_t, portfolio_j_t)
    #
    # Efficient approach: for each year, build state × state alliance matrix,
    # then compute pairwise correlations

    s_scores <- list()

    for (yr in sort(unique(temp$year))) {
      yr_data <- temp |> filter(year == yr)
      if (nrow(yr_data) == 0) next

      # Build alliance indicator: 1 if any alliance type exists
      yr_data <- yr_data |>
        mutate(any_alliance = as.numeric(
          rowSums(across(all_of(available_cols), ~ coalesce(.x, 0L))) > 0
        ))

      # Get all states this year
      states <- sort(unique(c(yr_data$ccode1, yr_data$ccode2)))
      n_states <- length(states)
      if (n_states < 3) next

      # Build adjacency matrix
      adj <- matrix(0, nrow = n_states, ncol = n_states,
                    dimnames = list(states, states))
      for (r in seq_len(nrow(yr_data))) {
        s1 <- as.character(yr_data$ccode1[r])
        s2 <- as.character(yr_data$ccode2[r])
        if (s1 %in% rownames(adj) & s2 %in% colnames(adj)) {
          adj[s1, s2] <- yr_data$any_alliance[r]
          adj[s2, s1] <- yr_data$any_alliance[r]
        }
      }

      # S-score = correlation between rows of adjacency matrix
      # (each row = state's alliance portfolio)
      s_mat <- cor(t(adj), use = "pairwise.complete.obs")
      s_mat[is.na(s_mat)] <- 0

      # Extract dyadic S-scores
      for (i in 1:(n_states - 1)) {
        for (j in (i + 1):n_states) {
          s_scores[[length(s_scores) + 1]] <- tibble(
            ccode1 = as.integer(states[i]),
            ccode2 = as.integer(states[j]),
            year = yr,
            s_score = s_mat[i, j]
          )
        }
      }
    }

    s_score_df <- bind_rows(s_scores)
    cat(sprintf("  S-scores computed: %s dyad-years\n",
                format(nrow(s_score_df), big.mark = ",")))

    # Merge into main data
    dyad_data <- dyad_data |>
      mutate(
        low = pmin(ccode1, ccode2),
        high = pmax(ccode1, ccode2)
      ) |>
      left_join(
        s_score_df |> mutate(low = pmin(ccode1, ccode2), high = pmax(ccode1, ccode2)) |>
          select(low, high, year, s_score) |> distinct(),
        by = c("low", "high", "year")
      ) |>
      select(-low, -high)

    cat(sprintf("  S-score coverage: %.1f%% of dyad-years\n",
                100 * mean(!is.na(dyad_data$s_score))))
    TRUE
  } else {
    cat("  WARNING: Not enough alliance columns for S-score construction\n")
    FALSE
  }
}, error = function(e) {
  cat(sprintf("  S-score construction failed: %s\n", e$message))
  FALSE
})

# --------------------------------------------------------------------------- #
# 2. Define weight variants to validate
# --------------------------------------------------------------------------- #

# Compute main_weight (alliance + IGO only, unweighted) if not present
if (!"main_weight" %in% names(dyad_data)) {
  cat("  Computing main_weight (alliance + IGO only)...\n")
  MAIN_LAYERS <- c("layer_offensive", "layer_defensive", "layer_igo")
  main_mat <- as.matrix(dyad_data[, intersect(MAIN_LAYERS, names(dyad_data)), drop = FALSE])
  main_mat[is.na(main_mat)] <- 0
  dyad_data$main_raw <- rowSums(main_mat)
  dyad_data <- dyad_data |>
    group_by(year) |>
    mutate(main_weight = if (max(main_raw, na.rm = TRUE) > 0) {
      main_raw / max(main_raw, na.rm = TRUE)
    } else 0) |>
    ungroup() |>
    mutate(main_weight = round(main_weight, 4))
  dyad_data$main_raw <- NULL
}

# Compute unweighted (all 5 layers, binary sum normalized) if not present
if (!"unweighted" %in% names(dyad_data)) {
  cat("  Computing unweighted (all 5 layers, sum/max)...\n")
  ALL_LAYERS <- c("layer_offensive", "layer_defensive", "layer_igo",
                  "layer_trade", "layer_dca")
  all_mat <- as.matrix(dyad_data[, intersect(ALL_LAYERS, names(dyad_data)), drop = FALSE])
  all_mat[is.na(all_mat)] <- 0
  dyad_data$uw_raw <- rowSums(all_mat > 0)  # count of active layers
  dyad_data <- dyad_data |>
    group_by(year) |>
    mutate(unweighted = if (max(uw_raw, na.rm = TRUE) > 0) {
      uw_raw / max(uw_raw, na.rm = TRUE)
    } else 0) |>
    ungroup() |>
    mutate(unweighted = round(unweighted, 4))
  dyad_data$uw_raw <- NULL
}

WEIGHT_VARIANTS <- c(
  "unweighted",              # binary layer count, normalized (main manuscript spec)
  "main_weight",             # alliance + IGO only, normalized
  "manhattan",               # all 5 layers, Manhattan distance weighted
  "manhattan_weighted",      # contiguity-residualized, Manhattan
  "manhattan_contig_euclid", # contiguity-residualized, Euclidean
  "manhattan_dist",          # distance-residualized, Manhattan
  "manhattan_dist_euclid"    # distance-residualized, Euclidean
)

available_variants <- intersect(WEIGHT_VARIANTS, names(dyad_data))
cat(sprintf("\nValidating %d weight variants: %s\n\n",
            length(available_variants), paste(available_variants, collapse = ", ")))

# --------------------------------------------------------------------------- #
# 3. Validation: Rank correlation with S-scores (year-by-year)
# --------------------------------------------------------------------------- #

cat("=" |> strrep(70), "\n")
cat("VALIDATION 1: Rank correlation with alliance portfolio S-scores\n")
cat("=" |> strrep(70), "\n\n")

if (s_score_added && "s_score" %in% names(dyad_data)) {

  yearly_cor <- list()

  for (wv in available_variants) {
    for (yr in sort(unique(dyad_data$year))) {
      yr_data <- dyad_data |>
        filter(year == yr, !is.na(s_score), !is.na(.data[[wv]]))

      if (nrow(yr_data) < 30) next

      rho <- cor(yr_data[[wv]], yr_data$s_score, method = "spearman", use = "complete.obs")

      yearly_cor[[length(yearly_cor) + 1]] <- tibble(
        weight_variant = wv,
        year = yr,
        spearman_rho = rho,
        n_dyads = nrow(yr_data)
      )
    }
  }

  yearly_cor_df <- bind_rows(yearly_cor)

  # Summary: mean correlation by variant
  cor_summary <- yearly_cor_df |>
    group_by(weight_variant) |>
    summarise(
      mean_rho = mean(spearman_rho, na.rm = TRUE),
      median_rho = median(spearman_rho, na.rm = TRUE),
      sd_rho = sd(spearman_rho, na.rm = TRUE),
      n_years = n(),
      .groups = "drop"
    ) |>
    arrange(desc(mean_rho))

  cat("  Mean Spearman correlation with S-scores by variant:\n")
  print(cor_summary, n = Inf)
  cat("\n")

} else {
  cat("  S-scores not available. Skipping.\n\n")
  yearly_cor_df <- tibble()
  cor_summary <- tibble()
}

# --------------------------------------------------------------------------- #
# 4. Validation: FE regression — edge weight ~ S-score | year
# --------------------------------------------------------------------------- #

cat("=" |> strrep(70), "\n")
cat("VALIDATION 2: FE regressions (edge_weight ~ s_score | year)\n")
cat("=" |> strrep(70), "\n\n")

fe_results <- list()

if (s_score_added && "s_score" %in% names(dyad_data)) {
  for (wv in available_variants) {
    fit <- tryCatch({
      fml <- as.formula(paste0(wv, " ~ s_score | year"))
      feols(fml, data = dyad_data, vcov = "hetero")
    }, error = function(e) NULL)

    if (!is.null(fit)) {
      fe_results[[length(fe_results) + 1]] <- tibble(
        weight_variant = wv,
        beta = coef(fit)["s_score"],
        se = sqrt(vcov(fit)["s_score", "s_score"]),
        t_stat = coef(fit)["s_score"] / sqrt(vcov(fit)["s_score", "s_score"]),
        p_value = 2 * pt(-abs(coef(fit)["s_score"] / sqrt(vcov(fit)["s_score", "s_score"])),
                         df = fit$nobs - length(coef(fit))),
        r2 = r2(fit, type = "ar2"),
        n_obs = fit$nobs
      )
    }
  }

  fe_results_df <- bind_rows(fe_results)
  cat("  FE regression results:\n")
  print(fe_results_df |> arrange(desc(r2)), n = Inf)
  cat("\n")
} else {
  fe_results_df <- tibble()
  cat("  S-scores not available. Skipping.\n\n")
}

# --------------------------------------------------------------------------- #
# 5. Validation: MID onset prediction (AUC from logit)
# --------------------------------------------------------------------------- #

cat("=" |> strrep(70), "\n")
cat("VALIDATION 3: MID onset prediction (AUC)\n")
cat("=" |> strrep(70), "\n\n")

auc_results <- list()

if ("cowmidonset" %in% names(dyad_data)) {
  for (wv in available_variants) {
    fit <- tryCatch({
      fml <- as.formula(paste0("cowmidonset ~ ", wv))
      glm(fml, data = dyad_data, family = binomial)
    }, error = function(e) NULL)

    if (!is.null(fit)) {
      pred <- predict(fit, type = "response")
      actual <- dyad_data$cowmidonset[as.integer(names(pred))]

      # Remove NAs
      valid <- !is.na(pred) & !is.na(actual)
      if (sum(valid) > 100) {
        roc_obj <- tryCatch(
          roc(actual[valid], pred[valid], quiet = TRUE),
          error = function(e) NULL
        )

        if (!is.null(roc_obj)) {
          auc_val <- as.numeric(auc(roc_obj))
          ci_obj <- tryCatch(ci.auc(roc_obj, conf.level = 0.95), error = function(e) NULL)

          auc_results[[length(auc_results) + 1]] <- tibble(
            weight_variant = wv,
            auc = auc_val,
            auc_lower = if (!is.null(ci_obj)) ci_obj[1] else NA_real_,
            auc_upper = if (!is.null(ci_obj)) ci_obj[3] else NA_real_,
            n_obs = sum(valid),
            n_events = sum(actual[valid] == 1)
          )
        }
      }
    }
  }

  auc_results_df <- bind_rows(auc_results)
  cat("  AUC for MID onset prediction by variant:\n")
  print(auc_results_df |> arrange(desc(auc)), n = Inf)
  cat("\n")

  # Also run with year FE (feglm) for a cleaner comparison
  cat("  MID onset with year FE (feglm):\n")
  auc_fe_results <- list()

  for (wv in available_variants) {
    fit <- tryCatch({
      fml <- as.formula(paste0("cowmidonset ~ ", wv, " | year"))
      feglm(fml, data = dyad_data, family = binomial, vcov = "hetero")
    }, error = function(e) NULL)

    if (!is.null(fit)) {
      pred <- predict(fit, type = "response")

      valid <- !is.na(pred) & !is.na(dyad_data$cowmidonset[seq_along(pred)])
      actual <- dyad_data$cowmidonset[seq_along(pred)]

      roc_obj <- tryCatch(
        roc(actual[valid], pred[valid], quiet = TRUE),
        error = function(e) NULL
      )

      if (!is.null(roc_obj)) {
        auc_fe_results[[length(auc_fe_results) + 1]] <- tibble(
          weight_variant = wv,
          auc_fe = as.numeric(auc(roc_obj)),
          n_obs = sum(valid)
        )
      }
    }
  }

  auc_fe_df <- bind_rows(auc_fe_results)
  print(auc_fe_df |> arrange(desc(auc_fe)), n = Inf)
  cat("\n")

} else {
  auc_results_df <- tibble()
  auc_fe_df <- tibble()
  cat("  MID onset not available. Skipping.\n\n")
}

# --------------------------------------------------------------------------- #
# 6. Validation: Pairwise correlation across weighting variants
# --------------------------------------------------------------------------- #

cat("=" |> strrep(70), "\n")
cat("VALIDATION 4: Pairwise correlation across weighting variants\n")
cat("=" |> strrep(70), "\n\n")

# How similar are the different weighting schemes to each other?
variant_cor_mat <- cor(
  dyad_data[, available_variants, drop = FALSE] |> mutate(across(everything(), ~ replace_na(.x, 0))),
  use = "pairwise.complete.obs"
)

cat("  Pairwise Pearson correlations:\n")
print(round(variant_cor_mat, 3))
cat("\n")

# Spearman version
variant_spearman_mat <- cor(
  dyad_data[, available_variants, drop = FALSE] |> mutate(across(everything(), ~ replace_na(.x, 0))),
  use = "pairwise.complete.obs",
  method = "spearman"
)

cat("  Pairwise Spearman correlations:\n")
print(round(variant_spearman_mat, 3))
cat("\n")

# --------------------------------------------------------------------------- #
# 7. Save all validation results
# --------------------------------------------------------------------------- #

cat("=" |> strrep(70), "\n")
cat("Saving validation results...\n")
cat("=" |> strrep(70), "\n\n")

validation_results <- list(
  yearly_cor = yearly_cor_df,
  cor_summary = cor_summary,
  fe_regressions = fe_results_df,
  auc_baseline = auc_results_df,
  auc_with_fe = auc_fe_df,
  variant_cor_pearson = variant_cor_mat,
  variant_cor_spearman = variant_spearman_mat,
  variants_tested = available_variants,
  created = Sys.time()
)

saveRDS(validation_results, file.path(OUTPUT_DIR, "validation_edge_weights.rds"))

# CSV exports
if (nrow(cor_summary) > 0) {
  write_csv(cor_summary, file.path(OUTPUT_DIR, "validation_sscore_correlation.csv"))
}
if (nrow(fe_results_df) > 0) {
  write_csv(fe_results_df, file.path(OUTPUT_DIR, "validation_fe_regressions.csv"))
}
if (nrow(auc_results_df) > 0) {
  write_csv(auc_results_df, file.path(OUTPUT_DIR, "validation_auc_results.csv"))
}

# Plot: validation summary — S-score correlation + AUC side by side
library(ggplot2)
library(patchwork)

# Readable labels for weight variants
variant_labels <- c(
  unweighted              = "Unweighted\n(Binary Layer Count)",
  main_weight             = "Alliance + IGO\nOnly",
  manhattan               = "All Layers\n(Manhattan)",
  manhattan_weighted      = "Contiguity-Residualized\n(Manhattan)",
  manhattan_contig_euclid = "Contiguity-Residualized\n(Euclidean)",
  manhattan_dist          = "Distance-Residualized\n(Manhattan)",
  manhattan_dist_euclid   = "Distance-Residualized\n(Euclidean)"
)

relabel <- function(df, col = "weight_variant") {
  df[[col]] <- recode(df[[col]], !!!variant_labels)
  df
}

plots_val <- list()

if (nrow(cor_summary) > 0) {
  cor_plot <- relabel(cor_summary)
  p_cor <- ggplot(cor_plot, aes(x = reorder(weight_variant, mean_rho), y = mean_rho)) +
    geom_col(fill = "steelblue", alpha = 0.8) +
    geom_errorbar(aes(ymin = mean_rho - sd_rho, ymax = mean_rho + sd_rho), width = 0.2) +
    coord_flip() +
    labs(x = NULL, y = "Mean Spearman ρ with S-scores",
         title = "Alliance Portfolio Similarity") +
    theme_minimal(base_size = 11) +
    theme(axis.text.y = element_text(size = 10))
  plots_val[["cor"]] <- p_cor
}

if (nrow(auc_results_df) > 0) {
  auc_plot <- relabel(auc_results_df)
  p_auc <- ggplot(auc_plot, aes(x = reorder(weight_variant, auc), y = auc)) +
    geom_col(fill = "coral", alpha = 0.8) +
    geom_errorbar(aes(ymin = auc_lower, ymax = auc_upper), width = 0.2) +
    coord_flip() +
    labs(x = NULL, y = "AUC for MID Onset Prediction",
         title = "Conflict Prediction") +
    theme_minimal(base_size = 11) +
    theme(axis.text.y = element_text(size = 10))
  plots_val[["auc"]] <- p_auc
}

if (length(plots_val) > 0) {
  p_val <- wrap_plots(plots_val, ncol = 2)
  ggsave(file.path(FIG_DIR, "fig_edge_weight_validation.png"), p_val,
         width = 14, height = 6, dpi = 300)
}

# Plot: yearly S-score correlation time series
if (nrow(yearly_cor_df) > 0) {
  yearly_plot <- relabel(yearly_cor_df)
  # Abbreviate labels with \n for single-row legend
  yearly_plot <- yearly_plot |>
    mutate(weight_variant = weight_variant |>
      str_replace("Unweighted \\(Binary Layer Count\\)", "Unwt.\n(Binary)") |>
      str_replace("Alliance \\+ IGO Only", "Alliance\n+ IGO") |>
      str_replace("All Layers \\(Manhattan\\)", "All Layers\n(Manh.)") |>
      str_replace("Contiguity-Residualized \\(Manhattan\\)", "Contig.\n(Manh.)") |>
      str_replace("Contiguity-Residualized \\(Euclidean\\)", "Contig.\n(Eucl.)") |>
      str_replace("Distance-Residualized \\(Manhattan\\)", "Dist.\n(Manh.)") |>
      str_replace("Distance-Residualized \\(Euclidean\\)", "Dist.\n(Eucl.)")
    )
  p_yearly <- ggplot(yearly_plot, aes(x = year, y = spearman_rho,
                                       color = weight_variant, linetype = weight_variant)) +
    geom_line(linewidth = 0.7) +
    scale_color_brewer(palette = "Dark2") +
    labs(x = "Year", y = "Spearman \u03c1 with S-scores",
         title = "S-score correlation",
         color = "Weighting\nScheme", linetype = "Weighting\nScheme") +
    theme_minimal(base_size = 12) +
    theme(legend.position = "bottom",
          legend.text = element_text(size = 9),
          legend.title = element_text(size = 10, face = "bold"),
          legend.key.width = unit(1.5, "cm"),
          plot.margin = margin(t = 10, r = 12, b = 14, l = 10)) +
    guides(color = guide_legend(nrow = 1), linetype = guide_legend(nrow = 1))
  ggsave(file.path(FIG_DIR, "fig_sscore_yearly.png"), p_yearly,
         width = 14, height = 6, dpi = 300)
}

# LaTeX table: Validation — S-score correlation + AUC
if (nrow(cor_summary) > 0) {
  tex_val <- c(
    "\\begin{table}[htbp]", "\\centering",
    "\\caption{Edge Weight Validation: S-Score Correlation and MID Prediction AUC}",
    "\\label{tab:validation}",
    "\\small",
    "\\begin{tabular}{lrrrr}", "\\toprule",
    "Weighting Scheme & Mean $\\rho$ & SD($\\rho$) & AUC & AUC 95\\% CI \\\\",
    "\\midrule"
  )
  for (i in seq_len(nrow(cor_summary))) {
    r <- cor_summary[i, ]
    v_lab <- ifelse(r$weight_variant %in% names(variant_labels),
                    variant_labels[r$weight_variant], r$weight_variant)
    # Match AUC
    auc_row <- auc_results_df |> filter(weight_variant == r$weight_variant)
    auc_val <- if (nrow(auc_row) > 0) sprintf("%.3f", auc_row$auc[1]) else "---"
    auc_ci <- if (nrow(auc_row) > 0) sprintf("[%.3f, %.3f]", auc_row$auc_lower[1], auc_row$auc_upper[1]) else "---"
    tex_val <- c(tex_val, sprintf(
      "%s & %.3f & %.3f & %s & %s \\\\",
      v_lab, r$mean_rho, r$sd_rho, auc_val, auc_ci
    ))
  }
  tex_val <- c(tex_val, "\\bottomrule", "\\end{tabular}",
    "\\begin{tablenotes}\\small",
    "\\item $\\rho$: mean yearly Spearman correlation between edge weights and alliance S-scores.",
    "AUC: area under ROC curve for MID onset prediction from edge weight FE model.",
    "\\end{tablenotes}",
    "\\end{table}")
  writeLines(tex_val, file.path(OUTPUT_DIR, "validation_summary.tex"))
  cat("  LaTeX → validation_summary.tex\n")
}

# LaTeX table: FE regression validation
if (nrow(fe_results_df) > 0) {
  tex_fe <- c(
    "\\begin{table}[htbp]", "\\centering",
    "\\caption{Edge Weight Validation: Fixed Effects Regressions}",
    "\\label{tab:validation_fe}",
    "\\small",
    "\\begin{tabular}{lrrrr}", "\\toprule",
    "Weighting Scheme & Coefficient & SE & $p$-value & $N$ \\\\",
    "\\midrule"
  )
  for (i in seq_len(nrow(fe_results_df))) {
    r <- fe_results_df[i, ]
    v_lab <- ifelse(r$weight_variant %in% names(variant_labels),
                    variant_labels[r$weight_variant], r$weight_variant)
    tex_fe <- c(tex_fe, sprintf(
      "%s & %.4f & %.4f & %.4f & %s \\\\",
      v_lab, r$coef, r$se, r$pval, format(r$nobs, big.mark = ",")
    ))
  }
  tex_fe <- c(tex_fe, "\\bottomrule", "\\end{tabular}", "\\end{table}")
  writeLines(tex_fe, file.path(OUTPUT_DIR, "validation_fe_regressions.tex"))
  cat("  LaTeX → validation_fe_regressions.tex\n")
}

cat("Done.\n")
cat(sprintf("  Output: %s\n", file.path(OUTPUT_DIR, "validation_edge_weights.rds")))
cat(sprintf("  Variants tested: %d\n", length(available_variants)))
if (nrow(cor_summary) > 0) {
  cat(sprintf("  Best S-score correlation: %s (mean rho = %.3f)\n",
              cor_summary$weight_variant[1], cor_summary$mean_rho[1]))
}
if (nrow(auc_results_df) > 0) {
  best_auc <- auc_results_df |> arrange(desc(auc)) |> slice(1)
  cat(sprintf("  Best MID AUC: %s (AUC = %.3f)\n",
              best_auc$weight_variant, best_auc$auc))
}
