###############################################################################
# 11_winner_loser_analysis.R
# War outcomes, capability change, and centrality reconfiguration
#
# Appendix analysis testing whether:
#   (a) War winners/losers experience differential centrality changes
#   (b) Delta CINC (start→end of war) predicts centrality change
#   (c) The capability-centrality link differs pre/post-WWII (cutoff: 1946)
#
# Three edge weightings: manhattan (primary), contiguity, distance
# Three designs: (a) winners vs losers, (b) full panel TWFE, (c) pre/post windows
###############################################################################

library(tidyverse)
library(igraph)
library(peacesciencer)
library(fixest)
library(texreg)
select <- dplyr::select

OUTPUT_DIR <- here::here("data", "processed")
FIG_DIR    <- here::here("code", "03_figures", "output")
APPX_DIR   <- here::here("results", "appendix")
dir.create(APPX_DIR, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# 0. Load data
# ============================================================

cat("Loading dyad data...\n")
dyad_data <- readRDS(file.path(OUTPUT_DIR, "dyad_year_network_layers.rds"))
years <- sort(unique(dyad_data$year))

# Compute main_weight (alliance + IGO only) if not already present
# This matches the fixed-panel main specification from 02_compute_network_measures.R
if (!"main_weight" %in% names(dyad_data)) {
  cat("  Computing main_weight (alliance + IGO only)...\n")
  MAIN_LAYERS <- c("layer_offensive", "layer_defensive", "layer_igo")
  main_mat <- as.matrix(dyad_data[, MAIN_LAYERS, drop = FALSE])
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

# ============================================================
# 1. War outcomes from COW interstate wars
# ============================================================

cat("Loading COW war data...\n")
data(cow_war_inter, package = "peacesciencer")

# State-level outcomes (ccode1 perspective)
state_war_outcomes <- cow_war_inter |>
  select(warnum, ccode = ccode1, year, outcome = outcome1, batdeath = batdeath1) |>
  distinct(warnum, ccode, year, .keep_all = TRUE)

# State-war level: start/end years and outcome
state_wars <- state_war_outcomes |>
  group_by(warnum, ccode) |>
  summarise(
    war_start = min(year),
    war_end = max(year),
    outcome = first(outcome),
    total_batdeath = sum(batdeath, na.rm = TRUE),
    .groups = "drop"
  ) |>
  mutate(
    war_winner = as.integer(outcome == 1),
    war_loser = as.integer(outcome == 2),
    war_draw = as.integer(outcome %in% c(3, 6)),
    war_duration = war_end - war_start + 1
  )

cat(sprintf("  State-war observations: %d (wars: %d, states: %d)\n",
            nrow(state_wars), n_distinct(state_wars$warnum), n_distinct(state_wars$ccode)))
cat(sprintf("  Winners: %d, Losers: %d, Draws: %d, Other: %d\n",
            sum(state_wars$war_winner), sum(state_wars$war_loser),
            sum(state_wars$war_draw),
            sum(!state_wars$war_winner & !state_wars$war_loser & !state_wars$war_draw)))

# War-end year indicators for panel merge
winner_years <- state_wars |>
  filter(war_winner == 1) |>
  select(ccode, year = war_end) |> distinct() |>
  mutate(war_winner = 1L)

loser_years <- state_wars |>
  filter(war_loser == 1) |>
  select(ccode, year = war_end) |> distinct() |>
  mutate(war_loser = 1L)

any_war_years <- state_war_outcomes |>
  distinct(ccode, year) |>
  mutate(in_cow_war = 1L)

# ============================================================
# 2. State-year CINC
# ============================================================

state_cinc <- dyad_data |>
  select(ccode = ccode1, year, cinc = cinc1) |>
  distinct(ccode, year, .keep_all = TRUE)

# State-level MID variables (for war IV × post_wwii interactions)
cat("Loading state-level MID variables...\n")
state_mid_vars <- readRDS(file.path(OUTPUT_DIR, "state_level_results.rds"))$data |>
  select(ccode, year, mid_onset_count, fatal_mid_onset_count,
         mid_count, fatal_mid_count) |>
  distinct(ccode, year, .keep_all = TRUE) |>
  mutate(
    log_mid_onset_count       = log1p(mid_onset_count),
    log_fatal_mid_onset_count = log1p(fatal_mid_onset_count),
    log_mid_count             = log1p(mid_count),
    log_fatal_mid_count       = log1p(fatal_mid_count)
  )

# Merge CINC at war start/end for delta CINC computation
state_wars <- state_wars |>
  left_join(state_cinc |> rename(cinc_start = cinc), by = c("ccode", "war_start" = "year")) |>
  left_join(state_cinc |> rename(cinc_end = cinc), by = c("ccode", "war_end" = "year")) |>
  mutate(
    delta_cinc_war = cinc_end - cinc_start,
    pct_cinc_chg_war = (cinc_end - cinc_start) / cinc_start
  )

cat(sprintf("  Delta CINC (war): mean=%.6f, sd=%.6f\n",
            mean(state_wars$delta_cinc_war, na.rm = TRUE),
            sd(state_wars$delta_cinc_war, na.rm = TRUE)))

# ============================================================
# 3. Compute centrality and run models for each weighting
# ============================================================

# Primary: main_weight (alliance + IGO only, fixed panel)
# Robustness: manhattan (all layers), contiguity-weighted, distance-weighted
weight_specs <- c("main_weight", "manhattan", "manhattan_contig_euclid", "manhattan_dist")
weight_labels <- c("Main (Alliance+IGO)", "All Layers Manhattan", "Contiguity-weighted", "Distance-weighted")

all_results <- list()

for (wi in seq_along(weight_specs)) {
  wt <- weight_specs[wi]
  wt_label <- weight_labels[wi]

  cat(sprintf("\n\n=== EDGE WEIGHT: %s (%s) ===\n", wt, wt_label))

  # --- 3a. Compute centrality per state-year ---
  centrality_list <- list()

  for (yr in years) {
    edges <- dyad_data |>
      filter(year == yr, !!sym(wt) > 0) |>
      select(ccode1, ccode2, w = !!sym(wt)) |>
      distinct(ccode1, ccode2, .keep_all = TRUE)

    if (nrow(edges) == 0) next

    g <- graph_from_data_frame(edges |> select(ccode1, ccode2), directed = FALSE)
    E(g)$weight <- edges$w
    nodes <- V(g)$name
    n <- vcount(g)
    inv_w <- 1 / E(g)$weight

    bt <- if (n > 2) betweenness(g, weights = inv_w, normalized = TRUE) else rep(0, n)
    ev <- if (n > 2) tryCatch(eigen_centrality(g, weights = E(g)$weight)$vector,
                               error = function(e) rep(NA, n)) else rep(NA, n)
    pr <- if (n > 2) page_rank(g, weights = E(g)$weight)$vector else rep(1/n, n)
    dg <- degree(g, normalized = TRUE)

    centrality_list[[as.character(yr)]] <- tibble(
      ccode = as.integer(nodes), year = yr,
      betweenness = bt, eigenvector = ev, pagerank = pr, degree_norm = dg
    )
  }

  state_cent <- bind_rows(centrality_list) |>
    group_by(ccode) |> arrange(year) |>
    mutate(
      betweenness_chg = betweenness - lag(betweenness),
      eigenvector_chg = eigenvector - lag(eigenvector),
      pagerank_chg = pagerank - lag(pagerank),
      degree_chg = degree_norm - lag(degree_norm)
    ) |> ungroup()

  cat(sprintf("  Centrality: %d state-years\n", nrow(state_cent)))

  # --- 3b. Build panel ---
  panel <- state_cent |>
    left_join(winner_years, by = c("ccode", "year")) |>
    left_join(loser_years, by = c("ccode", "year")) |>
    left_join(any_war_years, by = c("ccode", "year")) |>
    left_join(state_cinc, by = c("ccode", "year")) |>
    left_join(state_mid_vars, by = c("ccode", "year")) |>
    mutate(
      war_winner = replace_na(war_winner, 0L),
      war_loser = replace_na(war_loser, 0L),
      in_cow_war = replace_na(in_cow_war, 0L),
      post_wwii = as.integer(year >= 1946)
    ) |>
    group_by(ccode) |> arrange(year) |>
    mutate(delta_cinc = cinc - lag(cinc)) |>
    ungroup()

  pre_panel <- panel |> filter(year < 1946)
  post_panel <- panel |> filter(year >= 1946)

  dvs_chg <- c("betweenness_chg", "eigenvector_chg", "pagerank_chg", "degree_chg")
  dvs_diff <- c("betweenness_diff", "eigenvector_diff", "pagerank_diff", "degree_diff")

  wt_results <- list()

  # ============================================================
  # MODEL SET 1: Winner/Loser TWFE (full panel + split sample)
  # ============================================================
  cat("  Running winner/loser TWFE models...\n")

  wl_models <- list()

  # Full panel
  for (dv in dvs_chg) {
    f <- as.formula(paste(dv, "~ war_winner + war_loser | ccode + year"))
    wl_models[[paste0("full_", dv)]] <- tryCatch(
      feols(f, data = panel, vcov = ~ccode), error = function(e) NULL)

    f_cinc <- as.formula(paste(dv, "~ war_winner + war_loser + cinc | ccode + year"))
    wl_models[[paste0("full_cinc_", dv)]] <- tryCatch(
      feols(f_cinc, data = panel, vcov = ~ccode), error = function(e) NULL)
  }

  # Split sample
  for (period in c("pre", "post")) {
    d <- if (period == "pre") pre_panel else post_panel
    for (dv in dvs_chg) {
      f <- as.formula(paste(dv, "~ war_winner + war_loser | ccode + year"))
      wl_models[[paste0(period, "_", dv)]] <- tryCatch(
        feols(f, data = d, vcov = ~ccode), error = function(e) NULL)
    }
  }

  # Interaction
  for (dv in dvs_chg) {
    f <- as.formula(paste(dv, "~ war_winner * post_wwii + war_loser * post_wwii | ccode + year"))
    wl_models[[paste0("interact_", dv)]] <- tryCatch(
      feols(f, data = panel, vcov = ~ccode), error = function(e) NULL)
  }

  wt_results$wl_models <- wl_models

  # ============================================================
  # MODEL SET 1b: War IV × post_wwii interaction (TWFE)
  # Tests whether the effect of MID exposure on centrality change
  # differs pre vs post WWII
  # ============================================================
  cat("  Running war IV × post_wwii interaction models...\n")

  war_iv_vars <- c("mid_onset_count", "log_mid_onset_count",
                   "fatal_mid_onset_count", "log_fatal_mid_onset_count",
                   "mid_count", "log_mid_count",
                   "fatal_mid_count", "log_fatal_mid_count")

  wariv_models <- list()

  for (dv in dvs_chg) {
    for (iv in war_iv_vars) {
      # Full interaction: iv * post_wwii
      f_int <- as.formula(paste(dv, "~", iv, "* post_wwii | ccode + year"))
      wariv_models[[paste0("interact_", dv, "_", iv)]] <- tryCatch(
        feols(f_int, data = panel, vcov = ~ccode), error = function(e) NULL)

      # Split sample: pre-WWII
      f_base <- as.formula(paste(dv, "~", iv, "| ccode + year"))
      wariv_models[[paste0("pre_", dv, "_", iv)]] <- tryCatch(
        feols(f_base, data = pre_panel, vcov = ~ccode), error = function(e) NULL)

      # Split sample: post-WWII
      wariv_models[[paste0("post_", dv, "_", iv)]] <- tryCatch(
        feols(f_base, data = post_panel, vcov = ~ccode), error = function(e) NULL)
    }
  }

  wt_results$wariv_models <- wariv_models

  # ============================================================
  # MODEL SET 2: Delta CINC → centrality (TWFE, split sample)
  # ============================================================
  cat("  Running delta CINC TWFE models...\n")

  cinc_models <- list()

  # Full panel
  for (dv in dvs_chg) {
    f <- as.formula(paste(dv, "~ delta_cinc | ccode + year"))
    cinc_models[[paste0("full_", dv)]] <- tryCatch(
      feols(f, data = panel, vcov = ~ccode), error = function(e) NULL)
  }

  # Split sample
  for (period in c("pre", "post")) {
    d <- if (period == "pre") pre_panel else post_panel
    for (dv in dvs_chg) {
      f <- as.formula(paste(dv, "~ delta_cinc | ccode + year"))
      cinc_models[[paste0(period, "_", dv)]] <- tryCatch(
        feols(f, data = d, vcov = ~ccode), error = function(e) NULL)
    }
  }

  # War participants only (full + split)
  for (period in c("full", "pre", "post")) {
    d <- switch(period,
      full = panel |> filter(in_cow_war == 1),
      pre = pre_panel |> filter(in_cow_war == 1),
      post = post_panel |> filter(in_cow_war == 1)
    )
    for (dv in dvs_chg) {
      f <- as.formula(paste(dv, "~ delta_cinc | ccode + year"))
      cinc_models[[paste0("war_", period, "_", dv)]] <- tryCatch(
        feols(f, data = d, vcov = ~ccode), error = function(e) NULL)
    }
  }

  wt_results$cinc_models <- cinc_models

  # ============================================================
  # MODEL SET 3: Delta CINC (start→end of war) → pre/post centrality
  # ============================================================
  cat("  Running pre/post war centrality models...\n")

  war_events <- state_wars |> filter(outcome %in% c(1, 2))
  prepost_list <- list()

  for (i in 1:nrow(war_events)) {
    cc <- war_events$ccode[i]
    ws <- war_events$war_start[i]
    we <- war_events$war_end[i]

    pre_d <- state_cent |> filter(ccode == cc, year >= ws - 3, year < ws)
    post_d <- state_cent |> filter(ccode == cc, year > we, year <= we + 3)

    if (nrow(pre_d) > 0 & nrow(post_d) > 0) {
      prepost_list[[length(prepost_list) + 1]] <- tibble(
        ccode = cc, warnum = war_events$warnum[i],
        war_start = ws, war_end = we,
        war_duration = war_events$war_duration[i],
        outcome = war_events$outcome[i],
        war_winner = war_events$war_winner[i],
        war_loser = war_events$war_loser[i],
        delta_cinc_war = war_events$delta_cinc_war[i],
        cinc_start = war_events$cinc_start[i],
        betweenness_diff = mean(post_d$betweenness, na.rm = TRUE) - mean(pre_d$betweenness, na.rm = TRUE),
        eigenvector_diff = mean(post_d$eigenvector, na.rm = TRUE) - mean(pre_d$eigenvector, na.rm = TRUE),
        pagerank_diff = mean(post_d$pagerank, na.rm = TRUE) - mean(pre_d$pagerank, na.rm = TRUE),
        degree_diff = mean(post_d$degree_norm, na.rm = TRUE) - mean(pre_d$degree_norm, na.rm = TRUE)
      )
    }
  }

  war_cent <- bind_rows(prepost_list) |>
    mutate(post_wwii = as.integer(war_end >= 1946))

  cat(sprintf("  Pre/post observations: %d (winners: %d, losers: %d)\n",
              nrow(war_cent), sum(war_cent$war_winner), sum(war_cent$war_loser)))

  prepost_models <- list()

  # A. Delta CINC (war) → centrality change (OLS, full sample)
  for (dv in dvs_diff) {
    prepost_models[[paste0("ols_full_", dv)]] <- tryCatch(
      lm(as.formula(paste(dv, "~ delta_cinc_war")), data = war_cent),
      error = function(e) NULL)

    prepost_models[[paste0("ols_ctrl_", dv)]] <- tryCatch(
      lm(as.formula(paste(dv, "~ delta_cinc_war + cinc_start + war_duration")), data = war_cent),
      error = function(e) NULL)
  }

  # B. Split sample
  for (period in c("pre", "post")) {
    d <- if (period == "pre") war_cent |> filter(post_wwii == 0) else war_cent |> filter(post_wwii == 1)
    for (dv in dvs_diff) {
      prepost_models[[paste0("ols_", period, "_", dv)]] <- tryCatch(
        lm(as.formula(paste(dv, "~ delta_cinc_war")), data = d),
        error = function(e) NULL)
    }
  }

  # C. Interaction: delta_cinc_war × post_wwii
  for (dv in dvs_diff) {
    prepost_models[[paste0("interact_", dv)]] <- tryCatch(
      lm(as.formula(paste(dv, "~ delta_cinc_war * post_wwii")), data = war_cent),
      error = function(e) NULL)
  }

  # D. Within winners / within losers
  for (grp in c("winner", "loser")) {
    d <- if (grp == "winner") war_cent |> filter(war_winner == 1) else war_cent |> filter(war_loser == 1)
    for (dv in dvs_diff) {
      prepost_models[[paste0(grp, "_", dv)]] <- tryCatch(
        lm(as.formula(paste(dv, "~ delta_cinc_war")), data = d),
        error = function(e) NULL)
    }
  }

  wt_results$prepost_models <- prepost_models
  wt_results$war_cent_data <- war_cent
  wt_results$panel_summary <- tibble(
    n_panel = nrow(panel),
    n_pre = nrow(pre_panel),
    n_post = nrow(post_panel),
    n_winners = sum(panel$war_winner),
    n_losers = sum(panel$war_loser),
    n_prepost = nrow(war_cent)
  )

  all_results[[wt]] <- wt_results
}

# ============================================================
# 4. Print summary of key results
# ============================================================

cat("\n\n")
cat(paste(rep("=", 72), collapse = ""))
cat("\n  SUMMARY OF KEY RESULTS\n")
cat(paste(rep("=", 72), collapse = ""))

for (wi in seq_along(weight_specs)) {
  wt <- weight_specs[wi]
  res <- all_results[[wt]]

  cat(sprintf("\n\n--- %s (%s) ---\n", wt, weight_labels[wi]))

  # Winner/Loser interaction
  cat("  Winner/Loser × post_wwii interaction (eigenvector):\n")
  m <- res$wl_models[["interact_eigenvector_chg"]]
  if (!is.null(m)) {
    s <- summary(m)
    for (v in intersect(names(coef(m)), c("war_winner", "war_loser", "war_winner:post_wwii", "war_loser:post_wwii"))) {
      cat(sprintf("    %s: %.6f (p=%.4f)\n", v, coef(m)[v], s$coeftable[v, "Pr(>|t|)"]))
    }
  }

  # Delta CINC split
  cat("  Delta CINC (war) → eigenvector (pre/post 3yr windows):\n")
  for (period in c("pre", "post")) {
    m <- res$prepost_models[[paste0("ols_", period, "_eigenvector_diff")]]
    if (!is.null(m)) {
      s <- summary(m)
      cat(sprintf("    %s-WWII: beta=%.4f, p=%.4f\n", toupper(period),
                  coef(m)["delta_cinc_war"], s$coefficients["delta_cinc_war", "Pr(>|t|)"]))
    }
  }

  # War IV × post_wwii interaction (eigenvector, mid_onset_count)
  cat("  War IV × post_wwii interaction (eigenvector):\n")
  for (iv in c("mid_onset_count", "fatal_mid_onset_count")) {
    m <- res$wariv_models[[paste0("interact_eigenvector_chg_", iv)]]
    if (!is.null(m)) {
      s <- summary(m)
      int_term <- paste0(iv, ":post_wwii")
      if (int_term %in% rownames(s$coeftable)) {
        cat(sprintf("    %s × post_wwii: %.6f (p=%.4f)\n", iv,
                    s$coeftable[int_term, "Estimate"], s$coeftable[int_term, "Pr(>|t|)"]))
      }
    }
  }

  # Winner/loser decomposition
  cat("  Delta CINC within winners/losers (eigenvector, contiguity):\n")
  for (grp in c("winner", "loser")) {
    m <- res$prepost_models[[paste0(grp, "_eigenvector_diff")]]
    if (!is.null(m) && "delta_cinc_war" %in% rownames(summary(m)$coefficients)) {
      s <- summary(m)
      cat(sprintf("    %s: beta=%.4f, p=%.4f\n", tools::toTitleCase(grp),
                  coef(m)["delta_cinc_war"], s$coefficients["delta_cinc_war", "Pr(>|t|)"]))
    }
  }
}

# ============================================================
# 5. Save
# ============================================================

cat("\n\nSaving results...\n")
saveRDS(all_results, file.path(OUTPUT_DIR, "winner_loser_results.rds"))

# CSV: Extract key coefficients from winner/loser TWFE models
wl_coef_list <- list()
for (wt in names(all_results)) {
  res <- all_results[[wt]]
  for (mn in names(res$wl_models)) {
    m <- res$wl_models[[mn]]
    if (is.null(m)) next
    ct <- tryCatch(coeftable(m), error = function(e) NULL)
    if (is.null(ct)) next
    wl_coef_list[[length(wl_coef_list) + 1]] <- tibble(
      weighting = wt, model = mn,
      term = rownames(ct), estimate = ct[, "Estimate"],
      se = ct[, "Std. Error"], pval = ct[, "Pr(>|t|)"]
    )
  }
}
if (length(wl_coef_list) > 0) {
  write_csv(bind_rows(wl_coef_list), file.path(OUTPUT_DIR, "winner_loser_twfe_coefs.csv"))
}

# CSV: War IV × post_wwii interaction coefficients
wariv_coef_list <- list()
for (wt in names(all_results)) {
  res <- all_results[[wt]]
  if (is.null(res$wariv_models)) next
  for (mn in names(res$wariv_models)) {
    m <- res$wariv_models[[mn]]
    if (is.null(m)) next
    ct <- tryCatch(coeftable(m), error = function(e) NULL)
    if (is.null(ct)) next
    wariv_coef_list[[length(wariv_coef_list) + 1]] <- tibble(
      weighting = wt, model = mn,
      term = rownames(ct), estimate = ct[, "Estimate"],
      se = ct[, "Std. Error"], pval = ct[, "Pr(>|t|)"]
    )
  }
}
if (length(wariv_coef_list) > 0) {
  write_csv(bind_rows(wariv_coef_list), file.path(OUTPUT_DIR, "winner_loser_wariv_interact_coefs.csv"))
  cat(sprintf("  Saved %d war IV interaction coefficients\n", nrow(bind_rows(wariv_coef_list))))
}

# CSV: Pre/post war centrality differences
wc_list <- list()
for (wt in names(all_results)) {
  wc <- all_results[[wt]]$war_cent_data
  if (!is.null(wc)) {
    wc_list[[wt]] <- wc |> mutate(weighting = wt)
  }
}
if (length(wc_list) > 0) {
  write_csv(bind_rows(wc_list), file.path(OUTPUT_DIR, "winner_loser_centrality_diffs.csv"))
}

# Plot: Winner vs Loser centrality change (eigenvector) — main weighting
library(ggplot2)
main_res <- all_results[["main_weight"]]
if (!is.null(main_res$war_cent_data)) {
  wc <- main_res$war_cent_data |>
    pivot_longer(cols = c(betweenness_diff, eigenvector_diff, pagerank_diff, degree_diff),
                 names_to = "centrality", values_to = "diff") |>
    mutate(
      outcome = case_when(war_winner == 1 ~ "Winner", war_loser == 1 ~ "Loser", TRUE ~ "Other"),
      era = ifelse(post_wwii == 1, "Post-WWII", "Pre-WWII"),
      centrality = gsub("_diff", "", centrality) |> str_to_title()
    )

  p_wl <- ggplot(wc, aes(x = outcome, y = diff, fill = era)) +
    geom_boxplot(outlier.size = 0.8, alpha = 0.7) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
    facet_wrap(~ centrality, scales = "free_y", ncol = 2,
               labeller = label_wrap_gen(width = 22)) +
    scale_fill_manual(values = c("Pre-WWII" = "steelblue", "Post-WWII" = "coral")) +
    labs(x = NULL, y = "Pre/Post War Centrality Change",
         title = "Winners vs losers",
         fill = "Era") +
    theme_minimal(base_size = 13) +
    theme(legend.position = "top",
          strip.text = element_text(face = "bold", size = 9),
          panel.spacing = unit(1.2, "lines"),
          panel.spacing.x = unit(1.5, "lines"),
          plot.margin = margin(t = 10, r = 12, b = 14, l = 10))

  ggsave(file.path(FIG_DIR, "fig_winner_loser_centrality.png"), p_wl,
         width = 10, height = 8, dpi = 300)
}

# LaTeX table: Winner/loser TWFE (main weighting only)
main_res <- all_results[["main_weight"]]
if (!is.null(main_res$wl_models)) {
  # Select the full models for the 4 centrality DVs
  tex_dvs <- c("full_eigenvector_chg", "full_betweenness_chg",
                "full_pagerank_chg", "full_degree_chg")
  tex_models <- compact(main_res$wl_models[tex_dvs])
  if (length(tex_models) > 0) {
    dv_names <- c(
      full_eigenvector_chg = "Eigenvector", full_betweenness_chg = "Betweenness",
      full_pagerank_chg = "PageRank", full_degree_chg = "Degree"
    )
    texreg(tex_models,
           file = file.path(APPX_DIR, "winner_loser_twfe.tex"),
           custom.model.names = dv_names[names(tex_models)],
           caption = "Winner/Loser TWFE: War Outcome and Centrality Change",
           label = "tab:winner_loser_twfe",
           float.pos = "htbp",
           booktabs = TRUE, use.packages = FALSE,
           stars = c(0.001, 0.01, 0.05, 0.1),
           symbol = "$^{\\dagger}$",
           custom.note = "DV: Pre/post-war change in centrality. FE: state. Clustered SEs (state). $^{***}p<0.001$; $^{**}p<0.01$; $^{*}p<0.05$; $^{\\dagger}p<0.1$.")
    cat("  LaTeX → results/appendix/winner_loser_twfe.tex\n")
  }
}

# LaTeX table: Winner/loser centrality differences summary
if (length(wl_coef_list) > 0) {
  wl_coefs <- bind_rows(wl_coef_list) |>
    filter(weighting == "main_weight",
           grepl("^full_", model),
           grepl("war_winner|war_loser|post_wwii", term))

  if (nrow(wl_coefs) > 0) {
    tex_wl <- c(
      "\\begin{table}[htbp]", "\\centering",
      "\\caption{Winner/Loser Key Coefficients (Main Weighting)}",
      "\\label{tab:winner_loser_coefs}",
      "\\small",
      "\\begin{tabular}{llrrr}", "\\toprule",
      "Model & Term & Estimate & SE & $p$ \\\\",
      "\\midrule"
    )
    for (i in seq_len(nrow(wl_coefs))) {
      r <- wl_coefs[i, ]
      tex_wl <- c(tex_wl, sprintf(
        "%s & %s & %.4f & %.4f & %.3f \\\\",
        r$model, r$term, r$estimate, r$se, r$pval
      ))
    }
    tex_wl <- c(tex_wl, "\\bottomrule", "\\end{tabular}", "\\end{table}")
    writeLines(tex_wl, file.path(APPX_DIR, "winner_loser_coefs.tex"))
    cat("  LaTeX → results/appendix/winner_loser_coefs.tex\n")
  }
}

cat("Done. Winner/loser analysis complete.\n")
