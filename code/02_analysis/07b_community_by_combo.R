###############################################################################
# 07b_community_by_combo.R   [2026-11-09]
#
# Community fixed effects and community clustering with communities detected on EACH layer combination.
#
# Why: in 07 the state timing x FE x clustering grid ("timing_fe") uses one canonical community for every
#   layer combination and weighting: leiden_comm from 02 (two-stage DynMux on the main network, offensive +
#   defensive alliances + IGOs). For a combination such as trade + IGO, the community FE then comes from a
#   different network than the churn outcome. This script re-detects communities on each layer combination
#   and re-estimates the timing_fe grid with them.
#
# What it does
#   1. Runs 07's set-up code (everything before the system-level loop) unchanged, read from
#      07_robustness_multiverse.R at run time, so data, controls, combos, belligerent drops and the
#      coverage windows are identical to the main multiverse. Only the checkpoint directory and the reset
#      variable are renamed (two text substitutions, each asserted to hit exactly once / twice).
#   2. Detects communities on each of the 12 layer combinations with 02's algorithm and settings:
#      edge weight = unweighted sum of the combination's layers divided by its yearly maximum, rounded to 4
#      digits (identical to main_weight for the main combination); Leiden CPM with gamma = mean edge weight;
#      adjacent-year Jaccard tracking, tracking = "filtered", repair_chains = TRUE, repair_min_similarity = 0.2;
#      seed 123. Years = the state-panel years inside the combination's coverage window.
#      CHECK: for the main combination the detected communities must equal 02's leiden_comm exactly (stops
#      otherwise). Communities do not depend on the weighting method (02's community is unweighted too).
#   3. Re-estimates the timing_fe grid (07 block A2, copied verbatim) with the combination's community in
#      place of the canonical one, for every weighting method. Because the grid holds the sample fixed across
#      the two FE structures, the state + year FE cells can also change (through the singleton rule).
#   4. Writes the results and a comparison with 07's canonical-community grid.
#
# Scope: belligerent rule "none" only by default (the only rule the heatmap rows (i)-(j) and the timing-grid
#   meta-regression use). COMM07B_RULES=all runs every rule in BELLIG_RULES.
#   Nothing in 07's main grid uses the community (main-grid clustering is ccode+year, year, ccode; the
#   comm_spill extra control is not in EXTRA_CTRLS), so the main grid is unaffected.
#
# Inputs : 07_robustness_multiverse.R and its inputs (data/processed/...);
#          data/processed/robustness_specifications_timing_fe.csv.gz from the main 07 run (comparison only)
# Outputs: data/processed/combo_communities_07b.rds          (communities, cached; reused when settings match)
#          data/processed/combo_communities_07b_summary.csv  (per combination: communities per year, agreement
#                                                              with the canonical community)
#          data/processed/robustness_specifications_timing_fe_combocomm.csv.gz
#          results/appendix/community_by_combo_comparison.csv (share significant-positive, canonical vs combo)
#          results/appendix/community_by_combo_comparison_by_layers.csv
#          code/03_figures/output/fig_multiverse_{hypotheses,meta_by_hyp}.png (regenerated; needs the 2026-11-09 07 patch)
#          data/processed/robustness_specifications_timing_fe_combocomm.source_md5 (md5 of 07's timing file it matches)
#          checkpoints: data/processed/.robustness_checkpoints_07b/ (COMM07B_RESET=1 clears them; the set-up still
#          prints "ROBUSTNESS_RESET=1" in that message). To redo detection, delete combo_communities_07b.rds; the grid
#          checkpoints carry a hash of the communities and are re-estimated automatically.
# After ANY rerun of 07 (which rewrites its timing file), rerun 07b; 07 then falls back to the canonical community
#   with a warning until 07b has been rerun. Do not export MV_TFE_COMMUNITY=combo globally (a full 07 run would stop).
# Run   : ROBUSTNESS_WORKERS=12 caffeinate -i Rscript code/02_analysis/07b_community_by_combo.R
#          (do not run while 07 is running; it reads 07's source file and output)
###############################################################################
# [2026-11-09] SUPERSEDED: 07_robustness_multiverse.R now estimates the combo-community timing grid itself (same code)
# and writes both versions to robustness_specifications_timing_fe.csv.gz (column community). Kept for reference only;
# it always stops (its comparison and figure steps assume the pre-2026-11-09 07). The header above describes the
# version validated on 2026-11-09.
stop("07b is superseded by 07_robustness_multiverse.R (2026-11-09), which now writes both community versions to the timing file")
suppressPackageStartupMessages({library(dplyr); library(igraph); library(dynamicmultiplex)})
setwd(here::here())
stopifnot(packageVersion("dynamicmultiplex") >= "1.4.0")
if (nzchar(Sys.getenv("MV_PERM_DRAW"))) stop("07b does not run in permutation mode; unset MV_PERM_DRAW")
T0_07B <- Sys.time()

# ---- 1. 07 set-up code -------------------------------------------------------------------------------------
SRC_07 <- here::here("code", "02_analysis", "07_robustness_multiverse.R")
src <- readLines(SRC_07, warn = FALSE)
end_i <- which(src == 'cat("\\n=== SYSTEM-LEVEL ROBUSTNESS ===\\n")')
stopifnot(length(end_i) == 1)
pre <- src[seq_len(end_i - 1)]
old_ck <- 'CHECKPOINT_DIR <- here::here("data", "processed", ".robustness_checkpoints")'
stopifnot(sum(pre == old_ck) == 1)
pre[pre == old_ck] <- 'CHECKPOINT_DIR <- here::here("data", "processed", ".robustness_checkpoints_07b")'
reset_i <- grep('Sys.getenv("ROBUSTNESS_RESET"', pre, fixed = TRUE)
stopifnot(length(reset_i) == 2)
pre[reset_i] <- gsub('ROBUSTNESS_RESET', 'COMM07B_RESET', pre[reset_i], fixed = TRUE)
SRC_07_MD5 <- unname(tools::md5sum(SRC_07))
cat(sprintf("07b: running 07 set-up (lines 1-%d of %s, md5 %s)\n", end_i - 1, basename(SRC_07), SRC_07_MD5))
eval(parse(text = pre), envir = globalenv())
stopifnot(!PERM_MODE, basename(CHECKPOINT_DIR) == ".robustness_checkpoints_07b",
          exists("all_combos"), exists("MAIN_LAYERS"), exists("bellig_drop"), exists("canon_churn_rule"),
          exists("tfe_clustering"), "leiden_comm" %in% names(state_data), !anyNA(state_data$leiden_comm))

CANON_FILE <- file.path(OUTPUT_DIR, "robustness_specifications_timing_fe.csv.gz")
if (!file.exists(CANON_FILE)) stop("07b: ", basename(CANON_FILE), " not found; run 07 first")
RULES_07B <- if (Sys.getenv("COMM07B_RULES", "none") == "all") BELLIG_RULES else "none"
stopifnot(all(RULES_07B %in% BELLIG_RULES))
TFE_DVS_07B <- intersect(state_dvs, STATE_TFE_DVS)
stopifnot(length(TFE_DVS_07B) > 0, STATE_TRANSFORMS == "raw")

combo_key <- function(layers) paste(gsub("layer_", "", sort(layers)), collapse = "+")
MAIN_KEY <- combo_key(MAIN_LAYERS)
stopifnot(MAIN_KEY %in% sapply(all_combos, combo_key))
COMM_SETTINGS <- paste("02 two-stage DynMux: unweighted sum / yearly max, round 4; leiden cpm gamma = mean weight;",
                       "adjacent links weight 1; filtered; repair_chains TRUE min 0.2; seed 123; dynamicmultiplex",
                       as.character(packageVersion("dynamicmultiplex")), "igraph", as.character(packageVersion("igraph")),
                       "| years = state-panel years in coverage window")
COMBO_WINDOWS <- setNames(lapply(all_combos, get_valid_year_range), sapply(all_combos, combo_key))
WINDOWS_STR <- paste(names(COMBO_WINDOWS), sapply(COMBO_WINDOWS, paste, collapse = "-"), sep = ":", collapse = ";")
SCHEMA_07B <- paste(MV_SCHEMA_VERSION, "| 07b combocomm", COMM_SETTINGS, "| windows", WINDOWS_STR,
                    "| rules", paste(RULES_07B, collapse = ","))

# ---- 2. Communities per layer combination ------------------------------------------------------------------
COMM_YEARS <- sort(unique(state_data$year))
detect_comm <- function(dsub, layers, years) {
  m <- as.matrix(dsub[, layers, drop = FALSE]); m[is.na(m)] <- 0
  raw <- rowSums(m)
  mx  <- ave(raw, dsub$year, FUN = function(x) max(x, na.rm = TRUE))
  dsub$w__ <- round(ifelse(mx > 0, raw / mx, 0), 4)            # = main_weight for the main combination
  layer_list <- list(); year_map <- integer(0)
  for (yr in years) {
    e <- dsub[dsub$year == yr & !is.na(dsub$w__), ] |>
      dplyr::transmute(a = pmin(ccode1, ccode2), b = pmax(ccode1, ccode2), w = w__) |>
      dplyr::group_by(a, b) |> dplyr::summarise(w = max(w), .groups = "drop")
    stopifnot(all(e$a != e$b), !anyDuplicated(e[c("a", "b")]))
    if (nrow(e) < 2) next
    g <- igraph::graph_from_data_frame(e |> dplyr::select(a, b), directed = FALSE)
    igraph::E(g)$weight <- e$w
    stopifnot(!any(igraph::which_multiple(g)), !any(igraph::which_loop(g)), igraph::ecount(g) == nrow(e))
    layer_list[[length(layer_list) + 1]] <- g; year_map <- c(year_map, yr)
  }
  links <- data.frame(from = seq_len(length(layer_list) - 1), to = seq_len(length(layer_list) - 1) + 1, weight = 1)
  gam <- mean(unlist(lapply(layer_list, function(g) igraph::E(g)$weight)), na.rm = TRUE)
  # 02 runs with R's default RNG; furrr workers use L'Ecuyer-CMRG, so set the default kind before seeding
  # (otherwise set.seed(123) would start a different stream) and restore the caller's kind afterwards.
  old_kind <- RNGkind(); on.exit(do.call(RNGkind, as.list(old_kind)), add = TRUE)
  RNGkind("Mersenne-Twister", "Inversion", "Rejection")
  set.seed(123)
  fit <- dynamicmultiplex::fit_multilayer_jaccard(layer_list, algorithm = "leiden", layer_links = links,
           resolution_parameter = gam, objective = "cpm", directed = FALSE, seed = 123,
           allow_unequal_nodes = TRUE, tracking = "filtered", repair_chains = TRUE, repair_min_similarity = 0.2)
  meta <- dynamicmultiplex::extract_meta_membership(fit)
  stopifnot(length(meta) == length(layer_list))
  out <- dplyr::bind_rows(lapply(seq_along(meta), function(i) {
    nm <- igraph::V(layer_list[[i]])$name
    stopifnot(length(nm) == length(meta[[i]]), !anyNA(meta[[i]]))
    data.frame(ccode = as.numeric(nm), year = year_map[i], leiden_comm = as.integer(meta[[i]]))
  }))
  stopifnot(!anyDuplicated(out[c("ccode", "year")]))
  attr(out, "gamma") <- gam
  out
}

COMM_FILE <- file.path(OUTPUT_DIR, "combo_communities_07b.rds")
dyad_md5 <- unname(tools::md5sum(file.path(OUTPUT_DIR, "dyad_year_network_layers.rds")))
cached <- if (file.exists(COMM_FILE)) readRDS(COMM_FILE) else NULL
if (!is.null(cached) && identical(cached$settings, COMM_SETTINGS) && identical(cached$dyad_md5, dyad_md5) &&
    identical(cached$years, COMM_YEARS) && identical(cached$windows, COMBO_WINDOWS) &&
    setequal(names(cached$comm), sapply(all_combos, combo_key))) {
  combo_comm <- cached$comm
  cat(sprintf("07b: communities loaded from %s (%d combinations)\n", basename(COMM_FILE), length(combo_comm)))
} else {
  cat(sprintf("07b: detecting communities on %d layer combinations (%d workers)...\n", length(all_combos), .n_workers))
  t_c <- Sys.time()
  combo_comm <- furrr::future_map(all_combos, function(layers) {
    yr <- get_valid_year_range(layers)
    yrs <- COMM_YEARS[COMM_YEARS >= yr[1] & COMM_YEARS <= yr[2]]
    dsub <- dyad_data[dyad_data$year %in% yrs, c("ccode1", "ccode2", "year", layers)]
    detect_comm(dsub, layers, yrs)
  }, .options = furrr::furrr_options(seed = TRUE, globals = c("detect_comm", "dyad_data", "COMM_YEARS",
                                                                "get_valid_year_range", "layer_year_ranges")))
  names(combo_comm) <- sapply(all_combos, combo_key)
  cat(sprintf("07b: communities done in %.1f min\n", as.numeric(difftime(Sys.time(), t_c, units = "mins"))))
  saveRDS(list(comm = combo_comm, settings = COMM_SETTINGS, dyad_md5 = dyad_md5, years = COMM_YEARS, windows = COMBO_WINDOWS,
               src_07_md5 = SRC_07_MD5, created = Sys.time()), COMM_FILE)
}

# CHECK: the main combination reproduces 02's canonical community exactly.
chk <- state_data |> dplyr::select(ccode, year, canon = leiden_comm) |>
  dplyr::left_join(combo_comm[[MAIN_KEY]] |> dplyr::rename(combo = leiden_comm), by = c("ccode", "year"))
if (anyNA(chk$combo) || !all(chk$canon == chk$combo))
  stop(sprintf(paste("07b: main combination (%s) does not reproduce 02's leiden_comm (%d of %d state-years differ or missing).",
                     "Likely igraph/dynamicmultiplex version drift or changed inputs since 02 ran; rerun 02 or check versions."),
               MAIN_KEY, sum(is.na(chk$combo) | chk$canon != chk$combo), nrow(chk)))
cat(sprintf("07b: CHECK passed, main combination %s reproduces leiden_comm on all %d state-years\n", MAIN_KEY, nrow(chk)))

# Every state-year inside a combination's coverage window has a community.
for (cb in all_combos) {
  yr <- get_valid_year_range(cb)
  need <- state_data |> dplyr::filter(year >= yr[1], year <= yr[2]) |> dplyr::select(ccode, year)
  miss <- dplyr::anti_join(need, combo_comm[[combo_key(cb)]], by = c("ccode", "year"))
  if (nrow(miss) > 0) stop(sprintf("07b: %s has no community for %d state-years in its window", combo_key(cb), nrow(miss)))
}

# Summary per combination (agreement with the canonical community: mean NMI and ARI across years).
comm_summary <- dplyr::bind_rows(lapply(names(combo_comm), function(k) {
  cc <- combo_comm[[k]] |> dplyr::inner_join(state_data |> dplyr::select(ccode, year, canon = leiden_comm),
                                             by = c("ccode", "year"))
  by_yr <- split(cc, cc$year)
  nmi <- sapply(by_yr, function(z) igraph::compare(z$leiden_comm, z$canon, method = "nmi"))
  ari <- sapply(by_yr, function(z) igraph::compare(z$leiden_comm, z$canon, method = "adjusted.rand"))
  per_yr <- sapply(by_yr, function(z) dplyr::n_distinct(z$leiden_comm))
  sizes <- table(paste(cc$year, cc$leiden_comm))
  data.frame(layers = k, gamma = attr(combo_comm[[k]], "gamma"), year_start = min(cc$year), year_end = max(cc$year),
             state_years = nrow(cc), meta_ids = dplyr::n_distinct(cc$leiden_comm),
             median_comm_per_year = stats::median(per_yr), share_singleton_state_years = mean(sizes[paste(cc$year, cc$leiden_comm)] == 1),
             mean_nmi_vs_canonical = mean(nmi), mean_ari_vs_canonical = mean(ari, na.rm = TRUE),
             years_ari_undefined = sum(is.na(ari)))   # ARI is NaN when both partitions have one community
}))
readr::write_csv(comm_summary, file.path(OUTPUT_DIR, "combo_communities_07b_summary.csv"))
print(comm_summary, digits = 3)

# ---- 3. timing_fe grid with the combination's community ----------------------------------------------------
# Checkpoints are tied to the detected communities too, so re-detected communities always re-estimate the grid.
SCHEMA_07B <- paste(SCHEMA_07B, "| communities", rlang::hash(lapply(combo_comm, as.data.frame)))
read_ckpt_07b <- function(path) {
  x <- readRDS(path)
  if (!identical(attr(x, "schema_07b"), SCHEMA_07B)) stop("07b checkpoint ", basename(path), " has a different schema; rerun with COMM07B_RESET=1")
  x
}
for (tfm in STATE_TRANSFORMS) {
  dyad_work_state <- dyad_data   # raw
  for (method in weighting_methods) {
    ckpt_path <- file.path(CHECKPOINT_DIR, sprintf("state07b__%s__%s.rds", tfm, method))
    if (file.exists(ckpt_path)) {
      cat(sprintf("  [07b %s/%s] LOADED CHECKPOINT (%d rows)\n", tfm, method, nrow(read_ckpt_07b(ckpt_path)))); next
    }
    cat(sprintf("  [07b %s/%s] dispatching %d combos (%d workers)...\n", tfm, method, length(all_combos), .n_workers))
    t_m <- Sys.time()
    res_m <- furrr::future_map_dfr(seq_along(all_combos), function(combo_idx) {
      layers_subset <- all_combos[[combo_idx]]
      combo_name <- paste(gsub("layer_", "", layers_subset), collapse = "+")   # as in 07
      comm_cb <- combo_comm[[combo_key(layers_subset)]] |> dplyr::select(ccode, year, leiden_comm)
      manhattan_vec <- compute_manhattan_vec(dyad_work_state, layers_subset, method)
      churn_prep <- prep_state_churn(dyad_work_state, manhattan_vec)
      dplyr::bind_rows(lapply(RULES_07B, function(bellig_rule) {
        out_list <- list()
        state_dvs_df <- churn_from_prep(churn_prep, drop = bellig_drop[[bellig_rule]])
        if (nrow(state_dvs_df) < 100) return(dplyr::bind_rows(out_list))
        # Same merge as 07 (columns the timing grid uses), with the combination's community in place of 02's.
        churn_merged <- state_dvs_df |>
          dplyr::left_join(state_data |>
                             dplyr::select(ccode, year, churn, log_churn, mid_count, fatal_mid_count,
                                           mid_onset_count, fatal_mid_onset_count, log_mid_count, log_fatal_mid_count,
                                           log_mid_onset_count, log_fatal_mid_onset_count, is_isolated_t1,
                                           democracy, autocracy, cinc, major_power,
                                           dplyr::all_of(paste0(state_war_vars, "__lag1"))) |> dplyr::distinct(),
                           by = c("ccode", "year")) |>
          dplyr::left_join(comm_cb, by = c("ccode", "year"))
        if (bellig_rule != "none") {
          churn_merged <- churn_merged |> dplyr::select(-churn, -log_churn) |>
            dplyr::left_join(canon_churn_rule[[bellig_rule]], by = c("ccode", "year"))
        }
        stopifnot(!anyDuplicated(churn_merged[c("ccode", "year")]))
        has_offensive <- "layer_offensive" %in% layers_subset; has_defensive <- "layer_defensive" %in% layers_subset
        has_igo <- "layer_igo" %in% layers_subset; has_trade <- "layer_trade" %in% layers_subset
        has_dca <- "layer_dca" %in% layers_subset
        yr_range_st <- get_valid_year_range(layers_subset); year_start <- yr_range_st[1]; year_end <- yr_range_st[2]
        churn_combo <- churn_merged |>
          dplyr::filter(year >= year_start, year <= year_end) |>
          dplyr::mutate(churn_alt     = dplyr::if_else(year == year_start, NA_real_, churn_alt),
                        log_churn_alt = dplyr::if_else(year == year_start, NA_real_, log_churn_alt))
        stopifnot(!anyNA(churn_combo$leiden_comm))   # every state-year in the window has a combo community (checked above)

        for (war_var in state_war_vars) {
          iv_main <- paste0(war_var, "__lag1")
          is_onset <- grepl("onset", war_var)
          for (dv in TFE_DVS_07B) {
            dv_col <- dv
            analysis_sub <- churn_combo |> dplyr::filter(!is.na(.data[[dv_col]]), !is.na(.data[[iv_main]]))
            if (nrow(analysis_sub) < 100) next
            iv_control_sets <- get_control_sets(war_var)
            for (ctrl_name in names(iv_control_sets)) {
              controls_rhs <- iv_control_sets[[ctrl_name]]
              exposure_label <- ctrl_name
              # ---- 07 block A2 (verbatim from 07, 2026-11-08 version) ----
              lag_col  <- paste0(war_var, "__lag1")
              tfe_data <- analysis_sub |>
                dplyr::filter(!is.na(.data[[war_var]]), !is.na(.data[[lag_col]]), !is.na(leiden_comm))
              tfe_data$iv_diff__ <- tfe_data[[lag_col]] - tfe_data[[war_var]]
              for (it in 1:10) {
                n_before <- nrow(tfe_data)
                if (n_before < 100) break
                for (fe_rhs_s in c("ccode + year + leiden_comm", "ccode + year")) {
                  fit_samp <- tryCatch(feols(as.formula(paste(dv_col, "~", war_var, "+", lag_col, "+", controls_rhs, "|", fe_rhs_s)),
                                             data = tfe_data, vcov = "iid"), error = function(e) NULL)
                  if (is.null(fit_samp)) {
                    warning(sprintf("07b timing_fe: sample-fixing fit (%s) failed for %s / %s / %s; cell skipped",
                                    fe_rhs_s, combo_name, war_var, dv))
                    tfe_data <- tfe_data[0, ]; break
                  }
                  tfe_data <- tfe_data[obs(fit_samp), ]
                }
                if (nrow(tfe_data) == n_before) break
              }
              if (nrow(tfe_data) >= 100) {
                for (fe_struct in STATE_TFE_FE) {
                  fe_rhs <- if (fe_struct == "state_year") "ccode + year" else "ccode + year + leiden_comm"
                  fit_t  <- tryCatch(feols(as.formula(paste(dv_col, "~", war_var, "+", controls_rhs, "|", fe_rhs)),
                                           data = tfe_data, vcov = "iid"), error = function(e) NULL)
                  fit_t1 <- tryCatch(feols(as.formula(paste(dv_col, "~", lag_col, "+", controls_rhs, "|", fe_rhs)),
                                           data = tfe_data, vcov = "iid"), error = function(e) NULL)
                  fit_j  <- tryCatch(feols(as.formula(paste(dv_col, "~", war_var, "+", lag_col, "+", controls_rhs, "|", fe_rhs)),
                                           data = tfe_data, vcov = "iid"), error = function(e) NULL)
                  fit_s  <- tryCatch(feols(as.formula(paste(dv_col, "~", war_var, "+ iv_diff__ +", controls_rhs, "|", fe_rhs)),
                                           data = tfe_data, vcov = "iid"), error = function(e) NULL)
                  if (!is.null(fit_j) && !is.null(fit_s) && all(c(war_var, lag_col) %in% names(coef(fit_j))) &&
                      war_var %in% names(coef(fit_s))) {
                    if (!isTRUE(all.equal(unname(coef(fit_s)[[war_var]]), unname(sum(coef(fit_j)[c(war_var, lag_col)])),
                                          tolerance = 1e-5))) {
                      warning(sprintf("07b timing_fe: joint_sum != b_t + b_t1 for %s / %s / %s (%s); joint_sum row dropped",
                                      combo_name, war_var, dv, fe_struct))
                      fit_s <- NULL
                    }
                  }
                  for (f_chk in list(fit_t, fit_t1, fit_j, fit_s)) {
                    if (!is.null(f_chk)) stopifnot(nobs(f_chk) == nrow(tfe_data))
                  }
                  spec_terms <- list()
                  if ("t" %in% STATE_TFE_TIMINGS)  spec_terms[[length(spec_terms) + 1]] <- list(timing = "t",  term = "t",  fit = fit_t,  var = war_var)
                  if ("t1" %in% STATE_TFE_TIMINGS) spec_terms[[length(spec_terms) + 1]] <- list(timing = "t1", term = "t1", fit = fit_t1, var = lag_col)
                  if ("joint" %in% STATE_TFE_TIMINGS) {
                    spec_terms[[length(spec_terms) + 1]] <- list(timing = "joint", term = "joint_t",   fit = fit_j, var = war_var)
                    spec_terms[[length(spec_terms) + 1]] <- list(timing = "joint", term = "joint_t1",  fit = fit_j, var = lag_col)
                    spec_terms[[length(spec_terms) + 1]] <- list(timing = "joint", term = "joint_sum", fit = fit_s, var = war_var)
                  }
                  for (st in spec_terms) {
                    if (is.null(st$fit)) next
                    for (clust_name in STATE_TFE_CLUST) {
                      ct <- tryCatch(coeftable(st$fit, vcov = tfe_clustering[[clust_name]]), error = function(e) NULL)
                      if (is.null(ct) || !st$var %in% rownames(ct)) next
                      out_list[[length(out_list) + 1]] <- tibble::tibble(
                        level = "state", dv = dv, war_variable = war_var,
                        weighting = method, transform = tfm, layers = combo_name,
                        n_layers = length(layers_subset),
                        model = paste0(if (fe_struct == "state_year") "FE_full_" else "FE_commFE_", clust_name),
                        extra_control = "none", exposure = exposure_label,
                        grid = "timing_fe", iv_timing = st$timing, iv_term = st$term,
                        fe_structure = fe_struct, clustering = clust_name,
                        coef = ct[st$var, "Estimate"], se = ct[st$var, "Std. Error"], pval = ct[st$var, "Pr(>|t|)"],
                        n_obs = nobs(st$fit),
                        has_offensive = has_offensive, has_defensive = has_defensive,
                        has_igo = has_igo, has_trade = has_trade, has_dca = has_dca,
                        is_onset = is_onset, year_start = year_start, year_end = year_end)
                    }
                  }
                }
              }
              # ---- end A2 ----
            }
          }
        }
        res_r <- dplyr::bind_rows(out_list)
        if (nrow(res_r) > 0) res_r$bellig_rule <- bellig_rule
        res_r
      }))
    }, .options = furrr::furrr_options(seed = TRUE))
    stopifnot(nrow(res_m) > 0)
    res_m$community <- "combo"
    attr(res_m, "schema_07b") <- SCHEMA_07B
    saveRDS(res_m, ckpt_path)
    cat(sprintf("    -> %d rows in %.1f min [checkpoint saved]\n", nrow(res_m),
                as.numeric(difftime(Sys.time(), t_m, units = "mins"))))
    rm(res_m); gc(verbose = FALSE)
  }
}

# ---- 4. Combine, write, compare ------------------------------------------------------------------------------
ck_files <- file.path(CHECKPOINT_DIR, sprintf("state07b__%s__%s.rds", rep(STATE_TRANSFORMS, each = length(weighting_methods)),
                                              rep(weighting_methods, times = length(STATE_TRANSFORMS))))
stopifnot(all(file.exists(ck_files)))
cc_df <- dplyr::bind_rows(lapply(ck_files, read_ckpt_07b))
stopifnot(nrow(cc_df) > 0, all(cc_df$bellig_rule %in% RULES_07B), all(cc_df$grid == "timing_fe"))
OUT_07B <- file.path(OUTPUT_DIR, "robustness_specifications_timing_fe_combocomm.csv.gz")
readr::write_csv(cc_df, OUT_07B)
# Sidecar: md5 of the canonical timing file this grid belongs to. 07's figure section uses the 07b grid only when
# this matches the current canonical file (content check, not file times).
writeLines(unname(tools::md5sum(CANON_FILE)), sub("\\.csv\\.gz$", ".source_md5", OUT_07B))
cat(sprintf("07b: %d rows -> %s (+ .source_md5)\n", nrow(cc_df), basename(OUT_07B)))

{
  KEYS <- c("dv", "war_variable", "weighting", "transform", "layers", "exposure", "iv_timing", "iv_term",
            "fe_structure", "clustering", "bellig_rule")
  canon <- readr::read_csv(CANON_FILE, show_col_types = FALSE) |>
    dplyr::filter(bellig_rule %in% RULES_07B, dplyr::coalesce(extra_control, "none") == "none") |>
    dplyr::select(dplyr::all_of(KEYS), coef_canon = coef, pval_canon = pval, n_canon = n_obs)
  combo <- cc_df |> dplyr::select(dplyr::all_of(KEYS), coef_combo = coef, pval_combo = pval, n_combo = n_obs)
  stopifnot(!anyDuplicated(canon[KEYS]), !anyDuplicated(combo[KEYS]))
  cmp <- dplyr::full_join(canon, combo, by = KEYS)
  cat(sprintf("07b: cells canonical %d, combo %d, matched %d, canonical only %d, combo only %d\n",
              nrow(canon), nrow(combo), sum(!is.na(cmp$coef_canon) & !is.na(cmp$coef_combo)),
              sum(is.na(cmp$coef_combo)), sum(is.na(cmp$coef_canon))))
  sp <- function(b, p) !is.na(b) & !is.na(p) & b > 0 & (p / 2) < 0.10   # heatmap rule
  cmp <- cmp |> dplyr::mutate(sig_canon = sp(coef_canon, pval_canon), sig_combo = sp(coef_combo, pval_combo),
                              main_combo = layers == MAIN_KEY)
  # Cells that do not involve the community at all (state+year FE, no community clustering) can still differ
  # through the sample-fixing rule; for the main combination nothing may differ.
  if (any(cmp$main_combo & (is.na(cmp$coef_canon) != is.na(cmp$coef_combo))))
    stop("07b: main-combination cells present on one side only (canonical vs combo)")
  m_main <- cmp |> dplyr::filter(main_combo, !is.na(coef_canon), !is.na(coef_combo))
  if (nrow(m_main) > 0 && (any(m_main$n_canon != m_main$n_combo) ||
                           !isTRUE(all.equal(m_main$coef_canon, m_main$coef_combo, tolerance = 1e-8))))
    stop("07b: main-combination cells differ from 07's canonical grid (should be identical)")
  summ <- cmp |> dplyr::filter(!is.na(coef_canon), !is.na(coef_combo),
                               iv_term %in% c("t", "t1", "joint_sum")) |>
    dplyr::group_by(fe_structure, clustering) |>
    dplyr::summarise(cells = dplyr::n(), share_sig_canonical = mean(sig_canon), share_sig_combo = mean(sig_combo),
                     agree = mean(sig_canon == sig_combo), n_changed = mean(n_canon != n_combo),
                     median_rel_coef = stats::median(coef_combo / coef_canon, na.rm = TRUE), .groups = "drop")
  summ_tot <- cmp |> dplyr::filter(!is.na(coef_canon), !is.na(coef_combo), iv_term %in% c("t", "t1", "joint_sum")) |>
    dplyr::summarise(fe_structure = "all", clustering = "all", cells = dplyr::n(), share_sig_canonical = mean(sig_canon),
                     share_sig_combo = mean(sig_combo), agree = mean(sig_canon == sig_combo),
                     n_changed = mean(n_canon != n_combo), median_rel_coef = stats::median(coef_combo / coef_canon, na.rm = TRUE))
  summ <- dplyr::bind_rows(summ, summ_tot)
  readr::write_csv(summ, file.path(APPX_DIR, "community_by_combo_comparison.csv"))
  # Per combination: separates "different community" from "smaller sample" (sample fixing drops singletons).
  by_combo <- cmp |> dplyr::filter(!is.na(coef_canon), !is.na(coef_combo), iv_term %in% c("t", "t1", "joint_sum")) |>
    dplyr::group_by(layers, fe_structure) |>
    dplyr::summarise(cells = dplyr::n(), share_sig_canonical = mean(sig_canon), share_sig_combo = mean(sig_combo),
                     share_n_changed = mean(n_canon != n_combo), median_n_ratio = stats::median(n_combo / n_canon),
                     .groups = "drop") |>
    dplyr::left_join(comm_summary |> dplyr::select(layers, share_singleton_state_years, mean_nmi_vs_canonical), by = "layers")
  readr::write_csv(by_combo, file.path(APPX_DIR, "community_by_combo_comparison_by_layers.csv"))
  cat("\n07b: by layer combination\n"); print(as.data.frame(by_combo), digits = 3)
  cat("\n07b: canonical vs per-combination community (heatmap terms t, t1, joint_sum; share significant-positive)\n")
  print(as.data.frame(summ), digits = 3)
}
# ---- 5. Regenerate the multiverse figures with the per-combination community ---------------------------------
# Runs 07's standalone figure section (it reads the saved CSVs only); rows (i)-(j) and the timing meta-regression
# then use the file written above (MV_TFE_COMMUNITY = combo makes 07 stop if it cannot).
FIG_MARK <- "# tmp_multiverse_appendix_fig.R   [TEMPORARY / standalone]"
src_now <- readLines(SRC_07, warn = FALSE)
fig_i <- which(src_now == FIG_MARK)
stopifnot(length(fig_i) == 1, any(grepl("TFE_COMBO_FILE", src_now[fig_i:length(src_now)], fixed = TRUE)))
if (!identical(unname(tools::md5sum(SRC_07)), SRC_07_MD5))
  cat("07b: note, 07_robustness_multiverse.R changed since this run started; figures use the current file\n")
Sys.setenv(MV_TFE_COMMUNITY = "combo")
cat("07b: regenerating fig_multiverse_hypotheses.png and fig_multiverse_meta_by_hyp.png...\n")
eval(parse(text = src_now[(fig_i - 1):length(src_now)]), envir = new.env(parent = globalenv()))
Sys.unsetenv("MV_TFE_COMMUNITY")

cat(sprintf("\n07b_community_by_combo.R complete (%.1f min)\n", as.numeric(difftime(Sys.time(), T0_07B, units = "mins"))))
