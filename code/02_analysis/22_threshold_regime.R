# [2026-09-26] Headline system-level DV switched from KL to Jensen-Shannon (JS)
# divergence. KL results are retained in results/appendix/kl_canonical/.
###############################################################################
# 06_threshold_regime.R
# Nonlinear models: threshold VAR (TVAR), Markov-switching, smooth-transition
# VAR (STVAR).
#
# Motivation: network reconfiguration may respond differently to war
# depending on the intensity of conflict (high vs. low war regimes)
# or the current state of the network (stable vs. transitional).
#
# TRANSFORMATION: mid_count and js_divergence are log1p-transformed before
# entering the TVAR / STVAR / Markov-switching models. Both series are
# non-negative and right-skewed; the log handles the skewness, makes
# residuals closer to normal, and lets regime-conditional coefficients
# read as approximate elasticities. Matches the transform applied in
# 05_svar.R and inside the VAR loop of 01_system_level_analysis.R.
# The main-text OLS (single-equation, LHS-only) retains raw scales.
###############################################################################

library(tidyverse)
library(tsDyn)     # Threshold VAR
library(MSwM)      # Markov-switching
select <- dplyr::select

OUTPUT_DIR <- here::here("data", "processed")
FIG_DIR <- here::here("code", "03_figures", "output")

# --------------------------------------------------------------------------- #
# 0. Load + log1p-transform data
# --------------------------------------------------------------------------- #

system_data <- readRDS(file.path(OUTPUT_DIR, "system_analysis_data.rds")) |>
  filter(!is.na(js_divergence)) |>
  arrange(year)

# Apply log1p to non-negative skewed series. After this point every
# reference to mid_count, fatal_mid_count, and js_divergence in the
# script refers to the LOGGED quantity. Plot labels and tex captions
# are updated accordingly.
.LOG_VARS <- c("mid_count", "fatal_mid_count", "js_divergence")
.log_targets <- intersect(.LOG_VARS, names(system_data))
system_data <- system_data |>
  mutate(across(all_of(.log_targets), ~ log1p(.x)))
cat(sprintf("Applied log1p to: %s\n",
            paste(.log_targets, collapse = ", ")))

# --------------------------------------------------------------------------- #
# Spec grid: 2 IVs x 2 controls applied to each model class.
# Controls match 01's full-controls VAR setup:
#   exogenous: log1p(n_democracies) + log1p(n_autocracies) + log1p(n_states)
# --------------------------------------------------------------------------- #

EXOG_FULL <- intersect(c("n_democracies", "n_autocracies", "n_states"),
                        names(system_data))

make_exog_matrix <- function(data) {
  if (length(EXOG_FULL) == 0) return(NULL)
  data |>
    select(all_of(EXOG_FULL)) |>
    mutate(across(everything(), log1p)) |>
    as.matrix()
}

# The 4 cells we hit for each model class.
spec_grid <- expand.grid(
  iv       = c("mid", "fatal"),
  controls = c("none", "full"),
  stringsAsFactors = FALSE
)

iv_col_lookup <- c(mid = "mid_count", fatal = "fatal_mid_count")

# --------------------------------------------------------------------------- #
# 1. Threshold VAR (TVAR) — 4 specs
# --------------------------------------------------------------------------- #

cat("--- Threshold VAR ---\n")

fit_tvar_one <- function(iv_col, controls) {
  blk <- system_data |>
    select(all_of(c(iv_col, "js_divergence", EXOG_FULL))) |>
    na.omit()
  endo <- blk |> select(all_of(c(iv_col, "js_divergence")))
  exog <- if (controls == "full") make_exog_matrix(blk) else NULL
  args <- list(data = endo, lag = 2, nthresh = 1,
               thDelay = 1, trim = 0.15, mTh = 1, plot = FALSE)
  if (!is.null(exog)) args$exogen <- exog
  tryCatch(do.call(TVAR, args),
           error = function(e) {
             cat(sprintf("  TVAR error [%s/%s]: %s\n",
                         iv_col, controls, e$message))
             NULL
           })
}

tvar_results <- list()
for (i in seq_len(nrow(spec_grid))) {
  sp <- spec_grid[i, ]
  iv_col <- iv_col_lookup[[sp$iv]]
  key <- sprintf("tvar__%s__%s", sp$iv, sp$controls)
  cat(sprintf("\n[%s]\n", key))
  fit <- fit_tvar_one(iv_col, sp$controls)
  if (!is.null(fit)) {
    cat("  TVAR estimated.\n")
    tvar_results[[key]] <- list(iv = sp$iv, controls = sp$controls,
                                  iv_col = iv_col, fit = fit)
  }
}

# Linearity test on the headline cell (MID, no controls) — same as before.
lin_test <- tryCatch({
  blk <- system_data |>
    select(mid_count, js_divergence) |> na.omit()
  TVAR.LRtest(blk, lag = 2, nthresh = 1, mTh = 1, plot = FALSE)
}, error = function(e) {
  cat(sprintf("  Linearity test error: %s\n", e$message)); NULL
})
if (!is.null(lin_test)) {
  cat("\nLinearity test (H0: linear VAR; headline cell, MID/no ctrl):\n")
  print(lin_test)
}

# --------------------------------------------------------------------------- #
# 2. Markov-switching — 4 specs
# --------------------------------------------------------------------------- #

cat("\n--- Markov-Switching Model ---\n")

ms_data <- system_data |>
  mutate(
    kl_lag1   = lag(js_divergence, 1),
    mid_lag1  = lag(mid_count, 1),
    fatal_lag1= lag(fatal_mid_count, 1),
    log_n_dem = log1p(.data[["n_democracies"]]),
    log_n_aut = log1p(.data[["n_autocracies"]]),
    log_n_st  = log1p(.data[["n_states"]])
  ) |>
  filter(!is.na(kl_lag1), !is.na(mid_lag1), !is.na(fatal_lag1))

fit_ms_one <- function(iv_col, lag_col, controls) {
  ctrl_terms <- if (controls == "full") "+ log_n_dem + log_n_aut + log_n_st" else ""
  fml <- as.formula(sprintf(
    "js_divergence ~ %s + kl_lag1 + %s%s",
    iv_col, lag_col, ctrl_terms))
  lm_base <- lm(fml, data = ms_data)
  # Switch all coefficients + variance.
  n_coef <- length(coef(lm_base))
  sw_full <- rep(TRUE, n_coef + 1)  # +1 for variance
  fit <- tryCatch(msmFit(lm_base, k = 2, sw = sw_full),
                   error = function(e) {
                     # Fallback: switch only intercept + IV + variance
                     sw_simple <- c(TRUE, TRUE, rep(FALSE, n_coef - 2), TRUE)
                     tryCatch(msmFit(lm_base, k = 2, sw = sw_simple),
                              error = function(e2) {
                                cat(sprintf("  MS fit failed: %s\n",
                                            e2$message)); NULL
                              })
                   })
  list(fit = fit, formula = fml, lm_base = lm_base)
}

ms_results <- list()
for (i in seq_len(nrow(spec_grid))) {
  sp <- spec_grid[i, ]
  iv_col  <- iv_col_lookup[[sp$iv]]
  lag_col <- if (sp$iv == "mid") "mid_lag1" else "fatal_lag1"
  key <- sprintf("ms__%s__%s", sp$iv, sp$controls)
  cat(sprintf("\n[%s]\n", key))
  res <- fit_ms_one(iv_col, lag_col, sp$controls)
  if (!is.null(res$fit)) {
    cat("  Markov-switching estimated.\n")
    ms_results[[key]] <- list(
      iv = sp$iv, controls = sp$controls, iv_col = iv_col,
      lag_col = lag_col, fit = res$fit, formula = res$formula,
      lm_base = res$lm_base)
  }
}

# Plot regime probabilities for the headline cell (MID, no controls).
ms_headline <- ms_results[["ms__mid__none"]]
if (!is.null(ms_headline)) {
  rp <- ms_headline$fit@Fit@smoProb
  n_use <- min(nrow(rp), nrow(ms_data))
  ms_plot_df <- ms_data[seq_len(n_use), ]
  ms_plot_df$regime1_prob <- rp[seq_len(n_use), 1]
  ms_plot_df$regime2_prob <- rp[seq_len(n_use), 2]

  p_ms <- ggplot(ms_plot_df, aes(x = year)) +
    geom_area(aes(y = regime2_prob), fill = "darkred", alpha = 0.4) +
    geom_line(aes(y = js_divergence / max(js_divergence, na.rm = TRUE)),
              linewidth = 0.5, alpha = 0.6) +
    geom_vline(xintercept = c(1914, 1939, 1950, 1991),
               linetype = "dashed", alpha = 0.4) +
    labs(x = "Year",
         y = "Regime 2 Probability (shaded) / Normalized log(1+JS) (line)",
         title = "Markov-switching regimes") +
    theme_minimal(base_size = 12) +
    theme(plot.title = element_text(face = "bold"),
          plot.margin = margin(t = 10, r = 12, b = 14, l = 10))
  ggsave(file.path(FIG_DIR, "fig_markov_switching.png"),
         p_ms, width = 10, height = 5, dpi = 300, bg = "white")
}

# --------------------------------------------------------------------------- #
# 3. Smooth Transition VAR (STVAR) — 4 specs
# --------------------------------------------------------------------------- #

cat("\n--- Smooth Transition Model ---\n")

fit_stvar_one <- function(iv_col, controls) {
  blk <- system_data |>
    select(all_of(c(iv_col, "js_divergence", EXOG_FULL))) |>
    na.omit()
  endo <- blk |> select(all_of(c(iv_col, "js_divergence")))
  exog <- if (controls == "full") make_exog_matrix(blk) else NULL
  args <- list(data = endo, lag = 2, nthresh = 1,
               thDelay = 1, mTh = 1, model = "MTAR",
               trim = 0.15, plot = FALSE)
  if (!is.null(exog)) args$exogen <- exog
  tryCatch(do.call(TVAR, args),
           error = function(e) {
             cat(sprintf("  STVAR error [%s/%s]: %s\n",
                         iv_col, controls, e$message))
             NULL
           })
}

stvar_results <- list()
for (i in seq_len(nrow(spec_grid))) {
  sp <- spec_grid[i, ]
  iv_col <- iv_col_lookup[[sp$iv]]
  key <- sprintf("stvar__%s__%s", sp$iv, sp$controls)
  cat(sprintf("\n[%s]\n", key))
  fit <- fit_stvar_one(iv_col, sp$controls)
  if (!is.null(fit)) {
    cat("  STVAR estimated.\n")
    stvar_results[[key]] <- list(iv = sp$iv, controls = sp$controls,
                                   iv_col = iv_col, fit = fit)
  }
}

# Backward-compat aliases that downstream code may reference.
tvar_test  <- if ("tvar__mid__none" %in% names(tvar_results))
                tvar_results[["tvar__mid__none"]]$fit else NULL
ms_fit     <- if ("ms__mid__none"   %in% names(ms_results))
                ms_results[["ms__mid__none"]]$fit else NULL
stvar_fit  <- if ("stvar__mid__none" %in% names(stvar_results))
                stvar_results[["stvar__mid__none"]]$fit else NULL

# --------------------------------------------------------------------------- #
# 4. Quantile VAR
# --------------------------------------------------------------------------- #

cat("\n--- Quantile Regression ---\n")

library(quantreg)

quantiles <- c(0.10, 0.25, 0.50, 0.75, 0.90)

fit_qr_one <- function(iv_col, lag_col, controls, tau) {
  ctrl_terms <- if (controls == "full") "+ log_n_dem + log_n_aut + log_n_st" else ""
  fml <- as.formula(sprintf(
    "js_divergence ~ %s + kl_lag1 + %s%s",
    iv_col, lag_col, ctrl_terms))
  fit <- rq(fml, data = ms_data, tau = tau)
  s   <- summary(fit, se = "boot")
  ct  <- s$coefficients
  if (!iv_col %in% rownames(ct)) {
    return(tibble(iv = NA, controls = NA, tau = NA,
                  iv_coef = NA, iv_se = NA, iv_pval = NA))
  }
  tibble(
    iv_coef = ct[iv_col, "Value"],
    iv_se   = ct[iv_col, "Std. Error"],
    iv_pval = ct[iv_col, "Pr(>|t|)"]
  )
}

qr_grid <- expand.grid(
  iv       = c("mid", "fatal"),
  controls = c("none", "full"),
  tau      = quantiles,
  stringsAsFactors = FALSE
)

qr_rows <- list()
for (i in seq_len(nrow(qr_grid))) {
  sp <- qr_grid[i, ]
  iv_col  <- iv_col_lookup[[sp$iv]]
  lag_col <- if (sp$iv == "mid") "mid_lag1" else "fatal_lag1"
  res <- fit_qr_one(iv_col, lag_col, sp$controls, sp$tau)
  qr_rows[[length(qr_rows) + 1]] <- tibble(
    iv = sp$iv, controls = sp$controls, tau = sp$tau,
    iv_coef = res$iv_coef, iv_se = res$iv_se, iv_pval = res$iv_pval
  )
  cat(sprintf("  [qr %s/%s tau=%.2f] coef=%.4f (p=%.4f)\n",
              sp$iv, sp$controls, sp$tau, res$iv_coef, res$iv_pval))
}

qr_summary <- bind_rows(qr_rows)

# Plot quantile regression coefficients (split by controls, faceted by IV)
qr_plot_df <- qr_summary |>
  mutate(controls_label = factor(controls, levels = c("none", "full"),
                                   labels = c("No controls", "Full controls")),
         iv_label = factor(iv, levels = c("mid", "fatal"),
                            labels = c("log(1+MID Count)",
                                        "log(1+Fatal MID Count)")))

p_qr <- ggplot(qr_plot_df,
                aes(x = tau, y = iv_coef,
                    color = controls_label)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = iv_coef - 1.96 * iv_se,
                    ymax = iv_coef + 1.96 * iv_se),
                width = 0.03) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_color_manual(values = c("No controls"   = "#08519c",
                                  "Full controls" = "#B03366"),
                     name = NULL) +
  facet_wrap(~ iv_label, scales = "free_y",
             labeller = label_wrap_gen(width = 22)) +
  labs(x = "Quantile (tau)", y = "Coefficient",
       title = "Quantile regression") +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold"),
        legend.position = "bottom",
        strip.text = element_text(size = 9),
        panel.spacing = unit(1.2, "lines"),
        panel.spacing.x = unit(1.5, "lines"),
        plot.margin = margin(t = 10, r = 12, b = 14, l = 10))

ggsave(file.path(FIG_DIR, "fig_quantile_regression.png"),
       p_qr, width = 10, height = 5, dpi = 300, bg = "white")

# --------------------------------------------------------------------------- #
# 5. Save
# --------------------------------------------------------------------------- #

TABLE_DIR <- here::here("results", "main")
APPX_DIR  <- here::here("results", "appendix")
dir.create(TABLE_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(APPX_DIR,  recursive = TRUE, showWarnings = FALSE)

nonlinear_results <- list(
  # 4-spec grids: 2 IVs x 2 controls each
  tvar      = tvar_results,
  ms        = ms_results,
  stvar     = stvar_results,
  quantile  = qr_summary,
  # Backward-compat aliases (headline = MID, no controls)
  tvar_headline  = tvar_test,
  ms_headline    = ms_fit,
  stvar_headline = stvar_fit
)

saveRDS(nonlinear_results, file.path(OUTPUT_DIR, "nonlinear_results.rds"))
write_csv(qr_summary, file.path(OUTPUT_DIR, "quantile_regression_summary.csv"))

# --- Coefficient extraction for TVAR / STVAR (regime-dependent) -------------
# tsDyn::TVAR returns regime-specific coefficient blocks. Pull the IV
# coefficient for each regime in each spec for downstream tables.
extract_tvar_coefs <- function(tvar_obj, iv_col) {
  if (is.null(tvar_obj)) return(NULL)
  cf <- tryCatch(coef(tvar_obj), error = function(e) NULL)
  if (is.null(cf) || !is.matrix(cf)) return(NULL)
  # Look for IV columns in coef matrix names. tsDyn names them like
  # "mid_count -1" for lag 1 in regime 1 and so on.
  out <- list()
  rn <- rownames(cf); cn <- colnames(cf)
  for (eq_row in seq_len(nrow(cf))) {
    eq_name <- rn[eq_row]
    for (col in seq_len(ncol(cf))) {
      term <- cn[col]
      if (grepl(iv_col, term, fixed = TRUE)) {
        out[[length(out) + 1]] <- tibble(
          equation = eq_name, term = term,
          estimate = cf[eq_row, col]
        )
      }
    }
  }
  if (length(out) == 0) return(NULL)
  bind_rows(out)
}

tvar_coef_rows <- list()
for (key in names(tvar_results)) {
  r <- tvar_results[[key]]
  e <- extract_tvar_coefs(r$fit, r$iv_col)
  if (!is.null(e)) {
    e <- e |> mutate(spec_key = key, iv = r$iv,
                      controls = r$controls, model = "tvar")
    tvar_coef_rows[[key]] <- e
  }
}
stvar_coef_rows <- list()
for (key in names(stvar_results)) {
  r <- stvar_results[[key]]
  e <- extract_tvar_coefs(r$fit, r$iv_col)
  if (!is.null(e)) {
    e <- e |> mutate(spec_key = key, iv = r$iv,
                      controls = r$controls, model = "stvar")
    stvar_coef_rows[[key]] <- e
  }
}

tvar_coefs_df  <- bind_rows(tvar_coef_rows)
stvar_coefs_df <- bind_rows(stvar_coef_rows)

if (nrow(tvar_coefs_df) > 0) {
  write_csv(tvar_coefs_df, file.path(OUTPUT_DIR,
            "threshold_tvar_coefs.csv"))
}
if (nrow(stvar_coefs_df) > 0) {
  write_csv(stvar_coefs_df, file.path(OUTPUT_DIR,
            "threshold_stvar_coefs.csv"))
}

# --- LaTeX: Quantile regression (combined IV x controls table) -------------
qr_iv_label <- c(mid = "$\\log(1{+}\\text{MID})$",
                  fatal = "$\\log(1{+}\\text{Fatal MID})$")
qr_ctrl_lbl <- c(none = "no ctrls", full = "full ctrls")

tex_qr <- c(
  "\\begin{table}[htbp]", "\\centering",
  "\\caption{Quantile Regression: Conflict Effect on $\\log(1+\\text{JS Divergence})$ by Quantile.}",
  "\\label{tab:quantile_reg}",
  "\\small",
  "\\begin{tabular}{llrrrr}", "\\toprule",
  "IV & Controls & Quantile ($\\tau$) & Coefficient & SE & $p$-value \\\\",
  "\\midrule"
)
for (i in seq_len(nrow(qr_summary))) {
  r <- qr_summary[i, ]
  tex_qr <- c(tex_qr, sprintf(
    "%s & %s & %.2f & %.4f & %.4f & %.3f \\\\",
    qr_iv_label[r$iv], qr_ctrl_lbl[r$controls],
    r$tau, r$iv_coef, r$iv_se, r$iv_pval
  ))
}
tex_qr <- c(tex_qr, "\\bottomrule", "\\end{tabular}",
  "\\begin{tablenotes}\\small",
  "\\item DV: $\\log(1+\\text{JS Divergence})$. Lagged DV and lagged IV always included. ",
  "Full controls add $\\log(1+\\text{n\\_dem})$, $\\log(1+\\text{n\\_aut})$, $\\log(1+\\text{n\\_states})$. ",
  "SEs bootstrapped.",
  "\\end{tablenotes}",
  "\\end{table}")
writeLines(tex_qr, file.path(APPX_DIR, "quantile_regression_summary.tex"))
cat("  LaTeX (quantile reg) -> quantile_regression_summary.tex\n")

# --- LaTeX: Markov-switching (one table per spec, headline goes to main) ----
# MSwM stores @Coef and @seCoef as data.frames (rows = variables, cols =
# regimes). Earlier check used `is.matrix()` which rejects data.frames --
# resulting in the "empty coef" skip even though coefs were present.
# Coerce to matrix and compute t-stats / p-values manually.
write_ms_table <- function(spec_key, file_path, caption, label,
                            transform = "logged") {
  res <- ms_results[[spec_key]]
  if (is.null(res) || is.null(res$fit)) {
    cat(sprintf("  Skipping MS table %s (no fit)\n", spec_key)); return(NULL)
  }
  fit <- res$fit
  cf <- tryCatch(as.matrix(fit@Coef),   error = function(e) NULL)
  se <- tryCatch(as.matrix(fit@seCoef), error = function(e) NULL)
  if (is.null(cf) || is.null(se) ||
      !all(dim(cf) == dim(se)) || nrow(cf) == 0) {
    cat(sprintf("  Skipping MS table %s (coef/se extraction failed)\n",
                spec_key))
    return(NULL)
  }
  # Compute t-stats and p-values from coef/se ratios.
  tstat <- cf / se
  pval  <- 2 * (1 - pnorm(abs(tstat)))

  var_labels <- c(
    "(Intercept)"        = "Intercept",
    mid_count            = "$\\log(1{+}\\text{MID Count})$",
    fatal_mid_count      = "$\\log(1{+}\\text{Fatal MID Count})$",
    kl_lag1              = "$\\log(1{+}\\text{JS})$ ($t{-}1$)",
    mid_lag1             = "$\\log(1{+}\\text{MID})$ ($t{-}1$)",
    fatal_lag1           = "$\\log(1{+}\\text{Fatal MID})$ ($t{-}1$)",
    log_n_dem            = "$\\log(1{+}\\text{Democracies})$",
    log_n_aut            = "$\\log(1{+}\\text{Autocracies})$",
    log_n_st             = "$\\log(1{+}\\text{States})$"
  )
  tex <- c(
    "\\begin{table}[htbp]", "\\centering",
    sprintf("\\caption{%s}", caption),
    sprintf("\\label{%s}", label),
    "\\small",
    "\\begin{tabular}{lrrrr}", "\\toprule",
    " & \\multicolumn{2}{c}{Regime 1 (Low)} & \\multicolumn{2}{c}{Regime 2 (High)} \\\\",
    "\\cmidrule(lr){2-3} \\cmidrule(lr){4-5}",
    "Variable & Coef. & $p$ & Coef. & $p$ \\\\",
    "\\midrule"
  )
  # Use POSITION-based indexing rather than name-based to avoid
  # "subscript out of bounds" when MSwM strips or alters rownames in
  # @Coef across versions. Variable label lookup still uses the
  # row name when available.
  var_names <- rownames(cf)
  n_vars <- nrow(cf)
  n_reg  <- ncol(cf)
  for (i in seq_len(n_vars)) {
    v <- if (!is.null(var_names) && length(var_names) >= i) var_names[i] else as.character(i)
    v_lab   <- ifelse(v %in% names(var_labels), var_labels[v], v)
    r1_coef <- cf[i, 1]
    r1_p    <- pval[i, 1]
    r2_coef <- if (n_reg >= 2) cf[i, 2]   else NA_real_
    r2_p    <- if (n_reg >= 2) pval[i, 2] else NA_real_
    tex <- c(tex, sprintf(
      "%s & %.4f & %.3f & %.4f & %.3f \\\\",
      v_lab, r1_coef, r1_p, r2_coef, r2_p))
  }
  tex <- c(tex, "\\bottomrule", "\\end{tabular}",
    "\\begin{tablenotes}\\small",
    sprintf("\\item Two-regime Markov-switching model on $\\log(1+\\text{JS})$. Spec: %s.",
            spec_key),
    "\\end{tablenotes}",
    "\\end{table}")
  writeLines(tex, file_path)
  cat(sprintf("  LaTeX (MS) [%s] -> %s\n", spec_key, basename(file_path)))
}

# Headline (MID, no ctrl) to main; the rest to appendix.
write_ms_table("ms__mid__none",
                file.path(TABLE_DIR, "table_markov_switching_mid_none.tex"),
                "Markov-Switching: log(1+JS) on log(1+MID), no controls.",
                "tab:ms_mid_none")
write_ms_table("ms__mid__full",
                file.path(APPX_DIR,  "tableA_markov_switching_mid_full.tex"),
                "Appendix. Markov-Switching: log(1+JS) on log(1+MID), full controls.",
                "tab:ms_mid_full")
write_ms_table("ms__fatal__none",
                file.path(APPX_DIR,  "tableA_markov_switching_fatal_none.tex"),
                "Appendix. Markov-Switching: log(1+JS) on log(1+Fatal MID), no controls.",
                "tab:ms_fatal_none")
write_ms_table("ms__fatal__full",
                file.path(APPX_DIR,  "tableA_markov_switching_fatal_full.tex"),
                "Appendix. Markov-Switching: log(1+JS) on log(1+Fatal MID), full controls.",
                "tab:ms_fatal_full")

cat("\nDone. Nonlinear model analysis complete.\n")
