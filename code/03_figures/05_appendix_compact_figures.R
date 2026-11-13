###############################################################################
# 05_appendix_compact_figures.R   [2026-11-09]
#
# Redraws eight appendix figures as four combined figures to save space.
# Reads saved results only. No model is re-estimated.
#
#   B1  fig_app_compact_sweep_families.png
#       = fig_contamination_democracy + fig_contamination_capitalism
#       inputs: results/replication_sweep/est_<id>.rds, code/04_replication_sweep/registry.csv
#       data logic copied from 91_figures.R (family loop) unchanged: z, status with alpha_pub, ratios, clipping, axes.
#       presentation changes: shorter panel titles; "|z| < 1.96" label moved below the band; top-5 name labels
#       not drawn (91_figures.R already empties them for non-TIE families); figure note only on the bottom row;
#       smaller point/label sizes for a 9.5 in canvas; PNG only (no PDF); orphan-est check and panel-B n printout dropped.
#   B2  fig_app_compact_contam_spec.png
#       = fig_ms_contamination_specification + fig_app_contamination_nocontrols
#       inputs: data/processed/contamination_specification.csv, contamination_nocontrols_ladder.csv (21_*.R)
#   B3  fig_app_compact_war_samples.png
#       = log(1 + fatal MIDs) at t-1, system over state; system cut on fatal-MID years, state fatality threshold raised [2026-11-12]
#       inputs: data/processed/robustness_war_samples.csv (11_*.R)
#   B4  fig_app_compact_year_drops.png
#       = fig_ms_random_year_drops_state + fig_ms_sliding_window_state
#       inputs: data/processed/robustness_year_drops_random.csv, robustness_year_drops_sliding.csv (10_*.R)
#   B5  fig_app_compact_igo_share.png   [2026-11-12]
#       = IGO membership-share drops, log(1 + fatal MIDs) at t-1 only; rows system/state x cols lower/upper/both
#       inputs: data/processed/robustness_igo_share_results.csv (14_*.R)
#
# Output: code/03_figures/output/manuscript/  (canvases sized for \linewidth in the appendix)
# Run from the repository root:  Rscript code/03_figures/05_appendix_compact_figures.R
#
# Other presentation changes: long subtitles dropped (the LaTeX captions carry the sign convention and the
# gap-window explanation); B4 in-bar percentages dropped (shares printed to the log); B3 two rows, fatal MIDs only.
# Plot rules: horizontal x-axis text (long labels broken over two lines), facet strips on top
# (facet_grid row strips stay on the right), and no text drawn over data marks.
###############################################################################

suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(readr)
  library(patchwork); library(forcats); library(scales)
})
set.seed(123)   # jitter in B1 panel B (same seed as 91_figures.R)

ROOT    <- Sys.getenv("PIPE_ROOT", unset = here::here())   # same override as 91_figures.R
PROC    <- file.path(ROOT, "data", "processed")
FIG_DIR <- file.path(ROOT, "code", "03_figures", "output", "manuscript")
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

ACCENT  <- "#2C6CB0"
ACCENT2 <- "#B03366"

theme_app <- function(base_size = 12) {
  theme_bw(base_size = base_size) %+replace%
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "gray90", linewidth = 0.3),
      plot.title   = element_text(size = base_size + 1, face = "bold", hjust = 0,
                                  margin = margin(b = 6)),
      axis.text    = element_text(color = "black", size = base_size - 1),
      axis.text.x  = element_text(angle = 0, hjust = 0.5),
      axis.title   = element_text(color = "black", size = base_size),
      strip.background = element_rect(fill = "gray92", color = "gray70"),
      strip.text   = element_text(face = "bold", size = base_size - 1, color = "black",
                                  margin = margin(3, 3, 3, 3)),
      strip.placement = "outside",
      legend.position = "bottom",
      legend.background = element_blank(),
      legend.text  = element_text(size = base_size - 1),
      plot.margin  = margin(6, 10, 6, 6)
    )
}

save_fig <- function(filename, plot, width, height) {
  out <- file.path(FIG_DIR, filename)
  ggsave(out, plot, width = width, height = height, dpi = 300, bg = "white")
  stopifnot(file.exists(out), file.info(out)$size > 10000)
  cat(sprintf("  -> %s (%.1f x %.1f in)\n", filename, width, height))
}

read_chk <- function(f, cols) {
  p <- file.path(PROC, f)
  if (!file.exists(p)) stop("missing input: ", p)
  d <- read_csv(p, show_col_types = FALSE)
  miss <- setdiff(cols, names(d))
  if (length(miss)) stop(f, " lacks columns: ", paste(miss, collapse = ", "))
  cat(sprintf("  read %-40s %6d rows  (modified %s)\n", f, nrow(d), format(file.info(p)$mtime, "%Y-%m-%d %H:%M")))
  d
}

cat("\n=== Inputs ===\n")

###############################################################################
# B1. Replication sweep: democratic peace + capitalist peace (one figure)
###############################################################################
cat("\n=== B1: sweep families ===\n")
SW  <- file.path(ROOT, "code", "04_replication_sweep")
EST <- file.path(ROOT, "results", "replication_sweep")
reg <- read.csv(file.path(SW, "registry.csv"), colClasses = "character")
ids <- unique(reg$id)                                   # same id set as 91_figures.R with RUN_APPENDIX_FAMILIES = TRUE
f   <- file.path(EST, paste0("est_", ids, ".rds"))
if (!all(file.exists(f))) stop("missing estimates for: ", paste(ids[!file.exists(f)], collapse = ", "))
est <- bind_rows(lapply(f, function(x) { e <- readRDS(x); e$row_id <- as.character(e$row_id); e }))
stopifnot(!anyDuplicated(est$row_id))
stopifnot("pool_status" %in% names(reg),
          all((reg$pool_status == "main") == (reg$study_headline_term == "Y" & reg$in_headline == "Y")))
D <- reg %>% filter(study_headline_term == "Y", in_headline == "Y") %>%
  select(row_id, id, family, short_cite, named_illustration, citations_openalex, alpha_pub) %>%
  mutate(alpha_pub = as.numeric(alpha_pub)) %>%
  inner_join(est %>% select(-id, -key_term), by = "row_id") %>%
  mutate(citations = suppressWarnings(as.numeric(citations_openalex)))
stopifnot(!anyDuplicated(D$id))

ZC <- 8
zf <- function(b, p) pmin(abs(qnorm(pmax(p, 1e-300) / 2)), ZC)
TH <- theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold", size = 13),
        legend.position = "bottom", plot.caption = element_text(hjust = 0, size = 10),
        plot.caption.position = "plot")
FIG_NOTE <- "Shaded band: |z| < 1.96 (.05 reference). Colors use each study's own significance threshold."   # [2026-11-10]
STEPS <- c("Conflict\nhistory", "Benchmark\n(lagged CINC ratio)", "Benchmark\n(joint democracy)")

# expected counts from the appendix text (Section 3.1): studies, and lose + flip
EXPECT <- list(DEMOCRACY = c(n = 36, lose_flip = 7, med = 0.78), CAPITALISM = c(n = 10, lose_flip = 3, med = 0.83))   # [2026-11-13] final audit pool (126)
FAM_NAME <- c(DEMOCRACY = "Democratic peace", CAPITALISM = "Capitalist peace")

family_row <- function(FAM, letter_a, letter_b, show_caption) {
  d <- D %>% filter(family == FAM)
  stopifnot(nrow(d) > 0, all(is.finite(d$b0)), all(is.finite(d$b_hist)))
  d <- d %>% mutate(
    z0 = zf(b0, p0), zh = zf(b_hist, p_hist) * sign(b_hist) * sign(b0),
    status = case_when(p0 >= alpha_pub ~ "ns", p_hist >= alpha_pub ~ "lose",
                       sign(b_hist) != sign(b0) ~ "flip", TRUE ~ "stay"),
    ratio_h = abs(b_hist) / abs(b0) * sign(b_hist) * sign(b0),
    ratio_p = 1 + sign(b0) * placebo_pct / 100, ratio_d = 1 + sign(b0) * jointdem_pct / 100,
    cites = ifelse(is.na(citations), 0, citations))
  n <- table(factor(d$status, c("stay", "lose", "flip", "ns")))
  cat(FAM, "status counts:\n"); print(n)
  # fail loudly if the figure would disagree with the appendix text
  if (nrow(d) != EXPECT[[FAM]][["n"]] || (n[["lose"]] + n[["flip"]]) != EXPECT[[FAM]][["lose_flip"]])
    stop(FAM, ": counts differ from the appendix text (n = ", nrow(d), ", lose + flip = ",
         n[["lose"]] + n[["flip"]], "). Update the text or check the inputs.")
  LAB <- c(stay = sprintf("Stays significant (%d)", n[["stay"]]), lose = sprintf("Loses significance (%d)", n[["lose"]]),
           flip = sprintf("Significant, sign flips (%d)", n[["flip"]]), ns = sprintf("Not significant at baseline (%d)", n[["ns"]]))
  d$status <- factor(LAB[d$status], LAB)
  PAL <- setNames(c("#2a78d6", "#eb6834", "#8a2be2", "grey65"), LAB)

  pA <- ggplot(d, aes(z0, zh)) +
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = -1.96, ymax = 1.96, fill = "grey92") +
    annotate("text", x = 0.15, y = -2.6, label = "|z| < 1.96", hjust = 0, size = 3.6, colour = "grey45") +  # below the band, clear of points
    geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey40") +
    geom_hline(yintercept = 0, colour = "grey70") +
    geom_point(aes(size = cites, colour = status), alpha = .8) +
    scale_colour_manual(values = PAL, name = NULL, drop = FALSE,
                        guide = guide_legend(order = 1, nrow = 1, override.aes = list(size = 3))) +
    scale_size_area(max_size = 6, breaks = c(50, 150, 300), name = "Citations",
                    guide = guide_legend(order = 2, nrow = 1)) +
    scale_x_continuous(limits = c(0, ZC + 0.3), breaks = c(0, 2, 4, 6, 8), labels = c("0", "2", "4", "6", "8+")) +
    scale_y_continuous(limits = c(-3, ZC + 0.3), breaks = c(-2, 0, 2, 4, 6, 8), labels = c("-2", "0", "2", "4", "6", "8+")) +
    labs(title = sprintf("%s. %s: significance (%d studies)", letter_a, FAM_NAME[[FAM]], nrow(d)),
         x = "|z|, published specification", y = "z with conflict-history controls\n(sign relative to baseline)",
         caption = if (show_caption) FIG_NOTE else NULL) + TH

  sl <- bind_rows(
    d %>% transmute(step = STEPS[1], r = ratio_h, status),
    d %>% transmute(step = STEPS[2], r = ratio_p, status),
    d %>% transmute(step = STEPS[3], r = ratio_d, status)) %>% filter(!is.na(r))
  k <- table(sl$step)
  sl <- sl %>% mutate(step = factor(sprintf("%s\nn = %d", step, k[step]), rev(sprintf("%s\nn = %d", STEPS, k[STEPS]))),
                      r = pmin(pmax(r, -0.5), 2))
  med <- sl %>% group_by(step) %>% summarise(m = median(r), .groups = "drop")
  print(med)
  m_hist <- med$m[grepl("^Conflict", med$step)]
  # [2026-11-12] summary for the appendix table that replaces this figure (written after both families run)
  FAM_SUMMARY[[FAM]] <<- list(n = nrow(d), stay = n[["stay"]], lose = n[["lose"]], flip = n[["flip"]], ns = n[["ns"]],
                              med_hist = m_hist, med_cinc = med$m[grepl("CINC", med$step)],
                              med_jd = med$m[grepl("joint democracy", med$step)],
                              k_hist = k[[STEPS[1]]], k_cinc = k[[STEPS[2]]], k_jd = k[[STEPS[3]]])
  if (length(m_hist) != 1 || round(m_hist, 2) != EXPECT[[FAM]][["med"]])
    stop(FAM, ": conflict-history median ", round(m_hist, 3), " differs from the appendix text (", EXPECT[[FAM]][["med"]], ").")
  pB <- ggplot(sl, aes(r, step)) +
    geom_vline(xintercept = 1, linetype = 2, colour = "grey40") +
    geom_jitter(aes(colour = status), height = .17, width = 0, alpha = .7, size = 1.6) +
    geom_segment(data = med, aes(x = m, xend = m, y = as.numeric(step) - .3, yend = as.numeric(step) + .3), linewidth = 1.1) +
    geom_label(data = med, aes(x = m, y = as.numeric(step) + .42, label = sprintf("median %.2f", m)), size = 3.6,
               label.size = 0, fill = "white", label.padding = unit(1.5, "pt")) +
    scale_colour_manual(values = PAL, guide = "none") +
    scale_x_continuous(limits = c(-0.55, 2.05), breaks = c(-0.5, 0, .5, 1, 1.5, 2),
                       labels = c("<= -0.5", "0", "0.5", "1\n(no change)", "1.5", ">= 2")) +
    labs(title = sprintf("%s. %s: coefficient ratio", letter_b, FAM_NAME[[FAM]]),
         x = "Coefficient after / coefficient before", y = NULL) + TH

  # each family keeps its own legend because the status counts differ
  wrap_elements(full = (pA + pB + plot_layout(widths = c(1, 1.1), guides = "collect") &
                          theme(legend.position = "bottom", legend.box = "vertical",
                                legend.spacing.y = unit(1, "pt"), legend.margin = margin(0, 0, 0, 0))))
}

FAM_SUMMARY <- list()
fig_b1 <- family_row("DEMOCRACY", "A", "B", show_caption = FALSE) /
          family_row("CAPITALISM", "C", "D", show_caption = TRUE)
save_fig("fig_app_compact_sweep_families.png", fig_b1, width = 9.5, height = 10.5)

# [2026-11-12] Table version of B1 (replaces the figure in the appendix): results/appendix/tableA_sweep_families.tex
fs <- FAM_SUMMARY
stopifnot(setequal(names(fs), c("DEMOCRACY", "CAPITALISM")),
          all(sapply(fs, function(x) x$ns) == 0),                                  # every pooled study significant at baseline
          all(sapply(fs, function(x) x$stay + x$lose + x$flip) == sapply(fs, function(x) x$n)),
          all(sapply(fs, function(x) length(x$med_cinc) == 1 && length(x$med_jd) == 1)))
f2 <- function(x) sprintf("%.2f", x)
tab_row <- function(FAM) { x <- fs[[FAM]]
  sprintf("%s & %d & %d & %d & %d & %s & %s & %s \\\\", FAM_NAME[[FAM]], x$n, x$stay, x$lose, x$flip,
          f2(x$med_hist), f2(x$med_cinc), f2(x$med_jd)) }
# [2026-11-12] "all 11" rather than "11 of the 11"
bench_part <- function(k, n, lab) if (k == n) sprintf("all %d %s studies", n, lab) else sprintf("%d of the %d %s studies", k, n, lab)
bench_n <- paste0("Benchmarks are estimable for ",
                  bench_part(min(fs$DEMOCRACY$k_cinc, fs$DEMOCRACY$k_jd), fs$DEMOCRACY$n, "democratic peace"), " and ",
                  bench_part(min(fs$CAPITALISM$k_cinc, fs$CAPITALISM$k_jd), fs$CAPITALISM$n, "capitalist peace"), ".")
tab <- c(
  "% Generated by code/03_figures/05_appendix_compact_figures.R (B1 table) -- do not edit by hand.",
  "\\begin{table}[!htbp]",
  "\\centering",
  "\\footnotesize",
  "\\caption{Replication sweep results for the democratic peace and capitalist peace studies.}",
  "\\label{fig:contam_democracy}",
  "\\label{fig:contam_capitalism}",
  "\\setlength{\\tabcolsep}{3pt}",   # 4pt overflowed \\linewidth by 6pt
  "\\begin{tabular}{@{}lrrrrrrr@{}}",
  "\\toprule",
  " & & \\multicolumn{3}{c}{With conflict-history controls} & \\multicolumn{3}{c}{Median coefficient ratio} \\\\",
  "\\cmidrule(lr){3-5}\\cmidrule(lr){6-8}",
  "Literature & Studies & Stays sig. & Loses sig. & Sign flips & Conflict history & Lagged CINC ratio & Joint democracy \\\\",
  "\\midrule",
  tab_row("DEMOCRACY"), tab_row("CAPITALISM"),
  "\\bottomrule",
  "\\end{tabular}",
  "\\par\\smallskip",
  paste0("\\begin{minipage}{\\linewidth}\\scriptsize \\textit{Note:} Significance uses each study's own published threshold; every pooled study is significant at baseline. ",
         "The coefficient ratio is the key coefficient after adding the control divided by the coefficient before the control is added (1 = no change; below 1 = attenuation; negative = sign flip). ",
         "The benchmark columns add a lagged CINC ratio or joint democracy in place of conflict history. ", bench_n, "\\end{minipage}"),
  "\\end{table}")
tab_out <- file.path(ROOT, "results", "appendix", "tableA_sweep_families.tex")
dir.create(dirname(tab_out), recursive = TRUE, showWarnings = FALSE)
writeLines(tab, tab_out)
cat("  ->", tab_out, "\n"); cat(tab, sep = "\n")

###############################################################################
# B2. Contamination specification: samples/gapped window + no controls
###############################################################################
cat("\n=== B2: contamination specification ===\n")
spec <- read_chk("contamination_specification.csv",
                 c("sample_name", "layer", "process", "window", "attenuation_pct"))
nc   <- read_chk("contamination_nocontrols_ladder.csv", c("layer", "term", "attenuation_pct"))
stopifnot(setequal(unique(spec$window), c("main", "gap")),
          setequal(unique(spec$layer), c("Alliance", "IGO", "Multiplex")),
          setequal(unique(nc$layer),   c("Alliance", "IGO", "Multiplex")))
NC_TERMS <- c("Dyadic conflict history", "Monadic conflict history", "Between polity (regime gap)",
              "Monadic polity (own regime)", "Between CINC (asymmetry)", "Monadic CINC (own capability)")
stopifnot(setequal(unique(nc$term), NC_TERMS))
stopifnot(nrow(spec) == 24, nrow(nc) == 18, !anyNA(spec$attenuation_pct), !anyNA(nc$attenuation_pct))

p_spec <- spec |>
  filter(!is.na(attenuation_pct)) |>
  mutate(spec_label = paste0(sample_name, "\n",
                             ifelse(window == "main", "history t-2..t-6", "history t-3..t-7")),
         layer   = factor(layer, levels = c("Alliance", "IGO", "Multiplex")),
         process = factor(process, levels = c("Dyadic", "Monadic"),
                          labels = c("Dyadic conflict\nhistory (i-j MIDs)",
                                     "Monadic conflict\nhistory (third-party)"))) |>
  ggplot(aes(x = attenuation_pct, y = spec_label)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray55") +
  geom_col(fill = ACCENT, width = 0.62) +
  facet_grid(process ~ layer, scales = "free_x") +
  labs(title = NULL,
       x = "Attenuation (%)", y = NULL) +
  theme_app(11) +
  theme(panel.grid.major.y = element_blank(), axis.text.y = element_text(size = 9, lineheight = 0.9))

p_nc <- nc |>
  filter(!is.na(attenuation_pct)) |>
  mutate(kind  = if_else(grepl("conflict", term), "Conflict history", "Benchmark"),
         term  = factor(term, levels = rev(NC_TERMS)),
         layer = factor(layer, levels = c("Alliance", "IGO", "Multiplex"))) |>
  ggplot(aes(x = attenuation_pct, y = term, fill = kind)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray55") +
  geom_col(width = 0.62) +
  facet_wrap(~ layer, ncol = 3, scales = "free_x") +
  scale_fill_manual(name = NULL, values = c("Conflict history" = ACCENT, "Benchmark" = "#8C8C8C")) +
  scale_y_discrete(labels = label_wrap(18)) +
  labs(title = "B. No controls (dyad and year fixed effects only)", x = "Attenuation (%)", y = NULL) +
  theme_app(11) +
  theme(panel.grid.major.y = element_blank(), axis.text.y = element_text(size = 9))

# [2026-11-12] no-controls panel (p_nc) dropped from the figure per author; p_nc is still built and checked above.
fig_b2 <- p_spec
save_fig("fig_app_compact_contam_spec.png", fig_b2, width = 9, height = 5)

###############################################################################
# B3. Fatal-MID samples: system | state, log(1 + fatal MIDs) at t-1
###############################################################################
cat("\n=== B3: war samples ===\n")
# [2026-11-12] Two rows x one column, log(1 + fatal MIDs) at t-1 only (sits beside the sweep dial in the appendix).
# A = system (JS divergence): full sample and years with >= 3/5/10 fatal MIDs system-wide (11, samples (1), (3)).
# B = state (churn): fatal-onset IV with the fatality threshold raised from level 1 (main) to 2, 3, 4 (11, rows (1), (4)).
ws <- read_chk("robustness_war_samples.csv", c("level", "sample", "iv", "iv_label", "estimate", "ci_low", "ci_high"))
ws_pick <- function(lev, iv_col, samples) {
  out <- ws |>
    filter(level == lev, iv == iv_col, sample %in% samples) |>
    mutate(sample_lab = factor(names(samples)[match(sample, samples)], levels = names(samples)))
  stopifnot(nrow(out) == length(samples), !anyDuplicated(out$sample_lab),
            !anyNA(out$estimate), !anyNA(out$ci_low), !anyNA(out$ci_high))
  out
}
sys_df <- ws_pick("system", "log_fatal_mid_count_lag1",
                  c("Full sample"                 = "(1) Full sample",
                    "Years with >= 3 fatal MIDs"  = "(3) Fatal MIDs >= 3",
                    "Years with >= 5 fatal MIDs"  = "(3) Fatal MIDs >= 5",
                    "Years with >= 10 fatal MIDs" = "(3) Fatal MIDs >= 10"))
st_df  <- ws_pick("state", "log_fatal_onset_lag1",
                  c("Fatality >= 1\n(any deaths; main)" = "(1) Full sample",
                    "Fatality >= 2\n(26+ deaths)"       = "(4) Fatality >= 2",
                    "Fatality >= 3\n(101+ deaths)"      = "(4) Fatality >= 3",
                    "Fatality >= 4\n(251+ deaths)"      = "(4) Fatality >= 4"))

ws_panel <- function(dat, ttl) {
  ggplot(dat, aes(x = estimate, y = fct_rev(sample_lab))) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray55") +
    geom_errorbarh(aes(xmin = ci_low, xmax = ci_high), height = 0.15, linewidth = 0.6, color = ACCENT) +
    geom_point(size = 2.2, color = ACCENT) +
    scale_x_continuous(breaks = breaks_pretty(n = 3),
                       labels = function(x) format(x, scientific = FALSE, drop0trailing = TRUE, trim = TRUE)) +
    labs(title = ttl, x = "Coefficient (95% CI)", y = NULL) +
    theme_app(11) +
    theme(axis.text.y = element_text(size = 9, lineheight = 0.9))
}
# panel titles carry no letters: 06_appendix_grid_fig3_4.R places this figure as Panel B beside the sweep dial
fig_b3 <- ws_panel(sys_df, "System (JS divergence)") / ws_panel(st_df, "State (churn)")
save_fig("fig_app_compact_war_samples.png", fig_b3, width = 5, height = 6)
# [2026-11-12] ggplot objects for the appendix grid with the sweep dial (06_appendix_grid_fig3_4.R)
saveRDS(list(sys = ws_panel(sys_df, "System (JS divergence)"), st = ws_panel(st_df, "State (churn)")),
        file.path(FIG_DIR, "fig_app_compact_war_samples_plots.rds"))

###############################################################################
# B4. Temporal robustness (state level): random drops + sliding window
###############################################################################
cat("\n=== B4: year drops ===\n")
rnd <- read_chk("robustness_year_drops_random.csv", c("level", "iv_label", "estimate", "se", "drop_size", "rep"))
sld <- read_chk("robustness_year_drops_sliding.csv", c("level", "iv_label", "estimate", "ci_low", "window_size", "mid_year"))
DROP_SIZES <- c(5, 10, 15, 20, 25)

rnd_st <- rnd |>
  filter(level == "state", iv_label == "log Fatal MIDs") |>             # same restriction as 10_robustness_year_drops.R
  mutate(drop_n  = factor(drop_size, levels = DROP_SIZES),
         sig_pos = (abs(estimate / se) > qnorm(0.975)) & (estimate > 0),
         sig_cat = factor(if_else(sig_pos, "Positive & significant (p < 0.05)", "Null"),
                          levels = c("Positive & significant (p < 0.05)", "Null")))
reps <- rnd_st |> count(drop_n)
print(reps)
stopifnot(nrow(reps) == 5, all(reps$n == 500))

p_box <- ggplot(rnd_st, aes(x = drop_n, y = estimate)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray60") +
  geom_boxplot(fill = ACCENT, alpha = 0.7, width = 0.6, outlier.size = 0.4, linewidth = 0.35) +
  labs(x = "Years dropped (random)", y = "Coefficient", title = "A. Random drops: coefficient\n(state churn, log fatal MIDs t-1)") +
  theme_app(11)

share_df <- rnd_st |>
  count(drop_n, sig_cat, .drop = FALSE) |>
  group_by(drop_n) |> mutate(share = n / sum(n)) |> ungroup()
print(share_df |> filter(sig_cat == "Positive & significant (p < 0.05)"))
# The original figure printed the share inside each bar. That text sat on the bars, so it is dropped here;
# the shares are printed to the log above.
p_bar <- ggplot(share_df, aes(x = drop_n, y = share, fill = sig_cat)) +
  geom_col(width = 0.6) +
  scale_fill_manual(values = c("Positive & significant (p < 0.05)" = "#2E75B6", "Null" = "#C0392B"), name = NULL) +
  scale_y_continuous(labels = percent_format(accuracy = 1), expand = expansion(mult = c(0, 0.02))) +
  labs(x = "Years dropped (random)", y = "Share of 500 reps", title = "B. Random drops: share positive and\nsignificant (state churn, log fatal MIDs t-1)") +
  theme_app(11)

sld_st <- sld |>
  filter(level == "state", iv_label == "log Fatal MIDs") |>
  mutate(window_label = factor(paste0("Window = ", window_size, " yrs"),
                               levels = paste0("Window = ", DROP_SIZES, " yrs")))
stopifnot(nrow(sld_st) > 0, setequal(unique(sld_st$window_size), DROP_SIZES))

p_sld <- ggplot(sld_st, aes(x = mid_year, y = estimate, color = window_label, fill = window_label)) +
  geom_hline(yintercept = 0, color = "gray55", linewidth = 0.3) +
  geom_ribbon(aes(ymin = ci_low, ymax = estimate), alpha = 0.12, color = NA) +
  geom_line(linewidth = 0.6) +
  scale_color_brewer(palette = "Dark2", name = NULL) +
  scale_fill_brewer(palette = "Dark2", name = NULL) +
  labs(title = "C. Sliding-window exclusion (state churn, log fatal MIDs t-1)", x = "Midpoint year of excluded window",
       y = "Coefficient\n(95% CI, one-sided)") +
  theme_app(11)

fig_b4 <- (p_box | p_bar) / p_sld + plot_layout(heights = c(1, 1))
save_fig("fig_app_compact_year_drops.png", fig_b4, width = 9, height = 8.5)

###############################################################################
# B5. IGO membership-share drops: log(1 + fatal MIDs) at t-1 only   [2026-11-12]
#     rows = system (JS) / state (churn); columns = lower / upper / both tails dropped.
#     Replaces fig_ms_igo_share_both_system.png + fig_ms_igo_share_both_state.png in the appendix.
#     input: data/processed/robustness_igo_share_results.csv (14_igo_membership_sliding.R; plot only, no refit)
###############################################################################
cat("\n=== B5: IGO share drops ===\n")
igo <- read_chk("robustness_igo_share_results.csv",
                c("level", "iv", "L", "U", "grid", "estimate", "ci_low", "ci_high"))
igo_f <- igo |>
  filter((level == "system" & iv == "log_fatal_mid_count_lag1") |
         (level == "state"  & iv == "log_fatal_onset_lag1"))
cnt <- table(igo_f$level, igo_f$grid)
stopifnot(nrow(igo_f) == 48, setequal(colnames(cnt), c("both", "lower", "upper")), all(cnt[, "lower"] == 4), all(cnt[, "upper"] == 4), all(cnt[, "both"] == 16),
          !anyNA(igo_f$estimate))
fmt_pct <- function(x) sprintf("%.0f%%", x * 100)

igo_panel <- function(lev, grd, ttl = NULL, ylab = NULL) {
  dat <- igo_f |> filter(level == lev, grid == grd) |>
    mutate(L_label = factor(fmt_pct(L), levels = fmt_pct(sort(unique(L)))),
           U_label = factor(fmt_pct(U), levels = fmt_pct(sort(unique(U)))))
  if (grd == "both") {
    p <- ggplot(dat, aes(x = L_label, y = estimate, color = U_label, group = U_label)) +
      geom_hline(yintercept = 0, color = "gray55", linewidth = 0.3) +
      geom_errorbar(aes(ymin = ci_low, ymax = ci_high), width = 0.15, linewidth = 0.5,
                    position = position_dodge(width = 0.6)) +
      geom_point(size = 1.8, position = position_dodge(width = 0.6)) +
      scale_color_brewer(palette = "Dark2", name = "Upper threshold U") +
      labs(x = "Lower threshold L")
  } else {
    xv <- if (grd == "lower") "L_label" else "U_label"
    p <- ggplot(dat, aes(x = .data[[xv]], y = estimate)) +
      geom_hline(yintercept = 0, color = "gray55", linewidth = 0.3) +
      geom_errorbar(aes(ymin = ci_low, ymax = ci_high), width = 0.15, linewidth = 0.6, color = ACCENT) +
      geom_point(size = 2.2, color = ACCENT) +
      labs(x = if (grd == "lower") "Lower threshold L" else "Upper threshold U")
  }
  p + labs(title = ttl, y = ylab) + theme_app(11) +
    theme(axis.text.x = element_text(size = 9), plot.title = element_text(size = 11))
}
YS <- "System (JS)\nCoefficient (95% CI)"
YT <- "State (churn)\nCoefficient (95% CI)"
fig_b5 <- (igo_panel("system", "lower", "Lower tail dropped", YS) |
           igo_panel("system", "upper", "Upper tail dropped") |
           igo_panel("system", "both",  "Both tails dropped")) /
          (igo_panel("state",  "lower", NULL, YT) |
           igo_panel("state",  "upper") |
           igo_panel("state",  "both")) +
  plot_layout(guides = "collect") & theme(legend.position = "bottom")
save_fig("fig_app_compact_igo_share.png", fig_b5, width = 10, height = 6.5)

# [2026-11-12] Table version of B5 (replaces the figure in the appendix): results/appendix/tableA_igo_share.tex
# Rows: each lower-tail and upper-tail threshold, then the 16 both-tail combinations summarized by their range.
# [2026-11-13] stars match the main-text tables: * p < 0.05, ** p < 0.01, *** p < 0.001 (two-sided, normal approximation)
# old (2026-11-12): * p < 0.10, ** p < 0.05, *** p < 0.01
stars <- function(b, se) { p <- 2 * pnorm(-abs(b / se)); ifelse(p < 0.001, "$^{***}$", ifelse(p < 0.01, "$^{**}$", ifelse(p < 0.05, "$^{*}$", ""))) }
cell <- function(b, se) sprintf("%.4f%s (%.4f)", b, stars(b, se), se)
stopifnot("se" %in% names(igo), !anyNA(igo_f$se))
igo_t <- igo_f |> mutate(sig05 = abs(estimate / se) > 1.96)
one <- function(lev, grd, thr_col, thr) { r <- igo_t |> filter(level == lev, grid == grd, .data[[thr_col]] == thr)
  stopifnot(nrow(r) == 1); cell(r$estimate, r$se) }
rows <- c()
for (L in sort(unique(igo_t$L[igo_t$grid == "lower"])))
  rows <- c(rows, sprintf("%s & $L = %s$ & %s & %s \\\\", if (length(rows) == 0) "Lower tail" else "", sub("%", "\\\\%", fmt_pct(L)), one("system", "lower", "L", L), one("state", "lower", "L", L)))
rows <- c(rows, "\\addlinespace")
for (U in sort(unique(igo_t$U[igo_t$grid == "upper"])))
  rows <- c(rows, sprintf("%s & $U = %s$ & %s & %s \\\\", if (rows[length(rows)] == "\\addlinespace") "Upper tail" else "", sub("%", "\\\\%", fmt_pct(U)), one("system", "upper", "U", U), one("state", "upper", "U", U)))
# [2026-11-12] every both-tail (L, U) pair on its own row, with its standard error
rows <- c(rows, "\\addlinespace")
LU <- igo_t |> filter(grid == "both") |> distinct(L, U) |> arrange(L, U)
stopifnot(nrow(LU) == 16)
for (k in seq_len(nrow(LU))) {
  get_both <- function(lev) { r <- igo_t |> filter(level == lev, grid == "both", L == LU$L[k], U == LU$U[k])
    stopifnot(nrow(r) == 1); cell(r$estimate, r$se) }
  rows <- c(rows, sprintf("%s & $L = %s$, $U = %s$ & %s & %s \\\\", if (k == 1) "Both tails" else "",
                          sub("%", "\\\\%", fmt_pct(LU$L[k])), sub("%", "\\\\%", fmt_pct(LU$U[k])),
                          get_both("system"), get_both("state")))
}
n_sig <- sum(igo_t$sig05); n_all <- nrow(igo_t)
cat(sprintf("  IGO table: %d of %d estimates significant at 0.05\n", n_sig, n_all))
igo_tab <- c(
  "% Generated by code/03_figures/05_appendix_compact_figures.R (B5 table) -- do not edit by hand.",
  "\\begin{table}[!htbp]",
  "\\centering",
  "\\footnotesize",
  "\\caption{Conflict and network reconfiguration when IGOs are dropped by membership share.}",
  "\\label{fig:igo_sliding_both}",
  "\\begin{tabular}{@{}llcc@{}}",
  "\\toprule",
  "IGOs dropped & Threshold & System (JS divergence) & State (churn) \\\\",
  "\\midrule",
  rows,
  "\\bottomrule",
  "\\end{tabular}",
  "\\par\\smallskip",
  paste0("\\begin{minipage}{\\linewidth}\\scriptsize \\textit{Note:} Cells are coefficients (standard errors) on fatal MIDs$^{a}$ at $t-1$ from the main system and state models after dropping IGOs with membership share below $L$, above $U$, or both. ",
         sprintf("%s statistically significant at the 0.05 level. $^{*}p < 0.05$, $^{**}p < 0.01$, $^{***}p < 0.001$.", if (n_sig == n_all) sprintf("All %d estimates are", n_all) else sprintf("%d of the %d estimates are", n_sig, n_all)),
         "\\\\ $^{a}$ $\\log(1+x)$ of the fatal MID count (system) or fatal MID onset count (state).\\end{minipage}"),
  "\\end{table}")
igo_out <- file.path(ROOT, "results", "appendix", "tableA_igo_share.tex")
writeLines(igo_tab, igo_out)
cat("  ->", igo_out, "\n"); cat(igo_tab, sep = "\n")

###############################################################################
# B6. Main state specification under alternative measurement decisions (b), (e), (k), (f)   [2026-11-12]
#     Main = main-text Table 4 col. 3: TWFE (state + year FE, SEs clustered by state and year), churn on
#     log(1 + fatal MID onsets) at t-1, unweighted Manhattan on offensive + defensive alliances + IGO, all dyads,
#     no extra control. Each row changes ONE decision and keeps the rest at the main values. Churn is the
#     per-combination recomputation (churn_alt), which reproduces the main estimate exactly in the main cell.
#     Right column: wild cluster bootstrap p (fwildclusterboot, 999 Rademacher draws, null imposed, clustered by
#     state and year) from data/processed/robustness_specifications_wildboot.csv; run on all dyads only, so the
#     (k) rows have none.
#     input : data/processed/robustness_specifications_bellig.csv.gz (07_robustness_multiverse.R)
#     output: fig_app_main_spec_decisions.png
###############################################################################
cat("\n=== B6: main state specification under alternative decisions ===\n")
f_mv <- file.path(PROC, "robustness_specifications_bellig.csv.gz")
if (!file.exists(f_mv)) stop("missing input: ", f_mv)
mv <- read_csv(f_mv, show_col_types = FALSE, guess_max = 100000,
               col_select = c(level, model, war_variable, dv, transform, layers, weighting, bellig_rule,
                              extra_control, coef, se, n_obs))
MAIN <- list(model = "FE_full_ccode_year", war_variable = "log_fatal_mid_onset_count", dv = "churn_alt",
             transform = "raw", layers = "defensive+igo+offensive", weighting = "unweighted",
             bellig_rule = "none", extra_control = "none")
base <- mv |> filter(level == "state", model == MAIN$model, war_variable == MAIN$war_variable,
                     dv == MAIN$dv, transform == MAIN$transform)
rm(mv); invisible(gc())
pick <- function(...) { a <- list(...); m <- MAIN; m[names(a)] <- a
  base |> filter(layers == m$layers, weighting == m$weighting, bellig_rule == m$bellig_rule,
                 extra_control == m$extra_control) }
main_row <- pick()
stopifnot(nrow(main_row) == 1, abs(main_row$coef - 0.0191088) < 1e-5)   # reproduces main-text Table 4 col. 3

LAYER_LAB <- c(offensive = "Off", defensive = "Def", igo = "IGO", trade = "Trade")
lay_lab <- function(x) sapply(strsplit(x, "\\+"), function(p) paste(LAYER_LAB[c("offensive", "defensive", "igo", "trade")][
  c("offensive", "defensive", "igo", "trade") %in% p], collapse = " + "))
W_LAB <- c(contig_manhattan = "Contiguity-residualized Manhattan", dist_manhattan = "Distance-residualized Manhattan",
           contig_euclidean = "Contiguity-residualized Euclidean", dist_euclidean = "Distance-residualized Euclidean")
R_LAB <- c(onset_any_t1 = "MID onset at t-1", onset_any_t1t = "MID onset at t-1 or t",
           onset_fatal_t1 = "Fatal MID onset at t-1", onset_fatal_t1t = "Fatal MID onset at t-1 or t",
           ongoing_any_t1 = "Ongoing MID at t-1", ongoing_any_t1t = "Ongoing MID at t-1 or t",
           ongoing_fatal_t1 = "Ongoing fatal MID at t-1", ongoing_fatal_t1t = "Ongoing fatal MID at t-1 or t")

alt_layers <- setdiff(sort(unique(base$layers)), MAIN$layers)
rows <- bind_rows(
  main_row |> mutate(decision = "Main", label = "Main specification"),
  bind_rows(lapply(alt_layers, function(l) pick(layers = l) |> mutate(decision = "(b) Network layers", label = lay_lab(l)))),
  bind_rows(lapply(names(W_LAB), function(w) pick(weighting = w) |> mutate(decision = "(e) Edge weighting", label = W_LAB[[w]]))),
  bind_rows(lapply(names(R_LAB), function(r) pick(bellig_rule = r) |> mutate(decision = "(k) Dyads dropped", label = R_LAB[[r]]))),
  pick(extra_control = "gdp") |> mutate(decision = "(f) Extra control", label = "Adds log GDP"))
stopifnot(nrow(rows) == 1 + 11 + 4 + 8 + 1, !anyNA(rows$coef), !anyNA(rows$se))
wb <- read_csv(file.path(PROC, "robustness_specifications_wildboot.csv"), show_col_types = FALSE) |>
  filter(level == "state", model == "FE_full_wildboot", war_variable == MAIN$war_variable, dv == MAIN$dv,
         transform == MAIN$transform, bellig_rule == "none") |>
  select(layers, weighting, extra_control, wb_coef = coef, wb_p = pval)
stopifnot(nrow(wb) == 120, !anyDuplicated(wb[c("layers", "weighting", "extra_control")]))
rows <- rows |> left_join(wb, by = c("layers", "weighting", "extra_control")) |>
  mutate(wb_p = ifelse(bellig_rule == "none", wb_p, NA_real_))
# same coefficient as the clustered model wherever a bootstrap exists; main cell reproduces the reported p = 0.026
stopifnot(all(is.na(rows$wb_p) | abs(rows$wb_coef - rows$coef) < 1e-6),
          sum(!is.na(rows$wb_p)) == 1 + 11 + 4 + 1, abs(rows$wb_p[rows$decision == "Main"] - 0.026) < 0.001)
rows <- rows |>
  mutate(ci_lo = coef - 1.96 * se, ci_hi = coef + 1.96 * se, p = 2 * pnorm(-abs(coef / se)),
         sig = ifelse(p < 0.05, "p < 0.05", "p \u2265 0.05"),
         decision = factor(decision, levels = c("Main", "(b) Network layers", "(e) Edge weighting", "(k) Dyads dropped", "(f) Extra control"))) |>
  group_by(decision) |> arrange(coef, .by_group = TRUE) |> ungroup() |>
  mutate(label = factor(paste(decision, label, sep = "||"), levels = rev(unique(paste(decision, label, sep = "||")))))
cat(sprintf("  rows: %d; share p < 0.05: %.0f%%; coef range %.4f to %.4f\n", nrow(rows), 100 * mean(rows$p < 0.05),
            min(rows$coef), max(rows$coef)))
print(rows |> select(decision, label, coef, se, p, n_obs) |> mutate(label = sub(".*\\|\\|", "", label)), n = Inf)

p_coef <- ggplot(rows, aes(x = coef, y = label)) +
  annotate("rect", xmin = main_row$coef - 1.96 * main_row$se, xmax = main_row$coef + 1.96 * main_row$se,
           ymin = -Inf, ymax = Inf, fill = "grey92") +
  geom_vline(xintercept = 0, colour = "grey55") +
  geom_vline(xintercept = main_row$coef, linetype = "dashed", colour = ACCENT) +
  geom_errorbarh(aes(xmin = ci_lo, xmax = ci_hi, colour = sig), height = 0.25, linewidth = 0.6) +
  geom_point(aes(colour = sig), size = 2) +
  scale_colour_manual(values = c("p < 0.05" = ACCENT, "p \u2265 0.05" = "grey55"), name = NULL) +
  scale_y_discrete(labels = function(x) sub(".*\\|\\|", "", x)) +
  facet_grid(decision ~ ., scales = "free_y", space = "free_y") +
  labs(x = "Coefficient on log(1 + fatal MID onsets), t-1 (95% CI)", y = NULL) +
  theme_app(11) +
  theme(panel.grid.major.y = element_blank(), axis.text.y = element_text(size = 9),
        strip.text.y = element_blank(), legend.position = "bottom")
rows <- rows |> mutate(wb_lab = ifelse(is.na(wb_p), "\u2014", ifelse(wb_p < 0.001, "< 0.001", sprintf("%.3f", wb_p))))
p_wb <- ggplot(rows, aes(x = 0, y = label, label = wb_lab)) +
  geom_text(size = 3.2) +
  scale_x_continuous(limits = c(-1, 1), breaks = NULL) +
  facet_grid(decision ~ ., scales = "free_y", space = "free_y") +
  labs(x = "Wild bootstrap p", y = NULL) +
  theme_app(11) +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), panel.grid = element_blank(),
        strip.text.y = element_text(angle = 0, hjust = 0, size = 9))
fig_b6 <- p_coef + p_wb + plot_layout(widths = c(5, 1.6), guides = "collect") & theme(legend.position = "bottom")
save_fig("fig_app_main_spec_decisions.png", fig_b6, width = 9, height = 7.5)


print(sessionInfo())
cat("\nDone. Upload the five fig_app_compact_*.png files to Overleaf figs/.\n")
