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
#       = fig_ms_war_samples_system + fig_ms_hostlev_war_samples_state
#       inputs: data/processed/robustness_war_samples.csv (11_*.R), robustness_hostlev_results.csv (12_*.R)
#   B4  fig_app_compact_year_drops.png
#       = fig_ms_random_year_drops_state + fig_ms_sliding_window_state
#       inputs: data/processed/robustness_year_drops_random.csv, robustness_year_drops_sliding.csv (10_*.R)
#
# Output: code/03_figures/output/manuscript/  (canvases sized for \linewidth in the appendix)
# Run from the repository root:  Rscript code/03_figures/05_appendix_compact_figures.R
#
# Other presentation changes: long subtitles dropped (the LaTeX captions carry the sign convention and the
# gap-window explanation); B4 in-bar percentages dropped (shares printed to the log); B3 system facets 2 x 2.
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
EXPECT <- list(DEMOCRACY = c(n = 37, lose_flip = 5, med = 0.78), CAPITALISM = c(n = 11, lose_flip = 5, med = 0.75))
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

fig_b1 <- family_row("DEMOCRACY", "A", "B", show_caption = FALSE) /
          family_row("CAPITALISM", "C", "D", show_caption = TRUE)
save_fig("fig_app_compact_sweep_families.png", fig_b1, width = 9.5, height = 10.5)

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
  labs(title = "A. Alternative samples and gapped conflict-history window",
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

fig_b2 <- p_spec / p_nc + plot_layout(heights = c(1.15, 1))
save_fig("fig_app_compact_contam_spec.png", fig_b2, width = 9, height = 9.5)

###############################################################################
# B3. Fatal MIDs / battle-death subsets: system + state
###############################################################################
cat("\n=== B3: war samples ===\n")
ws <- read_chk("robustness_war_samples.csv", c("level", "sample", "iv_label", "estimate", "ci_low", "ci_high"))
hl <- read_chk("robustness_hostlev_results.csv", c("level", "iv_label", "estimate", "ci_low", "ci_high"))
stopifnot(all(hl$level == "state"))

stopifnot("(2) Drop WWI + WWII" %in% ws$sample)
sys_df <- ws |>
  filter(level == "system", sample != "(2) Drop WWI + WWII") |>      # same restriction as 11_robustness_war_samples.R
  mutate(sample   = factor(sample, levels = unique(sample)),
         iv_label = factor(iv_label, levels = unique(iv_label)))
stopifnot(nrow(sys_df) == 16)   # 4 samples x 4 IVs

p_sys <- ggplot(sys_df, aes(x = estimate, y = fct_rev(sample))) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray55") +
  geom_errorbarh(aes(xmin = ci_low, xmax = ci_high), height = 0.15, linewidth = 0.6, color = ACCENT) +
  geom_point(size = 2.2, color = ACCENT) +
  facet_wrap(~ iv_label, scales = "free_x", nrow = 2, labeller = label_wrap_gen(width = 32)) +
  scale_x_continuous(breaks = breaks_pretty(n = 3), labels = function(x) format(x, scientific = FALSE, drop0trailing = TRUE, trim = TRUE)) +
  labs(title = "A. System level (JS divergence): conflict-sample restrictions",
       x = "Coefficient (95% CI)", y = NULL) +
  theme_app(11) + theme(panel.spacing.x = unit(1.2, "lines"))

LEV_ORDER <- c("MID onset (t-1)", "log(1 + MID onset) (t-1)",
               "Fatal MID onset (t-1)", "log(1 + Fatal MID onset) (t-1)",
               "Mean HostLev (onsets, t-1)", "log(1 + Mean HostLev) (onsets, t-1)")
st_df <- bind_rows(
  ws |> filter(level == "state", sample == "(1) Full sample") |>
    transmute(iv_label, estimate, ci_low, ci_high, group = "MID sample (full)"),
  hl |> transmute(iv_label, estimate, ci_low, ci_high, group = "Hostility level"))
stopifnot(nrow(st_df) == 6, setequal(st_df$iv_label, LEV_ORDER))   # 4 war-sample IVs + 2 HostLev IVs
st_df <- st_df |> mutate(iv_label = factor(iv_label, levels = LEV_ORDER),
                         group = factor(group, levels = c("MID sample (full)", "Hostility level")))

p_st <- ggplot(st_df, aes(x = estimate, y = fct_rev(iv_label), color = group)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray55") +
  geom_errorbarh(aes(xmin = ci_low, xmax = ci_high), height = 0.18, linewidth = 0.7) +
  geom_point(size = 2.6) +
  scale_color_manual(values = c("MID sample (full)" = ACCENT, "Hostility level" = ACCENT2), name = NULL) +
  scale_y_discrete(labels = label_wrap(22)) +
  labs(title = "B. State level (churn): MID-sample and hostility-level variants",
       x = "Coefficient (95% CI)", y = NULL) +
  theme_app(11)

fig_b3 <- p_sys / p_st + plot_layout(heights = c(1.3, 1))
save_fig("fig_app_compact_war_samples.png", fig_b3, width = 9, height = 9)

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

print(sessionInfo())
cat("\nDone. Upload the four fig_app_compact_*.png files to Overleaf figs/.\n")
