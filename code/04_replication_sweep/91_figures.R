# 91_figures.R -- contamination figures, built in memory from est_<id>.rds + registry.csv (no study table written).
#   Main text: results/main/fig_contamination_ties.{pdf,png} (ties; labels the five pre-specified network papers,
#              the five most-cited network-measure studies, registry named_illustration == "Y"; four since 2026-11-08).
#   Appendix:  results/appendix/fig_contamination_{democracy,capitalism}.{pdf,png} (no study labels).
#   Appendix:  results/appendix/fig_contamination_all_attempts.{pdf,png} [2026-11-08]: pooled studies plus every attempt kept out
#              of the pool (failed rebuilds, not reproduced, null on the common sample, author nulls); no study labels.
# Panel A: |z| of the key coefficient in the published specification vs z with conflict-history controls.
# Panel B: coefficient ratio (after / before) for history, the lagged-CINC placebo and the joint-democracy benchmark.
suppressPackageStartupMessages({library(ggplot2); library(dplyr); library(patchwork)})
set.seed(123)   # jitter in panel B
PIPE <- Sys.getenv("PIPE_ROOT", unset = here::here())
SW <- file.path(PIPE, "code/04_replication_sweep"); EST <- file.path(PIPE, "results/replication_sweep")
reg <- read.csv(file.path(SW, "registry.csv"), colClasses = "character")
RUN_APPENDIX_FAMILIES <- as.logical(Sys.getenv("RUN_APPENDIX_FAMILIES", "TRUE"))
ids <- unique(reg$id[reg$family == "TIE" | RUN_APPENDIX_FAMILIES])
f <- file.path(EST, paste0("est_", ids, ".rds"))
if (!all(file.exists(f))) stop("missing estimates for: ", paste(ids[!file.exists(f)], collapse = ", "))
extra <- setdiff(sub("^est_(.*)\\.rds$", "\\1", list.files(EST, pattern = "^est_.*\\.rds$")), ids)
if (length(extra)) stop("est files for studies not in this run: ", paste(extra, collapse = ", "))
est <- bind_rows(lapply(f, function(x) { e <- readRDS(x); e$row_id <- as.character(e$row_id); e }))
stopifnot(!anyDuplicated(est$row_id))
stopifnot("pool_status" %in% names(reg),
          all((reg$pool_status == "main") == (reg$study_headline_term == "Y" & reg$in_headline == "Y")))
D <- reg %>% filter(id %in% ids, study_headline_term == "Y", in_headline == "Y") %>%
  select(row_id, id, family, short_cite, named_illustration, citations_openalex, alpha_pub) %>%
  mutate(alpha_pub = as.numeric(alpha_pub)) %>%
  inner_join(est %>% select(-id, -key_term), by = "row_id") %>%
  mutate(citations = suppressWarnings(as.numeric(citations_openalex)))
stopifnot(!anyDuplicated(D$id))
FAMS <- list(TIE = c("tie", "main", "fig_contamination_ties"))
if (RUN_APPENDIX_FAMILIES) FAMS <- c(FAMS, list(DEMOCRACY = c("democracy", "appendix", "fig_contamination_democracy"),
                                             CAPITALISM = c("capitalist-peace", "appendix", "fig_contamination_capitalism")))
for (FAM in names(FAMS)) {
NOUN <- FAMS[[FAM]][1]; OUTDIR <- file.path(PIPE, "results", FAMS[[FAM]][2]); STEM <- FAMS[[FAM]][3]
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)
d <- D %>% filter(family == FAM)
stopifnot(nrow(d) > 0, all(is.finite(d$b0)), all(is.finite(d$b_hist)))
ZC <- 8
zf <- function(b, p) pmin(abs(qnorm(pmax(p, 1e-300) / 2)), ZC)
d <- d %>% mutate(
  z0 = zf(b0, p0), zh = zf(b_hist, p_hist) * sign(b_hist) * sign(b0),
  status = case_when(p0 >= alpha_pub ~ "ns", p_hist >= alpha_pub ~ "lose",   # paper-specific threshold (registry alpha_pub)
                     sign(b_hist) != sign(b0) ~ "flip", TRUE ~ "stay"),
  ratio_h = abs(b_hist) / abs(b0) * sign(b_hist) * sign(b0),
  # placebo_pct / jointdem_pct = 100 * (b_bench - b0) / |b0|, so b_bench / b0 = 1 + sign(b0) * pct / 100 (Codex audit 2026-10-01)
  ratio_p = 1 + sign(b0) * placebo_pct / 100, ratio_d = 1 + sign(b0) * jointdem_pct / 100,
  cites = ifelse(is.na(citations), 0, citations))
n <- table(factor(d$status, c("stay", "lose", "flip", "ns"))); cat(FAM, "status counts:\n"); print(n)
LAB <- c(stay = sprintf("Stays significant (%d)", n[["stay"]]), lose = sprintf("Loses significance (%d)", n[["lose"]]),
         flip = sprintf("Significant, sign flips (%d)", n[["flip"]]), ns = sprintf("Not significant at baseline (%d)", n[["ns"]]))
d$status <- factor(LAB[d$status], LAB)
PAL <- setNames(c("#2a78d6", "#eb6834", "#8a2be2", "grey65"), LAB)
TH <- theme_minimal(base_size = 14) + theme(panel.grid.minor = element_blank(),   # [2026-10-27] larger text for the main-text figure
  plot.title = element_text(face = "bold", size = 16), legend.position = "bottom",
  plot.caption = element_text(hjust = 0, size = 11), plot.caption.position = "plot")

# named illustrations: top-5 network papers by citations (rule fixed in advance); label positions set by hand
top5 <- data.frame(short_cite = c("Hafner-Burton & Montgomery 2006", "Maoz et al. 2007", "Maoz 2009", "Kinne 2012", "Maoz et al. 2006"),
                   hj = c(0, 1, 0, 0, 0), lx = c(0.2, 7.6, 5.9, 0.2, 5.9), ly = c(6.4, 8, 1.35, 5.4, -0.9))
lab <- d %>% inner_join(top5, by = "short_cite")
# [2026-11-08] Maoz 2009 left the pool (not reproduced), so four of the five pre-specified papers are labelled; the rule is not refilled
if (FAM == "TIE") stopifnot(nrow(lab) == 4, all(lab$named_illustration == "Y"), setequal(lab$id, D$id[D$named_illustration == "Y"])) else lab <- lab[0, ]

# Figure note (Codex fix 2026-10-11; analysis unchanged): the shaded band marks |z| < 1.96, a fixed .05 visual reference only.
# Status colours are classified with each paper's own threshold (registry alpha_pub, .05 or .10), so a point can sit inside the
# band and still be coloured significant (alpha_pub = .10), or outside it without changing status.
FIG_NOTE <- "Shaded band: |z| < 1.96 (.05 reference). Colours use each paper's own threshold (alpha_pub)."
pA <- ggplot(d, aes(z0, zh)) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = -1.96, ymax = 1.96, fill = "grey92") +
  annotate("text", x = 0.15, y = -1.75, label = "|z| < 1.96", hjust = 0, size = 4.2, colour = "grey45") +
  geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey40") +
  geom_hline(yintercept = 0, colour = "grey70") +
  geom_point(aes(size = cites, colour = status), alpha = .8) +
  geom_segment(data = lab %>% filter(short_cite != "Maoz et al. 2007"), aes(x = lx + ifelse(hj == 0, 0, 0), y = ly - 0.2, xend = z0, yend = zh),
               colour = "grey35", linewidth = .3) +
  geom_label(data = lab, aes(x = lx, y = ly, label = short_cite, hjust = hj), size = 4.2, label.size = 0, fill = "white",
             label.padding = unit(1.5, "pt")) +
  scale_colour_manual(values = PAL, name = NULL, drop = FALSE, guide = guide_legend(order = 1, nrow = 1, override.aes = list(size = 3))) +
  scale_size_area(max_size = 7, breaks = c(50, 150, 300), name = "Citations", guide = guide_legend(order = 2, nrow = 1)) +
  scale_x_continuous(limits = c(0, ZC + 0.3), breaks = c(0, 2, 4, 6, 8), labels = c("0", "2", "4", "6", "8+")) +
  scale_y_continuous(limits = c(-3, ZC + 0.3), breaks = c(-2, 0, 2, 4, 6, 8), labels = c("-2", "0", "2", "4", "6", "8+")) +
  labs(title = sprintf("A. Significance of the key %s coefficient\n(%d studies)", NOUN, nrow(d)),
       x = "|z|, published specification", y = "z with conflict-history controls\n(sign relative to baseline)", caption = FIG_NOTE) + TH

sl <- bind_rows(
  d %>% transmute(step = "Conflict\nhistory", r = ratio_h, status),
  d %>% transmute(step = "Benchmark\n(lagged CINC ratio)", r = ratio_p, status),
  d %>% transmute(step = "Benchmark\n(joint democracy)", r = ratio_d, status)) %>%
  filter(!is.na(r))
# studies without an estimable benchmark are reported, not silently dropped: the axis label carries the plotted n
k <- table(sl$step); cat(FAM, "panel B studies plotted (of", nrow(d), "):\n"); print(k)
STEPS <- c("Conflict\nhistory", "Benchmark\n(lagged CINC ratio)", "Benchmark\n(joint democracy)")
sl <- sl %>% mutate(step = factor(sprintf("%s\nn = %d", step, k[step]), rev(sprintf("%s\nn = %d", STEPS, k[STEPS]))),
                    r = pmin(pmax(r, -0.5), 2))
med <- sl %>% group_by(step) %>% summarise(m = median(r), k = n(), .groups = "drop")
print(med)
pB <- ggplot(sl, aes(r, step)) +
  geom_vline(xintercept = 1, linetype = 2, colour = "grey40") +
  geom_jitter(aes(colour = status), height = .17, width = 0, alpha = .7, size = 1.8) +
  geom_segment(data = med, aes(x = m, xend = m, y = as.numeric(step) - .3, yend = as.numeric(step) + .3), linewidth = 1.1) +
  geom_label(data = med, aes(x = m, y = as.numeric(step) + .38, label = sprintf("median %.2f", m)), size = 4.4,
             label.size = 0, fill = "white", label.padding = unit(1.5, "pt")) +   # white box so the label does not sit on the x = 1 line
  scale_colour_manual(values = PAL, guide = "none") +
  scale_x_continuous(limits = c(-0.55, 2.05), breaks = c(-0.5, 0, .5, 1, 1.5, 2),
                     labels = c("<= -0.5", "0", "0.5", "1\n(no change)", "1.5", ">= 2")) +
  labs(title = "B. Change in the key coefficient\nwhen each control is added",
       x = "Coefficient after / coefficient before", y = NULL) + TH

fig <- pA + pB + plot_layout(widths = c(1, 1.1), guides = "collect") &
  theme(legend.position = "bottom", legend.box = "vertical", legend.spacing.y = unit(2, "pt"),
        legend.margin = margin(0, 0, 0, 0))

ggsave(file.path(OUTDIR, paste0(STEM, ".png")), fig, width = 12, height = 6.8, dpi = 200, bg = "white")
# cairo_pdf needs XQuartz on macOS; fall back to the base pdf device where cairo is not available [2026-10-05]
# capabilities("cairo") can report TRUE on macOS even when cairo.so cannot load (no XQuartz), so test-open a device [2026-11-02]
PDF_DEV <- tryCatch({ tf <- tempfile(fileext = ".pdf"); cairo_pdf(tf); grDevices::dev.off(); unlink(tf); cairo_pdf },
                    warning = function(w) grDevices::pdf, error = function(e) grDevices::pdf)
ggsave(file.path(OUTDIR, paste0(STEM, ".pdf")), fig, width = 12, height = 6.8, device = PDF_DEV)
}

# ---- [2026-11-08] appendix: all attempts, including those kept out of the pool ----
# Same panel A design as above, one facet per family. Colour = how the study entered (or did not enter) the pool; shape = what
# happens to the key term with conflict history at the paper's own threshold. Out of scope, duplicates, methods-paper
# reanalyses and non-comparable estimands are not plotted (they do not test the article's headline claim).
GRP <- c(main = "Pooled", rebuild_outside_rule = "Rebuilt, fails the round-2 match rule", not_reproduced = "Not reproduced",
         null_common_sample = "Null on the common sample\nbefore adding history", author_null = "Author argues a null")
GPAL <- setNames(c("#2a78d6", "#eb6834", "#8a2be2", "grey55", "#1a9e77"), GRP)
FLAB <- c(TIE = "Ties", DEMOCRACY = "Democratic peace", CAPITALISM = "Capitalist peace")
DA <- reg %>% filter(id %in% ids, study_headline_term == "Y", pool_status %in% names(GRP)) %>%
  select(row_id, id, family, pool_status, alpha_pub) %>% mutate(alpha_pub = as.numeric(alpha_pub)) %>%
  inner_join(est %>% select(row_id, b0, p0, b_hist, p_hist), by = "row_id")
stopifnot(!anyDuplicated(DA$id), all(is.finite(DA$b0)), all(is.finite(DA$b_hist)), !anyNA(DA$p0), !anyNA(DA$p_hist))
ZC <- 8
zf <- function(b, p) pmin(abs(qnorm(pmax(p, 1e-300) / 2)), ZC)
SH <- c(stay = "Significant before and after", lose = "Loses significance", flip = "Significant, sign flips",
        gain = "Becomes significant", none = "Not significant before or after")
DA <- DA %>% mutate(
  z0 = zf(b0, p0), zh = zf(b_hist, p_hist) * sign(b_hist) * sign(b0),
  sig0 = p0 < alpha_pub, sigh = p_hist < alpha_pub,
  out = case_when(sig0 & !sigh ~ "lose", sig0 & sign(b_hist) != sign(b0) ~ "flip", sig0 ~ "stay", sigh ~ "gain", TRUE ~ "none"),
  grp = factor(GRP[pool_status], GRP), out = factor(SH[out], SH),
  fam = factor(FLAB[family], FLAB[FLAB %in% FLAB[unique(family)]]),
  zh = pmax(zh, -3)) %>%                      # bottom of the axis reads "-3 or below" (same lower limit as panel A above)
  arrange(desc(pool_status == "main"))   # excluded attempts drawn on top
ca <- DA %>% count(fam, grp, .drop = FALSE); cat("ALL ATTEMPTS plotted:\n"); print(as.data.frame(ca))
# caption counts: attempts whose author claims an effect (author nulls left out), significant at baseline, changed with history
cc <- DA %>% filter(pool_status != "author_null") %>% group_by(fam) %>%
  summarise(k = n(), s = sum(sig0), ch = sum(sig0 & (!sigh | sign(b_hist) != sign(b0))), .groups = "drop")
pc <- DA %>% filter(family == "TIE", pool_status == "main")
CAP <- paste0("All attempts whose author claims an effect, significant at baseline, that lose significance or flip sign: ",
              paste(sprintf("%s %d of %d", cc$fam, cc$ch, cc$s), collapse = "; "),
              sprintf(" (pooled ties: %d of %d).", sum(pc$sig0 & (!pc$sigh | sign(pc$b_hist) != sign(pc$b0))), sum(pc$sig0)),
              "\nShaded band: |z| < 1.96 (.05 reference). Colours and shapes use each paper's own threshold (alpha_pub).")
pAll <- ggplot(DA, aes(z0, zh)) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = -1.96, ymax = 1.96, fill = "grey92") +
  annotate("rect", xmin = 0, xmax = 1.96, ymin = -Inf, ymax = Inf, fill = "grey92", alpha = .6) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey40") +
  geom_hline(yintercept = 0, colour = "grey70") +
  geom_point(aes(colour = grp, shape = out), size = 2.6, alpha = .85, stroke = .9) +
  facet_wrap(~ fam, nrow = 1) +
  scale_colour_manual(values = GPAL, name = NULL, drop = FALSE, guide = guide_legend(order = 1, nrow = 2, override.aes = list(size = 3, shape = 16))) +
  scale_shape_manual(values = setNames(c(16, 17, 4, 1, 0), SH), name = NULL, drop = FALSE,
                     guide = guide_legend(order = 2, nrow = 2, override.aes = list(size = 3))) +
  scale_x_continuous(limits = c(0, ZC + 0.3), breaks = c(0, 2, 4, 6, 8), labels = c("0", "2", "4", "6", "8+")) +
  scale_y_continuous(limits = c(-3, ZC + 0.3), breaks = c(-3, -2, 0, 2, 4, 6, 8), labels = c("<= -3", "-2", "0", "2", "4", "6", "8+")) +
  labs(title = sprintf("Significance of the key coefficient, all replication attempts (%d studies)", nrow(DA)),
       x = "|z|, published specification (common sample)", y = "z with conflict-history controls\n(sign relative to baseline)",
       caption = CAP) +
  theme_minimal(base_size = 13) +
  theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold", size = 15), legend.position = "bottom",
        legend.box = "vertical", plot.caption = element_text(hjust = 0, size = 10), plot.caption.position = "plot",
        strip.text = element_text(face = "bold", size = 13))
ADIR <- file.path(PIPE, "results", "appendix"); dir.create(ADIR, recursive = TRUE, showWarnings = FALSE)
ggsave(file.path(ADIR, "fig_contamination_all_attempts.png"), pAll, width = 13, height = 6.6, dpi = 200, bg = "white")
ggsave(file.path(ADIR, "fig_contamination_all_attempts.pdf"), pAll, width = 13, height = 6.6, device = PDF_DEV)
