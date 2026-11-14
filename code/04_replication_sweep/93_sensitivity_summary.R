# 93_sensitivity_summary.R -- pooled appendix output from results/replication_sweep/sens/sens_<id>.rds
# (written by run_sensitivity.R; definitions in legacy/hist/sens_helper.R). No study names are read or written here.
# [2026-11-14] Headline set = EXPANDED set: studies in sens_studies.csv whose headline row has pool_status "main" (pooled,
# in_headline == "Y") or "outcome_overlap" (outcome codes every year of a dispute). The main text reports the expanded set;
# the appendix subset table (tableA_sens_subsets.tex) reports the pooled-only and overlap-only subsets beside it.
#
# Writes
#   results/replication_sweep/sens_summary.rds          aggregates only
#   results/appendix/tableA_sens_split.tex              dyadic only / monadic only / both: coefficient and SE ratios,
#                                                       significance lost, and the matched placebo for each
#   results/appendix/tableA_sens_strength.tex           partial R2 of history and of the two benchmarks
#   results/appendix/fig_sweep_split.{pdf,png}           coefficient and standard-error ratios, three specifications
#   results/appendix/fig_sweep_strength.{pdf,png}        confounder strength: history next to the benchmarks
#   results/appendix/fig_sweep_placebo.{pdf,png}         real history against its matched placebo, three specifications
#   results/appendix/fig_sweep_dial.{pdf,png}            lambda dial (studies with dial == "Y")
# Ratios are after / before, signed relative to the baseline coefficient (1 = no change, < 0 = sign flip).
suppressPackageStartupMessages({library(ggplot2); library(dplyr); library(patchwork)})
set.seed(123)   # jitter
PIPE <- Sys.getenv("PIPE_ROOT", unset = here::here())
SW <- file.path(PIPE, "code/04_replication_sweep"); EST <- file.path(PIPE, "results/replication_sweep")
SENS <- file.path(EST, "sens"); APP <- file.path(PIPE, "results/appendix"); dir.create(APP, recursive = TRUE, showWarnings = FALSE)
# capabilities("cairo") can report TRUE on macOS even when cairo.so cannot load (no XQuartz), so test-open a device [2026-11-02]
PDF_DEV <- tryCatch({ tf <- tempfile(fileext = ".pdf"); cairo_pdf(tf); grDevices::dev.off(); unlink(tf); cairo_pdf },
                    warning = function(w) grDevices::pdf, error = function(e) grDevices::pdf)
reg <- read.csv(file.path(SW, "registry.csv"), colClasses = "character")
ss  <- read.csv(file.path(SW, "sens_studies.csv"), colClasses = "character")
f <- file.path(SENS, paste0("sens_", ss$id, ".rds"))
if (!all(file.exists(f))) stop("missing sensitivity results for: ", paste(ss$id[!file.exists(f)], collapse = ", "))
extra <- setdiff(sub("^sens_(.*)\\.rds$", "\\1", list.files(SENS, pattern = "^sens_.*\\.rds$")), ss$id)
if (length(extra)) stop("sens files for studies not in sens_studies.csv: ", paste(extra, collapse = ", "))
O <- lapply(f, readRDS); names(O) <- ss$id
# scope = full: single additive tie column, all analyses. scope = split: headline term is an interaction or product term,
# three-way split only (no placebo, strength or dial; see sens_helper.R).
stopifnot(all(ss$scope %in% c("full", "split")), identical(unname(sapply(O, function(o) isTRUE(o$split_only))), ss$scope == "split"))
FULL <- ss$id[ss$scope == "full"]; scope <- setNames(ss$scope, ss$id)
stopifnot(all(sapply(O[FULL], function(o) isTRUE(o$check_ok))), length(unique(sapply(O[FULL], `[[`, "draws"))) == 1)
hd <- reg[reg$study_headline_term == "Y" & reg$id %in% ss$id, c("id", "family", "in_headline", "alpha_pub", "pool_status")]
alp <- setNames(as.numeric(hd$alpha_pub), hd$id); stopifnot(!anyNA(alp))
stopifnot(!anyDuplicated(hd$id), setequal(hd$id, ss$id))
EXP_STATUS <- c("main", "outcome_overlap")   # [2026-11-14] expanded set
keep <- hd$id[hd$pool_status %in% EXP_STATUS]; fam <- setNames(hd$family, hd$id); keepF <- intersect(keep, FULL)
pstat <- setNames(hd$pool_status, hd$id)
stopifnot(all((hd$pool_status == "main") == (hd$in_headline == "Y")), sum(pstat[keep] == "outcome_overlap") > 0)
cat("sensitivity summary:", length(keep), "expanded-set studies (", sum(pstat[keep] == "main"), "pooled +", sum(pstat[keep] == "outcome_overlap"), "outcome overlap;", length(ss$id) - length(keep), "others not summarised );",
    length(keepF), "with placebo and strength analyses,", length(keep) - length(keepF), "split only;",
    O[[FULL[1]]]$draws, "placebo draws per specification\n")
# [2026-11-14] each expanded-set study must reproduce its main-sweep estimate with both history controls (est_<id>.rds b_hist),
# so the sensitivity tables describe the same coefficients as the headline counts. Stops on a gap above 1%.
rid <- setNames(reg$row_id[reg$study_headline_term == "Y"], reg$id[reg$study_headline_term == "Y"])
repro <- data.frame(id = keep, b_sens = sapply(keep, function(i) O[[i]]$split["both", "b"]),
                    b_est = sapply(keep, function(i) { e <- readRDS(file.path(EST, paste0("est_", i, ".rds"))); e$b_hist[as.character(e$row_id) == rid[[i]]] }))
repro$rel <- abs(repro$b_sens - repro$b_est) / abs(repro$b_est)
cat("reproduction check (sens 'both' vs est b_hist): max relative gap", signif(max(repro$rel), 3), "\n")
if (any(repro$rel > 1e-6)) print(repro[repro$rel > 1e-6, ], row.names = FALSE)
if (any(repro$rel > 0.01)) stop("sensitivity results do not reproduce the main sweep for: ", paste(repro$id[repro$rel > 0.01], collapse = ", "),
                                " -- rerun: ONLY_IDS=<id> Rscript code/04_replication_sweep/run_sensitivity.R")
SPEC <- c(dyad = "Dyadic history\nonly", mono = "Monadic history\nonly", both = "Dyadic and\nmonadic history")

# ---- study x specification (in memory only) ----
D <- bind_rows(lapply(keep, function(id) { o <- O[[id]]; s <- o$split
  pl <- if (scope[[id]] == "full") o$placebo %>% filter(lambda == 0) %>% group_by(spec) %>% summarise(b_plac = mean(b), .groups = "drop") else
    data.frame(spec = character(0), b_plac = numeric(0))
  data.frame(id = id, family = fam[[id]], set = pstat[[id]], alpha_pub = alp[[id]], scope = scope[[id]], spec = names(SPEC), b0 = s["none", "b"], se0 = s["none", "se"], p0 = s["none", "p"],
             b = s[names(SPEC), "b"], se = s[names(SPEC), "se"], p = s[names(SPEC), "p"],
             b_plac = pl$b_plac[match(names(SPEC), pl$spec)]) })) %>%
  mutate(r = b / b0, r_se = se / se0, r_plac = b_plac / b0, sig0 = p0 < alpha_pub, lose = sig0 & p >= alpha_pub,   # paper-specific threshold (registry alpha_pub)
         flip = sig0 & p < alpha_pub & sign(b) != sign(b0))
stopifnot(all(is.finite(D$r)), all(is.finite(D$r_se)), all(is.finite(D$r_plac[D$scope == "full"])), all(is.na(D$r_plac[D$scope == "split"])), !anyNA(D$p))

tab_split <- function(x) x %>% group_by(spec) %>% summarise(
  studies = n(), med_coef_ratio = median(r), share_shrink = mean(abs(r) < 1), med_se_ratio = median(r_se),
  share_se_up = mean(r_se > 1), share_se_up5 = mean(r_se > 1.05), share_se_down5 = mean(r_se < 0.95),
  sig_baseline = sum(sig0), lose_sig = sum(lose), sign_flips = sum(flip),
  placebo_n = sum(!is.na(r_plac)), med_placebo_ratio = median(r_plac, na.rm = TRUE),
  share_placebo_within5 = mean(abs(r_plac - 1) < .05, na.rm = TRUE), .groups = "drop") %>%
  slice(match(names(SPEC), spec))
T_all <- tab_split(D); T_tie <- tab_split(D %>% filter(family == "TIE"))
T_add <- tab_split(D %>% filter(scope == "full")); T_int <- tab_split(D %>% filter(scope == "split"))
cat("\nAll pooled studies\n"); print(as.data.frame(T_all), digits = 3, row.names = FALSE)
cat("\nAdditive tie term\n"); print(as.data.frame(T_add), digits = 3, row.names = FALSE)
cat("\nInteraction or product term\n"); print(as.data.frame(T_int), digits = 3, row.names = FALSE)
cat("\nTie studies\n"); print(as.data.frame(T_tie), digits = 3, row.names = FALSE)

# ---- strength: partial R2 with the tie variable and with the outcome ----
SV <- c(history_both = "Conflict history\n(both)", history_dyad = "Conflict history\n(dyadic)", history_mono = "Conflict history\n(monadic)",
        joint_democracy = "Joint democracy\n(benchmark)", capability_ratio = "Capability ratio\n(benchmark)")
G <- bind_rows(lapply(keepF, function(id) { s <- O[[id]]$strength
  data.frame(id = id, family = fam[[id]], set = pstat[[id]], var = rownames(s), r2_tie = s[, "r2_tie"], r2_out = s[, "r2_out"]) }))
HISTV <- c("history_both", "history_dyad", "history_mono")
stopifnot(all(G$var %in% names(SV)), all(is.finite(G$r2_tie[G$var %in% HISTV])), all(is.finite(G$r2_out[G$var %in% HISTV])))
G <- G %>% filter(is.finite(r2_tie), is.finite(r2_out))   # a benchmark that is constant in a study's sample has no strength there
tab_str <- function(x) x %>% group_by(var) %>% summarise(studies = n(), med_r2_tie = median(r2_tie), med_r2_out = median(r2_out),
  .groups = "drop") %>% slice(match(names(SV), var))
# within-study comparison: is history a stronger predictor than each benchmark?
W <- G %>% filter(var %in% c("history_both", "joint_democracy", "capability_ratio")) %>%
  tidyr::pivot_wider(id_cols = c(id, family), names_from = var, values_from = c(r2_tie, r2_out))
gt <- function(a, b) mean(a > b, na.rm = TRUE)
cmp <- function(x) data.frame(studies = nrow(x), democ_n = sum(!is.na(x$r2_out_joint_democracy)), cap_n = sum(!is.na(x$r2_out_capability_ratio)),
  out_gt_democ = gt(x$r2_out_history_both, x$r2_out_joint_democracy), out_gt_cap = gt(x$r2_out_history_both, x$r2_out_capability_ratio),
  tie_gt_democ = gt(x$r2_tie_history_both, x$r2_tie_joint_democracy), tie_gt_cap = gt(x$r2_tie_history_both, x$r2_tie_capability_ratio),
  med_out_ratio_democ = median(x$r2_out_history_both / x$r2_out_joint_democracy, na.rm = TRUE),
  med_out_ratio_cap = median(x$r2_out_history_both / x$r2_out_capability_ratio, na.rm = TRUE))
S_all <- tab_str(G); S_tie <- tab_str(G %>% filter(family == "TIE")); C_all <- cmp(W); C_tie <- cmp(W %>% filter(family == "TIE"))
cat("\nStrength, all pooled studies\n"); print(as.data.frame(S_all), digits = 3, row.names = FALSE); print(C_all, digits = 3, row.names = FALSE)

# ---- dial ----
dial_ids <- intersect(keep, ss$id[ss$dial == "Y"])
L <- bind_rows(lapply(dial_ids, function(id) { o <- O[[id]]; stopifnot(isTRUE(o$dial))
  o$placebo %>% filter(spec == "both") %>% group_by(lambda) %>% summarise(b = mean(b), .groups = "drop") %>%
    bind_rows(data.frame(lambda = 1, b = o$split["both", "b"])) %>% mutate(id = id, r = b / o$split["none", "b"]) }))
T_dial <- if (nrow(L)) L %>% group_by(lambda) %>% summarise(studies = n(), med_ratio = median(r), .groups = "drop") else NULL
if (nrow(L)) { stopifnot(all(table(L$id) == 5)); cat("\nDial (", length(dial_ids), "studies)\n"); print(as.data.frame(T_dial), digits = 3, row.names = FALSE) }

# ---- [2026-11-14] appendix subset table: pooled only vs outcome overlap only vs expanded (both history controls) ----
# Same quantities as the main tables, computed separately for each subset. Partial R2 rows use studies with an additive key term.
sub_row <- function(dd, gg, lab) {
  t <- tab_split(dd) %>% filter(spec == "both")
  h <- gg %>% filter(var == "history_both")
  # [2026-11-14] readable version (no \\resizebox): share shrinking and the placebo ratio (always 1.00) dropped; both are in Table sens_split
  sprintf("%s & %d & %s & %s & %d of %d & %s & %s & %s \\\\", lab, t$studies, f2s(t$med_coef_ratio),
          f2s(t$med_se_ratio), t$lose_sig + t$sign_flips, t$sig_baseline,
          ifelse(t$placebo_n > 0, pcs(t$share_placebo_within5), "--"),
          ifelse(nrow(h) > 0, e2s(median(h$r2_tie)), "--"), ifelse(nrow(h) > 0, e2s(median(h$r2_out)), "--"))
}
f2s <- function(x) sprintf("%.2f", x); pcs <- function(x) sprintf("%.0f\\%%", 100 * x)
e2s <- function(x) ifelse(x < 0.00005, "$<$0.0001", ifelse(x < 0.0005, sprintf("%.4f", x), sprintf("%.3f", x)))
SUBS <- list(c("Pooled studies", "main"), c("Outcome counts ongoing years", "outcome_overlap"), c("Expanded set", "both"))
sub_lines <- function(dd, gg) unlist(lapply(SUBS, function(z) {
  if (z[2] == "both") sub_row(dd, gg, z[1]) else sub_row(dd %>% filter(set == z[2]), gg %>% filter(set == z[2]), z[1]) }))
stopifnot(sum(D$set == "main") > 0, sum(D$set == "outcome_overlap") > 0)
writeLines(c("\\begin{tabular}{@{}lccccccc@{}}", "\\toprule",
  " & & Coefficient & Standard error & Lose significance & Placebo within & Partial $R^2$ & Partial $R^2$ \\\\",
  "Subset & Studies & ratio & ratio & or flip sign & 5\\% of baseline & key variable & outcome \\\\", "\\midrule",
  "\\multicolumn{8}{@{}l}{\\textit{All studies}} \\\\", sub_lines(D, G), "\\addlinespace",
  "\\multicolumn{8}{@{}l}{\\textit{Tie studies}} \\\\", sub_lines(D %>% filter(family == "TIE"), G %>% filter(family == "TIE")),
  "\\bottomrule", "\\end{tabular}"), file.path(APP, "tableA_sens_subsets.tex"))
cat("\nSubset table (both history controls)\n"); cat(readLines(file.path(APP, "tableA_sens_subsets.tex")), sep = "\n")

saveRDS(list(split_all = T_all, split_tie = T_tie, split_additive = T_add, split_interaction = T_int, n_full = length(keepF), strength_all = S_all, strength_tie = S_tie, compare_all = C_all, compare_tie = C_tie,
             dial = T_dial, n_pooled = length(keep), draws = O[[FULL[1]]]$draws), file.path(EST, "sens_summary.rds"))

# ---- tables ----
f2 <- function(x) sprintf("%.2f", x); pc <- function(x) sprintf("%.0f\\%%", 100 * x)
# [2026-11-14] placebo dial folded into the split table (replaces fig_sweep_dial in the appendix): one text row of medians by lambda
stopifnot(!is.null(T_dial), isTRUE(all.equal(T_dial$lambda, c(0, .25, .5, .75, 1))))
dial_vals <- sprintf("%.2f", T_dial$med_ratio)
dial_rows <- c("\\addlinespace",
  sprintf("\\multicolumn{10}{l}{\\textit{Placebo dial, both history controls (%d pre-selected tie studies)}} \\\\", length(dial_ids)),
  sprintf("\\multicolumn{10}{l}{Median coefficient ratio at $\\lambda = 0, 0.25, 0.5, 0.75, 1$ is %s, %s, %s, %s, and %s.} \\\\",
          dial_vals[1], dial_vals[2], dial_vals[3], dial_vals[4], dial_vals[5]))
# [2026-11-14] "Placebo within 5%" column added (replaces the appendix placebo scatter, fig_sweep_placebo)
row_split <- function(t, lab) c(sprintf("\\multicolumn{10}{l}{\\textit{%s}} \\\\", lab),
  sprintf("%s & %d & %s & %s & %s & %s & %d of %d & %s & %s & %s \\\\", gsub("\n", " ", SPEC[t$spec]), t$studies, f2(t$med_coef_ratio), pc(t$share_shrink),
          f2(t$med_se_ratio), pc(t$share_se_up), t$lose_sig, t$sig_baseline, ifelse(t$placebo_n > 0, f2(t$med_placebo_ratio), "--"),
          ifelse(t$placebo_n > 0, pc(t$share_placebo_within5), "--"), ifelse(t$placebo_n > 0, as.character(t$placebo_n), "--")))
writeLines(c("\\begin{tabular}{lccccccccc}", "\\toprule",
  " & Studies & Coefficient & Share & Standard error & Share with & Lose & Matched placebo & Placebo within & Placebo \\\\",
  "History controls added & & ratio & shrinking & ratio & larger SE & significance & coefficient ratio & 5\\% of baseline & studies \\\\", "\\midrule",
  row_split(T_all, "All studies"), "\\addlinespace", row_split(T_add, "Additive key term"), "\\addlinespace",
  row_split(T_int, "Interaction or product key term"), "\\addlinespace", row_split(T_tie, "Tie studies"), dial_rows, "\\bottomrule", "\\end{tabular}"),
  file.path(APP, "tableA_sens_split.tex"))   # placebo column: studies with an additive key term only
e2 <- function(x) ifelse(x < 0.00005, "$<$0.0001", ifelse(x < 0.0005, sprintf("%.4f", x), sprintf("%.3f", x)))
row_str <- function(t, lab) c(sprintf("\\multicolumn{4}{l}{\\textit{%s}} \\\\", lab),
  sprintf("%s & %d & %s & %s \\\\", gsub("\n", " ", SV[t$var]), t$studies, e2(t$med_r2_tie), e2(t$med_r2_out)))
writeLines(c("\\begin{tabular}{lccc}", "\\toprule", " & Studies & Partial $R^2$ with & Partial $R^2$ with \\\\",
  # [2026-11-14] audit: the first block is the studies with an additive key term (strength is not computed for interaction terms)
  "Variable & & the key variable & the outcome \\\\", "\\midrule", row_str(S_all, "Additive key term"), "\\addlinespace",
  row_str(S_tie, "Tie studies"), "\\bottomrule", "\\end{tabular}"), file.path(APP, "tableA_sens_strength.tex"))

# ---- figures ----
TH <- theme_minimal(base_size = 14) + theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold", size = 16),   # larger text [2026-11-05]
                                            legend.position = "bottom")
save2 <- function(p, stem, w, h) { ggsave(file.path(APP, paste0(stem, ".png")), p, width = w, height = h, dpi = 200, bg = "white")
  ggsave(file.path(APP, paste0(stem, ".pdf")), p, width = w, height = h, device = PDF_DEV) }
D$spec_f <- factor(SPEC[D$spec], rev(SPEC))
D$status <- factor(ifelse(!D$sig0, "Not significant at baseline", ifelse(D$lose, "Loses significance", ifelse(D$flip, "Significant, sign flips", "Stays significant"))),
                   c("Stays significant", "Loses significance", "Significant, sign flips", "Not significant at baseline"))
PAL <- c("Stays significant" = "#2a78d6", "Loses significance" = "#eb6834", "Significant, sign flips" = "#8a2be2", "Not significant at baseline" = "grey65")
strip <- function(d, v, title, xlab, lim, brk, labs) {
  d$x <- pmin(pmax(d[[v]], lim[1]), lim[2]); med <- d %>% group_by(spec_f) %>% summarise(m = median(.data[[v]]), .groups = "drop")
  ggplot(d, aes(x, spec_f)) + geom_vline(xintercept = 1, linetype = 2, colour = "grey40") +
    geom_jitter(aes(colour = status), height = .17, width = 0, alpha = .7, size = 1.8) +
    geom_segment(data = med, aes(x = m, xend = m, y = as.numeric(spec_f) - .3, yend = as.numeric(spec_f) + .3), linewidth = 1.1) +
    geom_text(data = med, aes(x = m, y = as.numeric(spec_f) + .38, label = sprintf("median %.2f", m)), size = 4.2) +
    scale_colour_manual(values = PAL, name = NULL, drop = TRUE, guide = guide_legend(nrow = 1, override.aes = list(size = 3))) +   # [2026-11-12] no empty keys
    scale_x_continuous(limits = lim + c(-.05, .05), breaks = brk, labels = labs) + labs(title = title, x = xlab, y = NULL) + TH
}
pS <- strip(D, "r", sprintf("A. Key coefficient (%d studies)", length(keep)), "Coefficient after / coefficient before", c(-0.5, 2),
            c(-0.5, 0, .5, 1, 1.5, 2), c("<= -0.5", "0", "0.5", "1\n(no change)", "1.5", ">= 2")) +
  strip(D, "r_se", "B. Standard error of the key coefficient", "Standard error after / standard error before", c(0.5, 1.5),
        c(.5, .75, 1, 1.25, 1.5), c("<= 0.5", "0.75", "1\n(no change)", "1.25", ">= 1.5")) + theme(axis.text.y = element_blank()) +
  plot_layout(guides = "collect") & theme(legend.position = "bottom")
save2(pS, "fig_sweep_split", 11, 5.2)

G$var_f <- factor(SV[G$var], SV); G3 <- G %>% filter(var %in% c("history_both", "joint_democracy", "capability_ratio"))
FL <- 1e-6; cl <- function(x) pmax(x, FL)   # partial R2 below 1e-6 drawn at the axis floor
mG <- G3 %>% group_by(var_f) %>% summarise(r2_tie = median(r2_tie), r2_out = median(r2_out), .groups = "drop")
BR <- 10^(-6:0); BL <- c("<= 0.0001%", "0.001%", "0.01%", "0.1%", "1%", "10%", "100%")
pG <- ggplot(G3, aes(cl(r2_tie), cl(r2_out), colour = var_f, shape = var_f)) + geom_point(alpha = .55, size = 1.9) +
  geom_point(data = mG, size = 5, stroke = 1.4, fill = "white", show.legend = FALSE) +
  scale_x_log10(limits = c(FL, 1), breaks = BR, labels = BL) + scale_y_log10(limits = c(FL, 1), breaks = BR, labels = BL) +
  scale_colour_manual(values = setNames(c("#eb6834", "#2a78d6", "grey45"), SV[c("history_both", "joint_democracy", "capability_ratio")]), name = NULL) +
  scale_shape_manual(values = setNames(c(16, 17, 15), SV[c("history_both", "joint_democracy", "capability_ratio")]), name = NULL) +
  labs(title = sprintf("Strength of conflict history as a confounder (%d studies)", length(keepF)),
       x = "Partial R2 with the tie variable, net of the study's controls",
       y = "Partial R2 with the outcome,\nnet of the tie variable and controls") + TH
save2(pG, "fig_sweep_strength", 8, 6)

P <- D %>% filter(scope == "full") %>% mutate(x = pmin(pmax(r_plac, -0.5), 2), y = pmin(pmax(r, -0.5), 2), spec_c = factor(gsub("\n", " ", SPEC[spec]), gsub("\n", " ", SPEC)))
pP <- ggplot(P, aes(x, y)) + geom_hline(yintercept = 1, linetype = 2, colour = "grey40") + geom_vline(xintercept = 1, linetype = 2, colour = "grey40") +
  geom_abline(slope = 1, intercept = 0, colour = "grey75") + geom_point(aes(colour = status), alpha = .75, size = 1.9) + facet_wrap(~spec_c) +
  scale_colour_manual(values = PAL, name = NULL, drop = TRUE, guide = guide_legend(nrow = 1, override.aes = list(size = 3))) +
  scale_x_continuous(limits = c(-0.55, 2.05), breaks = c(0, 1, 2), labels = c("0", "1\n(no change)", ">= 2")) +
  scale_y_continuous(limits = c(-0.55, 2.05), breaks = c(-0.5, 0, .5, 1, 1.5, 2), labels = c("<= -0.5", "0", "0.5", "1", "1.5", ">= 2")) +
  # [2026-11-13] in-plot title removed (study count reported in the caption/text): sprintf("... (%d studies)", length(keepF))
  labs(x = "Coefficient ratio with the matched placebo\n(same association with the key variable, no information about the outcome)",
       y = "Coefficient ratio with\nconflict history") + TH
save2(pP, "fig_sweep_placebo", 11, 4.8)

if (nrow(L)) {
  pL <- ggplot(L %>% mutate(r = pmin(pmax(r, -0.5), 2)), aes(lambda, r)) + geom_hline(yintercept = 1, linetype = 2, colour = "grey40") +
    geom_line(aes(group = id), colour = "grey60", alpha = .7) + geom_line(data = T_dial, aes(lambda, med_ratio), linewidth = 1.2, colour = "#eb6834") +
    geom_point(data = T_dial, aes(lambda, med_ratio), size = 2.5, colour = "#eb6834") +
    scale_x_continuous(breaks = c(0, .25, .5, .75, 1), labels = c("0\n(placebo)", "0.25", "0.5", "0.75", "1\n(conflict history)")) +
    scale_y_continuous(limits = c(-0.55, 2.05), breaks = c(-0.5, 0, .5, 1, 1.5, 2), labels = c("<= -0.5", "0", "0.5", "1", "1.5", ">= 2")) +
    # [2026-11-12] no in-plot title (caption carries it; the old title was clipped); wider right margin so "1 (conflict history)" fits
    labs(title = NULL, x = "Share of the history remainder that is real (lambda)", y = "Coefficient after / coefficient before") + TH +
    theme(plot.margin = margin(8, 40, 8, 8))
  save2(pL, "fig_sweep_dial", 7.5, 5)
  # [2026-11-12] ggplot object without the title, for the appendix grid with the war-sample figure (03_figures/06_*.R)
  saveRDS(pL + labs(title = NULL), file.path(APP, "fig_sweep_dial_plot.rds"))
}
cat("\n93_sensitivity_summary.R complete\n")
