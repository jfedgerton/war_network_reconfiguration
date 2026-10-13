###############################################################################
# 24_robustness_cinc_nonmobilization.R
#
# Responds to R3 (APSR, Sept 2026): Delta-CINC is a poor test of the
# "realized power redistribution" alternative because CINC includes military
# personnel and expenditure, which jump during wars (mobilization) and fall
# afterwards (demobilization). This script re-runs the state-level mediation
# from 05_robustness_alt_mechanisms.R with capability measures that strip out
# mobilization and/or bracket the war.
#
# Capability indices (built from COW NMC v6 via peacesciencer::cow_nmc):
#   cinc_full    : mean of the six component shares (replicates COW CINC)
#   cinc_mob     : mean of the two mobilization shares (milex, milper)
#   cinc_nonmob  : mean of the four latent-capability shares
#                  (irst, pec, tpop, upop) -- "realized redistribution"
#                  without wartime mobilization
#
# Mediators (state-year, TWFE design identical to 05):
#   (1) d_cinc_full        : cinc_full[t] - cinc_full[t-1]   (replicates 05)
#   (2) d_cinc_mob         : cinc_mob[t]  - cinc_mob[t-1]    (diagnostic)
#   (3) d_cinc_nonmob      : cinc_nonmob[t] - cinc_nonmob[t-1]
#   (4) win_cinc_nonmob    : mean(cinc_nonmob[t+3..t+5]) - mean(cinc_nonmob[t-3..t-1])
#   (5) abswin_cinc_nonmob : |win_cinc_nonmob|  (redistribution in either direction)
#   (6) win_cinc_full      : same window, full CINC (post-demobilization)
#
# The window mediators (4)-(6) compare pre-war capabilities with capabilities
# three to five years later, after demobilization. They are measured after the
# churn outcome at t, so they are a sensitivity check on whether conditioning
# on the eventual redistribution absorbs the conflict -> churn association,
# not a sequential mediation estimate.
#
# Outputs:
#   data/processed/robustness_cinc_nonmob_mediation.csv
#   results/appendix/tableA_cinc_nonmob.tex
#   code/03_figures/output/manuscript/fig_ms_cinc_nonmob_mediation.png
###############################################################################

library(tidyverse)
library(fixest)
library(peacesciencer)
select <- dplyr::select

set.seed(123)
N_BOOT <- 1000

OUTPUT_DIR <- here::here("data", "processed")
TABLE_DIR  <- here::here("results", "appendix")
FIG_DIR    <- here::here("code", "03_figures", "output", "manuscript")
dir.create(TABLE_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG_DIR,   recursive = TRUE, showWarnings = FALSE)

cat("R", R.version.string, "| fixest", as.character(packageVersion("fixest")),
    "| peacesciencer", as.character(packageVersion("peacesciencer")), "\n")

# ============================================================
# 1. Capability indices from COW NMC component shares
# ============================================================
cat("Building capability indices from cow_nmc...\n")
nmc <- peacesciencer::cow_nmc |>
  select(ccode, year, milex, milper, irst, pec, tpop, upop, cinc) |>
  mutate(across(c(milex, milper, irst, pec, tpop, upop),
                ~ if_else(.x < 0, NA_real_, as.numeric(.x))))   # COW codes missing as -9

nmc <- nmc |>
  group_by(year) |>
  mutate(across(c(milex, milper, irst, pec, tpop, upop),
                ~ .x / sum(.x, na.rm = TRUE), .names = "sh_{.col}")) |>
  ungroup()

share_mean <- function(df, cols) {
  m <- as.matrix(df[, cols]); out <- rowMeans(m, na.rm = TRUE)
  out[is.nan(out)] <- NA_real_; out
}
nmc$cinc_full   <- share_mean(nmc, c("sh_milex", "sh_milper", "sh_irst", "sh_pec", "sh_tpop", "sh_upop"))
nmc$cinc_mob    <- share_mean(nmc, c("sh_milex", "sh_milper"))
nmc$cinc_nonmob <- share_mean(nmc, c("sh_irst", "sh_pec", "sh_tpop", "sh_upop"))

# Validation: the rebuilt full index must reproduce COW's published CINC.
ok <- complete.cases(nmc$cinc_full, nmc$cinc)
r_full <- cor(nmc$cinc_full[ok], nmc$cinc[ok])
cat(sprintf("  cor(rebuilt CINC, COW CINC) = %.4f over %d state-years\n", r_full, sum(ok)))
stopifnot(r_full > 0.99)
r_nonmob <- cor(nmc$cinc_nonmob, nmc$cinc_full, use = "complete.obs")
cat(sprintf("  cor(non-mobilization CINC, full CINC) = %.4f\n", r_nonmob))

# Year-over-year and pre/post window changes, within state
roll_mean <- function(x, offsets) {
  m <- sapply(offsets, function(k) if (k < 0) dplyr::lag(x, -k) else dplyr::lead(x, k))
  n_ok <- rowSums(!is.na(m))
  out <- rowMeans(m, na.rm = TRUE); out[n_ok < 2] <- NA_real_   # need >= 2 of 3 years
  out
}
nmc <- nmc |>
  arrange(ccode, year) |>
  group_by(ccode) |>
  mutate(
    d_cinc_full     = cinc_full   - dplyr::lag(cinc_full),
    d_cinc_mob      = cinc_mob    - dplyr::lag(cinc_mob),
    d_cinc_nonmob   = cinc_nonmob - dplyr::lag(cinc_nonmob),
    win_cinc_nonmob = roll_mean(cinc_nonmob, 3:5) - roll_mean(cinc_nonmob, -3:-1),
    win_cinc_full   = roll_mean(cinc_full,   3:5) - roll_mean(cinc_full,   -3:-1)
  ) |>
  ungroup() |>
  mutate(abswin_cinc_nonmob = abs(win_cinc_nonmob))

stopifnot(!anyDuplicated(nmc[, c("ccode", "year")]))

# ============================================================
# 2. State-level panel (identical construction to 05)
# ============================================================
cat("Loading state panel...\n")
state_churn <- readRDS(file.path(OUTPUT_DIR, "state_network_churn.rds"))
dyad_data   <- readRDS(file.path(OUTPUT_DIR, "dyad_year_network_layers.rds"))

ML <- c("layer_offensive", "layer_defensive", "layer_igo")
mat <- as.matrix(dyad_data[, ML, drop = FALSE]); mat[is.na(mat)] <- 0
dyad_data$main_raw <- rowSums(mat)
dyad_data <- dyad_data |>
  group_by(year) |>
  mutate(main_weight = if (max(main_raw, na.rm = TRUE) > 0) main_raw / max(main_raw, na.rm = TRUE) else 0) |>
  ungroup()

state_attrs <- dyad_data |>
  select(ccode = ccode1, year, cinc = cinc1, major_power = cowmaj1) |>
  distinct(ccode, year, .keep_all = TRUE) |>
  mutate(major_power = as.integer(coalesce(major_power, 0L) == 1))
state_churn <- state_churn |> left_join(state_attrs, by = c("ccode", "year"))

state_years <- create_stateyears(subset_years = min(state_churn$year):max(state_churn$year))
attr(state_years, "ps_data_type") <- "state_year"
attr(state_years, "ps_system")    <- "cow"
state_years <- state_years |> add_democracy()
polity <- state_years |>
  select(ccode, year, polity2) |>
  mutate(democracy = as.integer(polity2 >= 6), autocracy = as.integer(polity2 <= -6))
state_churn <- state_churn |> left_join(polity, by = c("ccode", "year"))

partner_sets <- dyad_data |>
  filter(main_weight > 0) |>
  select(ccode1, ccode2, year, w = main_weight) |>
  bind_rows(dyad_data |> filter(main_weight > 0) |>
              select(ccode1 = ccode2, ccode2 = ccode1, year, w = main_weight))
at_war <- dyad_data |>
  filter(cowmidongoing == 1) |>
  select(ccode1, ccode2, year) |>
  pivot_longer(c(ccode1, ccode2), values_to = "ccode") |>
  distinct(ccode, year) |>
  mutate(at_war = 1L)
pwe <- partner_sets |>
  mutate(year_next = year + 1L) |>
  left_join(at_war, by = c("ccode2" = "ccode", "year_next" = "year")) |>
  mutate(at_war = replace_na(at_war, 0L)) |>
  group_by(ccode1, year_next) |>
  summarise(partner_war_frac_weighted = sum(at_war * w) / sum(w), .groups = "drop") |>
  rename(ccode = ccode1, year = year_next)
state_churn <- state_churn |>
  left_join(pwe, by = c("ccode", "year")) |>
  mutate(partner_war_frac_weighted = replace_na(partner_war_frac_weighted, 0),
         major_power = replace_na(major_power, 0),
         log_churn = log1p(churn),
         log_mid_onset_count = log1p(mid_onset_count))
# [2026-10-12] Headline state IVs are onset counts at t-1, taken from the full COW state-year universe in
# dyad_data (same construction as 01_data_cleaning/02_compute_network_measures.R and 02). The churn panel
# holds transition rows only, so lagging the previous ROW would drop every state's first transition year.
# A t-1 value is NA only if the state was not in the system in t-1 (after a COW membership gap).
onset_univ <- dplyr::bind_rows(dyad_data |> dplyr::distinct(ccode = ccode1, year),
                               dyad_data |> dplyr::distinct(ccode = ccode2, year)) |>
  dplyr::distinct() |>
  dplyr::left_join(
    dyad_data |>
      dplyr::filter(cowmidongoing == 1, ccode1 < ccode2) |>
      dplyr::select(ccode1, ccode2, year, cowmidfatal, cowmidonset) |>
      tidyr::pivot_longer(c(ccode1, ccode2), values_to = "ccode") |>
      dplyr::group_by(ccode, year) |>
      dplyr::summarise(u_onset = sum(cowmidonset == 1, na.rm = TRUE),
                       u_fatal_onset = sum(cowmidfatal == 1 & cowmidonset == 1, na.rm = TRUE),
                       .groups = "drop"),
    by = c("ccode", "year")) |>
  dplyr::mutate(u_onset = as.numeric(tidyr::replace_na(u_onset, 0L)),
                u_fatal_onset = as.numeric(tidyr::replace_na(u_fatal_onset, 0L)))
chk_univ <- dplyr::inner_join(state_churn |> dplyr::select(ccode, year, mid_onset_count, fatal_mid_onset_count),
                              onset_univ, by = c("ccode", "year"))
stopifnot(nrow(chk_univ) == nrow(state_churn), all(chk_univ$mid_onset_count == chk_univ$u_onset),
          all(chk_univ$fatal_mid_onset_count == chk_univ$u_fatal_onset))
rm(chk_univ)
state_churn <- state_churn |>
  dplyr::select(-dplyr::any_of(c("onset_lag1", "fatal_onset_lag1", "log_onset_lag1", "log_fatal_onset_lag1"))) |>
  dplyr::left_join(onset_univ |>
                     dplyr::transmute(ccode, year = year + 1L,
                                      onset_lag1 = u_onset, fatal_onset_lag1 = u_fatal_onset,
                                      log_onset_lag1 = log1p(u_onset), log_fatal_onset_lag1 = log1p(u_fatal_onset)),
                   by = c("ccode", "year"))
stopifnot(!anyDuplicated(state_churn[c("ccode", "year")]))

n_before <- nrow(state_churn)
state_panel <- state_churn |>
  left_join(nmc |> select(ccode, year, cinc_full, cinc_mob, cinc_nonmob,
                          d_cinc_full, d_cinc_mob, d_cinc_nonmob,
                          win_cinc_nonmob, abswin_cinc_nonmob, win_cinc_full),
            by = c("ccode", "year"))
stopifnot(nrow(state_panel) == n_before)   # join must not duplicate rows
cat(sprintf("  state-years: %d; with year-over-year non-mob change: %d; with window change: %d\n",
            nrow(state_panel), sum(!is.na(state_panel$d_cinc_nonmob)),
            sum(!is.na(state_panel$win_cinc_nonmob))))

# [2026-10-12] partner war exposure removed from the control set (H3 dropped from the paper).
CTRL <- "democracy + autocracy + cinc + major_power + is_isolated_t1"
# [2026-10-12] focal IV follows the fatal t-1 headline
IV   <- "log_fatal_onset_lag1"

# ============================================================
# 3. Baron-Kenny mediation, TWFE, cluster bootstrap over states (as in 05)
# ============================================================
bk_state <- function(data, iv, mediator, n_boot = N_BOOT) {
  f_tot <- as.formula(paste0("churn ~ ", iv, " + ", CTRL, " | ccode + year"))
  f_a   <- as.formula(paste0(mediator, " ~ ", iv, " + ", CTRL, " | ccode + year"))
  f_bc  <- as.formula(paste0("churn ~ ", iv, " + ", mediator, " + ", CTRL, " | ccode + year"))
  d <- data |> filter(!is.na(.data[[mediator]]), !is.na(.data[[iv]]))

  m_tot <- feols(f_tot, data = d, cluster = ~ ccode + year)
  m_a   <- feols(f_a,   data = d, cluster = ~ ccode + year)
  m_bc  <- feols(f_bc,  data = d, cluster = ~ ccode + year)

  c_tot <- coef(m_tot)[iv]; a <- coef(m_a)[iv]; a_se <- se(m_a)[iv]
  b <- coef(m_bc)[mediator]; b_se <- se(m_bc)[mediator]; c_prm <- coef(m_bc)[iv]
  ind <- a * b

  cl_idx <- split(seq_len(nrow(d)), d$ccode); n_cl <- length(cl_idx)
  ind_bt <- rep(NA_real_, n_boot)
  for (brep in seq_len(n_boot)) {
    draw <- sample.int(n_cl, n_cl, replace = TRUE)
    rows <- unlist(cl_idx[draw], use.names = FALSE)
    db <- d[rows, ]; db$ccode <- rep.int(seq_along(draw), lengths(cl_idx[draw]))
    ma <- tryCatch(feols(f_a,  data = db), error = function(e) NULL)
    mb <- tryCatch(feols(f_bc, data = db), error = function(e) NULL)
    if (is.null(ma) || is.null(mb)) next
    ind_bt[brep] <- coef(ma)[iv] * coef(mb)[mediator]
  }
  n_fail <- sum(is.na(ind_bt))
  if (n_fail > 0.05 * n_boot) stop(sprintf("%s: %d of %d bootstrap draws failed", mediator, n_fail, n_boot))
  qs <- quantile(ind_bt, c(0.025, 0.975), na.rm = TRUE)

  tibble(mediator = mediator, n_obs = nobs(m_tot),
         c_total = unname(c_tot), c_prime = unname(c_prm),
         a_path = unname(a), a_se = unname(a_se), b_path = unname(b), b_se = unname(b_se),
         indirect_ab = unname(ind), ci_low = unname(qs[1]), ci_high = unname(qs[2]),
         pct_mediated = unname(ind / c_tot), boot_failed = n_fail)
}

mediators <- tribble(
  ~col,                 ~label,                                                   ~label_tex,                                                        ~window,
  "d_cinc_full",        "Change in CINC, all six components (main text)",         "$\\Delta$ CINC, all six components (main text)",                 "$t$ vs.\\ $t-1$",
  "d_cinc_mob",         "Change in CINC, mobilization only (milex, milper)",      "$\\Delta$ CINC, mobilization only (milex, milper)",              "$t$ vs.\\ $t-1$",
  "d_cinc_nonmob",      "Change in CINC, non-mobilization (irst, pec, tpop, upop)", "$\\Delta$ CINC, non-mobilization (irst, pec, tpop, upop)",     "$t$ vs.\\ $t-1$",
  "win_cinc_nonmob",    "Change in non-mobilization CINC, pre/post window",       "$\\Delta$ non-mobilization CINC, pre/post window",               "$t{+}3..t{+}5$ vs.\\ $t{-}3..t{-}1$",
  "abswin_cinc_nonmob", "Absolute change, non-mobilization CINC, pre/post window", "$|\\Delta|$ non-mobilization CINC, pre/post window",           "$t{+}3..t{+}5$ vs.\\ $t{-}3..t{-}1$",
  "win_cinc_full",      "Change in CINC (all components), pre/post window",       "$\\Delta$ CINC (all components), pre/post window",               "$t{+}3..t{+}5$ vs.\\ $t{-}3..t{-}1$"
)

res <- list()
for (k in seq_len(nrow(mediators))) {
  cat(sprintf("  [%d/%d] %s ...\n", k, nrow(mediators), mediators$col[k]))
  res[[k]] <- bk_state(state_panel, IV, mediators$col[k]) |>
    mutate(label = mediators$label[k], label_tex = mediators$label_tex[k], window = mediators$window[k])
}
res <- bind_rows(res)
write_csv(res, file.path(OUTPUT_DIR, "robustness_cinc_nonmob_mediation.csv"))
print(res |> select(mediator, n_obs, c_total, a_path, a_se, indirect_ab, ci_low, ci_high, pct_mediated), width = 200)

# ============================================================
# 4. Table and figure
# ============================================================
f3 <- function(x) formatC(1e4 * x, format = "f", digits = 2)   # reported in units of 10^-4
fp <- function(x) sprintf("%.1f\\%%", 100 * x)
tex <- c(
  "\\begin{table}[htbp]", "\\centering", "\\footnotesize",
  "\\caption{Alternative capability measures for the realized-power mediation. State-level TWFE (state and year fixed effects), churn (not logged) on $\\log(1+\\text{fatal MID onset})_{t-1}$ with the full control set. Indirect effects are $a \\times b$ with 95\\% cluster-bootstrap intervals (1,000 draws, resampling states). The $a$ path, indirect effect, and interval are multiplied by $10^{4}$.}",
  "\\label{tab:cinc_nonmob}",
  "\\begin{tabular}{p{5.2cm}lrrrrr}", "\\toprule",
  "Mediator & Window & $N$ & $a$ path (SE) & Indirect $ab$ & 95\\% CI & \\% mediated \\\\", "\\midrule",
  sprintf("%s & %s & %s & %s (%s) & %s & [%s, %s] & %s \\\\",
          res$label_tex, res$window, format(res$n_obs, big.mark = ","),
          f3(res$a_path), f3(res$a_se), f3(res$indirect_ab), f3(res$ci_low), f3(res$ci_high),
          fp(res$pct_mediated)),
  "\\bottomrule", "\\end{tabular}",
  "\\par\\vspace{4pt}",
  "{\\scriptsize \\textit{Note:} Capability shares are rebuilt from COW NMC components. Non-mobilization CINC averages the iron and steel, energy, total population, and urban population shares; mobilization CINC averages military expenditure and personnel shares. Window measures compare mean capabilities three to five years after $t$ with the three years before $t$ (at least two of three years required). Because window measures are observed after the churn outcome at $t$, they test whether conditioning on the eventual redistribution absorbs the association rather than estimating sequential mediation.}",
  "\\end{table}")
writeLines(tex, file.path(TABLE_DIR, "tableA_cinc_nonmob.tex"))

plot_df <- res |> mutate(label = factor(label, levels = rev(mediators$label)))
p <- ggplot(plot_df, aes(x = indirect_ab, y = label)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
  geom_errorbar(aes(xmin = ci_low, xmax = ci_high), width = 0.2, orientation = "y", color = "#2C6CB0") +
  geom_point(size = 2.4, color = "#2C6CB0") +
  labs(x = "Indirect effect (a x b) with 95% cluster-bootstrap CI", y = NULL,
       title = "Realized-power mediators with and without mobilization") +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold", size = 11))
ggsave(file.path(FIG_DIR, "fig_ms_cinc_nonmob_mediation.png"), p, width = 8.5, height = 4, dpi = 400, bg = "white")
cat("Done: robustness_cinc_nonmob_mediation.csv, tableA_cinc_nonmob.tex, fig_ms_cinc_nonmob_mediation.png\n")
