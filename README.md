# War and Network Reconfiguration: Replication Code

Replication code for a manuscript on how interstate conflict reconfigures the alliance and IGO ties among states, and how that reconfiguration contaminates downstream analyses that treat network ties as exogenous predictors of conflict.

## Research design

1. **System level (1816–2014).** Year-over-year Jensen–Shannon (JS) divergence of the multiplex (alliance + IGO) edge distribution, regressed on the lagged count of fatal militarized interstate disputes (MIDs), over all dyads and over peace dyads. OLS with Newey–West HAC errors, VAR/Granger causality, and forecast-error variance decompositions.
2. **State level.** Weighted-Jaccard churn in each state's multilayer neighborhood, regressed on fatal MID onsets in t−1, over all dyads and over peace dyads (fatal-MID opponents in t−1 dropped). OLS and two-way fixed effects; community-clustered and community fixed-effects versions in the appendix.
3. **Contamination.** Dyad-year models of MID onset on lagged tie strength, conditioning on conflict history at t−2…t−6.

KL divergence (the earlier headline system-level DV) is reported as an appendix robustness check.

## Structure

```
code/
  01_data_cleaning/   build network layers; compute JS/KL, churn, Leiden communities
  02_analysis/        01–06 main text; 07–30 appendix robustness
  03_figures/         multiverse summaries, case-study map, mediation table
data/
  raw/                not tracked -- see data/raw/DATA_DOWNLOAD_GUIDE.md
  processed/          not tracked -- produced by code/01_data_cleaning
results/              not tracked -- LaTeX tables written by the pipeline
docs/codebook.docx
run_all.R             master pipeline
renv.lock             pinned R package versions
```

The replication sweep of published studies is in `code/04_replication_sweep/` (Stage 5 of `run_all.R`, off by default; `RUN_STAGE5=TRUE` runs it). The authors' replication data are not redistributed: `code/04_replication_sweep/data_sources.csv` lists each source, and `01_fetch_data.R` downloads them into `data/raw/replication_sweep/`. The appendix sensitivity analyses run afterwards with `Rscript code/04_replication_sweep/run_sensitivity.R`.

## Running

```r
renv::restore()
source("run_all.R")
```

Seeds are fixed at 123. Paths resolve from the repository root via `here` (the `.here` file marks the root).

Script notes:

| Script | Main-text / appendix output | Notes |
|---|---|---|
| `02_analysis/01` | System tables (JS headline; peace-dyad panel) | KL copies kept in `results/appendix/kl_canonical/` |
| `02_analysis/02` | State table (fatal MID onset, t−1), local projections | |
| `02_analysis/03` | Contamination ladder | needs `peacesciencer::download_extdata()` |
| `02_analysis/04–06` | Placebos, mediation, multiverse shares | |
| `02_analysis/07` | Multiverse | layer combos x weightings x conflict measures x churn measures x models x belligerent-dyad rules x extra control (log GDP); several hours |
| `02_analysis/10` | Random / sliding year drops | ~10 min |
| `02_analysis/14` | IGO membership-share windows | needs `igoR` |
| `02_analysis/23` | Wild cluster bootstrap | |
| `02_analysis/24` | CINC without mobilization; pre/post-war windows | ~7 min (1,000 bootstrap draws) |
| `02_analysis/32` | Multiverse permutation inference: reruns 07 on permuted conflict (500 draws, 3 fixed layer combos) and compares support shares with the observed draw | several hours; `RUN_MV_PERM=FALSE` skips |
| `02_analysis/25–30` | Strength and GDP controls, non-participant KL, structural time series, out-of-sample prediction, winners and losers | |
| `02_analysis/22` | Threshold / Markov-switching | needs `tsDyn`, `MSwM`; not in the manuscript |

## Script for each table and figure

Moved here from the appendix (Replication Materials) on 2026-11-09.

| Output | Script(s) |
|---|---|
| Main text, Table 1 | `02_analysis/31_table1_sumstats.R` |
| Main text, Tables 2–3 | `02_analysis/01_system_level_analysis.R` |
| Main text, Table 4 and Figure 4 | `02_analysis/02_state_level_analysis.R` |
| Main text, Figure 5 | `02_analysis/04_robustness_placebos.R`, `02_analysis/06_sensitivity_plots.R` |
| Main text, Figure 6 | `03_figures/02_case_study_map.R` |
| Main text, Figure 7 | `02_analysis/03_contamination_demo.R` |
| Main text, Figure 8 | `04_replication_sweep/run_sweep.R`, `04_replication_sweep/91_figures.R` |
| Main text, Table 5 | hand-coded from the sources in Appendix Section 7 |
| Main text, Figures 1–3 | TikZ diagrams in the manuscript source |
| Appendix §1 (multiverse) | `02_analysis/01`, `02` (tables), `02_analysis/07_robustness_multiverse.R` (multiverse figures), `02_analysis/32_multiverse_permutation.R` (permutation figure) |
| Appendix §2 (contamination robustness) | `02_analysis/19_contamination_multiverse.R`, `20_contamination_falsification.R`, `21_contamination_specification.R` |
| Appendix §3 (replication sweep) | `04_replication_sweep/run_sweep.R`, `run_sensitivity.R` (through `legacy/hist/sens_helper.R`), `91_figures.R`, `93_sensitivity_summary.R`, `94_study_flow_table.R`; `registry.csv` records in `pool_status` whether each study is pooled and, if not, why |
| Appendix §4 (fatal MIDs) | `02_analysis/11_robustness_war_samples.R`, `12_robustness_hostlev.R` |
| Appendix §5 (IGO windows) | `02_analysis/14_igo_membership_sliding.R` |
| Appendix §6 (temporal) | `02_analysis/10_robustness_year_drops.R` (random and sliding-window drops), `04_robustness_placebos.R` (permutation tests), `13_robustness_functional_form.R` (long differences) |
| Combined appendix figures (§§2, 3, 4, 6) | `03_figures/05_appendix_compact_figures.R` (redraws saved results; no re-estimation) |

Results in the replication archive only: full regression tables for the main-text figures (`02_analysis/export_figure_tables_appendix.R`), the wild cluster bootstrap figure (`02_analysis/23_wild_cluster_bootstrap.R`), Granger tests by outcome measure (`03_figures/04_granger_forward_reverse.R`), and the lower- and upper-only IGO windows (`02_analysis/14_igo_membership_sliding.R`).

## Change log

- **2026-11-05.** Multiverse (07): the network-strength control (post-treatment) and the two third-party war-exposure controls removed; log GDP is the only extra control. Granger forward/reverse figure rebuilt by `03_figures/04_granger_forward_reverse.R`.
- **2026-10-13.** Multiverse (07): third-party war-exposure controls (tie-weighted, actual ties at t−2, same-dispute partners excluded; community version), extra controls run on the headline cell, weighted x non-raw duplicate cells removed, state grid raw only, fit-status logging, per-control meta-regressions. 06 reports multiverse shares on baseline specifications.
- **2026-10-12.** Headline state measure is fatal MID onset at t−1 (full COW universe); partner war exposure and H3 removed; peace-dyad panels; belligerent-dyad rule dimension and IV-timing x FE grid in 07; scripts 17 and 18 retired.
- **2026-09-26.** Headline system-level DV switched from KL to JS divergence.
- **2026-09-26.** Added `24_robustness_cinc_nonmobilization.R`.
