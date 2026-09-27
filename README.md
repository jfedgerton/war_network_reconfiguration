# War and Network Reconfiguration: Replication Code

Replication code for a manuscript on how interstate conflict reconfigures the alliance and IGO ties among states, and how that reconfiguration contaminates downstream analyses that treat network ties as exogenous predictors of conflict.

## Research design

1. **System level (1816–2014).** Year-over-year Jensen–Shannon (JS) divergence of the multiplex (alliance + IGO) edge distribution, regressed on the lagged count of militarized interstate disputes (MIDs). OLS with Newey–West HAC errors, VAR/Granger causality, and forecast-error variance decompositions.
2. **State level.** Weighted-Jaccard churn in each state's multilayer neighborhood, regressed on MID onsets and tie-weighted partner war exposure. OLS and two-way fixed effects.
3. **Contamination.** Dyad-year models of MID onset on lagged tie strength, conditioning on conflict history at t−2…t−6.

KL divergence (the earlier headline system-level DV) is reported as an appendix robustness check.

## Structure

```
code/
  01_data_cleaning/   build network layers; compute JS/KL, churn, Leiden communities
  02_analysis/        01–06 main text; 07–23 appendix robustness
  03_figures/         multiverse summaries, case-study map, mediation table
  regen_appendix_figures.R
data/
  raw/                not tracked -- see data/raw/DATA_DOWNLOAD_GUIDE.md
  processed/          not tracked -- produced by code/01_data_cleaning
results/              not tracked -- LaTeX tables written by the pipeline
docs/codebook.docx
run_all.R             master pipeline
renv.lock             pinned R package versions (R 4.5.0)
```

## Running

```r
renv::restore()
source("run_all.R")
```

Seeds are fixed at 123. Paths resolve from the repository root via `here` (the `.here` file marks the root).

Script notes:

| Script | Main-text / appendix output | Notes |
|---|---|---|
| `02_analysis/01` | Tables 1–3 (system) | JS headline; KL copies kept in `results/appendix/kl_canonical/` |
| `02_analysis/02` | Table 4, Figure 4 (state) | |
| `02_analysis/03` | Figure 8 (contamination) | needs `peacesciencer::download_extdata()` |
| `02_analysis/04–06` | Figures 5–6 (placebos, mediation) | |
| `02_analysis/07` | multiverse | ~200k specifications; several hours |
| `02_analysis/10` | random / sliding year drops | ~10 min |
| `02_analysis/14` | IGO membership-share windows | needs `igoR` |
| `02_analysis/24` | Appendix: CINC without mobilization; pre/post-war windows | ~7 min (1,000 bootstrap draws) |
| `02_analysis/22` | threshold / Markov-switching | needs `tsDyn`, `MSwM`; not in the manuscript |

## Change log

- **2026-09-26.** Headline system-level DV switched from KL to JS divergence in scripts 01, 04, 05, 10, 11, 13, 14, 15, 16, 22 and `export_figure_tables_appendix.R`. Scripts 13 and 14, which compute the divergence directly, now compute JS. The multiverse (07) already contains both measures and was not rerun.
- **2026-09-26.** Added `24_robustness_cinc_nonmobilization.R` (realized-power mediation with non-mobilization CINC and pre/post-war capability windows).
