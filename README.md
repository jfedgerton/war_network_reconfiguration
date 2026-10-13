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

The replication sweep of published studies (Stage 5 of `run_all.R`, off by default) is not included in this repository.

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
| `02_analysis/07` | Multiverse | layer combos x weightings x conflict measures x churn measures x models x belligerent-dyad rules x extra controls (GDP, strength, third-party war exposure, community war exposure); several hours |
| `02_analysis/10` | Random / sliding year drops | ~10 min |
| `02_analysis/14` | IGO membership-share windows | needs `igoR` |
| `02_analysis/23` | Wild cluster bootstrap | |
| `02_analysis/24` | CINC without mobilization; pre/post-war windows | ~7 min (1,000 bootstrap draws) |
| `02_analysis/25–30` | Strength and GDP controls, non-participant KL, structural time series, out-of-sample prediction, winners and losers | |
| `02_analysis/22` | Threshold / Markov-switching | needs `tsDyn`, `MSwM`; not in the manuscript |

## Change log

- **2026-10-13.** Multiverse (07): third-party war-exposure controls (tie-weighted, actual ties at t−2, same-dispute partners excluded; community version), extra controls run on the headline cell, weighted x non-raw duplicate cells removed, state grid raw only, fit-status logging, per-control meta-regressions. 06 reports multiverse shares on baseline specifications.
- **2026-10-12.** Headline state measure is fatal MID onset at t−1 (full COW universe); partner war exposure and H3 removed; peace-dyad panels; belligerent-dyad rule dimension and IV-timing x FE grid in 07; scripts 17 and 18 retired.
- **2026-09-26.** Headline system-level DV switched from KL to JS divergence.
- **2026-09-26.** Added `24_robustness_cinc_nonmobilization.R`.
