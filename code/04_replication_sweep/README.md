# Stage 5 — Replication sweep (conflict history as an omitted confounder)

Re-estimates the headline model of each published study in `registry.csv` with the article's own estimator, then adds
conflict-history controls (dyadic and monadic MID onsets in t-1..t-5 and t-2..t-6) and two benchmarks (lagged CINC ratio
placebo; lagged joint democracy).

## Run
```r
Rscript code/04_replication_sweep/01_fetch_data.R   # once: authors' replication data -> data/raw/replication_sweep/
Rscript code/04_replication_sweep/run_sweep.R       # history tables, estimates, pooled statistics, figures
```
`01_fetch_data.R` stops with an error while any file is missing. The 37 files from PRIO (Journal of Peace Research
replication archives) cannot be scripted (the server redirects scripted requests to a login page) and are listed with
download links and the archive member to extract. Three replication packages the PI downloaded by hand (P9_045 UNL Digital Commons, P9_048 SAGE supplement, P9_064) are listed with `source_type = unknown` and are likewise reported for manual download. `run_sweep.R` checks every input file of the selected studies
(existence, size, md5) before estimating anything.

Environment flags for `run_sweep.R`: `RUN_APPENDIX_FAMILIES` (TRUE), `RERUN_SLOW` (FALSE), `ONLY_IDS` (comma list; skips
pooled statistics and figures), `RESUME` (FALSE), `CHECK_MD5` (TRUE).

## Files
| File | Role |
|---|---|
| `registry.csv` | Input. One row per study x key term: citation, outlet, DOI, family, tie type, key term, headline model, estimator, published b/SE, citations, flags (`study_headline_term`, `in_headline`, `named_illustration`, `slow`), script paths. **No results.** |
| `data_sources.csv`, `01_fetch_data.R` | Where every input file comes from (Dataverse DOI/file id, URL, PRIO archive, derived) with md5; downloader with md5 checks. `01_fetch_data.R` is the only step that writes under `data/raw/`. |
| `00_build_history.R` | Conflict-history lookup tables from GML MIDs (peacesciencer) -> `data/processed/replication_sweep/hist/history_tables.rds`. |
| `run_one.R` | Runs one study's ported script(s) in a fresh R process (or, with `--cached`, takes a slow study's cached estimates), extracts its numbers via `extract_map.csv`, checks them against `expected_est.rds`, and only then saves `results/replication_sweep/est_<id>.rds` (ID-keyed). |
| `run_sweep.R` | Preflight, loop over studies, then `90_summarise.R`, `91_figures.R` and `94_study_flow_table.R`. |
| `90_summarise.R` | Pooled statistics by family (aggregates only) -> `results/replication_sweep/summary_stats.rds`, including `$topline` (pooled headline terms significant at baseline that become null or flip sign with history, overall and by family, split by whether the authors already control for conflict history, per registry `author_history`). Also writes `results/replication_sweep/study_level_by_id.csv`, one row per pooled study keyed by study id only (no names or citations). |
| `run_sensitivity.R` | Appendix sensitivity analyses, run separately AFTER `run_sweep.R` (`Rscript code/04_replication_sweep/run_sensitivity.R`). Re-runs each study in `sens_studies.csv` with `SWEEP_SENS=1`, which activates the guarded hook in that study's script (`legacy/hist/sens_helper.R`): tie coefficient with dyadic history only, monadic history only and both; history's partial R2 next to the benchmarks; a matched placebo; the lambda dial for `dial = Y` studies. A study's result is kept only if the hook reproduces its verified baseline and with-history estimates. Flags: `SENS_WORKERS`, `SENS_DRAWS`, `ONLY_IDS`, `RESUME`. |
| `sens_studies.csv` | Studies with a hook. `scope = full` (84): single additive tie column with a logit, probit, GEE, rare-events logit or linear estimator; all analyses. `scope = split` (25): headline term is an interaction, product or squared term, or a constituent of one; three-way split and standard errors only (`sens_split()`), because the placebo and partial R2 need a tie column whose association with history is defined net of the controls. `dial = Y`: the 20 most-cited tie studies among the `full` ones, by `citations_openalex`. |
| `93_sensitivity_summary.R` | Pooled sensitivity tables and figures (aggregates only) -> `results/replication_sweep/sens_summary.rds`, `results/appendix/tableA_sens_*.tex`, `results/appendix/fig_sweep_*`. |
| `92_rebuild_expected.R` | Only after a deliberate change to the history controls: copies the with-history estimates of a run checked on the history-independent columns (`SWEEP_CHECK_COLS=b0,p0,placebo,jointdem`) into `expected_est.rds`. Steps are in the file header. |
| `91_figures.R` | `results/main/fig_contamination_ties.*`, `results/appendix/fig_contamination_{democracy,capitalism}.*`. |
| `94_study_flow_table.R` | `results/appendix/tab_study_flow.{tex,csv}` from `study_flow.csv` (article-level attrition; asserts the pooled articles match the registry). |
| `extract_map.csv` | Where each study's script writes each number (inside the temporary `_build/` folder). |
| `expected_est.rds` | Reference estimates (ID-keyed) used as a regression test: b0, b_hist, b_26 (relative tolerance `tol`), p0, p_hist, p_26 (1e-5), placebo and joint-democracy shifts (1e-2), n (exact). |
| `legacy/` | Ported study code. Paths come from `SWEEP_CODE`/`SWEEP_DATA`/`SWEEP_OUT`/`SWEEP_HIST`; estimation code is unchanged from the original reruns (see `PORTING_CONVENTIONS.md`). Folder names under `legacy/` are the original working names. Most scripts are ports that re-run the authors' replication data. `rebuild_*.R` scripts reconstruct the estimating data. 32 registered studies are protocol-B rebuilds (52, X074, X076, P8_004, P8_005, P8_012, P8_043, P8_048, P9_028-P9_031, P9_033-P9_036, P9_060, P9_061, P9_066-P9_068, P9_071, P9_074, P9_075, P9_078, P9_080, P9_082, P9_085, P9_087, P9_089, P9_093, P9_094; 18 of them pooled). They have no usable author analysis data, are rebuilt from public sources, and are judged by the round-2 rebuild rule (see Rules applied), not by exact reproduction. |
| `legacy/_slow/`, `SLOW_MODELS.md` | Entry scripts and requirements for the 16 slow studies (TERGM/SAOM/AME/latent space/ERGM network models, one Bayesian MCMC model, three Heckman probit models: one rebuild and two ports, one random-effects logit port, P9_092). Cached estimates in `results/replication_sweep/cache_slow/`, checked against `expected_est.rds` like any other study. |

## Disclosure design
Apart from the main-text figure (below), no file written by the pipeline pairs study names with results. Study IDs are neutral codes; estimates are saved by
study ID only; intermediate result files that carry citations live in `results/replication_sweep/_build/`, which is
cleared at the start of every run and deleted after each study, including after an error. Pooled statistics and figures
are computed in memory. A user who wants the study-level table must join `est_<id>.rds` to `registry.csv` deliberately.
The main-text figure `results/main/fig_contamination_ties.*` does name five studies: the five network-measure studies chosen
in advance by citation count (`named_illustration == "Y"`). The PI chose on 2026-10-10 to keep these five names on the figure.

## Rules applied
- One headline term per study (`study_headline_term`).
- `in_headline == "N"`: the row is estimated but not pooled. Reasons are in `methods_note`: headline term not significant
  at baseline at the paper's own threshold (`alpha_pub`), rebuilds that fail the round-2 rule, out-of-scope outcomes (coup and terrorism DVs, third-party intervention in civil wars), headline estimands that cannot be ported (e.g. a
  time-varying-coefficient model where only a constant-coefficient comparator could be run), methods papers' reanalyses of
  another study's variable, and ports that cannot use the article's estimator or data. The pooled set is therefore all
  significant at baseline (`90_summarise.R` asserts this).
- Monadic history of a dyad = both states' own MID onsets minus their onsets in disputes against each other (counted at
  each state's own start year). Monadic designs: dyadic history is replaced by history with the dispute's opponents where
  defined.
- `alpha_pub` (`registry.csv`, evidence in `alpha_pub_coding.csv`): 0.10 only when the paper uses one-tailed tests or the selected
  headline claim itself rests on p < .10; otherwise 0.05 (a .10 star in the legend alone does not qualify).
- Round-2 rebuild rule (PI, 2026-10-11): a protocol-B rebuild matches only if (1) the rebuilt baseline b0 has the published sign,
  (2) b0 is significant at the paper's `alpha_pub` whenever the published estimate is, and (3) b0 lies within 0.5x-2x of the
  published coefficient; otherwise `REBUILT_MISMATCH` and `in_headline = N`. Implemented once in
  `legacy/work_p8/round2_rule.R` (`rebuild_verdict()`, alpha read from `registry.csv`) and used by every helper and rebuild
  script. Ports of author data keep the REPRODUCED/CLOSE rule with no magnitude gate.
- Common-sample guard: the 1822-2011 common sample must keep more than 80% of complete cases. One named exception, P9_092
  (Millerd 2023) at 0.75, PI decision 2026-10-11, because its 2012-2014 rows lie outside the history window (see
  `PORTING_CONVENTIONS.md`).
- Author history control (`author_history`): any prior-conflict control in any stage or equation of the author's model counts.
- Benchmark shifts are stored as `100 * (b_bench - b0) / |b0|` for every study.
- Single-threaded BLAS during estimation for deterministic results.

## Requirements
R 4.5.3 (the reference values were produced with it; other R versions may differ in the last digits and fail the
1e-6 tolerance); fixest, haven, dplyr, geepack, lme4, sandwich, glmx (heteroskedastic probit, `hetglm`), peacesciencer 1.2.0, ggplot2, patchwork (pinned in
`renv.lock`); `run_sweep.R` calls `peacesciencer::download_extdata()` once if the package's external data are missing
(several ports read `cow_trade_ddy.rds`). Slow studies additionally need the package versions in `SLOW_MODELS.md`.
