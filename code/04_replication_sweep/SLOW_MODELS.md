# Slow studies (registry `slow == "Y"`)

These 16 studies (36 registry rows; network models, one Bayesian MCMC model, three Heckman probits and one random-effects
logit, P9_092) take between about 11 minutes and 11 hours each. P9_092 uses the named 0.75 common-sample guard exception (PI
2026-10-11; see PORTING_CONVENTIONS.md). Rebuild verdicts follow the round-2 rule in PORTING_CONVENTIONS.md. By default the pipeline uses the
cached estimates in `results/replication_sweep/cache_slow/est_<id>.rds`. With `RERUN_SLOW = TRUE` it runs the ported chain
below.

| id | Model | Rows | Entry script | Output read by `extract_map.csv` |
|---|---|---|---|---|
| X060 | SAOM (RSiena) | 2 | `_slow/X060.R` | `_slow/X060.csv` |
| X085 | SAOM (RSiena) | 1 | `_slow/X085.R` | `_slow/X085.csv` |
| X055 | TERGM (btergm bootstrap MPLE) | 3 | `_slow/X055.R` | `_slow/X055.csv` |
| X057 | TERGM (btergm bootstrap MPLE) | 1 | `_slow/X057.R` | `_slow/X057.csv` |
| N01 | TERGM (btergm bootstrap MPLE) | 5 | `_slow/N01.R` | `_slow/N01.csv` |
| N02 | Multilayer ERGM (2019 stack) | 2 | `_slow/N02.R` | `_slow/N02.csv` |
| N03 | Signed TERGM (ergm.sign) | 4 | `_slow/N03.R` | `_slow/N03.csv` |
| N04 | AME (amen 1.4 pa2018) | 2 | `_slow/N04.R` | `_slow/N04.csv` |
| N05 | AME (amen 1.4 pa2018) | 4 | `_slow/N05.R` | `_slow/N05.csv` |
| N07 | Longitudinal SRM probit | 4 | `_slow/N07.R` | `_slow/N07.csv` |
| N08 | Latent space (gbme) | 3 | `_slow/N08.R` | `_slow/N08.csv` |
| 101 | Bayesian split-population probit, custom MCMC | 1 | `_slow/101.R` | `results_estaudit/101__smigoabi.csv` |
| P9_035 | Heckman probit rebuild | 1 | `_slow/P9_035.R` | `_slow/P9_035.csv` |
| P9_045 | Heckman probit port (heckprob, cluster-robust) | 1 | `_slow/P9_045.R` | `_slow/P9_045.csv` |
| P9_057 | Heckman probit port (heckprob, OIM) | 1 | `_slow/P9_057.R` | `_slow/P9_057.csv` |
| P9_092 | Random-effects logit port (glmer, nAGQ 12) | 1 | `_slow/P9_092.R` | `_slow/P9_092.csv` |

## How a rerun works

One entry script per id, `legacy/_slow/<id>.R`, runs the recorded chain (estimation, then the parse or aggregate step) and
writes `SWEEP_OUT/_slow/<id>.csv` with columns `id,key_term,b0,p0,b_dm,p_dm,b_26,p_26,pct_bC,pct_bD,n_used`, one row per
master row (101 writes `results_estaudit/101__smigoabi.csv` instead). `extract_map.csv` points each master row at that file
(`row_filter` on `key_term`).

```sh
cd <PIPE>
export SWEEP_CODE=$PWD/code/04_replication_sweep/legacy SWEEP_DATA=$PWD/data/raw/replication_sweep \
       SWEEP_OUT=$PWD/results/replication_sweep/_build \
       SWEEP_HIST=$PWD/data/processed/replication_sweep/hist/history_tables.rds
export SWEEP_LIB_AMEN=<lib with amen 1.4 pa2018>   # MINHAS2022_* only
export SWEEP_LIB_CHEN=<lib with ergm 3.10.1 stack>  # N02 only
Rscript code/04_replication_sweep/legacy/_slow/<id>.R
```

- The entry scripts call `Rscript` for each step themselves and take no lock. Whoever launches `<id>.R` should hold the R
  lock (on this machine `flock <SWEEP>/R.lock`) for the whole chain. Nesting a flock inside would deadlock.
- `SLOW_SKIP_EST=1` skips every estimation step and re-runs only parse/aggregate on fits already under `SWEEP_OUT`.
  This mode was used to test the 11 network-model entry scripts (all except 101 and P9_035) against copies of the old fits/logs. All 31 rows reproduce `b0` and `b_hist`
  of `expected_est.rds` (tolerance 1e-6), and every other cached column matches too, with one exception. For X060 row 243,
  the hand-typed master `p0` is 0.00433 and the rule-based extraction gives 0.00434 (p = 0.0043365).
- `_slow/_common.R` holds the shared step runner (`run_r`, which fails loudly on a non-zero exit) and the final writer. The
  writer asserts exactly the registry's master rows for the id and finite `b0` and `b_dm`.
- The ported shell drivers (`work_saom/driver*.sh`, `work_net/queue*.sh`, `q2.sh`, `*/run_*.sh`) are kept for provenance.
  They use the same env vars and have flock removed. The pipeline does not call them.
- Python parsers were ported to R (`work_net/parse_cd.R`, `work_net/parse_fritz.R`), so there is no Python dependency.
  The long free-text columns (study, notes, convergence prose) were not ported. All numeric columns match the originals.
- `work_net/agg_mine.R` gained optional tag arguments so that each study can be aggregated on its own. The numbers are unchanged.
- Data construction scripts are not part of the rerun, because their outputs are staged as raw inputs. These are
  `work_net/CD2011/build_cd2011.R` (makes `cd2011_nets.rds`) and `work_net/WH2011/build_wh.R` + `final_wh.R` (make `wh_*.csv`).
- Newly staged input: `work_saom/X057/data/comm_detection/detected_comms/weak/*.rds` (6 files, hard links).
- Disk: X060 writes two full `sienaFit` objects per model (>100 MB each) to `_build/work_saom/X060/`.

## Environment and package versions (all runs: R 4.5.3, x86_64 conda build, 2 CPUs)

Main R library, pin these:
RSiena 1.6.6, btergm 1.11.1, ergm 4.12.0, tergm 4.2.2, network 1.20.0, sna 2.8, statnet.common 4.13.0, peacesciencer 1.2.0
(X055 uses its bundled `atop_alliance` data), igraph 2.3.3, dplyr 1.2.1, data.table 1.18.6.1, countrycode 1.9.0,
reshape2 1.4.5, ggplot2 4.0.3, plyr 1.8.9, fastglm 0.1.2, coda 0.19-4.1, MASS 7.3-66, mvtnorm 1.4-2, Matrix 1.7-5,
Rcpp 1.1.2, RcppArmadillo 15.6.0-1, stringr 1.6.0, lpSolveAPI 5.5.2.0-17.15, doParallel 1.0.17, ergm.sign 1.0.
ergm.sign 1.0 is not on CRAN. It was installed from the tarball in the authors' replication archive
(`SWEEP/work_net/FRITZ2025/src/ergm.sign_1.0.tar.gz`, packaged 2024-03-25). A `Makevars.local` with compiler overrides sits in `SWEEP/work_net/FRITZ2025/`.

`SWEEP_LIB_AMEN` (a separate library directory) holds one package:
- amen 1.4. It was built on 2026-09-29 from GitHub `s7minhas/amen`, branch `pa2018_version`, commit
  `76381be87407cd6aafd230f330d6c1231707ec25` (2017-09-30; source clone in `SWEEP/work_net/MINHAS2022/amen_pa2018`).
  It is not the CRAN amen, which is 1.4.5 in the main library. `run_minhas.R` asserts version 1.4 and the presence of
  `rUV_rep_fc_cpp`. Its imports (Rcpp, RcppArmadillo, reshape2, ggplot2, plyr) come from the main library.

`SWEEP_LIB_CHEN` (a separate library directory) holds the 2019 stack of the N02 author. All four packages were built
on 2026-09-30:
- ergm 3.10.1 (CRAN 2019-05-14)
- network 1.15 (CRAN 2019-04-02)
- statnet.common 4.3.0 (CRAN 2019-06-02)
- multilayer.ergm 0.1.6 (dated 2019-04-07; its DESCRIPTION has no Repository field, so it was installed from source)

`chen_run.R` asserts these versions. The remaining dependencies come from the main library (robustbase 0.99-7,
trust 0.1-9, lpSolve 5.6.23, purrr 1.2.2, rlang 1.3.0, tibble 3.3.1, magrittr 2.0.5, coda, Matrix, MASS). The current stack
(ergm 4.12 + multilayer.ergm 0.2.3) miscomputes the cross-layer statistics, so it must not be used.

## Per study

Runtimes are wall-clock from the old logs.

### X060

SAOM, RSiena. Ported files: `work_saom/X060/rerun_X060.R`, with the helpers `prep_data.R` and `prep_effects.R` in the same
folder, and `hist/add_history.R`. Originally run by `work_saom/driver.sh` and `driver2.sh`.

- **Runtime:** m0 5.28 h, mH 5.94 h (~11.2 h on 2 cores).
- **Convergence settings:** these are the authors' settings, `sienaAlgorithmCreate(nsub=4, n3=5000, firstg=.03, modelType=1 for both networks, useStdInits=F)`.
  Other settings are 2 cluster nodes, `thetaBound=1000`, and one `siena07` call. A continuation run with `prevAns` happens
  only if max|t| > 0.1 or the overall ratio > 0.25; it was not triggered. Logged results: m0 max|t| 0.052 and ratio 0.221,
  mH 0.047 and 0.208.
- **Seeds:** seed 123 in the algorithm. mH is warm-started from `res_m0.rds`.
- **Rerun:** `Rscript legacy/_slow/X060.R` runs m0, then mH.
- **Notes:** m26 and mP were not run, which is why the master has no `b_26` or placebo. The master's values were typed from
  `res_*.rds` (b to 4 dp, p to 5 dp).

### X085

SAOM, RSiena. Ported files: `work_saom/X085/rerun_X085.R`, with the helpers `prep_warren.R` and `RSiena_helper_functions.R`,
and `hist/add_history.R`. Originally run by `driver3.sh` and `driver4.sh`.

- **Runtime:** m0 part 1 1.65 h, part 2 0.93 h, mH 2.07 h.
- **Convergence settings:** the authors' est_siena ladder, re-implemented. n2start 20 -> 800, n3 500 -> 5000, diagonalize 0.8,
  firstg 0.1. It stops at tconv < .1 and tcmax < .25, with at most 100 reps, modelType 3, truncation 10, nsub 1, and 2 cores.
  A rep is discarded if its tcmax exceeds 1.5 x the last one or if it exceeds `REP_MIN` minutes. The whole ladder stops at `CAP_H` hours.
- **Seeds:** the seed is 123, 124, ... per ladder rep, and the counter restarts in every process.
- **Rerun:** `Rscript legacy/_slow/X085.R`.
- **Two-part m0:** part 1 was run with `CAP_H=4 REP_MIN=40` and stopped by hand after rep 6 (n3 switch). Part 2 was run with
  `RESUME=1 N2SW=1 N3SW=1 CAP_H=3 REP_MIN=90` and converged in one rep. mH was then run with `CAP_H=4.5 REP_MIN=90` (warm start
  from m0) and converged at rep 3. The entry script reproduces the manual stop with `X085_PART1_CAP_H`. The default is 1.5 h
  (rep 6 started at 1.29 h). The script also asserts the state that part 1 ended in.
- **Machine dependence:** the timing caps depend on machine speed, so bit-for-bit reproduction is not guaranteed.
- **Notes:** the master's values were typed from `res_*.rds` (b to 4 dp, p to 5 dp).

### X055

TERGM, btergm MPLE with bootstrap. Ported file: `work_saom/X055/rebuild_X055.R`. It was run by hand with no driver.

- **Runtime:** m0 19 min, mH 21 min, m26 2.3 min, mP 2.5 min (~45 min).
- **Convergence settings:** R=1000 bootstrap for m0 and mH, R=100 for m26 and mP. `parallel="multicore", ncpus=2`.
- **Seeds:** `set.seed(123)` before each model.
- **Rerun:** `Rscript legacy/_slow/X055.R` (calls `rebuild_X055.R 1000`).
- **Notes:** m26 is estimated but not in the master. The master's values were typed with b to 4 dp, p_boot as is, and
  pct_bC = 100(\|b0\|-\|bP\|)/\|b0\| to 1 dp.

### X057

TERGM, btergm MPLE with bootstrap. Ported file: `work_saom/X057/rerun_X057.R`. It was run by hand with no driver.

- **Runtime:** 24 fits, ~74 min.
- **Convergence settings:** R=500 (the value passed on the command line) for m0 and mH, R=100 for m26 and mP. There are 6
  embedding files, and the results are averaged over them.
- **Seeds:** `set.seed(1912)` per fit.
- **Rerun:** `Rscript legacy/_slow/X057.R` (calls `rerun_X057.R 500`).
- **Notes:** the run with the default R=2000 was aborted. The master uses the mean over the 6 fits (b to 5 dp, p to 4 dp).

### N03

Signed TERGM, ergm.sign 1.0. Ported files: `work_net/FRITZ2025/fritz_run.R`, `fritz_mple.R` and `run_q.sh`, plus
`work_net/q2.sh` and `work_net/parse_fritz.R` (port of the `.py` parser).

- **Runtime:** m0 42.9 min, hist 46.1 min (2 PSOCK workers).
- **Convergence settings:** tsergm Stepping, tol 0.2, max_it 30, Fisher variance. Variance sampler: 3000 networks, burn-in
  1e5, 1e4 proposals. Estimation sampler: 2000 networks, burn-in 1e4, 1e3 proposals. Stepping grid 1000. m0 took 25 iterations
  and hist took 28.
- **Seeds:** `set.seed(13)`, and sampler seeds 13.
- **Rerun:** `Rscript legacy/_slow/N03.R`.
- **Notes:** m26 and placebo were run as MPLE only, so `b_26` is blank.

### N05

AME, amen 1.4 pa2018. Ported files: `work_net/MINHAS2022/run_minhas.R`, the helpers `x/replArchive/helpers/ame_repL.R` and
`rZBinfc.R`, `work_net/agg_mine.R`, and `queue.sh`, `queue6.sh` and `queue7.sh`.

- **Runtime:** 29.5, 30.3, 32.5 and 32.8 min for m0, hist, m26 and placebo.
- **Convergence settings:** burn 500, nscan 2500, odens 10 (250 draws), R=2, bin. The run warm-starts from the authors'
  `ameFit$startVals` and uses the authors' `rZ_bin_fc` swapped in.
- **Seeds:** 6886.
- **Rerun:** `Rscript legacy/_slow/N05.R` (needs `SWEEP_LIB_AMEN`).

### N04

AME, amen 1.4 pa2018. Ported files: the same as N05.

- **Runtime:** 32.1, 33.0, 31.6 and 31.3 min for m0, hist, m26 and placebo.
- **Convergence settings:** burn 500, nscan 3000, odens 10 (300 draws), R=2. The run warm-starts from the authors' state.
- **Seeds:** 6886.
- **Rerun:** `Rscript legacy/_slow/N04.R` (needs `SWEEP_LIB_AMEN`).

### N08

Latent space, gbme. Ported files: `work_net/WSC2007/run_wsc.R`, `gbme_ws.r`, `work_net/agg_mine.R` and `queue7.sh`.

- **Runtime:** 31.4, 38.7, 33.9 and 33.9 min for m0, hist, m26 and placebo.
- **Convergence settings:** gbme binomial, k=3, 25,000 scans with odens 20. The first 5,000 scans are dropped in agg (1000 draws).
  The run warm-starts from the last draw of the authors' 497k-scan chain.
- **Seeds:** `set.seed(123)` and gbme seed 123.
- **Rerun:** `Rscript legacy/_slow/N08.R`.
- **Notes:** the cold-start runs from queue.sh and queue6.sh were discarded.

### N07

Longitudinal SRM probit. Ported files: `work_net/WH2011/probit_wh.R`, `work_net/agg_mine.R`, `queue.sh` and `queue7.sh`.

- **Runtime:** 39.4, 35.7, 42.3 and 44.0 min for m0, hist, m26 and placebo.
- **Convergence settings:** 12,000 iterations after a burn-in of 3,000, thin 12 (1000 draws). The authors' chain was 7e6
  iterations. The priors are the authors'.
- **Seeds:** `set.seed(123)`.
- **Rerun:** `Rscript legacy/_slow/N07.R`.
- **Notes:** the .ok markers for m0 and hist were originally touched by hand.

### N01

TERGM, btergm MPLE with bootstrap. Ported files: `work_net/CD2011/est_cd2011.R`, `run_all.sh` and `work_net/parse_cd.R`
(port of the `.py` parser).

- **Runtime:** 7 models x ~10-11 min (73.6 min).
- **Convergence settings:** R=1000 network-year bootstrap, fastglm.
- **Seeds:** `set.seed(123)`.
- **Rerun:** `Rscript legacy/_slow/N01.R`.
- **kstar2 and triangle rows:** these rows have no results CSV. They are read from the parsed log (`log_cd2011_parsed.rds`),
  which reproduces the master's hand-made conventions: a bootstrap p of 0 is reported as 0.0005, and
  pct_bC = 100(\|placebo\|-\|pbase\|)/\|pbase\| unrounded.

### N02

Multilayer ERGM on the author's 2019 stack. Ported files: `work_net/CHEN2021/chen_run.R`, `write_chen.R`, `run_mcmle.sh`,
`run_mcmle2.sh` and `run_diag.sh`.

- **Runtime:** MCMLE m0 33 min (4 iterations) and hist 40 min (5 iterations), 63 and 70 min wall. DIAG 9.4 min each.
  MPLE takes under 1 min each.
- **Convergence settings:** MCMLE with burn-in 1e5, interval 2000, samplesize 4000, maxit 8, 2 parallel chains, and ergm 3.10.1
  defaults (Hummel). The authors used interval 1e4, samplesize 2e4 and 4 chains. The post-hoc DIAG uses 4000 draws at interval 2000.
- **Seeds:** `set.seed(20842)` and control seed 30748.
- **Rerun:** `Rscript legacy/_slow/N02.R` (needs `SWEEP_LIB_CHEN`).
- **Notes:** m26 and placebo were run as MPLE only, so `b_26` is blank.

### 101

Bayesian split-population (Boolean) probit with a custom MCMC sampler (`xiang_mcmc.cpp`). Entry script `_slow/101.R` runs
the ported shell driver `work/101/est/run_chains.sh`.

- **Runtime:** about 7.4 h logged in the sweep (16 segments x 25,000 iterations per model, m0 and mDM interleaved; burn-in
  100,000 in the write step).
- **Rerun:** `Rscript legacy/_slow/101.R`; it writes `results_estaudit/101__smigoabi.csv` under `SWEEP_OUT`.

### P9_035

Rebuilt censored (Heckman) probit on directed COW dyad-years 1816-2000. Ported file: `work_p9/P9_035/rebuild_P9_035.R`; slow wrapper `_slow/P9_035.R`.

- **Runtime:** about 16 min for the six-step run (heckprob ML via `work_p8/heckprob_b3.R`, OIM SE).
- **Cache:** `results/replication_sweep/cache_slow/est_P9_035.rds`; `run_one.R P9_035 --cached` checks it.
- **Rerun:** `Rscript legacy/_slow/P9_035.R`.

### P9_045

Heckman probit (heckprob) port of the authors' Table 1 Model 1 on 907k directed dyad-years (selection = MID initiation, outcome =
target reciprocation). Ported file: `work_p9/P9_045/rerun_P9_045.R` (heckprob ML via `work_p8/heckprob_b3.R`, cluster-robust SE on
`ddyadid`); slow wrapper `_slow/P9_045.R`. The data have no year column; the script rebuilds it from row order and checks it against
the dated MID file (see the script header).

- **Runtime:** about 46 min on 1 BLAS thread (7 heckprob fits, about 6.5 min each).
- **Cache:** `results/replication_sweep/cache_slow/est_P9_045.rds`; `run_one.R P9_045 --cached` checks it.
- **Rerun:** `Rscript legacy/_slow/P9_045.R`.

### P9_057

Heckman probit (heckprob) port of the authors' Table 2 on 585k directed challenger-target dyad-years, OIM SE. Ported file:
`work_p9/P9_057/rerun_P9_057.R`; slow wrapper `_slow/P9_057.R`.

- **Runtime:** about 11 min under `run_one.R` on 1 BLAS thread (7 heckprob fits; 2026-10-10).
- **Cache:** `results/replication_sweep/cache_slow/est_P9_057.rds`; `run_one.R P9_057 --cached` checks it.
- **Rerun:** `Rscript legacy/_slow/P9_057.R`.

### P9_092

Random-effects logit port of Millerd (2023) Table 2 Model 3 (author's Dataverse archive doi:10.7910/DVN/I9ZSPR; `xtlogit, re` by
dyad on 50,216 PR dyad-years 1989-2014). Ported file: `work_p9/P9_092/rerun_P9_092.R` (lme4::glmer, random intercept by `dyad_id`,
nAGQ = 12 as Stata's default adaptive quadrature, `autoscale = TRUE`); slow wrapper `_slow/P9_092.R`. The common sample uses the
named 0.75 guard exception in `work_p8/common_b3.R` (PI decision 2026-10-11; see PORTING_CONVENTIONS.md).

- **Runtime:** about 20 min per glmer fit, about 2.5 h for the seven fits on 1 BLAS thread (2026-10-11).
- **Cache:** `results/replication_sweep/cache_slow/est_P9_092.rds`; `run_one.R P9_092 --cached` checks it.
- **Rerun:** `Rscript legacy/_slow/P9_092.R`.
