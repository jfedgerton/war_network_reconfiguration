# Porting conventions (internal build notes for code/04_replication_sweep)

Goal: move each study's FINAL estimating script (per `lineage_map.csv`) from the exploratory sweep into the
manuscript pipeline without changing any estimation logic, so the pipeline reproduces every number in the
current master table exactly.

## Four paths (set by `run_one.R` as environment variables; never hard-code absolute paths)
| Variable | Meaning | Pipeline location |
|---|---|---|
| `SWEEP_CODE` | ported R code | `code/04_replication_sweep/legacy/` |
| `SWEEP_DATA` | raw author data + history tables (read-only) | `data/raw/replication_sweep/` |
| `SWEEP_OUT`  | everything the scripts write (logs, fits, intermediate/result CSVs) | `results/replication_sweep/_build/` (deleted by `run_one.R` after extraction, also on error) |
| `SWEEP_HIST` | conflict-history lookup tables read by `hist/add_history.R` | `data/processed/replication_sweep/hist/history_tables.rds` (built by `00_build_history.R`) |

Every ported file keeps its ORIGINAL relative path under each root (e.g. `work/64/rerun_key_64.R` lives at
`legacy/work/64/rerun_key_64.R`; its data at `data/raw/replication_sweep/work/64/...`; its outputs at
`_build/work/64/...` and `_build/results_key/...`).

## Mechanical path rules for ported code
1. Any hard-coded `ROOT <- "<absolute path to the original sweep folder>"` becomes three lines:
   `ROOT <- Sys.getenv("SWEEP_DATA"); OUT <- Sys.getenv("SWEEP_OUT"); CODE <- Sys.getenv("SWEEP_CODE")`
   followed by `stopifnot(nzchar(ROOT), nzchar(OUT), nzchar(CODE))`.
2. `source(<abs or ROOT>/X.R)` -> `source(file.path(CODE, "X.R"))`.
3. READS of raw inputs (author data, `hist/history_tables.rds`, zips) -> `file.path(ROOT, ...)`.
4. WRITES (sink, write.csv, fwrite, saveRDS, save, ggsave, pdf/png, dir.create for outputs) -> `file.path(OUT, ...)`
   with the same relative subpath; create parent dirs with `dir.create(..., recursive = TRUE, showWarnings = FALSE)`.
5. READS of files the pipeline itself creates (fits, slim .rds, intermediate CSVs, prior result CSVs) -> `file.path(OUT, ...)`.
6. Files in the parent scratchpad read via `dirname(ROOT)` (e.g. `jcr2006replication.dta`) -> copy under
   `SWEEP_DATA/_parent/` and read `file.path(ROOT, "_parent", ...)`.
7. Data read from an installed package's extdata (`peacesciencer` cow_trade_ddy.rds, dyadic_fp_similarity.rds)
   -> `system.file("extdata", "<file>", package = "peacesciencer")`, and assert it exists.
8. Zips that were unzipped and then deleted: add an idempotent unzip into `OUT` (not into DATA) before reading,
   e.g. `unzip(file.path(ROOT, "work/64/web-appendix.zip"), files = "<member>", exdir = file.path(OUT, "work/64"))`,
   and read from OUT.
9. Do not change any estimation code, sample restriction, formula, seed, or SE computation.

## Names must not reach disk or logs with results
- Do not print study names/citations alongside estimates. Replace `cat()`/`print()` lines that print a citation
  or `meta$short_cite` with the study id. Header comments in code may keep citations (code is not output).
- Result CSVs written by legacy writers may still contain `short_cite`; they live only in `_build/` and are
  deleted after `run_one.R` extracts ID-keyed numbers.
- One deliberate exception: `91_figures.R` labels the five pre-chosen network-measure studies (`named_illustration == "Y"`)
  on `results/main/fig_contamination_ties.*`; the PI chose on 2026-10-10 to keep these five names on the figure.

## Extraction map
For every master row, add a line to `extract_map.csv` saying where the ported script's output lands:
`id, key_term, out_file (relative to SWEEP_OUT), row_filter (R expression over that file's columns, or blank),
 col_b0, col_p0, col_b_hist, col_p_hist, col_b_26, col_p_26, col_placebo_pct, col_jointdem_pct, col_n`
Use the file that the ESTIMATING script writes (if the master value was later copied by a Python chain
script, point to the estimating script's own output instead). Blank a column if that quantity is not produced.

## Verification
After porting, run the script through `run_one.R <id>` and check the extracted values against
`expected_est.rds` (ID-keyed master values) with relative tolerance 1e-6 on b0 and b_hist.

## Phase work_p9 additions (registered 2026-10-10)
- **Slow studies.** A study whose `run_one.R` run takes more than about 10 minutes gets `slow = Y`, an entry script
  `legacy/_slow/<id>.R` (pattern `_slow/P9_035.R`: run the estimating script with `run_r()`, then `write_slow()`), an
  `extract_map.csv` row pointing at `_slow/<id>.csv` with a `key_term` row filter, and a cached
  `results/replication_sweep/cache_slow/est_<id>.rds` copied from a passing uncached run (P9_035, P9_045, P9_057).
- **Missing year column.** If the author's file has no year, rebuild it only for attaching history, validate it against an
  independent dated source, and fail loudly below a stated agreement (P9_045: row-order blocks within target, anchored on the
  dated MID file, at least 97% agreement required).
- **Hand-downloaded packages.** Archives the PI downloaded by hand (Cloudflare or publisher 403) go into `data_sources.csv` with
  `source_type = unknown`, the landing page as `source` and the original member name in `archive_member`.
- **Sensitivity hooks.** A headline product term, or a constituent of one, uses `sens_split()` (scope `split`); a single
  additive term uses `sens_run()` (scope `full`).
- **alpha_pub.** Every new headline row gets `alpha_pub` with evidence in `alpha_pub_coding.csv`. Use 0.10 only when the paper
  uses one-tailed tests or the selected headline claim itself rests on p < .10; otherwise 0.05. A .10 star in the legend alone, or
  .10 language about other terms, does not qualify (P9_041, P9_053, P9_056 recoded to 0.05 on 2026-10-11; reason in
  `decision_note`). Mark a value coded from a star legend alone as low-confidence in `methods_note`.
- **Author history control.** `author_history` records any prior-conflict control (MID counts or lags, rivalry, ongoing or
  other MID involvement, third-party MID involvement) in ANY stage or equation of the author's model (A_dyadic / A_monadic);
  peace years and splines alone are B_duration. 2026-10-11 recodes under this rule: P9_089 (thirdL/thirdH, first stage),
  P9_094 (other MID involvement, second stage; MID involvement, first stage), 55 (sharerival, alliance equation).

## Round-2 rebuild rule (PI, 2026-10-11)
- A protocol-B rebuild (no usable author analysis data) MATCHES only if all three hold: (1) the rebuilt baseline b0 has the
  published sign; (2) b0 is significant at the paper's `alpha_pub` whenever the published estimate is significant at
  `alpha_pub`; (3) b0 lies within 0.5x-2x of the published coefficient. Otherwise `REBUILT_MISMATCH`, `in_headline = N`, and a
  `methods_note` saying why (e.g. "magnitude outside 0.5-2x of published").
- One implementation: `legacy/work_p8/round2_rule.R` (`rebuild_verdict(id, term, b0, p0, pub_b, pub_sig, pub_se)`), which reads
  `alpha_pub` from `registry.csv` (unregistered rebuilds fall back to 0.05 with a printed note). Every helper that writes a
  rebuild verdict (`work_p8/common_b3.R`, `common_b2.R`, `p8_common.R`, `vb3/vhelp.R`, `vb4/vhelp.R`, `vb5/verify_helpers.R`,
  `b1work/v/vhelp.R`) and every registered rebuild script with an inline verdict calls it; the old `.05` rebuild branches stop.
- Ports of author data keep the REPRODUCED / CLOSE / NOT_REPRODUCED rule (no magnitude gate).
- Verdict strings are informational: `run_one.R` checks numbers only, and pooling follows `registry.csv` `in_headline`.

## Common-sample guard exceptions (2026-10-11)
- `common_sample()` in `legacy/work_p8/common_b3.R` requires the common sample (1822-2011, complete on history terms) to keep more
  than 80% of the complete-case rows. The global threshold `CS_MIN_SHARE` stays 0.80.
- Exceptions are named per study in `CS_MIN_SHARE_EXCEPTIONS` and apply only when the study script passes its own `id`.
  Current list, P9_092 (Millerd 2023) at 0.75. PI decision of 2026-10-11. Its sample runs 1989-2014 and the 2012-2014 rows lie
  outside the history window (the GML MIDs behind `hist/` end in 2010, so history is defined through 2011), leaving 79.1% of rows.
  The exception is also recorded in the study's `methods_note`. Any new exception needs a PI decision and an entry here.
