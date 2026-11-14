# Tremors in the System — replication sweep and submission decision log (Nov 13–14, 2026)

This log records every decision made during the Nov 13–14 audits, so they do not need to be revisited. The per-study registry state as of the P3 fixes is appended at the end. Code: Dropbox ottoman_byzantine/code/04_replication_sweep. Backups: tmp/backup_add20261113/ and tmp/backup_restructure_20261114/.

## 1. State after round P3 (pending the P3 rerun `code/run_P3_20261114.sh`)

**Pool.** The main pool is 90 studies. The P3 rerun will finalize the exact family split, expected to be about 59 tie / 24 democratic peace / 7 capitalist peace.

**Appendix expanded set** (main + outcome_overlap, Table 10). It holds 32 overlap studies, about 122 studies in total.

**Last published numbers** (92-study pool, in Overleaf before the P3 rerun):

| | Main text | Expanded set |
|---|---|---|
| Tie | 13 of 61 change (21%); median −23.0% | 17 of 77 (22%); median −23.0% |
| Democratic peace | 3 of 24 (13%) | 5 of 35 (14%) |
| Capitalist peace | 1 of 7 (14%) | 3 of 10 (30%) |

- Test statistics: 74% move closer to zero.
- Tie authors with no history control: 47 of 61.
- Broader rule: 13 of 71 (18%), median −17.8%.
- All attempts: tie 19 of 85, DP 7 of 39, CP 3 of 12; 164 plotted.
- Sensitivity set: 65 studies (46 additive, 19 interaction; 45 tie / 14 DP / 6 CP).
  - Coefficient ratios: dyadic 0.87, monadic 0.95, both 0.81.
  - 12 of 65 lose significance.
  - Partial R² of history 3.1%, larger than either benchmark in every study.
  - Placebo within 5% in 93% of 46 studies.
- Dial (8 studies): 1.01, 0.95, 0.91, 0.86, 0.86.
- Word count 9,990.
- SI body ends on printed p. 26, so the SI-exception text should say 26.

## 2. Pooling rules (cumulative)

1. **Groups by the article's own result.** There are three:
   - The author claims an effect and the term is significant at the author's threshold. These enter the headline count.
   - The author argues a null. These are reported separately.
   - The term is not significant and the author does not argue a null. These are excluded.
2. **Changed.** A study counts as changed if it is significant at baseline and either becomes not significant at the article's own alpha or stays significant with a flipped sign. One headline term per study.
3. **Rebuilds.** They must land within 0.5–2× the published coefficient with the same sign and significance (the "round-2 rule").
4. **Studies that do not reproduce are not pooled.** This applies even when the key term matches, if the published specification itself (its controls and peace years) does not.
5. **Common sample.** It keeps rows with conflict history AND the lagged CINC/Polity benchmark controls, so every compared model runs on the same rows. This is by design.
6. **History.** Counts of MID onsets in t−1…t−5, dyadic and monadic, anchored on the year the outcome is measured. If the outcome is a t+1 lead, the anchor is year+1; if it is same-year, the anchor is year.
7. **Outcome-overlap rule (Nov 14, PI decision: most generous main text).** A study goes to pool_status outcome_overlap, listed only in appendix Table 10 with an expanded-set footnote, if either:
   - its outcome codes every year of a dispute; or
   - its outcome is conditional on an ongoing dispute or spell.

   **How studies were classified.** For each outcome I computed two shares of outcome-years against COW dyadic MID 4.02: the share that are COW onset years, and the share that are COW MID years. A study is classified as involvement when the onset share is well below the MID-year share (ratio below about 0.8) and the continuation share is high (about 0.4 or more). The script is tmp/datachk/chk.py on the Mac.
8. **Main-text figures name no studies.** The PI decided no one is called out. The named-illustration check in 91_figures.R was removed.

## 3. Studies moved to outcome_overlap (32)

| Group | Studies | Basis |
|---|---|---|
| TRIANGLE `dispute1` (involvement and t+1 lead) | P9_028, P9_029, P9_089, P9_061, X040, X092, 101 | dataset check |
| `maoznewl` (946 of 1,750 MID years are continuation years) | C1, C2, D2 | data check |
| Other COW-checked involvement outcomes | 62 `force_y1`, 118 `disp_l1`, X093 `mid`, P9_092 `mid`, 55 `conflict`, X037 `fatal_mid`, 115 `fatal_mid` (his JOP 2024), X031 `esc_s1`, P9_043 `ffatalmid`, P8_036 `MID`, inv34 `weakini`, P9_050 `midsidea` | COW onset-match vs MID-year-match shares, listed below |
| Involvement by construction or definition | 119 (MID under way), X059 (ongoing MID), X033 (gmlmidongoing), P8_048 (ongoing or onset), P8_042 (participation), P9_085 (involvement), X085 (15-year persistent conflict ties), P9_053 (involvement including joiners) | variable definition |
| Conditional on an ongoing dispute or spell | P8_030 (escalation within a confrontation), 100 (Cox spells that start with a MID) | design |

COW shares (onset match / MID-year match) for the checked involvement outcomes:

| Study | Variable | Onset match | MID-year match |
|---|---|---|---|
| 62 | `force_y1` (lead) | 0.57 | 0.93 |
| 118 | `disp_l1` (lead) | 0.56 | 0.69 |
| X093 | `mid` | 0.64 | 0.83 |
| P9_092 | `mid` | 0.78 | 1.00 |
| 55 | `conflict` | 0.52 | 0.93 |
| X037 | `fatal_mid` | 0.55 | 0.96 |
| 115 | `fatal_mid` | 0.45 | 0.84 |
| X031 | `esc_s1` | 0.73 | 1.00 |
| P9_043 | `ffatalmid` (lead) | 0.51 | 0.93 |
| P8_036 | `MID` | 0.50 | 0.72 |
| inv34 | `weakini` | 0.71 | 1.00 |
| P9_050 | `midsidea` | 0.75 | 1.00 |

## 4. Studies removed from the pool (not overlap)

- **Leeds & Mattes 2007 (ADD_T03) → not_reproduced.** The geeglm AR(1) port does not reproduce the printed table: key term −0.303 vs −0.268, NAP −0.466 vs −0.178. The gate had been hard-coded "REPRODUCED".
- **McDonald 2015 (P9_087) → rebuild_outside_rule.** The full rebuilt sample is 2.0× the published value. The GPDEM main term is 2.20 vs a published 0.897, so the variable is mis-built.
- **Powers 2006 (57) → not_reproduced.** With the printed control "potential target in REI alliance" (`lagptreiall`) added, the coefficient is −0.963 (p = .08) vs a printed −1.184 (.441). The script verdict is NOT_REPRODUCED.
- **Park & Colaresi 2014 (X020) → not_reproduced.** The key term matches (−0.0297 vs −0.030, same N), but the peace-year and spline terms and the controls do not match published Model 4 (peace years −0.468 vs −0.350; spline about 15×).
- **Nov 13 drops:**
  - The RULES list in drop_20261113.py: X060, X079, N03, X077, 50, P9_078, P8_037, C6, X025, P9_057, P9_045, X022 and Benson 52 (after its switch to DyadTotalDissim, dropped by PI decision).
  - Earlier Nov 13: Peterson 2014, Westveld & Hoff 2011, Barnhart et al. 2020.

## 5. Code fixes (rerun and rebaselined)

**Nov 13: history anchored on the outcome year for t+1 leads.**
- Studies: D2, P8_003, P9_028, P9_029, P9_061, P9_089, X040, X092, 12, 63, 11.
- The TRIANGLE rebuilds whose article defines the outcome in year t (P9_028, P9_089, P9_061) use the prior row's `dispute1` as `dispute_t`, anchored on year.
- The 63 assert was changed to 204,921.

**Nov 14 fixes:**
- **P8_031 (Gartzke & Weisiger).** `deadlyl` is a t+1 lead: COW MID-year match 0.97 at t+1 vs 0.59 at t. The peace-year reset test cannot detect a lead, so `lead_dv` is set to TRUE.
- **P9_029.** The affinity double lag was removed. `dispute1` is t+1, so affinity in the row year is the one-period lag. The baseline is now 0.976.
- **P8_048.** The printed Developed dummy was added (the study is now in overlap).
- **X040.** Uses an independence working correlation, because the AR(1) fits with history diverged (alpha about 0.92, coefficients about 1e8). Now b0 −1.301, b_hist −1.101 (in overlap).
- **P9_080 (P3).** A missing ally power change is now NA instead of 0 in `allypc`; 5,362 of 82,547 links have no CINC. Rerun pending.
- **Kinne 2012 (11).** The authors' `sample==1` cannot be applied because the variable is absent from the data. PI: leave as is.

**Script-editing lesson.** An inline `#` comment placed mid-line commented out the rest of an R vector or pipe (Powers base_ctrl, P9_029). Put comments on their own line.

## 6. Settled facts (do not re-audit)

**Authoritative results.**
- The est_<id>.rds files are the authoritative pooled numbers.
- The *_check.txt files are not rewritten after rebaselining, so agents' "stale pooled numbers" claims were all false. Verified by reading the rds files for D2, 11, 12, 63, P8_003, P9_028, P9_029, P9_061, P9_089, X040, X092 and 103.

**Outcome datasets verified against COW** (onset-coded at the anchored year):
- **Same-year:** X043, X001 (`sideaa`), X003, X087, X067, P9_040, X034, P8_046, P9_048, 10, 61, 1, P8_050, 91, P9_056, 122, P8_008, ADD_D11, 106, P8_027, 64, 88, P8_035, P8_015, X068, 72, 42, P8_020, ADD_T05, 68, X032, 54, 17, 109, P8_019, ADD_D10, X061, 14, 15, 65, X080, 103, X052, X078, X009 (fatal onset), 32.
- **Leads, anchored at year+1:** ADD_D01, 71, P8_025, X035, P9_047, P8_003, 63, 12, 86, 66, X014, P9_064, P8_031.

**Kenwick, Vasquez & Powers 2015 (13).** `year` is the alliance-formation year (the treated USA–Canada row is 1949, NATO). History t−1…t−5 is the pre-treatment window, so the anchor is correct.

## 7. Accepted issues — documented, not to be re-raised

**Specific studies:**
- **Malis 2021 (25).** The history window overlaps post-instrument years in the IV design; the coefficient moves −2.6%. Disclose only.
- **Barbieri & Peters (P9_071).** A weak rebuild (openness wrong sign, N 2×) that passes the gate on its headline term.
- **Sparse or fragile fits:**
  - P8_046 (McDonald 2010): sparse rare-events fit with SE inflation.
  - P8_012 (Choi rebuild): dyad-clustered vs robust SEs, baseline p = .045.
  - X036: random-effects gap, −0.458 vs −0.536.
  - X055 (SAOM): controls depart from the article, and the bootstrap SE ratio uses |b_hist|/|b0|.
  - P9_079: change uses the absolute-value convention.
  - X066: year reconstruction.
- **Rebuild closeness:** P9_093, inv8, P9_080, P9_082, P9_068, P8_003 (partial reconstruction), P8_049 (N +2.7%), P8_025 (N gap, SE 23% lower).
- **Alphas.** 32 uses alpha 0.10. One-tailed alphas would be 0.10 for X014, X001 and 64, with no classification change.
- **Bootstrap reps.** P8_008 uses 100 reps instead of 500 (disclosed).
- **Benchmark models fit on their own subsamples** (55, 61, 109, 68, 73, 69). The history comparison is unaffected.

**Table 9 issues the PI chose to skip:**
- the (t+1) marker is applied unevenly;
- constituent-of-interaction terms are not marked (no § marker);
- unit labels omit sample restrictions (case-control, target-has-ally, PRD, small power gap, exclusive PTA, etc.);
- blank published SEs;
- stale "unverified / PDF not accessible" methods notes (ADD_T05, ADD_T03, ADD_D01, P9_056, ADD_T01, P8_046, 32);
- headline-model choice questions for Gibler & Braithwaite (col 3 vs 2) and Li et al. 2017 (one-degree term);
- author-history recodes for X009, 42 and X069, not applied.

**Supplements not in the folder** (reproduction unverifiable): Renshon 2016 App. C, Morrow 2017 Table A1, Kenwick & Vasquez 2017 appendix, Kenwick et al. 2015 appendix, Dafoe 2011 SI Table 7, P9_062 Table SF-3, Beardsley et al. 2020 appendix.

## 8. Table 9 (appendix study table) changes applied Nov 13–14

- **Reiter & Stam 2003 (X001).** The registry published value was changed to the printed 1.05 (0.140); it had been the replicated 1.04494.
- **Outcome labels:**
  - Dafoe 2011: "MID involvement (t+1)".
  - Goldsmith 2013: "MID onset (t+1) (selection equation)".
  - Gartzke 2007 and Choi 2011: "MID involvement (t+1)".
  - P8_019: "MID onset", per the COW check.
  - Several overlap studies relabeled "involvement".
- **Locations:**
  - Hafner-Burton & Montgomery 2012: "Table 1 Model (3)".
  - Nieman 2016: "Table III (status-quo equation)".
- **Key-variable labels:**
  - McDonald 2015: "Great-power-supported democracy × democracy".
  - Clare 2007: "Redemocratization", estimator Logit.
  - Joshi et al. 2015: "Institutional democracy index".
- **Note.** The † (rebuilt) and ‡ (P9_093 contagion term) note was merged into the longtable.
- **New table.** 95_study_table.R now writes a second table, tableA_sweep_overlap.tex (Table 10).

## 9. Manuscript and submission edits

**Nov 13 fixes:**
- The Table 3 note says the Granger p-values are bootstrap p-values.
- The AI disclosure is an unnumbered note after \maketitle.

**Round R1 manuscript audit — applied:**
- word count;
- SI §3 "list that follows" clause deleted;
- "network-model search" phrase deleted;
- title-page word-count wording;
- abstract and title page set 12-pt double-spaced.

**Round R1 — skipped by the PI:**
- "dyadic" qualifier on IGO amplification;
- fn 2 "residualized" wording;
- SI Table 10 / Fig 8 systemic SE;
- Imber and İnalcık page cites (to check against his own copies).

**Outside (GPT) audit — applied quick fixes:**
- 12-pt double-spaced footnotes (footmisc);
- "doubling of one plus the count";
- FEVD Cholesky-ordering sentence (conflict first);
- the Table 9 labels in §8.

**Outside audit — deferred by the PI:**
- causal framing ("backdoor pathway");
- permutation null and exchangeability;
- the p < .10 "significant" wording vs AJPS .05;
- future-conflict diagnostic;
- placebo interpretation;
- normalization and JS invariance;
- earthquake "resolve" inference;
- novelty framing;
- median-SE footnote;
- "null" wording;
- Ruggiero citation;
- table and figure density.

**Overlap-rule prose (Nov 14).**
- New main-text footnote sentences on the set-aside studies and the expanded set.
- New SI paragraph explaining the overlap set.
- Main-text reference to the IGO table is now appendix Table 11.

**Upload and author details.**
- SI exception text for the AJPS Comments to the Editorial Office (change 25 → 26): "I request an exception to the 25-page SI limit. The appendix documents [N] replication studies and their source references…"
- Gmail on the title page is intended.

## 10. Open decisions

- **N06 (his PSRM piece).** Which version is accepted:
  - Replication files: Table 3 col 3, 0.94 (0.21), IID SEs.
  - Current draft: Table 2 col 6, 1.13 (0.53), SEs three-way clustered by sender, receiver and year.
- **SI Table 10 / Fig 8.** The systemic SE (0.0023) vs main Table 2 Newey-West SE (0.0048). Options: state the SE type, or rerun with Newey-West.
- **Pipeline audit loop.** P1 and P2 found real issues. P3 found 0 MUST items but three pool-changing fixes (P9_080, P8_030, 100). P4 is still needed. Stop after a round with no MUST items and no pool-changing findings.

## 11. Per-study registry state (generated from registry.csv and study_table_labels.csv after the P3 edits)

### pool_status = main (90)

| id | study | family | key term (label) | outcome label | Nov 13–14 notes |
|---|---|---|---|---|---|
| C5 | Mousseau 2013 ISQ | CAPITALISM | CIEl (Contract-intensive economy) | Fatal MID onset ($t+1$) |  |
| P8_031 | Gartzke & Weisiger 2014 | CAPITALISM | pcenerg (Systemic development (world energy per capita)) | Fatal MID onset |  |
| P8_046 | McDonald 2010 | CAPITALISM | lpublicb (Target public property) | War targeting |  |
| P8_049 | Chatagnier & Castelli 2016 | CAPITALISM | ind_hi (Industry share (higher)) | Fatal MID onset |  |
| P8_050 | Mousseau 2019 | CAPITALISM | baxecon (Both contractualist or axial economy) | Fatal MID |  |
| X032 | Boehmer & Sobek 2005 | CAPITALISM | logdev2 (Development squared) | MID initiation |  |
| X043 | Mousseau 2016 | CAPITALISM | CIEl (Contract-intensive economy) | Fatal MID onset |  |
| ADD_D01 | Lektzian & Souva 2009 | DEMOCRACY | sel_jointd6 (Joint democracy (selection equation)) | MID onset ($t+1$) |  |
| ADD_D10 | Zellman et al. 2026 | DEMOCRACY | pop_x_poly (Populism $\times$ polyarchy) | Hostility level | [2026-11-13 audit] verified against the published ISQ 70(4) sqag085 Table 1 Model 2: -1.535* (0.629). |
| ADD_D11 | Joo et al. 2025 | DEMOCRACY | rightparti2 (RWP leader $\times$ participatory democracy) | MID initiation |  |
| P8_017 | Huth & Allee 2002 | DEMOCRACY | accountability (Accountability) | Threat or use of force |  |
| P8_019 | Peceny et al. 2002 | DEMOCRACY | single (Single-party dyad) | MID onset | [2026-11-13 audit] published_b corrected to printed Table 3 M4 value (-0.691); was -0.8407 from the Reiter & Stam log. |
| P8_020 | Reiter & Tillman 2002 | DEMOCRACY | vanpart (Political participation) | MID initiation |  |
| P8_022 | Bueno de Mesquita et al. 2004 | DEMOCRACY | Coalition size x proportion of capabilities (Wratioany) (Coalition size $\times$ capabilities) | War involvement given dispute |  |
| P8_025 | Souva 2004 | DEMOCRACY | polsim (Political similarity) | MID onset ($t+1$) |  |
| P8_027 | Lai & Slater 2006 | DEMOCRACY | military (Military regime) | MID initiation count |  |
| P8_033 | Colgan & Weeks 2015 | DEMOCRACY | revpers10_lag1 (Revolutionary personalist regime) | MID initiation |  |
| P8_035 | Bak, Kenwick & Palmer 2016 | DEMOCRACY | relcap_targ_cgv_mil (Target capability $\times$ military regime) | MID initiation |  |
| P9_062 | Hegre et al. 2020 | DEMOCRACY | cspart_ij (Joint civil-society participation) | Fatal MID |  |
| P9_064 | Anderson & Souva 2010 | DEMOCRACY | w1power (Winning coalition $\times$ power) | MID initiation ($t+1$) |  |
| P9_068 | Henderson 2009 | DEMOCRACY | afr_x_dem (Africa $\times$ lower democracy) | MID onset |  |
| X001 | Reiter & Stam 2003 APSR | DEMOCRACY | pdemdtar (Personalist initiator, democratic target) | MID initiation | [2026-11-13 audit] published_b set to printed 1.05 (.140), Table 1 Model 4; earlier value 1.04494 was the replicated estimate. |
| X003 | Weeks 2012 | DEMOCRACY | strongmanjlw_1 (Strongman regime) | MID initiation |  |
| X009 | Gibler & Braithwaite 2013 BJPS | DEMOCRACY | lowdem (Lower democracy) | Fatal MID onset |  |
| X014 | Oneal & Tir 2006 | DEMOCRACY | demXgrow (Democracy $\times$ growth) | Fatal MID initiation ($t+1$) |  |
| X021 | Crescenzi & Kadera 2015 ISQ | DEMOCRACY | c_regstrengthsq (Democratic community strength (squared)) | Fatal MID onset |  |
| X034 | Clare 2007 | DEMOCRACY | redem1 (Redemocratization) | MID initiation |  |
| X035 | Goldsmith, Chalup & Quinlan 2008 | DEMOCRACY | polsim3 (Political similarity) | Fatal MID onset ($t+1$) |  |
| X036 | Joshi, Maloy & Peterson 2015 | DEMOCRACY | L_IDI_weighted (Institutional democracy index) | MID initiation |  |
| X052 | Daxecker 2007 | DEMOCRACY | abs_democ (Regime change $\times$ democratization) | MID onset | [2026-11-13 audit, PI] pooled headline: H2 interaction (regime change x democratization), printed -0.046** (0.019). |
| 1 | Leeds 2003 | TIE | ptargdef (Target has defensive ally) | MID initiation |  |
| 10 | Benson 2011 | TIE | uncondcompel (Unconditional compellent alliance) | MID initiation |  |
| 103 | Bak 2018 | TIE | cowcontigdum (Contiguous defensive ally (any)) | MID initiation | [2026-11-13] key term switched from relcincally (a control, p. 118) to the H1 headline term cowcontigdum. |
| 106 | Choi 2023 | TIE | lnbet1 (Betweenness centrality (log)) | MID initiation |  |
| 109 | Johnson & Leeds 2011 | TIE | ptargdef (Target has defense pact) | MID initiation |  |
| 11 | Kinne 2012 | TIE | cent50 (Network centrality) | MID initiation count ($t+1$) | [2026-11-13 audit] history anchored on the outcome year (DV is a t+1 lead; checked vs COW dyadic MID). |
| 114 | Beardsley et al. 2020 | TIE | samecomm_hierlow_close (Same community $\times$ low hierarchy) | MID onset ($t+1$) |  |
| 12 | Kleinberg, Robinson & French 2012 | TIE | lower_conc (Trade concentration (lower)) | MID onset | [2026-11-13 audit] history anchored on the outcome year (DV is a t+1 lead; checked vs COW dyadic MID). |
| 122 | Chatagnier & Kavakli 2017 | TIE | L_exports (Export similarity (lagged)) | MID onset |  |
| 13 | Kenwick, Vasquez & Powers 2015 | TIE | longitudinal_treatment (Alliance formation (before/after)) | MID after alliance formation |  |
| 14 | Morrow 2017 | TIE | ptargdef (Target has defensive ally) | MID initiation |  |
| 15 | Kenwick & Vasquez 2017 | TIE | nfbdef (Target pact formed after birth) | MID initiation |  |
| 17 | Chen 2021 | TIE | extdeporig (Extended dependence (trade with target's allies)) | MID initiation |  |
| 25 | Malis 2021 | TIE | turnover (Ambassador turnover) | MID onset |  |
| 32 | Bas & Schub 2016 | TIE | msod (Mutual optimism) | MID onset |  |
| 42 | Lee & Mitchell 2012 | TIE | dflow_avegdplag (Bilateral FDI/GDP) | Escalation of ICOW claim to fatal MID |  |
| 5 | Fuhrmann & Sechser 2014 | TIE | defense_pact_nuclear_cow (Defense pact with nuclear state) | Fatal MID initiation |  |
| 54 | Bearce & Omori 2005 | TIE | high_organs (PTA high-level organs) | MID onset |  |
| 61 | Gelpi & Grieco 2008 | TIE | ordm1dpa (Trade dependence $\times$ initiator democracy) | MID initiation |  |
| 63 | Bussmann 2010 | TIE | smlinflowsgdpim (FDI inflows/GDP (lower)) | Fatal MID onset | [2026-11-13 audit] history anchored on the outcome year (DV is a t+1 lead; checked vs COW dyadic MID). |
| 64 | Goenner 2010 | TIE | L.open1_lev (Energy share of bilateral trade) | MID onset | [2026-11-13 audit] key term is Energy level (energy share of bilateral trade), not trade openness. |
| 65 | Hegre, Oneal & Russett 2010 | TIE | lnrtrade (Trade (log)) | Fatal MID |  |
| 66 | Peterson 2011 | TIE | lnot1xswtglo (Third-party trade $\times$ similarity) | MID initiation ($t+1$) |  |
| 68 | Kinne 2013 | TIE | requivlag5 (Five-year change in IGO network structural similarity) | MID initiation |  |
| 69 | Mitchell & Zawahri 2015 | TIE | infoexchange2 (River treaty information exchange) | Militarized river claim |  |
| 71 | Lupu & Greenhill 2017 | TIE | samecomm (Same network community) | Hostile MID onset ($t+1$) |  |
| 72 | DiLorenzo & Rooney 2018 | TIE | sols_count_exclusive (Coalition (SOLS) changes in exclusive allies) | MID onset |  |
| 73 | Zeng 2024 | TIE | ltXext (Strategic trade $\times$ externality) | Hostile MID onset ($t+1$) |  |
| 86 | Hafner-Burton & Montgomery 2012 | TIE | DependD (Dependence (DependD)) | MID onset ($t+1$) |  |
| 88 | Kinne 2013 | TIE | signaltriads (Signal triads) | MID initiation |  |
| 91 | Kang 2015 | TIE | c1allcap (Initiator's allies' capability) | MID initiation |  |
| ADD_T01 | Mazumder 2017 | TIE | PTACent1_lag (PTA centrality (autocratic challenger)) | MID initiation |  |
| ADD_T05 | Choi 2026 | TIE | ccgap (Relational status decline) | MID initiation |  |
| N06 | Edgerton (forthcoming) PSRM | TIE | cinc_pr_dff (Capability rank difference) | Fatal MID initiation |  |
| P8_003 | Goldsmith 2013 | TIE | sel_GDPShareGledLow (Trade/GDP (lower)) | MID onset ($t+1$) (selection equation) | [2026-11-13 audit] script outcome fcwmid251 corresponds to Model 3 (key coefficient -14.994 vs Model 1 -14.960). [2026-11-13 audit, PI] history anchored on the outcome year (DV is t+1); verified against COW dyadic MID data. |
| P8_008 | Nieman 2016 | TIE | FUS_SecHierA (US security subordination (status-quo eq.; + = less challenge)) | Conflict initiation |  |
| P8_012 | Choi 2022 | TIE | trade_x_nationalist_leader (Trade $\times$ nationalist leader) | MID onset |  |
| P8_015 | Isgandarov 2026 | TIE | jblag (Joint-body strength) | MID onset |  |
| P9_040 | Renshon 2016 | TIE | L1std_comm_fg_pr_dissat_A (Community status deficit) | MID initiation |  |
| P9_047 | Peterson 2015 | TIE | max_dyadic_dist (PTA trade distortion) | MID onset ($t+1$) |  |
| P9_048 | Bareis 2024 | TIE | oil_dep_x_degree (Oil dependence $\times$ network degree) | MID initiation |  |
| P9_051 | Rahman 2026 | TIE | l1_poor_rep_t_def (Target's own reputation as unreliable ally) | MID targeting |  |
| P9_056 | Choi 2025 | TIE | fbdef (Target pact from birth) | War initiation |  |
| P9_071 | Barbieri & Peters 2003 | TIE | salience (Trade salience) | MID onset ($t+1$) |  |
| P9_079 | Peterson & Drury 2011 | TIE | anysanctB (Third-party sanctions on target) | MID initiation ($t+1$) |  |
| P9_080 | Johnson & Joiner 2019 | TIE | allypc (Ally power change ($\geq$30\%)) | MID initiation |  |
| P9_082 | Pevehouse & Russett 2006 | TIE | demigo (Democratic IGOs) | Fatal MID onset ($t+1$) |  |
| P9_093 | Gartzke & Gleditsch 2022 | TIE | tie3 (Third-order alliance tie$^{\ddagger}$) | MID onset | [2026-11-13 PI] kept; note: the third-order alliance tie is a lagged contagion term (alliance links to dyads with disputes), so it embeds prior conflict in linked dyads. REBUILT (no Dataverse deposit found; round-3 PI-supplied PDF). O-R all dyads 1950-92 (or_trade/DOR), DV Maoz dyadic MID onset in t+1 incl. joiners (mzmid1_all), continuing dispute years dropped using GML dyadic disputes; third-order alliance tie = dyad pairing a disputant with an ally of its opponent in a GML dispute ongoing in t (COW v4.1 alliances, any type); fourth-order tie, MID inbetweenness ratio and Max(A,B IBR) x logge |
| P9_094 | Krause 2004 | TIE | arms (Arms transfers) | MID initiation |  |
| X044 | Owsiak 2016 | TIE | settlem (Border settlement) | MID onset |  |
| X055 | Li, Bradshaw, Clary & Cranmer 2017 Sci Adv | TIE | Contiguity x one-degree alliance (Table 1 Model 1) (Contiguity $\times$ one-degree (direct) alliance) | Serious MID onset |  |
| X061 | Johnson, Leeds & Wu 2015 II | TIE | allies_s (Ally credibility) | MID initiation |  |
| X066 | Aydin 2010 | TIE | instrade (Trade through regional institutions) | MID initiation |  |
| X067 | Benson, Bentley & Ray 2013 | TIE | idgNcrevA (Unconditional commitment $\times$ revisionism) | MID initiation |  |
| X068 | Chatagnier & Lim 2020 | TIE | did.wto (Joint GATT/WTO membership $\times$ WTO era) | MID onset |  |
| X069 | Tomashevskiy 2025 | TIE | friendlassonotime (Foreign policy affinity (three-step)) | War onset | [2026-11-13 audit] outcome is war onset (Table 5, Dafoe A+B+C replication of Gartzke 2007 M7). |
| X078 | Peterson & Rudloff 2015 | TIE | F_signed_noforce (PTA signed, not in force) | MID onset |  |
| X080 | Ryou-Ellison & Gold 2020 | TIE | defense (Defense pact) | MID over claim |  |
| X087 | Owsiak 2012 | TIE | settlem (Border settlement) | MID onset |  |
| inv8 | Long et al. 2007 | TIE | nonaggression (Alliance non-aggression provisions) | MID with mutual force |  |

### pool_status = outcome_overlap (32)

| id | study | family | key term (label) | outcome label | Nov 13–14 notes |
|---|---|---|---|---|---|
| C1 | Gartzke 2007 AJPS | CAPITALISM | capopenl (Capital openness) | MID involvement ($t+1$) | [2026-11-14] outcome codes every year of a dispute (data check vs COW MID 4.02), so a dispute that began in the history window and is ongoing counts in both the outcome and history; moved to the appendix expanded set (pool_status outcome_overlap). |
| P8_048 | Ekmekci 2014 | CAPITALISM | development_low (Lower development) | MID involvement | [2026-11-14] outcome codes every year of a dispute (data check vs COW MID 4.02), so a dispute that began in the history window and is ongoing counts in both the outcome and history; moved to the appendix expanded set (pool_status outcome_overlap). |
| P9_085 | Mousseau 2000 | CAPITALISM | devdem (Developed democracy) | MID involvement ($t+1$) | [2026-11-14] outcome codes every year of a dispute (data check vs COW MID 4.02), so a dispute that began in the history window and is ongoing counts in both the outcome and history; moved to the appendix expanded set (pool_status outcome_overlap). |
| C2 | Choi 2011 ISQ | DEMOCRACY | demlo (Lower democracy) | MID involvement ($t+1$) | [2026-11-14] outcome codes every year of a dispute (data check vs COW MID 4.02), so a dispute that began in the history window and is ongoing counts in both the outcome and history; moved to the appendix expanded set (pool_status outcome_overlap). |
| D2 | Dafoe 2011 AJPS | DEMOCRACY | demlo (Lower democracy) | MID involvement ($t+1$) | [2026-11-13 audit, PI] history anchored on the outcome year (DV is t+1); verified against COW dyadic MID data. [2026-11-14] outcome codes every year of a dispute (data check vs COW MID 4.02), so a dispute that began in the history window and is ongoing counts in both the outcome and history; moved to the appendix expanded set (pool_status outcome_overlap). |
| P8_030 | Chiba & Fang 2014 | DEMOCRACY | xbi_netdem1 (Net democracy) | Escalation under IO opposition | [2026-11-14 P3] outcome is escalation within an ongoing confrontation (month-level, after a threat); history window can contain the confrontation onset; moved to the appendix expanded set (outcome_overlap). |
| P8_036 | McCallister 2016 | DEMOCRACY | interact_region_diss (Regime dissimilarity $\times$ democratic strength) | MID involvement | [2026-11-14 P2] outcome codes continuing dispute years (COW MID 4.02 check: share of outcome years that are COW onsets well below share that are COW MID-years; or TRIANGLE dispute1); moved to the appendix expanded set (outcome_overlap). |
| P8_042 | Boldt 2026 | DEMOCRACY | att_mid_pool3 (Democratic transition (pooled ATT)) | MID participation | [2026-11-14] outcome codes every year of a dispute (data check vs COW MID 4.02), so a dispute that began in the history window and is ongoing counts in both the outcome and history; moved to the appendix expanded set (pool_status outcome_overlap). |
| P9_061 | Cederman 2001 | DEMOCRACY | dem_year (Year $\times$ joint democracy) | MID involvement | [2026-11-13 audit, PI] history anchored on the outcome year (DV is t+1); verified against COW dyadic MID data. REBUILT (no archive). TRIANGLE PR dyads 1885-1991 minus world wars and 2-yr shadows; D = both Polity IV DEMOC >= 6; D x Year and (1-D) x Year; INTER, COLD, allies, weaker/stronger CINC, lower NMC energy per capita. Published values from the article table. [2026-11-14 P2] outcome codes continuing dispute years (COW MID 4.02 check: share of outcome years that are COW onsets well below share that are COW MID-years; or TRIANGLE dispute1); moved to the appendix expanded set (outcome_overla |
| P9_089 | Reuveny & Li 2003 | DEMOCRACY | DEML_hat (Lower democracy (instrumented)) | MID involvement | [2026-11-13 audit, PI] history anchored on the outcome year (DV is t+1); verified against COW dyadic MID data. REBUILT (no archive; round-2 estimator build, simultaneous equations). Hand-coded two-stage system: OLS DEML/DEMH on all exogenous variables, probit of MID involvement (TRIANGLE dispute1) on fitted DEML/DEMH plus contiguity, alliance, log capability ratio, lower dependence, lower growth; second-stage SEs clustered on dyad as in the article (fn. 25, no Maddala correction). TRIANGLE PR dyads 1950-92; Polity IV 2014; WDI CPI inflation (from 1960); peacesciencer GDP pc; COW national trade |
| X033 | Rasler & Thompson 2001 | DEMOCRACY | jointdem (Joint democracy) | MID involvement | [2026-11-14] outcome codes every year of a dispute (data check vs COW MID 4.02), so a dispute that began in the history window and is ongoing counts in both the outcome and history; moved to the appendix expanded set (pool_status outcome_overlap). |
| X037 | Altman, Rojas-de-Galarreta & Urdinez 2021 | DEMOCRACY | dmXds (Democracy mean $\times$ spread) | Fatal MID involvement | [2026-11-14 P2] outcome codes continuing dispute years (COW MID 4.02 check: share of outcome years that are COW onsets well below share that are COW MID-years; or TRIANGLE dispute1); moved to the appendix expanded set (outcome_overlap). |
| X040 | Bayer & Bernhard 2010 | DEMOCRACY | jointdem6 (Joint democracy (dichotomous Polity)) | MID involvement | [2026-11-13 audit, PI] history anchored on the outcome year (DV is t+1); verified against COW dyadic MID data. [2026-11-13 audit, PI] pooled headline: authors' preferred dichotomous joint democracy (Table 1 col 2). [2026-11-14] outcome codes every year of a dispute (data check vs COW MID 4.02), so a dispute that began in the history window and is ongoing counts in both the outcome and history; moved to the appendix expanded set (pool_status outcome_overlap). |
| X085 | Warren 2016 JPR | DEMOCRACY | con_est: same regime x democracy (sameX x altX; joint-democracy effect on conflict ties, Table I si6con) (Joint democracy) | MID ties | [2026-11-14] outcome codes every year of a dispute (data check vs COW MID 4.02), so a dispute that began in the history window and is ongoing counts in both the outcome and history; moved to the appendix expanded set (pool_status outcome_overlap). |
| X092 | Goldsmith 2006 | DEMOCRACY | smldmat (Lower democracy) | MID involvement | [2026-11-13 audit] removed from the sensitivity set: after anchoring history on the outcome year, the dyad-only GEE AR(1) split does not converge (b ~ -6e12); main battery unaffected. [2026-11-13 audit] history anchored on the outcome year (DV is a t+1 lead; checked vs COW dyadic MID). [2026-11-14] outcome codes every year of a dispute (data check vs COW MID 4.02), so a dispute that began in the history window and is ongoing counts in both the outcome and history; moved to the appendix expanded set (pool_status outcome_overlap). |
| 100 | Peterson & Quackenbush 2010 | TIE | imp_interdep (Recipient $\times$ imposer trade dependence) | MID recurrence | [2026-11-14 P3] Cox risk set is peace spells that begin with a MID; dyadic history counts the spell-opening MID when it began within 5 years; moved to the appendix expanded set (outcome_overlap). |
| 101 | Xiang 2016 | TIE | smigoabi (Shared IGO membership) | MID involvement ($t+1$) | [2026-11-14] outcome codes every year of a dispute (data check vs COW MID 4.02), so a dispute that began in the history window and is ongoing counts in both the outcome and history; moved to the appendix expanded set (pool_status outcome_overlap). |
| 115 | Edgerton 2024 | TIE | out_group [FE logit] (Out-group ties) | Fatal MID involvement | [2026-11-14 P2] outcome codes continuing dispute years (COW MID 4.02 check: share of outcome years that are COW onsets well below share that are COW MID-years; or TRIANGLE dispute1); moved to the appendix expanded set (outcome_overlap). |
| 118 | Hafner-Burton & Montgomery 2006 | TIE | CentDif (Prestige difference) | MID involvement ($t+1$) | [2026-11-14] outcome codes every year of a dispute (data check vs COW MID 4.02), so a dispute that began in the history window and is ongoing counts in both the outcome and history; moved to the appendix expanded set (pool_status outcome_overlap). |
| 119 | Maoz et al. 2006 | TIE | lgintegse (Integrated structural equivalence (lagged)) | MID under way | [2026-11-14] outcome codes every year of a dispute (data check vs COW MID 4.02), so a dispute that began in the history window and is ongoing counts in both the outcome and history; moved to the appendix expanded set (pool_status outcome_overlap). |
| 55 | Kimball 2006 | TIE | allform2 (Alliance formation) | Reciprocated use-of-force MID | [2026-11-14 P2] outcome codes continuing dispute years (COW MID 4.02 check: share of outcome years that are COW onsets well below share that are COW MID-years; or TRIANGLE dispute1); moved to the appendix expanded set (outcome_overlap). |
| 62 | Dorussen & Ward 2010 | TIE | lntradelinks3_index (Triadic trade (ln)) | MID involvement ($t+1$) | [2026-11-14] outcome codes every year of a dispute (data check vs COW MID 4.02), so a dispute that began in the history window and is ongoing counts in both the outcome and history; moved to the appendix expanded set (pool_status outcome_overlap). |
| P9_028 | Russett, Oneal & Davis 1998 | TIE | igo_l1 (Shared IGOs) | MID involvement | [2026-11-13 audit, PI] history anchored on the outcome year (DV is t+1); verified against COW dyadic MID data. REBUILT (no archive). TRIANGLE PR dyads 1950-85 + SDP growth; IGO and dependence lagged 1 yr from TRIANGLE; allies = direct or both allied with US. Published values from the article table. [2026-11-14] outcome codes every year of a dispute (data check vs COW MID 4.02), so a dispute that began in the history window and is ongoing counts in both the outcome and history; moved to the appendix expanded set (pool_status outcome_overlap). |
| P9_029 | Gartzke 1998 | TIE | aff_l1 (UN affinity) | MID involvement | [2026-11-13 audit, PI] history anchored on the outcome year (DV is t+1); verified against COW dyadic MID data. REBUILT (no archive) from Russett-Oneal TRIANGLE + Gartzke Affinity 2010 (S-score; paper used Spearman rho) + SDP growth. Family note: Oneal/Russett-type data (same base as X076). [2026-11-14] outcome codes every year of a dispute (data check vs COW MID 4.02), so a dispute that began in the history window and is ongoing counts in both the outcome and history; moved to the appendix expanded set (pool_status outcome_overlap). |
| P9_043 | Peterson & Thies 2012 | TIE | intra_index_final (Intra-industry trade) | Fatal MID involvement ($t+1$) | [2026-11-14 P2] outcome codes continuing dispute years (COW MID 4.02 check: share of outcome years that are COW onsets well below share that are COW MID-years; or TRIANGLE dispute1); moved to the appendix expanded set (outcome_overlap). |
| P9_050 | Li & Yang 2026 | TIE | ext_x_tradecent (Extended trade dependence $\times$ alliance centrality) | MID involvement (side A) | [2026-11-14 P2] outcome codes continuing dispute years (COW MID 4.02 check: share of outcome years that are COW onsets well below share that are COW MID-years; or TRIANGLE dispute1); moved to the appendix expanded set (outcome_overlap). |
| P9_053 | Smith & Vreeland 2023 | TIE | unscA (UN Security Council membership) | MID on initiating side (incl.\ joiners) | [2026-11-14 P2] outcome codes continuing dispute years (COW MID 4.02 check: share of outcome years that are COW onsets well below share that are COW MID-years; or TRIANGLE dispute1); moved to the appendix expanded set (outcome_overlap). |
| P9_092 | Millerd 2023 | TIE | gendense_lag (Shared gender-equal IOs) | MID involvement | [2026-11-14] outcome codes every year of a dispute (data check vs COW MID 4.02), so a dispute that began in the history window and is ongoing counts in both the outcome and history; moved to the appendix expanded set (pool_status outcome_overlap). |
| X031 | Kuo & Blankenship 2022 | TIE | jally_new2_target (New joint exercise with ally) | MID escalation | [2026-11-14 P2] outcome codes continuing dispute years (COW MID 4.02 check: share of outcome years that are COW onsets well below share that are COW MID-years; or TRIANGLE dispute1); moved to the appendix expanded set (outcome_overlap). |
| X059 | Maoz et al. 2007 | TIE | lag_imbal (Imbalanced relations (lagged)) | Ongoing MID | [2026-11-14] outcome codes every year of a dispute (data check vs COW MID 4.02), so a dispute that began in the history window and is ongoing counts in both the outcome and history; moved to the appendix expanded set (pool_status outcome_overlap). |
| X093 | Martin, Mayer & Thoenig 2008 RESTUD | TIE | lopenl (Multilateral openness ($t-4$)) | MID involvement (hostility level $\geq$ 3) | [2026-11-14] outcome codes every year of a dispute (data check vs COW MID 4.02), so a dispute that began in the history window and is ongoing counts in both the outcome and history; moved to the appendix expanded set (pool_status outcome_overlap). |
| inv34 | Suhas 2023 | TIE | onlyweakallied (Only the weaker rival allied) | Weaker rival initiates MID | [2026-11-13 audit] table cite corrected to Table 4 Model 1. [2026-11-14 P2] outcome codes continuing dispute years (COW MID 4.02 check: share of outcome years that are COW onsets well below share that are COW MID-years; or TRIANGLE dispute1); moved to the appendix expanded set (outcome_overlap). |

### pool_status = author_null (1)

| id | study | family | key term (label) | outcome label | Nov 13–14 notes |
|---|---|---|---|---|---|
| P8_032 | Choi 2015 | DEMOCRACY | sysdems (nan) | nan |  |

### pool_status = duplicate (4)

| id | study | family | key term (label) | outcome label | Nov 13–14 notes |
|---|---|---|---|---|---|
| C7 | Dafoe & Kelsey 2014 (JPR) | CAPITALISM | capopenl (nan) | nan |  |
| 99 | Gartzke & Hewitt 2010 | TIE | sun2cati (nan) | nan |  |
| P9_004 | Kang 2017 | TIE | c1allcap (nan) | nan |  |
| X062 | Leeds & Johnson 2017 | TIE | ptargdef (nan) | nan |  |

### pool_status = methods_reanalysis (5)

| id | study | family | key term (label) | outcome label | Nov 13–14 notes |
|---|---|---|---|---|---|
| N04 | Minhas et al. 2022 (Reiter & Stam app.) | DEMOCRACY | pers_democ_directed_dyad (nan) | nan |  |
| N05 | Minhas et al. 2022 (Weeks app.) | DEMOCRACY | boss (nan) | nan |  |
| 50 | Gartzke & Li 2003 | TIE | dependll (nan) | nan |  |
| N03 | Fritz, Mehrl, Thurner & Kauermann 2025 | TIE | GWESF- (MID tie between states sharing friends) (nan) | nan |  |
| N07 | Westveld & Hoff 2011 | TIE | alliance_level (nan) | nan |  |

### pool_status = not_comparable (5)

| id | study | family | key term (label) | outcome label | Nov 13–14 notes |
|---|---|---|---|---|---|
| P8_041 | Barnhart et al. 2020 | DEMOCRACY | lagjointMSdem (nan) | nan |  |
| X005 | Hensel, Goertz & Diehl 2000 JOP | DEMOCRACY | jointdem (nan) | nan |  |
| 52 | Benson 2004 | TIE | DyadTotalDissim (nan) | nan | [2026-11-13 audit, PI] key switched from Model 1 DyadEconDissim (base model) to Model 2 DyadTotalDissim, the author's preferred model. |
| P9_017 | Anderson, Mitchell & Schilling 2016 | TIE | i2_numigo (nan) | nan |  |
| X077 | Hwang 2010 II | TIE | Ipow (nan) | nan |  |

### pool_status = not_reproduced (4)

| id | study | family | key term (label) | outcome label | Nov 13–14 notes |
|---|---|---|---|---|---|
| X020 | Park & Colaresi 2014 | DEMOCRACY | lowdem (Lower democracy) | MID onset | [2026-11-14 P2] key term matches (-0.0297 vs -0.030, same N) but peace-year/spline terms and controls do not match published Table 1 Model 4 (peace years -0.468 vs -0.350; spline ~15x); baseline spec not reproduced. Not pooled. |
| 4 | Maoz 2009 | TIE | lgstrintrdep (nan) | nan |  |
| 57 | Powers 2006 | TIE | lagreiall (Shared REI alliance) | MID initiation | [2026-11-14] with the printed target-REI-alliance control added, lagreiall = -0.963 (SE 0.543, p = .08) vs printed -1.184 (.441); script verdict NOT_REPRODUCED. Not pooled. |
| ADD_T03 | Leeds & Mattes 2007 | TIE | notnonagg (Alliance beyond non-aggression) | MID involvement ($t+1$) | [2026-11-13 audit] published_b corrected to printed Table 6 col 3 value; was -0.303 (0.149), which no printed column shows. [2026-11-14] geeglm AR(1) port does not reproduce the printed Table (key -0.303 vs -0.268; NAP -0.466 vs -0.178); the gate had been hard-coded. Not pooled. |

### pool_status = null_common_sample (17)

| id | study | family | key term (label) | outcome label | Nov 13–14 notes |
|---|---|---|---|---|---|
| P8_044 | Mousseau 2009 | CAPITALISM | dml (nan) | nan |  |
| D1 | Dafoe, Oneal & Russett 2013 ISQ | DEMOCRACY | jtdemxpropdem (nan) | nan |  |
| N08 | Ward, Siverson & Cao 2007 | DEMOCRACY | joint_democracy_polity_i_x_polity_j (nan) | nan |  |
| X016 | Lektzian, Prins & Souva 2010 | DEMOCRACY | demnonternonriv (nan) | nan |  |
| X019 | Gibler 2014 | DEMOCRACY | jointdem7 (nan) | nan |  |
| 108 | Boehmer, Gartzke & Nordstrom 2004 | TIE | IGOmedst (nan) | nan |  |
| 43 | Lupu & Traag 2013 | TIE | stc_low (nan) | nan |  |
| 47 | Feldman, Eiran & Rubin 2021 | TIE | third_navA (nan) | nan |  |
| 53 | Kim & Rousseau 2005 | TIE | depend_m (nan) | nan |  |
| 58 | Goldsmith 2007 | TIE | lngleddeplowAsia (nan) | nan |  |
| 89 | Wright & Rider 2014 | TIE | terrXdef (nan) | nan |  |
| P9_002 | Conrad 2017 | TIE | interaction2 (nan) | nan |  |
| P9_054 | Peterson 2013 | TIE | ldep1_OR (nan) | nan |  |
| X057 | Pomeroy, Dasandi & Mikhaylov 2019 PLOS ONE | TIE | edgecov weak-signal community tie (Table 1 Model 2) (nan) | nan |  |
| X082 | Green, Kim & Yoon 2001 | TIE | lowdep (nan) | nan |  |
| X090 | Choi 2024 | TIE | extdeporig (nan) | nan |  |
| X091 | Mousseau, Hegre & Oneal 2003 | TIE | deplo (nan) | nan |  |

### pool_status = out_of_scope (22)

| id | study | family | key term (label) | outcome label | Nov 13–14 notes |
|---|---|---|---|---|---|
| C6 | Gartzke & Hewitt 2010 | CAPITALISM | capopnl2 (nan) | nan |  |
| P8_037 | Goldsmith et al. 2017 | DEMOCRACY | BGcomp1pXpolity22p (nan) | nan |  |
| X022 | Poznansky & Scroggs 2016 | DEMOCRACY | polyprod (nan) | nan |  |
| X025 | Gelpi & Grieco 2001 JCR | DEMOCRACY | defpol4_notenure (nan) | nan |  |
| 104 | Karreth 2018 | TIE | igo_lev3_count_use_bc (nan) | nan |  |
| N01 | Cranmer & Desmarais 2011 | TIE | Two-star (popularity; authors' own claim) (nan) | nan |  |
| N02 | Chen 2021 | TIE | MID-Civil Conflict alternating k-star(1): state's MID ties cluster with its civil-conflict ties (nan) | nan |  |
| P9_001 | Rahman 2026 | TIE | dca_dum (nan) | nan |  |
| P9_008 | Salehyan 2008 | TIE | logref1 (nan) | nan |  |
| P9_012 | Johnson, Leeds & Savun 2024 | TIE | gov_defense_pure (nan) | nan |  |
| P9_016 | Johnson 2016 | TIE | tgt_defense (nan) | nan |  |
| P9_023 | Krustev 2006 | TIE | kdeplo (nan) | nan |  |
| P9_036 | Norrevik & Sarwari 2025 | TIE | dca (nan) | nan |  |
| P9_041 | Fausett & Volgy 2010 | TIE | loglagTIGO (nan) | nan |  |
| P9_045 | Kim 2014 | TIE | Ldepkg1 (nan) | nan |  |
| P9_046 | Peterson 2014 | TIE | lnsg_strat_int (nan) | nan |  |
| P9_057 | Fang 2014 | TIE | ally_cost (nan) | nan |  |
| P9_065 | Zeng 2021 | TIE | cost_x_sal (nan) | nan |  |
| P9_078 | Hewitt 2003 | TIE | smigoabi (nan) | nan |  |
| X060 | Beardsley 2024 JOP | TIE | threat.net: j's support military spending (altX; H1, marginal effect at mean closeness) (nan) | nan |  |
| X079 | Bell & Long 2016 | TIE | l_interdep (nan) | nan |  |
| X088 | Atkinson, Owsiak, Jackson & Buechler 2024 | TIE | new_atop_ally (nan) | nan |  |

### pool_status = rebuild_outside_rule (20)

| id | study | family | key term (label) | outcome label | Nov 13–14 notes |
|---|---|---|---|---|---|
| P8_043 | Francis 2009 | CAPITALISM | hc_lo (nan) | nan |  |
| P9_067 | Mousseau 2005 | CAPITALISM | dem_x_dev (nan) | nan |  |
| P9_075 | Helfstein 2012 | CAPITALISM | debt (nan) | nan |  |
| P9_060 | Debs & Goemans 2010 | DEMOCRACY | civilian (nan) | nan |  |
| P9_066 | Kim & James 2010 | DEMOCRACY | demdyad (nan) | nan |  |
| X017 | Gartzke & Weisiger 2013 ISQ | DEMOCRACY | dxs (nan) | nan |  |
| X039 | Quackenbush & Rudy 2009 | DEMOCRACY | bothdem (nan) | nan |  |
| 20 | Oneal & Russett 2001 | TIE | smldep (nan) | nan |  |
| 78 | Benson 2005 | TIE | sal_x_size (nan) | nan |  |
| P8_004 | Aaronson et al. 2015 | TIE | wto_membership_x_trade (nan) | nan |  |
| P8_005 | Beyene 2015 | TIE | trade_dependence (nan) | nan |  |
| P9_030 | Bearce, Flanagan & Floros 2006 | TIE | lnprep_x_ally (nan) | nan |  |
| P9_031 | Mansfield & Pevehouse 2000 | TIE | PTA_traBilMax (nan) | nan |  |
| P9_033 | Hegre 2004 | TIE | e_s2 (nan) | nan |  |
| P9_034 | Barbieri 1996 | TIE | interdep (nan) | nan |  |
| P9_035 | Johnson & Wolford 2023 | TIE | out_allypc (nan) | nan |  |
| P9_074 | Souva & Prins 2006 | TIE | tdep (nan) | nan |  |
| P9_087 | McDonald 2015 | TIE | gpdem_x_dem (Great-power-supported democracy $\times$ democracy) | MID onset | [2026-11-14] full rebuilt sample gives 2.003x published (outside 0.5-2); GPDEM main term 2.20 vs 0.897. Not pooled. |
| X074 | Oneal & Russett 2005 | TIE | lowdep (nan) | nan |  |
| X076 | Gartzke 2000 | TIE | lowdem (nan) | nan |  |
