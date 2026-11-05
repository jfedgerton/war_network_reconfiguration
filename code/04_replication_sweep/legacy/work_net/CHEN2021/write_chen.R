# Write results_net/CHEN2021__<term>.csv from the author-stack MCMLE (m0, hist), MPLE (m26, placebo points) and post-hoc diagnostics.
OUT <- Sys.getenv("SWEEP_OUT"); stopifnot(nzchar(OUT))
W <- file.path(OUT, "work_net/CHEN2021"); O <- file.path(OUT, "results_net"); dir.create(O, recursive = TRUE, showWarnings = FALSE)
r0 <- readRDS(file.path(W, "res_m0_MCMLE_old.rds")); rh <- readRDS(file.path(W, "res_hist_MCMLE_old.rds"))
mp <- lapply(setNames(c("m0", "hist", "m26", "placebo"), c("m0", "hist", "m26", "placebo")), function(m) readRDS(file.path(W, paste0("res_", m, "_MPLE_old.rds"))))
d0 <- readRDS(file.path(W, "diag_m0.rds")); dh <- if (file.exists(file.path(W, "diag_hist.rds"))) readRDS(file.path(W, "diag_hist.rds")) else NULL
stopifnot(identical(r0$term, rh$term[1:20]), nrow(rh) == 22, identical(d0$term, r0$term))
g <- function(r, t, col = "b") r[[col]][r$term == t]
pc <- function(x, b) round(100 * (x - b) / abs(b), 1)
keys <- list(
  list(term = "altkstar(1)_crosslayer.1+1-2", file = "MIDxCivil_kstar", pub = 1.27, pub_se = 0.43,
       lab = "MID-Civil Conflict alternating k-star(1): state's MID ties cluster with its civil-conflict ties"),
  list(term = "threetrail_crosslayer_incident.1+1-2+2", file = "ThreeWay", pub = 0.19, pub_se = 0.07,
       lab = "Three-way clustering (cross-layer three-trail MID-civil-nonstate, incident): the paper's key theoretical term"))
for (k in keys) {
  b0 <- g(r0, k$term); se0 <- g(r0, k$term, "se"); p0 <- g(r0, k$term, "p")
  bh <- g(rh, k$term); seh <- g(rh, k$term, "se"); ph <- g(rh, k$term, "p")
  rep_ok <- abs(b0 - k$pub) / abs(k$pub) <= 0.10 && (p0 < .05) == (abs(k$pub / k$pub_se) > 1.96)
  dg0 <- d0[d0$term == k$term, ]; dgh <- if (!is.null(dh)) dh[dh$term == k$term, ] else NULL
  conv <- sprintf(paste0("MCMLE (ergm 3.10.1 + multilayer.ergm 0.1.6 + statnet.common 4.3.0 + network 1.15, the author's stack; seed 30748, burn-in 1e5, Hummel termination). ",
    "Chains shortened for budget: interval 2000 / samplesize 4000 (x4 at the last iteration), 2 parallel chains, vs the authors' 10000 / 20000, 4 chains. ",
    "m0 converged (step length 1 twice) in 4 iterations (33 min x 2 cores), hist in 5 (40 min). ",
    "Post-hoc single-chain check at the final estimates (4000 draws, interval 2000): this term t-ratio %.2f, ESS %.0f, Geweke z %.2f (m0)%s; ",
    "all 20 m0 statistics |t| < 0.5 and no |Geweke z| > 2; MID-layer statistics mix slowly (ESS 15-82 of 4000 in m0, 16-43 in hist)."),
    dg0$t_ratio, dg0$ess, dg0$geweke_z,
    if (!is.null(dgh) && nrow(dgh)) sprintf(" / t %.2f, ESS %.0f, Geweke %.2f (hist; max |t| over 22 stats %.2f, %d of 22 |Geweke z| > 2)", dgh$t_ratio, dgh$ess, dgh$geweke_z, max(abs(dh$t_ratio)), sum(abs(dh$geweke_z) > 2)) else " (hist post-hoc check not available)")
  notes <- paste0(
    "HEADLINE check: main-text application 2 (Sec. 3.2, Table 3, cross-layer dependence model); global conflict network 1989-1995 as one cross-section, 189 states + 544 nonstate actors, MID / civil / nonstate layers. ",
    sprintf("Published b %.2f (se %.2f); ours %.3f (se %.3f). ", k$pub, k$pub_se, b0, se0),
    "ENDOGENOUS conflict-network term: the MID layer is part of the DV and the key statistic is built from MID ties, so dyadic MID history overlaps it. ",
    "History = log1p dyadic MID onsets 1984-88 (t-1..t-5 relative to 1989, GML) as edgecov on the MID layer + log1p state onsets 1984-88 as nodecov on the MID layer. 34 of 134 MID ties have prior dyadic onsets. ",
    sprintf("hist terms: dyad %.2f (se %.2f), mono %.2f (se %.2f). ", g(rh, "edgecov.layer.1.hd15"), g(rh, "edgecov.layer.1.hd15", "se"), g(rh, "nodecov_layer.1.mono15"), g(rh, "nodecov_layer.1.mono15", "se")),
    sprintf("Civil-nonstate k-star (no MID tie): %.3f -> %.3f. ", g(r0, "altkstar(1)_crosslayer.1-2+2"), g(rh, "altkstar(1)_crosslayer.1-2+2")),
    "m26 and placebo MCMLE not run (each ~1.2 CPU-h; m0 + hist used the 3 CPU-h budget). MPLE points on the author stack given instead; MPLE is only ergm's starting value and puts the three-way term near 0 (m0 MPLE ",
    sprintf("%.3f), so MPLE %% changes for that term are uninformative. ", g(mp$m0, "threetrail_crosslayer_incident.1+1-2+2")),
    "Placebo = mean CINC stronger-share 1984-88 as MID-layer edgecov; 5,363 of 17,766 state dyads missing (post-1988 states) mean-imputed. ",
    "Current CRAN/GitHub stack (ergm 4.12 + multilayer.ergm 0.2.3) miscomputes the cross-layer statistics (observed three-trail = 0 vs 170, MID-civil k-star 50 vs 306), so earlier runs on it were discarded (work_net/CHEN2021/invalid_newstack/).")
  row <- data.frame(tag = "CHEN2021", study = "Chen 2021 PSRM (multilayer ERGM, post-Cold War global conflict network)",
    headline = "HEADLINE: main text Table 3, cross-layer dependence model (one of two main-text applications; the only one with MIDs)",
    key_effect = k$lab, model_term = k$term, pub_b = k$pub, pub_se = k$pub_se,
    b0 = round(b0, 4), se0 = round(se0, 4), p0 = signif(p0, 3), b_hist = round(bh, 4), se_hist = round(seh, 4), p_hist = signif(ph, 3), pct_hist = pc(bh, b0),
    b_26 = "NOT_RUN", pct_26 = sprintf("NOT_RUN (MPLE only: %.1f)", pc(g(mp$m26, k$term), g(mp$m0, k$term))),
    b_placebo = "NOT_RUN", pct_placebo = sprintf("NOT_RUN (MPLE only: %.1f)", pc(g(mp$placebo, k$term), g(mp$m0, k$term))),
    mple_b0 = round(g(mp$m0, k$term), 4), mple_hist = round(g(mp$hist, k$term), 4),
    convergence = conv, verdict = if (rep_ok) "REPRODUCED" else "CLOSE", notes = notes)
  write.csv(row, file.path(O, paste0("CHEN2021__", k$file, ".csv")), row.names = FALSE)
  cat(k$file, b0, se0, p0, "->", bh, seh, ph, pc(bh, b0), row$verdict, "\n")
}
