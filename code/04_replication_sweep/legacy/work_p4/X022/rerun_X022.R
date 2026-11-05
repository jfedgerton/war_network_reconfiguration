# X022 Poznansky & Scroggs 2016 ISQ, "Ballots and Blackmail: Coercive Diplomacy and the Democratic Peace"
# Key terms (decided from abstract BEFORE estimating): joint democracy. Headline = Table 3 Model 1 (polyprod, Polity product).
#   Abstract: "pairs of democracies are significantly less likely to engage in coercive diplomacy"
# The same Model 1 spec is re-run with the paper's two alternative joint-democracy measures (Table 4 M1 polity_low,
# Table 5 M1 polity_both), one row each. Outcome = militarized compellent threat (MCT, Sechser), not MIDs.
# Estimator: logit, SE clustered on dyad_id (main.do). No author log / table access -> UNVERIFIED.
ROOT <- Sys.getenv("SWEEP_DATA"); OUT <- Sys.getenv("SWEEP_OUT"); CODE <- Sys.getenv("SWEEP_CODE")
stopifnot(nzchar(ROOT), nzchar(OUT), nzchar(CODE))
source(file.path(CODE, "work_p4/b6/common_b6.R"))
WO <- file.path(OUT, "work_p4/X022"); dir.create(WO, recursive = TRUE, showWarnings = FALSE)  # port: outputs
W <- file.path(ROOT, "work_p4/X022"); sink(file.path(WO, "log_X022.txt"), split = TRUE)
m <- list(xid = "X022", short_cite = "Poznansky & Scroggs 2016", journal = "ISQ", year = 2016, doi = "10.1093/isq/sqw016")
d <- haven::read_dta(file.path(W, "PS.dta"), col_select = c(ccode1, ccode2, year, dyad_id, mct, polyprod, polity_low, polity_both,
  lowdpnd, cincratio, alliance1, absidealdiff, conttype2, time, time_sq, time_cu))
d <- as.data.frame(haven::zap_labels(d))
cat("rows:", nrow(d), " dup ccode1-ccode2-year:", sum(duplicated(d[, c("ccode1", "ccode2", "year")])), " years:", range(d$year), "\n")
stopifnot(all(d$mct %in% c(0, 1, NA)))

ctrl <- "lowdpnd + cincratio + alliance1 + absidealdiff + conttype2 + time + time_sq + time_cu"
spec <- list(polyprod = "Table 3 Model 1", polity_low = "Table 4 Model 1", polity_both = "Table 5 Model 1")
for (k in names(spec)) {
  f0 <- as.formula(paste("mct ~", k, "+", ctrl))
  Wk <- file.path(WO, k); dir.create(Wk, recursive = TRUE, showWarnings = FALSE)
  r <- fit_keys(f0, d, k, "~dyad_id", "logit", file.path(Wk, "cache_repro.rds"))
  ds <- common_sample(d, f0, "dyad_id", "ccode1", "ccode2", "year")
  M <- run_six(f0, ds, k, "~dyad_id", "logit", Wk); rm(ds); gc()
  # ---- sensitivity hook: runs only under run_sensitivity.R (SWEEP_SENS=1), never in the main sweep ----
  if (Sys.getenv("SWEEP_SENS") == "1" && k == "polity_low") {   # headline term only (the loop also runs polyprod and polity_both)
    source(file.path(CODE, "hist/sens_helper.R"))
    # ds is removed on the run_six() line above: rebuild the same common sample with the script's own call (same d, same f0)
    sens_run(id = "X022", dat = common_sample(d, f0, "dyad_id", "ccode1", "ccode2", "year"), y = "mct", tie = "polity_low",
             ctrl = c("lowdpnd", "cincratio", "alliance1", "absidealdiff", "conttype2", "time", "time_sq", "time_cu"), fe = NULL,
             fit = function(dat, add) { r <- fit_keys(if (length(add)) update(f0, as.formula(paste(". ~ . +", paste(add, collapse = " + ")))) else f0, dat, "polity_low", "~dyad_id", "logit")$k$polity_low
               c(b = unname(r[["b"]]), se = unname(r[["se"]]), p = unname(r[["p"]])) })
  }
  write_p4(m, "DEMOCRACY", "KEY", k, "abstract: 'pairs of democracies are significantly less likely to engage in coercive diplomacy'",
    "UNVERIFIED", NA, unname(r$k[[k]]["b"]), "logit, dyad-clustered SE (nondirected dyad-year, MCT onset)", r$n, M,
    "Y (time, time^2, time^3 since last MCT)",
    paste0(spec[[k]], " (main.do) full sample; no author log and paper table not accessible -> UNVERIFIED. DV is militarized compellent threats (MCT), not MIDs; history = GML MID onsets. ",
           "Years ", paste(range(d$year), collapse = "-"), " (inside GML range); dropped vs complete-case: ", r$n - M$m0$n,
           ". bD (lagged joint democracy) is not a placebo for a democracy key term."))
}
sink()
