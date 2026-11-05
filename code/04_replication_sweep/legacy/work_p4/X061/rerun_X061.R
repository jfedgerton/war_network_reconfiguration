# X061 Johnson, Leeds & Wu 2015 II "Capability, Credibility, and Extended General Deterrence".
# Directed dyad-years (challenger -> target with defense pact); DV dispute = challenger initiates MID against target.
# Headline = Table 1 hypothesis-specific probit models (authors compute substantive effects from these, fn 23-24):
#   H1: probit dispute allytarg_win pchaloff pchalneu ln_distance jdem s_un_glo peaceyrs peaceyrs2 peaceyrs3   (capability)
#   H2: same with allies_s (credibility);  H4: same with milinst (peacetime military coordination). Default (iid) SEs.
# Key terms fixed from abstract before estimating: allytarg_win, allies_s, milinst. CONTROL: jdem (joint democracy) from H1 model.
# b_pub from the authors' posted log JohnsonLeedsWuIIrep.smcl.
ROOT <- Sys.getenv("SWEEP_DATA"); OUT <- Sys.getenv("SWEEP_OUT"); CODE <- Sys.getenv("SWEEP_CODE")
stopifnot(nzchar(ROOT), nzchar(OUT), nzchar(CODE))
source(file.path(CODE, "work_p4/P4_common.R"))
WO <- file.path(OUT, "work_p4/X061"); dir.create(WO, recursive = TRUE, showWarnings = FALSE)  # port: outputs
W <- file.path(ROOT, "work_p4/X061"); sink(file.path(WO, "log_X061.txt"), split = TRUE)
V <- c("challenger","target","year","dispute","allytarg_win","allies_s","milinst","pchaloff","pchalneu","ln_distance","jdem","s_un_glo","peaceyrs","peaceyrs2","peaceyrs3")
slim <- file.path(WO, "slim.rds")  # port: derived cache -> OUT
if (!file.exists(slim)) saveRDS(as.data.frame(haven::zap_labels(haven::read_dta(file.path(W, "JohnsonLeedsWuIIrep.dta"), col_select = tidyselect::all_of(V)))), slim)
d <- readRDS(slim)
stopifnot(!anyDuplicated(d[, c("challenger","target","year")]), all(d$dispute %in% c(0, 1, NA)))
cat("years:", range(d$year), "\n")
PR <- binomial("probit")
ctl <- "pchaloff + pchalneu + ln_distance + jdem + s_un_glo + peaceyrs + peaceyrs2 + peaceyrs3"
spec <- list(H1 = list(k = "allytarg_win", pub = -.4901084, n = 585432),
             H2 = list(k = "allies_s",     pub = -.7307049, n = 585432),
             H4 = list(k = "milinst",      pub = -.0365115, n = 572178))
jdem_pub <- -.0370351  # H1 model, p=.306
sc <- "Johnson, Leeds & Wu 2015 II"; est <- "probit, default (iid) SEs as in authors' log (fixest feglm probit)"
hc <- "Y (peace years + squared + cubed only)"
ev <- c(allytarg_win = "abstract: 'defense pacts with more capability and more credibility reduce the probability' of being targeted",
        allies_s = "abstract: 'defense pacts with more capability and more credibility reduce the probability' of being targeted",
        milinst = "abstract: 'pacts that include higher levels of peacetime military coordination are less likely to be attacked'")
for (h in names(spec)) {
  s <- spec[[h]]; K <- if (h == "H1") c(s$k, "jdem") else s$k
  f0 <- as.formula(paste("dispute ~", s$k, "+", ctl))
  r <- fit_keys(f0, d, K, "iid", file.path(WO, paste0("cache_", h, "_repro.rds")), PR)
  cat(h, "N repro", r$n, "vs pub", s$n, "\n")
  ds <- common_sample(d, f0, character(0), "challenger", "target", "year")
  M <- run_six(f0, ds, K, "iid", WO, paste0(h, "_"), PR); rm(ds); gc()
  # ---- sensitivity hook: runs only under run_sensitivity.R (SWEEP_SENS=1), never in the main sweep ----
  if (Sys.getenv("SWEEP_SENS") == "1" && h == "H2") {   # headline term only (allies_s; the loop also runs H1 and H4)
    source(file.path(CODE, "hist/sens_helper.R"))
    # ds is removed on the run_six() line above: rebuild the same common sample with the script's own call (same d, same f0)
    sens_run(id = "X061", dat = common_sample(d, f0, character(0), "challenger", "target", "year"), y = "dispute", tie = "allies_s",
             ctrl = c("pchaloff", "pchalneu", "ln_distance", "jdem", "s_un_glo", "peaceyrs", "peaceyrs2", "peaceyrs3"), fe = NULL,
             fit = function(dat, add) { r <- fit_keys(if (length(add)) update(f0, as.formula(paste(". ~ . +", paste(add, collapse = " + ")))) else f0, dat, "allies_s", "iid", NULL, PR)$k$allies_s
               c(b = unname(r[["b"]]), se = unname(r[["se"]]), p = unname(r[["p"]])) })
  }
  nt <- paste0("Table 1 ", h, " model (", r$n, " directed dyad-years vs pub ", s$n, "; sample = targets with >=1 defense pact). b_pub from authors' log (p<.05). ",
               "Common sample drops ", r$n - M$m0$n, " obs (years <1822 / benchmark missingness). History is undirected-dyad and monadic MID onsets.")
  write_p4("X061", sc, "TIE", "KEY", s$k, ev[[s$k]], repro_status(s$pub, TRUE, r$k[[s$k]]), s$pub, r$k[[s$k]][["b"]], est, r$n, M, hc, nt)
  if (h == "H1") write_p4("X061", sc, "DEMOCRACY", "CONTROL", "jdem", "control in all Table 1 models (both polity2>=6)",
      repro_status(jdem_pub, FALSE, r$k$jdem), jdem_pub, r$k$jdem[["b"]], est, r$n, M, hc, paste(nt, "jdem row from H1 model; bD NOT a placebo for this term."))
}
sink()
