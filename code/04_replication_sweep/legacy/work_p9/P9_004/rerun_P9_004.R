# P9_004 Kang 2017 CMPS 34(5) "Capability Revisited: Ally's Capability and Dispute Initiation" (10.1177/0738894215604966)
# Data + log: Harvard Dataverse doi:10.7910/DVN/AMOREX ("Alliance Capability & MID Initiation Replication Data.dta",
#   "Kang CMPS2015 Replication STATA Log File.smcl", "Replication Commands.txt"); staged as kang_data.dta / kang_log.smcl / kang_commands.txt.
# Directed dyads, DV cwinit = cc1 initiates a MID vs cc2. All models: logit, vce(cluster ddyad), peaceyrs + square + cube.
# HEADLINE (PI decision; abstract's first claim, allies' capability raises initiation against non-members):
#   Table 1 Model 3 = "logit cwinit c1allcap c2allcap alleach kgdriv relcap jodem sreg condum peaceyrs peaceyrsq peaceyrcube".
#   KEY TERM c1allcap (strongest ally's CINC for cc1). Author's log: 1.481363 (robust SE .5355168), N = 164,391.
# SECOND ROW (non-headline; abstract's second claim, effect "strengthens as the allies' common interests grow"):
#   Table 3 = Table 1 Model 3 + capinter + c1alls2. KEY TERM capinter = c1allcap x c1alls2. Log: 2.602303 (1.210935), N = 164,390.
# Estimator: fixest feglm logit, vcov ~ddyad, all six steps. Each spec gets its own common sample.
CODE <- Sys.getenv("SWEEP_CODE"); stopifnot(nzchar(CODE))
source(file.path(CODE, "work_p8/common_b3.R"))
WD <- file.path(ROOT, "work_p9/P9_004"); W <- file.path(OUT, "work_p9/P9_004"); dir.create(W, recursive = TRUE, showWarnings = FALSE)
sink(file.path(W, "log_P9_004.txt"), split = TRUE)
d <- as.data.frame(haven::zap_labels(haven::read_dta(file.path(WD, "kang_data.dta"))))
stopifnot(nrow(d) == 167811, !anyDuplicated(d[, c("ccode1", "ccode2", "year")]))
stopifnot(max(abs(d$capinter - d$c1allcap * d$c1alls2), na.rm = TRUE) < 1e-6)
cat("year range:", range(d$year), "\n")
ctl <- c("c2allcap", "alleach", "kgdriv", "relcap", "jodem", "sreg", "condum", "peaceyrs", "peaceyrsq", "peaceyrcube")
specs <- list(
  T1M3 = list(term = "c1allcap", rhs = c("c1allcap", ctl), pub = c(b = 1.481363, se = 0.5355168, n = 164391)),
  T3   = list(term = "capinter", rhs = c("c1allcap", "capinter", "c1alls2", ctl), pub = c(b = 2.602303, se = 1.210935, n = 164390)))
out <- list()
for (s in names(specs)) {
  sp <- specs[[s]]; rhs0 <- sp$rhs; vars <- c("cwinit", rhs0, "ddyad")
  fitter <- function(r, data) fit_bin("cwinit", r, data, "logit", "~ddyad", terms = sp$term)
  cat("\n==== spec", s, "\n")
  d0 <- d[complete.cases(d[, vars]), ]
  full <- fitter(rhs0, d0)
  cat("\nREPRO vs log:", sp$term, "b", full[[sp$term]][["b"]], "(log", sp$pub[["b"]], ") se", full[[sp$term]][["se"]],
      "(log", sp$pub[["se"]], ") N", attr(full, "n"), "(log", sp$pub[["n"]], ")\n")
  stopifnot(attr(full, "n") == sp$pub[["n"]], abs(full[[sp$term]][["b"]] - sp$pub[["b"]]) / abs(sp$pub[["b"]]) <= 0.01)
  ds <- common_sample(d, vars, "ccode1", "ccode2", "year")
  out[[s]] <- list(M = run_six(rhs0, ds, fitter), full = full, n_used = nrow(ds))

  # ---- sensitivity hook: runs only under run_sensitivity.R (SWEEP_SENS=1), never in the main sweep ----
  if (Sys.getenv("SWEEP_SENS") == "1") {
    source(file.path(CODE, "hist/sens_helper.R"))
    if (s == "T1M3")   # headline: c1allcap enters additively -> full sensitivity
      sens_run(id = "P9_004", dat = ds, y = "cwinit", tie = "c1allcap", ctrl = setdiff(rhs0, "c1allcap"), fe = NULL,
               fit = function(dat, add) fitter(c(rhs0, add), dat)[["c1allcap"]])
  }
}
saveRDS(out, file.path(W, "M_P9_004.rds"))

meta <- list(pid = "P9_004", cite = "Kang 2017", journal = "CMPS", year = 2017, doi = "10.1177/0738894215604966", family = "TIE")
art <- "logit, SE clustered on directed dyad, peace years + square + cube"
used <- "feglm logit, vcov clustered ddyad, all six steps"
ahc <- "N (peace years + square + cube only)"
nt <- function(s, extra) { o <- out[[s]]; sp <- specs[[s]]
  sprintf(paste0("Directed dyads (Kang data, %d-%d). Published values from the author's Stata log. Full-sample repro b=%.6f se=%.6f N=%d ",
    "vs log %s (%s) N=%d. Common 1822-2011 sample N=%d. %s"), min(d$year), max(d$year), o$full[[sp$term]][["b"]],
    o$full[[sp$term]][["se"]], attr(o$full, "n"), format(sp$pub[["b"]]), format(sp$pub[["se"]]), as.integer(sp$pub[["n"]]), o$n_used, extra) }
write_p8(meta, "c1allcap",
  "abstract: 'the capability of allies increases' members' likelihood of initiating disputes",
  "Table 1 Model 3 (ally capability, full controls; author's log)", art, used,
  1.481363, 0.5355168, 164391, TRUE, FALSE, out$T1M3$M, ahc, nt("T1M3", "Headline (PI decision)."))
write_p8(meta, "capinter",
  "abstract: ally capability effect 'strengthens as the allies' common interests grow'",
  "Table 3 (ally capability x strongest ally's S-score; author's log)", art, used,
  2.602303, 1.210935, 164390, TRUE, FALSE, out$T3$M, ahc,
  nt("T3", "Non-headline second row. Interaction: constituent terms c1allcap and c1alls2 kept in all steps."))
# phase-9 output location: move the write_p8 row files into results_p9 (helper is shared and writes to results_p8)
dir.create(file.path(OUT, "results_p9"), recursive = TRUE, showWarnings = FALSE)
for (tm in c("c1allcap", "capinter")) {
  src <- file.path(OUT, "results_p8", paste0("P9_004__", tm, ".csv")); dst <- file.path(OUT, "results_p9", paste0("P9_004__", tm, ".csv"))
  stopifnot(file.exists(src), file.rename(src, dst))
}
sink()
