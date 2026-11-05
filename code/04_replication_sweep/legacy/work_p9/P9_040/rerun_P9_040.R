# P9_040 Renshon 2016 IO 70(3) "Status Deficits and War" (10.1017/S0020818316000163)
# Data + do-files: author's replication zip (Dropbox link on jonathanrenshon.com/publications, "Status Deficits and War Replication.zip");
#   staged SupportingFiles/Renshon2015StatusAndWar.dta as renshon_data.dta. Do-file: SupportingFiles/StatusDeficitsandWar/StatusDeficitsandWar.do.
# Directed dyads 1816-2005. DV mzinit (MID initiation by cc1 vs cc2, contemporaneous: peaceyears resets at t+1 after an event) /
#   init5 (war initiation). IV L1std_comm_fg_pr_dissat_A = status deficit within the detected community, lagged 1 year.
# Estimator (do-file): relogit dv iv controls, cluster(dyadid) pc(sample share of 1s). pc = sample share -> prior correction is a no-op
#   on the intercept; only the King-Zeng finite-sample bias correction and (n/(n+k))^2 variance scaling apply. The archive's data have
#   no "dyadid" column; "dyad" is the directed-dyad id (one ordered pair per id) and is used for clustering.
# HEADLINE: Appendix C Table 2 Model 2 (MID initiation, controls; the model behind the article's Figure 4A first differences).
#   KEY TERM L1std_comm_fg_pr_dissat_A. Published (appendix table = Table5.tex in the archive): .208 (.0357), N = 1,130,777.
# SECOND ROW: Model 5 (war initiation, same controls): .182 (.0867), N = 1,130,777.
# Estimator used: King-Zeng relogit (work/relogit_kz_ea1.R), vcov clustered dyad, all six steps. Each spec gets its own common sample.
CODE <- Sys.getenv("SWEEP_CODE"); stopifnot(nzchar(CODE))
source(file.path(CODE, "work_p8/common_b3.R"))
source(file.path(CODE, "work/relogit_kz_ea1.R"))
WD <- file.path(ROOT, "work_p9/P9_040"); W <- file.path(OUT, "work_p9/P9_040"); dir.create(W, recursive = TRUE, showWarnings = FALSE)
sink(file.path(W, "log_P9_040.txt"), split = TRUE)
IV <- "L1std_comm_fg_pr_dissat_A"
ctl <- c("capshareA", "contig", "distance", "ongoing", "polity21", "jointdemocracy", "defense", "offense", "neutral", "nonagg", "consul")
d <- as.data.frame(haven::zap_labels(haven::read_dta(file.path(WD, "renshon_data.dta"),
       col_select = all_of(c("ccode1", "ccode2", "year", "dyad", "mzinit", "init5", "peaceyears", IV, ctl)))))
stopifnot(nrow(d) == 1498924, !anyDuplicated(d[, c("ccode1", "ccode2", "year")]))
# dyad is a directed-dyad id: each id maps to exactly one (ccode1, ccode2) pair
stopifnot(all(tapply(paste(d$ccode1, d$ccode2), d$dyad, function(z) length(unique(z))) == 1))
# DV timing: contemporaneous (after an mzinit event in t, peaceyears is 0 at t+1)
o <- d[order(d$ccode1, d$ccode2, d$year), ]; no <- nrow(o)
cs <- o$ccode1[-1] == o$ccode1[-no] & o$ccode2[-1] == o$ccode2[-no] & o$year[-1] == o$year[-no] + 1
ev <- which(cs & o$mzinit[-no] %in% 1)
cat("P(peaceyears[t+1] == 0 | mzinit[t] == 1):", mean(o$peaceyears[ev + 1] == 0, na.rm = TRUE), "\n")
stopifnot(mean(o$peaceyears[ev + 1] == 0, na.rm = TRUE) > 0.99); rm(o, cs, ev, no)
cat("year range:", range(d$year), "\n")
rhs0 <- c(IV, ctl)
specs <- list(
  M2_mid = list(dv = "mzinit", pub = c(b = 0.208, se = 0.0357, n = 1130777)),
  M5_war = list(dv = "init5",  pub = c(b = 0.182, se = 0.0867, n = 1130777)))
out <- list()
for (s in names(specs)) {
  sp <- specs[[s]]; vars <- c(sp$dv, rhs0, "dyad")
  fitter <- function(r, data) {
    o <- relogit_kz(as.formula(paste(sp$dv, "~", paste(r, collapse = " + "))), data, cluster = "dyad")
    cat("\n---- relogit (KZ-corrected, dyad-clustered)", sp$dv, "N =", attr(o, "n"), "\n"); print(round(o[c(IV, paste0(IV, "__unc")), ], 6))
    res <- list(o[IV, c("b", "se", "p")]); names(res) <- IV; attr(res, "n") <- attr(o, "n"); res
  }
  cat("\n==== spec", s, "\n")
  d0 <- d[complete.cases(d[, vars]), ]
  full <- fitter(rhs0, d0)
  cat("\nREPRO vs published:", IV, "b", full[[IV]][["b"]], "(pub", sp$pub[["b"]], ") se", full[[IV]][["se"]],
      "(pub", sp$pub[["se"]], ") N", attr(full, "n"), "(pub", sp$pub[["n"]], ")\n")
  # published b rounded to 3 decimals: require |b - pub| <= 0.0005 and the exact N
  stopifnot(attr(full, "n") == sp$pub[["n"]], abs(full[[IV]][["b"]] - sp$pub[["b"]]) <= 0.0005)
  ds <- common_sample(d, vars, "ccode1", "ccode2", "year")
  out[[s]] <- list(M = run_six(rhs0, ds, fitter), full = full, n_used = nrow(ds))

  # ---- sensitivity hook: runs only under run_sensitivity.R (SWEEP_SENS=1), never in the main sweep ----
  if (Sys.getenv("SWEEP_SENS") == "1" && s == "M2_mid") {
    source(file.path(CODE, "hist/sens_helper.R"))
    # headline: status deficit enters additively -> full sensitivity
    sens_run(id = "P9_040", dat = ds, y = "mzinit", tie = IV, ctrl = ctl, fe = NULL,
             fit = function(dat, add) fitter(c(rhs0, add), dat)[[IV]])
  }
}
saveRDS(out, file.path(W, "M_P9_040.rds"))

meta <- list(pid = "P9_040", cite = "Renshon 2016", journal = "IO", year = 2016, doi = "10.1017/S0020818316000163", family = "TIE")
art <- "rare-events logit (relogit, pc = sample share), SE clustered on directed dyad, no temporal controls"
used <- "King-Zeng relogit (work/relogit_kz_ea1.R), vcov clustered dyad, all six steps"
ahc <- "N (no temporal-dependence controls in the headline model)"
nt <- function(s, extra) { o <- out[[s]]; sp <- specs[[s]]
  sprintf(paste0("Directed dyads 1816-2005 (author's data). Published values from online appendix C Table 2 (3 decimals). ",
    "Full-sample repro b=%.6f se=%.6f N=%d vs pub %s (%s) N=%d. Common 1822-2011 sample N=%d. Clustered on 'dyad' (do-file names dyadid, absent from data). %s"),
    o$full[[IV]][["b"]], o$full[[IV]][["se"]], attr(o$full, "n"), format(sp$pub[["b"]]), format(sp$pub[["se"]]),
    as.integer(sp$pub[["n"]]), o$n_used, extra) }
dir.create(file.path(OUT, "results_p9"), recursive = TRUE, showWarnings = FALSE)
mv <- function(dst) { src <- file.path(OUT, "results_p8", paste0("P9_040__", IV, ".csv"))
  stopifnot(file.exists(src), file.rename(src, file.path(OUT, "results_p9", dst))) }
ev <- "abstract: 'status deficits are significantly associated with an increased probability of war and MID initiation'"
write_p8(meta, IV, ev, "Appendix C Table 2 Model 2 (MID initiation, controls)", art, used, 0.208, 0.0357, 1130777, TRUE, FALSE,
  out$M2_mid$M, ahc, nt("M2_mid", "Headline (MID initiation)."))
mv(paste0("P9_040__", IV, ".csv"))
# second row: same term, war-initiation DV -> distinct output file (suffix __war)
write_p8(meta, IV, ev, "Appendix C Table 2 Model 5 (war initiation, controls)", art, used, 0.182, 0.0867, 1130777, TRUE, FALSE,
  out$M5_war$M, ahc, nt("M5_war", "Non-headline second row (war initiation DV; output file suffix __war)."))
mv(paste0("P9_040__", IV, "__war.csv"))
sink()
