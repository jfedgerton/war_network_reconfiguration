# P9_057 Fang, Johnson & Leeds 2014 IO 68(4) "To Concede or to Resist? The Restraining Effect of Military Alliances" (10.1017/S0020818314000137)
# Data + do-file + log: Harvard Dataverse doi:10.7910/DVN/S0020818314000137 (FangJohnson&LeedsIOrep.zip; staged as fjl.dta / fjl.do /
#   fjl.smcl / fjl_appendix.pdf). Earlier stop note (P9_009): selection estimator; now portable via work_p8/heckprob_b3.R (as P8_003).
# HEADLINE: Table 2 = heckprob targ_resist ally_cost pchaloff pchalneu c_win,
#   select(dispute = pchaloff pchalneu c_win ln_distance jdem s_un_glo peaceyrs peaceyrs2 peaceyrs3)   [OIM SEs, no vce option]
#   Directed challenger-target dyad-years 1816-2000; selection = challenger initiates a MID vs target; outcome = target resists.
# KEY TERM out:ally_cost (target's cost of damaging its alliances by resisting) in the OUTCOME equation.
#   NOTE: ally_cost is a property of the TARGET's alliance portfolio, not a challenger-target tie.
#   Author's log: ally_cost -.4043512 (.1365033), N = 585,467 (1,085 uncensored), athrho -.578801.
# Stata marks out outcome-equation regressors only on selected rows (ally_cost is defined only where dispute == 1); censored rows
#   contribute Phi(-Zg) only, so their X is unused and is set to 0 before fitting (no effect on the likelihood).
# Estimator: heckprob by ML (work_p8/heckprob_b3.R), OIM SE, all six steps; history/benchmark terms enter both equations (as P8_003).
CODE <- Sys.getenv("SWEEP_CODE"); stopifnot(nzchar(CODE))
source(file.path(CODE, "work_p8/common_b3.R"))
source(file.path(CODE, "work_p8/heckprob_b3.R"))
WD <- file.path(ROOT, "work_p9/P9_057"); W <- file.path(OUT, "work_p9/P9_057"); dir.create(W, recursive = TRUE, showWarnings = FALSE)
sink(file.path(W, "log_P9_057.txt"), split = TRUE)
outX <- c("ally_cost", "pchaloff", "pchalneu", "c_win")
selZ <- c("pchaloff", "pchalneu", "c_win", "ln_distance", "jdem", "s_un_glo", "peaceyrs", "peaceyrs2", "peaceyrs3")
d <- haven::read_dta(file.path(WD, "fjl.dta"), col_select = all_of(unique(c("challenger", "target", "year", "dispute", "targ_resist", outX, selZ))))
d <- as.data.frame(haven::zap_labels(d)); cat("rows", nrow(d), "years", range(d$year), "\n")
stopifnot(nrow(d) == 585467, !anyDuplicated(d[, c("challenger", "target", "year")]), all(d$dispute %in% 0:1))
stopifnot(all(is.na(d$ally_cost[d$dispute == 0])), all(!is.na(d$ally_cost[d$dispute == 1])), all(d$targ_resist[d$dispute == 1] %in% 0:1))
# outcome-side columns: zero on censored rows (unused there), checked complete on selected rows
d$ally_cost <- ifelse(d$dispute == 1, d$ally_cost, 0); d$targ_resist_chk <- ifelse(d$dispute == 1, d$targ_resist, 0)
KEY <- "out:ally_cost"
fit_hp <- function(extra, data) {
  X <- cbind(as.matrix(data[, c(outX, extra), drop = FALSE]), `_cons` = 1)
  Z <- cbind(as.matrix(data[, c(selZ, extra), drop = FALSE]), `_cons` = 1)
  s <- data$dispute; y <- data$targ_resist_chk
  f <- fit_heckprob(X, Z, y, s)
  z <- f$coef / f$se
  cat("---- heckprob | extra:", paste(extra, collapse = " + "), "| N =", f$n, "selected =", f$n_sel, "| LL =", round(f$ll, 3), "| rho =", round(f$rho, 4), "\n")
  print(round(cbind(b = f$coef, se = f$se, z = z), 6))
  out <- list(ally_cost = c(b = f$coef[[KEY]], se = f$se[[KEY]], p = 2 * pnorm(-abs(z[[KEY]]))))
  attr(out, "n") <- f$n; attr(out, "athrho") <- f$coef[["athrho"]]; out
}
# gradient check of the analytic heckprob scores on synthetic data (fail loudly if wrong; same check as P8_003)
set.seed(123); nn <- 4000; Xs0 <- cbind(a = rnorm(nn), c = 1); Zs0 <- cbind(a = rnorm(nn), b = rnorm(nn), c = 1)
ee <- MASS::mvrnorm(nn, c(0, 0), matrix(c(1, .4, .4, 1), 2)); ss <- as.numeric(Zs0 %*% c(.5, -.3, -.2) + ee[, 2] > 0)
yy <- as.numeric(Xs0 %*% c(.7, .1) + ee[, 1] > 0); yy[ss == 0] <- 0; th <- c(.6, .05, .45, -.25, -.15, .3)
ga <- colSums(hp_parts(th, Xs0, Zs0, yy, ss)$S); gn <- numDeriv::grad(function(t) sum(hp_parts(t, Xs0, Zs0, yy, ss)$ll), th)
cat("synthetic gradient check max abs diff:", max(abs(ga - gn)), "\n"); stopifnot(max(abs(ga - gn)) < 1e-5)
vars <- unique(c("dispute", "targ_resist_chk", outX, selZ))
full <- fit_hp(character(0), d[complete.cases(d[, vars]), ])
pub <- c(b = -0.4043512, se = 0.1365033); pub_n <- 585467
cat("\nREPRO vs log: ally_cost b", full$ally_cost[["b"]], "(log", pub[["b"]], ") se", full$ally_cost[["se"]], "(log", pub[["se"]], ") N", attr(full, "n"),
    "(log", pub_n, ") athrho", attr(full, "athrho"), "(log -.578801)\n")
stopifnot(attr(full, "n") == pub_n, abs(full$ally_cost[["b"]] - pub[["b"]]) / abs(pub[["b"]]) <= 0.01)
ds <- common_sample(d, vars, "challenger", "target", "year"); rm(d); invisible(gc())
cat("selected rows in common sample:", sum(ds$dispute), "\n")
fitter <- function(rhs, data) fit_hp(setdiff(rhs, vars), data)
M <- run_six(vars, ds, fitter)
saveRDS(list(M = M, full = full, n_used = nrow(ds)), file.path(W, "M_P9_057.rds"))

meta <- list(pid = "P9_057", cite = "Fang 2014", journal = "IO", year = 2014, doi = "10.1017/S0020818314000137", family = "TIE")
nt <- sprintf(paste0("Directed challenger-target dyad-years 1816-2000 (author's data). Published values from the author's Stata log (Table 2). ",
  "Full-sample repro b=%.6f se=%.6f N=%d vs log -0.4043512 (0.1365033) N=585467. Common 1822-2011 sample N=%d (%d selected). ",
  "OIM SEs as in the do-file. History terms in both equations. Key term is a property of the target's alliances (not a dyadic tie)."),
  full$ally_cost[["b"]], full$ally_cost[["se"]], attr(full, "n"), nrow(ds), sum(ds$dispute))
write_p8(meta, "ally_cost", "title: 'The Restraining Effect of Military Alliances' (targets' alliance costs restrain resistance)",
  "Table 2 (heckprob; outcome eq. target resists, selection eq. dispute; author's log)",
  "probit with sample selection (heckprob), OIM SE, peace years cubic in selection eq.",
  "heckprob ML (R port, analytic scores), OIM SE, all six steps",
  pub[["b"]], pub[["se"]], pub_n, TRUE, FALSE, M, "N (peace years cubic in selection eq. only)", nt)
dir.create(file.path(OUT, "results_p9"), recursive = TRUE, showWarnings = FALSE)
src <- file.path(OUT, "results_p8", "P9_057__ally_cost.csv"); dst <- file.path(OUT, "results_p9", "P9_057__ally_cost.csv")
stopifnot(file.exists(src), file.rename(src, dst))
sink()
