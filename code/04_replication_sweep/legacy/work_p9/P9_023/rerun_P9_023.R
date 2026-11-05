# P9_023 Krustev 2006 JPR "Interdependence and the Duration of Militarized Conflict" (10.1177/0022343306063930)
# Data + do-file + log: PRIO JPR replication archive "Replication Krustev.zip" (duration.dta, analysis_no_graphs.do, analysis.log).
# HEADLINE: Table I.2 (do-file "Table I.2", described there as "the main specification"):
#   stset calcdur if marker==1, failure(censored==.)
#   stcox ldistance lev4cont powdisparity allies jointdem kdeplo, tvc(allies jointdem) texp(ln(_t)) efron nohr
#   Unit = dyadic MID (first-year record, marker==1), 1950-92; outcome = MID termination (hazard), model-based (OIM) SE, no cluster.
#   Published (author's analysis.log): kdeplo 9.654224 (3.604127), N = 1,450 subjects, 1,439 failures.
# Non-headline row: Table I.1 (same model without the tvc terms): kdeplo 10.12026 (3.543817), N = 1,450.
# KEY TERM (a priori, abstract): "inverse relationship between interdependence and conflict duration" -> kdeplo
#   (lower dyadic trade/GDP ratio); positive b = higher hazard of termination = shorter MIDs.
# Estimator (precedent work/100/rerun_100.R, work_p4/X016/rerun_verify_X016.R): survival::coxph, same tie method as the author
#   (Efron), model-based SE (Stata stcox default). tvc(x) texp(ln(_t)) -> coxph tt() term x * log(t), all six steps.
#   History terms enter the rh (time-fixed) part additively. History keyed to the MID start year (year == strtyr for marker==1).
# No randomness (seed 123 set in common_b3.R).
CODE <- Sys.getenv("SWEEP_CODE"); stopifnot(nzchar(CODE))
source(file.path(CODE, "work_p8/common_b3.R"))
suppressPackageStartupMessages(library(survival))
WD <- file.path(ROOT, "work_p9/P9_023"); W <- file.path(OUT, "work_p9/P9_023"); dir.create(W, recursive = TRUE, showWarnings = FALSE)
sink(file.path(W, "log_P9_023.txt"), split = TRUE)
# raw archive is a zip: idempotent unzip of the one member into OUT (convention 8)
unzip(file.path(ROOT, "work_p9/P9_023/krustev.zip"), files = "duration.dta", exdir = W, overwrite = TRUE)
d <- as.data.frame(haven::zap_labels(haven::read_dta(file.path(W, "duration.dta"))))
stopifnot(nrow(d) == 2205)
d <- d[d$marker == 1, ]                                   # stset ... if marker==1
stopifnot(nrow(d) == 1524, !anyDuplicated(d$caseid), all(d$year == d$strtyr), all(d$calcdur > 0))
d$fail <- as.numeric(is.na(d$censored))                   # failure(censored==.)
stopifnot(sum(d$fail) == 1509)

ctrl <- c("ldistance", "lev4cont", "powdisparity", "allies", "jointdem")
TVC <- c("allies", "jointdem")
fit_cox <- function(term, tvc) function(rhs, data) {
  rhs_f <- if (tvc) c(rhs, paste0("tt(", TVC, ")")) else rhs
  f <- as.formula(paste("Surv(calcdur, fail) ~", paste(rhs_f, collapse = " + ")))
  m <- if (tvc) coxph(f, data = data, ties = "efron", tt = rep(list(function(x, t, ...) x * log(t)), length(TVC)))
       else coxph(f, data = data, ties = "efron")
  stopifnot(m$iter < 30, all(is.finite(coef(m))))
  ct <- summary(m)$coefficients
  cat("\n---- coxph (Efron)", deparse(f, width.cutoff = 500), if (tvc) "| tt = x*log(t)" else "", "| N =", m$n, "| failures =", m$nevent,
      "| logLik =", round(m$loglik[2], 4), "\n")
  print(signif(ct[, c("coef", "se(coef)", "z", "Pr(>|z|)")], 7))
  stopifnot("kdeplo" %in% rownames(ct))
  out <- setNames(list(c(b = ct["kdeplo", "coef"], se = ct["kdeplo", "se(coef)"], p = ct["kdeplo", "Pr(>|z|)"])), term)
  attr(out, "n") <- m$n; out
}
SPEC <- list(kdeplo    = list(tvc = TRUE,  pub = c(b = 9.654224, se = 3.604127, n = 1450), ll = -9012.0944),
             kdeplo_I1 = list(tvc = FALSE, pub = c(b = 10.12026, se = 3.543817, n = 1450), ll = -9019.0167))
rhs0 <- c(ctrl, "kdeplo"); vars <- c("calcdur", "fail", rhs0, "statea", "stateb", "year")
full <- d[complete.cases(d[, vars]), ]
R <- list()
for (term in names(SPEC)) {
  fitter <- fit_cox(term, SPEC[[term]]$tvc)
  cat("\n==== term", term, "(", if (SPEC[[term]]$tvc) "Table I.2" else "Table I.1", ") complete-case N:", nrow(full), "\n")
  r_full <- fitter(rhs0, full); n_orig <- attr(r_full, "n")
  pb <- SPEC[[term]]$pub
  cat("REPRO vs author log:", term, "b", r_full[[term]][["b"]], "(pub", pb[["b"]], ") se", r_full[[term]][["se"]], "(pub", pb[["se"]],
      ") N", n_orig, "(pub", pb[["n"]], ") rel.diff b", abs(r_full[[term]][["b"]] - pb[["b"]]) / abs(pb[["b"]]), "\n")
  stopifnot(n_orig == pb[["n"]], abs(r_full[[term]][["b"]] - pb[["b"]]) / abs(pb[["b"]]) <= 0.01)
  ds <- common_sample(d, vars, "statea", "stateb", "year")
  M <- run_six(rhs0, ds, fitter)
  R[[term]] <- list(M = M, r_full = r_full, n_orig = n_orig, n_used = nrow(ds))
}
saveRDS(R, file.path(W, "M_P9_023.rds"))
# ---- result-CSV write step (write_p8 writes to results_p8; the file is moved to results_p9 unchanged) ----
meta <- list(pid = "P9_023", cite = "Krustev 2006", journal = "JPR", year = 2006, doi = "10.1177/0022343306063930", family = "TIE")
ev <- "abstract: 'inverse relationship between interdependence and conflict duration'"
art <- "Cox PH (stcox, Efron ties), OIM SE; Table I.2 adds tvc(allies jointdem) x ln(t); dyadic MIDs 1950-92"
est <- "survival::coxph, Efron ties, model-based SE; tvc via tt() = x*log(t) (Table I.2 only), all six steps"
ahc <- "N"
mv <- function(term) {
  src <- file.path(OUT, "results_p8", paste0("P9_023__", term, ".csv")); dst <- file.path(OUT, "results_p9", paste0("P9_023__", term, ".csv"))
  dir.create(dirname(dst), recursive = TRUE, showWarnings = FALSE); stopifnot(file.exists(src), file.rename(src, dst))
}
nt <- function(term) sprintf(paste0("Dyadic MID duration (Cox, outcome = termination hazard; positive b = shorter MIDs). History keyed to MID start year. ",
  "Authors' PRIO data + do-file; published values from author's analysis.log. Full-sample repro b=%.4f se=%.4f N=%d. Common 1822-2011 sample N=%d."),
  R[[term]]$r_full[[term]][["b"]], R[[term]]$r_full[[term]][["se"]], R[[term]]$n_orig, R[[term]]$n_used)
write_p8(meta, "kdeplo", ev, "Table I.2 (Cox PH with tvc(allies jointdem) x ln(t); main specification)", art, est,
  9.654224, 3.604127, 1450, TRUE, FALSE, R$kdeplo$M, ahc, nt("kdeplo")); mv("kdeplo")
write_p8(meta, "kdeplo_I1", ev, "Table I.1 (Cox PH, time-fixed covariates)", art, est,
  10.12026, 3.543817, 1450, TRUE, FALSE, R$kdeplo_I1$M, ahc, nt("kdeplo_I1")); mv("kdeplo_I1")
sink()
