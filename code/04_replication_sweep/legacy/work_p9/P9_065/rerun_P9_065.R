# P9_065 Zeng 2021 JPR 58(4) "Biding time versus timely retreat: Asymmetric dependence, issue salience, and conflict duration"
#   (10.1177/0022343320918917). TIE family (relative economic cost of a crisis to each side = asymmetric trade dependence).
# Data + code: PRIO JPR replication archive Zeng.zip (https://files.prio.org/journals/JPR/2021/58/4/Zeng.zip ->
#   cdn.cloud.prio.org/files/65bf06f4-.../Zeng.zip); members used: JPR/Data/DurationTimeCov_ICB.RData (object durB),
#   JPR/Code/02_MainModel.R (author's R code; read as text, not executed). The article PDF was not available here.
# Unit: ICB crisis actor-opponent directed dyads 1917-2014 in counting-process form (tstart, tstop], DV quit = actor leaves the
#   crisis (termination hazard; dispute duration is in scope per PI 2026-10-10).
# Author's models (02_MainModel.R), reproduced line by line:
#   joint_demo = demo1 > 5 & demo2 > 5; costratio = econcost1 / (econcost2 + .Machine$double.eps); issuesalience = gravty in 3:6;
#   defenseratio = defense1 / (defense2 + defense1)
#   fit1bc ("original model"): coxph(Surv(tstart, tstop, quit) ~ costratio + issuesalience + costratio:issuesalience + joint_demo +
#            contbinary + powerratio + defenseratio + cluster(crisno))
#   cmod2 ("adjusted model", time interactions): fit1bc + tcontbinary + tcostissue, where t<x> = x * log(tstop + lamhat) and lamhat
#            minimises BIC over lam = 0:99 (author's grid).
# HEADLINE (article read 2026-10-11, Zeng (2021) JPR.pdf Table I): Model 1 (adjusted model, time interaction; BIC-selected, BIC 765.99),
#   term Cost x Salience (cost_x_sal, H1: cost ratio conditional on issue salience). Second row: Cost ratio (constituent), same model.
#   Table I reports b with 95% CIs; published SE = (upper - lower) / (2 * 1.96). N = 667. Legend: dagger p<0.1; * p<0.05; ** p<0.01.
# Estimator: survival::coxph (Efron ties, as the author), robust SE clustered on crisno (cluster()); lamhat chosen once on the baseline
#   specification and held fixed in the six steps. History keyed to the crisis-dyad start year (first record's year; precedent P9_023).
CODE <- Sys.getenv("SWEEP_CODE"); stopifnot(nzchar(CODE))
source(file.path(CODE, "work_p8/common_b3.R"))
suppressPackageStartupMessages(library(survival))
WD <- file.path(ROOT, "work_p9/P9_065"); W <- file.path(OUT, "work_p9/P9_065"); dir.create(W, recursive = TRUE, showWarnings = FALSE)
sink(file.path(W, "log_P9_065.txt"), split = TRUE)
e <- new.env(); load(file.path(WD, "DurationTimeCov_ICB.RData"), envir = e); stopifnot(identical(ls(e), "durB"))
d <- as.data.frame(e$durB); rm(e)
stopifnot(nrow(d) == 1597, all(d$tstop > d$tstart), all(d$quit %in% 0:1))
d$joint_demo <- ifelse(d$demo1 > 5 & d$demo2 > 5, 1, 0)
d$costratio <- d$econcost1 / (d$econcost2 + .Machine$double.eps)
d$issuesalience <- as.numeric(d$gravty %in% c(3, 4, 5, 6))
d$defenseratio <- d$defense1 / (d$defense2 + d$defense1)
d$cost_x_sal <- d$costratio * d$issuesalience
d$time <- d$tstop
d$start_year <- ave(d$year, d$crisno, d$ccode1, d$ccode2, FUN = min)
base <- c("costratio", "issuesalience", "cost_x_sal", "joint_demo", "contbinary", "powerratio", "defenseratio")
mk_t <- function(dat, lam) { dat$tcontbinary <- dat$contbinary * log(dat$time + lam); dat$tcostissue <- dat$cost_x_sal * log(dat$time + lam); dat }
cfit <- function(rhs, data) coxph(as.formula(paste("Surv(tstart, tstop, quit) ~", paste(rhs, collapse = " + "), "+ cluster(crisno)")),
                                 data = data, ties = "efron")
# author's lambda search (BIC) on the full data
bic <- sapply(0:99, function(l) BIC(cfit(c(base, "tcontbinary", "tcostissue"), mk_t(d, l))))
lamhat <- (0:99)[which.min(bic)]; cat("lamhat (BIC over 0:99):", lamhat, "\n")
d <- mk_t(d, lamhat)
rhs0 <- c(base, "tcontbinary", "tcostissue")
TERMS <- c("costratio", "cost_x_sal")
fitter <- function(r, data) {
  m <- cfit(r, data); ct <- summary(m)$coefficients
  cat("\n---- coxph (Efron, cluster crisno)", paste(r, collapse = " + "), "| N =", m$n, "| events =", m$nevent, "\n"); print(round(ct, 6))
  out <- lapply(setNames(TERMS, TERMS), function(t) c(b = ct[t, "coef"], se = ct[t, "robust se"], p = ct[t, "Pr(>|z|)"]))
  attr(out, "n") <- m$n; out
}
vars <- c("tstart", "tstop", "quit", rhs0, "crisno")
orig <- cfit(base, d); cat("\nORIGINAL MODEL (fit1bc):\n"); print(round(summary(orig)$coefficients, 6))
full <- fitter(rhs0, d[complete.cases(d[, vars]), ])
for (t in TERMS) cat("FULL-SAMPLE cmod2:", t, "b", full[[t]][["b"]], "se", full[[t]][["se"]], "N", attr(full, "n"), "\n")
stopifnot(attr(full, "n") == 667, abs(full$cost_x_sal[["b"]] - 3.008) < 0.001, abs(full$costratio[["b"]] - 0.236) < 0.001)  # Table I Model 1
ds <- common_sample(d, vars, "ccode1", "ccode2", "start_year")
M <- run_six(rhs0, ds, fitter)
saveRDS(list(M = M, full = full, lamhat = lamhat, n_used = nrow(ds)), file.path(W, "M_P9_065.rds"))

meta <- list(pid = "P9_065", cite = "Zeng 2021", journal = "JPR", year = 2021, doi = "10.1177/0022343320918917", family = "TIE")
art <- "Cox PH (Efron), robust SE clustered on crisis, time interactions x*log(t+lambda) for contiguity and cost x salience"
used <- "survival::coxph Efron, cluster(crisno), lambda fixed at the baseline BIC choice, all six steps"
dir.create(file.path(OUT, "results_p9"), recursive = TRUE, showWarnings = FALSE)
PUB <- list(cost_x_sal = list(b = 3.008, se = (5.726 - 0.290) / 3.92, sig = TRUE,
                ev = "Table I Model 1: Cost x Salience 3.008* (95% CI 0.290, 5.726); H1 (cost ratio conditional on salience); fn. 24"),
            costratio  = list(b = 0.236, se = (1.073 + 0.601) / 3.92, sig = FALSE,
                ev = "Table I Model 1: Cost ratio 0.236 (95% CI -0.601, 1.073); constituent of the headline product"))
for (t in TERMS) {
  write_p8(meta, t, PUB[[t]]$ev,
    "Table I Model 1 (adjusted Cox model with time interactions; ICB crisis actor-dyads; days before quitting)", art, used,
    PUB[[t]]$b, PUB[[t]]$se, 667, PUB[[t]]$sig, FALSE, M, "N (no prior-conflict term)",
    sprintf("ICB crisis actor-dyads 1917-2014, counting process. lamhat=%d. Full-sample b=%.6f se=%.6f N=%d (published N 667). Common sample N=%d (history at crisis-dyad start year). Published SE derived from the 95%% CI.",
            lamhat, full[[t]][["b"]], full[[t]][["se"]], attr(full, "n"), nrow(ds)))
  src <- file.path(OUT, "results_p8", paste0("P9_065__", t, ".csv")); stopifnot(file.exists(src), file.rename(src, file.path(OUT, "results_p9", basename(src))))
}
sink()
