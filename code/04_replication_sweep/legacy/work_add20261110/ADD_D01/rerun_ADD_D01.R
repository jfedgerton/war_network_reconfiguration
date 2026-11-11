# ADD_D01 Lektzian & Souva 2009 JPR "A Comparative Theory Test of Democratic Peace Arguments, 1946-2000" (10.1177/0022343308098402).
# PORT of the authors' PRIO/JPR replication archive ("Lektzian Souva JPR Democratic Peace Replication Files.zip": do-file + 4 .dta).
# ARTICLE PDF NOT ACCESSIBLE (SAGE paywall; not in literature/pdfs_for_sweep) -> no printed coefficients to check against. As for
#   P8_003, the published value is our run of the authors' own do-file specification on their own file (published_check = unverified).
# Design (do-file, Table 3 Model 3A): heckprob, outcome = reciprocation (recipt1) | selection = MID onset (midjt1), no joiners,
#   nondirected dyad-years 1946-2000, cluster(dyadid). midjt1/recipt1 are measured in year t+1 (checked below against the GML MID
#   onsets in the history tables), so history is keyed to the OUTCOME year (hy = year + 1), as in work_p8/p8_common.R.
# HEADLINE (task: joint regime categories -> MID onset; onset equation of the censored probit): selection-equation jointd6 (joint
#   democracy, both Polity >= 6) in Model 3A. Joint autocratic subtypes (monarchy, personal, socialist, military) are kept as
#   non-headline rows. Model 3B (fatal escalation) has the identical selection specification.
# Estimator: heckprob by ML (work_p8/heckprob_b3.R: analytic scores), cluster-robust on dyadid, all six steps; history/benchmark
#   terms enter both equations (P8_003 convention).
# Author history control: peace years + 3 splines only in the onset equation -> B_duration.
ROOT <- Sys.getenv("SWEEP_DATA"); OUT <- Sys.getenv("SWEEP_OUT"); CODE <- Sys.getenv("SWEEP_CODE")
stopifnot(nzchar(ROOT), nzchar(OUT), nzchar(CODE))
source(file.path(CODE, "work_p8/common_b3.R"))     # dplyr, fixest, add_history, common_sample, run_six, write_p8; set.seed(123)
source(file.path(CODE, "work_p8/heckprob_b3.R"))   # fit_heckprob()
set.seed(123)
suppressMessages(library(haven))
W <- file.path(OUT, "work_add20261110/ADD_D01"); dir.create(W, recursive = TRUE, showWarnings = FALSE)
sink(file.path(W, "log_ADD_D01.txt"), split = TRUE)

# ---- 1. authors' file 1 (MIDs, no joiners) ----
outX <- c("monarchy", "personal", "socialist", "military", "jointd6", "powerratio", "lndist", "majdyad", "allies", "tradeweak", "devweak")
selZ <- c("monarchy", "personal", "socialist", "military", "jointd6", "powerratio", "majdyad", "tradeweak", "devweak", "lndist", "allies",
          "midjpyr", "midspline1", "midspline2", "midspline3")
d <- read_dta(file.path(ROOT, "_incoming/Lektzian2009/JPR replication file 1 no joiners mids.dta"),
              col_select = all_of(unique(c("statea", "stateb", "year", "dyadid", "midjt1", "recipt1", outX, selZ)))) %>% zap_labels() %>% as.data.frame()
stopifnot(nrow(d) == 536522, anyDuplicated(d[, c("statea", "stateb", "year")]) == 0, all(d$statea < d$stateb),
          all(d$midjt1 %in% c(0, 1, NA)), all(d$recipt1[d$midjt1 %in% 1] %in% 0:1))
# timing check: midjt1 = 1 should coincide with a GML onset in year + 1 far more often than in year (onset(y) ~ d15(y+1) - d26(y+1))
k <- function(yy) { x <- data.frame(c1 = d$statea, c2 = d$stateb, year = yy) %>% left_join(.H$dyh, by = c("c1", "c2", "year"))
  coalesce(x$d15, 0) - coalesce(x$d26, 0) }
i1 <- which(d$midjt1 == 1); on_t <- mean(k(d$year + 1)[i1] > 0); on_t1 <- mean(k(d$year + 2)[i1] > 0)
cat("share of midjt1=1 with a GML onset in year t:", round(on_t, 3), " in t+1:", round(on_t1, 3), "\n"); stopifnot(on_t1 > 0.5, on_t1 > 2 * on_t)
d$hy <- d$year + 1

KEYS <- c("sel:jointd6", "sel:monarchy", "sel:personal", "sel:socialist", "sel:military")
fit_hp <- function(extra, data) {
  X <- cbind(as.matrix(data[, c(outX, extra), drop = FALSE]), `_cons` = 1)
  Z <- cbind(as.matrix(data[, c(selZ, extra), drop = FALSE]), `_cons` = 1)
  s <- data$midjt1; y <- ifelse(s == 1, data$recipt1, 0)
  f <- fit_heckprob(X, Z, y, s, cl = data$dyadid)
  z <- f$coef / f$se
  cat("---- heckprob | extra:", paste(extra, collapse = " + "), "| N =", f$n, "selected =", f$n_sel, "| LL =", round(f$ll, 3), "| rho =", round(f$rho, 4), "\n")
  print(round(cbind(b = f$coef, se = f$se, z = z), 5))
  out <- lapply(setNames(KEYS, KEYS), function(k) c(b = f$coef[[k]], se = f$se[[k]], p = 2 * pnorm(-abs(z[[k]]))))
  attr(out, "n") <- f$n; out
}
vars <- unique(c("midjt1", outX, selZ, "dyadid"))
d$recipt1_chk <- ifelse(d$midjt1 == 1, d$recipt1, 0)

# ---- 2. authors' specification on the full sample (Stata heckprob sample: complete on both equations' regressors and midjt1) ----
full <- d[complete.cases(d[, c(vars, "recipt1_chk")]), ]
cat("full-sample N:", nrow(full), " selected (MID onsets):", sum(full$midjt1), "\n")
r_full <- fit_hp(character(0), full); n_orig <- attr(r_full, "n")
b_rep <- r_full[["sel:jointd6"]][["b"]]; se_rep <- r_full[["sel:jointd6"]][["se"]]
stopifnot(is.finite(b_rep), se_rep > 0)

# ---- 3. history window (outcome year 1822-2011) + six models on one common sample ----
ds <- common_sample(d, c(vars, "recipt1_chk"), "statea", "stateb", "hy"); rm(full); invisible(gc())
fitter <- function(rhs, data) fit_hp(setdiff(rhs, vars), data)
M <- run_six(vars, ds, fitter)
saveRDS(list(M = M, r_full = r_full, n_orig = n_orig, n_used = nrow(ds)), file.path(W, "M_ADD_D01.rds"))

# ---- sensitivity hook: runs only under run_sensitivity.R (SWEEP_SENS=1), never in the main sweep ----
# heckprob (two equations; history enters both) is not one of sens_run()'s estimator classes and the run is slow -> split only.
if (Sys.getenv("SWEEP_SENS") == "1") {
  source(file.path(CODE, "hist/sens_helper.R"))
  sens_split(id = "ADD_D01", dat = ds, fit = function(dat, add) fitter(c(vars, add), dat)[["sel:jointd6"]])
}

meta <- list(pid = "ADD_D01", cite = "ADD_D01", journal = "JPR", year = 2009, doi = "10.1177/0022343308098402", family = "DEMOCRACY")
nt <- sprintf(paste0("PORT of the PRIO/JPR archive (authors' do-file Table 3 Model 3A: heckprob recipt1 | midjt1, cluster dyadid, no joiners). ",
  "Article PDF inaccessible: 'published' = our full-sample run of the authors' specification, b=%.6f se=%.6f N=%d (UNVERIFIED against print). ",
  "DV in t+1 -> history keyed to outcome year. Common sample N=%d. History terms in both equations."), b_rep, se_rep, n_orig, nrow(ds))
for (k in KEYS) write_p8(meta, gsub(":", "_", k), "Task/abstract: joint regime categories and MID onset (onset = selection equation)",
  paste("Table 3 Model 3A (heckprob, no joiners), onset equation,", k), "probit with sample selection (heckprob), onset = selection eq, reciprocation = outcome eq, SE clustered on dyad",
  "heckprob ML (R port, analytic scores), cluster-robust on dyadid, all six steps", r_full[[k]][["b"]], r_full[[k]][["se"]], n_orig, NA, FALSE,
  lapply(M, function(m) { o <- setNames(list(m[[k]]), gsub(":", "_", k)); attr(o, "n") <- attr(m, "n"); o }), "N (peace years + 3 splines)", nt, verdict = "UNVERIFIED")
dir.create(file.path(OUT, "results_add"), recursive = TRUE, showWarnings = FALSE)
for (k in KEYS) { f <- paste0("ADD_D01__", gsub(":", "_", k), ".csv")
  stopifnot(file.rename(file.path(OUT, "results_p8", f), file.path(OUT, "results_add", f))) }
sink()
