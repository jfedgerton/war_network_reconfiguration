# X034 Clare 2007 JPR "Democratization and International Conflict: The Impact of Institutional Legacies" (10.1177/0022343307076628)
# Data + do-file: PRIO JPR archive (file49642_clare-replication_materials.zip).
# Headline: Table III dyadic model (the paper's only dyad-level model; Table II is a monadic ZINB count model):
#   logit initiate stable_dem1 stable_aut1 redem1 redem1_pautdur1 redem1_pdemdur1 democratizing2 relcap alliance01
#         prop_demsregion1 growth riots1 peaceyears spline1 spline2 spline3, cluster(dyad)
# KEY terms (fixed from abstract BEFORE estimating): redem1 (redemocratizing state), redem1_pdemdur1 (redem x prior democratic duration),
#   redem1_pautdur1 (redem x prior authoritarian duration). Abstract: "Leaders in redemocratizing states, and those with a longer legacy
#   of previous democratic rule, should be more likely to initiate disputes"; "a stronger legacy of authoritarian rule ... more cautious".
# No joint-democracy / capitalist terms in the model (regime terms are monadic, initiator side) -> no CONTROL rows.
# dyad = ccode1*1000 + ccode2 (directed: initiator = ccode1); decoded for the history merge.
ROOT <- Sys.getenv("SWEEP_DATA"); OUT <- Sys.getenv("SWEEP_OUT"); CODE <- Sys.getenv("SWEEP_CODE")
stopifnot(nzchar(ROOT), nzchar(OUT), nzchar(CODE))
source(file.path(CODE, "work_p4/common_p4_b02.R"))
suppressPackageStartupMessages(library(haven))
WO <- file.path(OUT, "work_p4/X034"); dir.create(WO, recursive = TRUE, showWarnings = FALSE)  # port: outputs
W <- file.path(ROOT, "work_p4/X034"); sink(file.path(WO, "log_X034.txt"), split = TRUE)

d <- as.data.frame(zap_labels(read_dta(file.path(W, "Replication Materials/Clare-Dyadic Replication.dta"))))
d <- d[!is.na(d$dyad), ]
d$c1 <- d$dyad %/% 1000; d$c2 <- d$dyad %% 1000
stopifnot(all(d$c1 >= 2 & d$c1 <= 990), all(d$c2 >= 2 & d$c2 <= 990), all(d$c1 != d$c2), !anyDuplicated(d[, c("dyad", "year")]))
cat("rows:", nrow(d), " years:", range(d$year), " initiate mean:", mean(d$initiate, na.rm = TRUE), "\n")

f0 <- initiate ~ stable_dem1 + stable_aut1 + redem1 + redem1_pautdur1 + redem1_pdemdur1 + democratizing2 + relcap + alliance01 +
  prop_demsregion1 + growth + riots1 + peaceyears + spline1 + spline2 + spline3
K <- c("redem1", "redem1_pdemdur1", "redem1_pautdur1")
FF <- ff_logit_cl("dyad")
r0 <- FF(f0, d[complete.cases(d[, c(all.vars(f0), "dyad")]), ]); print(round(r0, 5)); n_orig <- attr(r0, "n")
cat("repro N:", n_orig, "\n")

ds <- common_sample(d, f0, "dyad", "c1", "c2", "year")
M <- run_six(f0, ds, K, FF, WO)
# ---- sensitivity hook (split only): runs only under run_sensitivity.R (SWEEP_SENS=1), never in the main sweep ----
if (Sys.getenv("SWEEP_SENS") == "1") {
  source(file.path(CODE, "hist/sens_helper.R"))
  # fit_keys without a cache file (run_six caches each model by file name)
  sens_split(id = "X034", dat = ds,
             fit = function(dat, add) { r <- fit_keys(if (length(add)) update(f0, as.formula(paste(". ~ . +", paste(add, collapse = " + ")))) else f0, dat, K, FF)$k$redem1
               c(b = unname(r[["b"]]), se = unname(r[["se"]]), p = unname(r[["p"]])) })
}

meta <- list(xid = "X034", short_cite = "Clare 2007", journal = "JPR", year = 2007, doi = "10.1177/0022343307076628")
est <- "logit, dyad-clustered SE (fixest feglm; Stata logit, cluster(dyad))"
hc <- "Y (peace years + 3 cubic splines only)"
nt <- paste0("Authors' PRIO archive data + do-file, Table III dyadic initiation model, directed dyads of democratizing-era states 1950-1990. ",
  "Published table not accessible (Sage PDF) and no log in archive -> UNVERIFIED; b_repro = estimate from authors' data/do-file. ",
  "repro N ", n_orig, "; common N ", M$m0$n, ". Table II (monadic ZINB) not used.")
write_p4(meta, "redem1", "DEMOCRACY", "KEY", "\"Leaders in redemocratizing states ... should be more likely to initiate disputes\"",
  "UNVERIFIED", NA, r0["redem1", 1], est, n_orig, M, hc, paste(nt, "redem1 = redemocratizing (constituent of the two legacy interactions: effect at 0 prior duration)."))
write_p4(meta, "redem1_pdemdur1", "DEMOCRACY", "KEY", "\"those with a longer legacy of previous democratic rule, should be more likely to initiate\"",
  "UNVERIFIED", NA, r0["redem1_pdemdur1", 1], est, n_orig, M, hc, paste(nt, "redem x prior democratic duration (product term)."))
write_p4(meta, "redem1_pautdur1", "DEMOCRACY", "KEY", "\"a stronger legacy of authoritarian rule is expected to lead to a more cautious foreign policy\"",
  "UNVERIFIED", NA, r0["redem1_pautdur1", 1], est, n_orig, M, hc, paste(nt, "redem x prior authoritarian duration (product term)."))
sink()
