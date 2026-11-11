# X078 Peterson & Rudloff 2015 II, "Preferential Trade Agreements and Trade Expectations Theory"
# Key term (decided from abstract BEFORE estimating): signed-but-not-in-force PTA (F.signed_noforce).
#   Abstract: "signed PTAs are pacifying, while in-force agreements have no statistically significant impact"
# Headline = model m2 of ii_replication_3_11.do (the signed-vs-in-force split, all MIDs):
#   logit F.cwmid F.signed_noforce F.entered_force lowerdep minchdep mingrowth F.contigdi F.lndist lowerpolity atopally caprat
#         cubic(F.py_mid) coldwar if F.cwongo==0, cluster(dyadid)
# Stata F. = value at year+1 within dyad (xtset dyadid year); built exactly that way here.
# lowerpolity (weak-link democracy) recorded as a CONTROL row. No author log -> UNVERIFIED.
ROOT <- Sys.getenv("SWEEP_DATA"); OUT <- Sys.getenv("SWEEP_OUT"); CODE <- Sys.getenv("SWEEP_CODE")
stopifnot(nzchar(ROOT), nzchar(OUT), nzchar(CODE))
source(file.path(CODE, "work_p4/b6/common_b6.R"))
WO <- file.path(OUT, "work_p4/X078"); dir.create(WO, recursive = TRUE, showWarnings = FALSE)  # port: outputs
W <- file.path(ROOT, "work_p4/X078"); sink(file.path(WO, "log_X078.txt"), split = TRUE)
m <- list(xid = "X078", short_cite = "Peterson & Rudloff 2015", journal = "II", year = 2015, doi = "10.1080/03050629.2014.948153")
d <- as.data.frame(haven::zap_labels(haven::read_dta(file.path(W, "ii_replication_3_11.dta"))))
d$dyadid <- d$ccode1 * 1000 + d$ccode2
stopifnot(!anyDuplicated(d[, c("dyadid", "year")]))
cat("rows:", nrow(d), " years:", range(d$year), "\n")
# leads: F.x = x in (dyadid, year+1)
lv <- c("cwmid", "signed_noforce", "entered_force", "contigdi", "lndist", "py_mid", "cwongo")
nx <- d[, c("dyadid", "year", lv)]; nx$year <- nx$year - 1; names(nx)[-(1:2)] <- paste0("F_", lv)
d <- left_join(d, nx, by = c("dyadid", "year"))
d$F_py2 <- d$F_py_mid^2; d$F_py3 <- d$F_py_mid^3
d <- d[!is.na(d$F_cwongo) & d$F_cwongo == 0, ]
stopifnot(all(d$F_cwmid %in% c(0, 1, NA)))
f0 <- F_cwmid ~ F_signed_noforce + F_entered_force + lowerdep + minchdep + mingrowth + F_contigdi + F_lndist + lowerpolity +
  atopally + caprat + F_py_mid + F_py2 + F_py3 + coldwar
K <- c("F_signed_noforce", "lowerpolity")
r <- fit_keys(f0, d, K, "~dyadid", "logit", file.path(WO, "cache_repro.rds"))
d$dvyear <- d$year + 1   # history anchored on the DV year (t+1)
ds <- common_sample(d, f0, "dyadid", "ccode1", "ccode2", "dvyear"); rm(d); gc()
M <- run_six(f0, ds, K, "~dyadid", "logit", WO)
# ---- sensitivity hook: runs only under run_sensitivity.R (SWEEP_SENS=1), never in the main sweep ----
if (Sys.getenv("SWEEP_SENS") == "1") {
  source(file.path(CODE, "hist/sens_helper.R"))
  sens_run(id = "X078", dat = ds, y = "F_cwmid", tie = "F_signed_noforce",   # [2026-11-10 audit2] hook moved to the new headline term (F_signed_noforce)
           ctrl = c("lowerpolity", "F_entered_force", "lowerdep", "minchdep", "mingrowth", "F_contigdi", "F_lndist", "atopally", "caprat",
                    "F_py_mid", "F_py2", "F_py3", "coldwar"), fe = NULL,
           fit = function(dat, add) { r <- fit_keys(if (length(add)) update(f0, as.formula(paste(". ~ . +", paste(add, collapse = " + ")))) else f0, dat, "F_signed_noforce", "~dyadid", "logit")$k$F_signed_noforce
             c(b = unname(r[["b"]]), se = unname(r[["se"]]), p = unname(r[["p"]])) })
}
nt <- paste0("Model 2 of authors' do-file (table not accessible; no author log) -> UNVERIFIED. DV F.cwmid (MID onset t+1), 1957-2000 per paper. ",
  "History anchored on DV year t+1; dropped vs complete-case: ", r$n - M$m0$n, ".")
write_p4(m, "TIE", "KEY", "F_signed_noforce", "abstract: 'signed PTAs are pacifying, while in-force agreements have no statistically significant impact'",
  "UNVERIFIED", NA, unname(r$k$F_signed_noforce["b"]), "logit, dyad-clustered SE (nondirected dyad-year MID onset t+1)", r$n, M,
  "Y (cubic peace-years F.py_mid)", nt)
write_p4(m, "DEMOCRACY", "CONTROL", "lowerpolity", "control in headline model (weak-link Polity)",
  "UNVERIFIED", NA, unname(r$k$lowerpolity["b"]), "logit, dyad-clustered SE (nondirected dyad-year MID onset t+1)", r$n, M,
  "Y (cubic peace-years F.py_mid)", paste(nt, "bD (lagged joint democracy) is not a placebo for this term."))
sink()
