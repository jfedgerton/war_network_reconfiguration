# X003 Weeks 2012 APSR "Strongmen and Straw Men" (10.1017/s0003055412000111)
# Data + do-file: jessicalpweeks.com (weeksapsr2012.zip). Published coefficients: web appendix Table 1 (full results).
# Headline: Table 1 Model 2 (cross-sectional logit with full controls), SE clustered on directed dyad (dirdyadid),
# sample democracy_1 != . ; directed dyad-years, Maoz MID initiation (mzinit); peace years + 3 splines.
# KEY terms (a priori; abstract: "substantial variation in conflict initiation occurs among authoritarian regimes";
#   hypotheses on juntas, bosses, strongmen (more belligerent) and machines (no more than democracies)):
#   strongmanjlw_1, juntajlw_1, bossjlw_1, machinejlw_1 (omitted category: democracy). No joint-democracy/capitalist controls.
ROOT <- Sys.getenv("SWEEP_DATA"); OUT <- Sys.getenv("SWEEP_OUT"); CODE <- Sys.getenv("SWEEP_CODE")
stopifnot(nzchar(ROOT), nzchar(OUT), nzchar(CODE))
source(file.path(CODE, "work_p4/b3/common_b3.R"))
WO <- file.path(OUT, "work_p4/X003"); dir.create(WO, recursive = TRUE, showWarnings = FALSE)  # port: outputs
W <- file.path(ROOT, "work_p4/X003"); sink(file.path(WO, "log_X003.txt"), split = TRUE)
m <- list(xid = "X003", short_cite = "Weeks 2012", journal = "APSR", year = 2012, doi = "10.1017/s0003055412000111")

rhsv <- c("machinejlw_1","juntajlw_1","bossjlw_1","strongmanjlw_1","allotherauts_1","newregime_1","democracy_2","cap_1","cap_2",
          "initshare","dependlow","majmaj","minmaj","majmin","contigdum","logdist","s_wt_glo","s_lead_1","s_lead_2",
          "pcyrsmzinit","pcyrsmzinits1","pcyrsmzinits2","pcyrsmzinits3")
vars <- c("ccode1","ccode2","year","dirdyadid","mzinit","democracy_1", rhsv)
d <- haven::read_dta(file.path(W, "WeeksAPSR2012/WeeksAPSR2012.dta"), col_select = all_of(vars))
d <- as.data.frame(lapply(d, as.numeric)); gc()
d <- d[!is.na(d$democracy_1), ]                         # "if democracy_1!=."
cat("years:", range(d$year), " n:", nrow(d), "\n")
f0 <- as.formula(paste("mzinit ~", paste(rhsv, collapse = " + ")))
terms <- c("strongmanjlw_1","juntajlw_1","bossjlw_1","machinejlw_1")
rep <- fit_keys(f0, d, terms, "~dirdyadid", "logit", file.path(WO, "cache_repro.rds"))
b_pub <- c(strongmanjlw_1 = 0.832, juntajlw_1 = 0.515, bossjlw_1 = 0.649, machinejlw_1 = -0.459)  # web appendix Table 1 col 2
sig_pub <- c(strongmanjlw_1 = TRUE, juntajlw_1 = TRUE, bossjlw_1 = TRUE, machinejlw_1 = TRUE)       # all ***
st <- sapply(terms, function(t) repro_status(b_pub[t], sig_pub[t], rep$k[[t]])); print(st)

ds <- common_sample(d, f0, "dirdyadid", "ccode1", "ccode2", "year"); rm(d); gc()
M <- run_six(f0, ds, terms, "~dirdyadid", "logit", WO)
# ---- sensitivity hook: runs only under run_sensitivity.R (SWEEP_SENS=1), never in the main sweep ----
if (Sys.getenv("SWEEP_SENS") == "1") {
  source(file.path(CODE, "hist/sens_helper.R"))
  sens_run(id = "X003", dat = ds, y = "mzinit", tie = "strongmanjlw_1",
           ctrl = setdiff(rhsv, "strongmanjlw_1"), fe = NULL,
           fit = function(dat, add) { r <- fit_keys(if (length(add)) update(f0, as.formula(paste(". ~ . +", paste(add, collapse = " + ")))) else f0, dat, "strongmanjlw_1", "~dirdyadid", "logit")$k$strongmanjlw_1
             c(b = unname(r[["b"]]), se = unname(r[["se"]]), p = unname(r[["p"]])) })
}
est <- "logit, SE clustered on directed dyad (fixest feglm)"
hc <- "Y (peace years since last initiation + 3 cubic splines)"
nt <- paste0("Table 1 Model 2 (logit, full controls); b_pub from web appendix Table 1 col 2 (3 decimals); repro N=", rep$n,
             "; common sample drops ", attr(ds, "n_cc") - nrow(ds), " obs (benchmark/history missingness). Directed dyad-years (Weeks' leader-regime data), Maoz MID initiation. ",
             "Autocratic subtype dummies vs democracy baseline (initiator side). ")
ev <- "abstract: 'substantial variation in conflict initiation occurs among authoritarian regimes'"
for (t in terms)
  write_p4(m, "DEMOCRACY", "KEY", t, ev, st[t], b_pub[t], rep$k[[t]]["b"], est, rep$n, M, hc,
           paste0(nt, "bD (lagged joint democracy) is a benchmark, not an exact placebo for regime-type terms."))
sink()
