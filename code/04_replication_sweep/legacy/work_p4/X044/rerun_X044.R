# X044 Owsiak 2016 CMPS "Foundations for integrating the democratic and territorial peace arguments" (10.1177/0738894216650635)
# Data/do-file/log: Dataverse doi:10.7910/DVN/6GKLIR ("MID Onset.dta", Do.do, All Published Results.smcl).
# Headline = Table 2 Model 2 (main model, 1816-2001, contiguous dyads):
#   logit midonset jtdem6 settlem interaction majinvolve atopally lncincratio peaceyears _midspline1-3 if conttype==1, cluster(cdyad)
# KEY terms (a priori, abstract): interaction = jtdem6 x settlem ("the most peaceful, contiguous dyads require both democracy and border
#   settlement"); jtdem6 ("Democratic dyads ... experience no such effect in the absence of border settlement"); settlem ("border settlement
#   significantly reduces conflict, even for non-democratic dyads"; territorial-peace tie term).
# Published (authors' smcl log): jtdem6 .2018305 (p .556); settlem -.4582362 (p<.001); interaction -.7458922 (p .082); N 16,845.
ROOT <- Sys.getenv("SWEEP_DATA"); OUT <- Sys.getenv("SWEEP_OUT"); CODE <- Sys.getenv("SWEEP_CODE")
stopifnot(nzchar(ROOT), nzchar(OUT), nzchar(CODE))
source(file.path(CODE, "work_p4/common_p4_b04.R"))
suppressPackageStartupMessages(library(haven))
WO <- file.path(OUT, "work_p4/X044"); dir.create(WO, recursive = TRUE, showWarnings = FALSE)  # port: outputs
W <- file.path(ROOT, "work_p4/X044"); sink(file.path(WO, "log_X044.txt"), split = TRUE)
fn <- file.path(W, "midonset.dta")
nm <- names(read_dta(fn, n_max = 1)); cat("columns:", length(nm), "\n"); print(nm[1:min(80, length(nm))])
idv <- intersect(nm, c("ccode1", "ccode2", "statea", "stateb", "state1", "state2", "year"))
keep <- c("midonset","jtdem6","settlem","interaction","majinvolve","atopally","lncincratio","peaceyears","_midspline1","_midspline2","_midspline3","conttype","cdyad", idv)
stopifnot(all(setdiff(keep, idv) %in% nm), "year" %in% idv)
d <- zap_labels(read_dta(fn, col_select = all_of(keep))) %>% as.data.frame()
d <- d[!is.na(d$conttype) & d$conttype == 1, ]
names(d) <- sub("^_", "", names(d))
c1 <- intersect(c("ccode1", "statea", "state1"), names(d))[1]; c2 <- intersect(c("ccode2", "stateb", "state2"), names(d))[1]
cat("id cols:", c1, c2, " rows (contiguous):", nrow(d), " years:", range(d$year), " directed dup check:", anyDuplicated(d[, c(c1, c2, "year")]), "\n")
f0 <- midonset ~ jtdem6 + settlem + interaction + majinvolve + atopally + lncincratio + peaceyears + midspline1 + midspline2 + midspline3
terms <- c("interaction", "jtdem6", "settlem")
fit <- ff_logit_cl("cdyad")
r0 <- fit(f0, d[complete.cases(d[, c(all.vars(f0), "cdyad")]), ]); print(round(r0, 5)); n_orig <- attr(r0, "n")
b_pub <- c(interaction = -.7458922, jtdem6 = .2018305, settlem = -.4582362); sig_pub <- c(interaction = FALSE, jtdem6 = FALSE, settlem = TRUE)
rep <- sapply(terms, function(t) if (abs(r0[t, 1] - b_pub[t]) / abs(b_pub[t]) <= .10 && (r0[t, 3] < .05) == sig_pub[t]) "REPRODUCED" else "NOT_REPRODUCED")
print(rep)
ds <- common_sample(d, f0, "cdyad", c1, c2, "year")
M <- run_six(f0, ds, terms, fit, WO)
# ---- sensitivity hook (split only): runs only under run_sensitivity.R (SWEEP_SENS=1), never in the main sweep ----
if (Sys.getenv("SWEEP_SENS") == "1") {
  source(file.path(CODE, "hist/sens_helper.R"))
  # fit_keys is called without a cache file (run_six caches by step name)
  sens_split(id = "X044", dat = ds,
             fit = function(dat, add) { r <- fit_keys(if (length(add)) update(f0, as.formula(paste(". ~ . +", paste(add, collapse = " + ")))) else f0, dat, terms, fit)$k$settlem
               c(b = unname(r[["b"]]), se = unname(r[["se"]]), p = unname(r[["p"]])) })
}
meta <- list(xid = "X044", short_cite = "Owsiak 2016", journal = "CMPS", year = 2016, doi = "10.1177/0738894216650635")
est <- "logit, SE clustered on cdyad (fixest, Stata-style small-sample adj.)"
nt <- paste0("Authors' data + do-file, Table 2 Model 2 (contiguous dyads 1816-2001). Common sample restricted to 1822-2011 (drops 1816-1821 and rows missing benchmarks). ",
  "Published values from authors' smcl log (N 16,845).")
hc <- "N (peace years + 3 splines only)"
write_p4(meta, "interaction", "DEMOCRACY", "KEY", "\"even the most peaceful, contiguous dyads require both democracy and border settlement\"", rep["interaction"],
  b_pub["interaction"], r0["interaction", 1], est, n_orig, M, hc, paste(nt, "Key = jtdem6 x settlem product term (published p=.082, n.s. at .05)."))
write_p4(meta, "jtdem6", "DEMOCRACY", "KEY", "\"Democratic dyads, however, experience no such effect in the absence of border settlement\"", rep["jtdem6"],
  b_pub["jtdem6"], r0["jtdem6", 1], est, n_orig, M, hc, paste(nt, "Joint democracy (both polity>=6) constituent = effect among unsettled contiguous dyads."))
write_p4(meta, "settlem", "TIE", "KEY", "\"Within contiguous dyads, border settlement significantly reduces conflict\"", rep["settlem"],
  b_pub["settlem"], r0["settlem", 1], est, n_orig, M, hc, paste(nt, "Border-settlement agreement (dyadic institutional tie; territorial peace) constituent = effect among non-joint-democracies."))
sink()
