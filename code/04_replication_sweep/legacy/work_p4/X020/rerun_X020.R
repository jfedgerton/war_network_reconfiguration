# X020 Park & Colaresi 2014 ISQ "Safe Across the Border" (10.1111/isqu.12114). Data + code: github.com/colaresi/SafeAcrosstheBorderReplication
# Headline = Table 1 ("TableLowDem") Model 4 (full: omitted main effects, splines, contiguity x all, + prev.mid.count, PostColdWar, majdyad):
#   zelig logit (= glm, ML SEs) cwmid ~ lowdem + lowgdp + landcontig + nlparity + beckpeace + splines + cvwonset + dyaddur + samecolony + ethborder
#     + terrsim + landcontig x {nlparity, beckpeace, spline1-3, cvwonset, dyaddur, samecolony, ethborder, terrsim} + prev.mid.count + PostColdWar + majdyad
# KEY (a priori, abstract: "even when controlling for stable border variables, democracy continues to be a consistently useful predictor"):
#   lowdem (weak-link democracy, Table 1 M4); also bothdem (joint-democracy dummy, parallel table model4zb) as a second KEY row.
# Author already controls prior-dispute history (prev.mid.count) -> author_has_history_control = Y.
ROOT <- Sys.getenv("SWEEP_DATA"); OUT <- Sys.getenv("SWEEP_OUT"); CODE <- Sys.getenv("SWEEP_CODE")
stopifnot(nzchar(ROOT), nzchar(OUT), nzchar(CODE))
source(file.path(CODE, "work_p4/common_p4_b04.R"))
WO <- file.path(OUT, "work_p4/X020"); dir.create(WO, recursive = TRUE, showWarnings = FALSE)  # port: outputs
W <- file.path(ROOT, "work_p4/X020"); sink(file.path(WO, "log_X020.txt"), split = TRUE)
e <- new.env(); load(file.path(W, "repo/SafeAcrossReplicationData"), envir = e); print(ls(e))
nm <- names(e$data1.subset); print(nm)
c1 <- intersect(c("ccode1", "statea", "ccodea", "state1", "ccode_1"), nm)[1]; c2 <- intersect(c("ccode2", "stateb", "ccodeb", "state2", "ccode_2"), nm)[1]
cat("id cols:", c1, c2, "\n"); stopifnot(!is.na(c1), !is.na(c2), "year" %in% nm)
V <- c("cwmid","lowdem","bothdem","lowgdp","landcontig","nlparity","beckpeace",".spline1",".spline2",".spline3","cvwonset","dyaddur","samecolony",
       "ethborder","terrsim","prev.mid.count","PostColdWar","majdyad")
stopifnot(all(V %in% nm))
d <- as.data.frame(e$data1.subset[, unique(c(V, c1, c2, "year"))]); rm(e); gc()
names(d) <- sub("^\\.", "", names(d))
cat("rows:", nrow(d), " years:", range(d$year), "\n")
R <- "lowgdp + landcontig + nlparity + beckpeace + spline1 + spline2 + spline3 + cvwonset + dyaddur + samecolony + ethborder + terrsim + landcontig:nlparity + landcontig:beckpeace + landcontig:spline1 + landcontig:spline2 + landcontig:spline3 + landcontig:cvwonset + landcontig:dyaddur + landcontig:samecolony + landcontig:ethborder + landcontig:terrsim + prev.mid.count + PostColdWar + majdyad"
fL <- as.formula(paste("cwmid ~ lowdem +", R)); fB <- as.formula(paste("cwmid ~ bothdem +", R))
meta <- list(xid = "X020", short_cite = "Park & Colaresi 2014", journal = "ISQ", year = 2014, doi = "10.1111/isqu.12114")
est <- "logit, conventional ML SEs (Zelig 'logit' = glm), as in authors' code"
nt <- paste0("Authors' GitHub data + R code (Gibler-2007 replication/extension data). Published tables not accessible (OUP paywall; authors' .out files not in repo) ",
  "-> UNVERIFIED. Common sample restricted to 1822-2011 and rows with history/benchmarks. Author model already includes prev.mid.count (prior MIDs).")
for (k in list(list(f = fL, t = "lowdem", ev = "\"democracy continues to be a consistently useful predictor of international conflict\" (Table 1 M4 lowdem)"),
               list(f = fB, t = "bothdem", ev = "same claim, joint-democracy dummy version (model4zb)"))) {
  Wk <- file.path(WO, k$t); dir.create(Wk, recursive = TRUE, showWarnings = FALSE)
  r0 <- ff_logit_iid(k$f, d); print(round(r0[k$t, , drop = FALSE], 5)); n_orig <- attr(r0, "n")
  ds <- common_sample(d, k$f, character(0), c1, c2, "year")
  M <- run_six(k$f, ds, k$t, ff_logit_iid, Wk)
  # ---- sensitivity hook: runs only under run_sensitivity.R (SWEEP_SENS=1), never in the main sweep ----
  # (headline row = bothdem model fB; controls = every term of the script's R string, i.e. all terms of fB except bothdem)
  if (Sys.getenv("SWEEP_SENS") == "1" && k$t == "bothdem") {
    source(file.path(CODE, "hist/sens_helper.R"))
    sens_run(id = "X020", dat = ds, y = "cwmid", tie = "bothdem", ctrl = setdiff(attr(terms(fB), "term.labels"), "bothdem"), fe = NULL,
             fit = function(dat, add) { ct <- ff_logit_iid(if (length(add)) update(fB, as.formula(paste(". ~ . +", paste(add, collapse = " + ")))) else fB, dat)
               c(b = unname(ct["bothdem", 1]), se = unname(ct["bothdem", 2]), p = unname(ct["bothdem", 3])) })
  }
  write_p4(meta, k$t, "DEMOCRACY", "KEY", k$ev, "UNVERIFIED", NA, r0[k$t, 1], est, n_orig, M, "Y (prev.mid.count: count of previous MIDs)",
    paste(nt, if (k$t == "lowdem") "Headline Table 1 Model 4." else "Parallel bothdem table, Model 4 spec."))
}
sink()
