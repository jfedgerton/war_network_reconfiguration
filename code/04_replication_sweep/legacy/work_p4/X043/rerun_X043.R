# X043 Mousseau 2016/2018 CMPS "Grasping the scientific evidence: The contractualist peace supersedes the democratic peace"
# (10.1177/0738894215616408). Data + do-file + outputs: UCF author page (Replication-Data.zip: CMPS2016.dta, Do File.do, Outputs.pdf).
# Headline: Table 2 Model 2 (listwise): logit fatal MID onset (mzfmid) on D (weak-link democracy) + CIEl (contract-intensive
# economy, lagged) + controls + fatal-MID peace years & 3 splines; SE clustered on dyadid. Nondirected dyad-years 1951-2001.
# KEY terms (a priori from abstract: "contractualist economy accounts for the democratic peace"; "no support for democracy"):
#   CIEl (CAPITALISM) and D (DEMOCRACY; claim = null once CIE included).
ROOT <- Sys.getenv("SWEEP_DATA"); OUT <- Sys.getenv("SWEEP_OUT"); CODE <- Sys.getenv("SWEEP_CODE")
stopifnot(nzchar(ROOT), nzchar(OUT), nzchar(CODE))
source(file.path(CODE, "work_p4/b3/common_b3.R"))
WO <- file.path(OUT, "work_p4/X043"); dir.create(WO, recursive = TRUE, showWarnings = FALSE)  # port: outputs
W <- file.path(ROOT, "work_p4/X043"); sink(file.path(WO, "log_X043.txt"), split = TRUE)
m <- list(xid = "X043", short_cite = "Mousseau 2016", journal = "CMPS", year = 2016, doi = "10.1177/0738894215616408")

vars <- c("ccode1","ccode2","year","dyadid","mzfmid","D","CIEl","lncprt","mjpw","cntg5","nmste","dist","fpceyrs","fspl1","fspl2","fspl3")
d <- haven::read_dta(file.path(W, "CMPS2016.dta"), col_select = all_of(vars))
d <- as.data.frame(lapply(d, as.numeric))
f0 <- mzfmid ~ D + CIEl + lncprt + mjpw + cntg5 + nmste + dist + fpceyrs + fspl1 + fspl2 + fspl3
terms <- c("CIEl", "D")
rep <- fit_keys(f0, d, terms, "~dyadid", "logit", file.path(WO, "cache_repro.rds"))
cat("repro N =", rep$n, "(author: 349635)\n")
b_pub <- c(CIEl = log(0.4075603), D = log(1.018974))   # Outputs.pdf, T2 model 2 (odds ratios -> logit coef)
sig_pub <- c(CIEl = TRUE, D = FALSE)                     # p = .000, .430
st <- sapply(terms, function(t) repro_status(b_pub[t], sig_pub[t], rep$k[[t]])); print(st)

ds <- common_sample(d, f0, "dyadid", "ccode1", "ccode2", "year")
M <- run_six(f0, ds, terms, "~dyadid", "logit", WO)
# ---- sensitivity hook: runs only under run_sensitivity.R (SWEEP_SENS=1), never in the main sweep ----
if (Sys.getenv("SWEEP_SENS") == "1") {
  source(file.path(CODE, "hist/sens_helper.R"))
  sens_run(id = "X043", dat = ds, y = "mzfmid", tie = "CIEl",
           ctrl = c("D", "lncprt", "mjpw", "cntg5", "nmste", "dist", "fpceyrs", "fspl1", "fspl2", "fspl3"), fe = NULL,
           fit = function(dat, add) { r <- fit_keys(if (length(add)) update(f0, as.formula(paste(". ~ . +", paste(add, collapse = " + ")))) else f0, dat, "CIEl", "~dyadid", "logit")$k$CIEl
             c(b = unname(r[["b"]]), se = unname(r[["se"]]), p = unname(r[["p"]])) })
}
est <- "logit, SE clustered on dyad (fixest feglm)"
hc <- "Y (fatal-MID peace years + 3 cubic splines)"
nt <- paste0("Table 2 Model 2 (listwise; the MI version is model 4); b_pub from author Outputs.pdf (odds ratios logged); repro N=", rep$n,
             " (author 349,635); common sample drops ", attr(ds, "n_cc") - nrow(ds), " obs (benchmark/history missingness). DV fatal MID onset 1951-2001. ")
write_p4(m, "CAPITALISM", "KEY", "CIEl", "abstract: 'The democratic peace is a statistical artifact explained by contractualist economy'",
         st["CIEl"], b_pub["CIEl"], rep$k$CIEl["b"], est, rep$n, M, hc, nt)
write_p4(m, "DEMOCRACY", "KEY", "D", "abstract: 'Results are consistent across all tests: there is no support for democracy'",
         st["D"], b_pub["D"], rep$k$D["b"], est, rep$n, M, hc,
         paste0(nt, "D = weak-link (lower) Polity score; author's claim is that it is NOT significant with CIE. bD (lagged joint democracy) is NOT a placebo for this term."))
sink()
