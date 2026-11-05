# ESTAUDIT X019 Gibler 2014 ISQ -- copy of rerun_X019.R with the ARTICLE's SE type.
# Article (ISQ 2014 Table 1 note, p.128): "Logit regression models with standard errors in parentheses" -- no robust/cluster
# correction stated; the replicated original (Gibler 2007 ISQ Tables 1-3) also reports plain SEs (Park & Colaresi 2014 fn 13:
# Gibler 2007 did not use robust/clustered SEs). Temporal terms: peace years + square + cube (p.127) -- already in the model.
# Earlier rerun used heteroskedasticity-robust SEs (from the book ch.7 log). Change here: vcov "iid" (conventional ML SEs) in
# ALL six steps; identical common sample and point estimates. Model spec is still the book ch.7 model (ISQ Table 1 data not public).
ROOT <- Sys.getenv("SWEEP_DATA"); OUT <- Sys.getenv("SWEEP_OUT"); CODE <- Sys.getenv("SWEEP_CODE")
stopifnot(nzchar(ROOT), nzchar(OUT), nzchar(CODE))
source(file.path(CODE, "work_p4/b3/common_b3.R"))
source(file.path(CODE, "ea4/est_helpers.R"))
WO <- file.path(OUT, "work_p4/X019"); dir.create(WO, recursive = TRUE, showWarnings = FALSE)  # port: outputs
W <- file.path(ROOT, "work_p4/X019"); WC <- file.path(WO, "est_cache"); dir.create(WC, recursive = TRUE, showWarnings = FALSE)
sink(file.path(WO, "log_est_X019.txt"), split = TRUE)
rhsv <- c("jointdem7","caprat","allied","mzpceyrs","mzpceyrs2","mzpceyrs3","lndyaddur1","samemaster","peacetrans","violenttrans",
          "low_peace_defpact","either_cowintra","high_threat_militarization","high_threat_terrMID","high_threat_terrMIDxmilpe")
vars <- c("ccode1","ccode2","year","contig","include1816","mzfatald", rhsv)
d <- haven::read_dta(file.path(W, "Chapter 7 replication data/Chapter 7 testing data.dta"), col_select = all_of(vars))
d <- as.data.frame(lapply(d, as.numeric)); gc()
d <- d[!is.na(d$contig) & d$contig == 1 & !is.na(d$include1816) & d$include1816 == 1, ]
d$fatal <- as.numeric(d$mzfatald != 0); d$fatal[is.na(d$mzfatald)] <- NA
f0 <- as.formula(paste("fatal ~", paste(rhsv, collapse = " + ")))
rep <- fit_keys(f0, d, "jointdem7", "iid", "logit", file.path(WC, "cache_repro.rds"))
stopifnot(rep$n == 19026, abs(rep$k$jointdem7["b"] - (-0.0604678)) < 1e-4)   # book log model 4 point estimate
ds <- common_sample(d, f0, character(0), "ccode1", "ccode2", "year")
M <- run_six(f0, ds, "jointdem7", "iid", "logit", WC)
stopifnot(all(sapply(M, function(z) z$n) == M$m0$n), M$m0$n == 18416)

# ---- sensitivity hook: runs only under run_sensitivity.R (SWEEP_SENS=1), never in the main sweep ----
if (Sys.getenv("SWEEP_SENS") == "1") {
  source(file.path(CODE, "hist/sens_helper.R"))
  sens_run(id = "X019", dat = ds, y = "fatal", tie = "jointdem7", ctrl = setdiff(rhsv, "jointdem7"), fe = NULL,
           fit = function(dat, add) { r <- fit_keys(if (length(add)) update(f0, as.formula(paste(". ~ . +", paste(add, collapse = " + ")))) else f0, dat, "jointdem7", "iid", "logit")$k$jointdem7
             c(b = unname(r[["b"]]), se = unname(r[["se"]]), p = unname(r[["p"]])) })
}

R <- lapply(M, function(z) z$k$jointdem7)
f_old <- file.path(OUT, "results_verify/X019__jointdem7.csv")  # port: prior result (check + notes); absent in pipeline
if (file.exists(f_old)) { old <- read.csv(f_old); stopifnot(abs(R$m0["b"] - old$b0) < 1e-8)   # same point estimate as before; only SEs change
} else old <- list(b0 = NA_real_, p0 = NA_real_, p_dm = NA_real_)
st <- "logit (fixest feglm), conventional ML SE, peace yrs+sq+cube"
write_est("X019", "Gibler 2014", "jointdem7",
  "logit, conventional (non-robust) SE, peace years + square + cube; contiguous dyads 1946-2001",
  "pdf p.127-128 Table 1 note ('Logit regression models with standard errors'); Gibler 2007 ISQ same convention",
  rep(st, 6), "N",
  "all six steps used heteroskedasticity-robust (HC) SEs (book ch.7 log); article reports plain ML SEs. Family/link and temporal terms already matched.",
  "Y", R, M$m0$n,
  paste0("SE-only fix: point estimates identical to results_verify (b0 ", signif(old$b0, 5), "); old robust p0 ", signif(old$p0, 3),
         ". Model is still the book ch.7 model 4 (fatal MIDs, jointdem7), not ISQ Table 1 (data not public; key flagged in verify: ISQ uses lowest democracy score). ",
         "Sig. change vs robust: m0 p ", signif(old$p0, 3), " -> ", signif(R$m0["p"], 3), "; mDM p ", signif(old$p_dm, 3), " -> ", signif(R$mDM["p"], 3), "."))
sink()
