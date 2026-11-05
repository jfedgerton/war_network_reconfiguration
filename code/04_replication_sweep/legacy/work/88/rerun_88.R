# Rerun: Kinne 2014 CMPS, Table 1 Model 3 (logit, SE clustered on dyadid, BKT peace-years splines)
# Headline tie: signaltriads (3rd-party signaling ties). b_pub = -9.136* (4.383), N = 943,152 (paper Table 1)
# Steps are cached (cache_*.rds) because the shared box OOM-kills large R processes; rerun until complete.
ROOT <- Sys.getenv("SWEEP_DATA"); OUT <- Sys.getenv("SWEEP_OUT"); CODE <- Sys.getenv("SWEEP_CODE")
stopifnot(nzchar(ROOT), nzchar(OUT), nzchar(CODE))
source(file.path(CODE, "work/rerun_helpers_b6.R"))
suppressPackageStartupMessages({library(haven); library(fixest)})
W <- file.path(OUT, "work/88"); dir.create(W, recursive = TRUE, showWarnings = FALSE); sink(file.path(W, "log_88.txt"), split = TRUE)   # port: W = step caches/log (OUT)
TIE <- "signaltriads"
X <- c("costlytriads","signaltriads","demlow","dep_ab","sum_lncap","allies","contig","lndist","systsize",
       "peaceyrs_mid","x_spline1_mid","x_spline2_mid","x_spline3_mid")   # Stata _spline* renamed x_spline*
SSC <- ssc(adj = FALSE, cluster.adj = TRUE)  # Stata logit vce(cluster): G/(G-1) only
fit <- function(dat, add) feglm(reformulate(c(X, add), "mid"), family = binomial, data = dat, mem.clean = TRUE,
                                cluster = ~dyadid, ssc = SSC, nthreads = 1)

# 1. load
base <- cached(W, "base", {
  d <- read_dta(file.path(ROOT, "work/88", "CMPS-Replication/data.dta"),
                col_select = any_of(c("mid", sub("^x_", "_", X), "dyadid","ccode1","ccode2","year"))) |>
    zap_labels() |> as.data.frame()
  names(d) <- sub("^_", "x_", names(d))
  d[complete.cases(d[, c("mid", X, "dyadid")]), ]
})
n_orig <- nrow(base); cat("n_orig =", n_orig, "\n"); stopifnot(n_orig == 943152)

# 2. reproduce
rep <- cached(W, "rep", { m <- fit(base, character(0)); list(ct = coeftable(m), n = nobs(m)) })
print(rep$ct)
b_pub <- -9.136; b_repro <- unname(rep$ct[TIE, 1]); p_repro <- rep$ct[TIE, 4]
repro_ok <- abs(b_repro - b_pub) / abs(b_pub) <= .10 && p_repro < .05
cat("b_repro =", b_repro, " p =", p_repro, " ok =", repro_ok, "\n")

# 3. history (directed dyads; history is symmetric)
s <- cached(W, "s", {
  d <- add_history(base, "ccode1", "ccode2", "year")
  d <- d[d$year >= 1822 & d$year <= 2011, ]
  d[complete.cases(d[, c("hist_dyad_w15","hist_mono_w15","hist_dyad_w26","hist_mono_w26","bm_cinc_w15","bm_democ_w15")]), ]
})
rm(base); invisible(gc())
n_used <- nrow(s); cat("n_used =", n_used, " dropped =", n_orig - n_used, "\n")
print(summary(s[, c("hist_dyad_w15","hist_mono_w15","bm_cinc_w15","bm_democ_w15")]))

# 4. model sequence
r <- list()
for (nm in names(ADD_SETS)) {
  res <- cached(W, paste0("m_", nm), { m <- fit(s, ADD_SETS[[nm]]); list(ct = coeftable(m), n = nobs(m)) })
  cat("\n====", nm, " N =", res$n, "\n"); print(res$ct); stopifnot(res$n == nrow(s))
  r[[nm]] <- getc(res$ct, TIE)
}

# ---- sensitivity hook: runs only under run_sensitivity.R (SWEEP_SENS=1), never in the main sweep ----
if (Sys.getenv("SWEEP_SENS") == "1") {
  source(file.path(CODE, "hist/sens_helper.R"))
  sens_run(id = "88", dat = s, y = "mid", tie = TIE, ctrl = setdiff(X, TIE), fe = NULL,
           fit = function(dat, add) getc(coeftable(fit(dat, add)), TIE))
}

# 4b. sensitivity: benchmark (CINC/polity) missingness drops ~5% of obs (mostly new/micro states).
#     Re-run m0 and mDM on the year-restricted sample WITHOUT requiring benchmarks (hist vars are never missing).
rm(s); invisible(gc())
sens <- cached(W, "sens", {
  b <- readRDS(file.path(W, "cache_base.rds")); b <- b[b$year >= 1822 & b$year <= 2011, ]
  b <- add_history(b, "ccode1", "ccode2", "year"); b$bm_cinc_w15 <- NULL; b$bm_democ_w15 <- NULL; invisible(gc())
  m0 <- fit(b, character(0)); c0 <- getc(coeftable(m0), TIE); n <- nobs(m0); rm(m0); invisible(gc())
  mdm <- fit(b, ADD_SETS$mDM); cdm <- getc(coeftable(mdm), TIE)
  list(n = n, c0 = c0, cdm = cdm)
})
cat("\n==== sensitivity (no benchmark requirement): N =", sens$n, "\n"); print(rbind(m0 = sens$c0, mDM = sens$cdm))
sens_note <- sprintf(" Sensitivity w/o benchmark-missingness drop (N=%d): b0=%.3f (p=%.3f), b_dm=%.3f (p=%.3f), pct_dm=%.1f.",
                     sens$n, sens$c0["b"], sens$c0["p"], sens$cdm["b"], sens$cdm["p"], pct(sens$cdm["b"], sens$c0["b"]))

row <- list(cid = "88", short_cite = "Kinne 2013", reproduced = ifelse(repro_ok, "REPRODUCED", "NOT_REPRODUCED"),
  b_pub = b_pub, b_repro = b_repro, estimator = "logit, cluster dyadid (fixest)", tie_var = "signaltriads",
  n_orig = n_orig, n_used = n_used, author_has_history_control = "Y",
  notes = paste0("Table 1 M3 (paper: -9.136*, SE 4.383; exact repro). Directed dyads 1950-2000; dropped ", n_orig - n_used,
                 " obs (missing CINC/polity benchmarks, mostly new states). Author ch_control=PEACE_YEARS (BKT splines, kept). costlytriads also in model.", sens_note))
row <- fill_seq(row, r)
write_result(row)
sink()
