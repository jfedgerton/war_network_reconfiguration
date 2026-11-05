# X087 Owsiak 2012 ISQ "Signing Up for Peace" (10.1111/j.1468-2478.2011.00699.x)
# KEY (fixed from abstract before estimating): settlem (settled border via boundary agreement).
#   Evidence: "once neighboring states settle their borders, they are less likely to go to war or experience MIDs"
# CONTROL (DEMOCRACY): jtdem6 (joint democracy, polity>=6).
# Headline: Table 3 Model 3 (author do-file): logit midonset settlem jtdem6 majinvolve atopally lncincratio peaceyears _spline1-3,
#   default (non-robust) SEs. Data: Dataverse doi:10.7910/DVN/40LBZP 'Onset' file (original .dta). Settlem Patch.do NOT applied
#   (published results predate it; newer deposit WJUNCD omits it).
ROOT <- Sys.getenv("SWEEP_DATA"); OUT <- Sys.getenv("SWEEP_OUT"); CODE <- Sys.getenv("SWEEP_CODE")
stopifnot(nzchar(ROOT), nzchar(OUT), nzchar(CODE))
source(file.path(CODE, "work_p4/b7/common.R"))
suppressMessages(library(readstata13))
WO <- file.path(OUT, "work_p4/X087"); dir.create(WO, recursive = TRUE, showWarnings = FALSE)  # port: outputs
W <- file.path(ROOT, "work_p4/X087"); logf <- file.path(WO, "log_X087.txt"); cat("", file = logf)
meta <- list(xid = "X087", cite = "Owsiak 2012", journal = "ISQ", year = 2012, doi = "10.1111/j.1468-2478.2011.00699.x")
d <- read.dta13(file.path(W, "onset.dta"), convert.factors = FALSE)
names(d) <- sub("^_", "x_", names(d))
rhs <- c("settlem","jtdem6","majinvolve","atopally","lncincratio","peaceyears","x_spline1","x_spline2","x_spline3")
vars <- c("midonset", rhs)
cat("raw N:", nrow(d), " years:", range(d$year), "\n"); print(table(d$conttype, useNA = "ifany"))
fitter <- function(r, dat) { m <- fit_bin(r, "midonset", dat); list(ct = coeftable(m), n = nobs(m)) }
# 1. reproduce on full data
d0 <- d[complete.cases(d[, vars]), ]; n_orig <- nrow(d0)
rep <- fitter(rhs, d0); print(round(rep$ct, 4)); capture.output(cat("REPRO\n"), print(rep$ct), file = logf, append = TRUE)
cat("n_orig =", n_orig, "\n")
# 2-4. history + six models
s <- common_sample(d, vars, "ccode1", "ccode2", "year")
R <- run_six(s, rhs, fitter, c("settlem","jtdem6"), logf)
# ---- sensitivity hook: runs only under run_sensitivity.R (SWEEP_SENS=1), never in the main sweep ----
if (Sys.getenv("SWEEP_SENS") == "1") {
  source(file.path(CODE, "hist/sens_helper.R"))
  sens_run(id = "X087", dat = s, y = "midonset", tie = "settlem", ctrl = setdiff(rhs, "settlem"), fe = NULL,
           fit = function(dat, add) kstat(fitter(c(rhs, add), dat)$ct, "settlem"))
}
est <- "logit, default (iid) SEs as in author do-file (fixest feglm)"
ahc <- "Y (peace years + 3 BTSCS cubic splines)"
nt <- sprintf("Table 3 Model 3 (all years 1816-2001, dyad-years of author onset file). No published table accessible (OUP PDF blocked) and no author log -> UNVERIFIED; repro N=%d. Common sample drops %d obs (years<1822 / benchmark missingness). Settlem patch not applied.", n_orig, n_orig - attr(R, "n"))
write_p4(meta, "settlem", "TIE", "KEY", "abstract: 'once neighboring states settle their borders, they are less likely to...experience MIDs'",
         "UNVERIFIED", NA, unname(kstat(rep$ct, "settlem")["b"]), est, n_orig, R, ahc, nt)
write_p4(meta, "jtdem6", "DEMOCRACY", "CONTROL", "Table 3 control; abstract: effects persist 'even after controlling for joint democracy'",
         "UNVERIFIED", NA, unname(kstat(rep$ct, "jtdem6")["b"]), est, n_orig, R, ahc, paste(nt, "bD (lagged joint democracy) is NOT a placebo for this term."))
