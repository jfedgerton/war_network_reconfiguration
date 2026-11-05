# Rerun 32: Bas & Schub 2016 ISQ -- Table 1 Model 1:
#   logit mid msod relcap cont jdem2 ally pceyrs pceyrs2 pceyrs3 if ongomid==0, cluster(dyad)
ROOT <- Sys.getenv("SWEEP_DATA"); OUT <- Sys.getenv("SWEEP_OUT"); CODE <- Sys.getenv("SWEEP_CODE")
stopifnot(nzchar(ROOT), nzchar(OUT), nzchar(CODE))
source(file.path(CODE, "work/rerun_common_g7.R"))
suppressPackageStartupMessages(library(haven))
W <- file.path(ROOT, "work/32")
WO <- file.path(OUT, "work/32"); dir.create(WO, recursive = TRUE, showWarnings = FALSE)  # port: outputs
sink(file.path(WO, "log_32.txt"), split = TRUE)

# 1. Authors' analysis data (Dataverse 10.7910/DVN/4W1JHE, MOmain_BasSchub.dta original format)
d <- as.data.frame(zap_labels(read_dta(file.path(W, "MOmain.dta"))))
stopifnot(nrow(d) == 25387, all(c("ccode1","ccode2","year","dyad","mid","msod","ongomid") %in% names(d)))
d <- d[d$ongomid == 0, ]
BASE <- "mid ~ msod + relcap + cont + jdem2 + ally + pceyrs + pceyrs2 + pceyrs3"
MV <- c("mid","msod","relcap","cont","jdem2","ally","pceyrs","pceyrs2","pceyrs3","dyad")
fitter <- function(add, dat) fit_glm_tie(as.formula(paste0(BASE, add)), dat, binomial(), "msod", cluster = "dyad")

# 2. Reproduce on the full published sample
cat("\n==== Reproduction: Table 1 Model 1 (full sample) ====\n")
d_full <- d[complete.cases(d[, MV]), ]
rp <- fitter("", d_full)
b_pub <- NA   # paper PDF not retrievable (DASH behind captcha); no author log in deposit
n_orig <- rp$n

# 3. Merge history; restrict to 1822-2011 (sample is 1816-1923, so 1816-1821 dropped)
d <- add_history(d, "ccode1", "ccode2", "year")
d_all <- d
d <- restrict_sample(d, MV)

# 4. Six models on one sample
R <- run_six(fitter, d)
SENS <- sens_hist_only(fitter, d_all, MV); rm(d_all)

# ---- sensitivity hook: runs only under run_sensitivity.R (SWEEP_SENS=1), never in the main sweep ----
if (Sys.getenv("SWEEP_SENS") == "1") {
  source(file.path(CODE, "hist/sens_helper.R"))
  sens_run(id = "32", dat = d, y = "mid", tie = "msod",
           ctrl = c("relcap", "cont", "jdem2", "ally", "pceyrs", "pceyrs2", "pceyrs3"), fe = NULL,
           fit = function(dat, add) { r <- fitter(if (length(add)) paste0(" + ", paste(add, collapse = " + ")) else "", dat); c(b = r$b, se = r$se, p = r$p) })
}

# 5. Results row
write_row(32, "Bas & Schub 2016", "UNVERIFIED", b_pub, rp$b, "logit, SE clustered on dyad (glm + sandwich::vcovCL)", "msod",
  n_orig, R$m0$n, R, ahc = "Y",
  notes = sprintf("Table 1 M1 (msod, 1816-1923, ongoing MIDs dropped); author controls cubic peace years (Y=peace years only). b_pub not obtainable (paper behind captcha, no log); code run verbatim from author .do. Dropped %d obs total vs orig N=%d; %s. Sensitivity %s.", n_orig - R$m0$n, n_orig, RS_INFO, SENS))
sink()
