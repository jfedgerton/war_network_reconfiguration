# Rerun 91: Kang 2015 CMPS -- Table 1 col 3:
#   logit cwinit c1allcap c2allcap alleach kgdriv relcap jodem sreg condum peaceyrs peaceyrsq peaceyrcube, vce(cluster ddyad)
# Headline tie: c1allcap (initiator's allies' capability). Published (author log): 1.481363, N=164,391.
ROOT <- Sys.getenv("SWEEP_DATA"); OUT <- Sys.getenv("SWEEP_OUT"); CODE <- Sys.getenv("SWEEP_CODE")
stopifnot(nzchar(ROOT), nzchar(OUT), nzchar(CODE))
source(file.path(CODE, "work/rerun_common_g7.R"))
suppressPackageStartupMessages(library(haven))
W <- file.path(ROOT, "work/91")
WO <- file.path(OUT, "work/91"); dir.create(WO, recursive = TRUE, showWarnings = FALSE)
sink(file.path(WO, "log_91.txt"), split = TRUE)

# 1. Authors' analysis data (Dataverse 10.7910/DVN/AMOREX, original .dta)
d <- as.data.frame(zap_labels(read_dta(file.path(W, "kang.dta"))))
stopifnot(nrow(d) == 167811, all(c("ddyad","ccode1","ccode2","year","cwinit","c1allcap") %in% names(d)))
BASE <- "cwinit ~ c1allcap + c2allcap + alleach + kgdriv + relcap + jodem + sreg + condum + peaceyrs + peaceyrsq + peaceyrcube"
MV <- c("cwinit","c1allcap","c2allcap","alleach","kgdriv","relcap","jodem","sreg","condum","peaceyrs","peaceyrsq","peaceyrcube","ddyad")
fitter <- function(add, dat) fit_glm_tie(as.formula(paste0(BASE, add)), dat, binomial(), "c1allcap", cluster = "ddyad")

# 2. Reproduce
cat("\n==== Reproduction: Table 1 col 3 (full sample) ====\n")
d_full <- d[complete.cases(d[, MV]), ]
rp <- fitter("", d_full)
b_pub <- 1.481363
stopifnot(rp$n == 164391)
repro <- if (abs(rp$b - b_pub) / abs(b_pub) <= 0.10 && (rp$p < .05) == TRUE) "REPRODUCED" else "NOT_REPRODUCED"
cat("\nb_pub =", b_pub, " b_repro =", rp$b, " status:", repro, "\n")
n_orig <- rp$n

# 3. History merge (directed dyad: ccode1 = potential initiator) and 1822-2011 restriction
d <- add_history(d, "ccode1", "ccode2", "year")
d_all <- d
d <- restrict_sample(d, MV)

# 4. Six models
R <- run_six(fitter, d)
SENS <- sens_hist_only(fitter, d_all, MV); rm(d_all)

# ---- sensitivity hook: runs only under run_sensitivity.R (SWEEP_SENS=1), never in the main sweep ----
if (Sys.getenv("SWEEP_SENS") == "1") {
  source(file.path(CODE, "hist/sens_helper.R"))
  sens_run(id = "91", dat = d, y = "cwinit", tie = "c1allcap",
           ctrl = c("c2allcap", "alleach", "kgdriv", "relcap", "jodem", "sreg", "condum", "peaceyrs", "peaceyrsq", "peaceyrcube"), fe = NULL,
           fit = function(dat, add) { r <- fitter(if (length(add)) paste0(" + ", paste(add, collapse = " + ")) else "", dat); c(b = r$b, se = r$se, p = r$p) })
}

# 5. Results
write_row(91, "Kang 2015", repro, b_pub, rp$b, "logit, SE clustered on directed dyad (glm + sandwich::vcovCL)", "c1allcap",
  n_orig, R$m0$n, R, ahc = "Y",
  notes = sprintf("Table 1 col 3; b_pub from author Stata log. Author controls cubic peace years + kgdriv rivalry (Y). Directed dyads 1816-2001; dropped %d obs total vs orig; %s. Sensitivity %s.", n_orig - R$m0$n, RS_INFO, SENS))
sink()
