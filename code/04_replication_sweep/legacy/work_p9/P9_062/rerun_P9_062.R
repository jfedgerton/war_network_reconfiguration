# P9_062 Hegre, Bernhard & Teorell 2020 JCR 64(1) "Civil Society and the Democratic Peace" (10.1177/0022002719850620). DEMOCRACY.
# HEADLINE (supplement read 2026-10-11, "Hegre et al. (2020) JCR supplement.pdf" Table SF-3 col 3 "Social"; abstract: social
#   accountability "the strongest and most consistent predictor of nonbelligerence in dyads"): Social accountability, int.
#   (cspart_ij = cspart_i * cspart_j, built in AssembleData2) = -6.092***, t = -7.72 -> SE = 6.092/7.72; N 541,560; AIC 5854.5.
#   Legend "* p < 0.05, ** p < 0.01, *** p < 0.001"; "Logit model, with robust standard errors" (do-file: cluster(dyadid)).
# Step 1: prep_P9_062.R rebuilds the model-7 sample from the authors' inputs (validated against their log; saves the model data).
# Step 2 (here): full-sample check against SF-3, then the six steps. DV = fatal MID in year t (contemporaneous) -> history keyed
#   to year. Undirected dyads (statea < stateb). Author controls include peace years (py1) -> B_duration.
CODE <- Sys.getenv("SWEEP_CODE"); stopifnot(nzchar(CODE))
source(file.path(CODE, "work_p9/P9_062/prep_P9_062.R"))     # writes <OUT>/work_p9/P9_062/P9_062_model7_data.rds
rm(list = intersect(c("d", "v", "m", "um", "ds"), ls())); invisible(gc())
source(file.path(CODE, "work_p8/common_b3.R"))
W <- file.path(OUT, "work_p9/P9_062")
ds0 <- readRDS(file.path(W, "P9_062_model7_data.rds"))
sink(file.path(W, "log_P9_062.txt"), split = TRUE)
stopifnot(all(ds0$statea < ds0$stateb), !anyDuplicated(ds0[, c("statea", "stateb", "year")]))
rhs0 <- c("cspart_i", "cspart_j", "cspart_ij", "lncapi", "lncapj", "lnpopi", "lnpopj", "dircont", "lndstab", "lnNt", "py1")
vars <- c("mzfatal_new", rhs0, "dyadid")
TERMS <- "cspart_ij"
fitter <- function(r, data) fit_bin("mzfatal_new", r, data, "logit", "~dyadid", terms = TERMS)
full <- fitter(rhs0, ds0[complete.cases(ds0[, vars]), ])
cat("FULL-SAMPLE model 7: b", full$cspart_ij[["b"]], "se", full$cspart_ij[["se"]], "N", attr(full, "n"), "(SF-3 col 3: -6.092, t -7.72, N 541,560)\n")
stopifnot(attr(full, "n") == 541560, abs(full$cspart_ij[["b"]] - (-6.092)) < 0.001)
ds <- common_sample(ds0, vars, "statea", "stateb", "year")
M <- run_six(rhs0, ds, fitter)
if (Sys.getenv("SWEEP_SENS") == "1") {   # sensitivity hook (split: product term); runs only under run_sensitivity.R
  source(file.path(CODE, "hist/sens_helper.R"))
  sens_split(id = "P9_062", dat = ds, fit = function(dat, add) fitter(c(rhs0, add), dat)[["cspart_ij"]])
}
meta <- list(pid = "P9_062", cite = "Hegre et al. 2020", journal = "JCR", year = 2020, doi = "10.1177/0022002719850620", family = "DEMOCRACY")
dir.create(file.path(OUT, "results_p9"), recursive = TRUE, showWarnings = FALSE)
write_p8(meta, "cspart_ij",
  "Supplement Table SF-3 col 3: Social accountability, int. -6.092*** (t -7.72); abstract: social accountability strongest predictor of dyadic nonbelligerence",
  "Supplement Table SF-3 col 3 (Social; fatal MID, undirected dyad-years 1900-2010; = article Figure 3)",
  "logit, SE clustered on dyad, peace-years decay (py1)", "fixest feglm logit, cluster dyadid, all six steps",
  -6.092, 6.092 / 7.72, 541560, TRUE, FALSE, M, "Y (peace years py1)",
  sprintf("Sample rebuilt from the authors' V-Dem v6.2 inputs (prep_P9_062.R, log counts matched). Full-sample b=%.6f se=%.6f N=%d (published N 541,560). Common sample N=%d. Published SE = |b|/t.",
          full$cspart_ij[["b"]], full$cspart_ij[["se"]], attr(full, "n"), nrow(ds)))
src <- file.path(OUT, "results_p8", "P9_062__cspart_ij.csv"); stopifnot(file.exists(src), file.rename(src, file.path(OUT, "results_p9", basename(src))))
sink()
