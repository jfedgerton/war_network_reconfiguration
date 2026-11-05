# ESTAUDIT 115: Edgerton 2024 JOP. Headline specs (user): LPM, FE logit and the latent-factor model (LFM), all co-equal.
# rerun_115.R ran only the FE logit in the six steps. This copy keeps the identical data prep / training split (author seed 1) /
# 1822-2011 common sample and runs the author's LPM (7_analyze_conflict_data.r l.51):
#   feols(fatal_mid ~ out_group + democ_joint + autoc_joint + major_power_joint + cinc_ratio | ccode1 + ccode2 + dyad + year),
#   SEs as in the author's summary(): cluster = c("year","dyad","ccode1","ccode2")   (same CL object as rerun_115.R)
# for all six steps. FE-logit numbers are re-estimated here too (unchanged spec) as a check against results/115.csv.
# The LFM (brms) is handled in rerun_est_115_lfm.R.
ROOT <- Sys.getenv("SWEEP_DATA"); OUT <- Sys.getenv("SWEEP_OUT"); CODE <- Sys.getenv("SWEEP_CODE")
stopifnot(nzchar(ROOT), nzchar(OUT), nzchar(CODE))
suppressPackageStartupMessages({library(fixest); library(dplyr)})
source(file.path(CODE, "hist/add_history.R")); source(file.path(CODE, "work/est_write_ea1.R"))
WO <- file.path(OUT, "work/115"); dir.create(WO, recursive = TRUE, showWarnings = FALSE)
sink(file.path(WO, "log_est_115.txt"), split = TRUE)
war.data <- readRDS(file.path(ROOT, "_parent", "coop/war_data_replication.rds")); stopifnot(nrow(war.data) == 993975)
CL <- ~year + dyad + ccode1 + ccode2
RHS <- "fatal_mid ~ out_group + democ_joint + autoc_joint + major_power_joint + cinc_ratio"
fe <- "| ccode1 + ccode2 + dyad + year"
fit <- function(add, dat, est) {
  fml <- as.formula(paste(RHS, add, fe))
  m <- if (est == "lpm") feols(fml, data = dat, cluster = CL) else feglm(fml, data = dat, family = "binomial", cluster = CL)
  ct <- coeftable(m); cat("\n----", est, add, "| N =", nobs(m), "|", attr(ct, "type"), "\n"); print(ct["out_group", ])
  c(b = unname(ct["out_group", 1]), se = unname(ct["out_group", 2]), p = unname(ct["out_group", 4]), n = nobs(m))
}
set.seed(1); tr <- sample(nrow(war.data), round(0.75 * nrow(war.data)))
d <- filter(war.data[tr, ], !is.na(fatal_mid))
rep_lpm <- fit("", d, "lpm"); cat("LPM full training sample out_group:", rep_lpm["b"], "\n")
d <- add_history(d, "ccode1", "ccode2", "year")
d <- filter(d, year >= 1822, year <= 2011)
d <- d[complete.cases(d[, c("fatal_mid","out_group","democ_joint","autoc_joint","major_power_joint","cinc_ratio",
                            "hist_dyad_w15","hist_mono_w15","hist_dyad_w26","hist_mono_w26","bm_cinc_w15","bm_democ_w15")]), ]
cat("common sample rows:", nrow(d), " (rerun_115: 582286)\n"); stopifnot(nrow(d) == 582286)
KEEPV <- c("fatal_mid","out_group","democ_joint","autoc_joint","major_power_joint","cinc_ratio","contig","state_count","year",
              "ccode1","ccode2","dyad","hist_dyad_w15","hist_mono_w15","hist_dyad_w26","hist_mono_w26","bm_cinc_w15","bm_democ_w15")
cat("LFM vars missing from data:", setdiff(KEEPV, names(d)), "\n")
saveRDS(d[, intersect(KEEPV, names(d))], file.path(WO, "common_sample_115.rds"))            # for the LFM script
ADD <- list(m0 = "", mD = "+ hist_dyad_w15", mDM = "+ hist_dyad_w15 + hist_mono_w15", m26 = "+ hist_dyad_w26 + hist_mono_w26",
            bC = "+ bm_cinc_w15", bD = "+ bm_democ_w15")
L <- lapply(ADD, fit, dat = d, est = "lpm")
G <- lapply(ADD, fit, dat = d, est = "logit")
stopifnot(abs(G$m0["b"] - 0.858522701419272) < 1e-6)   # = results/115.csv (the mDM check was dropped 2026-10-01: monadic history redefined)
print(rbind(lpm = sapply(L, `[`, "b"), logit = sapply(G, `[`, "b")))
saveRDS(list(L = L, G = G, rep_lpm = rep_lpm), file.path(WO, "est_115.rds"))
write_est_ea1("115", "Edgerton 2024", "out_group [LPM]",
  "co-equal headline specs: LPM with ccode1+ccode2+dyad+year FE; FE logit (same FE); Bayesian latent-factor/crossed-RE logit (brms). This row: LPM, SE clustered year+dyad+ccode1+ccode2",
  "code (7_analyze_conflict_data.r l.51 full_lpm + summary(cluster=c('year','dyad','ccode1','ccode2'))); user instruction (three co-equal headline specs); no PDF in folder",
  "LPM (fixest feols), FE ccode1+ccode2+dyad+year, 4-way clustered SE",
  "N", "LPM headline spec never run in the six steps (only FE logit)", "Y", lapply(L, `[`, c("b","se","p")), unname(L$m0["n"]),
  sprintf("LPM full training-sample repro out_group b=%.5g (SE %.3g). LPM uses all %d common-sample rows (no separation drop; FE logit keeps %d). Same fixest cluster handling as rerun_115.R.",
          rep_lpm["b"], rep_lpm["se"], unname(L$m0["n"]), unname(G$m0["n"])))
sink()
