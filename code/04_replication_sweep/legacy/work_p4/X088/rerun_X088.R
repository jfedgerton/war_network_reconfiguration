# X088 Atkinson, Owsiak, Jackson & Buechler 2024 II "Sequencing the Steps to War". Dataverse doi:10.7910/DVN/7IJTNQ
# (replication.tab + do-file + Stata log). Estimator: Stata -sem- path analysis, observed variables only, recursive,
# uncorrelated structural errors -> ML coefficients of each equation == equation-by-equation OLS on the SEM listwise sample.
# We reproduce the war_onset equation (LPM) by OLS. Key tie term (decided before estimating): new_atop_ally
# ("secure defensive allies" step) direct path to war_onset; CONTROL: jtdem (joint democracy) in the same equation.
# Scope caveat: the paper's headline claim is the modal SEQUENCE (indirect paths), not one coefficient.
ROOT <- Sys.getenv("SWEEP_DATA"); OUT <- Sys.getenv("SWEEP_OUT"); CODE <- Sys.getenv("SWEEP_CODE")
stopifnot(nzchar(ROOT), nzchar(OUT), nzchar(CODE))
source(file.path(CODE, "work_p4/common_p4_x8.R"))
WO <- file.path(OUT, "work_p4/X088"); dir.create(WO, recursive = TRUE, showWarnings = FALSE)  # port: outputs
W <- file.path(ROOT, "work_p4/X088")
sink(file.path(WO, "log_X088.txt"), split = TRUE)

d <- data.table::fread(file.path(W, "replication.tab"), data.table = FALSE)
SEMV <- c("jtdem","new_icow","rel_pow_rat","war_onset","new_spatial_riv","icb_onset_any","majinvolve","new_atop_ally",
          "new_counter_alliance","new_ghr_war_adj_darmsrace","lag_total_outside_atop_ally")
d <- d[complete.cases(d[, SEMV]), ]
cat("SEM listwise N:", nrow(d), "(Stata log: 539,068)\n"); stopifnot(nrow(d) == 539068)
d$ccode1 <- d$dyad %/% 1000; d$ccode2 <- d$dyad %% 1000
cat("directed? share of dyads with ccode1>ccode2:", mean(d$ccode1 > d$ccode2), "\n")
f0 <- "war_onset ~ new_icow + new_spatial_riv + icb_onset_any + new_atop_ally + new_counter_alliance + new_ghr_war_adj_darmsrace + jtdem + majinvolve + rel_pow_rat"
KEY <- c("new_atop_ally", "jtdem")
fit <- function(extra, s) fit_fx(as.formula(paste(f0, extra)), s, KEY, "ols", "iid")

r_full <- fit("", d)
B_PUB <- c(new_atop_ally = .0010748, jtdem = -.0004668)
n_orig <- r_full[[1]]$n

d <- add_history(d, "ccode1", "ccode2", "year")
s <- subset(d, year >= 1822 & year <= 2011)
s <- s[complete.cases(s[, HISTV]), ]
cat("common N", nrow(s), "dropped", n_orig - nrow(s), "\n")
M <- run_six(fit, s)

meta <- list(xid = "X088", short_cite = "Atkinson, Owsiak, Jackson & Buechler 2024", journal = "II", year = 2024, doi = "10.1080/03050629.2024.2378779")
ev <- "abstract: 'war ... rises as ... states ... secure defensive allies' (steps-to-war sequence)"
for (k in KEY) {
  rep <- if (abs(r_full[[k]]$b - B_PUB[[k]]) / abs(B_PUB[[k]]) <= 0.10 && r_full[[k]]$p < .05) "REPRODUCED"  # both published paths p<.001 else "NOT_REPRODUCED"
  nt <- sprintf(paste0("SEM path analysis (Stata sem, ML); war_onset equation reproduced exactly as OLS LPM (recursive observed-variable SEM). ",
    "b_pub from authors' log (Table 1 SEM): %s = %.7f; repro %.7f (p=%.3g). Borderline scope: paper's claim is the modal step sequence, key tie = defensive-ally step's direct path. ",
    "Data 1816-2022 dyad-years (dyad = ccode1*1000+ccode2); common sample 1822-2011 N=%d (orig %d). SE = OLS iid (SEM OIM analogue). ",
    "Authors include icb_onset_any and rivalry (new_spatial_riv) as steps but no lagged MID-count history."),
    k, B_PUB[[k]], r_full[[k]]$b, r_full[[k]]$p, M$m0[[k]]$n, n_orig)
  write_p4(meta, k, if (k == "jtdem") "DEMOCRACY" else "TIE", if (k == "jtdem") "CONTROL" else "KEY", ev, rep, B_PUB[[k]], r_full[[k]]$b,
           "SEM path (war eq. = OLS LPM, iid SE)", n_orig, M, "N (rivalry/ICB steps, no MID history)", nt)
}
sink()
