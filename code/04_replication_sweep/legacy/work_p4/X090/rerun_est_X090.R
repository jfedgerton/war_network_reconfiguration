# ESTAUDIT X090: Choi 2024 PLOS ONE, Table 1 Model 3 (war initiation). Copy of rerun_verify_X090.R (same DyMID 4.02 war DV,
# Chen Model 8 RHS, fatal-MID sample w/o 1968, common sample) with the ESTIMATOR fixed. Article: "the study chooses a rare event
# logistic regression" (c24.txt l.313, l.459 "Chen and I use rare event logit") + year fixed effects (Table 1 'Year Fixed Effects Yes'),
# robust SE (Chen's code: relogit ... _Iyear_*, cluster(ddyad)). All previous steps used ordinary logit + year FE.
# Now every step: relogit_kz() = ML logit with year FE (fixest; years without a war initiation drop out as perfect predictors, as
# in Stata) + King-Zeng bias correction computed with the full year-dummy design; SE cluster(ddyad) x n/(n+k). No tau (full
# directed-dyad population). Caches -> est_cache/.
# VERIFY X090 Choi 2024 PLOS ONE. Copy of rerun_X090.R with ONE change justified by the paper text:
#   DV source. Paper p.4 ("Replication and extended empirical analysis"): war DV built from "Maoz et al.'s dyadic militarized disputes
#   dataset ... dyadic war collection", >=1000 battle deaths, State A initiates a NEW war vs B, joiners excluded.
#   Earlier rebuild used GML (hostlev 5). Now: Maoz Dyadic MID v4.02 (COW site; v3.1 not posted). Joiner rule calibrated on the
#   MID DV only: Chen's midsidea == DyMID4.02 {first dyadic-MID year, sideaa==1, orignata==1, orignatb==1} in 1599/1603 onsets
#   (4 extra) -> same rule applied to war: first year with war==1 in the dispute-dyad, A on side A, both states originators.
# Everything else unchanged: Chen Model 8 RHS, Chen fatal-MID sample (non-missing midsideafatal, year!=1968; N 1,362,097 = Choi N),
# logit with year FE (paper: relogit + year dummies), SE clustered on ddyad. Table 1 Model 3: extdeporig -0.072 (0.045), n.s.
ROOT <- Sys.getenv("SWEEP_DATA"); OUT <- Sys.getenv("SWEEP_OUT"); CODE <- Sys.getenv("SWEEP_CODE")
stopifnot(nzchar(ROOT), nzchar(OUT), nzchar(CODE))
source(file.path(CODE, "vb2/vcommon.R"))
source(file.path(CODE, "work/relogit_kz_ea1.R")); source(file.path(CODE, "work/est_write_ea1.R"))
suppressPackageStartupMessages(library(haven))
WO <- file.path(OUT, "work_p4/X090"); dir.create(WO, recursive = TRUE, showWarnings = FALSE)  # port: outputs
W <- file.path(ROOT, "work_p4/X090"); sink(file.path(WO, "log_est_X090.txt"), split = TRUE)
V <- c("ccode1", "ccode2", "year", "ddyad", "midsideafatal", "extdeporig", "dyaddep0", "dpcapshare", "comdpno", "dyaddp",
       "c2majpow", "c1dem", "c2dem", "cont", "dist", "pcyrmidsideafatal1", "pcyrmidsideafatal2", "pcyrmidsideafatal3")
d <- as.data.frame(zap_labels(read_dta(file.path(W, "replication.dta"), col_select = all_of(V))))
stopifnot(!anyDuplicated(d[, c("ccode1", "ccode2", "year")]))

y <- read.csv(file.path(W, "dymid/dyadic_mid_4.02.csv"))
y <- y[y$sideaa == 1 & y$orignata == 1 & y$orignatb == 1 & y$war == 1, c("disno", "statea", "stateb", "year")]
w <- y %>% group_by(disno, statea, stateb) %>% summarise(year = min(year), .groups = "drop") %>%
  distinct(ccode1 = statea, ccode2 = stateb, year) %>% mutate(war = 1)
cat("DyMID4.02 directed war initiations (all years):", nrow(w), "; 1951-2010:", sum(w$year >= 1951 & w$year <= 2010), "\n")
d <- d %>% left_join(w, by = c("ccode1", "ccode2", "year")) %>% mutate(war = coalesce(war, 0))
d <- d[!is.na(d$midsideafatal) & d$year != 1968, ]
cat("war sample rows:", nrow(d), "; war initiations in sample:", sum(d$war), "\n")

f0 <- war ~ extdeporig + dyaddep0 + dpcapshare + comdpno + dyaddp + c2majpow + c1dem + c2dem + cont + dist +
  pcyrmidsideafatal1 + pcyrmidsideafatal2 + pcyrmidsideafatal3
terms <- c("extdeporig", "extdeporig__unc")
FF <- function(f, data) relogit_kz(f, data, cluster = "ddyad", fe = "year")
CD <- file.path(WO, "est_cache"); dir.create(CD, recursive = TRUE, showWarnings = FALSE)
rep <- vfit(f0, d[complete.cases(d[, c(all.vars(f0), "ddyad")]), ], terms, FF, file.path(CD, "e_repro.rds"))
ds <- vcommon_sample(d, f0, "ddyad", "ccode1", "ccode2", "year"); rm(d); gc()
C <- function(k) file.path(CD, paste0("e_", k, ".rds"))
M <- list(m0  = vfit(f0, ds, terms, FF, C("m0")),
          mD  = vfit(update(f0, . ~ . + hist_dyad_w15), ds, terms, FF, C("mD")),
          mDM = vfit(update(f0, . ~ . + hist_dyad_w15 + hist_mono_w15), ds, terms, FF, C("mDM")),
          m26 = vfit(update(f0, . ~ . + hist_dyad_w26 + hist_mono_w26), ds, terms, FF, C("m26")),
          bC  = vfit(update(f0, . ~ . + bm_cinc_w15), ds, terms, FF, C("bC")),
          bD  = vfit(update(f0, . ~ . + bm_democ_w15), ds, terms, FF, C("bD")))
stopifnot(all(sapply(M, function(m) m$n) == M$m0$n), M$m0$n == 472614)   # = results_verify n_used

# ---- sensitivity hook: runs only under run_sensitivity.R (SWEEP_SENS=1), never in the main sweep ----
if (Sys.getenv("SWEEP_SENS") == "1") {
  source(file.path(CODE, "hist/sens_helper.R"))
  # the year-FE logit drops years without variation in the outcome (same rows in every step); pass exactly the sample the model uses
  sens_ds <- ds[ave(ds$war, ds$year, FUN = function(z) as.numeric(min(z) != max(z))) == 1, ]
  sens_run(id = "X090", dat = sens_ds, y = "war", tie = "extdeporig",
           ctrl = c("dyaddep0", "dpcapshare", "comdpno", "dyaddp", "c2majpow", "c1dem", "c2dem", "cont", "dist",
                    "pcyrmidsideafatal1", "pcyrmidsideafatal2", "pcyrmidsideafatal3"), fe = "year",
           fit = function(dat, add) { ct <- FF(if (length(add)) update(f0, as.formula(paste(". ~ . +", paste(add, collapse = " + ")))) else f0, dat)
             c(b = unname(ct["extdeporig", "b"]), se = unname(ct["extdeporig", "se"]), p = unname(ct["extdeporig", "p"])) })
}

RK <- lapply(M, function(m) m$k$extdeporig); RU <- lapply(M, function(m) m$k$extdeporig__unc)
print(rbind(KZ = sapply(RK, `[`, "b"), logit = sapply(RU, `[`, "b")))
pc <- function(a, b) 100 * (a - b) / abs(b)
write_est_ea1("X090", "Choi 2024", "extdeporig",
  "rare-events logit (relogit, King-Zeng) with year fixed effects, cubic peace years, SE clustered on directed dyad; 1951-2010",
  "pdf (Choi (2024) PLOSONE: 'chooses a rare event logistic regression'; Table 1 'Year Fixed Effects: Yes', robust SE) + Chen's code (relogit ..., cluster(ddyad))",
  "relogit: ML logit + year FE (fixest) + King-Zeng bias correction with full year-dummy design (work/relogit_kz_ea1.R), SE cluster ddyad x n/(n+k)",
  "N", "all six steps (rerun_verify_X090.R) used ordinary logit with year FE; King-Zeng rare-events bias correction not applied",
  "Y", RK, M$m0$n,
  sprintf(paste0("REBUILT_MATCH war DV (DyMID 4.02) from rerun_verify_X090.R unchanged; 31 war initiations -> KZ correction matters more here. ",
    "KZ b0=%.6g (p=%.3g); uncorrected logit b0=%.6g (p=%.3g) = results_verify; uncorrected mDM %.6g (p=%.3g, pct_dm %.1f). ",
    "Full-sample relogit repro %.6g (pub -0.072, n.s.). Years with no war drop out of the FE fit (as Stata)."),
    RK$m0["b"], RK$m0["p"], RU$m0["b"], RU$m0["p"], RU$mDM["b"], RU$mDM["p"], pc(RU$mDM["b"], RU$m0["b"]), rep$k$extdeporig["b"]))
sink()
