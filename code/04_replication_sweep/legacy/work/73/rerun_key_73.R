# KEY rerun 73 (key term = ltXext = logtrade_strategic x externality, Table II M1 interaction; pub 0.189).
# Copy of rerun_73.R with TIE switched to the interaction; results saved to key_73.rds (results/73.csv untouched).
# Rerun 73: Zeng 2024 JPR "Microchips and sneakers" -- strategic trade x relative externality -> hostile MID onset (t+1)
# Headline: Table II Model 1 (small power gap sample, commitcut==1): logit, ML SE (authors' 02-Estimation.R fits1).
# Headline interaction is logtrade_strategic x externality (0.189**); per protocol the tie = constituent term
# logtrade_strategic (pub -0.029, CI -0.099..0.040). Interaction term tracked in notes.
ROOT <- Sys.getenv("SWEEP_DATA"); OUT <- Sys.getenv("SWEEP_OUT"); CODE <- Sys.getenv("SWEEP_CODE")
stopifnot(nzchar(ROOT), nzchar(OUT), nzchar(CODE))
source(file.path(CODE, "work/rerun_helpers_s68.R"))
WO <- file.path(OUT, "work/73"); dir.create(WO, recursive = TRUE, showWarnings = FALSE)  # port: outputs
W <- file.path(ROOT, "work/73"); sink(file.path(WO, "log_key_73.txt"), split = TRUE)
suppressPackageStartupMessages({library(dplyr); library(tidyr); library(zoo)})

# 1. Authors' data prep, copied from 02-Estimation.R (RcppRoll::roll_mean(na.rm=TRUE) -> data.table::frollmean(na.rm=TRUE), RcppRoll not installable)
load(file.path(W, "DA-RT_JPR_2022/Data/commit_dd.RData"))
# lead check (fix 2026-10-10): midlead1_onset[t] == midonset[t+1] on consecutive directed-dyad records (before columns are dropped)
lc <- commit_dat_dd[order(commit_dat_dd$ccode1, commit_dat_dd$ccode2, commit_dat_dd$year), c("ccode1", "ccode2", "year", "midonset", "midlead1_onset")]
nl <- nrow(lc); cs <- which(lc$ccode1[-1] == lc$ccode1[-nl] & lc$ccode2[-1] == lc$ccode2[-nl] & lc$year[-1] == lc$year[-nl] + 1)
la <- lc$midlead1_onset[cs]; lb <- lc$midonset[cs + 1]; ok <- !is.na(la) & !is.na(lb)
stopifnot(length(cs) > 1000, sum(ok) > 1000, all(la[ok] == lb[ok]))
cat("lead check: midlead1_onset[t] == midonset[t+1] on", sum(ok), "consecutive dyad pairs with both non-missing\n"); rm(lc, cs, la, lb, ok, nl)
commit_dat_dd <- commit_dat_dd[, c("ccode1","ccode2","year","midlead1_onset","externality_manufacexp_1","externality_manufacexp_2","tradetotal_1","tradetotal_2","cinc_en_md_1","cinc_en_md_2","defense","neutrality","nonaggression","dyadtrade_strategic","mindist","polity2_1","polity2_2","conttype","major1","major2","peaceyear")]; invisible(gc())
rm_na <- function(x) { r <- data.table::frollmean(x, 5, na.rm = TRUE, align = "right", algo = "exact"); r[is.nan(r)] <- NA; r }  # = RcppRoll::roll_mean(5, align="right", fill=NA, na.rm=TRUE)
commit_dat <- commit_dat_dd %>%
  mutate(externality_1 = externality_manufacexp_1, externality_2 = externality_manufacexp_2) %>%
  mutate(externality_1 = ifelse(externality_1 == 0 & tradetotal_1 == 0, NA, externality_1),
         externality_2 = ifelse(externality_2 == 0 & tradetotal_2 == 0, NA, externality_2)) %>%
  mutate(externality = (externality_1 - externality_2) / (externality_1 + externality_2)) %>%
  arrange(ccode1, ccode2, year) %>% group_by(ccode1, ccode2) %>%
  mutate(externality = if (n() >= 5) rm_na(externality) else NA_real_) %>% ungroup() %>%
  filter(year >= 1962 & year <= 2009)
rm(commit_dat_dd); invisible(gc())
commit_dat <- commit_dat %>%
  mutate(cinc1 = cinc_en_md_1, cinc2 = cinc_en_md_2) %>%
  filter(cinc1 < cinc2) %>%
  mutate(powerratio = cinc1 / (cinc1 + cinc2),
         ally = as.numeric(defense + neutrality + nonaggression > 0),
         nationaltrade_1 = log(tradetotal_1 + 1), nationaltrade_2 = log(tradetotal_2 + 1),
         logtrade_strategic = log(dyadtrade_strategic + 1),
         distance = log(mindist + 1),
         jointdemo = (polity2_1 > 5) & (polity2_2 > 5),
         coldwar = as.factor(year < 1992),
         dyad = ifelse(ccode1 < ccode2, ccode1 * 1000 + ccode2, ccode2 * 1000 + ccode1)) %>%
  replace_na(list(ally = 0)) %>% mutate(ally = as.factor(ally))
relv_dat <- commit_dat %>% filter(conttype < 6 | (major1 + major2 > 0)) %>%
  group_by(year) %>% mutate(commitcut = as.numeric(powerratio > quantile(powerratio, 0.75, na.rm = TRUE))) %>% ungroup()
d <- as.data.frame(relv_dat %>% filter(commitcut == 1))
d$pY2 <- d$peaceyear^2; d$pY3 <- d$peaceyear^3; d$ltXext <- d$logtrade_strategic * d$externality
rm(commit_dat, relv_dat); invisible(gc())

# 2. Reproduce Table II M1
TIE <- "ltXext"
f0 <- midlead1_onset ~ logtrade_strategic + externality + ltXext + nationaltrade_1 + nationaltrade_2 +
  distance + ally + jointdemo + coldwar + peaceyear + pY2 + pY3
b_pub <- 0.189
rep <- fit_glm(f0, d, TIE, label = "repro")
repI <- fit_glm(f0, d, "ltXext", label = "repro interaction")
n_orig <- rep$n
cat("b_pub", b_pub, "b_repro", rep$b, "N", n_orig, "(pub 12521); interaction pub 0.189 repro", repI$b, "\n")
# pub CI (-0.099, 0.040) -> not significant; interaction significant at .01
status <- if (abs(rep$b - b_pub) / abs(b_pub) <= 0.10 && rep$p < .05) "Y" else "N"

# 3. Merge history (one row per dyad-year, ccode1 weaker). DV is onset at t+1, so history and benchmark windows are keyed on the
#    outcome year t+1 (fix 2026-10-10; previously the data year t). Years 1962-2009: the 1822-2011 filter (row year) is unaffected.
d$outcome_year <- d$year + 1
d <- add_history(d, "ccode1", "ccode2", "outcome_year")
s <- d[d$year >= 1822 & d$year <= 2011, ]
s <- s[complete.cases(s[, c(all.vars(f0), "hist_dyad_w15", "hist_mono_w15", "hist_dyad_w26", "hist_mono_w26")]), ]
cat("Estimation sample", nrow(s), "; dropped vs original", n_orig - nrow(s), "\n")

# 4. Models (tie + interaction each time)
run <- function(f, dat, lab) fit_glm(f, dat, TIE, label = lab)
m0  <- run(f0, s, "m0")
mD  <- run(update(f0, . ~ . + hist_dyad_w15), s, "mD")
mDM <- run(update(f0, . ~ . + hist_dyad_w15 + hist_mono_w15), s, "mDM")
m26 <- run(update(f0, . ~ . + hist_dyad_w26 + hist_mono_w26), s, "m26")
sC <- s[!is.na(s$bm_cinc_w15), ]; sD <- s[!is.na(s$bm_democ_w15), ]
cat("bC sample", nrow(sC), " bD sample", nrow(sD), "\n")
b0C <- if (nrow(sC) < nrow(s)) fit_glm(f0, sC, TIE, label = "m0 on bC sample")$b else NULL
bC  <- run(update(f0, . ~ . + bm_cinc_w15), sC, "bC")
b0D <- if (nrow(sD) < nrow(s)) fit_glm(f0, sD, TIE, label = "m0 on bD sample")$b else NULL
bD  <- run(update(f0, . ~ . + bm_democ_w15), sD, "bD")
# ---- sensitivity hook (split only): runs only under run_sensitivity.R (SWEEP_SENS=1), never in the main sweep ----
if (Sys.getenv("SWEEP_SENS") == "1") {
  source(file.path(CODE, "hist/sens_helper.R"))
  sens_split(id = "73", dat = s,
             fit = function(dat, add) { r <- run(if (length(add)) update(f0, as.formula(paste(". ~ . +", paste(add, collapse = " + ")))) else f0, dat, "sens")
               c(b = r$b, se = r$se, p = r$p) })
}

saveRDS(list(rep = rep, b_pub = b_pub, status = status, n_orig = n_orig, R = list(m0 = m0, mD = mD, mDM = mDM, m26 = m26, bC = bC, bD = bD), b0C = b0C, b0D = b0D), file.path(WO, "key_73.rds"))
print(sapply(list(m0 = m0, mD = mD, mDM = mDM, m26 = m26, bC = bC, bD = bD), unlist)); cat("b0C", if (is.null(b0C)) NA else b0C, "b0D", if (is.null(b0D)) NA else b0D, "status", status, "\n")
# port: one-row wide copy of the same numbers (b_hist = mDM; benchmarks vs m0 on the benchmark sample when it differs, as in
# the chain step that filled results_key/73.csv) so run_one.R can extract by column; no new estimation
wide <- data.frame(id = "73", key_term = TIE, b0 = m0$b, p0 = m0$p, b_dm = mDM$b, p_dm = mDM$p, b_26 = m26$b, p_26 = m26$p,
                   pct_bC = pct(bC$b, if (is.null(b0C)) m0$b else b0C), pct_bD = pct(bD$b, if (is.null(b0D)) m0$b else b0D), n_used = m0$n)
write.csv(wide, file.path(WO, "key_73_wide.csv"), row.names = FALSE)
sink()
