# Rerun 14: Morrow 2017 JOP -- Fig. 2 CEM-matched logit, initiators WITHOUT a defensive ally (pchalally==0),
# full period, with controls:
#   cem jdem(#0) pchaloff(#0) pchalneu(#0) contig(#0) c_majpow(#0) t_majpow(#0) rivalyears(#0) midcount(#0) ptargoff(#0) ptargneu(#0), tr(ptargdef)
#   logit dispute ptargdef jdem pchaloff pchalneu contig c_majpow t_majpow rivalyears midcount ptargoff ptargneu [iweight=cem_weights]
ROOT <- Sys.getenv("SWEEP_DATA"); OUT <- Sys.getenv("SWEEP_OUT"); CODE <- Sys.getenv("SWEEP_CODE")
stopifnot(nzchar(ROOT), nzchar(OUT), nzchar(CODE))
source(file.path(CODE, "work/rerun_common_b8.R"))
library(haven)
WO <- file.path(OUT, "work/14"); dir.create(WO, recursive = TRUE, showWarnings = FALSE)  # port: outputs
sink(file.path(WO, "log_14.txt"), split = TRUE)

# 1. Load authors' data (Leeds-Johnson directed dyad-years; challenger/target COW codes)
CEMV <- c("jdem","pchaloff","pchalneu","contig","c_majpow","t_majpow","rivalyears","midcount","ptargoff","ptargneu")
v <- c("challenger","target","year","dispute","ptargdef","pchalally", CEMV)
d0 <- as.data.frame(zap_labels(read_dta(file.path(ROOT, "work/14/LeedsJohnsonJOPrep.dta"), col_select = all_of(v))))
cat("rows:", nrow(d0), " years:", range(d0$year), "\n")

# CEM with (#0) on every variable == exact matching on the full covariate profile.
# Weights (Iacus-King-Porro): treated = 1; control in matched stratum s = (mT_s/mC_s) * (mC/mT); unmatched = 0.
cem_exact <- function(dat, treat) {
  key <- do.call(paste, c(dat[, CEMV], sep = "|"))
  ok <- complete.cases(dat[, c(CEMV, treat)])
  s <- data.frame(key = key, t = dat[[treat]], ok = ok)
  tab <- s[ok, ] %>% group_by(key) %>% summarise(nT = sum(t == 1), nC = sum(t == 0), .groups = "drop") %>% filter(nT > 0, nC > 0)
  s <- left_join(s, tab, by = "key")
  matched <- s$ok & !is.na(s$nT)
  mT <- sum(matched & s$t == 1); mC <- sum(matched & s$t == 0)
  w <- ifelse(!matched, 0, ifelse(s$t == 1, 1, (s$nT / s$nC) * (mC / mT)))
  cat("CEM exact: matched treated", mT, " matched controls", mC, " strata", nrow(tab), "\n")
  w
}

d <- d0[d0$pchalally == 0, ]
d1 <- d0[d0$pchalally == 1, ]; rm(d0); gc(verbose = FALSE)   # memory
d$w <- cem_exact(d, "ptargdef"); d <- d[d$w > 0, ]   # only matched obs enter the iweight logit
# NOTE: the file's stored cem_weights_1..4 come from undocumented full-sample CEM runs (matched both pchalally groups,
# 3.6k-11k strata => other covariates); they do not correspond to the do-file's cem calls, so they are not used.

TIE <- "ptargdef"
RHS <- paste("dispute ~ ptargdef +", paste(CEMV, collapse = " + "))
# iweight logit, conventional (non-robust) ML SEs
fitter <- function(add, dat) {
  dat <- dat[dat$w > 0, ]
  f <- as.formula(paste(RHS, add))
  # fixest feglm with weights + iid vcov == Stata logit [iweight] default (inverse weighted information) SEs
  m <- fixest::feglm(f, data = dat, family = binomial(), weights = ~w, vcov = "iid", notes = FALSE, mem.clean = TRUE)
  stopifnot(m$convStatus)
  ct <- fixest::coeftable(m)
  cat("\n---- ", paste(deparse(f, width.cutoff = 500), collapse = ""), " | N =", nobs(m), "\n"); print(signif(ct, 5))
  out <- list(b = ct[TIE, 1], se = ct[TIE, 2], p = ct[TIE, 4], n = nobs(m)); rm(m); gc(verbose = FALSE); out
}

# 2. Reproduce. Fig. 2 reports risk ratios graphically; published coefficient not obtainable (paywalled, no OA) and
#    no author log -> UNVERIFIED. Also print the authors' own no-control model for context.
cat("==== Reproduction (no controls, pchalally==0) ====\n")
RHS_full <- RHS; RHS <- "dispute ~ ptargdef"; r_nc <- fitter("", d); RHS <- RHS_full
cat("==== Reproduction (controls, pchalally==0) ====\n")
r_pub <- fitter("", d)
n_orig <- r_pub$n

# 3. Merge history (challenger, target); restrict to 1822-2011. CEM weights kept as matched on the original data.
source(file.path(CODE, "hist/add_history.R"))
d <- add_history(d, "challenger", "target", "year")
n_pre <- sum(d$w > 0)
d <- filter(d, year >= 1822, year <= 2011)
cat("Weighted-sample obs dropped by 1822-2011 restriction:", n_pre - sum(d$w > 0), "\n")
cat("\nNA counts in history/benchmark vars (w>0):\n"); print(colSums(is.na(d[d$w > 0, HIST_VARS])))
d <- d[complete.cases(d[, c("dispute","ptargdef", CEMV, HIST_VARS)]), ]
cat("Rows (w>0) after year restriction + complete cases:", sum(d$w > 0), " (orig model N", n_orig, ")\n")

# 4. Models on the same sample
R <- run_six(fitter, d)

# ---- sensitivity hook (split only): runs only under run_sensitivity.R (SWEEP_SENS=1), never in the main sweep ----
# split only: the model is CEM-weighted, and the placebo / partial R2 steps of sens_run() are unweighted
if (Sys.getenv("SWEEP_SENS") == "1") {
  source(file.path(CODE, "hist/sens_helper.R"))
  sens_split(id = "14", dat = d,
             fit = function(dat, add) { r <- fitter(if (length(add)) paste("+", paste(add, collapse = " + ")) else "", dat); c(b = r$b, se = r$se, p = r$p) })
}

# Context: same exercise for initiators WITH a defensive ally (pchalally==1)
cat("\n==== Context: pchalally==1 ====\n")
d1$w <- cem_exact(d1, "ptargdef"); d1 <- d1[d1$w > 0, ]
d1 <- add_history(d1, "challenger", "target", "year"); d1 <- filter(d1, year >= 1822, year <= 2011)
d1 <- d1[complete.cases(d1[, c("dispute","ptargdef", CEMV, HIST_VARS)]), ]
c0 <- fitter("", d1); cDM <- fitter(" + hist_dyad_w15 + hist_mono_w15", d1)

# 5. Write row
write_row(14, "Morrow 2017", "UNVERIFIED", NA, r_pub$b, "CEM (exact) weighted logit (fixest feglm, iweights), ML SE",
  "ptargdef (target has defensive ally)", n_orig, R$m0$n, R, ahc = "Y",
  notes = sprintf(paste0("Fig.2 pchalally==0 w/ controls; published coef not obtainable; CEM rebuilt as exact matching (#0) with IKP weights; ",
    "author already matches+controls on rivalyears & midcount (MULTIPLE); pchalally==1 analog: b0 %.3f p %.3f -> b_dm %.3f p %.3f; dropped %d obs vs orig"),
    c0$b, c0$p, cDM$b, cDM$p, n_orig - R$m0$n))
sink()
