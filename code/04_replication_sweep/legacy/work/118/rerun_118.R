# Rerun 118: Hafner-Burton & Montgomery 2006 JCR "Power Positions"
# Headline: Table 1 col 2, logit disp_l1 (MID involvement t+1), cluster(dyadid).
# Headline tie: CentDif (difference in IGO-network prestige); ClusSame and ClusSizeMax also logged.
ROOT <- Sys.getenv("SWEEP_DATA"); OUT <- Sys.getenv("SWEEP_OUT"); CODE <- Sys.getenv("SWEEP_CODE")
stopifnot(nzchar(ROOT), nzchar(OUT), nzchar(CODE))
source(file.path(CODE, "work/rerun_helpers_b118.R"))
library(haven)
WO <- file.path(OUT, "work/118"); dir.create(WO, recursive = TRUE, showWarnings = FALSE)  # port: outputs
W <- file.path(ROOT, "work/118"); logf <- file.path(WO, "log_118.txt"); cat("", file = logf)
sink(logf, split = TRUE, append = TRUE)

# ---- 1. Load authors' data and apply IGODataSetup.do exactly (Stata missing semantics) ----
d <- as.data.frame(read_dta(file.path(W, "IGOData.dta")))
n_orig <- nrow(d); stopifnot(n_orig == 149403)
d$CentDif <- abs(d$Cent1 - d$Cent2)
# Stata: gen ClusSame=1 if Clus1==Clus2 (missing==missing is TRUE in Stata); else 0
d$ClusSame <- as.numeric((is.na(d$Clus1) & is.na(d$Clus2)) | (!is.na(d$Clus1) & !is.na(d$Clus2) & d$Clus1 == d$Clus2))
# Stata max() ignores missing
d$ClusSizeMax <- suppressWarnings(pmax(d$ClusSize1, d$ClusSize2, na.rm = TRUE))
d$smldmat <- pmin(d$demauta, d$demautb)   # NA if either missing (matches Stata if-conditions)
d$smldep  <- pmin(d$dependa, d$dependb)
sa <- d$statea; sb <- d$stateb; y <- d$year
mp <- rep(0, nrow(d))
mp[(sa == 2 | sb == 2) & y > 1898] <- 1
mp[(sa == 200 | sb == 200)] <- 1
mp[(sa == 220 | sb == 220) & (y < 1941 | y > 1944)] <- 1
mp[(sa == 255 | sb == 255) & (y < 1919 | (y > 1924 & y < 1946) | y > 1990)] <- 1
mp[(sa == 300 | sb == 300) & y < 1919] <- 1
mp[(sa == 325 | sb == 325) & (y > 1859 & y < 1944)] <- 1
mp[(sa == 365 | sb == 365) & (y < 1918 | y > 1922)] <- 1
mp[(sa == 710 | sb == 710) & y > 1949] <- 1
mp[(sa == 740 | sb == 740) & ((y > 1894 & y < 1946) | y > 1990)] <- 1
d$majpower <- mp
d$dyadid <- 1000 * d$statea + d$stateb

base_rhs <- c("IGOSame","ClusSame","CentDif","ClusSizeMax","smldmat","smldep","lcaprat2","allies","hegdefb",
              "contigkb","logdstab","majpower","disp_spl0","disp_spl1","disp_spl2","disp_spl3")
dv <- "disp_l1"; tie <- "CentDif"

# ---- 2. Reproduce Table 1 col 2 on the authors' full sample ----
m_rep <- fit_glm(base_rhs, dv, d, "logit", "dyadid")
cat("\n==== REPRODUCTION (full sample) ====\n"); print(coeftable(m_rep)); cat("N =", nobs(m_rep), "\n")
b_pub <- -.0003548  # IGOTables.log Table 1 col 2 (se .0001031, p .001); N=149403
b_repro <- unname(coef(m_rep)[tie]); p_repro <- coeftable(m_rep)[tie, 4]
stopifnot(nobs(m_rep) == 149403)
rel <- abs(b_repro - b_pub) / abs(b_pub)
reproduced <- ifelse(rel <= 0.10 & p_repro < .05, "REPRODUCED", "NOT_REPRODUCED")
cat(sprintf("b_pub=%g b_repro=%g rel.diff=%.4f p_repro=%.4g -> %s\n", b_pub, b_repro, rel, p_repro, reproduced))
cat("ClusSame pub -.1605983 repro", coef(m_rep)["ClusSame"], "; ClusSizeMax pub .0116839 repro", coef(m_rep)["ClusSizeMax"], "\n")

# ---- 3. Merge history. DV disp_l1 is MID involvement at t+1, so history and benchmark windows are keyed on the outcome year t+1
#      (fix 2026-10-10; previously the row year t). year_h (1822-2011 filter) stays the row year: 1885-1991, unaffected. ----
# lead check: the authors' file has no same-year involvement variable, so check against COW MID 5.0 dyadic involvement
# (peacesciencer cow_mid_dirdisps): dyad-years with a COW MID in t+1 must almost all carry disp_l1 = 1, and far more than in t.
mk <- as.data.frame(peacesciencer::cow_mid_dirdisps)
mk <- unique(paste(pmin(mk$ccode1, mk$ccode2), pmax(mk$ccode1, mk$ccode2), mk$year))
ua <- pmin(d$statea, d$stateb); ub <- pmax(d$statea, d$stateb)
i1 <- paste(ua, ub, d$year + 1) %in% mk & d$year < max(d$year); i0 <- paste(ua, ub, d$year) %in% mk
s1 <- mean(d$disp_l1[i1] == 1); s0 <- mean(d$disp_l1[i0] == 1)
cat(sprintf("lead check: share of COW MID dyad-years with disp_l1 = 1: MID in t+1 %.4f (n %d), MID in t %.4f (n %d)\n", s1, sum(i1), s0, sum(i0)))
stopifnot(sum(i1) > 500, s1 > 0.95, s1 > s0 + 0.3); rm(mk, ua, ub, i1, i0, s1, s0)
d$outcome_year <- d$year + 1
d <- add_history(d, "statea", "stateb", "outcome_year"); d$year_h <- d$year
# ---- 4. Battery on common restricted sample ----
ds <- battery_sample(d, c(dv, base_rhs, "dyadid"))
fitter <- function(rhs, data) fit_glm(rhs, dv, data, "logit", "dyadid")
R <- run_battery(ds, base_rhs, fitter, function(m) tie_stats(m, tie), logf)
# also log the other network terms under m0 and mDM
m0 <- fitter(base_rhs, ds); mdm <- fitter(c(base_rhs, "hist_dyad_w15", "hist_mono_w15"), ds)
cat("\nOther network terms m0 vs mDM:\n"); print(cbind(m0 = coef(m0)[c("ClusSame","ClusSizeMax","IGOSame")],
                                                  mDM = coef(mdm)[c("ClusSame","ClusSizeMax","IGOSame")]))
print(coeftable(mdm)[c("ClusSame","ClusSizeMax","IGOSame","hist_dyad_w15","hist_mono_w15"), ])

# ---- sensitivity hook: runs only under run_sensitivity.R (SWEEP_SENS=1), never in the main sweep ----
if (Sys.getenv("SWEEP_SENS") == "1") {
  source(file.path(CODE, "hist/sens_helper.R"))
  sens_run(id = "118", dat = ds, y = dv, tie = tie, ctrl = setdiff(base_rhs, tie), fe = NULL,
           fit = function(dat, add) tie_stats(fitter(c(base_rhs, add), dat), tie))
}

# ---- 5. Write result ----
notes <- sprintf("Headline tie CentDif (IGO prestige difference) from Tbl1 col2; exact reproduction (rel.diff %.4f). Years 1885-1991 all within GML range; drops only from benchmark missingness. DV is disp at t+1; history and benchmark windows keyed on outcome year t+1. mDM: ClusSame %.3f->%.3f, ClusSizeMax %.4f->%.4f, IGOSame %.4f->%.4f (sign flips, ns).",
                 rel, coef(m0)["ClusSame"], coef(mdm)["ClusSame"], coef(m0)["ClusSizeMax"], coef(mdm)["ClusSizeMax"], coef(m0)["IGOSame"], coef(mdm)["IGOSame"])
write_result(118, "Hafner-Burton & Montgomery 2006", reproduced, b_pub, b_repro, "logit, cluster dyadid", tie,
             n_orig, nrow(ds), R, "N", notes)
sink()
