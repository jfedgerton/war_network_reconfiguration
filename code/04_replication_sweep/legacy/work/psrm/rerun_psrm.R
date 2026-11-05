# Edgerton (PSRM, forthcoming) hierarchy misalignment: Table 3 col 3 + conflict history
suppressPackageStartupMessages({library(dplyr); library(fixest)})
ROOT <- Sys.getenv("SWEEP_DATA"); OUT <- Sys.getenv("SWEEP_OUT"); CODE <- Sys.getenv("SWEEP_CODE")
stopifnot(nzchar(ROOT), nzchar(OUT), nzchar(CODE))
source(file.path(CODE, "hist/add_history.R"))
t <- readRDS(file.path(ROOT, "work/psrm", "target_dyads_for_analysis.rds"))
t <- t %>% select(ccode1, ccode2, year, cinc1, cinc2, fatal_init, joint_democ, joint_autoc, capdist, cinc_pr_dff)
lag_iv <- t %>% select(ccode1, ccode2, year, cinc_pr_dff) %>% mutate(year = year + 1)
d <- t %>% select(-cinc_pr_dff) %>% left_join(lag_iv, by = c("ccode1","ccode2","year")) %>%
  mutate(cinc_ratio = ifelse(cinc1 > cinc2, cinc1/(cinc1+cinc2), cinc2/(cinc1+cinc2)))
rm(t, lag_iv); gc()
f0 <- fatal_init ~ cinc_pr_dff + cinc_ratio + joint_democ + joint_autoc + capdist | year + ccode1 + ccode2
m_pub <- feglm(f0, data = d, family = "binomial")
cat("REPRO:", coef(m_pub)["cinc_pr_dff"], se(m_pub)["cinc_pr_dff"], nobs(m_pub), "(pub 0.94, 0.21, 789510)\n")
d <- add_history(d, "ccode1", "ccode2", "year")
s <- d %>% filter(year >= 1822, year <= 2011, !is.na(bm_cinc_w15), !is.na(bm_democ_w15))
fx <- function(x) as.formula(paste("fatal_init ~ cinc_pr_dff + cinc_ratio + joint_democ + joint_autoc + capdist", x, "| year + ccode1 + ccode2"))
m0 <- feglm(fx(""), s, family="binomial"); mD <- feglm(fx("+ hist_dyad_w15"), s, family="binomial")
mDM <- feglm(fx("+ hist_dyad_w15 + hist_mono_w15"), s, family="binomial"); m26 <- feglm(fx("+ hist_dyad_w26 + hist_mono_w26"), s, family="binomial")
bC <- feglm(fx("+ bm_cinc_w15"), s, family="binomial"); bD <- feglm(fx("+ bm_democ_w15"), s, family="binomial")
etable(m0, mD, mDM, m26, bC, bD, keep = c("cinc_pr_dff","hist","bm_"))
# ---- sensitivity hook: runs only under run_sensitivity.R (SWEEP_SENS=1), never in the main sweep ----
if (Sys.getenv("SWEEP_SENS") == "1") {
  source(file.path(CODE, "hist/sens_helper.R"))
  sens_run(id = "N06", dat = s[complete.cases(s[, c("fatal_init", "cinc_pr_dff", "cinc_ratio", "joint_democ", "joint_autoc", "capdist", "year", "ccode1", "ccode2")]), ],
           y = "fatal_init", tie = "cinc_pr_dff", ctrl = c("cinc_ratio", "joint_democ", "joint_autoc", "capdist"), fe = c("year", "ccode1", "ccode2"), fe_sep = TRUE,
           fit = function(dat, add) { m <- feglm(fx(if (length(add)) paste("+", paste(add, collapse = " + ")) else ""), dat, family="binomial")
             c(b = unname(coef(m)["cinc_pr_dff"]), se = unname(se(m)["cinc_pr_dff"]), p = unname(pvalue(m)["cinc_pr_dff"])) })
}
k <- "cinc_pr_dff"; b0 <- coef(m0)[k]; pc <- function(m) round(100*(coef(m)[k]-b0)/abs(b0),1)
out <- data.frame(cid="PSRM", short_cite="Edgerton (forthcoming PSRM)", reproduced=ifelse(abs(coef(m_pub)[k]-0.94)/0.94<=.1,"REPRODUCED","NOT_REPRODUCED"),
  b_pub=0.94, b_repro=coef(m_pub)[k], estimator="FE logit (year, ccode1, ccode2), fixest default SE", tie_var=k, n_orig=nobs(m_pub), n_used=nobs(m0),
  b0=b0, se0=se(m0)[k], p0=pvalue(m0)[k], b_dm=coef(mDM)[k], se_dm=se(mDM)[k], p_dm=pvalue(mDM)[k], pct_dm=pc(mDM), b_d=coef(mD)[k], pct_d=pc(mD),
  b_26=coef(m26)[k], se_26=se(m26)[k], p_26=pvalue(m26)[k], pct_26=pc(m26), pct_bC=pc(bC), pct_bD=pc(bD), author_has_history_control="N", notes="Table 3 col 3; IV lagged one year as in author code")
dir.create(file.path(OUT, "results"), recursive = TRUE, showWarnings = FALSE)
write.csv(out, file.path(OUT, "results/PSRM.csv"), row.names=FALSE); print(out[, names(out) != "short_cite", drop = FALSE])
