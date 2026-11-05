# ESTAUDIT X082 Green, Kim & Yoon 2001 IO. Copy of rerun_verify_X082.R (latest; MID-presence DV) with ONE change:
#   estimator unconditional dyad-FE logit (fixest feglm, ff_felogit_iid) -> conditional (fixed-effects) logit = Stata clogit,
#   survival::clogit(method = "exact", strata = dyadid), model-based (OIM) SEs, as in the paper's Table 3 col 2 ("fixed-effects logit").
#   Same data prep, same >=20-obs rule, same common sample for all six steps. New cache dir work_p4/X082/est/.
ROOT <- Sys.getenv("SWEEP_DATA"); OUT <- Sys.getenv("SWEEP_OUT"); CODE <- Sys.getenv("SWEEP_CODE")
stopifnot(nzchar(ROOT), nzchar(OUT), nzchar(CODE))
source(file.path(CODE, "work_p4/common_p4_b02.R"))
source(file.path(CODE, "ea2/est_helpers.R"))
suppressPackageStartupMessages({library(peacesciencer); library(survival)})
WO <- file.path(OUT, "work_p4/X082/est"); dir.create(WO, recursive = TRUE, showWarnings = FALSE)  # port: outputs
W <- file.path(ROOT, "work_p4/X082/est"); sink(file.path(OUT, "work_p4/X082/log_est_X082.txt"), split = TRUE)
d <- readRDS(file.path(ROOT, "work_p4/X074/or_trade.rds"))
d <- d[d$year >= 1950 & d$year <= 1991, c("statea", "stateb", "year", "dyadid", "mzmid1_all", "contigland", "lcaprat", "lowgrowth",
                                         "anyally", "lowdem", "lowdep")]
stopifnot(all(d$statea < d$stateb))
cm <- as.data.frame(cow_mid_ddydisps)[, c("ccode1", "ccode2", "year", "cowmidongoing")]
cm$c1 <- pmin(cm$ccode1, cm$ccode2); cm$c2 <- pmax(cm$ccode1, cm$ccode2)
pres <- aggregate(cowmidongoing ~ c1 + c2 + year, data = cm, FUN = max)
names(pres) <- c("statea", "stateb", "yr1", "mid_pres1")
d$yr1 <- d$year + 1
d <- merge(d, pres, by = c("statea", "stateb", "yr1"), all.x = TRUE)
d$mid_pres1[is.na(d$mid_pres1)] <- 0
f0 <- mid_pres1 ~ contigland + lcaprat + lowgrowth + anyally + lowdem + lowdep
d <- d[complete.cases(d[, all.vars(f0)]), ]
d$nobs <- ave(d$year, d$dyadid, FUN = length); d <- d[d$nobs >= 20, ]
cat("complete-case dyad-years in dyads with >=20 obs:", nrow(d), "(published 93,755) dyads:", length(unique(d$dyadid)), "\n")
K <- c("lowdep", "lowdem")
# conditional logit (Stata clogit): exact conditional likelihood within dyad strata; OIM SEs. N = informative obs (strata with outcome variation).
ff_clogit <- function(f, data) {
  y <- all.vars(f)[1]
  v <- ave(data[[y]], data$dyadid, FUN = function(z) as.numeric(min(z) != max(z)))
  data <- data[v == 1, ]
  f2 <- update(f, paste("~ . + strata(dyadid)"))
  m <- clogit(f2, data = data, method = "exact")
  ct <- summary(m)$coefficients[, c(1, 3, 5), drop = FALSE]      # coef, se(coef), Pr(>|z|)
  stopifnot(all(is.finite(ct[, 1])))
  attr(ct, "n") <- nrow(data); cat("   [clogit informative N =", nrow(data), " strata =", length(unique(data$dyadid)), " loglik =", round(m$loglik[2], 3), "]\n"); ct
}
r0 <- ff_clogit(f0, d); print(round(r0, 5)); n_orig <- attr(r0, "n")
rF <- ff_felogit_iid("dyadid")(f0, d); cat("for comparison, previous unconditional FE logit on same data:\n"); print(round(rF, 5))
ds <- common_sample(d, f0, "dyadid", "statea", "stateb", "year"); rm(d); gc()
M <- run_six(f0, ds, K, ff_clogit, WO)
ns <- sapply(M, `[[`, "n"); print(ns)

# ---- sensitivity hook: runs only under run_sensitivity.R (SWEEP_SENS=1), never in the main sweep ----
if (Sys.getenv("SWEEP_SENS") == "1") {
  source(file.path(CODE, "hist/sens_helper.R"))
  # ff_clogit() keeps only dyads with outcome variation (the informative sample, same for every step); pass exactly that sample
  sens_ds <- ds[ave(ds$mid_pres1, ds$dyadid, FUN = function(z) as.numeric(min(z) != max(z))) == 1, ]
  sens_run(id = "X082", dat = sens_ds, y = "mid_pres1", tie = "lowdep",
           ctrl = c("contigland", "lcaprat", "lowgrowth", "anyally", "lowdem"), fe = "dyadid",
           fit = function(dat, add) { ct <- ff_clogit(if (length(add)) update(f0, as.formula(paste(". ~ . +", paste(add, collapse = " + ")))) else f0, dat)
             c(b = unname(ct["lowdep", 1]), se = unname(ct["lowdep", 2]), p = unname(ct["lowdep", 3])) })
}

AE <- "conditional (dyad fixed-effects) logit, Stata clogit, OIM SEs; dyads with >=20 obs; no temporal terms (Table 3 col 2)"
ES <- "conditional logit (survival::clogit exact, strata dyadid), OIM SEs"
MM <- "unconditional dyad-FE logit (fixest feglm dummies) instead of conditional logit in all six steps"
for (k in K) {
  R <- lapply(M, function(z) z$k[[k]])
  write_est("X082", "Green, Kim & Yoon 2001", k, AE, "pdf (paper.txt Table 3 col 2 'fixed effects' logit; clogit)", ES, "N", MM, "Y", R, M$m0$n,
    paste0(if (k == "lowdep") "KEY" else "CONTROL (DEMOCRACY)", ". REBUILT data (verify spec: DV = COW MID presence t+1). Full rebuilt sample clogit b=", signif(r0[k, 1], 5),
      " (p=", signif(r0[k, 3], 3), "), informative N=", n_orig, "; previous unconditional FE-logit b=", signif(rF[k, 1], 5),
      ". Common-sample N (informative) per step: ", paste(ns, collapse = "/"), ". lowdep scale not comparable to paper (fraction vs unstated units)."))
}
sink()
