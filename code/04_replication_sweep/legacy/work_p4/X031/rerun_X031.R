# X031 Kuo & Blankenship 2022 JCR "Deterrence and Restraint" (10.1177/00220027211023147)
# KEY (fixed from abstract before estimating): jally_new2 = JME conducted with an ally (vs jnally2 = non-ally JME).
#   Evidence: "JMEs conducted with allies in particular reduce the probability of conflict escalation"
# Headline: Table 2 Model 2 (target escalation) and Model 4 (participant escalation), authors' R code:
#   stage 1 logit onset ~ minpow + rival + contig_land; resid = glm working residuals; stage 2 logit esc ~ jnally2 + jally_new2 +
#   atop_sy2 + jointdemo + cinc + unga + dy_trade + escalate_l + resid + factor(year), default glm SEs. Directed dyad-years 1973-2003.
# CONTROL (DEMOCRACY): jointdemo (target model). History added to stage 2 only; stage-1 residual kept from the authors' stage 1.
ROOT <- Sys.getenv("SWEEP_DATA"); OUT <- Sys.getenv("SWEEP_OUT"); CODE <- Sys.getenv("SWEEP_CODE")
stopifnot(nzchar(ROOT), nzchar(OUT), nzchar(CODE))
source(file.path(CODE, "work_p4/b7/common.R"))
WO <- file.path(OUT, "work_p4/X031"); dir.create(WO, recursive = TRUE, showWarnings = FALSE)  # port: outputs
W <- file.path(ROOT, "work_p4/X031"); logf <- file.path(WO, "log_X031.txt"); cat("", file = logf)
meta <- list(xid = "X031", cite = "Kuo & Blankenship 2022", journal = "JCR", year = 2022, doi = "10.1177/00220027211023147")
est <- "two-stage: logit 2nd stage with 1st-stage working residual + year dummies, iid SEs (fixest feglm)"
ahc <- "Y (lagged escalation escalate_l; 1st-stage rivalry)"
run_one <- function(file, dv, tag, terms) {
  keep <- c("onset","minpow","rival","contig_land", dv, "jnally2","jally_new2","atop_sy2","jointdemo","cinc","unga","dy_trade","escalate_l","year","ccode1","ccode2")
  d <- read.csv(file.path(W, file))[, keep]
  first <- glm(onset ~ minpow + rival + contig_land, data = d, family = binomial("logit"), control = glm.control(maxit = 200))
  stopifnot(length(first$residuals) == nrow(d)); d$resid <- first$residuals; rm(first); gc()
  d$fyear <- factor(d$year)
  rhs <- c("jnally2","jally_new2","atop_sy2","jointdemo","cinc","unga","dy_trade","escalate_l","resid","fyear")
  fitter <- function(r, dat) { m <- fit_bin(r, dv, dat); list(ct = coeftable(m), n = nobs(m)) }
  d0 <- d[complete.cases(d[, c(dv, rhs)]), ]
  rep <- fitter(rhs, d0); cat("\n####", tag, "REPRO N =", rep$n, "\n"); print(round(rep$ct[!grepl("fyear", rownames(rep$ct)), ], 4))
  cat("\n####", tag, "REPRO\n", file = logf, append = TRUE); capture.output(print(rep$ct), file = logf, append = TRUE)
  s <- common_sample(d, c(dv, rhs), "ccode1", "ccode2", "year")
  s$fyear <- droplevels(s$fyear)
  cat("\n####", tag, "SIX MODELS\n", file = logf, append = TRUE)
  R <- run_six(s, rhs, fitter, terms, logf)
  list(rep = rep, R = R, n_orig = nrow(d0))
}
T <- run_one("tgt/main target data.csv", "esc_s1", "TARGET", c("jally_new2","jnally2","jointdemo"))
nt <- function(x, lab) sprintf("%s escalation, Table 2 %s. Published tables not accessible (only 6-page preprint front matter online) -> UNVERIFIED. repro N=%d; common sample N=%d. jnally2 (non-ally JME) in model: m0 b=%.3f. Table 3 Heckman version INFEASIBLE (no sampleSelection).",
  lab, ifelse(lab == "Target", "Model 2", "Model 4"), x$n_orig, attr(x$R, "n"), x$R$m0$jnally2["b"])
write_p4(meta, "jally_new2", "TIE", "KEY", "abstract: 'JMEs conducted with allies in particular reduce the probability of conflict escalation'",
         "UNVERIFIED", NA, unname(kstat(T$rep$ct, "jally_new2")["b"]), est, T$n_orig, T$R, ahc, nt(T, "Target"), fname_term = "jally_new2_target")
write_p4(meta, "jointdemo", "DEMOCRACY", "CONTROL", "Table 2 Model 2 control (joint democracy)",
         "UNVERIFIED", NA, unname(kstat(T$rep$ct, "jointdemo")["b"]), est, T$n_orig, T$R, ahc, paste(nt(T, "Target"), "bD NOT a placebo for this term."))
rm(T); gc()
P <- run_one("par/main participant data.csv", "esc_s2", "PARTICIPANT", c("jally_new2","jnally2"))
write_p4(meta, "jally_new2", "TIE", "KEY", "abstract: 'JMEs conducted with allies in particular reduce the probability of conflict escalation'",
         "UNVERIFIED", NA, unname(kstat(P$rep$ct, "jally_new2")["b"]), est, P$n_orig, P$R, ahc, nt(P, "Participant"), fname_term = "jally_new2_participant")
