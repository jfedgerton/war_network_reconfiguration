# Rerun 57: Powers 2006 JPR -- REI (African regional economic institution) alliance membership -> MID initiation
# Headline: Model 1: xtlogit (random-effects logit) cwinit lagdyadreimem lagdyadreiall lagdemdyad caprat anyally
#   lagdtrade distance swtreg cwpceyrs _spline1-3 ; directed dyad-years (Africa), 1950-2000.
# Tie = lagreiall (dyad REI-with-alliance membership). b_pub = -1.184 (authors' do-file, predicted-probability
# code reproduces Model 1 estimates: reiall -1.184, reimem -.0606, lagdemdyad 1.108, cons -4.223).
# Data names differ from the do-file (lagdyadreiall -> lagreiall, swtreg -> s_wt_reg, anyally -> laganyally);
# 'caprat' is ambiguous (lagcaprat vs capratio) and the xtlogit panel id is not in the do-file, so the
# candidates are tried and the one matching the author's estimates is used (identification, not tuning).
# Estimator: lme4::glmer random intercept, adaptive Gauss-Hermite nAGQ = 12 (Stata xtlogit re default 12 pts).
ROOT <- Sys.getenv("SWEEP_DATA"); OUT <- Sys.getenv("SWEEP_OUT"); CODE <- Sys.getenv("SWEEP_CODE")
stopifnot(nzchar(ROOT), nzchar(OUT), nzchar(CODE))
source(file.path(CODE, "work/rerun_helpers_s68.R"))
WO <- file.path(OUT, "work/57"); dir.create(WO, recursive = TRUE, showWarnings = FALSE)  # port: outputs
W <- file.path(ROOT, "work/57"); sink(file.path(WO, "log_57.txt"), split = TRUE)
suppressPackageStartupMessages({library(haven); library(lme4)})

d <- as.data.frame(zap_labels(read_dta(file.path(W, "powersjprdata/powersjprdata.dta"))))
names(d) <- sub("^_", "x_", names(d))
stopifnot(nrow(d) == 90722)
d$ddyad <- d$ccode1 * 1000 + d$ccode2                                   # directed dyad
d$udyad <- pmin(d$ccode1, d$ccode2) * 1000 + pmax(d$ccode1, d$ccode2)   # undirected dyad
TIE <- "lagreiall"

fit_re <- function(rhs, dat, grp, label) {
  f <- as.formula(paste("cwinit ~", paste(rhs, collapse = " + "), "+ (1|", grp, ")"))
  dat <- dat[complete.cases(dat[, c("cwinit", rhs, grp)]), ]
  # rescale nothing (keep author's scale); use bobyqa with more evaluations
  m <- glmer(f, data = dat, family = binomial, nAGQ = 12,
             control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1e5)))
  ct <- summary(m)$coefficients
  cat("\n---- [", label, "] ", deparse(f, width.cutoff = 500), "| N =", nobs(m), " groups =", ngrps(m),
      " sd_u =", round(sqrt(unlist(VarCorr(m))), 4), " ll =", round(as.numeric(logLik(m)), 3), "\n")
  print(round(ct, 5))
  msgs <- m@optinfo$conv$lme4$messages; if (length(msgs)) cat("CONVERGENCE NOTE:", msgs, "\n")
  list(b = unname(ct[TIE, 1]), se = unname(ct[TIE, 2]), p = unname(ct[TIE, 4]), n = nobs(m), ct = ct)
}

base_ctrl <- function(cap) c("lagreimem", "lagdemdyad", cap, "laganyally", "lagdtrade", "distance", "s_wt_reg",
                             "cwpceyrs", "x_spline1", "x_spline2", "x_spline3")
b_pub <- -1.184
pub_other <- c(lagreimem = -0.0606, lagdemdyad = 1.108, laganyally = 1.012, lagdtrade = 0.0059, s_wt_reg = -0.5895, cwpceyrs = -0.1393)

# 2. Reproduce: try candidate specs, pick the one closest to the author's full coefficient vector
cands <- expand.grid(cap = c("lagcaprat", "capratio"), grp = c("ddyad", "udyad"), stringsAsFactors = FALSE)
R <- list(); dist <- numeric(nrow(cands))
for (i in seq_len(nrow(cands))) {
  R[[i]] <- fit_re(c(TIE, base_ctrl(cands$cap[i])), d, cands$grp[i], paste("repro cand", i, cands$cap[i], cands$grp[i]))
  est <- c(R[[i]]$b, R[[i]]$ct[names(pub_other), 1])
  dist[i] <- sum(abs(est - c(b_pub, pub_other)) / abs(c(b_pub, pub_other)))
  cat("cand", i, "tie b =", R[[i]]$b, " sum rel dev over 7 coefs =", dist[i], "\n")
}
best <- which.min(dist); cap <- cands$cap[best]; grp <- cands$grp[best]
rep <- R[[best]]; n_orig <- rep$n
cat("Chosen:", cap, grp, " b_repro =", rep$b, " rel diff =", abs(rep$b - b_pub) / abs(b_pub), "\n")
# significance of pub unknown from do-file; assume published as significant (key finding) -> require p<.05
status <- if (abs(rep$b - b_pub) / abs(b_pub) <= 0.10) "REPRODUCED" else "NOT_REPRODUCED"
rhs0 <- c(TIE, base_ctrl(cap))

# 3. History (directed dyads)
d <- add_history(d, "ccode1", "ccode2", "year")
s <- d[d$year >= 1822 & d$year <= 2011, ]
s <- s[complete.cases(s[, c("cwinit", rhs0, grp, "hist_dyad_w15","hist_mono_w15","hist_dyad_w26","hist_mono_w26")]), ]
cat("Estimation sample", nrow(s), "; dropped vs original", n_orig - nrow(s), "\n")

# 4. Models
m0  <- fit_re(rhs0, s, grp, "m0")
mD  <- fit_re(c(rhs0, "hist_dyad_w15"), s, grp, "mD")
mDM <- fit_re(c(rhs0, "hist_dyad_w15", "hist_mono_w15"), s, grp, "mDM")
m26 <- fit_re(c(rhs0, "hist_dyad_w26", "hist_mono_w26"), s, grp, "m26")
sC <- s[!is.na(s$bm_cinc_w15), ]; sD <- s[!is.na(s$bm_democ_w15), ]
cat("bC sample", nrow(sC), " bD sample", nrow(sD), "\n")
b0C <- if (nrow(sC) < nrow(s)) fit_re(rhs0, sC, grp, "m0 on bC sample")$b else NULL
bC  <- fit_re(c(rhs0, "bm_cinc_w15"), sC, grp, "bC")
b0D <- if (nrow(sD) < nrow(s)) fit_re(rhs0, sD, grp, "m0 on bD sample")$b else NULL
bD  <- fit_re(c(rhs0, "bm_democ_w15"), sD, grp, "bD")

write_row(57, "Powers 2006", status, b_pub, rep$b, sprintf("RE logit (glmer, AGHQ 12), random intercept %s", grp), "lagreiall",
  n_orig, m0$n, m0, mD, mDM, m26, bC, bD, b0C, b0D, ahc = "Y",
  notes = sprintf("Model 1; b_pub from author's do-file (pp code). Spec identification: caprat=%s, panel=%s (4 candidates tried; see log). Dropped %d obs. Rare DV (~0.26%%).",
                  cap, grp, n_orig - m0$n))
sink()
